#!/usr/bin/env python3
"""ClawBio skill: rank pathogen proteins as targets by host dissimilarity.

Takes a pathogen proteome, a host proteome and optional panel proteomes, and
returns a ranked target list scored on divergence from the host, published
essentiality, and conservation across the panel.

Runs offline. Search is done in-process with phmmer, or read from a
BLAST/DIAMOND tabular file for full-proteome runs done elsewhere.
"""

import argparse
import csv
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)
from homology_search import (  # noqa: E402
    DEFAULT_EVALUE,
    HomologyHit,
    best_hit_per_query,
    filter_by_coverage,
    parse_tabular_hits,
    parse_uniprot_header,
    read_fasta,
    read_fasta_with_headers,
    search_pyhmmer,
)
from target_scoring import TargetVerdict, score_target  # noqa: E402

DEMO_DIR = Path(__file__).resolve().parent / "demo"
DEFAULT_MIN_COVERAGE = 50.0
DISCLAIMER = (
    "*ClawBio is a research tool. Not a medical device. Scores rank the evidence "
    "available for each protein, not the quality of the biology; a low score on a "
    "protein with blank annotations means the work has not been done.*"
)
ANNOTATION_COLUMNS = ("protein_id", "gene_name", "description", "essentiality", "essentiality_source")


def read_annotations(path: Path) -> dict[str, dict[str, str]]:
    """Read the per-protein annotation table.

    Args:
        path: CSV with columns protein_id, gene_name, description,
            essentiality, essentiality_source.

    Returns:
        Mapping of protein identifier to its annotation fields.

    Raises:
        SystemExit: If the file is missing or lacks required columns.
    """
    if not path.exists():
        raise SystemExit(f"Annotation file not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    missing = [column for column in ANNOTATION_COLUMNS if column not in fieldnames]
    if missing:
        raise SystemExit(
            f"{path.name} is missing required columns: {', '.join(missing)}. "
            f"Found: {', '.join(fieldnames)}."
        )
    return {row["protein_id"]: row for row in rows if row.get("protein_id")}


def annotations_from_fasta(path: Path) -> dict[str, dict[str, str]]:
    """Derive gene names and descriptions from FASTA headers.

    Used when no annotation CSV is supplied, so that a run against a real
    UniProt proteome shows gene names rather than bare accessions.
    Essentiality is left empty, because a FASTA header carries no such
    evidence and inventing one would defeat the scoring.

    Args:
        path: Path to the pathogen proteome FASTA.

    Returns:
        Mapping of protein identifier to annotation fields.
    """
    annotations: dict[str, dict[str, str]] = {}
    for identifier, _, header in read_fasta_with_headers(path):
        parsed = parse_uniprot_header(header)
        annotations[identifier] = {
            "protein_id": identifier,
            "gene_name": parsed["gene_name"],
            "description": parsed["description"],
            "essentiality": "",
            "essentiality_source": "",
        }
    return annotations


def merge_annotations(
    base: dict[str, dict[str, str]], override: dict[str, dict[str, str]]
) -> dict[str, dict[str, str]]:
    """Merge curated annotations over header-derived ones.

    A non-empty field in the override wins. An empty field falls back to the
    header-derived value, so a CSV carrying only essentiality calls does not
    wipe out gene names parsed from the FASTA.

    Args:
        base: Header-derived annotations.
        override: Annotations read from a user-supplied CSV.

    Returns:
        The merged mapping.
    """
    merged = {identifier: dict(fields) for identifier, fields in base.items()}
    for identifier, fields in override.items():
        target = merged.setdefault(identifier, {"protein_id": identifier})
        for key, value in fields.items():
            if value.strip():
                target[key] = value
            else:
                target.setdefault(key, "")
    return merged


def build_demo_inputs() -> dict[str, object]:
    """Load the committed synthetic demo proteomes and annotations.

    Returns:
        Keyword arguments for `run`: pathogen, host, panels, annotations and
        a source label.
    """
    panels = {
        path.stem.replace("demo_panel_", ""): read_fasta(path)
        for path in sorted(DEMO_DIR.glob("demo_panel_*.faa"))
    }
    return {
        "pathogen": read_fasta(DEMO_DIR / "demo_pathogen.faa"),
        "host": read_fasta(DEMO_DIR / "demo_host.faa"),
        "panels": panels,
        "annotations": read_annotations(DEMO_DIR / "demo_essentiality.csv"),
        "source_label": "built-in synthetic demo proteomes (10 pathogen proteins, 3 panel species)",
    }


def count_panel_presence(
    pathogen: list[tuple[str, str]],
    panels: dict[str, list[tuple[str, str]]],
    evalue: float,
    min_coverage: float,
) -> dict[str, int]:
    """Count how many panel species carry a homologue of each pathogen protein.

    Args:
        pathogen: Pathogen proteins.
        panels: Mapping of panel species name to its proteins.
        evalue: E-value threshold.
        min_coverage: Minimum query coverage for a hit to count.

    Returns:
        Mapping of pathogen protein identifier to the number of panel species
        with a homologue.
    """
    counts = {identifier: 0 for identifier, _ in pathogen}
    for proteins in panels.values():
        hits = filter_by_coverage(search_pyhmmer(pathogen, proteins, evalue), min_coverage)
        for identifier in best_hit_per_query(hits):
            counts[identifier] += 1
    return counts


def resolve_host_hits(
    pathogen: list[tuple[str, str]],
    host: list[tuple[str, str]],
    hits_file: Path | None,
    evalue: float,
    min_coverage: float,
) -> tuple[dict[str, HomologyHit], bool]:
    """Obtain the best host hit per pathogen protein.

    Args:
        pathogen: Pathogen proteins.
        host: Host proteins; may be empty when a hits file is supplied.
        hits_file: Precomputed BLAST/DIAMOND tabular file, or None.
        evalue: E-value threshold.
        min_coverage: Minimum query coverage.

    Returns:
        Tuple of the best-hit mapping and whether a host search was available
        at all. The flag matters: no hits from a search that ran means no host
        homologue, while no search means no information.
    """
    if hits_file is not None:
        hits = filter_by_coverage(parse_tabular_hits(hits_file), min_coverage)
        return best_hit_per_query(hits), True
    if not host:
        return {}, False
    hits = filter_by_coverage(search_pyhmmer(pathogen, host, evalue), min_coverage)
    return best_hit_per_query(hits), True


def write_target_table(path: Path, verdicts: list[TargetVerdict]) -> None:
    """Write the scored target table.

    Args:
        path: Destination CSV path.
        verdicts: Scored targets in output order.
    """
    columns = [
        "protein_id", "gene_name", "description", "score", "decision",
        "host_identity_pct", "panel_present", "panel_total",
        "flags", "missing_annotations", "reasons",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for verdict in verdicts:
            writer.writerow([
                verdict.protein_id,
                verdict.gene_name,
                verdict.description,
                verdict.score,
                verdict.decision,
                "" if verdict.host_identity_pct is None else verdict.host_identity_pct,
                verdict.panel_present,
                verdict.panel_total,
                "|".join(verdict.flags),
                "|".join(verdict.missing_annotations),
                "; ".join(verdict.reasons),
            ])


def write_host_hits_table(path: Path, hits: dict[str, HomologyHit]) -> None:
    """Write the best host hit per pathogen protein.

    Args:
        path: Destination CSV path.
        hits: Best host hit per pathogen protein; empty writes a header only.
    """
    columns = ["query_id", "host_subject_id", "identity_pct", "coverage_pct", "evalue"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for query_id in sorted(hits):
            hit = hits[query_id]
            writer.writerow([
                hit.query_id,
                hit.subject_id,
                hit.identity_pct,
                "" if hit.coverage_pct is None else hit.coverage_pct,
                hit.evalue,
            ])


def build_report(
    verdicts: list[TargetVerdict],
    source_label: str,
    panel_names: list[str],
    host_searched: bool,
) -> str:
    """Assemble the markdown report.

    Args:
        verdicts: Scored targets in rank order.
        source_label: Description of the inputs.
        panel_names: Panel species searched.
        host_searched: Whether a host search was available.

    Returns:
        Complete report text.
    """
    pursued = [verdict for verdict in verdicts if verdict.decision == "pursue"]
    no_host = [verdict for verdict in verdicts if verdict.host_identity_pct is None]
    incomplete = [verdict for verdict in verdicts if verdict.missing_annotations]

    lines = [
        "# Pathogen Target Dissimilarity Report",
        "",
        f"**Input**: {source_label}  ",
        f"**Panel**: {', '.join(panel_names) if panel_names else 'none supplied'}  ",
        f"**Date**: {date.today().isoformat()}",
        "",
        f"Scored **{len(verdicts)} proteins**; **{len(pursued)}** reach the pursue threshold; "
        f"**{len(no_host)}** have no detectable host homologue.",
        "",
        "| Target | Description | Score | Decision | Host identity | Panel | Unannotated |",
        "|---|---|---|---|---|---|---|",
    ]
    for verdict in verdicts:
        identity = "none" if verdict.host_identity_pct is None else f"{verdict.host_identity_pct:.0f}%"
        if not host_searched:
            identity = "not assessed"
        panel = f"{verdict.panel_present}/{verdict.panel_total}" if verdict.panel_total else "n/a"
        missing = ", ".join(verdict.missing_annotations) if verdict.missing_annotations else "none"
        lines.append(
            f"| {verdict.gene_name or verdict.protein_id} | {verdict.description} | {verdict.score} | "
            f"{verdict.decision} | {identity} | {panel} | {missing} |"
        )

    orphans = [verdict for verdict in verdicts if "orphan_candidate" in verdict.flags]
    top_score = verdicts[0].score if verdicts else 0
    tied_at_top = [verdict for verdict in verdicts if verdict.score == top_score]

    lines.extend(["", "## Interpretation", ""])
    if orphans:
        lines.append(
            f"{len(orphans)} proteins have no host homologue and none in any panel species. "
            "They forfeit the host-divergence points: absence everywhere means the protein "
            "is not conserved, which is not the same as being selectively targetable."
        )
        lines.append("")
    if len(tied_at_top) > 1:
        lines.append(
            f"{len(tied_at_top)} proteins tie at the top score of {top_score}. A tie this wide "
            "means the available axes cannot separate them; the order within the tie is "
            "arbitrary and should not be read as a ranking. Supplying essentiality data or a "
            "wider panel is what breaks it."
        )
        lines.append("")
    if not host_searched:
        lines.append(
            "No host proteome was searched, so every selectivity score is zero and the "
            "ranking reflects essentiality and conservation alone. Supply a host proteome "
            "or a precomputed hits file before using this ranking to choose a target."
        )
    else:
        lines.append(
            f"{len(pursued)} of {len(verdicts)} proteins clear the pursue threshold. "
            f"{len(incomplete)} carry at least one unannotated field and were scored "
            "conservatively; read the unannotated column before treating a low score as a "
            "verdict on the protein."
        )
    lines.extend([
        "",
        "Percent identity to the closest host protein is a first filter and not a "
        "selectivity guarantee. Two proteins at 40% identity can share a near-identical "
        "binding site, which this score does not see. Confirm at the pocket before "
        "committing to a target.",
        "",
        DISCLAIMER,
        "",
    ])
    return "\n".join(lines)


def run(
    pathogen: list[tuple[str, str]],
    host: list[tuple[str, str]],
    panels: dict[str, list[tuple[str, str]]],
    annotations: dict[str, dict[str, str]],
    output_dir: Path,
    argv: list[str],
    source_label: str = "user-supplied proteomes",
    hits_file: Path | None = None,
    evalue: float = DEFAULT_EVALUE,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
) -> dict[str, object]:
    """Score every pathogen protein and write all outputs.

    Args:
        pathogen: Pathogen proteins as (identifier, sequence).
        host: Host proteins; may be empty when hits_file is given.
        panels: Mapping of panel species name to its proteins.
        annotations: Per-protein annotation fields keyed by identifier.
        output_dir: Directory to write into.
        argv: Command-line arguments, recorded for reproducibility.
        source_label: Description of the inputs for the report header.
        hits_file: Precomputed BLAST/DIAMOND tabular file, or None.
        evalue: E-value threshold.
        min_coverage: Minimum query coverage for a hit to count.

    Returns:
        The result dictionary also written to result.json.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    host_hits, host_searched = resolve_host_hits(pathogen, host, hits_file, evalue, min_coverage)
    panel_counts = count_panel_presence(pathogen, panels, evalue, min_coverage)
    panel_total = len(panels)

    verdicts: list[TargetVerdict] = []
    for identifier, _ in pathogen:
        annotation = annotations.get(identifier, {})
        verdicts.append(
            score_target(
                protein_id=identifier,
                gene_name=annotation.get("gene_name", ""),
                description=annotation.get("description", ""),
                host_hit=host_hits.get(identifier),
                essentiality=annotation.get("essentiality", ""),
                essentiality_source=annotation.get("essentiality_source", ""),
                panel_present=panel_counts.get(identifier, 0),
                panel_total=panel_total,
                host_searched=host_searched,
            )
        )
    verdicts.sort(key=lambda item: (-item.score, item.protein_id))

    write_target_table(tables_dir / "target_scores.csv", verdicts)
    write_host_hits_table(tables_dir / "host_hits.csv", host_hits)

    panel_names = sorted(panels)
    report = build_report(verdicts, source_label, panel_names, host_searched)
    (output_dir / "report.md").write_text(report, encoding="utf-8")

    result: dict[str, object] = {
        "date": date.today().isoformat(),
        "source": source_label,
        "host_searched": host_searched,
        "panel_species": panel_names,
        "parameters": {"evalue": evalue, "min_coverage": min_coverage},
        "targets": [
            {
                "protein_id": verdict.protein_id,
                "gene_name": verdict.gene_name,
                "description": verdict.description,
                "score": verdict.score,
                "decision": verdict.decision,
                "host_identity_pct": verdict.host_identity_pct,
                "panel_present": verdict.panel_present,
                "panel_total": verdict.panel_total,
                "flags": verdict.flags,
                "missing_annotations": verdict.missing_annotations,
                "reasons": verdict.reasons,
            }
            for verdict in verdicts
        ],
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    command = " ".join(["python", "skills/pathogen-target-dissimilarity/pathogen_target_dissimilarity.py", *argv])
    write_commands_sh(output_dir, command)
    write_environment_yml(
        output_dir,
        env_name="clawbio-pathogen-target-dissimilarity",
        pip_deps=["pyhmmer>=0.10"],
        python_version="3.11",
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            tables_dir / "target_scores.csv",
            tables_dir / "host_hits.csv",
        ],
        output_dir,
        anchor=output_dir,
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser.

    Returns:
        Configured parser.
    """
    parser = argparse.ArgumentParser(
        description="Rank pathogen proteins as targets by host dissimilarity, essentiality and conservation."
    )
    parser.add_argument("--pathogen", type=Path, help="Pathogen proteome FASTA")
    parser.add_argument("--host", type=Path, help="Host proteome FASTA")
    parser.add_argument("--panel", type=Path, action="append", default=[], help="Panel proteome FASTA; repeatable")
    parser.add_argument("--annotations", type=Path, help="CSV of gene names, descriptions and essentiality calls")
    parser.add_argument("--hits", type=Path, help="Precomputed BLAST/DIAMOND outfmt 6 hits against the host")
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
    parser.add_argument("--evalue", type=float, default=DEFAULT_EVALUE, help="E-value threshold")
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE, help="Minimum query coverage percent")
    parser.add_argument("--demo", action="store_true", help="Run on built-in synthetic data")
    return parser


def resolve_inputs(args: argparse.Namespace) -> dict[str, object]:
    """Turn parsed arguments into keyword arguments for `run`.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Keyword arguments for `run`.

    Raises:
        SystemExit: If required inputs are missing.
    """
    if args.demo:
        return build_demo_inputs()
    if args.pathogen is None:
        raise SystemExit("A pathogen proteome is required: pass --pathogen, or run with --demo.")
    if args.host is None and args.hits is None:
        raise SystemExit(
            "Host evidence is required: pass --host with a host proteome, or --hits with a "
            "precomputed BLAST/DIAMOND table. Running without either scores every protein "
            "zero on selectivity."
        )
    annotations = annotations_from_fasta(args.pathogen)
    if args.annotations:
        annotations = merge_annotations(annotations, read_annotations(args.annotations))
    panels = {path.stem: read_fasta(path) for path in args.panel}
    labels = [f"{args.pathogen.name}"]
    if args.host is not None:
        labels.append(f"host {args.host.name}")
    if args.hits is not None:
        labels.append(f"hits {args.hits.name}")
    return {
        "pathogen": read_fasta(args.pathogen),
        "host": read_fasta(args.host) if args.host is not None else [],
        "panels": panels,
        "annotations": annotations,
        "source_label": "; ".join(labels),
        "hits_file": args.hits,
    }


def main() -> int:
    """Entry point for the command-line interface.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args()
    inputs = resolve_inputs(args)
    result = run(
        output_dir=args.output,
        argv=sys.argv[1:],
        evalue=args.evalue,
        min_coverage=args.min_coverage,
        **inputs,
    )
    pursued = [target for target in result["targets"] if target["decision"] == "pursue"]
    print(f"Report written to {args.output / 'report.md'}")
    print(f"Scored {len(result['targets'])} proteins; {len(pursued)} at pursue threshold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
