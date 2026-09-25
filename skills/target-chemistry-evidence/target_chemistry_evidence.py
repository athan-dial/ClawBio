#!/usr/bin/env python3
"""ClawBio skill: find existing chemistry for a ranked target list.

Takes a ranked target list and a bioactivity export, and reports, for each
target, what compounds already have measured activity against it. Evidence is
split into four tiers that are never combined into one number: assays against
the protein itself, against a homologue in another organism, against another
member of the same domain family, and whole-organism activity with no protein
assigned.

Runs offline. The skill never queries ChEMBL, PubChem or BindingDB; it reads
the table the user exported from them. SKILL.md documents the exact queries
that produce that table.
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
from bioactivity_records import BioactivityRecord, read_bioactivity  # noqa: E402
from evidence_tiers import (  # noqa: E402
    EXACT_TARGET,
    NO_EVIDENCE_PROFILE,
    ORTHOLOG,
    PHENOTYPIC,
    TIER_DEFINITIONS,
    TIER_ORDER,
    TargetEvidence,
    compound_count,
    gather_evidence,
)
from target_inputs import (  # noqa: E402
    ACCESSION_FROM_PROTEIN_ID,
    TargetSpec,
    attach_orthologs,
    read_orthologs,
    read_targets,
)

DEMO_DIR = Path(__file__).resolve().parent / "demo"
DEMO_ORGANISM = "Candida albicans"
NOT_ASSESSED = "n/a"

DISCLAIMER = (
    "*ClawBio is a research tool. Not a medical device. Evidence tiers describe "
    "what has been measured, not how good a target is; an empty exact-target tier "
    "usually means nobody has run that assay, not that the protein resists "
    "chemistry.*"
)

DESIGN_NOTE = (
    "Tiers are reported separately and never summed, averaged or weighted into a "
    "single druggability score. A measured IC50 against this protein and a "
    "whole-organism MIC are different claims, and one number covering both would "
    "hide which claim is in hand."
)


def build_demo_inputs() -> dict[str, object]:
    """Load the committed synthetic demo files.

    Returns:
        Keyword arguments for `run`: targets, records and a source label.
    """
    targets = read_targets(DEMO_DIR / "demo_targets.csv", DEMO_ORGANISM)
    attach_orthologs(targets, read_orthologs(DEMO_DIR / "demo_orthologs.csv"))
    return {
        "targets": targets,
        "records": read_bioactivity(DEMO_DIR / "demo_bioactivity.tsv"),
        "source_label": (
            "built-in synthetic demo data (7 targets, 24 fabricated bioactivity "
            "records, no real compound or assay identifiers)"
        ),
        "organism": DEMO_ORGANISM,
    }


def format_count(value: int | None) -> str:
    """Render a per-tier compound count for a table cell.

    Args:
        value: Compound count, or None when the tier was not assessable.

    Returns:
        The count as text, or 'n/a'. A tier that could not be evaluated is
        never shown as 0.
    """
    return NOT_ASSESSED if value is None else str(value)


def format_activity(evidence: TargetEvidence, tier: str) -> str:
    """Render a tier's best activity per unit.

    Args:
        evidence: The target's evidence.
        tier: Tier name.

    Returns:
        A text summary such as '12 nM (IC50)', or an empty string when the
        tier holds no convertible measurement.
    """
    summary = evidence.tiers[tier]
    parts = [
        f"{summary.best[unit].value:g} {unit} ({summary.best[unit].activity_type})"
        for unit in sorted(summary.best)
    ]
    return "; ".join(parts)


def write_summary_table(path: Path, evidences: list[TargetEvidence]) -> None:
    """Write one row per target with the per-tier compound counts.

    Args:
        path: Destination CSV path.
        evidences: Target evidence in output order.
    """
    columns = [
        "protein_id", "gene_name", "description", "accession", "accession_source",
        "organism", "pfam_accession", "upstream_score", "upstream_decision", "profile",
        *[f"{tier}_compounds" for tier in TIER_ORDER],
        "not_assessable",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for evidence in evidences:
            writer.writerow([
                evidence.protein_id,
                evidence.gene_name,
                evidence.description,
                evidence.accession,
                evidence.accession_source,
                evidence.organism,
                evidence.pfam_accession,
                evidence.upstream_score,
                evidence.upstream_decision,
                evidence.profile,
                *[format_count(compound_count(evidence, tier)) for tier in TIER_ORDER],
                "|".join(evidence.not_assessable),
            ])


def write_tier_table(path: Path, evidences: list[TargetEvidence]) -> None:
    """Write one row per target and tier.

    The long shape is deliberate: every row names its tier, so no consumer of
    this table can add two tiers together without noticing.

    Args:
        path: Destination CSV path.
        evidences: Target evidence in output order.
    """
    columns = [
        "protein_id", "gene_name", "tier", "tier_definition", "assessable",
        "unassessable_reason", "compound_count", "record_count",
        "best_value", "best_unit", "best_activity_type", "best_compound_id",
        "activity_types", "uncomparable_records", "source_ids",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for evidence in evidences:
            for tier in TIER_ORDER:
                summary = evidence.tiers[tier]
                units = sorted(summary.best)
                writer.writerow([
                    evidence.protein_id,
                    evidence.gene_name,
                    tier,
                    TIER_DEFINITIONS[tier],
                    "yes" if summary.assessable else "no",
                    summary.unassessable_reason,
                    format_count(compound_count(evidence, tier)),
                    summary.record_count if summary.assessable else NOT_ASSESSED,
                    "|".join(f"{summary.best[unit].value:g}" for unit in units),
                    "|".join(units),
                    "|".join(summary.best[unit].activity_type for unit in units),
                    "|".join(summary.best[unit].compound_id for unit in units),
                    "|".join(summary.activity_types),
                    summary.uncomparable_records,
                    "|".join(summary.source_ids),
                ])


def write_compound_table(
    path: Path, attributions: list[tuple[str, BioactivityRecord, str]]
) -> None:
    """Write one row per attributed bioactivity measurement.

    Args:
        path: Destination CSV path.
        attributions: (protein identifier, record, tier) triples.
    """
    columns = [
        "protein_id", "tier", "compound_id", "compound_name", "activity_type",
        "activity_value", "activity_units", "comparable_value", "comparable_unit",
        "assay_type", "assay_organism", "target_accession", "target_pfam", "source_id",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for protein_id, record, tier in attributions:
            writer.writerow([
                protein_id,
                tier,
                record.compound_id,
                record.compound_name,
                record.activity_type,
                "" if record.activity_value is None else f"{record.activity_value:g}",
                record.activity_units,
                "" if record.comparable_value is None else f"{record.comparable_value:g}",
                record.comparable_unit,
                record.assay_type,
                record.assay_organism,
                record.target_accession,
                record.target_pfam,
                record.source_id,
            ])


def count_with_evidence(evidences: list[TargetEvidence], tier: str) -> int:
    """Count targets carrying at least one compound in a tier.

    Args:
        evidences: Target evidence.
        tier: Tier name.

    Returns:
        The number of targets with evidence in that tier.
    """
    return sum(1 for evidence in evidences if evidence.tiers[tier].compound_ids)


def count_unassessable(evidences: list[TargetEvidence], tier: str) -> int:
    """Count targets whose tier could not be evaluated.

    Args:
        evidences: Target evidence.
        tier: Tier name.

    Returns:
        The number of targets for which the tier was not assessable.
    """
    return sum(1 for evidence in evidences if not evidence.tiers[tier].assessable)


def first_unassessable_reason(evidences: list[TargetEvidence], tier: str) -> str:
    """Return the reason the first target could not be assessed in a tier.

    Args:
        evidences: Target evidence.
        tier: Tier name.

    Returns:
        The reason text, or an empty string when every target was assessable.
    """
    for evidence in evidences:
        if not evidence.tiers[tier].assessable:
            return evidence.tiers[tier].unassessable_reason
    return ""


def build_summary_lines(evidences: list[TargetEvidence]) -> list[str]:
    """Build the per-target summary table.

    Args:
        evidences: Target evidence in output order.

    Returns:
        Markdown lines for the table.
    """
    lines = [
        "| Target | Description | Exact | Ortholog | Family | Phenotypic | Profile |",
        "|---|---|---|---|---|---|---|",
    ]
    for evidence in evidences:
        counts = [format_count(compound_count(evidence, tier)) for tier in TIER_ORDER]
        lines.append(
            f"| {evidence.gene_name or evidence.protein_id} | {evidence.description} | "
            f"{' | '.join(counts)} | {evidence.profile} |"
        )
    return lines


def build_activity_lines(evidences: list[TargetEvidence]) -> list[str]:
    """Build the best-activity-per-tier table.

    Args:
        evidences: Target evidence in output order.

    Returns:
        Markdown lines, or a single explanatory line when nothing converted.
    """
    rows: list[str] = []
    for evidence in evidences:
        for tier in TIER_ORDER:
            activity = format_activity(evidence, tier)
            if activity:
                rows.append(
                    f"| {evidence.gene_name or evidence.protein_id} | {tier} | {activity} |"
                )
    if not rows:
        return ["No measurement in any tier carried convertible units."]
    return ["| Target | Tier | Best activity |", "|---|---|---|", *rows]


def build_interpretation(
    evidences: list[TargetEvidence], unattributed: int, total_records: int
) -> list[str]:
    """Build the interpretation section.

    Args:
        evidences: Target evidence in output order.
        unattributed: Records that matched no target in any tier.
        total_records: Records read from the bioactivity file.

    Returns:
        Markdown lines.
    """
    exact = count_with_evidence(evidences, EXACT_TARGET)
    phenotypic_only = [
        evidence for evidence in evidences if evidence.profile == PHENOTYPIC + "_only"
    ]
    ortholog_carrying = count_with_evidence(evidences, ORTHOLOG)
    empty = [evidence for evidence in evidences if evidence.profile == NO_EVIDENCE_PROFILE]

    lines = ["## Interpretation", "", DESIGN_NOTE, ""]
    lines.append(
        f"**Exact-target evidence**: {exact} of {len(evidences)} targets. "
        "Target-level annotation for fungal proteins is sparse in public bioactivity "
        "databases, so an empty exact-target tier is the common case and is not a "
        "failure of the search; it is a statement about what has been deposited."
    )
    lines.append("")
    fallback = [
        evidence
        for evidence in evidences
        if evidence.accession_source == ACCESSION_FROM_PROTEIN_ID
        and not evidence.tiers[EXACT_TARGET].compound_ids
    ]
    if fallback:
        lines.append(
            f"**Accessions taken from `protein_id`**: {len(fallback)} of {len(evidences)} "
            "targets have an empty exact-target tier and an accession that was derived "
            "from their identifier rather than supplied. An internal identifier cannot "
            "match a database export, so those zeros carry no information at all. Add a "
            "`uniprot_accession` column before reading them as an absence of chemistry."
        )
        lines.append("")
    if phenotypic_only:
        names = ", ".join(
            evidence.gene_name or evidence.protein_id for evidence in phenotypic_only
        )
        lines.append(
            f"**Phenotypic evidence and nothing else**: {len(phenotypic_only)} of "
            f"{len(evidences)} targets ({names}). Whole-organism activity says a compound "
            "kills the organism. It says nothing about which protein it hits, so it is "
            "not evidence that this target is engaged and must not be reported as "
            "chemistry against it."
        )
        lines.append("")
    if ortholog_carrying:
        lines.append(
            f"**Ortholog evidence**: {ortholog_carrying} of {len(evidences)} targets. "
            "Read it with the tension it carries: a target with rich ortholog chemistry "
            "is by definition similar to a well-studied protein in another organism, and "
            "that same similarity is what creates selectivity risk. Ortholog chemistry is "
            "a starting point for a chemical series and a warning about the host "
            "counterpart at the same time."
        )
        lines.append("")
    if empty:
        lines.append(
            f"**No evidence in any assessable tier**: {len(empty)} of {len(evidences)} "
            "targets. That is a statement about the export in hand, not about the "
            "protein; widen the export or add an ortholog mapping before concluding the "
            "chemistry does not exist."
        )
        lines.append("")
    unassessed_any = False
    for tier in TIER_ORDER:
        unassessed = count_unassessable(evidences, tier)
        if not unassessed:
            continue
        unassessed_any = True
        lines.append(
            f"- `{tier}` could not be evaluated for {unassessed} of {len(evidences)} "
            f"targets: {first_unassessable_reason(evidences, tier)}"
        )
    if unassessed_any:
        lines.append("")
    lines.append(
        f"**Unattributed records**: {unattributed} of {total_records} records in the "
        "export matched no target in any tier. A high proportion usually means the "
        "export was pulled for a wider organism or target set than this list covers."
    )
    lines.append("")
    return lines


def build_report(
    evidences: list[TargetEvidence],
    source_label: str,
    organism: str,
    unattributed: int,
    total_records: int,
) -> str:
    """Assemble the markdown report.

    Args:
        evidences: Target evidence in output order.
        source_label: Description of the inputs.
        organism: Organism used for the phenotypic tier.
        unattributed: Records that matched no target.
        total_records: Records read from the bioactivity file.

    Returns:
        Complete report text.
    """
    with_any = sum(1 for evidence in evidences if evidence.profile != NO_EVIDENCE_PROFILE)
    lines = [
        "# Target Chemistry Evidence Report",
        "",
        f"**Input**: {source_label}  ",
        f"**Organism**: {organism or 'not supplied'}  ",
        f"**Date**: {date.today().isoformat()}",
        "",
        f"Checked **{len(evidences)} targets** against **{total_records} bioactivity "
        f"records**; **{with_any}** carry evidence in at least one tier.",
        "",
        "## Evidence Tiers",
        "",
        "| Tier | Meaning |",
        "|---|---|",
    ]
    lines.extend(f"| `{tier}` | {TIER_DEFINITIONS[tier]} |" for tier in TIER_ORDER)
    lines.extend([
        "",
        "Counts below are distinct compounds per tier. `n/a` means the tier could not "
        "be evaluated for that target, which is not the same as no compounds found.",
        "",
    ])
    lines.extend(build_summary_lines(evidences))
    lines.extend(["", "## Best Activity Per Tier", ""])
    lines.extend(build_activity_lines(evidences))
    lines.extend(["", ""])
    lines.extend(build_interpretation(evidences, unattributed, total_records))
    lines.extend([DISCLAIMER, ""])
    return "\n".join(lines)


def tier_payload(evidence: TargetEvidence, tier: str) -> dict[str, object]:
    """Build the machine-readable payload for one tier.

    Args:
        evidence: The target's evidence.
        tier: Tier name.

    Returns:
        A dictionary for result.json.
    """
    summary = evidence.tiers[tier]
    return {
        "definition": TIER_DEFINITIONS[tier],
        "assessable": summary.assessable,
        "unassessable_reason": summary.unassessable_reason,
        "compound_count": compound_count(evidence, tier),
        "compound_ids": summary.compound_ids,
        "record_count": summary.record_count if summary.assessable else None,
        "activity_types": summary.activity_types,
        "uncomparable_records": summary.uncomparable_records,
        "source_ids": summary.source_ids,
        "best": {
            unit: {
                "value": summary.best[unit].value,
                "unit": summary.best[unit].unit,
                "activity_type": summary.best[unit].activity_type,
                "compound_id": summary.best[unit].compound_id,
                "source_id": summary.best[unit].source_id,
            }
            for unit in sorted(summary.best)
        },
    }


def build_result(
    evidences: list[TargetEvidence],
    source_label: str,
    organism: str,
    unattributed: int,
    total_records: int,
) -> dict[str, object]:
    """Build the machine-readable result.

    Args:
        evidences: Target evidence in output order.
        source_label: Description of the inputs.
        organism: Organism used for the phenotypic tier.
        unattributed: Records that matched no target.
        total_records: Records read from the bioactivity file.

    Returns:
        The result dictionary written to result.json.
    """
    return {
        "date": date.today().isoformat(),
        "source": source_label,
        "organism": organism,
        "records_read": total_records,
        "records_unattributed": unattributed,
        "tier_definitions": dict(TIER_DEFINITIONS),
        "design": {"tiers_combined": False, "note": DESIGN_NOTE},
        "targets": [
            {
                "protein_id": evidence.protein_id,
                "gene_name": evidence.gene_name,
                "description": evidence.description,
                "accession": evidence.accession,
                "accession_source": evidence.accession_source,
                "organism": evidence.organism,
                "pfam_accession": evidence.pfam_accession,
                "upstream_score": evidence.upstream_score,
                "upstream_decision": evidence.upstream_decision,
                "profile": evidence.profile,
                "not_assessable": evidence.not_assessable,
                "tiers": {tier: tier_payload(evidence, tier) for tier in TIER_ORDER},
            }
            for evidence in evidences
        ],
    }


def run(
    targets: list[TargetSpec],
    records: list[BioactivityRecord],
    output_dir: Path,
    argv: list[str],
    source_label: str = "user-supplied target list and bioactivity export",
    organism: str = "",
) -> dict[str, object]:
    """Gather evidence for every target and write all outputs.

    Target order is preserved from the input file, because that file is
    already a ranking and re-ordering it here would silently replace the
    upstream priority with a chemistry-driven one.

    Args:
        targets: Targets to gather evidence for.
        records: Bioactivity records read from the export.
        output_dir: Directory to write into.
        argv: Command-line arguments, recorded for reproducibility.
        source_label: Description of the inputs for the report header.
        organism: Organism used for the phenotypic tier.

    Returns:
        The result dictionary also written to result.json.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    evidences: list[TargetEvidence] = []
    attributions: list[tuple[str, BioactivityRecord, str]] = []
    matched: set[int] = set()
    for target in targets:
        evidence, attributed = gather_evidence(target, records)
        evidences.append(evidence)
        for index, tier in attributed:
            attributions.append((target.protein_id, records[index], tier))
            matched.add(index)
    unattributed = len(records) - len(matched)

    write_summary_table(tables_dir / "target_evidence_summary.csv", evidences)
    write_tier_table(tables_dir / "evidence_by_tier.csv", evidences)
    write_compound_table(tables_dir / "compound_evidence.csv", attributions)

    report = build_report(evidences, source_label, organism, unattributed, len(records))
    (output_dir / "report.md").write_text(report, encoding="utf-8")

    result = build_result(evidences, source_label, organism, unattributed, len(records))
    (output_dir / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    command = " ".join(
        ["python", "skills/target-chemistry-evidence/target_chemistry_evidence.py", *argv]
    )
    write_commands_sh(output_dir, command)
    write_environment_yml(
        output_dir,
        env_name="clawbio-target-chemistry-evidence",
        pip_deps=[],
        python_version="3.11",
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            tables_dir / "target_evidence_summary.csv",
            tables_dir / "evidence_by_tier.csv",
            tables_dir / "compound_evidence.csv",
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
        description=(
            "Find existing compounds with measured activity against a ranked target "
            "list, keeping exact-target, ortholog, family and phenotypic evidence "
            "strictly separate."
        )
    )
    parser.add_argument("--targets", type=Path, help="Ranked target list CSV")
    parser.add_argument("--bioactivity", type=Path, help="Bioactivity export CSV or TSV")
    parser.add_argument("--orthologs", type=Path, help="Ortholog mapping CSV or TSV")
    parser.add_argument(
        "--organism",
        default="",
        help="Organism for the phenotypic tier when the target list carries no organism column",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
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
    if args.targets is None:
        raise SystemExit("A target list is required: pass --targets, or run with --demo.")
    if args.bioactivity is None:
        raise SystemExit(
            "A bioactivity export is required: pass --bioactivity. This skill has no "
            "network access to ChEMBL, PubChem or BindingDB; SKILL.md documents the "
            "queries that produce the file."
        )
    targets = read_targets(args.targets, args.organism)
    if args.orthologs is not None:
        attach_orthologs(targets, read_orthologs(args.orthologs))
    labels = [args.targets.name, f"bioactivity {args.bioactivity.name}"]
    if args.orthologs is not None:
        labels.append(f"orthologs {args.orthologs.name}")
    return {
        "targets": targets,
        "records": read_bioactivity(args.bioactivity),
        "source_label": "; ".join(labels),
        "organism": args.organism,
    }


def main() -> int:
    """Entry point for the command-line interface.

    Returns:
        Process exit code.
    """
    args = build_parser().parse_args()
    inputs = resolve_inputs(args)
    result = run(output_dir=args.output, argv=sys.argv[1:], **inputs)
    targets = result["targets"]
    with_exact = sum(
        1 for target in targets if target["tiers"][EXACT_TARGET]["compound_ids"]
    )
    print(f"Report written to {args.output / 'report.md'}")
    print(f"Checked {len(targets)} targets; {with_exact} carry exact-target evidence")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
