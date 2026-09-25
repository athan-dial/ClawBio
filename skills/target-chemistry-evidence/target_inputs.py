"""Read the ranked target list and the optional ortholog mapping.

The target list is the `tables/target_scores.csv` produced by a target-triage
skill, or any CSV carrying a `protein_id` column. Everything else is optional
and every optional field that is absent switches a whole evidence tier from
`0 compounds` to `not assessable`, which is the distinction this skill exists
to preserve.
"""

from dataclasses import dataclass, field
from pathlib import Path

from bioactivity_records import (
    accession_from_protein_id,
    is_uniprot_header_id,
    normalise_accession,
    read_delimited,
    require_columns,
)

# Where a target's accession came from. The fallback matters: an accession
# taken from an internal identifier will match a bioactivity export only if
# that export uses the same internal identifiers, which a database export
# never does.
ACCESSION_FROM_COLUMN = "uniprot_accession column"
ACCESSION_FROM_HEADER = "UniProt FASTA header"
ACCESSION_FROM_PROTEIN_ID = "protein_id fallback"

TARGET_REQUIRED_COLUMNS: tuple[str, ...] = ("protein_id",)
TARGET_OPTIONAL_COLUMNS: tuple[str, ...] = (
    "gene_name",
    "description",
    "score",
    "decision",
    "uniprot_accession",
    "pfam_accession",
    "organism",
)

ORTHOLOG_REQUIRED_COLUMNS: tuple[str, ...] = ("protein_id", "ortholog_accession")
ORTHOLOG_OPTIONAL_COLUMNS: tuple[str, ...] = (
    "ortholog_organism",
    "identity_pct",
    "source",
)


@dataclass(frozen=True)
class Ortholog:
    """One homologous protein in another organism.

    Attributes:
        accession: UniProt accession of the homologue, normalised.
        organism: Organism the homologue belongs to; may be empty.
        identity_pct: Percent identity to the target, or None when the
            mapping file did not carry it.
        source: Where the mapping came from, e.g. a search tool or a
            database release; may be empty.
    """

    accession: str
    organism: str
    identity_pct: float | None
    source: str


@dataclass
class TargetSpec:
    """One target to gather chemistry evidence for.

    Attributes:
        protein_id: Identifier as written in the target list.
        gene_name: Gene symbol; may be empty.
        description: Free-text description; may be empty.
        accession: UniProt accession used for the exact-target tier. Empty
            when none could be derived.
        accession_source: Where that accession came from, which decides how
            much a zero in the exact-target tier is worth.
        organism: Organism used for the phenotypic tier; may be empty.
        pfam_accession: Pfam accession used for the family tier; may be
            empty.
        upstream_score: Score carried over from the target list, verbatim.
        upstream_decision: Decision carried over from the target list,
            verbatim.
        orthologs: Homologues used for the ortholog tier.
    """

    protein_id: str
    gene_name: str
    description: str
    accession: str
    organism: str
    pfam_accession: str
    accession_source: str = ACCESSION_FROM_COLUMN
    upstream_score: str = ""
    upstream_decision: str = ""
    orthologs: list[Ortholog] = field(default_factory=list)


def ortholog_accessions(target: TargetSpec) -> frozenset[str]:
    """Collect the accessions of a target's orthologs.

    Args:
        target: The target.

    Returns:
        The set of normalised ortholog accessions.
    """
    return frozenset(ortholog.accession for ortholog in target.orthologs if ortholog.accession)


def parse_identity(raw: str) -> float | None:
    """Read a percent identity as a number.

    Args:
        raw: Value as written in the mapping file.

    Returns:
        The value as a float, or None when blank or not a number.
    """
    trimmed = raw.strip()
    if not trimmed:
        return None
    try:
        return float(trimmed)
    except ValueError:
        return None


def accession_source_for(protein_id: str, explicit: str) -> str:
    """Say where a target's accession came from.

    Args:
        protein_id: Identifier from the target list.
        explicit: Value of the `uniprot_accession` column, normalised.

    Returns:
        One of the ACCESSION_FROM_* labels.
    """
    if explicit:
        return ACCESSION_FROM_COLUMN
    if is_uniprot_header_id(protein_id):
        return ACCESSION_FROM_HEADER
    return ACCESSION_FROM_PROTEIN_ID


def target_from_row(row: dict[str, str], default_organism: str) -> TargetSpec:
    """Build one target from a row of the target list.

    The accession comes from an explicit `uniprot_accession` column when
    present, then from a UniProt FASTA header identifier, and otherwise from
    `protein_id` itself. That last case is recorded rather than hidden: a
    target list keyed on internal identifiers can never match a database
    export, so an empty exact-target tier there says nothing about the
    chemistry.

    Args:
        row: Row from the target list.
        default_organism: Organism applied when the row carries none, from
            the `--organism` option.

    Returns:
        The target specification.
    """
    protein_id = (row.get("protein_id", "") or "").strip()
    explicit = normalise_accession(row.get("uniprot_accession", "") or "")
    return TargetSpec(
        protein_id=protein_id,
        gene_name=(row.get("gene_name", "") or "").strip(),
        description=(row.get("description", "") or "").strip(),
        accession=explicit or accession_from_protein_id(protein_id),
        organism=(row.get("organism", "") or "").strip() or default_organism,
        pfam_accession=normalise_accession(row.get("pfam_accession", "") or ""),
        accession_source=accession_source_for(protein_id, explicit),
        upstream_score=(row.get("score", "") or "").strip(),
        upstream_decision=(row.get("decision", "") or "").strip(),
    )


def read_targets(path: Path | str, default_organism: str = "") -> list[TargetSpec]:
    """Read the ranked target list.

    Args:
        path: Path to the target CSV, typically `tables/target_scores.csv`
            from an upstream target-triage run.
        default_organism: Organism applied to rows with no `organism` column.

    Returns:
        Targets in file order.

    Raises:
        SystemExit: If the file is missing, lacks `protein_id`, or contains
            no usable rows.
    """
    path = Path(path)
    rows, fieldnames = read_delimited(path)
    require_columns(path, fieldnames, TARGET_REQUIRED_COLUMNS)
    targets = [
        target_from_row(row, default_organism)
        for row in rows
        if (row.get("protein_id") or "").strip()
    ]
    if not targets:
        raise SystemExit(f"{path.name} contains no rows with a protein_id.")
    return targets


def read_orthologs(path: Path | str) -> dict[str, list[Ortholog]]:
    """Read the ortholog mapping.

    Args:
        path: Path to the mapping CSV or TSV.

    Returns:
        Mapping of target protein identifier to its orthologs. A target with
        no row is absent from the mapping, which the caller reads as an
        unassessable ortholog tier rather than as an absence of chemistry.

    Raises:
        SystemExit: If the file is missing or lacks a required column.
    """
    path = Path(path)
    rows, fieldnames = read_delimited(path)
    require_columns(path, fieldnames, ORTHOLOG_REQUIRED_COLUMNS)
    mapping: dict[str, list[Ortholog]] = {}
    for row in rows:
        protein_id = (row.get("protein_id", "") or "").strip()
        accession = normalise_accession(row.get("ortholog_accession", "") or "")
        if not protein_id or not accession:
            continue
        mapping.setdefault(protein_id, []).append(
            Ortholog(
                accession=accession,
                organism=(row.get("ortholog_organism", "") or "").strip(),
                identity_pct=parse_identity(row.get("identity_pct", "") or ""),
                source=(row.get("source", "") or "").strip(),
            )
        )
    return mapping


def attach_orthologs(targets: list[TargetSpec], mapping: dict[str, list[Ortholog]]) -> None:
    """Attach the ortholog mapping to each target in place.

    Args:
        targets: Targets to update.
        mapping: Mapping of protein identifier to orthologs.
    """
    for target in targets:
        target.orthologs = list(mapping.get(target.protein_id, []))
