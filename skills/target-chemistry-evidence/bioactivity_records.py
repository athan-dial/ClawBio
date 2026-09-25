"""Read and normalise bioactivity records exported from a chemistry database.

The skill never queries ChEMBL, PubChem or BindingDB. It reads a table the
user exported, so the same file gives the same answer on any machine and the
run works with no network at all.

This module does two jobs and nothing else: turn a delimited file into
`BioactivityRecord` objects, and convert an activity value into a comparable
number where the units allow it. It assigns no evidence tier and makes no
judgement about potency; that is `evidence_tiers`.
"""

import csv
from dataclasses import dataclass
from pathlib import Path

# Columns every bioactivity export must carry. A missing one is refused rather
# than defaulted, because a blank target accession is the difference between
# exact-target evidence and organism-level evidence.
REQUIRED_COLUMNS: tuple[str, ...] = (
    "compound_id",
    "target_accession",
    "target_organism",
    "activity_type",
    "activity_value",
    "activity_units",
    "source_id",
)

# Columns that improve attribution but are allowed to be absent entirely.
OPTIONAL_COLUMNS: tuple[str, ...] = (
    "compound_name",
    "target_id",
    "target_pfam",
    "assay_type",
    "assay_organism",
)

# Molar units, converted to nM. Matched exactly, never case-folded: folding
# would map 'mM' and 'M' onto the same key and silently change a value by six
# orders of magnitude.
MOLAR_UNITS_TO_NM: dict[str, float] = {
    "M": 1e9,
    "mM": 1e6,
    "uM": 1e3,
    "µM": 1e3,
    "nM": 1.0,
    "pM": 1e-3,
    "fM": 1e-6,
}

# Mass-concentration units, converted to ug.mL-1. mg.L-1 and ug.mL-1 are the
# same quantity, which is why MIC tables written either way compare directly.
MASS_UNITS_TO_UG_PER_ML: dict[str, float] = {
    "g.L-1": 1000.0,
    "g/L": 1000.0,
    "mg.mL-1": 1000.0,
    "mg/mL": 1000.0,
    "mg.L-1": 1.0,
    "mg/L": 1.0,
    "ug.mL-1": 1.0,
    "ug/mL": 1.0,
    "µg/mL": 1.0,
    "µg.mL-1": 1.0,
    "ng.mL-1": 0.001,
    "ng/mL": 0.001,
}

NANOMOLAR = "nM"
MICROGRAM_PER_ML = "ug.mL-1"

TAB_SUFFIXES: frozenset[str] = frozenset({".tsv", ".tab", ".txt"})


@dataclass(frozen=True)
class BioactivityRecord:
    """One measured compound-target activity, as exported by the user.

    Attributes:
        compound_id: Database identifier of the compound, e.g. a ChEMBL
            molecule identifier.
        compound_name: Human-readable compound name; may be empty.
        target_id: Database identifier of the assay target; may be empty.
        target_accession: UniProt accession of the assayed protein. Empty for
            an organism-level assay, which is the whole point of the
            phenotypic tier.
        target_organism: Organism the assay target belongs to; may be empty.
        target_pfam: Pfam accession of the assayed protein's domain; may be
            empty.
        activity_type: Measurement type, e.g. IC50, Ki, Kd, MIC.
        activity_value: Reported value, or None when it was blank or not a
            number.
        activity_units: Units as reported, verbatim.
        comparable_value: activity_value converted into `comparable_unit`, or
            None when the units are not convertible.
        comparable_unit: 'nM', 'ug.mL-1', or '' when no conversion was made.
        assay_type: Assay classification as exported, e.g. B, F, A.
        assay_organism: Organism the assay was run in; may be empty.
        source_id: Assay or document identifier the measurement came from.
    """

    compound_id: str
    compound_name: str
    target_id: str
    target_accession: str
    target_organism: str
    target_pfam: str
    activity_type: str
    activity_value: float | None
    activity_units: str
    comparable_value: float | None
    comparable_unit: str
    assay_type: str
    assay_organism: str
    source_id: str


def normalise_accession(value: str) -> str:
    """Put a protein accession into a comparable form.

    Upper-cases and drops an isoform or version suffix, so that `P10635-2`
    and `P10635.3` both compare equal to `P10635`. The suffix is dropped
    rather than kept because a bioactivity export and a proteome FASTA rarely
    agree on it, and disagreeing on a suffix would silently demote
    exact-target evidence to no evidence.

    Args:
        value: Accession as written in the input file.

    Returns:
        The normalised accession, or an empty string when the input is blank.
    """
    trimmed = value.strip().upper()
    if not trimmed:
        return ""
    for separator in ("-", "."):
        head, found, _ = trimmed.partition(separator)
        if found:
            trimmed = head
    return trimmed


def accession_from_protein_id(protein_id: str) -> str:
    """Extract a UniProt accession from a proteome identifier.

    Handles the UniProt FASTA convention ``db|accession|entry_name`` and
    leaves anything else alone, so a target list keyed on bare accessions and
    one keyed on full UniProt headers both work.

    Args:
        protein_id: Protein identifier from the target list.

    Returns:
        The accession in normalised form, or the normalised identifier itself
        when it carries no UniProt prefix.
    """
    fields = protein_id.strip().split("|")
    if len(fields) == 3 and fields[0].lower() in {"sp", "tr"}:
        return normalise_accession(fields[1])
    return normalise_accession(protein_id)


def is_uniprot_header_id(protein_id: str) -> bool:
    """Report whether a protein identifier is a UniProt FASTA header identifier.

    Args:
        protein_id: Protein identifier from the target list.

    Returns:
        True for the ``db|accession|entry_name`` form, where the accession is
        stated rather than assumed.
    """
    fields = protein_id.strip().split("|")
    return len(fields) == 3 and fields[0].lower() in {"sp", "tr"} and bool(fields[1].strip())


def organism_matches(left: str, right: str) -> bool:
    """Decide whether two organism labels refer to the same organism.

    Comparison is case-insensitive and tolerates a strain suffix, so that
    `Candida albicans` matches `Candida albicans SC5314`. The suffix must
    start at a space, so `Candida albicans` never matches `Candida
    albicansomething`. Two blank labels never match: absence of a label is
    not evidence of agreement.

    Args:
        left: First organism label.
        right: Second organism label.

    Returns:
        True when the labels refer to the same organism.
    """
    first = left.strip().casefold()
    second = right.strip().casefold()
    if not first or not second:
        return False
    if first == second:
        return True
    return first.startswith(second + " ") or second.startswith(first + " ")


def parse_activity_value(raw: str) -> float | None:
    """Read a reported activity value as a number.

    Args:
        raw: Value as written in the input file.

    Returns:
        The value as a float, or None when it is blank or not a number. A
        qualified value such as `>100` returns None rather than 100, because
        treating a censored measurement as an exact one would report a
        potency that was never observed.
    """
    trimmed = raw.strip()
    if not trimmed:
        return None
    try:
        return float(trimmed)
    except ValueError:
        return None


def convert_activity(value: float | None, units: str) -> tuple[float | None, str]:
    """Convert an activity value into a comparable unit.

    Molar units become nM and mass-concentration units become ug.mL-1. The
    two families are deliberately not interconverted: that needs a molecular
    weight, which a bioactivity export does not carry, and guessing one would
    make a MIC look like a binding constant.

    Args:
        value: Parsed activity value, or None.
        units: Units as reported.

    Returns:
        Tuple of the converted value and its unit, or (None, '') when no
        conversion applies.
    """
    if value is None:
        return None, ""
    trimmed = units.strip()
    molar_factor = MOLAR_UNITS_TO_NM.get(trimmed)
    if molar_factor is not None:
        return value * molar_factor, NANOMOLAR
    mass_factor = MASS_UNITS_TO_UG_PER_ML.get(trimmed)
    if mass_factor is not None:
        return value * mass_factor, MICROGRAM_PER_ML
    return None, ""


def delimiter_for(path: Path) -> str:
    """Choose the field delimiter from a file extension.

    Args:
        path: Path to the data file.

    Returns:
        A tab for `.tsv`, `.tab` and `.txt`, otherwise a comma.
    """
    return "\t" if path.suffix.lower() in TAB_SUFFIXES else ","


def read_delimited(path: Path | str) -> tuple[list[dict[str, str]], list[str]]:
    """Read a CSV or TSV file into rows and its header.

    Args:
        path: Path to the file.

    Returns:
        Tuple of the rows as dictionaries and the header field names.

    Raises:
        SystemExit: If the file does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"File not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter_for(path))
        rows = [row for row in reader]
        fieldnames = list(reader.fieldnames or [])
    return rows, fieldnames


def require_columns(path: Path, fieldnames: list[str], required: tuple[str, ...]) -> None:
    """Refuse a file that is missing a required column.

    Args:
        path: Path the columns came from, used in the message.
        fieldnames: Header fields found in the file.
        required: Columns that must be present.

    Raises:
        SystemExit: If any required column is absent.
    """
    missing = [column for column in required if column not in fieldnames]
    if missing:
        raise SystemExit(
            f"{path.name} is missing required columns: {', '.join(missing)}. "
            f"Found: {', '.join(fieldnames) if fieldnames else 'no header'}. "
            "Re-export the file rather than renaming columns by hand."
        )


def record_from_row(row: dict[str, str]) -> BioactivityRecord:
    """Build one record from a parsed row.

    Args:
        row: Row as read from the bioactivity file.

    Returns:
        The record, with the activity value parsed and converted.
    """
    value = parse_activity_value(row.get("activity_value", "") or "")
    units = (row.get("activity_units", "") or "").strip()
    comparable_value, comparable_unit = convert_activity(value, units)
    return BioactivityRecord(
        compound_id=(row.get("compound_id", "") or "").strip(),
        compound_name=(row.get("compound_name", "") or "").strip(),
        target_id=(row.get("target_id", "") or "").strip(),
        target_accession=normalise_accession(row.get("target_accession", "") or ""),
        target_organism=(row.get("target_organism", "") or "").strip(),
        target_pfam=normalise_accession(row.get("target_pfam", "") or ""),
        activity_type=(row.get("activity_type", "") or "").strip(),
        activity_value=value,
        activity_units=units,
        comparable_value=comparable_value,
        comparable_unit=comparable_unit,
        assay_type=(row.get("assay_type", "") or "").strip(),
        assay_organism=(row.get("assay_organism", "") or "").strip(),
        source_id=(row.get("source_id", "") or "").strip(),
    )


def read_bioactivity(path: Path | str) -> list[BioactivityRecord]:
    """Read a bioactivity export into records.

    Rows with no compound identifier are dropped, because a measurement that
    cannot be traced to a compound cannot be acted on and counting it would
    inflate every tier.

    Args:
        path: Path to the bioactivity CSV or TSV.

    Returns:
        Records in file order.

    Raises:
        SystemExit: If the file is missing or lacks a required column.
    """
    path = Path(path)
    rows, fieldnames = read_delimited(path)
    require_columns(path, fieldnames, REQUIRED_COLUMNS)
    return [record_from_row(row) for row in rows if (row.get("compound_id") or "").strip()]


def record_organism(record: BioactivityRecord) -> str:
    """Return the organism a record's activity was measured in.

    The assay organism wins when present. For an organism-level assay it is
    the only organism field that is filled, and for a protein assay the two
    agree.

    Args:
        record: The bioactivity record.

    Returns:
        The organism label, possibly empty.
    """
    return record.assay_organism or record.target_organism
