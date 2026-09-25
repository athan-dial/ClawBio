#!/usr/bin/env python3
"""Regenerate the synthetic demo files for target-chemistry-evidence.

Everything written here is fabricated. Gene names and Pfam accessions are
real so the report reads like a real one, but every compound identifier,
assay identifier and activity value is invented, and the accessions carry a
`DEMO` prefix so they can never be mistaken for UniProt entries. Nothing in
`demo/` is evidence about any compound or any protein.

Run this only to regenerate the committed files in demo/; the skill itself
reads those files and never calls this module.
"""

import csv
from pathlib import Path

DEMO_DIR = Path(__file__).resolve().parent / "demo"
ORGANISM = "Candida albicans"
HOST_ORGANISM = "Homo sapiens"

# protein_id, gene, description, demo accession, Pfam, organism, score, decision.
# A blank organism falls back to the --organism option, which is how a target
# list exported straight from a target-triage run behaves.
TARGETS: tuple[tuple[str, str, str, str, str, str, str, str], ...] = (
    ("CAL_ERG11", "ERG11", "lanosterol 14-alpha-demethylase", "DEMO0001", "PF00067", "", "80", "pursue"),
    ("CAL_FKS1", "FKS1", "beta-1,3-glucan synthase catalytic subunit", "DEMO0002", "PF02364", "", "100", "pursue"),
    ("CAL_HSP90", "HSP90", "molecular chaperone", "DEMO0003", "PF02518", "", "55", "hold"),
    ("CAL_CHS1", "CHS1", "chitin synthase", "DEMO0004", "PF01644", "", "72", "pursue"),
    ("CAL_GWT1", "GWT1", "GPI inositol acyltransferase", "DEMO0005", "PF06423", "", "75", "pursue"),
    ("CAL_PYRE", "PYRE", "dihydroorotate dehydrogenase", "DEMO0006", "PF01180", "", "47", "hold"),
    ("CGL_ERG3", "ERG3", "C-5 sterol desaturase", "DEMO0007", "PF04116", "Candida glabrata", "68", "hold"),
)

# protein_id, ortholog accession, organism, identity, source
ORTHOLOGS: tuple[tuple[str, str, str, str, str], ...] = (
    ("CAL_HSP90", "DEMOH003", HOST_ORGANISM, "58.0", "synthetic demo mapping"),
    ("CAL_PYRE", "DEMOH006", HOST_ORGANISM, "31.0", "synthetic demo mapping"),
    ("CAL_ERG11", "DEMOX001", "Aspergillus fumigatus", "62.0", "synthetic demo mapping"),
)

BIOACTIVITY_COLUMNS: tuple[str, ...] = (
    "compound_id",
    "compound_name",
    "target_id",
    "target_accession",
    "target_organism",
    "target_pfam",
    "activity_type",
    "activity_value",
    "activity_units",
    "assay_type",
    "assay_organism",
    "source_id",
)

# Each row is fabricated. The set is chosen to exercise every tier, both unit
# families, an uncensored-but-unconvertible measurement, a blank value, and a
# record belonging to an organism outside this target list.
BIOACTIVITY: tuple[tuple[str, ...], ...] = (
    # ERG11: assayed against the protein itself, plus whole-organism activity.
    ("DEMO-CPD-001", "demo-azole-A", "DEMOTGT01", "DEMO0001", ORGANISM, "PF00067", "IC50", "42", "nM", "B", ORGANISM, "DEMO-ASSAY-001"),
    ("DEMO-CPD-002", "demo-azole-B", "DEMOTGT01", "DEMO0001", ORGANISM, "PF00067", "IC50", "310", "nM", "B", ORGANISM, "DEMO-ASSAY-002"),
    ("DEMO-CPD-003", "demo-azole-C", "DEMOTGT01", "DEMO0001", ORGANISM, "PF00067", "Ki", "0.9", "uM", "B", ORGANISM, "DEMO-ASSAY-003"),
    ("DEMO-CPD-001", "demo-azole-A", "", "", ORGANISM, "", "MIC", "0.25", "ug.mL-1", "F", "Candida albicans SC5314", "DEMO-ASSAY-004"),
    # ERG11 family member in another organism, and a same-family protein.
    ("DEMO-CPD-004", "demo-azole-D", "DEMOTGT02", "DEMOX001", "Aspergillus fumigatus", "PF00067", "IC50", "88", "nM", "B", "Aspergillus fumigatus", "DEMO-ASSAY-005"),
    ("DEMO-CPD-005", "demo-P450-probe", "DEMOTGT03", "DEMOF001", "Trypanosoma cruzi", "PF00067", "IC50", "1200", "nM", "B", "Trypanosoma cruzi", "DEMO-ASSAY-006"),
    # FKS1: whole-organism only. This is the case the report must make visible.
    ("DEMO-CPD-010", "demo-echinocandin-A", "", "", ORGANISM, "", "MIC", "0.06", "ug.mL-1", "F", ORGANISM, "DEMO-ASSAY-010"),
    ("DEMO-CPD-011", "demo-echinocandin-B", "", "", ORGANISM, "", "MIC", "0.12", "ug.mL-1", "F", ORGANISM, "DEMO-ASSAY-011"),
    ("DEMO-CPD-012", "demo-echinocandin-C", "", "", ORGANISM, "", "MIC", "0.5", "mg.L-1", "F", ORGANISM, "DEMO-ASSAY-012"),
    ("DEMO-CPD-013", "demo-extract-1", "", "", ORGANISM, "", "Inhibition", "62", "%", "F", ORGANISM, "DEMO-ASSAY-013"),
    ("DEMO-CPD-014", "demo-extract-2", "", "", ORGANISM, "", "MIC", "", "ug.mL-1", "F", ORGANISM, "DEMO-ASSAY-014"),
    # HSP90: no fungal protein assay, but a human ortholog with rich chemistry.
    ("DEMO-CPD-020", "demo-hsp90-inhibitor-A", "DEMOTGT10", "DEMOH003", HOST_ORGANISM, "PF02518", "Kd", "7", "nM", "B", HOST_ORGANISM, "DEMO-ASSAY-020"),
    ("DEMO-CPD-021", "demo-hsp90-inhibitor-B", "DEMOTGT10", "DEMOH003", HOST_ORGANISM, "PF02518", "IC50", "35", "nM", "B", HOST_ORGANISM, "DEMO-ASSAY-021"),
    ("DEMO-CPD-022", "demo-hsp90-inhibitor-C", "DEMOTGT11", "DEMOF003", "Mus musculus", "PF02518", "IC50", "140", "nM", "B", "Mus musculus", "DEMO-ASSAY-022"),
    # CHS1: nothing against the protein or a homologue, only same-family chemistry.
    ("DEMO-CPD-030", "demo-chitin-probe-A", "DEMOTGT20", "DEMOF004", "Spodoptera frugiperda", "PF01644", "IC50", "4.5", "uM", "B", "Spodoptera frugiperda", "DEMO-ASSAY-030"),
    ("DEMO-CPD-031", "demo-chitin-probe-B", "DEMOTGT20", "DEMOF004", "Spodoptera frugiperda", "PF01644", "IC50", "", "", "B", "Spodoptera frugiperda", "DEMO-ASSAY-031"),
    # PYRE: a human ortholog and whole-organism activity, but no fungal protein assay.
    ("DEMO-CPD-040", "demo-dhodh-inhibitor-A", "DEMOTGT30", "DEMOH006", HOST_ORGANISM, "PF01180", "IC50", "18", "nM", "B", HOST_ORGANISM, "DEMO-ASSAY-040"),
    ("DEMO-CPD-041", "demo-dhodh-inhibitor-B", "DEMOTGT30", "DEMOH006", HOST_ORGANISM, "PF01180", "IC50", "0.25", "uM", "B", HOST_ORGANISM, "DEMO-ASSAY-041"),
    ("DEMO-CPD-042", "demo-dhodh-inhibitor-C", "", "", ORGANISM, "", "MIC", "8", "ug.mL-1", "F", ORGANISM, "DEMO-ASSAY-042"),
    # Records that belong to no target in this list.
    ("DEMO-CPD-050", "demo-antibacterial-A", "", "", "Escherichia coli", "", "MIC", "2", "ug.mL-1", "F", "Escherichia coli", "DEMO-ASSAY-050"),
    ("DEMO-CPD-051", "demo-kinase-inhibitor", "DEMOTGT40", "DEMOF050", HOST_ORGANISM, "PF00069", "IC50", "3", "nM", "B", HOST_ORGANISM, "DEMO-ASSAY-051"),
    ("DEMO-CPD-052", "demo-unrelated-A", "DEMOTGT41", "DEMOF051", HOST_ORGANISM, "PF00089", "Ki", "64", "nM", "B", HOST_ORGANISM, "DEMO-ASSAY-052"),
    ("DEMO-CPD-053", "demo-unrelated-B", "", "", "Aspergillus fumigatus", "", "MIC", "1", "ug.mL-1", "F", "Aspergillus fumigatus", "DEMO-ASSAY-053"),
    ("DEMO-CPD-054", "demo-unrelated-C", "DEMOTGT42", "DEMOF052", "Rattus norvegicus", "PF00001", "EC50", "500", "nM", "F", "Rattus norvegicus", "DEMO-ASSAY-054"),
)


def write_targets(path: Path) -> int:
    """Write the demo target list.

    Args:
        path: Destination CSV path.

    Returns:
        Number of targets written.
    """
    columns = [
        "protein_id", "gene_name", "description", "score", "decision",
        "uniprot_accession", "pfam_accession", "organism",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for protein_id, gene, description, accession, pfam, organism, score, decision in TARGETS:
            writer.writerow(
                [protein_id, gene, description, score, decision, accession, pfam, organism]
            )
    return len(TARGETS)


def write_orthologs(path: Path) -> int:
    """Write the demo ortholog mapping.

    Args:
        path: Destination CSV path.

    Returns:
        Number of mappings written.
    """
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["protein_id", "ortholog_accession", "ortholog_organism", "identity_pct", "source"]
        )
        writer.writerows(ORTHOLOGS)
    return len(ORTHOLOGS)


def write_bioactivity(path: Path) -> int:
    """Write the demo bioactivity export.

    Args:
        path: Destination TSV path.

    Returns:
        Number of records written.
    """
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(BIOACTIVITY_COLUMNS)
        writer.writerows(BIOACTIVITY)
    return len(BIOACTIVITY)


def main() -> int:
    """Generate every demo file.

    Returns:
        Process exit code.
    """
    DEMO_DIR.mkdir(parents=True, exist_ok=True)
    targets = write_targets(DEMO_DIR / "demo_targets.csv")
    orthologs = write_orthologs(DEMO_DIR / "demo_orthologs.csv")
    records = write_bioactivity(DEMO_DIR / "demo_bioactivity.tsv")
    print(f"Wrote {targets} targets, {orthologs} ortholog mappings, {records} bioactivity records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
