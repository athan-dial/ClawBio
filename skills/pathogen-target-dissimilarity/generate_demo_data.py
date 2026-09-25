#!/usr/bin/env python3
"""Regenerate the synthetic demo proteomes for pathogen-target-dissimilarity.

The demo data is synthetic. Gene names are real so the report reads like a
real one, but the sequences are generated, not the actual proteins, and the
essentiality calls are placeholders. Nothing here is evidence about any gene.

Run this only to regenerate the committed files in demo/; the skill itself
reads those files and never calls this module.
"""

import csv
import random
from pathlib import Path

DEMO_DIR = Path(__file__).resolve().parent / "demo"
SEED = 20260925

# Amino acid frequencies roughly matching observed proteome composition.
RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
WEIGHTS = [8.3, 1.4, 5.4, 6.2, 3.9, 7.1, 2.3, 5.9, 5.8, 9.7, 2.4, 4.1, 4.7, 3.9, 5.5, 6.6, 5.4, 6.9, 1.1, 2.9]

# Substitution groups keep a mutated copy alignable while lowering identity.
SIMILAR = {
    "A": "GSTV", "C": "SA", "D": "ENQ", "E": "DQK", "F": "YWL", "G": "ASN",
    "H": "NQY", "I": "LVMF", "K": "REQ", "L": "IVMF", "M": "LIV", "N": "DQSH",
    "P": "AGS", "Q": "ENKR", "R": "KQH", "S": "TAGN", "T": "SAVN", "V": "ILAM",
    "W": "FY", "Y": "FWH",
}

# gene, description, length, host identity target (None = no host counterpart),
# essentiality call, essentiality source, panel species carrying an orthologue
SPEC: list[tuple[str, str, int, float | None, str, str, int]] = [
    ("FKS1", "beta-1,3-glucan synthase catalytic subunit", 300, None, "essential", "GRACE (synthetic)", 3),
    ("GWT1", "GPI inositol acyltransferase", 260, None, "essential", "GRACE (synthetic)", 2),
    ("CHS1", "chitin synthase", 280, None, "unknown", "", 3),
    ("ERG6", "sterol 24-C-methyltransferase", 240, None, "unknown", "", 2),
    ("ERG11", "lanosterol 14-alpha-demethylase", 290, 0.36, "essential", "GRACE (synthetic)", 3),
    ("PYRE", "dihydroorotate dehydrogenase", 250, 0.31, "unknown", "", 1),
    ("HSP90", "molecular chaperone", 300, 0.55, "essential", "GRACE (synthetic)", 3),
    ("RPL3", "60S ribosomal protein L3", 270, 0.66, "essential", "GRACE (synthetic)", 3),
    ("TUB2", "beta-tubulin", 280, 0.76, "essential", "GRACE (synthetic)", 3),
    ("ACT1", "actin", 290, 0.90, "essential", "GRACE (synthetic)", 3),
]

PANEL_SPECIES = ("Candida_auris", "Candida_glabrata", "Candida_tropicalis")
PANEL_IDENTITY = 0.86
HOST_DECOY_COUNT = 6


def random_protein(rng: random.Random, length: int) -> str:
    """Generate a random protein sequence with plausible composition.

    Args:
        rng: Seeded random generator.
        length: Number of residues.

    Returns:
        A protein sequence string starting with methionine.
    """
    body = rng.choices(RESIDUES, weights=WEIGHTS, k=length - 1)
    return "M" + "".join(body)


def mutate(rng: random.Random, sequence: str, identity: float) -> str:
    """Produce a homologue of a sequence at approximately a target identity.

    Substitutions are drawn from a similarity group so the pair stays
    alignable; identity falls because the residue still differs.

    Args:
        rng: Seeded random generator.
        sequence: Source sequence.
        identity: Fraction of positions to leave unchanged, 0 to 1.

    Returns:
        The mutated sequence, same length as the source.
    """
    residues: list[str] = []
    for residue in sequence:
        if rng.random() < identity:
            residues.append(residue)
        else:
            residues.append(rng.choice(SIMILAR.get(residue, RESIDUES)))
    return "".join(residues)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write records to a FASTA file, wrapped at 60 columns.

    Args:
        path: Destination path.
        records: (identifier with description, sequence) pairs.
    """
    lines: list[str] = []
    for header, sequence in records:
        lines.append(f">{header}")
        lines.extend(sequence[index : index + 60] for index in range(0, len(sequence), 60))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    """Generate every demo file.

    Returns:
        Process exit code.
    """
    rng = random.Random(SEED)
    DEMO_DIR.mkdir(parents=True, exist_ok=True)

    pathogen: list[tuple[str, str]] = []
    host: list[tuple[str, str]] = []
    panels: dict[str, list[tuple[str, str]]] = {species: [] for species in PANEL_SPECIES}
    annotations: list[list[str]] = [
        ["protein_id", "gene_name", "description", "essentiality", "essentiality_source"]
    ]

    for gene, description, length, host_identity, call, source, panel_count in SPEC:
        sequence = random_protein(rng, length)
        protein_id = f"CAL_{gene}"
        pathogen.append((f"{protein_id} {description}", sequence))
        annotations.append([protein_id, gene, description, call, source])
        if host_identity is not None:
            host.append((f"HSA_{gene}_like human homologue of {gene}", mutate(rng, sequence, host_identity)))
        for species in PANEL_SPECIES[:panel_count]:
            prefix = "".join(part[0] for part in species.split("_")).upper()
            panels[species].append((f"{prefix}_{gene}", mutate(rng, sequence, PANEL_IDENTITY)))

    for index in range(HOST_DECOY_COUNT):
        host.append((f"HSA_DECOY_{index + 1} unrelated human protein", random_protein(rng, 250)))

    write_fasta(DEMO_DIR / "demo_pathogen.faa", pathogen)
    write_fasta(DEMO_DIR / "demo_host.faa", host)
    for species, records in panels.items():
        write_fasta(DEMO_DIR / f"demo_panel_{species.lower()}.faa", records)
    with (DEMO_DIR / "demo_essentiality.csv").open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows(annotations)

    print(f"Wrote {len(pathogen)} pathogen proteins, {len(host)} host proteins")
    for species, records in panels.items():
        print(f"  panel {species}: {len(records)} proteins")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
