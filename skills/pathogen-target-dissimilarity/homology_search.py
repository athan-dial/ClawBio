"""Homology search between a pathogen proteome and a host or panel proteome.

Two routes to the same record type. `search_pyhmmer` runs phmmer in-process,
which needs no external binary and suits proteome subsets. `parse_tabular_hits`
reads a BLAST or DIAMOND tabular-6 file produced elsewhere, which is what a
real run against a full host proteome uses.

Neither route decides anything. They report identity, coverage and e-value so
the scoring module can band them; a search that was never run is distinguished
from a search that found nothing.
"""

from dataclasses import dataclass
from pathlib import Path

import pyhmmer

DEFAULT_EVALUE = 1e-3
TABULAR_MIN_COLUMNS = 12


@dataclass(frozen=True)
class HomologyHit:
    """One protein-protein match.

    Attributes:
        query_id: Pathogen protein identifier.
        subject_id: Matched host or panel protein identifier.
        identity_pct: Percent identical residues over aligned columns.
        coverage_pct: Percent of the query covered by the alignment, or None
            when the source did not carry query length.
        evalue: Search e-value.
    """

    query_id: str
    subject_id: str
    identity_pct: float
    coverage_pct: float | None
    evalue: float


def read_fasta(path: Path | str) -> list[tuple[str, str]]:
    """Read a protein FASTA file into (identifier, sequence) pairs.

    The identifier is the header up to the first whitespace, kept verbatim so
    that UniProt-style compound identifiers survive intact.

    Args:
        path: Path to the FASTA file.

    Returns:
        List of (identifier, sequence) pairs in file order.

    Raises:
        SystemExit: If the file does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"FASTA file not found: {path}")
    records: list[tuple[str, str]] = []
    identifier = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if identifier:
                records.append((identifier, "".join(chunks)))
            identifier = line[1:].split()[0] if line[1:].split() else ""
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if identifier:
        records.append((identifier, "".join(chunks)))
    return records


def read_fasta_with_headers(path: Path | str) -> list[tuple[str, str, str]]:
    """Read a FASTA file keeping the full header line.

    Args:
        path: Path to the FASTA file.

    Returns:
        List of (identifier, sequence, full header) triples in file order.

    Raises:
        SystemExit: If the file does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"FASTA file not found: {path}")
    records: list[tuple[str, str, str]] = []
    header = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header.split()[0], "".join(chunks), header))
            header = line[1:].strip()
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header:
        records.append((header.split()[0], "".join(chunks), header))
    return records


def parse_uniprot_header(header: str) -> dict[str, str]:
    """Extract gene name and description from a FASTA header.

    Handles the UniProt convention
    ``db|accession|entry_name description OS=... OX=... GN=... PE=... SV=...``
    and leaves anything else alone. When no ``GN=`` field is present the entry
    name is used, so a real-proteome report shows something readable rather
    than a bare accession.

    Args:
        header: Header line without the leading '>'.

    Returns:
        Mapping with 'gene_name' and 'description'. Both may be empty.
    """
    tokens = header.split()
    if not tokens:
        return {"gene_name": "", "description": ""}

    identifier = tokens[0]
    gene_name = ""
    description_parts: list[str] = []
    in_description = True
    for token in tokens[1:]:
        key, separator, value = token.partition("=")
        if separator and len(key) == 2 and key.isupper():
            in_description = False
            if key == "GN":
                gene_name = value
            continue
        if in_description:
            description_parts.append(token)

    fields = identifier.split("|")
    if not gene_name and len(fields) == 3 and fields[0] in {"sp", "tr"}:
        gene_name = fields[2]
    return {"gene_name": gene_name, "description": " ".join(description_parts)}


def parse_tabular_hits(path: Path | str) -> list[HomologyHit]:
    """Read BLAST or DIAMOND tabular output (outfmt 6).

    Expects the standard twelve columns: qseqid, sseqid, pident, length,
    mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore. A
    thirteenth column is read as query length and used for coverage, which is
    what `--outfmt 6 ... qlen` produces.

    Args:
        path: Path to the tabular file.

    Returns:
        List of hits in file order.

    Raises:
        SystemExit: If the file is missing or any data line has too few
            columns. A short line is rejected rather than padded, because
            guessing which column is missing silently corrupts every score
            downstream.
    """
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"Hits file not found: {path}")
    hits: list[HomologyHit] = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < TABULAR_MIN_COLUMNS:
            raise SystemExit(
                f"{path.name} line {number}: expected at least {TABULAR_MIN_COLUMNS} "
                f"tab-separated columns (BLAST/DIAMOND outfmt 6), found {len(fields)}. "
                "Re-export the file rather than letting columns be inferred."
            )
        coverage = None
        if len(fields) >= 13:
            coverage = coverage_from_span(fields[6], fields[7], fields[12])
        hits.append(
            HomologyHit(
                query_id=fields[0],
                subject_id=fields[1],
                identity_pct=float(fields[2]),
                coverage_pct=coverage,
                evalue=float(fields[10]),
            )
        )
    return hits


def coverage_from_span(start: str, end: str, query_length: str) -> float | None:
    """Compute query coverage from alignment span and query length.

    Args:
        start: Alignment start on the query, 1-based.
        end: Alignment end on the query, 1-based.
        query_length: Total query length.

    Returns:
        Percent coverage, or None if any field is not numeric or the length
        is zero.
    """
    try:
        span = abs(int(end) - int(start)) + 1
        total = int(query_length)
    except ValueError:
        return None
    if total <= 0:
        return None
    return round(100.0 * span / total, 1)


def best_hit_per_query(hits: list[HomologyHit]) -> dict[str, HomologyHit]:
    """Select the single strongest hit for each query.

    Ranking is by e-value ascending, with percent identity descending as the
    tie-break, so that equally significant matches resolve to the closer one.

    Args:
        hits: All hits across queries.

    Returns:
        Mapping of query identifier to its best hit.
    """
    best: dict[str, HomologyHit] = {}
    for hit in hits:
        current = best.get(hit.query_id)
        if current is None:
            best[hit.query_id] = hit
            continue
        if (hit.evalue, -hit.identity_pct) < (current.evalue, -current.identity_pct):
            best[hit.query_id] = hit
    return best


def alignment_identity(query_aligned: str, target_aligned: str) -> float:
    """Compute percent identity over aligned columns.

    HMMER lowercases low-posterior residues, so comparison is case-folded.
    Columns where either sequence has a gap are counted in the denominator,
    which is the conservative convention.

    Args:
        query_aligned: Aligned query string.
        target_aligned: Aligned target string.

    Returns:
        Percent identity, or 0.0 for an empty alignment.
    """
    columns = min(len(query_aligned), len(target_aligned))
    if columns == 0:
        return 0.0
    identical = sum(
        1
        for index in range(columns)
        if query_aligned[index].upper() == target_aligned[index].upper()
        and query_aligned[index] not in {"-", "."}
    )
    return round(100.0 * identical / columns, 1)


def decode_name(name: bytes | str) -> str:
    """Return a hit name as text.

    pyhmmer returns sequence names as bytes in some releases and as str in
    others, so both are handled rather than pinning a version.

    Args:
        name: Name as returned by pyhmmer.

    Returns:
        The name as a string.
    """
    return name.decode() if isinstance(name, bytes) else name


def search_pyhmmer(
    queries: list[tuple[str, str]],
    targets: list[tuple[str, str]],
    evalue_threshold: float = DEFAULT_EVALUE,
) -> list[HomologyHit]:
    """Search query proteins against target proteins with phmmer.

    Runs single-threaded so that repeated runs on the same input produce
    identical output.

    Args:
        queries: Pathogen proteins as (identifier, sequence).
        targets: Host or panel proteins as (identifier, sequence).
        evalue_threshold: Hits weaker than this are discarded.

    Returns:
        List of hits, one per query-target domain match that passes the
        threshold. An empty list means no detectable homology, which is a
        result; it is the caller's job not to confuse it with a search that
        never ran.
    """
    if not queries or not targets:
        return []
    alphabet = pyhmmer.easel.Alphabet.amino()
    query_block = [
        pyhmmer.easel.TextSequence(name=identifier.encode(), sequence=sequence).digitize(alphabet)
        for identifier, sequence in queries
    ]
    target_block = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            pyhmmer.easel.TextSequence(name=identifier.encode(), sequence=sequence).digitize(alphabet)
            for identifier, sequence in targets
        ],
    )
    lengths = {identifier: len(sequence) for identifier, sequence in queries}
    hits: list[HomologyHit] = []
    results = pyhmmer.hmmer.phmmer(query_block, target_block, cpus=1, E=evalue_threshold)
    for (query_id, _), top_hits in zip(queries, results):
        for hit in top_hits:
            if hit.evalue > evalue_threshold:
                continue
            domain = hit.best_domain
            alignment = domain.alignment
            identity = alignment_identity(alignment.hmm_sequence, alignment.target_sequence)
            span = abs(alignment.hmm_to - alignment.hmm_from) + 1
            coverage = round(100.0 * span / lengths[query_id], 1) if lengths[query_id] else None
            hits.append(
                HomologyHit(
                    query_id=query_id,
                    subject_id=decode_name(hit.name),
                    identity_pct=identity,
                    coverage_pct=coverage,
                    evalue=float(hit.evalue),
                )
            )
    return hits


def filter_by_coverage(hits: list[HomologyHit], min_coverage: float) -> list[HomologyHit]:
    """Drop hits whose alignment covers too little of the query.

    A high-identity match over a short segment is a shared domain, not an
    orthologue, and treating it as one wrongly disqualifies a target. Hits
    with unknown coverage are kept, since absence of the field is not
    evidence of a short alignment.

    Args:
        hits: Hits to filter.
        min_coverage: Minimum percent of the query that must be aligned.

    Returns:
        The hits that pass.
    """
    return [hit for hit in hits if hit.coverage_pct is None or hit.coverage_pct >= min_coverage]
