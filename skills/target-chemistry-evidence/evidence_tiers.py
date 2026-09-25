"""Assign bioactivity records to evidence tiers and summarise each tier alone.

Four tiers, reported separately and never combined:

1. `exact_target`  bioactivity assayed against that protein in that organism
2. `ortholog`      bioactivity against a homologue in another organism
3. `family`        bioactivity against the same domain family
4. `phenotypic`    whole-organism activity with no protein assigned

There is deliberately no total, no weighted sum and no druggability score. A
compound with a measured IC50 against the actual protein and a compound with
a whole-organism MIC are different kinds of claim, and one number covering
both would hide which kind was in hand. A target with only phenotypic
evidence must read differently from one with exact-target evidence, and the
only way to guarantee that is to never add them together.

A tier that could not be evaluated is marked not assessable and carries the
reason. That is not the same as a tier that was evaluated and found nothing.
"""

from dataclasses import dataclass, field

from bioactivity_records import BioactivityRecord, organism_matches, record_organism
from target_inputs import TargetSpec, ortholog_accessions

EXACT_TARGET = "exact_target"
ORTHOLOG = "ortholog"
FAMILY = "family"
PHENOTYPIC = "phenotypic"

# Precedence order. A record is attributed to exactly one tier per target,
# the strongest one it qualifies for, so the same measurement is never
# counted twice for the same target.
TIER_ORDER: tuple[str, ...] = (EXACT_TARGET, ORTHOLOG, FAMILY, PHENOTYPIC)

TIER_DEFINITIONS: dict[str, str] = {
    EXACT_TARGET: "bioactivity assayed against this protein in this organism",
    ORTHOLOG: "bioactivity against a homologous protein in another organism",
    FAMILY: "bioactivity against another protein in the same domain family",
    PHENOTYPIC: "whole-organism activity in this organism, with no protein assigned",
}

UNASSESSABLE_REASONS: dict[str, str] = {
    EXACT_TARGET: (
        "no UniProt accession for this target; add a uniprot_accession column "
        "to the target list"
    ),
    ORTHOLOG: "no ortholog mapping for this target; supply --orthologs",
    FAMILY: (
        "no Pfam accession for this target; add a pfam_accession column to the "
        "target list"
    ),
    PHENOTYPIC: (
        "no organism for this target; add an organism column to the target list "
        "or pass --organism"
    ),
}

NO_EVIDENCE_PROFILE = "no_evidence"
ONLY_SUFFIX = "_only"


@dataclass(frozen=True)
class BestActivity:
    """The most potent measurement in one tier, in one unit.

    Attributes:
        value: Activity converted into `unit`.
        unit: 'nM' or 'ug.mL-1'.
        activity_type: Measurement type behind the value, e.g. IC50 or MIC.
        compound_id: Compound the value belongs to.
        source_id: Assay or document the value came from.
    """

    value: float
    unit: str
    activity_type: str
    compound_id: str
    source_id: str


@dataclass
class TierEvidence:
    """Everything known about one target in one tier.

    Attributes:
        tier: Tier name.
        assessable: False when the target lacked the field this tier needs.
        unassessable_reason: Why the tier could not be evaluated; empty when
            it was.
        compound_ids: Distinct compounds with evidence in this tier, sorted.
        record_count: Measurements attributed to this tier.
        best: Most potent measurement per unit, keyed by unit. Molar and
            mass-concentration measurements are kept apart because converting
            between them needs a molecular weight the export does not carry.
        activity_types: Distinct measurement types seen, sorted.
        source_ids: Distinct assay or document identifiers, sorted.
        uncomparable_records: Measurements counted in the tier whose units
            could not be converted, so they contribute to the count but not
            to `best`.
    """

    tier: str
    assessable: bool
    unassessable_reason: str = ""
    compound_ids: list[str] = field(default_factory=list)
    record_count: int = 0
    best: dict[str, BestActivity] = field(default_factory=dict)
    activity_types: list[str] = field(default_factory=list)
    source_ids: list[str] = field(default_factory=list)
    uncomparable_records: int = 0


@dataclass
class TargetEvidence:
    """The four tiers for one target, plus what identified it.

    Attributes:
        protein_id: Identifier from the target list.
        gene_name: Gene symbol; may be empty.
        description: Free-text description; may be empty.
        accession: Accession used for the exact-target tier.
        accession_source: Where that accession came from. A zero in the
            exact-target tier means much less when the accession was a
            fallback from an internal identifier.
        organism: Organism used for the phenotypic tier.
        pfam_accession: Pfam accession used for the family tier.
        upstream_score: Score carried over from the target list, verbatim.
        upstream_decision: Decision carried over from the target list.
        tiers: One TierEvidence per tier, keyed by tier name.
        profile: Which tiers carry evidence, as a label. Never a number.
        not_assessable: Tiers that could not be evaluated, in tier order.
    """

    protein_id: str
    gene_name: str
    description: str
    accession: str
    accession_source: str
    organism: str
    pfam_accession: str
    upstream_score: str
    upstream_decision: str
    tiers: dict[str, TierEvidence]
    profile: str
    not_assessable: list[str] = field(default_factory=list)


def tier_is_assessable(tier: str, target: TargetSpec) -> bool:
    """Report whether a tier can be evaluated for a target.

    Args:
        tier: Tier name.
        target: The target.

    Returns:
        True when the target carries the field the tier needs.
    """
    if tier == EXACT_TARGET:
        return bool(target.accession)
    if tier == ORTHOLOG:
        return bool(ortholog_accessions(target))
    if tier == FAMILY:
        return bool(target.pfam_accession)
    if tier == PHENOTYPIC:
        return bool(target.organism)
    return False


def classify_record(record: BioactivityRecord, target: TargetSpec) -> str | None:
    """Decide which tier a record belongs to for one target.

    Tiers are tested in precedence order and the first match wins, so a
    measurement against the target's own protein is never also counted as
    family evidence.

    Args:
        record: The bioactivity record.
        target: The target.

    Returns:
        The tier name, or None when the record says nothing about this
        target.
    """
    if target.accession and record.target_accession == target.accession:
        return EXACT_TARGET
    if record.target_accession and record.target_accession in ortholog_accessions(target):
        return ORTHOLOG
    if target.pfam_accession and record.target_pfam == target.pfam_accession:
        return FAMILY
    if not record.target_accession and organism_matches(record_organism(record), target.organism):
        return PHENOTYPIC
    return None


def better_activity(candidate: BestActivity, incumbent: BestActivity | None) -> bool:
    """Decide whether a candidate measurement replaces the current best.

    Lower is more potent for every measurement type this skill converts
    (IC50, Ki, Kd, EC50, MIC), so the smallest value wins. Ties resolve on
    compound then source identifier, which keeps repeated runs identical.

    Args:
        candidate: Measurement under consideration.
        incumbent: Current best, or None.

    Returns:
        True when the candidate should replace the incumbent.
    """
    if incumbent is None:
        return True
    return (candidate.value, candidate.compound_id, candidate.source_id) < (
        incumbent.value,
        incumbent.compound_id,
        incumbent.source_id,
    )


def summarise_tier(tier: str, records: list[BioactivityRecord]) -> TierEvidence:
    """Summarise the records attributed to one assessable tier.

    Args:
        tier: Tier name.
        records: Records attributed to this tier for one target.

    Returns:
        The tier summary. Counts and best values describe this tier alone and
        are never merged with another tier's.
    """
    evidence = TierEvidence(tier=tier, assessable=True, record_count=len(records))
    compounds: set[str] = set()
    activity_types: set[str] = set()
    sources: set[str] = set()
    for record in records:
        if record.compound_id:
            compounds.add(record.compound_id)
        if record.activity_type:
            activity_types.add(record.activity_type)
        if record.source_id:
            sources.add(record.source_id)
        if record.comparable_value is None or not record.comparable_unit:
            evidence.uncomparable_records += 1
            continue
        candidate = BestActivity(
            value=record.comparable_value,
            unit=record.comparable_unit,
            activity_type=record.activity_type,
            compound_id=record.compound_id,
            source_id=record.source_id,
        )
        if better_activity(candidate, evidence.best.get(record.comparable_unit)):
            evidence.best[record.comparable_unit] = candidate
    evidence.compound_ids = sorted(compounds)
    evidence.activity_types = sorted(activity_types)
    evidence.source_ids = sorted(sources)
    return evidence


def unassessable_tier(tier: str) -> TierEvidence:
    """Build the summary for a tier that could not be evaluated.

    Args:
        tier: Tier name.

    Returns:
        A TierEvidence marked not assessable, carrying the reason. Its zero
        counts must be rendered as 'n/a', never as 'no compounds found'.
    """
    return TierEvidence(
        tier=tier,
        assessable=False,
        unassessable_reason=UNASSESSABLE_REASONS[tier],
    )


def evidence_profile(tiers: dict[str, TierEvidence]) -> str:
    """Label which tiers carry evidence.

    Args:
        tiers: The four tier summaries.

    Returns:
        A label such as 'exact_target_only', 'phenotypic_only',
        'ortholog+family' or 'no_evidence'. It is a label, not a rank: the
        point is that a phenotypic-only target reads differently at a glance
        from one with exact-target chemistry.
    """
    populated = [
        tier for tier in TIER_ORDER if tiers[tier].assessable and tiers[tier].compound_ids
    ]
    if not populated:
        return NO_EVIDENCE_PROFILE
    if len(populated) == 1:
        return populated[0] + ONLY_SUFFIX
    return "+".join(populated)


def gather_evidence(
    target: TargetSpec, records: list[BioactivityRecord]
) -> tuple[TargetEvidence, list[tuple[int, str]]]:
    """Assign every record to a tier for one target and summarise each tier.

    Args:
        target: The target.
        records: All bioactivity records in the export.

    Returns:
        Tuple of the target's evidence and the (record index, tier) pairs
        attributed to it, for the per-record output table. Indices refer to
        `records` and let the caller tell which records matched nothing at
        all.
    """
    attributed: list[tuple[int, str]] = []
    by_tier: dict[str, list[BioactivityRecord]] = {tier: [] for tier in TIER_ORDER}
    for index, record in enumerate(records):
        tier = classify_record(record, target)
        if tier is None:
            continue
        by_tier[tier].append(record)
        attributed.append((index, tier))

    tiers: dict[str, TierEvidence] = {}
    not_assessable: list[str] = []
    for tier in TIER_ORDER:
        if tier_is_assessable(tier, target):
            tiers[tier] = summarise_tier(tier, by_tier[tier])
        else:
            tiers[tier] = unassessable_tier(tier)
            not_assessable.append(tier)

    evidence = TargetEvidence(
        protein_id=target.protein_id,
        gene_name=target.gene_name,
        description=target.description,
        accession=target.accession,
        accession_source=target.accession_source,
        organism=target.organism,
        pfam_accession=target.pfam_accession,
        upstream_score=target.upstream_score,
        upstream_decision=target.upstream_decision,
        tiers=tiers,
        profile=evidence_profile(tiers),
        not_assessable=not_assessable,
    )
    return evidence, attributed


def compound_count(evidence: TargetEvidence, tier: str) -> int | None:
    """Return a tier's compound count, or None when the tier was not assessed.

    Args:
        evidence: The target's evidence.
        tier: Tier name.

    Returns:
        The number of distinct compounds, or None when the tier could not be
        evaluated. None must not be rendered as zero.
    """
    summary = evidence.tiers[tier]
    if not summary.assessable:
        return None
    return len(summary.compound_ids)
