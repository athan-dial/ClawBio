"""Score pathogen proteins as antifungal or antimicrobial target candidates.

Three components, 100 points: divergence from the host proteome (40),
published essentiality (35), and conservation across a pathogen panel (25).
A close host orthologue additionally carries a penalty, because a selectivity
liability is not offset by good biology elsewhere.

Every component scores an absent annotation as zero and records it. An
unrun search is not a negative result, and a score of 40 for a protein with
three blank fields means the work has not been done, not that the target is
poor.
"""

from dataclasses import dataclass, field

from homology_search import HomologyHit

# Host divergence bands, applied to the best host hit's percent identity.
SELECTIVITY_NO_HIT_POINTS = 40
IDENTITY_BANDS: tuple[tuple[float, int, str], ...] = (
    (25.0, 32, "distant host homologue"),
    (40.0, 20, "moderate host homologue"),
    (60.0, 8, "close host homologue; selectivity work required"),
)
SELECTIVITY_FLOOR_POINTS = 0
SAFETY_IDENTITY_CUTOFF = 60.0
SAFETY_PENALTY = 25

ESSENTIAL_POINTS = 35
ESSENTIAL_UNKNOWN_POINTS = 7
ESSENTIAL_CALLS = {"essential", "essential_in_vitro", "conditionally_essential"}
NON_ESSENTIAL_CALLS = {"non-essential", "non_essential", "nonessential", "no", "none", "dispensable"}
UNKNOWN_CALLS = {"", "unknown", "untested", "na", "n/a"}

CONSERVATION_ALL_POINTS = 25
CONSERVATION_MAJORITY_POINTS = 13
CONSERVATION_MAJORITY_THRESHOLD = 0.5

PURSUE_THRESHOLD = 70
HOLD_THRESHOLD = 45


@dataclass
class TargetVerdict:
    """Scored result for one candidate target.

    Attributes:
        protein_id: Pathogen protein identifier.
        gene_name: Gene symbol, as supplied.
        description: Free-text description, as supplied.
        score: Total score, 0-100 after clamping.
        decision: One of 'pursue', 'hold', 'drop'.
        host_identity_pct: Percent identity of the best host hit, or None.
        panel_present: Panel species carrying a homologue.
        panel_total: Panel species searched.
        reasons: Human-readable score components.
        missing_annotations: Inputs that were absent and scored as zero.
    """

    protein_id: str
    gene_name: str
    description: str
    score: int
    decision: str
    host_identity_pct: float | None
    panel_present: int
    panel_total: int
    reasons: list[str] = field(default_factory=list)
    missing_annotations: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)


def score_selectivity(host_hit: HomologyHit | None) -> tuple[int, str]:
    """Score divergence from the closest host protein.

    Args:
        host_hit: Best host hit, or None when the search found nothing.

    Returns:
        Tuple of points awarded and a reason string.
    """
    if host_hit is None:
        return SELECTIVITY_NO_HIT_POINTS, (
            f"no detectable host homologue: +{SELECTIVITY_NO_HIT_POINTS}"
        )
    for ceiling, points, label in IDENTITY_BANDS:
        if host_hit.identity_pct < ceiling:
            return points, f"{label} ({host_hit.identity_pct:.0f}% identity): +{points}"
    return SELECTIVITY_FLOOR_POINTS, (
        f"host orthologue at {host_hit.identity_pct:.0f}% identity: +{SELECTIVITY_FLOOR_POINTS}"
    )


def score_essentiality(call: str, source: str) -> tuple[int, str]:
    """Score published essentiality evidence.

    An essentiality call with no source is treated as unknown. A claim that a
    gene is essential is only as good as the screen behind it, and an
    unsourced call cannot be checked.

    Args:
        call: Essentiality call, e.g. 'essential', 'non-essential', 'unknown'.
        source: Screen or publication the call came from.

    Returns:
        Tuple of points awarded and a reason string.
    """
    normalised = call.strip().lower()
    if normalised in UNKNOWN_CALLS:
        return ESSENTIAL_UNKNOWN_POINTS, f"essentiality unknown: +{ESSENTIAL_UNKNOWN_POINTS}"
    if normalised in NON_ESSENTIAL_CALLS:
        return 0, "reported non-essential: +0"
    if normalised in ESSENTIAL_CALLS:
        if not source.strip():
            return ESSENTIAL_UNKNOWN_POINTS, (
                f"essential but no source given, scored as unknown: +{ESSENTIAL_UNKNOWN_POINTS}"
            )
        return ESSENTIAL_POINTS, f"essential ({source.strip()}): +{ESSENTIAL_POINTS}"
    return ESSENTIAL_UNKNOWN_POINTS, (
        f"unrecognised essentiality call '{call.strip()}', scored as unknown: "
        f"+{ESSENTIAL_UNKNOWN_POINTS}"
    )


def score_conservation(present: int, total: int) -> tuple[int, str]:
    """Score conservation across the pathogen panel.

    Args:
        present: Panel species in which a homologue was found.
        total: Panel species searched.

    Returns:
        Tuple of points awarded and a reason string.
    """
    if total <= 0:
        return 0, "no panel supplied, spectrum not assessed: +0"
    if present >= total:
        return CONSERVATION_ALL_POINTS, (
            f"present in all {total} panel species: +{CONSERVATION_ALL_POINTS}"
        )
    if present / total >= CONSERVATION_MAJORITY_THRESHOLD:
        return CONSERVATION_MAJORITY_POINTS, (
            f"present in {present}/{total} panel species: +{CONSERVATION_MAJORITY_POINTS}"
        )
    return 0, f"present in {present}/{total} panel species: +0"


def decide(score: int) -> str:
    """Convert a score into a decision.

    Args:
        score: Total target score.

    Returns:
        'pursue', 'hold', or 'drop'.
    """
    if score >= PURSUE_THRESHOLD:
        return "pursue"
    if score >= HOLD_THRESHOLD:
        return "hold"
    return "drop"


def score_target(
    protein_id: str,
    gene_name: str,
    description: str,
    host_hit: HomologyHit | None,
    essentiality: str,
    essentiality_source: str,
    panel_present: int,
    panel_total: int,
    host_searched: bool = True,
) -> TargetVerdict:
    """Score one pathogen protein as a target candidate.

    Args:
        protein_id: Pathogen protein identifier.
        gene_name: Gene symbol, may be empty.
        description: Free-text description, may be empty.
        host_hit: Best host hit, or None.
        essentiality: Essentiality call.
        essentiality_source: Screen or publication behind the call.
        panel_present: Panel species carrying a homologue.
        panel_total: Panel species searched.
        host_searched: False when no host search was run for this protein,
            which scores zero and is recorded, rather than being read as
            absence of a host homologue.

    Returns:
        A TargetVerdict.
    """
    reasons: list[str] = []
    missing: list[str] = []
    flags: list[str] = []
    score = 0

    if host_searched:
        selectivity_points, selectivity_reason = score_selectivity(host_hit)
    else:
        selectivity_points, selectivity_reason = 0, "host homology not assessed: +0"
        missing.append("host_homology")

    orphan = (
        host_searched
        and host_hit is None
        and panel_total > 0
        and panel_present == 0
    )
    if orphan:
        selectivity_points = 0
        selectivity_reason = (
            "no host homologue and none in any panel species: +0 "
            "(absence here means the protein is not conserved anywhere, "
            "not that it is selectively targetable)"
        )
        flags.append("orphan_candidate")
    reasons.append(selectivity_reason)
    score += selectivity_points

    essentiality_points, essentiality_reason = score_essentiality(essentiality, essentiality_source)
    reasons.append(essentiality_reason)
    score += essentiality_points
    if essentiality_points == ESSENTIAL_UNKNOWN_POINTS:
        missing.append("essentiality")

    conservation_points, conservation_reason = score_conservation(panel_present, panel_total)
    reasons.append(conservation_reason)
    score += conservation_points
    if panel_total <= 0:
        missing.append("conservation_panel")

    if host_searched and host_hit is not None and host_hit.identity_pct >= SAFETY_IDENTITY_CUTOFF:
        score -= SAFETY_PENALTY
        reasons.append(
            f"selectivity penalty, host identity at or above "
            f"{SAFETY_IDENTITY_CUTOFF:.0f}%: -{SAFETY_PENALTY}"
        )

    score = max(0, min(100, score))
    return TargetVerdict(
        protein_id=protein_id,
        gene_name=gene_name,
        description=description,
        score=score,
        decision=decide(score),
        host_identity_pct=host_hit.identity_pct if host_hit is not None else None,
        panel_present=panel_present,
        panel_total=panel_total,
        reasons=reasons,
        missing_annotations=missing,
        flags=flags,
    )
