"""Tests for target-chemistry-evidence.

Written before the implementation, per the repo's red/green rule.

The suite is organised around the one property the skill exists to hold: the
four evidence tiers are never collapsed into one number, and a tier that
could not be evaluated never looks like a tier that was evaluated and found
nothing.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bioactivity_records import (  # noqa: E402
    BioactivityRecord,
    accession_from_protein_id,
    is_uniprot_header_id,
    convert_activity,
    delimiter_for,
    normalise_accession,
    organism_matches,
    parse_activity_value,
    read_bioactivity,
    read_delimited,
    record_organism,
)
from evidence_tiers import (  # noqa: E402
    EXACT_TARGET,
    BestActivity,
    FAMILY,
    NO_EVIDENCE_PROFILE,
    ORTHOLOG,
    PHENOTYPIC,
    TIER_ORDER,
    better_activity,
    classify_record,
    compound_count,
    evidence_profile,
    gather_evidence,
    summarise_tier,
    tier_is_assessable,
)
from target_chemistry_evidence import (  # noqa: E402
    DEMO_DIR,
    build_demo_inputs,
    build_parser,
    format_count,
    resolve_inputs,
    run,
)
from target_inputs import (  # noqa: E402
    ACCESSION_FROM_COLUMN,
    ACCESSION_FROM_HEADER,
    ACCESSION_FROM_PROTEIN_ID,
    Ortholog,
    TargetSpec,
    attach_orthologs,
    ortholog_accessions,
    parse_identity,
    read_orthologs,
    read_targets,
)

BIOACTIVITY_HEADER = (
    "compound_id\tcompound_name\ttarget_id\ttarget_accession\ttarget_organism\t"
    "target_pfam\tactivity_type\tactivity_value\tactivity_units\tassay_type\t"
    "assay_organism\tsource_id\n"
)


def make_record(**kwargs) -> BioactivityRecord:
    """Build a record with sensible defaults for a single-field test.

    Args:
        **kwargs: Fields to override.

    Returns:
        The record.
    """
    base = {
        "compound_id": "CPD1",
        "compound_name": "compound one",
        "target_id": "TGT1",
        "target_accession": "P00001",
        "target_organism": "Candida albicans",
        "target_pfam": "PF00067",
        "activity_type": "IC50",
        "activity_value": 10.0,
        "activity_units": "nM",
        "comparable_value": 10.0,
        "comparable_unit": "nM",
        "assay_type": "B",
        "assay_organism": "Candida albicans",
        "source_id": "ASSAY1",
    }
    base.update(kwargs)
    return BioactivityRecord(**base)


def make_target(**kwargs) -> TargetSpec:
    """Build a target with sensible defaults for a single-field test.

    Args:
        **kwargs: Fields to override.

    Returns:
        The target.
    """
    base = {
        "protein_id": "CAL_X",
        "gene_name": "X",
        "description": "a protein",
        "accession": "P00001",
        "organism": "Candida albicans",
        "pfam_accession": "PF00067",
    }
    base.update(kwargs)
    return TargetSpec(**base)


# ---------------------------------------------------------------------------
# Accessions, organisms and units
# ---------------------------------------------------------------------------

class TestNormaliseAccession:
    def test_upper_cases(self):
        assert normalise_accession("p10635") == "P10635"

    def test_strips_whitespace(self):
        assert normalise_accession("  P10635 ") == "P10635"

    def test_drops_isoform_suffix(self):
        assert normalise_accession("P10635-2") == "P10635"

    def test_drops_version_suffix(self):
        assert normalise_accession("P10635.3") == "P10635"

    def test_blank_stays_blank(self):
        assert normalise_accession("   ") == ""


class TestAccessionFromProteinId:
    def test_uniprot_header_identifier(self):
        assert accession_from_protein_id("sp|P0CY33|ERG11_CANAL") == "P0CY33"

    def test_trembl_identifier(self):
        assert accession_from_protein_id("tr|A0A1D8PI61|A0A1D8PI61_CANAL") == "A0A1D8PI61"

    def test_bare_accession_passes_through(self):
        assert accession_from_protein_id("P0CY33") == "P0CY33"

    def test_internal_identifier_is_not_invented_into_an_accession(self):
        assert accession_from_protein_id("CAL_FKS1") == "CAL_FKS1"

    def test_blank(self):
        assert accession_from_protein_id("") == ""


class TestOrganismMatches:
    def test_exact_match(self):
        assert organism_matches("Candida albicans", "Candida albicans")

    def test_case_insensitive(self):
        assert organism_matches("candida albicans", "Candida Albicans")

    def test_strain_suffix_matches(self):
        assert organism_matches("Candida albicans SC5314", "Candida albicans")

    def test_different_species_do_not_match(self):
        assert not organism_matches("Candida albicans", "Candida glabrata")

    def test_prefix_without_word_boundary_does_not_match(self):
        assert not organism_matches("Candida albicansoid", "Candida albicans")

    def test_blank_never_matches(self):
        assert not organism_matches("", "Candida albicans")
        assert not organism_matches("Candida albicans", "")


class TestParseActivityValue:
    def test_number(self):
        assert parse_activity_value("42.5") == 42.5

    def test_blank_is_none(self):
        assert parse_activity_value("") is None

    def test_censored_value_is_not_read_as_exact(self):
        assert parse_activity_value(">100") is None

    def test_text_is_none(self):
        assert parse_activity_value("not determined") is None


class TestConvertActivity:
    @pytest.mark.parametrize("units,expected", [
        ("nM", 5.0), ("uM", 5000.0), ("mM", 5_000_000.0), ("M", 5e9), ("pM", 0.005),
    ])
    def test_molar_family_converts_to_nanomolar(self, units, expected):
        value, unit = convert_activity(5.0, units)
        assert value == pytest.approx(expected)
        assert unit == "nM"

    @pytest.mark.parametrize("units,expected", [
        ("ug.mL-1", 5.0), ("ug/mL", 5.0), ("mg.L-1", 5.0), ("mg/L", 5.0),
        ("ng.mL-1", 0.005), ("mg/mL", 5000.0),
    ])
    def test_mass_family_converts_to_microgram_per_ml(self, units, expected):
        value, unit = convert_activity(5.0, units)
        assert value == pytest.approx(expected)
        assert unit == "ug.mL-1"

    def test_molar_and_mass_are_not_interconverted(self):
        """A MIC in ug.mL-1 has no molar equivalent without a molecular weight."""
        assert convert_activity(1.0, "ug.mL-1")[1] != convert_activity(1.0, "nM")[1]

    def test_unknown_units_are_not_guessed(self):
        assert convert_activity(62.0, "%") == (None, "")

    def test_missing_value_is_not_converted(self):
        assert convert_activity(None, "nM") == (None, "")

    def test_case_is_significant_for_molar_units(self):
        """Folding case would map mM onto M and move a value a millionfold."""
        assert convert_activity(1.0, "mM")[0] != convert_activity(1.0, "M")[0]


class TestDelimiterFor:
    @pytest.mark.parametrize("name,expected", [
        ("a.tsv", "\t"), ("a.tab", "\t"), ("a.txt", "\t"), ("a.csv", ","), ("a.dat", ","),
    ])
    def test_extension_chooses_delimiter(self, name, expected):
        assert delimiter_for(Path(name)) == expected


class TestRecordOrganism:
    def test_assay_organism_wins(self):
        record = make_record(assay_organism="Candida albicans SC5314", target_organism="x")
        assert record_organism(record) == "Candida albicans SC5314"

    def test_falls_back_to_target_organism(self):
        assert record_organism(make_record(assay_organism="")) == "Candida albicans"


# ---------------------------------------------------------------------------
# Reading input files
# ---------------------------------------------------------------------------

class TestReadDelimited:
    def test_reads_csv(self, tmp_path):
        path = tmp_path / "a.csv"
        path.write_text("a,b\n1,2\n")
        rows, fieldnames = read_delimited(path)
        assert rows == [{"a": "1", "b": "2"}]
        assert fieldnames == ["a", "b"]

    def test_reads_tsv(self, tmp_path):
        path = tmp_path / "a.tsv"
        path.write_text("a\tb\n1\t2\n")
        rows, _ = read_delimited(path)
        assert rows == [{"a": "1", "b": "2"}]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(SystemExit):
            read_delimited(tmp_path / "nope.csv")


class TestReadBioactivity:
    def test_parses_records(self, tmp_path):
        path = tmp_path / "b.tsv"
        path.write_text(
            BIOACTIVITY_HEADER
            + "CPD1\tone\tTGT1\tP00001\tCandida albicans\tPF00067\tIC50\t42\tnM\tB\t"
              "Candida albicans\tASSAY1\n"
        )
        records = read_bioactivity(path)
        assert len(records) == 1
        assert records[0].compound_id == "CPD1"
        assert records[0].comparable_value == 42.0
        assert records[0].comparable_unit == "nM"

    def test_missing_required_column_is_refused(self, tmp_path):
        path = tmp_path / "b.tsv"
        path.write_text("compound_id\ttarget_accession\n" "CPD1\tP00001\n")
        with pytest.raises(SystemExit) as error:
            read_bioactivity(path)
        assert "missing required columns" in str(error.value)

    def test_rows_without_a_compound_are_dropped(self, tmp_path):
        path = tmp_path / "b.tsv"
        path.write_text(
            BIOACTIVITY_HEADER
            + "\t\t\tP00001\tCandida albicans\t\tIC50\t42\tnM\tB\tCandida albicans\tA1\n"
        )
        assert read_bioactivity(path) == []

    def test_blank_activity_value_does_not_crash(self, tmp_path):
        path = tmp_path / "b.tsv"
        path.write_text(
            BIOACTIVITY_HEADER
            + "CPD1\tone\t\t\tCandida albicans\t\tMIC\t\tug.mL-1\tF\tCandida albicans\tA1\n"
        )
        record = read_bioactivity(path)[0]
        assert record.activity_value is None
        assert record.comparable_value is None

    def test_accession_is_normalised_on_read(self, tmp_path):
        path = tmp_path / "b.tsv"
        path.write_text(
            BIOACTIVITY_HEADER
            + "CPD1\tone\tTGT1\tp00001-2\tCandida albicans\t\tIC50\t1\tnM\tB\t"
              "Candida albicans\tA1\n"
        )
        assert read_bioactivity(path)[0].target_accession == "P00001"

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(SystemExit):
            read_bioactivity(tmp_path / "nope.tsv")


class TestReadTargets:
    def test_reads_minimal_target_list(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("protein_id\nP00001\n")
        targets = read_targets(path)
        assert targets[0].protein_id == "P00001"
        assert targets[0].accession == "P00001"

    def test_target_scores_csv_from_upstream_is_accepted(self, tmp_path):
        path = tmp_path / "target_scores.csv"
        path.write_text(
            "protein_id,gene_name,description,score,decision,host_identity_pct,"
            "panel_present,panel_total,flags,missing_annotations,reasons\n"
            'CAL_FKS1,FKS1,"beta-1,3-glucan synthase",100,pursue,,3,3,,,"a; b"\n'
        )
        target = read_targets(path)[0]
        assert target.gene_name == "FKS1"
        assert target.upstream_score == "100"
        assert target.upstream_decision == "pursue"
        assert target.description == "beta-1,3-glucan synthase"

    def test_explicit_accession_wins_over_derived(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("protein_id,uniprot_accession\nCAL_FKS1,P12345\n")
        assert read_targets(path)[0].accession == "P12345"

    def test_default_organism_fills_blank_rows(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("protein_id,organism\nA,\nB,Candida glabrata\n")
        targets = read_targets(path, "Candida albicans")
        assert targets[0].organism == "Candida albicans"
        assert targets[1].organism == "Candida glabrata"

    def test_missing_protein_id_column_is_refused(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("gene_name\nFKS1\n")
        with pytest.raises(SystemExit):
            read_targets(path)

    def test_file_with_no_usable_rows_is_refused(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("protein_id\n\n")
        with pytest.raises(SystemExit):
            read_targets(path)


class TestReadOrthologs:
    def test_groups_by_target(self, tmp_path):
        path = tmp_path / "o.csv"
        path.write_text(
            "protein_id,ortholog_accession,ortholog_organism,identity_pct,source\n"
            "A,P1,Homo sapiens,58,search\n"
            "A,P2,Mus musculus,55,search\n"
            "B,P3,Homo sapiens,31,search\n"
        )
        mapping = read_orthologs(path)
        assert len(mapping["A"]) == 2
        assert mapping["B"][0].identity_pct == 31.0

    def test_rows_without_an_accession_are_skipped(self, tmp_path):
        path = tmp_path / "o.csv"
        path.write_text("protein_id,ortholog_accession\nA,\n")
        assert read_orthologs(path) == {}

    def test_missing_column_is_refused(self, tmp_path):
        path = tmp_path / "o.csv"
        path.write_text("protein_id\nA\n")
        with pytest.raises(SystemExit):
            read_orthologs(path)

    def test_attach_puts_orthologs_on_the_target(self):
        targets = [make_target(protein_id="A")]
        attach_orthologs(targets, {"A": [Ortholog("P1", "Homo sapiens", 58.0, "s")]})
        assert ortholog_accessions(targets[0]) == {"P1"}

    def test_attach_leaves_unmapped_targets_empty(self):
        targets = [make_target(protein_id="A")]
        attach_orthologs(targets, {})
        assert ortholog_accessions(targets[0]) == frozenset()


class TestParseIdentity:
    def test_number(self):
        assert parse_identity("58.5") == 58.5

    def test_blank_is_none(self):
        assert parse_identity("") is None

    def test_text_is_none(self):
        assert parse_identity("high") is None


# ---------------------------------------------------------------------------
# Tier assignment
# ---------------------------------------------------------------------------

class TestTierIsAssessable:
    def test_exact_needs_an_accession(self):
        assert tier_is_assessable(EXACT_TARGET, make_target())
        assert not tier_is_assessable(EXACT_TARGET, make_target(accession=""))

    def test_ortholog_needs_a_mapping(self):
        assert not tier_is_assessable(ORTHOLOG, make_target())
        target = make_target()
        target.orthologs = [Ortholog("P9", "Homo sapiens", 50.0, "s")]
        assert tier_is_assessable(ORTHOLOG, target)

    def test_family_needs_a_pfam_accession(self):
        assert tier_is_assessable(FAMILY, make_target())
        assert not tier_is_assessable(FAMILY, make_target(pfam_accession=""))

    def test_phenotypic_needs_an_organism(self):
        assert tier_is_assessable(PHENOTYPIC, make_target())
        assert not tier_is_assessable(PHENOTYPIC, make_target(organism=""))


class TestClassifyRecord:
    def test_same_accession_is_exact_target(self):
        assert classify_record(make_record(), make_target()) == EXACT_TARGET

    def test_ortholog_accession_is_ortholog(self):
        target = make_target()
        target.orthologs = [Ortholog("P00009", "Homo sapiens", 58.0, "s")]
        record = make_record(target_accession="P00009", target_organism="Homo sapiens")
        assert classify_record(record, target) == ORTHOLOG

    def test_same_pfam_different_protein_is_family(self):
        record = make_record(target_accession="P99999", target_organism="Mus musculus")
        assert classify_record(record, make_target()) == FAMILY

    def test_organism_assay_without_an_accession_is_phenotypic(self):
        record = make_record(target_accession="", target_pfam="")
        assert classify_record(record, make_target()) == PHENOTYPIC

    def test_unrelated_record_matches_nothing(self):
        record = make_record(
            target_accession="P99999", target_pfam="PF99999", target_organism="Homo sapiens",
            assay_organism="Homo sapiens",
        )
        assert classify_record(record, make_target()) is None

    def test_exact_beats_family_so_nothing_is_double_counted(self):
        """The record shares the target's accession and its Pfam."""
        assert classify_record(make_record(), make_target()) == EXACT_TARGET

    def test_ortholog_beats_family(self):
        target = make_target()
        target.orthologs = [Ortholog("P00009", "Homo sapiens", 58.0, "s")]
        record = make_record(target_accession="P00009", target_organism="Homo sapiens")
        assert classify_record(record, target) == ORTHOLOG

    def test_protein_assay_in_the_same_organism_is_not_phenotypic(self):
        record = make_record(target_accession="P99999", target_pfam="")
        assert classify_record(record, make_target()) is None

    def test_phenotypic_requires_the_organism_to_match(self):
        record = make_record(
            target_accession="", target_pfam="", target_organism="Candida glabrata",
            assay_organism="Candida glabrata",
        )
        assert classify_record(record, make_target()) is None

    def test_phenotypic_tolerates_a_strain_suffix(self):
        record = make_record(
            target_accession="", target_pfam="", assay_organism="Candida albicans SC5314"
        )
        assert classify_record(record, make_target()) == PHENOTYPIC

    def test_target_without_an_accession_has_no_exact_tier(self):
        record = make_record(target_accession="")
        assert classify_record(record, make_target(accession="")) != EXACT_TARGET


class TestSummariseTier:
    def test_counts_distinct_compounds_not_records(self):
        records = [make_record(source_id="A1"), make_record(source_id="A2")]
        summary = summarise_tier(EXACT_TARGET, records)
        assert len(summary.compound_ids) == 1
        assert summary.record_count == 2

    def test_best_is_the_most_potent(self):
        records = [
            make_record(compound_id="C1", comparable_value=100.0),
            make_record(compound_id="C2", comparable_value=4.0),
        ]
        summary = summarise_tier(EXACT_TARGET, records)
        assert summary.best["nM"].compound_id == "C2"
        assert summary.best["nM"].value == 4.0

    def test_units_are_kept_apart(self):
        records = [
            make_record(compound_id="C1", comparable_value=4.0, comparable_unit="nM"),
            make_record(
                compound_id="C2", comparable_value=0.5, comparable_unit="ug.mL-1",
                activity_type="MIC", activity_units="ug.mL-1",
            ),
        ]
        summary = summarise_tier(PHENOTYPIC, records)
        assert set(summary.best) == {"nM", "ug.mL-1"}
        assert summary.best["ug.mL-1"].value == 0.5

    def test_unconvertible_records_count_but_do_not_set_best(self):
        records = [
            make_record(compound_id="C1", comparable_value=None, comparable_unit=""),
        ]
        summary = summarise_tier(PHENOTYPIC, records)
        assert summary.uncomparable_records == 1
        assert summary.best == {}
        assert summary.compound_ids == ["C1"]

    def test_activity_types_and_sources_are_listed(self):
        records = [
            make_record(activity_type="IC50", source_id="A1"),
            make_record(activity_type="Ki", source_id="A2"),
        ]
        summary = summarise_tier(EXACT_TARGET, records)
        assert summary.activity_types == ["IC50", "Ki"]
        assert summary.source_ids == ["A1", "A2"]

    def test_empty_tier_is_assessable_and_empty(self):
        summary = summarise_tier(EXACT_TARGET, [])
        assert summary.assessable
        assert summary.compound_ids == []
        assert summary.record_count == 0


class TestBetterActivity:
    def test_lower_value_wins(self):
        incumbent = BestActivity(10.0, "nM", "IC50", "C1", "A1")
        candidate = BestActivity(2.0, "nM", "IC50", "C2", "A2")
        assert better_activity(candidate, incumbent)

    def test_ties_resolve_deterministically(self):
        incumbent = BestActivity(10.0, "nM", "IC50", "C2", "A1")
        candidate = BestActivity(10.0, "nM", "IC50", "C1", "A1")
        assert better_activity(candidate, incumbent)
        assert not better_activity(incumbent, candidate)

    def test_any_candidate_beats_nothing(self):
        assert better_activity(BestActivity(10.0, "nM", "IC50", "C1", "A1"), None)


# ---------------------------------------------------------------------------
# Tier separation — the property the skill exists to hold
# ---------------------------------------------------------------------------

class TestTierSeparation:
    def test_profile_names_the_tiers_with_evidence(self):
        target = make_target()
        target.orthologs = [Ortholog("P00009", "Homo sapiens", 58.0, "s")]
        records = [make_record(), make_record(compound_id="C9", target_accession="P00009")]
        evidence, _ = gather_evidence(target, records)
        assert evidence.profile == "exact_target+ortholog"

    def test_phenotypic_only_target_is_labelled_as_such(self):
        records = [make_record(compound_id="C1", target_accession="", target_pfam="")]
        evidence, _ = gather_evidence(make_target(pfam_accession="PF99999"), records)
        assert evidence.profile == "phenotypic_only"

    def test_phenotypic_only_reads_differently_from_exact_target(self):
        """The requirement in one assertion: the two must not look alike."""
        phenotypic_records = [make_record(compound_id="C1", target_accession="", target_pfam="")]
        exact_records = [make_record(compound_id="C1")]
        phenotypic, _ = gather_evidence(make_target(pfam_accession="PF99999"), phenotypic_records)
        exact, _ = gather_evidence(make_target(pfam_accession="PF99999"), exact_records)
        assert phenotypic.profile != exact.profile
        assert phenotypic.tiers[EXACT_TARGET].compound_ids == []
        assert exact.tiers[EXACT_TARGET].compound_ids == ["C1"]

    def test_no_evidence_profile(self):
        record = make_record(
            target_accession="P99999", target_pfam="PF99999",
            target_organism="Homo sapiens", assay_organism="Homo sapiens",
        )
        evidence, _ = gather_evidence(make_target(), [record])
        assert evidence.profile == NO_EVIDENCE_PROFILE

    def test_evidence_profile_of_empty_tiers(self):
        evidence, _ = gather_evidence(make_target(), [])
        assert evidence_profile(evidence.tiers) == NO_EVIDENCE_PROFILE

    def test_a_record_is_attributed_to_exactly_one_tier(self):
        target = make_target()
        target.orthologs = [Ortholog("P00009", "Homo sapiens", 58.0, "s")]
        evidence, attributed = gather_evidence(target, [make_record()])
        assert len(attributed) == 1
        counts = [len(evidence.tiers[tier].compound_ids) for tier in TIER_ORDER]
        assert sum(counts) == 1

    def test_target_evidence_carries_no_total_or_score_field(self):
        evidence, _ = gather_evidence(make_target(), [make_record()])
        fields = set(vars(evidence))
        assert "total" not in fields
        assert "druggability_score" not in fields
        assert "combined_score" not in fields

    def test_unassessable_tier_is_not_reported_as_zero(self):
        evidence, _ = gather_evidence(make_target(), [make_record()])
        assert compound_count(evidence, ORTHOLOG) is None
        assert compound_count(evidence, EXACT_TARGET) == 1

    def test_unassessable_tier_carries_a_reason(self):
        evidence, _ = gather_evidence(make_target(), [])
        assert "--orthologs" in evidence.tiers[ORTHOLOG].unassessable_reason
        assert ORTHOLOG in evidence.not_assessable

    def test_unassessable_tier_never_carries_evidence(self):
        evidence, _ = gather_evidence(make_target(pfam_accession=""), [make_record()])
        assert not evidence.tiers[FAMILY].assessable
        assert evidence.tiers[FAMILY].compound_ids == []

    def test_format_count_renders_not_assessed_rather_than_zero(self):
        assert format_count(None) == "n/a"
        assert format_count(0) == "0"


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def demo_output(tmp_path_factory):
    """Run the skill once on the built-in demo inputs."""
    output = tmp_path_factory.mktemp("tce")
    run(output_dir=output, argv=["--demo"], **build_demo_inputs())
    return output


@pytest.fixture(scope="module")
def demo_result(demo_output):
    """Parse the demo run's result.json."""
    return json.loads((demo_output / "result.json").read_text())


def targets_by_gene(payload: dict) -> dict:
    """Index a result payload's targets by gene name.

    Args:
        payload: Parsed result.json.

    Returns:
        Mapping of gene name to the target's payload.
    """
    return {target["gene_name"]: target for target in payload["targets"]}


class TestOutputContract:
    @pytest.mark.parametrize("relative", [
        "report.md",
        "result.json",
        "tables/target_evidence_summary.csv",
        "tables/evidence_by_tier.csv",
        "tables/compound_evidence.csv",
        "reproducibility/commands.sh",
        "reproducibility/environment.yml",
        "reproducibility/checksums.sha256",
    ])
    def test_artifact_written(self, demo_output, relative):
        assert (demo_output / relative).exists()

    def test_report_carries_disclaimer(self, demo_output):
        assert "Not a medical device" in (demo_output / "report.md").read_text()

    def test_demo_files_exist_on_disk(self):
        assert (DEMO_DIR / "demo_targets.csv").exists()
        assert (DEMO_DIR / "demo_bioactivity.tsv").exists()
        assert (DEMO_DIR / "demo_orthologs.csv").exists()

    def test_every_target_reports_all_four_tiers(self, demo_result):
        for target in demo_result["targets"]:
            assert set(target["tiers"]) == set(TIER_ORDER)

    def test_result_declares_that_tiers_are_not_combined(self, demo_result):
        assert demo_result["design"]["tiers_combined"] is False

    def test_result_carries_no_blended_score_anywhere(self, demo_result):
        """No key in a target payload may be a number summarising all tiers."""
        forbidden = {"druggability_score", "combined_score", "evidence_score", "total"}
        for target in demo_result["targets"]:
            assert forbidden.isdisjoint(target)
            assert forbidden.isdisjoint(target["tiers"])

    def test_target_order_is_preserved_from_the_input(self, demo_result):
        assert [target["protein_id"] for target in demo_result["targets"]] == [
            target.protein_id for target in build_demo_inputs()["targets"]
        ]

    def test_run_is_deterministic(self, tmp_path):
        first = tmp_path / "a"
        second = tmp_path / "b"
        run(output_dir=first, argv=["--demo"], **build_demo_inputs())
        run(output_dir=second, argv=["--demo"], **build_demo_inputs())
        assert (first / "result.json").read_text() == (second / "result.json").read_text()


class TestDemoExpectations:
    def test_target_with_protein_assays_has_exact_evidence(self, demo_result):
        erg11 = targets_by_gene(demo_result)["ERG11"]
        assert erg11["tiers"][EXACT_TARGET]["compound_count"] == 3
        assert erg11["tiers"][EXACT_TARGET]["best"]["nM"]["value"] == 42.0

    def test_target_with_only_organism_data_is_phenotypic_only(self, demo_result):
        fks1 = targets_by_gene(demo_result)["FKS1"]
        assert fks1["profile"] == "phenotypic_only"
        assert fks1["tiers"][EXACT_TARGET]["compound_count"] == 0
        assert fks1["tiers"][PHENOTYPIC]["compound_count"] > 0

    def test_phenotypic_only_and_exact_target_are_distinguishable_in_the_report(
        self, demo_output
    ):
        report = (demo_output / "report.md").read_text()
        assert "phenotypic_only" in report
        assert "exact_target+ortholog+family+phenotypic" in report

    def test_target_with_only_a_human_homologue_shows_ortholog_not_exact(self, demo_result):
        hsp90 = targets_by_gene(demo_result)["HSP90"]
        assert hsp90["tiers"][EXACT_TARGET]["compound_count"] == 0
        assert hsp90["tiers"][ORTHOLOG]["compound_count"] == 2

    def test_target_with_only_family_chemistry(self, demo_result):
        chs1 = targets_by_gene(demo_result)["CHS1"]
        assert chs1["tiers"][FAMILY]["compound_count"] == 2
        assert chs1["tiers"][ORTHOLOG]["compound_count"] is None

    def test_target_in_another_organism_has_no_evidence(self, demo_result):
        erg3 = targets_by_gene(demo_result)["ERG3"]
        assert erg3["profile"] == NO_EVIDENCE_PROFILE
        assert erg3["tiers"][PHENOTYPIC]["compound_count"] == 0

    def test_unassessable_ortholog_tier_reports_none_not_zero(self, demo_result):
        fks1 = targets_by_gene(demo_result)["FKS1"]
        assert fks1["tiers"][ORTHOLOG]["assessable"] is False
        assert fks1["tiers"][ORTHOLOG]["compound_count"] is None

    def test_mic_stays_in_mass_units(self, demo_result):
        fks1 = targets_by_gene(demo_result)["FKS1"]
        assert "ug.mL-1" in fks1["tiers"][PHENOTYPIC]["best"]
        assert "nM" not in fks1["tiers"][PHENOTYPIC]["best"]

    def test_unconvertible_records_are_counted_not_dropped(self, demo_result):
        fks1 = targets_by_gene(demo_result)["FKS1"]
        assert fks1["tiers"][PHENOTYPIC]["uncomparable_records"] == 2

    def test_records_outside_the_target_list_are_reported(self, demo_result):
        assert demo_result["records_unattributed"] == 5
        assert demo_result["records_read"] == 24

    def test_descriptions_with_commas_survive_parsing(self, demo_result):
        fks1 = targets_by_gene(demo_result)["FKS1"]
        assert fks1["description"] == "beta-1,3-glucan synthase catalytic subunit"

    def test_tier_table_names_its_tier_on_every_row(self, demo_output):
        lines = (demo_output / "tables" / "evidence_by_tier.csv").read_text().splitlines()
        assert lines[0].split(",")[2] == "tier"
        assert len(lines) == 1 + 7 * len(TIER_ORDER)

    def test_summary_table_uses_n_a_for_unassessable_tiers(self, demo_output):
        text = (demo_output / "tables" / "target_evidence_summary.csv").read_text()
        assert "n/a" in text

    def test_compound_table_records_the_tier_of_every_measurement(self, demo_output):
        lines = (demo_output / "tables" / "compound_evidence.csv").read_text().splitlines()
        tiers = {line.split(",")[1] for line in lines[1:]}
        assert tiers.issubset(set(TIER_ORDER))


class TestReportHonesty:
    def test_report_states_that_sparse_annotation_is_expected(self, demo_output):
        report = (demo_output / "report.md").read_text()
        assert "sparse" in report
        assert "not a failure of the search" in report

    def test_report_states_the_ortholog_selectivity_tension(self, demo_output):
        report = (demo_output / "report.md").read_text()
        assert "selectivity risk" in report

    def test_report_states_that_tiers_are_never_combined(self, demo_output):
        report = (demo_output / "report.md").read_text()
        assert "never summed" in report

    def test_report_explains_what_n_a_means(self, demo_output):
        report = (demo_output / "report.md").read_text()
        assert "not the same as no compounds found" in report

    def test_report_names_the_unassessable_tier_reason(self, demo_output):
        report = (demo_output / "report.md").read_text()
        assert "--orthologs" in report


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

class TestCommandLine:
    def test_demo_needs_no_other_input(self, tmp_path):
        args = build_parser().parse_args(["--demo", "--output", str(tmp_path)])
        inputs = resolve_inputs(args)
        assert inputs["targets"]
        assert inputs["records"]

    def test_missing_targets_is_refused(self, tmp_path):
        args = build_parser().parse_args(["--output", str(tmp_path)])
        with pytest.raises(SystemExit) as error:
            resolve_inputs(args)
        assert "--targets" in str(error.value)

    def test_missing_bioactivity_is_refused_with_the_offline_reason(self, tmp_path):
        targets = tmp_path / "t.csv"
        targets.write_text("protein_id\nP00001\n")
        args = build_parser().parse_args(
            ["--targets", str(targets), "--output", str(tmp_path)]
        )
        with pytest.raises(SystemExit) as error:
            resolve_inputs(args)
        assert "no network access" in str(error.value)

    def test_orthologs_are_optional(self, tmp_path):
        targets = tmp_path / "t.csv"
        targets.write_text("protein_id\nP00001\n")
        bioactivity = tmp_path / "b.tsv"
        bioactivity.write_text(BIOACTIVITY_HEADER)
        args = build_parser().parse_args([
            "--targets", str(targets), "--bioactivity", str(bioactivity),
            "--output", str(tmp_path),
        ])
        assert resolve_inputs(args)["records"] == []

    def test_organism_option_reaches_the_phenotypic_tier(self, tmp_path):
        targets = tmp_path / "t.csv"
        targets.write_text("protein_id\nP00001\n")
        bioactivity = tmp_path / "b.tsv"
        bioactivity.write_text(
            BIOACTIVITY_HEADER
            + "CPD1\tone\t\t\tCandida albicans\t\tMIC\t1\tug.mL-1\tF\tCandida albicans\tA1\n"
        )
        args = build_parser().parse_args([
            "--targets", str(targets), "--bioactivity", str(bioactivity),
            "--organism", "Candida albicans", "--output", str(tmp_path / "out"),
        ])
        result = run(output_dir=tmp_path / "out", argv=[], **resolve_inputs(args))
        assert result["targets"][0]["tiers"][PHENOTYPIC]["compound_count"] == 1

    def test_run_with_an_empty_export_still_writes_every_artifact(self, tmp_path):
        targets = tmp_path / "t.csv"
        targets.write_text("protein_id\nP00001\n")
        bioactivity = tmp_path / "b.tsv"
        bioactivity.write_text(BIOACTIVITY_HEADER)
        args = build_parser().parse_args([
            "--targets", str(targets), "--bioactivity", str(bioactivity),
            "--output", str(tmp_path / "out"),
        ])
        run(output_dir=tmp_path / "out", argv=[], **resolve_inputs(args))
        assert (tmp_path / "out" / "report.md").exists()
        assert (tmp_path / "out" / "tables" / "compound_evidence.csv").exists()

    def test_no_organism_makes_the_phenotypic_tier_unassessable(self, tmp_path):
        targets = tmp_path / "t.csv"
        targets.write_text("protein_id\nP00001\n")
        bioactivity = tmp_path / "b.tsv"
        bioactivity.write_text(
            BIOACTIVITY_HEADER
            + "CPD1\tone\t\t\tCandida albicans\t\tMIC\t1\tug.mL-1\tF\tCandida albicans\tA1\n"
        )
        args = build_parser().parse_args([
            "--targets", str(targets), "--bioactivity", str(bioactivity),
            "--output", str(tmp_path / "out"),
        ])
        result = run(output_dir=tmp_path / "out", argv=[], **resolve_inputs(args))
        phenotypic = result["targets"][0]["tiers"][PHENOTYPIC]
        assert phenotypic["assessable"] is False
        assert phenotypic["compound_count"] is None


# ---------------------------------------------------------------------------
# Accession provenance
#
# An empty exact-target tier only says something when the accession could
# have matched in the first place, so where it came from is reported too.
# ---------------------------------------------------------------------------

class TestIsUniprotHeaderId:
    def test_swissprot_header(self):
        assert is_uniprot_header_id("sp|P0CY33|ERG11_CANAL")

    def test_trembl_header(self):
        assert is_uniprot_header_id("tr|A0A1D8PI61|A0A1D8PI61_CANAL")

    def test_bare_accession_is_not_a_header(self):
        assert not is_uniprot_header_id("P0CY33")

    def test_internal_identifier_is_not_a_header(self):
        assert not is_uniprot_header_id("CAL_FKS1")

    def test_header_with_a_blank_accession_is_rejected(self):
        assert not is_uniprot_header_id("sp||ERG11_CANAL")


class TestAccessionProvenance:
    def test_explicit_column_is_recorded(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("protein_id,uniprot_accession\nCAL_FKS1,P12345\n")
        assert read_targets(path)[0].accession_source == ACCESSION_FROM_COLUMN

    def test_uniprot_header_is_recorded(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("protein_id\nsp|P0CY33|ERG11_CANAL\n")
        assert read_targets(path)[0].accession_source == ACCESSION_FROM_HEADER

    def test_internal_identifier_is_recorded_as_a_fallback(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("protein_id\nCAL_FKS1\n")
        assert read_targets(path)[0].accession_source == ACCESSION_FROM_PROTEIN_ID

    def test_provenance_reaches_result_json(self, demo_result):
        for target in demo_result["targets"]:
            assert target["accession_source"] == ACCESSION_FROM_COLUMN

    def test_provenance_reaches_the_summary_table(self, demo_output):
        path = demo_output / "tables" / "target_evidence_summary.csv"
        assert "accession_source" in path.read_text().splitlines()[0]

    def test_report_warns_when_a_fallback_accession_found_nothing(self, tmp_path):
        targets = tmp_path / "t.csv"
        targets.write_text("protein_id\nCAL_FKS1\n")
        bioactivity = tmp_path / "b.tsv"
        bioactivity.write_text(BIOACTIVITY_HEADER)
        args = build_parser().parse_args([
            "--targets", str(targets), "--bioactivity", str(bioactivity),
            "--output", str(tmp_path / "out"),
        ])
        run(output_dir=tmp_path / "out", argv=[], **resolve_inputs(args))
        report = (tmp_path / "out" / "report.md").read_text()
        assert "those zeros carry no information at all" in report

    def test_report_is_silent_when_the_accession_was_supplied(self, demo_output):
        report = (demo_output / "report.md").read_text()
        assert "those zeros carry no information" not in report
