"""DA009 pfhrp2/pfhrp3 joint deletion: both evidence models, joint states, prevalence and trend."""

from __future__ import annotations

import random
from pathlib import Path

import pysam
import pytest

from pfvarx.constants import HrpEvidence, JointDeletion
from pfvarx.downstream import hrp

LENGTH = 20000
PFHRP2 = ("ctg", 5000, 5300)   # panel target with no reads: the gene's amplicon dropped out
PFHRP3 = ("ctg", 9000, 9300)   # panel target with reads: present
CONTROLS = [("ctg", 1000, 1300), ("ctg", 2000, 2300), ("ctg", 3000, 3300)]
PARAMS = {"da009_dropout_max_ratio": 0.05, "da009_present_min_ratio": 0.2, "da009_min_control_depth": 20,
          "da009_min_control_amplicons": 2, "da009_min_samples_per_unit": 3, "da009_interval_level": 0.95}
DELETED, PRESENT, UNKNOWN = hrp.DELETED, hrp.PRESENT, hrp.UNKNOWN


@pytest.fixture(scope="module")
def panel_bam(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Amplicon-like coverage: reads over pfhrp3 and the controls, none over pfhrp2."""
    d = tmp_path_factory.mktemp("hrp")
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"},
                                              "SQ": [{"SN": "ctg", "LN": LENGTH}]})
    seq = "".join(random.Random(3).choices("ACGT", k=100))
    reads = []
    n = 0
    for chrom, start, end in [PFHRP3, *CONTROLS]:
        for offset in range(-100, end - start, 4):  # even coverage across the whole target
            n += 1
            a = pysam.AlignedSegment(header)
            a.query_name, a.reference_id, a.reference_start = f"r{n}", 0, start + offset
            a.query_sequence, a.cigarstring, a.mapping_quality = seq, "100M", 60
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            reads.append(a)
    raw = d / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    path = d / "panel.bam"
    pysam.sort("-o", str(path), str(raw))
    pysam.index(str(path))
    return path


def _joint(sample: str, state: str, population: str = "P1", country: str = "C1", month: str = "2026-01",
           evidence: str = HrpEvidence.DEPTH_RATIO.value) -> dict[str, object]:
    return {"sample_id": sample, "joint_state": state, "population": population, "country": country,
            "month": month, "evidence_class": evidence}


def test_joint_state_never_assumes_a_gene_is_present() -> None:
    assert hrp.joint_state(DELETED, DELETED) == JointDeletion.DOUBLE_DELETION.value
    assert hrp.joint_state(DELETED, PRESENT) == JointDeletion.PFHRP2_ONLY.value
    assert hrp.joint_state(PRESENT, DELETED) == JointDeletion.PFHRP3_ONLY.value
    assert hrp.joint_state(PRESENT, PRESENT) == JointDeletion.NO_DELETION.value
    # An unjudged gene is indeterminate, whatever the other gene shows.
    for other in (DELETED, PRESENT, UNKNOWN):
        assert hrp.joint_state(UNKNOWN, other) == JointDeletion.INDETERMINATE.value
        assert hrp.joint_state(other, UNKNOWN) == JointDeletion.INDETERMINATE.value


def test_measures_follow_the_joint_state() -> None:
    assert hrp.carries(JointDeletion.DOUBLE_DELETION.value, "double_deletion") is True
    assert hrp.carries(JointDeletion.PFHRP2_ONLY.value, "double_deletion") is False
    assert hrp.carries(JointDeletion.PFHRP2_ONLY.value, "pfhrp2_deleted") is True
    assert hrp.carries(JointDeletion.PFHRP3_ONLY.value, "pfhrp2_deleted") is False
    assert hrp.carries(JointDeletion.PFHRP3_ONLY.value, "any_deletion") is True
    assert hrp.carries(JointDeletion.NO_DELETION.value, "any_deletion") is False
    assert hrp.carries(JointDeletion.INDETERMINATE.value, "any_deletion") is None


def test_amplicon_dropout_reads_the_panel(panel_bam: Path) -> None:
    with pysam.AlignmentFile(str(panel_bam)) as bam:
        control = hrp.median_target_depth(bam, CONTROLS, 20)
        gone = hrp.median_target_depth(bam, [PFHRP2], 20)
        kept = hrp.median_target_depth(bam, [PFHRP3], 20)
    assert control is not None and control >= 20
    assert gone == 0.0 and kept is not None and kept >= 20
    assert hrp.dropout_call(gone, control, len(CONTROLS), PARAMS)[0] == DELETED
    assert hrp.dropout_call(kept, control, len(CONTROLS), PARAMS)[0] == PRESENT


def test_dropout_without_working_controls_is_unknown_not_a_deletion() -> None:
    # The assay did not work for this sample: a missing amplicon says nothing.
    state, reason, _ = hrp.dropout_call(0.0, 4.0, 3, PARAMS)
    assert (state, reason) == (UNKNOWN, "CONTROL_AMPLICON_DEPTH_LOW")
    # Too few controls in the panel to tell a deletion from a failed PCR.
    assert hrp.dropout_call(0.0, 100.0, 1, PARAMS)[:2] == (UNKNOWN, "TOO_FEW_CONTROL_AMPLICONS")
    # The panel carries no target for this gene.
    assert hrp.dropout_call(None, 100.0, 3, PARAMS)[:2] == (UNKNOWN, "GENE_TARGET_NOT_IN_PANEL")
    # Between the two ratios nothing is claimed.
    assert hrp.dropout_call(10.0, 100.0, 3, PARAMS)[:2] == (UNKNOWN, "AMBIGUOUS_COVERAGE")


def test_panel_targets_split_into_gene_and_control_targets() -> None:
    places = {"pfhrp2_deletion": PFHRP2, "pfhrp3_deletion": PFHRP3}
    gene_targets, controls = hrp.split_targets([PFHRP2, PFHRP3, *CONTROLS], places)
    assert gene_targets["pfhrp2_deletion"] == [PFHRP2] and gene_targets["pfhrp3_deletion"] == [PFHRP3]
    assert controls == CONTROLS
    # A gene that could not be placed in the GFF claims no target.
    gene_targets, controls = hrp.split_targets([PFHRP2, *CONTROLS], {"pfhrp2_deletion": None})
    assert gene_targets["pfhrp2_deletion"] == [] and len(controls) == 4


def test_prevalence_excludes_indeterminate_from_the_denominator() -> None:
    rows = [_joint("A", JointDeletion.DOUBLE_DELETION.value), _joint("B", JointDeletion.PFHRP2_ONLY.value),
            _joint("C", JointDeletion.NO_DELETION.value), _joint("D", JointDeletion.INDETERMINATE.value)]
    pooled = {r["measure"]: r for r in hrp.prevalence(rows, 0.95, by_evidence=False) if r["group"] == "ALL"}
    assert pooled["pfhrp2_deleted"]["judged"] == 3 and pooled["pfhrp2_deleted"]["carriers"] == 2
    assert pooled["double_deletion"]["carriers"] == 1
    assert pooled["any_deletion"]["prevalence"] == pytest.approx(2 / 3)
    assert pooled["any_deletion"]["indeterminate"] == 1
    lo, hi = pooled["any_deletion"]["ci_lower"], pooled["any_deletion"]["ci_upper"]
    assert lo is not None and hi is not None and lo < 2 / 3 < hi


def test_prevalence_records_and_can_stratify_the_evidence_classes() -> None:
    rows = [_joint("A", JointDeletion.DOUBLE_DELETION.value, evidence=HrpEvidence.DEPTH_RATIO.value),
            _joint("B", JointDeletion.NO_DELETION.value, evidence=HrpEvidence.AMPLICON_DROPOUT.value)]
    pooled = [r for r in hrp.prevalence(rows, 0.95, by_evidence=False)
              if r["group"] == "ALL" and r["measure"] == "double_deletion"][0]
    # The pooled row is the headline figure, but it always says what it pooled.
    assert pooled["judged"] == 2 and pooled["pooled_evidence_classes"] == 2
    assert pooled["evidence_classes"] == "AMPLICON_DROPOUT,DEPTH_RATIO"
    split = {r["evidence_class"]: r for r in hrp.prevalence(rows, 0.95, by_evidence=True)
             if r["group"] == "ALL" and r["measure"] == "double_deletion"}
    assert split[HrpEvidence.DEPTH_RATIO.value]["prevalence"] == 1.0
    assert split[HrpEvidence.AMPLICON_DROPOUT.value]["prevalence"] == 0.0


def test_cooccurrence_against_independence() -> None:
    # Deletions always together: observed double deletions exceed the independent expectation.
    rows = [_joint(f"D{i}", JointDeletion.DOUBLE_DELETION.value) for i in range(6)] + \
           [_joint(f"N{i}", JointDeletion.NO_DELETION.value) for i in range(6)]
    row = [r for r in hrp.cooccurrence(rows, 0.95) if r["group"] == "ALL"][0]
    assert row["both_deleted"] == 6 and row["judged"] == 12
    assert row["expected_both_if_independent"] == pytest.approx(3.0)
    assert row["p_value"] is not None and row["p_value"] < 0.01
    # Without variation in one gene no association can be estimated, and none is invented.
    flat = [_joint(f"P{i}", JointDeletion.NO_DELETION.value) for i in range(5)]
    only = [r for r in hrp.cooccurrence(flat, 0.95) if r["group"] == "ALL"][0]
    assert only["status"] == "NO_VARIATION" and only["p_value"] is None


def test_trend_withholds_prevalence_below_the_minimum() -> None:
    rows = [_joint("A", JointDeletion.DOUBLE_DELETION.value, month="2026-01"),
            _joint("B", JointDeletion.NO_DELETION.value, month="2026-01"),
            _joint("C", JointDeletion.NO_DELETION.value, month="2026-01"),
            _joint("D", JointDeletion.DOUBLE_DELETION.value, month="2026-02")]
    by_unit = {(r["month"], r["measure"]): r for r in hrp.trend(rows, 3, 0.95)}
    full = by_unit[("2026-01", "double_deletion")]
    assert full["judged"] == 3 and full["prevalence"] == pytest.approx(1 / 3) and full["status"] == "COMPUTED"
    thin = by_unit[("2026-02", "double_deletion")]
    assert thin["judged"] == 1 and thin["prevalence"] is None and thin["status"] == "INSUFFICIENT_SAMPLES"
    assert thin["carriers"] == 1  # the counts are still reported


def test_joint_rows_pair_the_two_genes_per_sample() -> None:
    genes = [{"locus_id": "pfhrp2_deletion", "gene_name": "pfhrp2", "gene_id": "G2"},
             {"locus_id": "pfhrp3_deletion", "gene_name": "pfhrp3", "gene_id": "G3"}]
    calls = [{"sample_id": "S1", "locus_id": "pfhrp2_deletion", "state": DELETED, "reason": "",
              "evidence_class": HrpEvidence.AMPLICON_DROPOUT.value, "design": "targeted_umi",
              "population": "P1", "country": "C1", "month": "2026-01"},
             {"sample_id": "S1", "locus_id": "pfhrp3_deletion", "state": UNKNOWN, "reason": "AMBIGUOUS_COVERAGE",
              "evidence_class": HrpEvidence.AMPLICON_DROPOUT.value, "design": "targeted_umi",
              "population": "P1", "country": "C1", "month": "2026-01"}]
    (row,) = hrp.joint_rows(calls, genes)
    assert row["pfhrp2_state"] == DELETED and row["pfhrp3_state"] == UNKNOWN
    assert row["joint_state"] == JointDeletion.INDETERMINATE.value
    assert row["indeterminate_reason"] == "AMBIGUOUS_COVERAGE"


def test_catalogue_and_contract_are_checked() -> None:
    catalogue = [{"locus_id": "pfhrp2_deletion", "event": "DELETION", "gene_id": "G2", "gene_name": "pfhrp2"},
                 {"locus_id": "pfhrp3_deletion", "event": "DELETION", "gene_id": "G3", "gene_name": "pfhrp3"},
                 {"locus_id": "pfmdr1_amplification", "event": "AMPLIFICATION", "gene_id": "G1", "gene_name": "pfmdr1"}]
    assert [g["gene_name"] for g in hrp.gene_rows(catalogue, ["pfhrp2_deletion", "pfhrp3_deletion"])] == \
        ["pfhrp2", "pfhrp3"]
    with pytest.raises(ValueError, match="no DELETION locus"):
        hrp.gene_rows(catalogue, ["pfhrp2_deletion", "pfmdr1_amplification"])
    with pytest.raises(ValueError, match="exactly two"):
        hrp.contract_loci({"needs": {"deletion_loci": ["pfhrp2_deletion"]}})


def test_shipped_registry_and_catalogue_agree_on_the_two_genes() -> None:
    from pfvarx.downstream.resistance import load_catalog
    from pfvarx.eligibility.engine import load_registry

    root = Path(__file__).resolve().parents[2]
    contract = load_registry(root)["analyses"]["DA009"]
    _, catalogue = load_catalog(root / "config/registry/cnv_catalog_pfvarx_v1.tsv")
    genes = hrp.gene_rows(catalogue, hrp.contract_loci(contract))
    assert [g["gene_name"] for g in genes] == ["pfhrp2", "pfhrp3"]
    assert sorted(contract["designs"]) == ["targeted_no_umi", "targeted_umi", "wgs"]


def test_interval_level_and_month_parsing() -> None:
    assert hrp.z_score(0.95) == pytest.approx(1.959963984540054)
    assert hrp.z_score(0.99) == pytest.approx(2.5758293035489004)
    assert hrp._month("2026-01-17") == "2026-01"
    assert hrp._month("") == "" and hrp._month("not-a-date") == ""


# --- evidence that needs no calibrated cut-off, and how much a call depends on one

def test_background_percentile_places_a_gene_among_the_samples_own_windows() -> None:
    from pfvarx.downstream.cnv import background_percentile

    # The sample's own single-copy core windows, centred on 1.
    windows = sorted([0.8, 0.9, 0.95, 1.0, 1.0, 1.05, 1.1, 1.2, 1.3, 1.4])
    assert background_percentile(windows, 0.0) == 0.0      # a deleted gene sits below every window
    assert background_percentile(windows, 1.0) == pytest.approx(0.5)   # a single copy sits in the middle
    assert background_percentile(windows, 2.0) == 1.0      # a duplication sits above them all
    assert background_percentile(windows, None) is None    # nothing measured, nothing claimed
    assert background_percentile([], 0.5) is None          # no background, no percentile


def test_background_ratios_come_from_the_samples_own_model() -> None:
    from pfvarx.downstream.cnv import GcModel, background_ratios

    model = GcModel.__new__(GcModel)        # the points a sample's windows produced
    model.points = [(0.2, 20.0), (0.25, 22.0), (0.3, 18.0), (0.35, 21.0), (0.4, 19.0)]
    model.min_per_bin = 3
    ratios = background_ratios(model)
    assert len(ratios) == 5 and ratios == sorted(ratios)
    # Single-copy windows against their own expectation: ratios cluster about 1.
    assert all(0.5 < r < 2.0 for r in ratios)


def _call(state: str, **extra: object) -> dict[str, object]:
    base = {"sample_id": "S1", "locus_id": "pfhrp2_deletion", "gene": "pfhrp2", "state": state,
            "reason": "", "evidence_class": HrpEvidence.AMPLICON_DROPOUT.value}
    return {**base, **extra}


def test_sensitivity_marks_a_call_that_depends_on_the_cut_off() -> None:
    params = {**PARAMS, "da009_dropout_ratio_grid": [0.02, 0.05, 0.10, 0.20]}
    # 8% of the control depth: deleted at a 0.10 or 0.20 cut-off, present-or-unknown below it.
    borderline = _call(UNKNOWN, reason="AMBIGUOUS_COVERAGE", gene_median_depth="8.0",
                       control_median_depth="100.0", control_amplicons=3)
    rows, unstable = hrp.sensitivity([borderline], HrpEvidence.AMPLICON_DROPOUT, params)
    assert unstable == 1, "a call that moves within the grid must be counted"
    assert {r["threshold"] for r in rows} == {0.02, 0.05, 0.10, 0.20}
    assert {r["state"] for r in rows} == {DELETED, UNKNOWN}
    assert all(r["reported_state"] == UNKNOWN for r in rows)   # nothing reported changed
    # A gene with no reads at all is deleted at every cut-off in the grid: the call stands alone.
    clear = _call(DELETED, gene_median_depth="0.0", control_median_depth="100.0", control_amplicons=3)
    rows, unstable = hrp.sensitivity([clear], HrpEvidence.AMPLICON_DROPOUT, params)
    assert unstable == 0 and {r["state"] for r in rows} == {DELETED}


def test_sensitivity_for_whole_genome_uses_the_depth_ratio_grid() -> None:
    params = {**PARAMS, "da009_deletion_ratio_grid": [0.05, 0.10, 0.20, 0.30],
              "da008_deletion_ratio": 0.1, "da008_present_ratio": 0.3, "da008_flank_min_ratio": 0.3,
              "da008_amplification_ratio": 1.5, "da008_amplification_lower": 1.25}
    call = {"sample_id": "S1", "locus_id": "pfhrp2_deletion", "gene": "pfhrp2", "state": UNKNOWN,
            "reason": "AMBIGUOUS_DEPTH", "evidence_class": HrpEvidence.DEPTH_RATIO.value,
            "depth_ratio": "0.1500", "ratio_lo": "0.1000", "ratio_hi": "0.2000", "flank_ratio": "0.9000"}
    rows, unstable = hrp.sensitivity([call], HrpEvidence.DEPTH_RATIO, params)
    assert len(rows) == 4 and unstable == 1
    by_threshold = {r["threshold"]: r["state"] for r in rows}
    assert by_threshold[0.05] == UNKNOWN and by_threshold[0.20] == DELETED


def test_sensitivity_leaves_an_unmeasured_locus_alone() -> None:
    params = {**PARAMS, "da009_deletion_ratio_grid": [0.05, 0.2]}
    missing = {"sample_id": "S1", "locus_id": "pfhrp3_deletion", "gene": "pfhrp3", "state": UNKNOWN,
               "reason": "GENE_NOT_IN_GFF", "evidence_class": HrpEvidence.DEPTH_RATIO.value}
    rows, unstable = hrp.sensitivity([missing], HrpEvidence.DEPTH_RATIO, params)
    assert unstable == 0 and {r["state"] for r in rows} == {UNKNOWN}


def test_background_spread_shows_whether_the_states_separate() -> None:
    calls = [
        {"gene": "pfhrp2", "state": DELETED, "background_percentile": "0.0010"},
        {"gene": "pfhrp2", "state": DELETED, "background_percentile": "0.0000"},
        {"gene": "pfhrp2", "state": PRESENT, "background_percentile": "0.4800"},
        {"gene": "pfhrp2", "state": PRESENT, "background_percentile": "0.5200"},
        {"gene": "pfhrp3", "state": UNKNOWN, "background_percentile": ""},   # nothing measured
    ]
    spread = hrp.background_spread(calls)
    assert spread["pfhrp2:DELETED"]["samples"] == 2 and spread["pfhrp2:DELETED"]["max"] == 0.001
    assert spread["pfhrp2:PRESENT"]["median"] == pytest.approx(0.5)
    # The deleted samples sit far below the present ones: the locus separates on the data itself.
    assert spread["pfhrp2:DELETED"]["max"] < spread["pfhrp2:PRESENT"]["min"]
    assert "pfhrp3:UNKNOWN" not in spread, "a locus with no measurement claims no spread"
