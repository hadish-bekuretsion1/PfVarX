"""Cohort QC (gate G5): completeness in the cohort, and its effect on downstream eligibility."""
import gzip
from pathlib import Path

from pfvarx.eligibility.engine import SampleState, load_registry, sample_eligibility
from pfvarx.quality.cohort import evaluate, sample_metrics, site_summary

ROOT = Path(__file__).resolve().parents[2]
FROZEN = {"config": {"assay_design": "wgs"}, "quality": {"policy_id": "pfvarx_qc_v1"}}


def _cohort(tmp_path: Path) -> tuple[Path, Path]:
    # 10 sites. GOOD misses 1 (10 %), WARN misses 3 (30 %), BAD misses 6 (60 %).
    missing = tmp_path / "missingness.tsv"
    missing.write_text(
        "sample_id\tALT_DOMINANT\tMIXED\tREF_DOMINANT\tMISSING_LOW_DEPTH\tUNCALLABLE\tmissing_fraction\n"
        "GOOD\t2\t1\t6\t1\t0\t0.1000\n"
        "WARN\t2\t2\t3\t2\t1\t0.3000\n"
        "BAD\t1\t0\t3\t2\t4\t0.6000\n")
    matrix = tmp_path / "evidence_matrix.tsv.gz"
    states = {"GOOD": ["REF_DOMINANT"] * 9 + ["UNCALLABLE"],
              "WARN": ["REF_DOMINANT"] * 7 + ["UNCALLABLE"] * 3,
              "BAD": ["REF_DOMINANT"] * 4 + ["UNCALLABLE"] * 6}
    with gzip.open(matrix, "wt") as handle:
        handle.write("site_id\tsample_id\tstate\n")
        for sid, column in states.items():
            for i, state in enumerate(column):
                handle.write(f"S{i}\t{sid}\t{state}\n")
    return missing, matrix


def test_g5_verdicts_follow_missingness_thresholds(tmp_path: Path) -> None:
    missing, matrix = _cohort(tmp_path)
    result = evaluate(ROOT, FROZEN, {"design": "wgs"}, missing, matrix)
    verdicts = {sid: q["verdict"] for sid, q in result["samples"].items()}
    assert verdicts == {"GOOD": "PASS", "WARN": "WARNING", "BAD": "FAIL"}  # warn above 0.2, fail above 0.5
    assert result["samples"]["BAD"]["permitted"]["DA002"] is False
    assert result["samples"]["WARN"]["permitted"]["DA002"] is True  # a warning may continue downstream
    assert {d["gate_id"] for d in result["samples"]["GOOD"]["decisions"]} == {"G5.cohort_missingness",
                                                                           "G5.cohort_mixed_fraction"}
    assert "G5.cohort_missingness" in result["samples"]["GOOD"]["uncalibrated_gates"]  # provisional limits


def test_metrics_and_site_call_rate(tmp_path: Path) -> None:
    missing, matrix = _cohort(tmp_path)
    metrics = sample_metrics(missing)
    assert metrics["WARN"]["cohort_sites"] == 10 and metrics["WARN"]["cohort_mixed_fraction"] == 0.2
    sites = site_summary(matrix, 3)
    # Sites S0-S3 are called in all 3 samples, S4-S6 in 2 of 3, S7-S9 in at most 1 of 3.
    assert sites["sites"] == 10 and sites["sites_meeting_call_rate"] == 4


def test_cohort_qc_fail_excludes_from_downstream_analyses(tmp_path: Path) -> None:
    missing, matrix = _cohort(tmp_path)
    result = evaluate(ROOT, FROZEN, {"design": "wgs"}, missing, matrix)["samples"]
    contract = load_registry(ROOT)["analyses"]["DA003"]
    layer1 = {"verdict": "PASS", "permitted": {}}
    ctx = {"run_id": "R1", "project_id": "P", "policy_version": "pfvarx_qc_v1"}
    records = {sid: sample_eligibility("DA003", contract, SampleState(sample_id=sid, design="wgs", qc=layer1,
                                                                      cohort_qc=result[sid]), "SAMPLE", ctx)
               for sid in result}
    assert records["BAD"].eligibility == "INELIGIBLE" and "COHORT_QC_FAIL" in records["BAD"].reason_codes
    assert records["GOOD"].eligibility == "ELIGIBLE" and records["WARN"].eligibility == "ELIGIBLE"
    assert records["BAD"].evidence["cohort_qc_verdict"] == "FAIL"
