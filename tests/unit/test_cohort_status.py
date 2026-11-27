"""Layer 3 outcome when no sample is admitted: an explicit skip, never a crashed merge."""

import json
from pathlib import Path

from pfvarx.downstream.common import load_context
from pfvarx.downstream.run import write_cohort_status

REPO = Path(__file__).resolve().parents[2]


def _run(tmp: Path, verdicts: dict[str, str], cohort: list[str] | None) -> Path:
    root = tmp / "ws"
    root.mkdir()
    (root / "config").symlink_to(REPO / "config")
    (root / "inputs/metadata").mkdir(parents=True)
    (root / "inputs/metadata/samples.tsv").write_text("sample_id\n" + "\n".join(verdicts) + "\n")
    run = root / "runs" / "r1"
    (run / "results/qc/samples").mkdir(parents=True)
    eff = {"config": {"project_id": "P", "samples_tsv": "inputs/metadata/samples.tsv"}, "parameters": {}}
    plan = {"samples": [{"sample_id": s} for s in verdicts if s != "BADINPUT"],
            "samples_failed_validation": ["BADINPUT"] if "BADINPUT" in verdicts else [],
            "layers": [1, 2, 3, 4, 5], "analyses_requested": ["DA002", "DA004"], "analyses_internal": ["DA001"]}
    (run / "frozen_config.json").write_text(json.dumps(eff))
    (run / "frozen_plan.json").write_text(json.dumps(plan))
    for sid, verdict in verdicts.items():
        if sid == "BADINPUT":
            continue
        d = run / "results/qc/samples" / sid
        d.mkdir()
        (d / "qc_decision.json").write_text(json.dumps({"verdict": verdict, "permitted": {"L3": verdict != "FAIL"}}))
    if cohort is not None:
        (run / "results/cohort").mkdir(parents=True)
        (run / "results/cohort/cohort_samples.txt").write_text("\n".join(cohort) + "\n")
    return run


def test_empty_cohort_skips_every_requested_analysis(tmp_path: Path) -> None:
    run = _run(tmp_path, {"S1": "FAIL", "S2": "FAIL", "BADINPUT": "FAIL"}, None)
    write_cohort_status(load_context(run.parents[1], run, with_matrix=False), empty=True)
    status = json.loads((run / "results/cohort/cohort_status.json").read_text())
    assert status["status"] == "SKIPPED_NO_ELIGIBLE_SAMPLES"
    assert status["blocked_samples"] == {"S1": "QC_GATE_FAIL", "S2": "QC_GATE_FAIL",
                                         "BADINPUT": "INPUT_VALIDATION_FAILED"}
    assert status["layers_skipped"] == [3, 4, 5] and status["analyses_skipped"] == ["DA002", "DA004", "DA001"]
    for da in ("DA001", "DA002", "DA004"):
        manifest = json.loads((run / f"results/downstream/{da}/status_manifest.json").read_text())
        assert manifest["terminal_state"] == "SKIPPED_INELIGIBLE"
        assert manifest["reason_codes"] == ["NO_ELIGIBLE_SAMPLES", "NO_COHORT"]
        # Every registered sample stays in the denominator trail.
        assert manifest["denominators"]["requested"] == 3 and manifest["denominators"]["analysed"] == 0


def test_formed_cohort_is_recorded_without_touching_analyses(tmp_path: Path) -> None:
    run = _run(tmp_path, {"S1": "PASS", "S2": "FAIL"}, ["S1"])
    write_cohort_status(load_context(run.parents[1], run, with_matrix=False), empty=False)
    status = json.loads((run / "results/cohort/cohort_status.json").read_text())
    assert status["status"] == "COMPLETED"
    assert status["cohort_samples"] == ["S1"] and status["blocked_samples"] == {"S2": "QC_GATE_FAIL"}
    assert status["analyses_skipped"] == [] and not (run / "results/downstream").exists()
