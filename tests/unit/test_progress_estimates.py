"""Run progress page, resource history and estimates (advice only; never changes results)."""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

from pfvarx.execution import progress
from pfvarx.planning import history

RUN = "20261123T093234Z-TEST-0000abcd"
LOG = """
[Mon Nov 23 09:33:00 2026]
localrule alignment_entry:
    jobid: 4
    wildcards: sample=S1
[Mon Nov 23 09:40:00 2026]
Finished jobid: 4 (Rule: alignment_entry)
1 of 6 steps (17%) done

[Mon Nov 23 09:41:00 2026]
localrule pfvarx_call:
    jobid: 7
    wildcards: sample=S1
"""


def _hist_run(run_id: str, gb: float, call_s: float, design: str = "WGS", input_type: str = "BAM") -> dict:
    return {"run_id": run_id, "design": design, "input_type": input_type, "input_gb": gb,
            "jobs": [{"job": "pfvarx_call", "rule": "pfvarx_call", "sample": "A", "input_gb": gb,
                      "wall_s": call_s, "cpu_s": call_s, "rss_mb": 4000.0},
                     {"job": "cohort_merge", "rule": "cohort_merge", "sample": None, "input_gb": gb,
                      "wall_s": 60.0, "cpu_s": 60.0, "rss_mb": 200.0}],
            "storage_mb": {"variants": 10.0 * gb, "alignments": 900.0 * gb}}


def test_parse_log_tracks_running_finished_and_step_count() -> None:
    parsed = progress.parse_log(LOG)
    assert parsed["steps"] == {"done": 1, "total": 6, "percent": 17}
    assert parsed["jobs"]["4"]["status"] == "DONE"
    assert parsed["jobs"]["7"] == {**parsed["jobs"]["7"], "rule": "pfvarx_call", "sample": "S1", "status": "RUNNING"}


def test_parse_log_marks_failed_jobs() -> None:
    parsed = progress.parse_log(LOG + "\n[Mon Nov 23 10:00:00 2026]\nError in rule pfvarx_call:\n    jobid: 7\n")
    assert parsed["jobs"]["7"]["status"] == "FAILED"


def test_estimate_scales_by_input_size_and_packs_samples() -> None:
    runs = [_hist_run("r1", 2.0, 1200.0), _hist_run("r2", 4.0, 2400.0)]  # 600 s per GB
    est = history.estimate(["pfvarx_call", "cohort_merge", "dashboard"], {"S1": 1.0, "S2": 3.0}, 1, "WGS", "BAM",
                           history=runs, layers=[0, 1, 2, 3])
    call = next(s for s in est["steps"] if s["job"] == "pfvarx_call")
    assert call["minutes_per_sample"] == {"S1": 10.0, "S2": 30.0}
    assert call["peak_ram_mb"] == 4000 and call["confidence"] == "MEASURED_HISTORY"
    assert est["wall_minutes"] is not None and est["wall_minutes"] > 40  # one job at a time
    assert est["storage_mb"]["variants"] == 40 and "annotation" not in est["storage_mb"]
    two = history.estimate(["pfvarx_call"], {"S1": 1.0, "S2": 3.0}, 2, "WGS", "BAM", history=runs)
    assert two["wall_minutes"] == 30.0  # two lanes: the longest sample decides


def test_estimate_is_unknown_without_history() -> None:
    est = history.estimate(["pfvarx_call"], {"S1": 1.0}, 1, history=[])
    assert est["wall_minutes"] is None and est["confidence"] == "UNKNOWN"
    assert est["steps"][0]["minutes"] is None and "unknown" in history.estimate_table(est)


def test_history_file_follows_environment(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("PFVARX_HISTORY_DIR", str(tmp_path))
    assert history.history_file().parent == tmp_path
    assert history.load_history() == []


def _run_dir(tmp_path: Path, status: str, est: dict | None) -> Path:
    run = tmp_path / RUN
    (run / "logs").mkdir(parents=True)
    (run / "logs" / "snakemake.pass1.log").write_text(LOG)
    (run / "status.json").write_text(json.dumps({"status": status, "updated_utc": "2026-11-23T10:30:00+00:00"}))
    (run / "frozen_plan.json").write_text(json.dumps({"resources": {"estimate": est}}))
    return run


def test_progress_snapshot_running_with_estimate(tmp_path: Path) -> None:
    est = history.estimate(["alignment_entry", "pfvarx_call"], {"S1": 1.0}, 1, history=[
        {**_hist_run("r1", 1.0, 1800.0), "jobs": _hist_run("r1", 1.0, 1800.0)["jobs"]
         + [{"job": "alignment_entry", "rule": "alignment_entry", "sample": "A", "input_gb": 1.0,
             "wall_s": 420.0, "cpu_s": 1.0, "rss_mb": 100.0}]}])
    run = _run_dir(tmp_path, "RUNNING", est)
    # Snakemake writes local clock times; the caller started at 09:41 local time.
    now = progress._parse_time("Mon Nov 23 09:41:00 2026") + dt.timedelta(minutes=10)
    p = progress.snapshot(run, now)
    assert p["started_utc"].startswith("2026-11-23T09:32:34")
    assert p["running"][0]["rule"] == "pfvarx_call" and p["running"][0]["running_min"] == 10.0
    assert p["running"][0]["expected_min"] == 30.0
    assert p["remaining_min"] == 20.0  # the caller's estimate less the 10 minutes it has run
    page = progress.render(p)
    assert 'http-equiv="refresh" content="15"' in page and "Variant calling" in page


def test_progress_page_finished_links_report_and_never_opens_in_tests(tmp_path: Path) -> None:
    run = _run_dir(tmp_path, "COMPLETED", None)
    p = progress.write(run)
    assert p["finished"] and p["elapsed_min"] == 57.4 and p["remaining_min"] == 0.0
    page = (run / "progress.html").read_text()
    assert "results/report/pfvarx_dashboard.html" in page and "Open the report dashboard" in page
    assert json.loads((run / "progress.json").read_text())["status"] == "COMPLETED"
    assert progress.open_in_browser(run / "progress.html") is False  # PFVARX_NO_BROWSER is set by conftest


def test_gui_progress_page_waits_for_the_new_run(tmp_path: Path) -> None:
    from pfvarx.configuration.web import WAITING_PAGE, run_page

    assert run_page(tmp_path, "progress") == WAITING_PAGE
    run = _run_dir(tmp_path / "runs", "RUNNING", None)
    assert run_page(tmp_path, "progress", after=RUN) == WAITING_PAGE  # the previous run: keep waiting
    page = run_page(tmp_path, "progress", after="20250101T000000Z-OLD-0000abcd")
    assert page != WAITING_PAGE and RUN.encode() in page
    assert (run / "progress.html").is_file()
