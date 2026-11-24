"""Run progress page, resource history and estimates (advice only; never changes results)."""

from __future__ import annotations

import datetime as dt
import json
import os
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
    (run / "results/report").mkdir(parents=True)
    (run / "results/report/pfvarx_dashboard.html").write_text("report")
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


def test_stopped_run_is_reported_interrupted_without_a_report_link(tmp_path: Path) -> None:
    from pfvarx.runs.registry import run_status

    run = _run_dir(tmp_path, "RUNNING", None)
    status = json.loads((run / "status.json").read_text())
    import socket
    status.update(pid=2**22 + 12345, host=socket.gethostname())  # no such process
    (run / "status.json").write_text(json.dumps(status))
    assert run_status(run)["status"] == "INTERRUPTED" and "RESUME" in run_status(run)["interruption"]
    p = progress.write(run)
    assert p["finished"] and not p["report_ready"] and p["running"] == []
    page = (run / "progress.html").read_text()
    assert "The run did not finish" in page and "Open the report dashboard" not in page
    status.update(pid=os.getpid())  # this process is alive and is PfVarX's test run
    (run / "status.json").write_text(json.dumps(status))
    assert run_status(run)["status"] in ("RUNNING", "INTERRUPTED")  # INTERRUPTED only if not a pfvarx process


def test_measured_peak_raises_the_plan_request(monkeypatch) -> None:
    from pfvarx.planning import plan

    runs = [_hist_run("r1", 2.0, 1200.0)]
    monkeypatch.setattr(history, "load_history", lambda: runs)
    assert history.measured_peak_mb("pfvarx_call", "WGS", "BAM") == 4000.0
    assert history.measured_peak_mb("snpeff") is None
    from types import SimpleNamespace as NS

    from pfvarx.constants import Design, FinalOutput, InputType

    def bundle(per_job_mem: int, max_jobs: int) -> NS:
        return NS(compute=NS(max_jobs=max_jobs, max_cores=8, max_mem_mb=12000, tmp_disk_mb=10**6,
                             per_job=NS(max_threads=8, max_mem_mb=per_job_mem)),
                  config=NS(final_output=FinalOutput.FO002, input_type=InputType.BAM, assay_design=Design.WGS,
                            analyses=NS(selection="all"), reference=NS(gff=None, host_fasta_used=None)))

    design = bundle(1, 1).config.assay_design.value
    runs[0]["design"] = design
    monkeypatch.setattr(plan, "_available_mem_mb", lambda: 100000)
    rows = {r["job"]: r for r in plan.resource_plan(Path("."), bundle(6500, 1), 4.0, 2)["jobs"]}
    assert rows["pfvarx_call"]["mem_mb"] == 5000 and rows["pfvarx_call"]["memory_basis"] == "FRAMEWORK_POLICY_REQUEST"
    runs[0]["jobs"][0]["rss_mb"] = 5000.0  # 5000 MB measured: 6000 MB with headroom
    rows = {r["job"]: r for r in plan.resource_plan(Path("."), bundle(6500, 1), 4.0, 2)["jobs"]}
    assert rows["pfvarx_call"]["mem_mb"] == 6000 and rows["pfvarx_call"]["memory_basis"] == "MEASURED_HISTORY"
    tight = plan.resource_plan(Path("."), bundle(5500, 1), 4.0, 2)
    assert tight["feasible"] and any(w["reason"] == "RESOURCE_MEMORY_TIGHT" for w in tight["warnings"])
    monkeypatch.setattr(plan, "_available_mem_mb", lambda: 7000)
    busy = plan.resource_plan(Path("."), bundle(6500, 2), 4.0, 2)
    assert any(w["reason"] == "RESOURCE_MEMORY_IN_USE" and "max_jobs to 1" in w["detail"] for w in busy["warnings"])


def test_damaged_sample_sheet_header_fails_cleanly(tmp_path: Path) -> None:
    from pfvarx.constants import Design, InputType
    from pfvarx.inputs.validate import validate_samples

    sheet = tmp_path / "samples.tsv"
    sheet.write_text("\x1b]633;C\x07sample_id\tlane\tlibrary_id\tbam\nS1\tL1\tlib\tx.bam\n")
    rep, rows = validate_samples(sheet, tmp_path, InputType.BAM, Design.WGS, set())
    assert rep.status == "FAIL" and rows == []
    assert "hidden or extra characters" in rep.findings[0].detail
    sheet.write_text("﻿sample_id\tlane\tlibrary_id\tbam\nS1\tL1\tlib\tx.bam\n")  # byte-order mark is fine
    rep, rows = validate_samples(sheet, tmp_path, InputType.BAM, Design.WGS, set())
    assert rows and rows[0]["sample_id"] == "S1"
