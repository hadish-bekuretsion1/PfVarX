"""Live run progress: what is running now, how far the run is, and when it should finish.

While Snakemake runs, the runner calls :func:`write` every few seconds. It reads the run's
Snakemake log (jobs started, finished and failed, and "N of M steps done") and the estimate
frozen with the plan (learned from earlier runs on this computer, see
``pfvarx.planning.history``), and writes ``runs/<RUN_ID>/progress.json`` and a self-refreshing
``runs/<RUN_ID>/progress.html``. When the run ends the page links to, and then opens, the
report dashboard. The page only reads run files; it never changes the run.
"""

from __future__ import annotations

import datetime as dt
import html
import json
import os
import re
import sys
import webbrowser
from pathlib import Path
from typing import Any

from pfvarx.planning.history import JOB_NAMES, PER_SAMPLE, RULE_TO_JOB

# Rule names in the Snakemake log that differ from their benchmark folder names.
RULE_ALIASES = {"bwa_mem2_lane": "bwa_mem2", "mark_duplicates": "markdup", "mosdepth_targets": "mosdepth",
                "entry_trim_primers": "trim_primers", "da_execute": "downstream", "da001_clonality": "downstream",
                "all": None}
STEP_LABELS = {
    "fastp": "Read QC and trimming", "bwa_mem2": "Read alignment", "merge_lanes": "Merge lanes",
    "markdup": "Mark duplicates", "alignment_entry": "Imported BAM checks", "entry_finalize": "Imported BAM finalise",
    "trim_primers": "Primer trimming", "flagstat": "Alignment statistics", "mosdepth_core": "Core-genome depth",
    "mosdepth": "Target depth", "l1_gate": "QC gate", "pfvarx_call": "Variant calling",
    "normalise_sample": "Normalise sample VCF", "cohort_merge": "Cohort merge", "cohort_status": "Cohort status",
    "cohort_evidence": "Cohort evidence", "cohort_qc": "Cohort QC", "snpeff_build": "Annotation database",
    "snpeff_annotate": "Annotation", "da_eligibility": "Analysis eligibility", "downstream": "Downstream analysis",
}
TIME = re.compile(r"^\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]$")
HEADER = re.compile(r"^(?:local)?(?:rule|checkpoint) (\w+):$")
FINISHED = re.compile(r"^Finished jobid: (\d+)")
ERROR = re.compile(r"^Error in rule (\w+):")
STEPS = re.compile(r"^(\d+) of (\d+) steps \((\d+)%\) done")


def _parse_time(text: str) -> dt.datetime:
    return dt.datetime.strptime(" ".join(text.split()), "%a %b %d %H:%M:%S %Y").astimezone(dt.UTC)


def parse_log(text: str) -> dict[str, Any]:
    """Jobs started, finished and failed, and the latest step count, from a Snakemake log."""
    jobs: dict[str, dict[str, Any]] = {}
    when: dt.datetime | None = None
    current: dict[str, Any] | None = None
    error_rule: str | None = None
    steps = None
    for raw in text.splitlines():
        line = raw.strip()
        if m := TIME.match(line):
            when = _parse_time(m.group(1))
            current = None
        elif m := HEADER.match(line):
            current = {"rule": m.group(1), "started": when, "sample": None, "status": "RUNNING"}
        elif line.startswith("jobid:") and (current is not None or error_rule):
            jobid = line.split(":", 1)[1].strip()
            if error_rule:
                jobs.setdefault(jobid, {"rule": error_rule, "started": None, "sample": None})["status"] = "FAILED"
                error_rule = None
            elif current is not None:
                jobs[jobid] = current
        elif line.startswith("wildcards:") and current is not None:
            found = re.search(r"\bsample=([^,\s]+)", line)
            da = re.search(r"\bda=([^,\s]+)", line)
            current["sample"] = found.group(1) if found else (da.group(1) if da else None)
        elif m := FINISHED.match(line):
            job = jobs.get(m.group(1))
            if job is not None:
                job["status"], job["finished"] = "DONE", when
        elif m := ERROR.match(line):
            error_rule = m.group(1)
        elif m := STEPS.match(line):
            steps = {"done": int(m.group(1)), "total": int(m.group(2)), "percent": int(m.group(3))}
    return {"jobs": jobs, "steps": steps}


def _job_class(rule: str) -> str | None:
    folder = RULE_ALIASES.get(rule, rule)
    return RULE_TO_JOB.get(folder) if folder else None


def snapshot(run: Path, now: dt.datetime | None = None) -> dict[str, Any]:
    now = now or dt.datetime.now(dt.UTC)
    status = json.loads((run / "status.json").read_text()) if (run / "status.json").is_file() else {}
    plan = json.loads((run / "frozen_plan.json").read_text()) if (run / "frozen_plan.json").is_file() else {}
    estimate = (plan.get("resources") or {}).get("estimate") or {}
    started = dt.datetime.strptime(run.name[:16], "%Y%m%dT%H%M%SZ").replace(tzinfo=dt.UTC)
    logs = sorted((run / "logs").glob("snakemake.pass*.log"))
    parsed = parse_log("\n".join(p.read_text(errors="replace") for p in logs)) if logs else {"jobs": {}, "steps": None}
    jobs = parsed["jobs"].values()
    running = [j for j in jobs if j.get("status") == "RUNNING"]
    failed = [j for j in jobs if j.get("status") == "FAILED"]
    state = status.get("status", "RUNNING")
    finished = state not in ("RUNNING", "PENDING", "READY", None)
    if finished and running:  # stopped without finishing these jobs
        for j in running:
            j["status"] = "STOPPED"
        running = []
    # Expected minutes per (job class, sample) from the frozen estimate; cohort steps use sample None.
    expected: dict[tuple[str, str | None], float] = {}
    for step in estimate.get("steps", []):
        if step.get("minutes") is None:
            continue
        if step.get("per_sample"):
            for sample, minutes in (step.get("minutes_per_sample") or {}).items():
                expected[(step["job"], sample)] = minutes
        else:
            expected[(step["job"], None)] = step["minutes"]
    known = bool(expected) and not estimate.get("unknown_steps")
    # Remaining work: untouched (class, sample) pairs count in full, running ones the rest of their estimate.
    touched: dict[tuple[str, str | None], str] = {}
    running_rows = []
    for j in jobs:
        job = _job_class(j["rule"])
        if job is None:
            continue
        key = (job, j["sample"] if job in PER_SAMPLE else None)
        if j.get("status") == "RUNNING":
            touched[key] = "RUNNING"
        else:
            touched.setdefault(key, "DONE")
    for j in running:
        job = _job_class(j["rule"])
        run_key = (job, j["sample"] if job in PER_SAMPLE else None) if job else None
        elapsed = (now - j["started"]).total_seconds() / 60 if j.get("started") else None
        running_rows.append({"step": STEP_LABELS.get(RULE_ALIASES.get(j["rule"], j["rule"]) or "", j["rule"]),
                             "rule": j["rule"], "sample": j["sample"], "tool_group": JOB_NAMES.get(job or "", ""),
                             "started": j["started"].isoformat() if j.get("started") else None,
                             "running_min": None if elapsed is None else round(elapsed, 1),
                             "expected_min": expected.get(run_key) if run_key else None})
    remaining_sample: dict[str, float] = {}
    remaining_cohort = 0.0
    for (job, sample), minutes in expected.items():
        state_of = touched.get((job, sample))
        if state_of == "DONE":
            continue
        left = minutes
        if state_of == "RUNNING":
            spent = max((r["running_min"] or 0) for r in running_rows
                        if _job_class(r["rule"]) == job and (sample is None or r["sample"] == sample))
            left = max(0.0, minutes - spent)
        if sample is None:
            remaining_cohort += left
        else:
            remaining_sample[sample] = remaining_sample.get(sample, 0.0) + left
    lanes = [0.0] * max(1, min(int(estimate.get("max_jobs") or 1), max(len(remaining_sample), 1)))
    for minutes in sorted(remaining_sample.values(), reverse=True):
        lanes[lanes.index(min(lanes))] += minutes
    remaining = max(lanes) + remaining_cohort if known and not finished else (0.0 if finished else None)
    end = dt.datetime.fromisoformat(status["updated_utc"]) if finished and status.get("updated_utc") else now
    elapsed_min = (end - started).total_seconds() / 60
    total_expected = estimate.get("wall_minutes")
    by_step: dict[str, dict[str, Any]] = {}
    for step in estimate.get("steps", []):
        by_step[step["job"]] = {"name": step["name"], "expected_min": step.get("minutes"), "peak_ram_mb": step.get("peak_ram_mb"),
                                "done": 0, "running": 0, "failed": 0}
    for j in jobs:
        job = _job_class(j["rule"])
        if job is None:
            continue
        row = by_step.setdefault(job, {"name": JOB_NAMES.get(job, job), "expected_min": None, "peak_ram_mb": None,
                                       "done": 0, "running": 0, "failed": 0})
        row[{"DONE": "done", "RUNNING": "running", "FAILED": "failed"}.get(j.get("status"), "failed")] += 1
    return {"schema": "pfvarx.progress/1", "run_id": run.name, "status": state, "finished": finished,
            "started_utc": started.isoformat(), "updated_utc": now.isoformat(),
            "finished_utc": end.isoformat() if finished else None, "elapsed_min": round(elapsed_min, 1),
            "steps": parsed["steps"], "running": running_rows,
            "failed": [{"rule": j["rule"], "sample": j["sample"]} for j in failed],
            "expected_total_min": total_expected, "remaining_min": None if remaining is None else round(remaining, 1),
            "expected_finish_utc": None if remaining is None else (now + dt.timedelta(minutes=remaining)).isoformat(),
            "estimate_basis": estimate.get("basis") or "no estimate in the frozen plan",
            "estimate_confidence": estimate.get("confidence", "UNKNOWN"),
            "steps_table": list(by_step.values()), "report": "results/report/pfvarx_dashboard.html"}


def _clock(iso: str | None) -> str:
    if not iso:
        return "—"
    return dt.datetime.fromisoformat(iso).astimezone().strftime("%a %d %b %Y, %H:%M")


def _duration(minutes: float | None) -> str:
    if minutes is None:
        return "unknown"
    if minutes < 1:
        return "under 1 min"
    whole = int(round(minutes))
    return f"{whole // 60} h {whole % 60:02d} min" if whole >= 60 else f"{whole} min"


def _over(r: dict[str, Any]) -> str:
    if r.get("expected_min") and r.get("running_min") and r["running_min"] > 1.5 * max(r["expected_min"], 1):
        return "<br><small class='warn'>longer than expected</small>"
    return ""


def render(p: dict[str, Any]) -> str:
    e = html.escape
    steps = p.get("steps") or {}
    percent = steps.get("percent")
    if p["finished"]:
        percent = 100
    elif p.get("expected_total_min") and p.get("remaining_min") is not None:
        percent = max(percent or 0, int(100 * (1 - p["remaining_min"] / max(p["expected_total_min"], 0.1))))
    tone = {"COMPLETED": "good", "RUNNING": "info"}.get(p["status"], "warning" if "EXCLUSIONS" in p["status"] or "WITH_FAILURES" in p["status"] else "critical")
    refresh = "" if p["finished"] else '<meta http-equiv="refresh" content="15">'
    if p["finished"]:
        refresh = f'<meta http-equiv="refresh" content="8;url={e(p["report"])}">'
    running = "".join(
        f"<tr><td><strong>{e(r['step'])}</strong><br><small>{e(r['tool_group'])}</small></td><td>{e(r['sample'] or 'cohort')}</td>"
        f"<td>{_clock(r['started'])}</td><td class='num'>{_duration(r['running_min'])}</td>"
        f"<td class='num'>{_duration(r['expected_min'])}{_over(r)}</td></tr>" for r in p["running"]) or \
        "<tr><td colspan='5' class='muted'>No step is running right now.</td></tr>"
    table = "".join(
        f"<tr><td>{e(s['name'])}</td><td class='num'>{s['done']}</td><td class='num'>{s['running']}</td>"
        f"<td class='num{' bad' if s['failed'] else ''}'>{s['failed']}</td><td class='num'>{_duration(s['expected_min'])}</td>"
        f"<td class='num'>{'—' if not s['peak_ram_mb'] else str(round(s['peak_ram_mb'] / 1024, 1)) + ' GB'}</td></tr>"
        for s in p["steps_table"])
    failed = "".join(f"<li>{e(f['rule'])} · {e(f['sample'] or 'cohort')}</li>" for f in p["failed"])
    finish = _clock(p["expected_finish_utc"]) if p["expected_finish_utc"] else "unknown"
    late = (not p["finished"] and p.get("expected_total_min") and p["elapsed_min"] > p["expected_total_min"])
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
{refresh}<title>PfVarX run progress</title><style>
:root{{--navy:#0f2a44;--teal:#16a394;--ink:#14263a;--muted:#5b6f82;--line:#dde5ee;--bg:#eef2f7;--good:#13703a;--good-bg:#e6f4ea;--info:#1c5cab;--info-bg:#e8f1fc;--warning:#8a5300;--warning-bg:#fff3d6;--critical:#b3261e;--critical-bg:#fdecea}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 "Inter","Segoe UI",system-ui,sans-serif}}
header{{background:var(--navy);color:#fff;padding:16px 20px;text-align:center}}header h1{{margin:0;font-size:24px;font-weight:800}}
main{{max-width:1100px;margin:0 auto;padding:18px 16px 40px}}h2{{color:var(--teal);text-transform:uppercase;letter-spacing:.08em;font-size:18px;margin:6px 0 4px}}
.panel{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin:14px 0;box-shadow:0 1px 2px #0f2a440d,0 4px 14px #0f2a440f}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}}.kpi{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:12px 14px}}
.kpi b{{display:block;font-size:20px}}.kpi span{{color:var(--muted);font-size:12px;font-weight:600}}
.bar{{height:14px;border-radius:999px;background:#e3e9f0;overflow:hidden;margin:10px 0 4px}}.bar i{{display:block;height:100%;background:linear-gradient(90deg,#1468a8,var(--teal))}}
.badge{{display:inline-block;padding:3px 12px;border-radius:999px;font-weight:700;font-size:12px}}.good{{color:var(--good);background:var(--good-bg)}}.info{{color:var(--info);background:var(--info-bg)}}
.warning{{color:var(--warning);background:var(--warning-bg)}}.critical{{color:var(--critical);background:var(--critical-bg)}}
table{{width:100%;border-collapse:collapse;font-size:13px}}th,td{{text-align:left;padding:7px 10px;border-bottom:1px solid #edf1f6;vertical-align:top}}
th{{background:#f3f6fa;font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:#3f5568}}td.num{{text-align:right;white-space:nowrap}}td.bad{{color:var(--critical);font-weight:700}}small.warn{{color:var(--warning);font-weight:700}}
.muted,small{{color:var(--muted)}}.note{{color:var(--muted);font-size:12px}}a.btn{{display:inline-block;background:var(--teal);color:#fff;padding:8px 16px;border-radius:8px;font-weight:700;text-decoration:none}}
.wrap{{overflow-x:auto}}</style></head><body>
<header><h1>PfVarX Genomic Analysis Dashboard</h1></header><main>
<h2>PfVarX · Run progress</h2><p class="note">Run {e(p['run_id'])} · updated {_clock(p['updated_utc'])}{'' if p['finished'] else ' · this page refreshes every 15 seconds'}</p>
<div class="panel"><span class="badge {tone}">{e(p['status'])}</span>
<div class="bar"><i style="width:{percent or 0}%"></i></div><div class="note">{percent if percent is not None else '—'} % · {steps.get('done', 0)} of {steps.get('total', '?')} workflow steps done (more steps are added after the QC and eligibility checkpoints)</div>
{f'<p><a class="btn" href="{e(p["report"])}">Open the report dashboard</a> <span class="note">opens by itself in a few seconds</span></p>' if p['finished'] else ''}
{'<p class="warning badge">Taking longer than expected</p>' if late else ''}</div>
<div class="kpis"><div class="kpi"><span>Started</span><b>{_clock(p['started_utc'])}</b></div>
<div class="kpi"><span>Elapsed</span><b>{_duration(p['elapsed_min'])}</b></div>
<div class="kpi"><span>{'Finished' if p['finished'] else 'Expected finish'}</span><b>{_clock(p['finished_utc']) if p['finished'] else finish}</b></div>
<div class="kpi"><span>{'Total time' if p['finished'] else 'Remaining'}</span><b>{_duration(p['elapsed_min']) if p['finished'] else _duration(p['remaining_min'])}</b></div></div>
<p class="note">Estimate: {e(p['estimate_basis'])} ({e(p['estimate_confidence'])}). Times are learned from earlier runs on this computer and scaled by each sample's input size; they are a guide, not a promise.</p>
<div class="panel"><h3>Running now</h3><div class="wrap"><table><tr><th>Step</th><th>Sample</th><th>Started</th><th>Running for</th><th>Expected</th></tr>{running}</table></div></div>
<div class="panel"><h3>Steps</h3><div class="wrap"><table><tr><th>Step</th><th>Done</th><th>Running</th><th>Failed</th><th>Expected time</th><th>Peak RAM</th></tr>{table}</table></div>
{f'<p class="critical badge">Failed jobs</p><ul>{failed}</ul>' if failed else ''}</div>
</main></body></html>"""


def write(run: Path) -> dict[str, Any]:
    p = snapshot(run)
    (run / "progress.json").write_text(json.dumps(p, indent=2, default=str))
    tmp = run / "progress.html.tmp"
    tmp.write_text(render(p))
    tmp.replace(run / "progress.html")
    return p


def open_in_browser(path: Path) -> bool:
    """Open the progress page for a person watching this terminal; never in scripts, tests or headless servers."""
    if os.environ.get("PFVARX_NO_BROWSER") or not sys.stdout.isatty():
        return False
    if sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        return False
    try:
        return webbrowser.open(path.resolve().as_uri())
    except Exception:  # noqa: BLE001 - a missing browser must never stop a run
        return False
