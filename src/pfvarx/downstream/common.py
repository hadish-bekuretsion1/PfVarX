"""Shared Layer 5 data adapters and the terminal status manifest (Appendix F.5)."""

from __future__ import annotations

import datetime as dt
import gzip
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pfvarx.constants import Design
from pfvarx.io import read_tsv, write_json

MANIFEST_SCHEMA = "pfvarx.analysis_status/1"
PANEL_CAVEAT = ("Targeted panel: only the panel's target regions are analysed, and sites within one "
                "amplicon or target region are linked, so the result describes the panel loci, not the genome.")


@dataclass
class Context:
    root: Path
    run: Path
    eff: dict[str, Any]
    plan: dict[str, Any]
    params: dict[str, Any]
    samples: list[str]
    metadata: dict[str, dict[str, str]]
    sites: list[dict[str, str]] = field(default_factory=list)
    matrix: dict[str, dict[str, dict[str, Any]]] = field(default_factory=dict)  # site -> sample -> cell

    @property
    def results(self) -> Path:
        return self.run / "results"

    def out(self, da: str) -> Path:
        d = self.results / "downstream" / da
        d.mkdir(parents=True, exist_ok=True)
        return d

    def gate(self, sample: str) -> dict[str, Any]:
        return json.loads((self.results / "qc" / "samples" / sample / "qc_decision.json").read_text())

    @property
    def design(self) -> Design:
        """The run's assay design (runs frozen before config schema 4 recorded assay.design)."""
        cfg = self.eff["config"]
        return Design.parse(self.plan.get("design") or cfg.get("assay_design") or cfg["assay"]["design"])


def load_context(root: Path, run: Path, with_matrix: bool = True) -> Context:
    eff = json.loads((run / "frozen_config.json").read_text())
    plan = json.loads((run / "frozen_plan.json").read_text())
    params = {k.split(".", 1)[1]: v["value"] for k, v in eff["parameters"].items()}
    cohort = run / "results" / "cohort" / "cohort_samples.txt"
    samples = cohort.read_text().split() if cohort.is_file() else []
    meta: dict[str, dict[str, str]] = {}
    for row in read_tsv(root / eff["config"]["samples_tsv"]):
        meta.setdefault(row["sample_id"], row)
    ctx = Context(root=root, run=run, eff=eff, plan=plan, params=params, samples=samples, metadata=meta)
    if with_matrix and (run / "results" / "cohort" / "sites.tsv").is_file():
        ctx.sites = read_tsv(run / "results" / "cohort" / "sites.tsv")
        with gzip.open(run / "results" / "cohort" / "evidence_matrix.tsv.gz", "rt") as h:
            header = h.readline().rstrip("\n").split("\t")
            for line in h:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                cell = {"state": row["state"], "ref": int(row["ref_count"]), "alt": int(row["alt_count"]),
                        "depth": int(row["depth"]), "wsaf": float(row["wsaf"]) if row["wsaf"] else None,
                        "caller_status": row["caller_status"]}
                ctx.matrix.setdefault(row["site_id"], {})[row["sample_id"]] = cell
    return ctx


def decisions(ctx: Context) -> dict[str, Any]:
    return json.loads((ctx.results / "eligibility" / "cohort_decisions.json").read_text())


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float | None, float | None]:
    """Wilson score interval (simple random sampling only)."""
    if n == 0:
        return None, None
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def write_manifest(ctx: Context, da: str, terminal: str, *, reason_codes: list[str], requested: list[str],
                   eligible: list[str], excluded: dict[str, list[str]], analysed: list[str],
                   outputs: list[Path], caveats: list[str], method: dict[str, Any],
                   results_summary: dict[str, Any] | None = None, scale: str | None = None) -> Path:
    from pfvarx.io import sha256

    rel = [p for p in outputs if p.is_file()]
    # Samples blocked before the cohort (input validation, QC gate, failed jobs) stay in the
    # denominator trail: they are never silently removed from membership reports.
    registered = sorted({s["sample_id"] for s in ctx.plan["samples"]} | set(ctx.plan.get("samples_failed_validation", [])))
    upstream = [s for s in registered if s not in requested]
    if upstream:
        excluded = {**excluded, "BLOCKED_UPSTREAM": sorted(set(excluded.get("BLOCKED_UPSTREAM", [])) | set(upstream))}
        requested = registered
    manifest = {
        "schema": MANIFEST_SCHEMA, "analysis_id": da, "terminal_state": terminal, "scale": scale,
        "reason_codes": reason_codes,
        "decision_evidence": {
            "sample_qc_directory": "results/qc/samples/",
            "input_validation": "metadata/validation_report.json",
            "cohort_decisions": "results/eligibility/cohort_decisions.json",
            "method_log": f"logs/downstream/{da}.log",
        },
        "recommended_action": (
            "Inspect the method log and dependencies, then resume after correcting the execution error."
            if terminal == "FAILED_EXECUTION" else
            "Review the reason codes, sample QC and method prerequisites before starting a new run."
            if terminal in {"SKIPPED_INELIGIBLE", "BLOCKED_MISSING_PREREQUISITE"} else None
        ),
        "denominators": {"requested": len(requested), "eligible": len(eligible),
                         "excluded": {k: len(v) for k, v in excluded.items()}, "analysed": len(analysed)},
        "samples": {"requested": requested, "eligible": eligible, "excluded": excluded, "analysed": analysed},
        "outputs": [{"path": str(p.relative_to(ctx.run)), "sha256": sha256(p)} for p in rel],
        "caveats": caveats, "method": method, "results_summary": results_summary or {},
        "run_id": ctx.run.name, "project_id": ctx.eff["config"]["project_id"],
        "decided_at": dt.datetime.now(dt.UTC).isoformat(),
    }
    path = ctx.out(da) / "status_manifest.json"
    write_json(path, manifest)
    return path
