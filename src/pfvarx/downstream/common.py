"""Shared Layer 5 data adapters and the terminal status manifest (Appendix F.5)."""

from __future__ import annotations

import datetime as dt
import gzip
import json
import math
from array import array
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pfvarx.constants import CallStatus, Design
from pfvarx.io import read_tsv, write_json

# The vocabularies the evidence matrix uses, stored as one byte per cell.
STATE_VOCABULARY = ("ALT_DOMINANT", "MIXED", "REF_DOMINANT", "MISSING_LOW_DEPTH", "UNCALLABLE")
STATUS_VOCABULARY = ("", *(s.value for s in CallStatus))
_STATE_CODES = {s: i for i, s in enumerate(STATE_VOCABULARY)}
_STATUS_CODES = {s: i for i, s in enumerate(STATUS_VOCABULARY)}

MANIFEST_SCHEMA = "pfvarx.analysis_status/1"
PANEL_CAVEAT = ("Targeted panel: only the panel's target regions are analysed, and sites within one "
                "amplicon or target region are linked, so the result describes the panel loci, not the genome.")


class Cell:
    """One sample's evidence at one cohort site. Reads like the mapping it replaces
    (``cell["wsaf"]``), so every consumer and hand-built test fixture keeps working."""

    __slots__ = ("state", "ref", "alt", "depth", "wsaf", "caller_status")

    def __init__(self, state: str, ref: int, alt: int, depth: int, wsaf: float | None,
                 caller_status: str) -> None:
        self.state, self.ref, self.alt = state, ref, alt
        self.depth, self.wsaf, self.caller_status = depth, wsaf, caller_status

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, default)

    def __repr__(self) -> str:
        return (f"Cell(state={self.state!r}, ref={self.ref}, alt={self.alt}, depth={self.depth}, "
                f"wsaf={self.wsaf!r}, caller_status={self.caller_status!r})")


class SiteRow:
    """Every sample's evidence at one site, held as arrays indexed by the cohort's sample order.

    A cohort of a few hundred samples over tens of thousands of sites is tens of millions of
    cells, and one dict per cell costs about 474 bytes of them (measured). Columns cost about 30
    and are addressed directly by the sample's position, so lookup stays constant-time. Counts
    are integers and ``wsaf`` is a float64, exactly the values the file holds: nothing is
    rounded, so no threshold comparison can move.
    """

    __slots__ = ("_index", "_state", "_status", "_ref", "_alt", "_depth", "_wsaf")
    ABSENT = 255  # no record for this sample at this site: .get() gives None, as a dict did

    def __init__(self, index: dict[str, int]) -> None:
        n = len(index)
        self._index = index
        self._state = array("B", [self.ABSENT]) * n
        self._status = array("B", [0]) * n
        self._ref = array("i", [0]) * n
        self._alt = array("i", [0]) * n
        self._depth = array("i", [0]) * n
        self._wsaf = array("d", [-1.0]) * n

    def set(self, position: int, state: int, status: int, ref: int, alt: int, depth: int,
            wsaf: float | None) -> None:
        self._state[position] = state
        self._status[position] = status
        self._ref[position] = ref
        self._alt[position] = alt
        self._depth[position] = depth
        self._wsaf[position] = -1.0 if wsaf is None else wsaf

    def get(self, sample: str, default: Any = None) -> Cell | Any:
        position = self._index.get(sample)
        if position is None or self._state[position] == self.ABSENT:
            return default
        wsaf = self._wsaf[position]
        return Cell(STATE_VOCABULARY[self._state[position]], self._ref[position], self._alt[position],
                    self._depth[position], None if wsaf < 0 else wsaf,
                    STATUS_VOCABULARY[self._status[position]])

    def __getitem__(self, sample: str) -> Cell:
        cell = self.get(sample)
        if cell is None:
            raise KeyError(sample)
        return cell

    def __contains__(self, sample: str) -> bool:
        return self.get(sample) is not None

    def __iter__(self) -> Iterator[str]:
        return (s for s, position in self._index.items() if self._state[position] != self.ABSENT)

    def __len__(self) -> int:
        return sum(1 for _ in self)


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
    matrix: dict[str, Any] = field(default_factory=dict)  # site -> SiteRow (or, in fixtures, sample -> cell)

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
        load_matrix(ctx, run / "results" / "cohort" / "evidence_matrix.tsv.gz")
    return ctx


def load_matrix(ctx: Context, path: Path) -> dict[str, int]:
    """Read the cohort evidence matrix into *ctx*; returns what was skipped, for diagnostics.

    The cohort's sample order is known before the file is read (``cohort_samples.txt``), so each
    site's columns are allocated once and addressed by position: one pass, constant-time lookup
    and no intermediate copy of the matrix.
    """
    index = {sample: position for position, sample in enumerate(ctx.samples)}
    skipped = {"unknown_sample": 0, "unknown_state": 0}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        columns = {name: i for i, name in enumerate(header)}
        site_at, sample_at = columns["site_id"], columns["sample_id"]
        state_at, status_at = columns["state"], columns["caller_status"]
        ref_at, alt_at, depth_at, wsaf_at = (columns["ref_count"], columns["alt_count"],
                                             columns["depth"], columns["wsaf"])
        for line in handle:
            field_values = line.rstrip("\n").split("\t")
            position = index.get(field_values[sample_at])
            if position is None:  # not a cohort member: the matrix is written from the cohort
                skipped["unknown_sample"] += 1
                continue
            state = _STATE_CODES.get(field_values[state_at])
            if state is None:
                skipped["unknown_state"] += 1
                continue
            site = field_values[site_at]
            row = ctx.matrix.get(site)
            if row is None:
                row = ctx.matrix[site] = SiteRow(index)
            wsaf = field_values[wsaf_at]
            row.set(position, state, _STATUS_CODES.get(field_values[status_at], 0),
                    int(field_values[ref_at]), int(field_values[alt_at]), int(field_values[depth_at]),
                    float(wsaf) if wsaf else None)
    return skipped


def cohort_sites(run: Path) -> int:
    """Cohort sites, counted from sites.tsv without reading the evidence matrix."""
    path = run / "results" / "cohort" / "sites.tsv"
    if not path.is_file():
        return 0
    return max(0, sum(1 for _ in path.open()) - 1)


def matrix_memory_mb(samples: int, sites: int, bytes_per_cell: int) -> int:
    """Memory the cohort evidence matrix needs in an analysis that reads it.

    Every cohort site holds a column for every sample (the evidence is site-forced), so the cost
    grows with the product. ``bytes_per_cell`` is measured, not guessed: see
    ``planning.cohort_matrix_bytes_per_cell`` in config/registry/job_resources.yaml.
    """
    return int(samples * sites * bytes_per_cell / 1e6)


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
