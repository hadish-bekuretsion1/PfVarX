"""Layer 3 cohort QC (gate G5): how complete each sample is in the cohort.

The cohort evidence step classifies every sample at every cohort site as ALT_DOMINANT, MIXED,
REF_DOMINANT, MISSING_LOW_DEPTH or UNCALLABLE (``cohort/missingness.tsv``). This step applies
the frozen QC policy's G5 gates to those counts, the same way Layer 1 applies G1-G3:

* ``G5.cohort_missingness``: the share of cohort sites where the sample is MISSING_LOW_DEPTH or
  UNCALLABLE. A FAIL excludes the sample from the downstream analyses (FO004); a WARNING lets
  it continue into the stages the gate lists.
* ``G5.cohort_mixed_fraction``: the share of cohort sites with MIXED evidence (descriptive).

It also summarises site completeness across the cohort (how many sites are called in at least
``site_call_rate_min`` of the samples), which is reported, not gated.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Any

from pfvarx.constants import Design
from pfvarx.quality.policy import evaluate_sample, load_policy, stage_permitted

COHORT_GATES = {"G5"}
CALLED = {"ALT_DOMINANT", "MIXED", "REF_DOMINANT"}
SITE_CALL_RATE_MIN = 0.8
DOWNSTREAM_STAGES = ["L4", "DA001", "DA002", "DA003", "DA004", "DA005", "DA006", "DA007", "DA008"]


def sample_metrics(missingness: Path) -> dict[str, dict[str, Any]]:
    """Per-sample G5 metrics from ``cohort/missingness.tsv``."""
    out: dict[str, dict[str, Any]] = {}
    with missingness.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts = {k: int(row[k]) for k in ("ALT_DOMINANT", "MIXED", "REF_DOMINANT", "MISSING_LOW_DEPTH", "UNCALLABLE")}
            sites = sum(counts.values())
            out[row["sample_id"]] = {"cohort_sites": sites, **counts,
                                     "cohort_missing_fraction": float(row["missing_fraction"]),
                                     "cohort_mixed_fraction": counts["MIXED"] / sites if sites else None}
    return out


def site_summary(matrix: Path, n_samples: int, call_rate_min: float = SITE_CALL_RATE_MIN) -> dict[str, Any]:
    """How many cohort sites are called (ALT/MIXED/REF) in at least *call_rate_min* of the samples."""
    called: Counter[str] = Counter()
    sites: set[str] = set()
    with gzip.open(matrix, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sites.add(row["site_id"])
            if row["state"] in CALLED:
                called[row["site_id"]] += 1
    n = len(sites)
    good = sum(1 for s in sites if n_samples and called[s] / n_samples >= call_rate_min)
    return {"sites": n, "samples": n_samples, "site_call_rate_min": call_rate_min,
            "sites_meeting_call_rate": good, "fraction_meeting_call_rate": good / n if n else None}


def evaluate(root: Path, frozen: dict[str, Any], plan: dict[str, Any], missingness: Path,
             matrix: Path | None) -> dict[str, Any]:
    design = Design.parse(plan.get("design") or frozen["config"]["assay_design"]).value
    policy = load_policy(root, frozen["quality"]["policy_id"])
    metrics = sample_metrics(missingness)
    samples: dict[str, Any] = {}
    for sid, values in sorted(metrics.items()):
        qc = evaluate_sample(policy, values, design, stages=COHORT_GATES)
        qc["permitted"] = {st: stage_permitted(qc, st) for st in DOWNSTREAM_STAGES}
        qc["metrics"] = values
        samples[sid] = qc
    return {"schema": "pfvarx.cohort_qc/1", "policy_id": policy["policy_id"], "policy_version": policy["version"],
            "gates": sorted(COHORT_GATES), "samples": samples,
            "sites": site_summary(matrix, len(metrics)) if matrix is not None and matrix.is_file() else None}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--frozen-config", required=True)
    ap.add_argument("--missingness", required=True)
    ap.add_argument("--matrix")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    frozen_path = Path(a.frozen_config)
    frozen = json.loads(frozen_path.read_text())
    plan = json.loads((frozen_path.parent / "frozen_plan.json").read_text())
    result = evaluate(Path(a.root), frozen, plan, Path(a.missingness), Path(a.matrix) if a.matrix else None)
    Path(a.out).write_text(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
