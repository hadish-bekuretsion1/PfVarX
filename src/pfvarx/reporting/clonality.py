"""Per-sample clonality report: monoclonal, or polyclonal with the number of clones and each
clone's share of the parasite DNA.

Sources, per sample:
  * the PfVarX-Call mixture model (``qc/samples/S/S.caller_qc.json``): the number of DNA
    components it selects and their weights. This gives the count and the percentages.
  * DA001 (``downstream/DA001/clonality.tsv``), when MOIRE ran (cohorts of at least
    ``population_frequency_min_samples``): the MOIRE monoclonal/polyclonal call. In smaller
    cohorts MOIRE is not run, and DA001's call comes from the same mixture test, so it is
    not counted again. The MOIRE COI count is never used here.

Written to ``results/report/clonality_report.txt`` (one readable line per sample) and
``clonality_report.tsv``, and shown in the dashboard overview.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pfvarx.io import read_tsv, write_tsv

COLUMNS = ["sample_id", "clonality", "number_of_clones", "clone_proportions", "basis", "note"]


def _mixture(results: Path, sample: str) -> dict[str, Any] | None:
    path = results / "qc" / "samples" / sample / f"{sample}.caller_qc.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text()).get("mixture")


def _da001(results: Path) -> dict[str, str]:
    path = results / "downstream" / "DA001" / "clonality.tsv"
    if not path.is_file():
        return {}
    # A call DA001 took from the caller's mixture test (small cohorts) is the same evidence
    # as the mixture model used here, so only MOIRE calls count as DA001's own.
    return {row["sample_id"]: row.get("clonality", "") for row in read_tsv(path)
            if row.get("basis", "MOIRE") in ("MOIRE", "")}


def _components(mix: dict[str, Any]) -> tuple[int, list[float], bool, float]:
    """Selected number of components, their weights (largest first), whether the largest
    tested number was reached, and the smallest weight the model treats as identifiable."""
    k = int(mix["k_selected"])
    fit = mix.get("fits", {}).get(str(k)) or mix.get("fit") or {}
    weights = sorted((float(w) for w in fit.get("weights", [])), reverse=True)
    tests = mix.get("tests") or []
    at_limit = bool(tests) and "K_MAX" in str(tests[-1].get("note", "")) and k >= int(tests[-1].get("k", k + 1))
    return k, weights, at_limit, float(mix.get("min_component_weight", 0.02))


def classify(sample: str, mix: dict[str, Any] | None, da001_call: str | None) -> dict[str, str]:
    """One report row for one sample."""
    row = {"sample_id": sample, "clonality": "NOT_DETERMINED", "number_of_clones": "", "clone_proportions": "",
           "basis": "", "note": ""}
    estimated = mix is not None and mix.get("status") == "ESTIMATED"
    if not estimated:
        reason = "no caller mixture estimate" if mix is None else f"caller mixture {mix.get('status', 'unavailable')}"
        if da001_call in ("MONOCLONAL", "POLYCLONAL"):
            row.update(clonality=da001_call, number_of_clones="1" if da001_call == "MONOCLONAL" else "not resolved",
                       basis="DA001 (MOIRE)", note=reason + "; clone proportions unavailable")
        else:
            row["note"] = reason + ("; DA001 " + da001_call if da001_call else "")
        return row
    assert mix is not None
    k, weights, at_limit, floor = _components(mix)
    shares = "; ".join(f"Clone {i}: {w * 100:.1f}%" for i, w in enumerate(weights, start=1))
    row["basis"] = "PfVarX-Call mixture" + (" + DA001 (MOIRE)" if da001_call else "")
    if k == 1:
        if da001_call == "POLYCLONAL":
            row.update(clonality="POLYCLONAL", number_of_clones="not resolved",
                       note="DA001 detects more than one clone, but no minor clone reached "
                            f"{floor:.0%} of the DNA in the caller's mixture model")
        else:
            row.update(clonality="MONOCLONAL", number_of_clones="1", clone_proportions="Clone 1: 100.0%")
        return row
    if da001_call == "MONOCLONAL":
        row.update(clonality="INDETERMINATE", number_of_clones="", clone_proportions=shares,
                   note=f"the caller's mixture model finds {k} components but DA001 calls the sample monoclonal")
        return row
    row.update(clonality="POLYCLONAL", number_of_clones=f"at least {k}" if at_limit else str(k),
               clone_proportions=shares,
               note=(f"the largest number of clones tested ({k}) was selected; more may be present. " if at_limit else "")
               + f"Further clones below {floor:.0%} of the DNA cannot be resolved")
    return row


def sentence(row: dict[str, str]) -> str:
    sample, call = row["sample_id"], row["clonality"]
    if call == "MONOCLONAL":
        text = f"Sample {sample} is monoclonal"
    elif call == "POLYCLONAL":
        n = row["number_of_clones"]
        text = f"Sample {sample} is polyclonal"
        text += f" with {n} clones" if n and n != "not resolved" else " (number of clones not resolved)"
        if row["clone_proportions"]:
            text += ": " + row["clone_proportions"].replace("; ", ", ")
    elif call == "INDETERMINATE":
        text = f"Sample {sample} is indeterminate"
    else:
        text = f"Sample {sample}: clonality not determined"
    return text + (f" ({row['note']})" if row["note"] and call != "MONOCLONAL" else "")


def write_report(results: Path, samples: list[str], run_id: str) -> list[dict[str, str]]:
    """Build, write and return the per-sample clonality report."""
    da001 = _da001(results)
    rows = [classify(s, _mixture(results, s), da001.get(s)) for s in samples]
    out = results / "report"
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "clonality_report.tsv", rows, COLUMNS)
    lines = [f"PfVarX clonality report, run {run_id}", ""]
    lines += [sentence(r) for r in rows] or ["No samples reached variant calling."]
    lines += ["", "Clone percentages are shares of the parasite DNA estimated by the PfVarX-Call mixture model",
              "(uncalibrated until the caller is benchmarked). They approximate, but are not exactly, the share",
              "of parasites. MOIRE (DA001) is used only in cohorts large enough to estimate allele frequencies; its clone count is not used here."]
    (out / "clonality_report.txt").write_text("\n".join(lines) + "\n")
    return rows
