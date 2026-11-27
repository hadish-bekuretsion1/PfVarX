"""Benchmark evaluation of a completed run against declared truth (benchmarking workspace).

* Monoclonal truth samples: SNV and indel sensitivity, precision, F1 and false positives
  per callable Mb within confident regions, after normalisation (bcftools norm) of both
  call and truth sets; exact (chrom, pos, ref, alt) matching after normalisation. RTG
  vcfeval availability is reported; its cross-check is run manually.
* Clustered uncertainty: wild cluster bootstrap over samples (isolates) of the
  linearised ratio estimator, Rademacher weights at >= 10 clusters and Webb six-point
  weights below, all weight vectors enumerated when there are <= 9999 of them.
* Mixtures: WSAF calibration (slope, intercept, RMSE, interval coverage) against
  expected fractions, detection by expected-frequency bin, technical K versus true K.
* Held-out samples are scored only with --holdout, once; the attempt is logged.
Scores are counts on synthetic or declared truth only; they are not field accuracy.
"""

from __future__ import annotations

import gzip
import itertools
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from pfvarx.benchmarking.partitions import load as load_partition
from pfvarx.configuration.loader import frozen_param
from pfvarx.genomic_evidence import in_intervals, read_bed
from pfvarx.io import read_json, write_json, write_tsv

WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5), math.sqrt(0.5), 1.0, math.sqrt(1.5)])


def _vcf_records(path: Path, pass_only: bool) -> set[tuple[str, int, str, str]]:
    out = set()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as h:  # type: ignore[operator]
        for line in h:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if pass_only and c[6] != "PASS":
                continue
            for alt in c[4].split(","):
                out.add((c[0], int(c[1]), c[3].upper(), alt.upper()))
    return out


def _norm(vcf: Path, ref: Path, out: Path) -> Path:
    subprocess.run(["bcftools", "norm", "-f", str(ref), "-m", "-any", "-c", "w", "-Oz", "-o", str(out), str(vcf)],
                   check=True, capture_output=True)
    return out


def wild_cluster_interval(num: np.ndarray, den: np.ndarray, level: float, seed: int, max_vectors: int,
                          rademacher_min_clusters: int) -> dict[str, Any]:
    """Interval for R = sum(num)/sum(den) by wild cluster bootstrap of the linearised estimator.

    Settings are layer6.benchmark_interval_level, benchmark_bootstrap_seed,
    benchmark_bootstrap_max_vectors and benchmark_rademacher_min_clusters."""
    G = len(num)
    if G == 0 or den.sum() == 0:
        return {"estimate": None, "lo": None, "hi": None, "clusters": G}
    R = num.sum() / den.sum()
    u = (num - R * den) / den.sum()  # influence of each cluster
    rademacher = G >= rademacher_min_clusters
    weights = np.array([-1.0, 1.0]) if rademacher else WEBB
    total = len(weights) ** G
    if total <= max_vectors:
        W = np.array(list(itertools.product(weights, repeat=G)))
        method = "enumerated"
    else:
        rng = np.random.default_rng(seed)
        W = rng.choice(weights, size=(max_vectors, G))
        method = f"{max_vectors} random draws"
    stats_ = R + W @ u
    a = (1 - level) / 2
    # design effect: clustered variance / binomial variance
    var_cl = float(G / max(G - 1, 1) * np.sum(u ** 2))
    var_bin = R * (1 - R) / den.sum()
    return {"estimate": float(R), "lo": float(max(0.0, np.quantile(stats_, a))),
            "hi": float(min(1.0, np.quantile(stats_, 1 - a))), "clusters": G,
            "weights": "Rademacher" if rademacher else "Webb six-point", "vectors": method,
            "design_effect": (var_cl / var_bin) if var_bin > 0 else None}


def evaluate(root: Path, run: Path, out: Path, holdout: bool = False) -> dict[str, Any]:
    truth = read_json(root / "inputs" / "truth" / "truth.json")
    part = load_partition(root)
    use = part["held_out"] if holdout else part["benchmarking"]
    log = out / "holdout_attempts.jsonl"
    if holdout:
        if log.is_file() and log.read_text().strip():
            raise ValueError("held-out evaluation already performed once; record a new experiment instead")
        out.mkdir(parents=True, exist_ok=True)
        log.write_text(json.dumps({"run_id": run.name, "partition_sha256": part["sha256"]}) + "\n")
    res = run / "results"
    ref = root / read_json(run / "frozen_plan.json")["reference"]["parasite_fasta"]
    # Settings come from the evaluated run's frozen parameters (layer6.benchmark_*); the interval
    # coverage is judged against the level the caller's WSAF intervals were computed at.
    params = {k.split(".", 1)[1]: v["value"] for k, v in read_json(run / "frozen_config.json")["parameters"].items()}
    level = float(frozen_param(params, "benchmark_interval_level"))
    seed = int(frozen_param(params, "benchmark_bootstrap_seed"))
    max_vectors = int(frozen_param(params, "benchmark_bootstrap_max_vectors"))
    rademacher_min = int(frozen_param(params, "benchmark_rademacher_min_clusters"))
    clonal = float(frozen_param(params, "benchmark_clonal_min_expected_wsaf"))
    wsaf_level = float(frozen_param(params, "wsaf_interval_level"))
    confident = read_bed(root / "inputs" / "truth" / "confident_regions.bed")
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for sid in use:
        info = truth["samples"][sid]
        if info["clonality"] != "MONOCLONAL":
            continue
        calls_vcf = res / "variants" / "samples" / f"{sid}.vcf.gz"
        tvcf = root / "inputs" / "truth" / f"{sid}.truth.vcf"
        if not calls_vcf.is_file():
            rows.append({"sample_id": sid, "status": "NO_CALLS"})
            continue
        n_calls = _norm(calls_vcf, ref, out / f"{sid}.calls.norm.vcf.gz")
        n_truth = _norm(tvcf, ref, out / f"{sid}.truth.norm.vcf.gz")
        called = {r for r in _vcf_records(n_calls, True) if in_intervals(confident.get(r[0], []), r[1] - 1)}
        tru = {r for r in _vcf_records(n_truth, False) if in_intervals(confident.get(r[0], []), r[1] - 1)}
        callable_iv = read_bed(res / "variants" / "samples" / f"{sid}.callable.bed")
        callable_mb = sum(e - s for iv in callable_iv.values() for s, e in iv) / 1e6
        for vt in ("SNV", "INDEL"):
            def is_t(r: tuple[str, int, str, str]) -> bool:
                return (len(r[2]) == 1 and len(r[3]) == 1) == (vt == "SNV")
            c, t = {r for r in called if is_t(r)}, {r for r in tru if is_t(r)}
            tp, fp, fn = len(c & t), len(c - t), len(t - c)
            rows.append({"sample_id": sid, "variant_type": vt, "status": "SCORED", "tp": tp, "fp": fp, "fn": fn,
                         "sensitivity": tp / (tp + fn) if tp + fn else None,
                         "precision": tp / (tp + fp) if tp + fp else None,
                         "fp_per_callable_mb": fp / callable_mb if callable_mb else None})
    summary = {}
    for vt in ("SNV", "INDEL"):
        sc = [r for r in rows if r.get("variant_type") == vt]
        true_positives = np.array([r["tp"] for r in sc], float)
        false_negatives = np.array([r["fn"] for r in sc], float)
        false_positives = np.array([r["fp"] for r in sc], float)
        sens = wild_cluster_interval(true_positives, true_positives + false_negatives, level, seed, max_vectors,
                                     rademacher_min)
        prec = wild_cluster_interval(true_positives, true_positives + false_positives, level, seed, max_vectors,
                                     rademacher_min)
        denominator = float((2 * true_positives + false_positives + false_negatives).sum())
        f1 = float(2 * true_positives.sum() / denominator) if denominator else None
        summary[vt] = {"sensitivity": sens, "precision": prec, "f1": f1, "samples": len(sc)}
    vcfeval = None
    if shutil.which("rtg"):
        vcfeval = "RTG vcfeval available; run 'rtg format' on the reference and 'rtg vcfeval' for cross-check (not automated)"
    # mixtures
    mix_rows, k_rows = [], []
    clones = truth["clones"]
    for sid in use:
        info = truth["samples"][sid]
        qc_p = res / "qc" / "samples" / sid / f"{sid}.caller_qc.json"
        if qc_p.is_file():
            q = read_json(qc_p)["mixture"]
            k_rows.append({"sample_id": sid, "true_k": info["true_k"], "technical_k": q.get("k_selected"),
                           "status": q["status"], "weights": q.get("fit", {}).get("weights")})
        if info["clonality"] == "MONOCLONAL":
            continue
        expect: dict[tuple[str, int, str, str], float] = {}
        for comp in info["components"]:
            for v in clones[comp["clone"]]:
                if len(v["ref"]) == 1 and len(v["alt"]) == 1:
                    key = (v["chrom"], v["pos"], v["ref"], v["alt"])
                    expect[key] = expect.get(key, 0.0) + comp["proportion"]
        ev = res / "variants" / "samples" / f"{sid}.evidence.tsv"
        obs = {}
        if ev.is_file():
            lines = ev.read_text().splitlines()
            col = {name: i for i, name in enumerate(lines[0].split("\t"))} if lines else {}
            for line in lines[1:]:
                fields = line.split("\t")
                if fields[col["kind"]] == "SNV":
                    obs[(fields[col["chrom"]], int(fields[col["pos"]]), fields[col["ref"]], fields[col["alt"]])] = (
                        float(fields[col["wsaf"]]), float(fields[col["wsaf_lo"]]), float(fields[col["wsaf_hi"]]),
                        fields[col["status"]])
        for key, e in expect.items():
            if e >= clonal:
                continue
            o = obs.get(key)
            mix_rows.append({"sample_id": sid, "chrom": key[0], "pos": key[1], "expected": round(e, 4),
                             "observed": None if o is None else o[0], "lo": None if o is None else o[1],
                             "hi": None if o is None else o[2], "status": "NOT_DETECTED" if o is None else o[3],
                             "covered": None if o is None else (o[1] <= e <= o[2])})
    det = [r for r in mix_rows if r["observed"] is not None]
    calib = {}
    if len(det) >= int(frozen_param(params, "benchmark_calibration_min_sites")):
        x = np.array([r["expected"] for r in det])
        y = np.array([r["observed"] for r in det])
        slope, intercept = np.polyfit(x, y, 1)
        calib = {"sites": len(det), "slope": float(slope), "intercept": float(intercept),
                 "rmse": float(np.sqrt(np.mean((y - x) ** 2))),
                 "interval_coverage": float(np.mean([bool(r["covered"]) for r in det])), "nominal": wsaf_level}
    bins = {}
    for lo, hi in ((float(a), float(b)) for a, b in frozen_param(params, "benchmark_wsaf_bins")):
        sel = [r for r in mix_rows if lo <= r["expected"] < hi]
        if sel:
            bins[f"{lo}-{hi}"] = {"sites": len(sel), "detected_pass": sum(1 for r in sel if r["status"] == "PASS"),
                                  "detection_rate": sum(1 for r in sel if r["status"] == "PASS") / len(sel)}
    write_tsv(out / "per_sample_accuracy.tsv", rows, ["sample_id", "variant_type", "status", "tp", "fp", "fn",
                                                       "sensitivity", "precision", "fp_per_callable_mb"])
    write_tsv(out / "mixture_sites.tsv", mix_rows, ["sample_id", "chrom", "pos", "expected", "observed", "lo", "hi",
                                                    "status", "covered"])
    write_tsv(out / "k_selection.tsv", k_rows, ["sample_id", "true_k", "technical_k", "status", "weights"])
    report = {"schema": "pfvarx.benchmark_report/1", "run_id": run.name, "partition": "held_out" if holdout else "benchmarking",
              "partition_sha256": part["sha256"], "accuracy": summary, "vcfeval": vcfeval,
              "wsaf_calibration": calib, "detection_by_expected_fraction": bins,
              "k_selection": {"correct": sum(1 for r in k_rows if r["technical_k"] == r["true_k"]), "total": len(k_rows)},
              "truth_note": "scores describe agreement with the declared truth set only (synthetic fixtures: software testing, not field accuracy)",
              "matching": "exact after bcftools norm (left-aligned, decomposed); within confident regions"}
    write_json(out / "benchmark_report.json", report)
    return report
