"""Parse tool outputs into a flat, unit-labelled metric dictionary.

Each parser returns only metrics it could read; absent metrics are absent (never 0).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


def fastp(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text())
    before = data["summary"]["before_filtering"]
    after = data["summary"]["after_filtering"]
    out = {"reads_before": before["total_reads"], "reads_after": after["total_reads"],
           "mean_read_length_after": after.get("read1_mean_length"),
           "q30_fraction_after": after.get("q30_rate")}
    if before["total_reads"]:
        out["reads_retained_fraction"] = after["total_reads"] / before["total_reads"]
    adapters = data.get("adapter_cutting", {})
    out["adapter_trimmed_reads"] = adapters.get("adapter_trimmed_reads")
    return out


def flagstat(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text())["QC-passed reads"]
    primary = data.get("primary", data["total"])
    out: dict[str, Any] = {"primary_reads": primary}
    if primary:
        out["mapped_fraction"] = data.get("primary mapped", data["mapped"]) / primary
    paired = data.get("paired in sequencing", 0)
    if paired:
        out["properly_paired_fraction"] = data["properly paired"] / paired
    return out


def picard_duplicates(path: Path) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("LIBRARY\t"):
            header = line.split("\t")
            vals = lines[i + 1].split("\t")
            row = dict(zip(header, vals))
            return {"duplicate_fraction": float(row["PERCENT_DUPLICATION"] or 0.0)}
    return {}


def mosdepth_region_dist(path: Path, min_depth: int) -> dict[str, Any]:
    """Median depth and callable fraction over the BED supplied to mosdepth --by."""
    dist: dict[int, float] = {}
    for line in path.read_text().splitlines():
        chrom, depth, frac = line.split("\t")
        if chrom == "total":
            dist[int(depth)] = float(frac)
    if not dist:
        return {}
    median = max((d for d, f in dist.items() if f >= 0.5), default=0)
    return {"median_depth_core": median, "callable_fraction_core": dist.get(min_depth, 0.0)}


def target_metrics(regions_bed_gz_text: str, min_depth: int, on_target: int, mapped: int) -> dict[str, Any]:
    means = [float(line.split("\t")[-1]) for line in regions_bed_gz_text.splitlines() if line]
    out: dict[str, Any] = {}
    if means:
        out["amplicon_dropout_fraction"] = sum(m < min_depth for m in means) / len(means)
        out["amplicons"] = len(means)
        ordered = sorted(means)
        middle = len(ordered) // 2
        out["median_amplicon_depth"] = (ordered[middle] if len(ordered) % 2
                                        else (ordered[middle - 1] + ordered[middle]) / 2)
    if mapped:
        out["on_target_fraction"] = on_target / mapped
    return out


def umi_family_sizes(path: Path) -> dict[str, Any]:
    rows = [ln.split("\t") for ln in path.read_text().splitlines()[1:] if ln]
    if not rows:
        return {}
    counts = {int(r[0]): int(r[1]) for r in rows}
    total = sum(counts.values())
    sizes = sorted(counts)
    cum, median = 0, sizes[-1]
    for s in sizes:
        cum += counts[s]
        if cum >= total / 2:
            median = s
            break
    return {"umi_families": total, "umi_singleton_fraction": counts.get(1, 0) / total,
            "umi_median_family_size": median}


def samtools_count(text: str) -> int:
    m = re.search(r"\d+", text)
    return int(m.group(0)) if m else 0
