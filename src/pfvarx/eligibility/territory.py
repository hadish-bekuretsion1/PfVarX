"""The genomic territory a run can speak about, and each sample's callable evidence within it.

Eligibility for downstream analyses follows from what was sequenced and what is callable,
not from the design label alone:

* the **territory** is the core genome for ``wgs`` and the panel's targeted intervals for
  targeted designs (``panel_id`` + ``panel_version``);
* a panel supplies as many **independent loci** as it has separate target regions (overlapping
  or touching intervals are one locus; sites inside one amplicon are linked); the core genome
  is genome-wide and is not limited this way;
* a **catalogue marker** is in the territory when its whole codon lies inside it, and is
  **callable** for a sample when its whole codon lies inside that sample's callable regions.

The facts are recorded in the eligibility decisions so each verdict can be traced.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

Intervals = dict[str, list[tuple[int, int]]]


def read_intervals(path: Path) -> Intervals:
    """Merged half-open intervals of a BED file, by contig."""
    raw: dict[str, list[tuple[int, int]]] = {}
    if not path.is_file():
        return {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        c = line.split("\t")
        raw.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    merged: Intervals = {}
    for chrom, spans in raw.items():
        out: list[tuple[int, int]] = []
        for s, e in sorted(spans):
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def covers(intervals: Intervals, chrom: str, positions: list[int]) -> bool:
    spans = intervals.get(chrom, [])
    return bool(positions) and all(any(s <= p < e for s, e in spans) for p in positions)


def describe(ctx: Any) -> dict[str, Any]:
    """Territory facts for the run in *ctx* (a downstream Context)."""
    from pfvarx.downstream.resistance import catalog_path, codon_positions, load_catalog

    design = ctx.design
    if design.targeted:
        panel = ctx.plan.get("panel") or {}
        bed = panel.get("targets_bed") or ctx.eff["config"].get("assay", {}).get("targets_bed")
        kind = "PANEL"
    else:
        panel, bed, kind = {}, ctx.eff["config"]["reference"]["core_mask_bed"], "CORE_GENOME"
    intervals = read_intervals(ctx.root / bed) if bed else {}
    facts: dict[str, Any] = {
        "kind": kind, "bed": bed, "panel_id": panel.get("panel_id"), "panel_version": panel.get("panel_version"),
        "regions": sum(len(v) for v in intervals.values()),
        "bp": sum(e - s for v in intervals.values() for s, e in v),
        # Separate panel regions are the independent loci a panel can contribute.
        "independent_loci": sum(len(v) for v in intervals.values()) if kind == "PANEL" else None,
    }
    markers: list[dict[str, str]] = []
    positions: list[list[int]] = []
    path = catalog_path(ctx)
    if path.is_file():
        _, markers = load_catalog(path)
        positions = codon_positions(ctx, markers)
    inside = [m["marker_id"] for m, pos in zip(markers, positions) if covers(intervals, m.get("chrom", ""), pos)]
    facts.update(catalogue_markers=len(markers), catalogue_markers_in_territory=len(inside),
                 markers_in_territory=inside)
    facts["_markers"] = [(m.get("chrom", ""), pos) for m, pos in zip(markers, positions)
                         if m["marker_id"] in set(inside)]
    return facts


def callable_markers(ctx: Any, sample: str, territory: dict[str, Any]) -> int:
    """Territory markers whose whole codon is callable in *sample*."""
    callable_ = read_intervals(ctx.results / "variants" / "samples" / f"{sample}.callable.bed")
    return sum(covers(callable_, chrom, pos) for chrom, pos in territory.get("_markers", []))
