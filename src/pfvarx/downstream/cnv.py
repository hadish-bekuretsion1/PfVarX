"""DA008 copy-number markers: gene amplification and deletion at catalogue loci (WGS only).

For each locus of the copy-number catalogue (default: pfmdr1, plasmepsin II and III and gch1
amplification; pfhrp2 and pfhrp3 deletion) the gene's read depth is compared with the
sample's median core-genome depth (both excluding duplicates and low mapping quality):

* depth ratio  = observed / expected depth over the gene (about 1 = one copy, about 2 = two
  copies, about 0 = absent), with a bootstrap interval over fixed windows. The expected
  depth of each window comes from core-genome windows of the same sample with the same GC
  content: P. falciparum is about 80 % A+T and libraries cover AT-rich stretches unevenly,
  differently in each library, so a plain ratio to the genome median would mistake that
  bias for copy-number change;
* flank ratio  = the same for the regions either side of the gene, so a deleted gene is
  told apart from a region that is simply not covered;
* read pairs   = outward-facing pairs near the gene (the signature of a tandem duplication
  junction) and pairs spanning the whole gene (the signature of a deletion), reported as
  supporting evidence.

A locus that cannot be judged is UNKNOWN with its reason, never SINGLE_COPY or PRESENT.
Ratios are DNA fractions: in a mixed infection a duplication carried by one clone gives a
ratio between 1 and 2. Thresholds are UNCALIBRATED_PROVISIONAL until benchmarked.
The catalogue must be MAINTAINER_VERIFIED (or TEST_FIXTURE), as for DA002.
"""

from __future__ import annotations

import random
import statistics
from pathlib import Path
from typing import Any

import pysam

from pfvarx.configuration.loader import frozen_param
from pfvarx.constants import Reason, Terminal
from pfvarx.downstream.common import Context, decisions, wilson, write_manifest
from pfvarx.downstream.resistance import load_catalog
from pfvarx.genomic_evidence import FLAG_FILTER
from pfvarx.io import write_json, write_tsv
from pfvarx.quality.metrics import mosdepth_region_dist

# Layer 5 memory: the cohort evidence matrix is loaded for this analysis only when it reads it
# (measures read depth from the alignments, not the matrix). tests/unit/test_matrix_declarations.py checks the declaration against the code.
USES_COHORT_MATRIX = False
CATALOG = "config/registry/cnv_catalog_pfvarx_v1.tsv"
CATALOG_COLUMNS = ("locus_id", "gene_id", "gene_name", "event", "meaning", "evidence_class", "source")
EVENTS = ("AMPLIFICATION", "DELETION")


def catalog_path(ctx: Context) -> Path:
    return ctx.root / CATALOG


def gene_location(gff: Path, gene_id: str) -> tuple[str, int, int] | None:
    """(contig, 0-based start, end) of a gene from its GFF3 gene feature, else from its CDS."""
    cds: list[tuple[str, int, int]] = []
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split("\t")
        if len(c) < 9:
            continue
        attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        if c[2].endswith("gene") and attrs.get("ID") == gene_id:
            return c[0], int(c[3]) - 1, int(c[4])
        if c[2] == "CDS" and attrs.get("Parent", "").split(".")[0] == gene_id:
            cds.append((c[0], int(c[3]) - 1, int(c[4])))
    if cds:
        return cds[0][0], min(x[1] for x in cds), max(x[2] for x in cds)
    return None


def depth_profile(bam: pysam.AlignmentFile, chrom: str, start: int, end: int, min_mapq: int) -> list[int]:
    """Per-base read depth (duplicates, secondary, supplementary and low-MAPQ reads excluded)."""
    if end <= start:
        return []

    def keep(read: pysam.AlignedSegment) -> bool:
        return not (read.flag & FLAG_FILTER) and read.mapping_quality >= min_mapq

    counts = bam.count_coverage(chrom, start, end, quality_threshold=0, read_callback=keep)
    return [sum(base[i] for base in counts) for i in range(end - start)]


def gc_fraction(seq: str) -> float:
    acgt = sum(seq.count(b) for b in "ACGTacgt")
    return (sum(seq.count(b) for b in "GCgc") / acgt) if acgt else 0.0


class GcModel:
    """Expected depth of a window from its GC content, learned from random core-genome windows of one sample."""

    def __init__(self, bam: pysam.AlignmentFile, fasta: pysam.FastaFile, core: list[tuple[str, int, int]],
                 window: int, n: int, min_mapq: int, seed: int, min_per_bin: int):
        rng = random.Random(seed)
        spans = [(c, a, b) for c, a, b in core if b - a > window and c in bam.references]
        weights = [b - a for _, a, b in spans]
        self.points: list[tuple[float, float]] = []
        for _ in range(n if spans else 0):
            c, a, b = rng.choices(spans, weights=weights)[0]
            st = rng.randrange(a, b - window)
            depth = depth_profile(bam, c, st, st + window, min_mapq)
            if depth:
                self.points.append((gc_fraction(fasta.fetch(c, st, st + window)), float(statistics.median(depth))))
        self.points.sort()
        self.min_per_bin = min_per_bin

    def expected(self, gc: float) -> float | None:
        """Median depth of the core windows closest in GC (at least ``min_per_bin`` of them)."""
        if len(self.points) < self.min_per_bin:
            return None
        nearest = sorted(self.points, key=lambda x: abs(x[0] - gc))[: self.min_per_bin]
        return float(statistics.median(d for _, d in nearest))


def window_ratios(bam: pysam.AlignmentFile, fasta: pysam.FastaFile, chrom: str, start: int, end: int,
                  window: int, min_mapq: int, model: GcModel) -> list[float]:
    """Observed / GC-expected depth of each full window of a region."""
    out = []
    depth = depth_profile(bam, chrom, start, end, min_mapq)
    for i in range(0, len(depth), window):
        w = depth[i:i + window]
        if len(w) < max(1, window // 2):
            continue
        exp = model.expected(gc_fraction(fasta.fetch(chrom, start + i, start + i + len(w))))
        if exp and exp > 0:
            out.append(float(statistics.median(w)) / exp)
    return out


def window_medians(depth: list[int], window: int) -> list[float]:
    """Median depth of each full window (a last partial window of at least half a window is kept)."""
    out = []
    for i in range(0, len(depth), window):
        w = depth[i:i + window]
        if len(w) >= max(1, window // 2):
            out.append(float(statistics.median(w)))
    return out


def ratio_interval(windows: list[float], core: float, reps: int, seed: int,
                   level: float) -> tuple[float | None, float | None]:
    """Percentile bootstrap interval at *level* (da008_interval_level; 0.95 gives 2.5 % and
    97.5 %) of the median window value divided by ``core`` (1 when the windows are already ratios)."""
    if not windows or core <= 0:
        return None, None
    rng = random.Random(seed)
    stats = sorted(statistics.median(rng.choices(windows, k=len(windows))) / core for _ in range(reps))
    tail = (1 - level) / 2
    return stats[int(tail * (reps - 1))], stats[int((1 - tail) * (reps - 1))]


def pair_support(bam: pysam.AlignmentFile, chrom: str, start: int, end: int, flank: int,
                 min_mapq: int) -> dict[str, int]:
    """Read pairs that support a duplication junction (outward-facing) or a deletion (spanning the gene)."""
    lo, hi = max(0, start - flank), min(bam.get_reference_length(chrom), end + flank)
    outward = spanning = 0
    for read in bam.fetch(chrom, lo, hi):
        if (read.flag & FLAG_FILTER or read.mapping_quality < min_mapq or not read.is_paired or not read.is_read1
                or read.mate_is_unmapped or read.next_reference_id != read.reference_id):
            continue
        left_is_read = read.reference_start <= read.next_reference_start
        left_reverse = read.is_reverse if left_is_read else read.mate_is_reverse
        right_reverse = read.mate_is_reverse if left_is_read else read.is_reverse
        if left_reverse and not right_reverse:
            outward += 1
        left_pos = min(read.reference_start, read.next_reference_start)
        right_pos = max(read.reference_start, read.next_reference_start)
        if not left_reverse and right_reverse and left_pos < start and right_pos >= end:
            spanning += 1
    return {"outward_pairs": outward, "spanning_pairs": spanning}


def call_locus(event: str, ratio: float | None, ratio_lo: float | None, flank_ratio: float | None,
               p: dict[str, Any]) -> tuple[str, str]:
    """(state, reason) from the depth ratios; reason is empty for a confident call."""
    if ratio is None:
        return "UNKNOWN", "NO_DEPTH"
    if event == "AMPLIFICATION":
        if ratio >= float(p["da008_amplification_ratio"]) and (ratio_lo or 0) >= float(p["da008_amplification_lower"]):
            return "AMPLIFIED", ""
        if ratio < float(p["da008_present_ratio"]):
            return "UNKNOWN", "GENE_DEPTH_LOW"
        return "SINGLE_COPY", ""
    if ratio >= float(p["da008_present_ratio"]):
        return "PRESENT", ""
    if flank_ratio is None or flank_ratio < float(p["da008_flank_min_ratio"]):
        return "UNKNOWN", "FLANKS_NOT_COVERED"
    if ratio <= float(p["da008_deletion_ratio"]):
        return "DELETED", ""
    return "UNKNOWN", "AMBIGUOUS_DEPTH"


def measure(bam: pysam.AlignmentFile, loc: tuple[str, int, int], core: float, p: dict[str, Any],
            min_mapq: int, seed: int, fasta: pysam.FastaFile | None = None, model: GcModel | None = None
            ) -> dict[str, Any]:
    """Depth ratios of a gene and its flanks: GC-corrected with a model, else relative to ``core``."""
    chrom, start, end = loc
    flank, window = int(p["da008_flank_bp"]), int(p["da008_window_bp"])
    length = bam.get_reference_length(chrom)
    gene_depth = window_medians(depth_profile(bam, chrom, start, end, min_mapq), window)
    if model is not None and fasta is not None:
        gene = window_ratios(bam, fasta, chrom, start, end, window, min_mapq, model)
        flanks = (window_ratios(bam, fasta, chrom, max(0, start - flank), start, window, min_mapq, model)
                  + window_ratios(bam, fasta, chrom, end, min(length, end + flank), window, min_mapq, model))
        scale = 1.0
    else:
        gene = gene_depth
        flanks = (window_medians(depth_profile(bam, chrom, max(0, start - flank), start, min_mapq), window)
                  + window_medians(depth_profile(bam, chrom, end, min(length, end + flank), min_mapq), window))
        scale = core
    ratio = statistics.median(gene) / scale if gene and scale > 0 else None
    lo, hi = ratio_interval(gene, scale, int(p["da008_bootstrap"]), seed,
                            float(frozen_param(p, "da008_interval_level")))
    flank_ratio = statistics.median(flanks) / scale if flanks and scale > 0 else None
    return {"gene_depth": None if not gene_depth else statistics.median(gene_depth), "depth_ratio": ratio, "ratio_lo": lo,
            "ratio_hi": hi, "flank_ratio": flank_ratio, **pair_support(bam, chrom, start, end, flank, min_mapq)}


def fmt(x: float | None) -> str:
    return "" if x is None else f"{x:.3f}"


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA008")
    dec = decisions(ctx)["analyses"]["DA008"]
    eligible = sorted({s for c in dec["scales"].values() for s in c["eligible"]})
    excluded = {"INELIGIBLE": sorted(set(ctx.samples) - set(eligible))}
    cat_path = catalog_path(ctx)
    meta, loci = load_catalog(cat_path)
    status = meta.get("status", "UNKNOWN")
    if status not in ("MAINTAINER_VERIFIED", "TEST_FIXTURE"):
        write_manifest(ctx, "DA008", Terminal.BLOCKED_MISSING_PREREQUISITE.value,
                       reason_codes=[Reason.CATALOG_NOT_VERIFIED.value], requested=ctx.samples, eligible=eligible,
                       excluded=excluded, analysed=[], outputs=[],
                       caveats=[f"catalogue {cat_path.name} has status {status}; verify its loci against the "
                                "literature (gene IDs and the meaning of each event) before relying on DA008"],
                       method={"catalog": str(cat_path), "catalog_meta": meta})
        return
    p = ctx.params
    min_mapq = int(ctx.eff["parameters"]["layer2.min_mapping_quality"]["value"])
    min_depth = int(ctx.eff["parameters"]["layer2.min_usable_depth"]["value"])
    min_core = float(p["da008_min_core_depth"])
    seed = int(ctx.eff["config"].get("random_seed") or 1)
    gff = ctx.root / ctx.eff["config"]["reference"]["gff"] if ctx.eff["config"]["reference"].get("gff") else None
    places = {lc["locus_id"]: (gene_location(gff, lc["gene_id"]) if gff and gff.is_file() else None) for lc in loci}
    pop_col = ctx.eff["config"]["analyses"].get("population_column") or "population"
    fasta_path = ctx.root / ctx.plan["reference"]["parasite_fasta"]
    core_bed = ctx.root / ctx.eff["config"]["reference"]["core_mask_bed"]
    core_spans = [(c[0], int(c[1]), int(c[2])) for c in (ln.split("\t") for ln in core_bed.read_text().splitlines()
                                                         if ln.strip() and not ln.startswith(("#", "track")))]
    calls: list[dict[str, Any]] = []
    for sid in eligible:
        dist = ctx.results / "qc" / "samples" / sid / "core.mosdepth.region.dist.txt"
        core = float(mosdepth_region_dist(dist, min_depth).get("median_depth_core", 0)) if dist.is_file() else 0.0
        with pysam.AlignmentFile(str(ctx.results / "alignments" / f"{sid}.bam")) as bam, \
                pysam.FastaFile(str(fasta_path)) as fasta:
            model = GcModel(bam, fasta, core_spans, int(p["da008_window_bp"]), int(p["da008_gc_windows"]),
                            min_mapq, seed, int(frozen_param(p, "da008_gc_min_windows_per_bin"))) \
                if core >= min_core else None
            for i, lc in enumerate(loci):
                loc = places[lc["locus_id"]]
                row: dict[str, Any] = {"sample_id": sid, "locus_id": lc["locus_id"], "gene_id": lc["gene_id"],
                                       "gene": lc["gene_name"], "event": lc["event"],
                                       "location": "" if loc is None else f"{loc[0]}:{loc[1] + 1}-{loc[2]}",
                                       "core_median_depth": fmt(core)}
                m: dict[str, Any] = {}
                if loc is None:
                    state, reason = "UNKNOWN", "GENE_NOT_IN_GFF"
                elif core < min_core:
                    state, reason = "UNKNOWN", "LOW_CORE_DEPTH"
                else:
                    m = measure(bam, loc, core, p, min_mapq, seed + i, fasta, model)
                    state, reason = call_locus(lc["event"], m["depth_ratio"], m["ratio_lo"], m["flank_ratio"], p)
                row.update({"state": state, "reason": reason, "gene_median_depth": fmt(m.get("gene_depth")),
                            "depth_ratio": fmt(m.get("depth_ratio")), "ratio_lo": fmt(m.get("ratio_lo")),
                            "ratio_hi": fmt(m.get("ratio_hi")), "flank_ratio": fmt(m.get("flank_ratio")),
                            "outward_pairs": m.get("outward_pairs", ""), "spanning_pairs": m.get("spanning_pairs", ""),
                            "population": ctx.metadata.get(sid, {}).get(pop_col, ""),
                            "country": ctx.metadata.get(sid, {}).get("country", ""),
                            "meaning": lc.get("meaning", ""), "evidence_class": lc.get("evidence_class", "")})
                calls.append(row)
    cols = list(calls[0]) if calls else ["sample_id", "locus_id", "state"]
    write_tsv(out / "cnv_calls.tsv", calls, cols)
    prev = prevalence(calls, loci)
    write_tsv(out / "prevalence.tsv", prev, list(prev[0]) if prev else ["locus_id"])
    write_json(out / "catalog_provenance.json", {"path": str(cat_path), "meta": meta, "loci": len(loci)})
    caveats = ["Copy number is estimated from read depth relative to core-genome windows of the same sample with "
               "similar GC content; it is a DNA-level estimate, not a measured copy count.",
               "In a mixed infection a change carried by one clone gives an intermediate ratio (for example 1.3).",
               "GC bias is corrected from the sample's own core genome; mapping bias is not. pfhrp2 and pfhrp3 lie "
               "in subtelomeric, repetitive regions where depth is less reliable.",
               "Thresholds (amplification ratio, deletion ratio, minimum depth) are UNCALIBRATED_PROVISIONAL.",
               "An UNKNOWN locus could not be judged; it is never counted as single copy or present.",
               "pfhrp2/pfhrp3 deletion status supports RDT surveillance; it is not a diagnosis of an individual infection."]
    if status == "TEST_FIXTURE":
        caveats.insert(0, "TEST FIXTURE CATALOGUE: results are for software testing, not for surveillance.")
    unknown = sum(1 for c in calls if c["state"] == "UNKNOWN")
    write_manifest(ctx, "DA008", (Terminal.COMPLETED_WITH_EXCLUSIONS if excluded["INELIGIBLE"] or unknown
                                  else Terminal.COMPLETED).value,
                   reason_codes=[], requested=ctx.samples, eligible=eligible, excluded=excluded, analysed=eligible,
                   outputs=[out / "cnv_calls.tsv", out / "prevalence.tsv", out / "catalog_provenance.json"],
                   caveats=caveats,
                   method={"catalog": meta, "normalisation": f"expected depth from {int(p['da008_gc_windows'])} random "
                                                              "core-genome windows of the same sample with the nearest GC content",
                           "interval": f"percentile bootstrap over {int(p['da008_window_bp'])} bp windows",
                           "amplification_ratio": float(p["da008_amplification_ratio"]),
                           "amplification_lower": float(p["da008_amplification_lower"]),
                           "deletion_ratio": float(p["da008_deletion_ratio"]),
                           "present_ratio": float(p["da008_present_ratio"]), "min_core_depth": min_core},
                   results_summary={"locus_calls": len(calls), "unknown_calls": unknown,
                                    "amplified": sum(1 for c in calls if c["state"] == "AMPLIFIED"),
                                    "deleted": sum(1 for c in calls if c["state"] == "DELETED")})


def prevalence(calls: list[dict[str, Any]], loci: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Share of judged infections carrying each locus's event, per population and overall."""
    rows = []
    groups = sorted({c["population"] for c in calls if c["population"]}) + ["ALL"]
    for lc in loci:
        for g in groups:
            sel = [c for c in calls if c["locus_id"] == lc["locus_id"] and (g == "ALL" or c["population"] == g)]
            judged = [c for c in sel if c["state"] != "UNKNOWN"]
            k = sum(1 for c in judged if c["state"] in ("AMPLIFIED", "DELETED"))
            n = len(judged)
            lo, hi = wilson(k, n)
            rows.append({"locus_id": lc["locus_id"], "gene": lc["gene_name"], "event": lc["event"], "group": g,
                         "samples": len(sel), "judged": n, "unknown": len(sel) - n, "with_event": k,
                         "prevalence": "" if not n else f"{k / n:.4f}",
                         "wilson_lo": "" if lo is None else f"{lo:.4f}", "wilson_hi": "" if hi is None else f"{hi:.4f}"})
    return rows
