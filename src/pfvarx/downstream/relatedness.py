"""DA004 relatedness and identity by descent.

Method is matched to clonal status (from DA001): pairs of MONOCLONAL samples are
analysed with hmmIBD on dominant-allele haploid genotypes; pairs including a POLYCLONAL
sample with Dcifer on distance-thinned allele-presence data, using DA001 COI medians.
INDETERMINATE samples were excluded at the eligibility checkpoint. Network thresholds in
the dashboard are display settings, not relatedness classifications.

Marker requirements depend on the assay. ``da004_min_informative_markers`` (about 200
biallelic markers keep the relatedness error below 0.1; Taylor et al. 2019) is the precision
level for every design and the minimum for WGS. Targeted panels have far fewer, partly linked
SNPs, so they use ``da004_min_informative_markers_targeted`` as the minimum. A pair between
the two is ESTIMATED_LOW_PRECISION: the estimate is reported without a relatedness label.
"""

from __future__ import annotations

import shutil
import subprocess
from itertools import combinations

from pfvarx.constants import Reason, Terminal
from pfvarx.downstream.coi import informative_sites
from pfvarx.downstream.common import Context, decisions, write_manifest
from pfvarx.downstream.diversity import dominant
from pfvarx.io import read_tsv, write_tsv


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA004")
    dec = decisions(ctx)["analyses"]["DA004"]
    eligible = sorted({s for c in dec["scales"].values() for s in c["eligible"]})
    excluded = {"INELIGIBLE": sorted(set(ctx.samples) - set(eligible))}
    clon = {r["sample_id"]: r for r in read_tsv(ctx.results / "downstream" / "DA001" / "clonality.tsv")}
    mono = [s for s in eligible if clon.get(s, {}).get("clonality") == "MONOCLONAL"]
    precise = int(ctx.params["da004_min_informative_markers"])
    targeted = ctx.design.targeted
    min_markers = int(ctx.params["da004_min_informative_markers_targeted"]) if targeted else precise

    def marker_status(n: int) -> str:
        if n < min_markers:
            return "EXCLUDED_FEW_MARKERS"
        return "ESTIMATED" if n >= precise else "ESTIMATED_LOW_PRECISION"
    floor = float(ctx.params["da001_allele_presence_floor"])
    related, clonal = float(ctx.params["da004_related_ibd"]), float(ctx.params["da004_clonal_ibd"])

    def label(estimate: str, status: str) -> str:
        """Reporting convention only; relatedness is continuous and the labels are not tests.
        Low-precision estimates get no label: their error can exceed the label boundaries."""
        if status != "ESTIMATED" or estimate in ("", "NA"):
            return ""
        value = float(estimate)
        return "NEAR_IDENTICAL" if value >= clonal else "CLOSELY_RELATED" if value >= related else "BELOW_RELATED_LEVEL"
    pairs_out = []
    # hmmIBD on monoclonal samples
    if len(mono) >= 2 and shutil.which("hmmIBD"):
        chrom_ids: dict[str, int] = {}
        lines = ["chrom\tpos\t" + "\t".join(mono)]
        for s in ctx.sites:
            if s["type"] != "SNV":
                continue
            g = [dominant(ctx, s["site_id"], x) for x in mono]
            cid = chrom_ids.setdefault(s["chrom"], len(chrom_ids) + 1)
            lines.append(f"{cid}\t{s['pos']}\t" + "\t".join("-1" if v is None else str(v) for v in g))
        geno = out / "hmmibd_input.txt"
        geno.write_text("\n".join(lines) + "\n")
        r = subprocess.run(["hmmIBD", "-i", str(geno), "-o", str(out / "hmmibd")], capture_output=True, text=True)
        (out / "hmmibd.log").write_text(r.stdout + r.stderr)
        if r.returncode != 0:
            raise SystemExit(f"hmmIBD failed with exit code {r.returncode}")
        for row in read_tsv(out / "hmmibd.hmm_fract.txt"):
            n = int(row["N_informative_sites"])
            pairs_out.append({"sample_1": row["sample1"], "sample_2": row["sample2"], "method": "hmmIBD",
                              "estimate": row["fract_sites_IBD"], "ci_lower": "", "ci_upper": "",
                              "informative_markers": n,
                              "status": marker_status(n)})
    # Dcifer for pairs with at least one polyclonal sample
    poly_pairs = [(a, b) for a, b in combinations(eligible, 2)
                  if "POLYCLONAL" in (clon.get(a, {}).get("clonality"), clon.get(b, {}).get("clonality"))]
    if poly_pairs and shutil.which("Rscript"):
        sites = informative_sites(ctx)
        rows = []
        for s in sites:
            for sid in eligible:
                c = ctx.matrix[s["site_id"]].get(sid)
                if not c or c["wsaf"] is None:
                    continue
                if c["wsaf"] <= 1 - floor:
                    rows.append({"sample_id": sid, "locus": s["site_id"], "allele": "R"})
                if c["wsaf"] >= floor:
                    rows.append({"sample_id": sid, "locus": s["site_id"], "allele": "A"})
        lf = out / "dcifer_input.tsv"
        write_tsv(lf, rows, ["sample_id", "locus", "allele"])
        # MOIRE's COI count where DA001 reported one; monoclonal samples are 1; otherwise "NA",
        # and the R script uses Dcifer's own naive estimate (getCOI) rather than assuming 1.
        coi_rows = [{"sample_id": s, "coi": clon.get(s, {}).get("coi_median")
                     or ("1" if clon.get(s, {}).get("clonality") == "MONOCLONAL" else "NA")} for s in eligible]
        write_tsv(out / "dcifer_coi.tsv", coi_rows, ["sample_id", "coi"])
        r = subprocess.run(["Rscript", str(ctx.root / "workflow" / "scripts" / "da004_dcifer.R"), str(lf),
                            str(out / "dcifer_coi.tsv"), str(out / "dcifer_pairs.tsv"),
                            str(ctx.params.get("da004_dcifer_alpha", 0.05)),
                            str(ctx.params.get("da004_dcifer_null_grid", 1000))], capture_output=True, text=True)
        (out / "dcifer.log").write_text(r.stdout + r.stderr)
        if r.returncode != 0:
            raise SystemExit(f"Dcifer failed with exit code {r.returncode}")
        wanted = {frozenset(p) for p in poly_pairs}
        for row in read_tsv(out / "dcifer_pairs.tsv"):
            if frozenset((row["sample_1"], row["sample_2"])) in wanted:
                pairs_out.append({"sample_1": row["sample_1"], "sample_2": row["sample_2"], "method": "Dcifer",
                                  "estimate": row["estimate"], "ci_lower": row["ci_lower"], "ci_upper": row["ci_upper"],
                                  "informative_markers": len(sites), "status": marker_status(len(sites))})
    for pair in pairs_out:
        pair["relatedness_label"] = label(str(pair["estimate"]), str(pair["status"]))
    cols = ["sample_1", "sample_2", "method", "estimate", "ci_lower", "ci_upper", "informative_markers", "status",
            "relatedness_label"]
    write_tsv(out / "pairwise_relatedness.tsv", pairs_out, cols)
    if not pairs_out:
        write_manifest(ctx, "DA004", Terminal.SKIPPED_INELIGIBLE.value, reason_codes=[Reason.NO_ELIGIBLE_PAIR.value],
                       requested=ctx.samples, eligible=eligible, excluded=excluded, analysed=[], outputs=[],
                       caveats=[], method={"tools": ["hmmIBD", "Dcifer"]})
        return
    counts = {s: sum(p["status"] == s for p in pairs_out)
              for s in ("ESTIMATED", "ESTIMATED_LOW_PRECISION", "EXCLUDED_FEW_MARKERS")}
    estimated = counts["ESTIMATED"] + counts["ESTIMATED_LOW_PRECISION"]
    caveats = ["hmmIBD estimates for monoclonal pairs have no interval in this release.",
               "Dcifer assumes approximately unlinked loci; distance thinning is a proxy."]
    if counts["ESTIMATED_LOW_PRECISION"]:
        caveats.append(f"{counts['ESTIMATED_LOW_PRECISION']} pair(s) have fewer than {precise} informative markers: "
                       "estimates are reported without relatedness labels because their error can exceed 0.1.")
    if not estimated:
        terminal, reasons, caveats = Terminal.SKIPPED_INELIGIBLE, [Reason.INSUFFICIENT_INFORMATIVE_SITES.value], [
            f"Every pair had fewer than {min_markers} informative markers; no relatedness was estimated."]
    elif excluded["INELIGIBLE"] or counts["EXCLUDED_FEW_MARKERS"]:
        terminal, reasons = Terminal.COMPLETED_WITH_EXCLUSIONS, []
    else:
        terminal, reasons = Terminal.COMPLETED, []
    write_manifest(ctx, "DA004", terminal.value,
                   reason_codes=reasons, requested=ctx.samples, eligible=eligible, excluded=excluded,
                   analysed=eligible if estimated else [], outputs=[out / "pairwise_relatedness.tsv"],
                   caveats=caveats,
                   method={"monoclonal": "hmmIBD", "polyclonal": "Dcifer", "min_markers": min_markers,
                           "precise_markers": precise, "marker_rule": "TARGETED" if targeted else "WGS",
                           "closely_related_ibd": related, "near_identical_ibd": clonal},
                   results_summary={"pairs": len(pairs_out), **{f"pairs_{k.lower()}": v for k, v in counts.items()}})
