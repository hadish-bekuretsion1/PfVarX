"""Layer 3 cohort evidence matrix and missingness semantics.

The merged VCF is variant-only: an absent record is NOT a reference call. For
downstream estimators that need a genotype at every cohort site, each sample's own
reads are re-queried at the union of cohort sites (same admission filters as the
caller). The resulting per-sample allele counts are *site-forced per-sample evidence*,
not joint genotyping; positions outside a sample's callable region are UNCALLABLE and
positions below the minimum usable depth are MISSING.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import TypedDict

import pysam

from pfvarx.genomic_evidence import FLAG_FILTER, in_intervals, read_bed

STATES = ("ALT_DOMINANT", "MIXED", "REF_DOMINANT", "MISSING_LOW_DEPTH", "UNCALLABLE")


class CohortSite(TypedDict):
    chrom: str
    pos: int
    ref: str
    alt: str
    type: object
    n_pass: int
    statuses: dict[str, str | None]


def site_list(vcf: Path) -> list[CohortSite]:
    sites: list[CohortSite] = []
    with pysam.VariantFile(str(vcf)) as f:
        for rec in f:
            if rec.alts is None or rec.ref is None:
                continue
            statuses = {str(s): rec.samples[s].get("PST") for s in rec.samples}
            sites.append({"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0],
                          "type": rec.info.get("TYPE"), "n_pass": sum(1 for v in statuses.values() if v == "PASS"),
                          "statuses": statuses})
    return sites


def allele_counts(bam: pysam.AlignmentFile, site: CohortSite, min_mapq: int, min_bq: int) -> tuple[int, int, int]:
    """(ref, alt, depth) at one site; indels counted from CIGAR at the anchor."""
    chrom, pos = str(site["chrom"]), int(site["pos"]) - 1
    ref, alt = str(site["ref"]), str(site["alt"])
    kind = "SNV" if len(ref) == len(alt) == 1 else ("INS" if len(alt) > len(ref) else "DEL")
    r = a = d = 0
    for col in bam.pileup(chrom, pos, pos + 1, truncate=True, stepper="nofilter", flag_filter=FLAG_FILTER,
                          min_mapping_quality=min_mapq, min_base_quality=0, ignore_overlaps=True,
                          ignore_orphans=False, max_depth=1_000_000, compute_baq=False):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            sequence = pr.alignment.query_sequence
            qualities = pr.alignment.query_qualities
            if sequence is None or qualities is None:
                continue
            q = qualities[pr.query_position]
            if q < min_bq:
                continue
            d += 1
            base = sequence[pr.query_position].upper()
            if kind == "SNV":
                if base == alt:
                    a += 1
                elif base == ref:
                    r += 1
            else:
                ind = pr.indel
                want = len(alt) - len(ref)
                if ind == want and (kind == "DEL" or sequence[
                        pr.query_position + 1: pr.query_position + 1 + ind].upper() == alt[1:]):
                    a += 1
                elif ind == 0:
                    r += 1
    return r, a, d


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--bams", nargs="+", required=True)
    ap.add_argument("--callable", nargs="+", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--frozen-config", required=True)
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-missing", required=True)
    ap.add_argument("--out-sites", required=True)
    a = ap.parse_args(argv)
    eff = json.loads(Path(a.frozen_config).read_text())
    prm = {k: v["value"] for k, v in eff["parameters"].items()}
    min_mapq, min_bq = int(prm["layer2.min_mapping_quality"]), int(prm["layer2.min_base_quality"])
    min_depth = int(prm["layer2.min_usable_depth"])
    dom = float(prm["layer3.dominant_allele_fraction"])
    samples = Path(a.samples).read_text().split()
    sites = site_list(Path(a.vcf))
    with Path(a.out_sites).open("w") as h:
        h.write("site_id\tchrom\tpos\tref\talt\ttype\tn_pass\n")
        for i, s in enumerate(sites):
            h.write(f"S{i}\t{s['chrom']}\t{s['pos']}\t{s['ref']}\t{s['alt']}\t{s['type']}\t{s['n_pass']}\n")
    summary = []
    with gzip.open(a.out_matrix, "wt") as out:
        out.write("site_id\tsample_id\tstate\tref_count\talt_count\tdepth\twsaf\tcaller_status\n")
        for sid, bam_path, cal in zip(samples, a.bams, a.callable):
            callable_iv = read_bed(Path(cal))
            counts = {k: 0 for k in STATES}
            with pysam.AlignmentFile(bam_path) as bam:
                for i, s in enumerate(sites):
                    status = s["statuses"].get(sid) or "."
                    if not in_intervals(callable_iv.get(str(s["chrom"]), []), int(s["pos"]) - 1):
                        state, r, al, d, w = "UNCALLABLE", 0, 0, 0, ""
                    else:
                        r, al, d = allele_counts(bam, s, min_mapq, min_bq)
                        tot = r + al
                        if tot < min_depth:
                            state, w = "MISSING_LOW_DEPTH", ""
                        else:
                            f = al / tot
                            w = f"{f:.4f}"
                            state = "ALT_DOMINANT" if f >= dom else ("REF_DOMINANT" if f <= 1 - dom else "MIXED")
                    counts[state] += 1
                    out.write(f"S{i}\t{sid}\t{state}\t{r}\t{al}\t{d}\t{w}\t{status}\n")
            n = len(sites) or 1
            summary.append({"sample_id": sid, **counts,
                            "missing_fraction": (counts["MISSING_LOW_DEPTH"] + counts["UNCALLABLE"]) / n})
    cols = ["sample_id", *STATES, "missing_fraction"]
    with Path(a.out_missing).open("w") as h:
        h.write("\t".join(cols) + "\n")
        for row in summary:
            h.write("\t".join(f"{row[c]:.4f}" if c == "missing_fraction" else str(row[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
