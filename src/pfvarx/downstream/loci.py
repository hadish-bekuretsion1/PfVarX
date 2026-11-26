"""The loci MOIRE (DA001) and Dcifer (DA004) are given, and which kind they are.

Both tools take long-form allele presence (``sample_id``, ``locus``, ``allele``) and both are
built for multiallelic loci. Two locus sets can supply that:

* ``MICROHAPLOTYPE`` — the Layer 3 microhaplotype artifact (``cohort/microhaplotypes.tsv.gz``):
  one locus per catalogue locus, one allele per haplotype the sample carries. This is what the
  tools expect, and a locus with several haplotypes separates clones far better than its sites
  taken one at a time.
* ``BIALLELIC_SNV`` — distance-thinned PASS SNVs, each contributing the alleles present in the
  sample (``R``, ``A``). This is what runs without a microhaplotype panel use.

Microhaplotypes are used when the cohort has at least ``layer3.microhap_min_loci_for_basis`` of
them; with fewer, the larger biallelic set is kept, because a handful of loci is worse than many
sites even when each locus carries more information. The choice is never silent: the basis, the
locus count and the reason are recorded in the analysis manifest.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from typing import Any

from pfvarx.configuration.loader import frozen_param
from pfvarx.downstream.common import Context

MICROHAPLOTYPE = "MICROHAPLOTYPE"
BIALLELIC_SNV = "BIALLELIC_SNV"
ARTIFACT = ("cohort", "microhaplotypes.tsv.gz")


@dataclass
class AlleleSet:
    """Long-form allele presence for one analysis, with how it was obtained."""

    rows: list[dict[str, str]]
    loci_per_sample: dict[str, int]
    basis: str
    loci: int
    reason: str = ""
    detail: dict[str, Any] = field(default_factory=dict)

    @property
    def provenance(self) -> dict[str, Any]:
        """What to record in the analysis manifest's method block."""
        return {"locus_basis": self.basis, "loci": self.loci, "basis_reason": self.reason, **self.detail}


def microhaplotype_rows(ctx: Context, eligible: list[str]) -> tuple[list[dict[str, str]], dict[str, int], int]:
    """Called haplotypes of the eligible samples, as long-form rows."""
    path = ctx.results.joinpath(*ARTIFACT)
    rows: list[dict[str, str]] = []
    per_sample = {sample: 0 for sample in eligible}
    keep = set(eligible)
    loci: set[str] = set()
    if not path.is_file():
        return rows, per_sample, 0
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row.get("state") != "CALLED" or row["sample_id"] not in keep or not row.get("haplotype"):
                continue
            rows.append({"sample_id": row["sample_id"], "locus": row["locus_id"], "allele": row["haplotype"]})
            loci.add(row["locus_id"])
    for locus_sample in {(r["locus"], r["sample_id"]) for r in rows}:
        per_sample[locus_sample[1]] += 1
    return rows, per_sample, len(loci)


def biallelic_rows(ctx: Context, eligible: list[str], floor: float
                   ) -> tuple[list[dict[str, str]], dict[str, int], int]:
    """Allele presence at distance-thinned PASS SNVs: the reference and alternative alleles
    whose within-sample fraction clears the presence floor."""
    from pfvarx.downstream.coi import informative_sites

    rows: list[dict[str, str]] = []
    per_sample = {sample: 0 for sample in eligible}
    sites = informative_sites(ctx)
    for site in sites:
        for sample in eligible:
            cell = ctx.matrix[site["site_id"]].get(sample)
            if not cell or cell["wsaf"] is None:
                continue
            per_sample[sample] += 1
            if cell["wsaf"] <= 1 - floor:
                rows.append({"sample_id": sample, "locus": site["site_id"], "allele": "R"})
            if cell["wsaf"] >= floor:
                rows.append({"sample_id": sample, "locus": site["site_id"], "allele": "A"})
    return rows, per_sample, len(sites)


def allele_set(ctx: Context, eligible: list[str], floor: float) -> AlleleSet:
    """Microhaplotypes when the cohort has enough of them, else biallelic SNVs; never silently."""
    minimum = int(frozen_param(ctx.params, "microhap_min_loci_for_basis"))
    rows, per_sample, loci = microhaplotype_rows(ctx, eligible)
    if loci >= minimum and rows:
        haplotypes = len({(r["locus"], r["allele"]) for r in rows})
        return AlleleSet(rows, per_sample, MICROHAPLOTYPE, loci,
                         f"{loci} microhaplotype loci called (minimum {minimum})",
                         {"haplotypes": haplotypes,
                          "mean_haplotypes_per_locus": round(haplotypes / loci, 2)})
    reason = ("no microhaplotype locus was called, so biallelic SNVs are used" if not loci else
              f"only {loci} microhaplotype loci were called (minimum {minimum}), so the larger "
              "biallelic SNV set is used")
    snv_rows, snv_per_sample, sites = biallelic_rows(ctx, eligible, floor)
    return AlleleSet(snv_rows, snv_per_sample, BIALLELIC_SNV, sites, reason,
                     {"microhaplotype_loci_available": loci})
