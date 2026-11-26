"""DA003 genetic diversity.

* Within-host: FWS per sample against the analysed samples' allele frequencies (descriptive).
* Population: nucleotide diversity per callable base (dominant-allele representation,
  sites callable in every analysed sample of the group) and mean expected heterozygosity,
  with a bootstrap over samples. SNP-only and SNP+indel versions are both produced.
"""

from __future__ import annotations

import random
from typing import Any

from pfvarx.configuration.loader import frozen_param
from pfvarx.constants import Terminal
from pfvarx.downstream.coi import fws
from pfvarx.downstream.common import Context, decisions, write_manifest
from pfvarx.genomic_evidence import read_bed
from pfvarx.io import write_tsv


def callable_intersection_bases(ctx: Context, samples: list[str]) -> int:
    """Number of bases callable in every sample (sweep over per-sample BEDs)."""
    per = [read_bed(ctx.results / "variants" / "samples" / f"{s}.callable.bed") for s in samples]
    total = 0
    chroms = set.intersection(*[set(p) for p in per]) if per else set()
    for c in chroms:
        events = []
        for p in per:
            for a, b in p[c]:
                events += [(a, 1), (b, -1)]
        events.sort()
        depth, last = 0, None
        for pos, d in events:
            if depth == len(per) and last is not None:
                total += pos - last
            depth += d
            last = pos
    return total


def dominant(ctx: Context, site: str, sid: str) -> int | None:
    c = ctx.matrix[site].get(sid)
    if not c or c["wsaf"] is None:
        return None
    return 1 if c["wsaf"] >= 0.5 else 0


def group_stats(ctx: Context, samples: list[str], include_indels: bool) -> dict[str, Any]:
    sites = [s for s in ctx.sites if include_indels or s["type"] == "SNV"]
    L = callable_intersection_bases(ctx, samples)
    pi_sum, he = 0.0, []
    used = 0
    for s in sites:
        alleles = [dominant(ctx, s["site_id"], x) for x in samples]
        if any(a is None for a in alleles):
            continue  # restricted to sites observed in every sample of the group
        n = len(alleles)
        if n < 2:
            continue
        p = sum(alleles) / n  # type: ignore[arg-type]
        h = n / (n - 1) * 2 * p * (1 - p)
        pi_sum += h
        he.append(h)
        used += 1
    return {"samples": len(samples), "sites_used": used, "callable_bases_intersection": L,
            "nucleotide_diversity_per_callable_base": (pi_sum / L) if L else None,
            "mean_expected_heterozygosity": (sum(he) / len(he)) if he else None}


def bootstrap(ctx: Context, samples: list[str], include_indels: bool, reps: int, seed: int,
              level: float, min_values: int) -> tuple[float | None, float | None]:
    """Percentile bootstrap interval at *level* (da003_interval_level) over samples; none with
    fewer than *min_values* replicates with a value (da003_min_bootstrap_values)."""
    rng = random.Random(seed)
    vals = []
    for _ in range(reps):
        draw = [rng.choice(samples) for _ in samples]
        if len(set(draw)) < 2:
            continue
        v = group_stats(ctx, draw, include_indels)["mean_expected_heterozygosity"]
        if v is not None:
            vals.append(v)
    if len(vals) < min_values:
        return None, None
    vals.sort()
    tail = (1 - level) / 2
    return vals[int(tail * len(vals))], vals[int((1 - tail) * len(vals)) - 1]


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA003")
    dec = decisions(ctx)["analyses"]["DA003"]
    eligible = sorted({s for c in dec["scales"].values() for s in c["eligible"]})
    excluded = {"INELIGIBLE": sorted(set(ctx.samples) - set(eligible))}
    col = ctx.eff["config"]["analyses"].get("population_column") or "population"
    snv_sites = [s for s in ctx.sites if s["type"] == "SNV"]
    # FWS needs population allele frequencies; with few samples each sets much of its own
    # reference (one sample: FWS is 0 by construction), so it is withheld below the minimum.
    population_ok = len(eligible) >= int(ctx.params["population_frequency_min_samples"])
    within = fws(ctx, snv_sites, eligible) if population_ok else {s: None for s in eligible}
    write_tsv(out / "within_host_fws.tsv", [{"sample_id": s, "fws": "" if within[s] is None else f"{within[s]:.4f}"}
                                            for s in eligible], ["sample_id", "fws"])
    groups: dict[str, list[str]] = {"ALL": eligible}
    for s in eligible:
        g = ctx.metadata.get(s, {}).get(col, "")
        if g:
            groups.setdefault(g, []).append(s)
    minimum = int(ctx.params["da003_min_group_size"])
    rows: list[dict[str, Any]] = []
    for g, members in sorted(groups.items()):
        for indels in (False, True):
            if len(members) < minimum:
                rows.append({"group": g, "variant_set": "SNP+indel" if indels else "SNP", "samples": len(members),
                             "status": "GROUP_TOO_SMALL"})
                continue
            st = group_stats(ctx, members, indels)
            lo, hi = bootstrap(ctx, members, indels, int(frozen_param(ctx.params, "da003_bootstrap_replicates")),
                               int(ctx.eff["config"]["random_seed"]),
                               float(frozen_param(ctx.params, "da003_interval_level")),
                               int(frozen_param(ctx.params, "da003_min_bootstrap_values")))
            rows.append({"group": g, "variant_set": "SNP+indel" if indels else "SNP", "status": "COMPUTED", **st,
                         "he_boot_lo": lo, "he_boot_hi": hi})
    cols = ["group", "variant_set", "status", "samples", "sites_used", "callable_bases_intersection",
            "nucleotide_diversity_per_callable_base", "mean_expected_heterozygosity", "he_boot_lo", "he_boot_hi"]
    write_tsv(out / "population_diversity.tsv", rows, cols)
    write_manifest(ctx, "DA003", (Terminal.COMPLETED_WITH_EXCLUSIONS if excluded["INELIGIBLE"] else Terminal.COMPLETED).value,
                   reason_codes=[], requested=ctx.samples, eligible=eligible, excluded=excluded, analysed=eligible,
                   outputs=[out / "within_host_fws.tsv", out / "population_diversity.tsv"],
                   caveats=["Population statistics use the dominant allele per sample (polyclonal samples simplified).",
                            "Sites restricted to those observed in every sample of the group; callable-base "
                            "denominator is the intersection of the samples' callability masks.",
                            "FWS is descriptive and depends on the reference group used."]
                   + ([] if population_ok else
                      [f"FWS withheld: fewer than {int(ctx.params['population_frequency_min_samples'])} samples "
                       "to estimate population allele frequencies."]),
                   method={"representation": "dominant allele", "bootstrap": "200 resamples of samples"},
                   results_summary={"groups": len(groups)})
