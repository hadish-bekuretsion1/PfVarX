"""DA001 complexity of infection and clonal status.

Estimator: MOIRE run on allele-presence calls at distance-thinned biallelic PASS SNV
sites. Classification rule (Section 3.10): POLYCLONAL if P(COI > 1) >= credibility,
MONOCLONAL if P(COI > 1) <= 1 - credibility, otherwise INDETERMINATE; samples with too
few informative loci are INDETERMINATE. FWS is reported as a descriptive index only, and
the PfVarX-Call technical mixture K is shown as context, never as COI.

MOIRE needs a population: with fewer than ``population_frequency_min_samples`` informative
samples it is not run, the monoclonal/polyclonal call comes from the caller's per-sample
mixture test (``basis`` CALLER_MIXTURE_LRT), and the COI count and FWS are withheld.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from typing import Any

from pfvarx.constants import Clonality, Reason, Terminal
from pfvarx.downstream.common import Context, write_manifest
from pfvarx.io import read_tsv, write_tsv


def informative_sites(ctx: Context) -> list[dict[str, str]]:
    thin = int(ctx.params["da001_thin_bp"])
    keep: list[dict[str, str]] = []
    for s in ctx.sites:
        if s["type"] != "SNV" or int(s["n_pass"]) == 0:
            continue
        if keep and keep[-1]["chrom"] == s["chrom"] and int(s["pos"]) - int(keep[-1]["pos"]) < thin:
            continue
        keep.append(s)
    return keep


def fws(ctx: Context, sites: list[dict[str, str]], samples: list[str]) -> dict[str, float | None]:
    """FWS = 1 - sum(Hw)/sum(Hs) across sites, population frequencies from the analysed samples."""
    out: dict[str, float | None] = {}
    pop: dict[str, float] = {}
    for s in sites:
        vals = [c["wsaf"] for sid in samples if (c := ctx.matrix[s["site_id"]].get(sid)) and c["wsaf"] is not None]
        if vals:
            pop[s["site_id"]] = sum(vals) / len(vals)
    for sid in samples:
        hw = hs = 0.0
        for s in sites:
            c = ctx.matrix[s["site_id"]].get(sid)
            if not c or c["wsaf"] is None or s["site_id"] not in pop:
                continue
            q, p = c["wsaf"], pop[s["site_id"]]
            hw += 2 * q * (1 - q)
            hs += 2 * p * (1 - p)
        out[sid] = None if hs == 0 else 1 - hw / hs
    return out


def _coi(row: dict[str, str] | None, key: str) -> str:
    """MOIRE COI summaries are posterior quantiles that R may interpolate (2.02499999999998);
    report them to two decimals, without trailing zeros (2.02, 1)."""
    if not row or row.get(key) in (None, "", "NA"):
        return ""
    try:
        return f"{float(row[key]):.2f}".rstrip("0").rstrip(".")
    except ValueError:
        return row[key]


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA001")
    floor = float(ctx.params["da001_allele_presence_floor"])
    cred = float(ctx.params["da001_credibility"])
    min_loci = int(ctx.params["da001_min_informative_loci"])
    min_population = int(ctx.params["population_frequency_min_samples"])
    eligible = [s for s in ctx.samples if ctx.gate(s)["permitted"].get("DA001")]
    excluded = {"QC_GATE": [s for s in ctx.samples if s not in eligible]}
    sites = informative_sites(ctx)
    rows = []
    loci_per_sample = {s: 0 for s in eligible}
    for s in sites:
        for sid in eligible:
            c = ctx.matrix[s["site_id"]].get(sid)
            if not c or c["wsaf"] is None:
                continue
            loci_per_sample[sid] += 1
            if c["wsaf"] <= 1 - floor:
                rows.append({"sample_id": sid, "locus": s["site_id"], "allele": "R"})
            if c["wsaf"] >= floor:
                rows.append({"sample_id": sid, "locus": s["site_id"], "allele": "A"})
    low_information = [sid for sid in eligible if loci_per_sample[sid] < min_loci]
    low_set = set(low_information)
    rows = [row for row in rows if row["sample_id"] not in low_set]
    excluded["INSUFFICIENT_INFORMATIVE_SITES"] = low_information
    long_form = out / "moire_input.tsv"
    write_tsv(long_form, rows, ["sample_id", "locus", "allele"])
    if not eligible or not rows:
        write_manifest(ctx, "DA001", Terminal.SKIPPED_INELIGIBLE.value,
                       reason_codes=[Reason.INSUFFICIENT_INFORMATIVE_SITES.value], requested=ctx.samples,
                       eligible=eligible, excluded=excluded, analysed=[], outputs=[], caveats=[],
                       method={"estimator": "MOIRE"})
        write_tsv(out / "clonality.tsv", [], ["sample_id", "clonality"])
        return
    caller_mix: dict[str, dict[str, Any]] = {}
    for sid in eligible:
        qc = ctx.results / "qc" / "samples" / sid / f"{sid}.caller_qc.json"
        if qc.is_file():
            caller_mix[sid] = json.loads(qc.read_text()).get("mixture") or {}
    informative = [sid for sid in eligible if loci_per_sample[sid] >= min_loci]
    # MOIRE estimates population allele frequencies jointly with each sample's COI. With few
    # samples neither is identified: on a 9-sample synthetic cohort the posterior P(COI>1) of
    # truly single-clone samples ranged from 0.00 to 1.00 across random seeds. Below the
    # minimum, MOIRE is not run and the per-sample call comes from the caller's mixture test
    # (bootstrap likelihood-ratio test of one against several genotypes, which needs no
    # population); the COI count and FWS are then withheld.
    population_ok = len(informative) >= min_population
    coi_tsv, af_tsv = out / "moire_coi.tsv", out / "moire_allele_freqs.tsv"
    summary: dict[str, dict[str, str]] = {}
    if population_ok:
        rscript = shutil.which("Rscript")
        if rscript is None:
            raise SystemExit("METHOD_NOT_INSTALLED: Rscript with the moire package is required for DA001")
        script = ctx.root / "workflow" / "scripts" / "da001_moire.R"
        r = subprocess.run([rscript, str(script), str(long_form), str(coi_tsv), str(af_tsv),
                            str(ctx.params["da001_moire_burnin"]), str(ctx.params["da001_moire_samples"]),
                            str(ctx.eff["config"]["random_seed"])], capture_output=True, text=True)
        (out / "moire.log").write_text(r.stdout + r.stderr)
        if r.returncode != 0:
            raise SystemExit(f"MOIRE failed (exit {r.returncode}); see {out / 'moire.log'}")
        summary = {row["sample_id"]: row for row in read_tsv(coi_tsv)}
    else:
        for stale in (coi_tsv, af_tsv, out / "moire.log"):
            stale.unlink(missing_ok=True)
    fw = fws(ctx, sites, eligible) if population_ok else {}
    out_rows = []
    for sid in eligible:
        mix = caller_mix.get(sid, {})
        ppoly: float | None = None
        srow = None
        if loci_per_sample[sid] < min_loci:
            clon, reason, basis = Clonality.INDETERMINATE, Reason.INSUFFICIENT_INFORMATIVE_SITES.value, ""
        elif population_ok:
            srow = summary.get(sid)
            basis = "MOIRE"
            if srow is None:
                clon, reason = Clonality.INDETERMINATE, "NOT_IN_MOIRE_OUTPUT"
            else:
                ppoly = float(srow["prob_polyclonal"])
                reason = ""
                if ppoly >= cred:
                    clon = Clonality.POLYCLONAL
                elif ppoly <= 1 - cred:
                    clon = Clonality.MONOCLONAL
                else:
                    clon, reason = Clonality.INDETERMINATE, "CREDIBLE_INTERVAL_SPANS_ONE"
        elif mix.get("status") == "ESTIMATED" and mix.get("k_selected") is not None:
            basis = "CALLER_MIXTURE_LRT"
            clon = Clonality.MONOCLONAL if int(mix["k_selected"]) == 1 else Clonality.POLYCLONAL
            reason = "MOIRE_NEEDS_POPULATION"
        else:
            clon, reason, basis = Clonality.INDETERMINATE, "MOIRE_NEEDS_POPULATION;CALLER_MIXTURE_UNAVAILABLE", ""
        sample_fws = fw.get(sid)
        out_rows.append({
            "sample_id": sid, "clonality": clon.value, "basis": basis, "reason": reason,
            "prob_polyclonal": "" if ppoly is None else f"{ppoly:.4f}",
            "coi_median": _coi(srow, "post_coi_med"), "coi_lower": _coi(srow, "post_coi_lower"),
            "coi_upper": _coi(srow, "post_coi_upper"), "informative_loci": loci_per_sample[sid],
            "fws": "" if sample_fws is None else f"{sample_fws:.4f}",
            "fws_ge_0_95_reference_only": "" if sample_fws is None else str(sample_fws >= float(ctx.params["da001_fws_reference_cutoff"])),
            "caller_technical_k_context_only": "" if mix.get("k_selected") is None else mix["k_selected"]})
    cols = list(out_rows[0])
    write_tsv(out / "clonality.tsv", out_rows, cols)
    counts: dict[str, int] = {}
    for r_ in out_rows:
        clonality = str(r_["clonality"])
        counts[clonality] = counts.get(clonality, 0) + 1
    assessed = [str(row["sample_id"]) for row in out_rows if row["basis"]]
    excluded["INSUFFICIENT_INFORMATIVE_SITES"] = [sid for sid in eligible if sid not in assessed]
    terminal = (Terminal.SKIPPED_INELIGIBLE if not assessed else
                Terminal.COMPLETED_WITH_EXCLUSIONS if any(excluded.values()) else Terminal.COMPLETED)
    outputs = [out / "clonality.tsv"] + ([coi_tsv, af_tsv] if population_ok else [])
    caveats = (["MOIRE parameters are UNCALIBRATED_PROVISIONAL; convergence not assessed.",
                "Classification rule uses P(COI>1) at the declared credibility level.",
                "Caller technical K is context only and is not COI."] if population_ok else
               [f"Fewer than {min_population} samples with enough informative loci: MOIRE was not run, because "
                "its allele frequencies and COI are not identified in so few samples. Monoclonal/polyclonal calls "
                "come from the PfVarX-Call mixture test (one genotype against several, per sample; basis "
                "CALLER_MIXTURE_LRT). The COI count and FWS are withheld."])
    write_manifest(ctx, "DA001", terminal.value, reason_codes=[] if assessed else [Reason.INSUFFICIENT_INFORMATIVE_SITES.value], requested=ctx.samples, eligible=eligible,
                   excluded=excluded, analysed=assessed, outputs=outputs if assessed else [], caveats=caveats,
                   method={"estimator": "MOIRE" if population_ok else "CALLER_MIXTURE_LRT", "credibility": cred,
                           "sites": len(sites), "population_min_samples": min_population,
                           "moire_run": population_ok, "presence_floor": floor, "thin_bp": ctx.params["da001_thin_bp"]},
                   results_summary={"classification_counts": counts})
