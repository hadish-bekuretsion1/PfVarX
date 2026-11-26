"""DA009 pfhrp2/pfhrp3 joint deletion: single and dual deletion in every assay design.

HRP2-based rapid diagnostic tests fail when a parasite has lost *pfhrp2*; loss of *pfhrp3*
matters because HRP3 cross-reacts with the same test. What a surveillance programme acts on is
therefore the **joint** state of the two genes, not either gene alone, and how that joint state
changes over place and time. This analysis reports both.

Two evidence models give each gene a per-sample state, chosen by the run's assay design:

* ``DEPTH_RATIO`` (``wgs``): the same measurement DA008 makes — the gene's read depth against
  core-genome windows of the same sample with similar GC content, with a bootstrap interval,
  the flanking depth and spanning read pairs. DA008 and DA009 call WGS deletions identically.
* ``AMPLICON_DROPOUT`` (``targeted_umi``, ``targeted_no_umi``): amplicon depth does not track
  copy number, so copy number cannot be estimated. What a panel shows is whether the
  pfhrp2/pfhrp3 amplicons produced reads while the panel's other amplicons (the controls) did.
  The gene's targets are compared with the control targets of the same sample, and the result
  is a presence/absence statement, not a copy number.

The two models answer related but different questions and fail differently: an amplicon can
drop out because the gene is deleted, because a primer site carries a mutation, or because that
PCR simply failed. Every row records its ``evidence_class`` for that reason, prevalence is
reported both pooled and stratified by evidence class, and the caveats say so. A run has one
assay design, so within a single run the two tables coincide; they differ only when results
from several runs are combined, which is a deliberate act outside this pipeline.

Joint state per sample (``INDETERMINATE`` whenever either gene is UNKNOWN; a gene is never
assumed present):

| pfhrp2 | pfhrp3 | joint state |
| --- | --- | --- |
| DELETED | DELETED | ``DOUBLE_DELETION`` |
| DELETED | PRESENT | ``PFHRP2_ONLY`` |
| PRESENT | DELETED | ``PFHRP3_ONLY`` |
| PRESENT | PRESENT | ``NO_DELETION`` |

Outputs: per-gene calls, per-sample joint calls, joint prevalence (pooled and by evidence
class), the co-occurrence of the two deletions against independence, and the joint prevalence
by country and month. Thresholds are UNCALIBRATED_PROVISIONAL and the copy-number catalogue
must be MAINTAINER_VERIFIED, as for DA008.
"""

from __future__ import annotations

import math
import statistics
from statistics import NormalDist
from typing import Any

import pysam

from pfvarx.configuration.loader import frozen_param
from pfvarx.constants import HrpEvidence, JointDeletion, Reason, Terminal
from pfvarx.downstream.cnv import (
    GcModel,
    call_locus,
    catalog_path,
    depth_profile,
    gene_location,
)
from pfvarx.downstream.cnv import measure as depth_measure
from pfvarx.downstream.common import PANEL_CAVEAT, Context, decisions, wilson, write_manifest
from pfvarx.downstream.resistance import load_catalog
from pfvarx.eligibility.territory import read_intervals
from pfvarx.io import write_json, write_tsv

MEASURES = ("pfhrp2_deleted", "pfhrp3_deleted", "double_deletion", "any_deletion")
DELETED, PRESENT, UNKNOWN = "DELETED", "PRESENT", "UNKNOWN"


def z_score(level: float) -> float:
    """Two-sided normal quantile of a confidence level (0.95 -> 1.96)."""
    return float(NormalDist().inv_cdf(1 - (1 - level) / 2))


def contract_loci(contract: dict[str, Any]) -> list[str]:
    """The catalogue locus ids this analysis reports, from the maintainer registry."""
    loci = [str(x) for x in ((contract.get("needs") or {}).get("deletion_loci") or [])]
    if len(loci) != 2:
        raise ValueError("DA009 needs exactly two deletion_loci in config/registry/downstream_analyses.yaml")
    return loci


def gene_rows(catalogue: list[dict[str, str]], loci: list[str]) -> list[dict[str, str]]:
    """Catalogue rows for *loci*, in that order; a missing or non-deletion locus is an error."""
    by_id = {row["locus_id"]: row for row in catalogue}
    out = []
    for locus in loci:
        row = by_id.get(locus)
        if row is None or row.get("event") != "DELETION":
            raise ValueError(f"copy-number catalogue has no DELETION locus {locus!r} for DA009")
        out.append(row)
    return out


def panel_target_spans(ctx: Context) -> list[tuple[str, int, int]]:
    """Merged target intervals of the run's panel (empty for whole-genome sequencing)."""
    panel = ctx.plan.get("panel") or {}
    bed = panel.get("targets_bed") or ctx.eff["config"].get("assay", {}).get("targets_bed")
    if not bed:
        return []
    return [(chrom, s, e) for chrom, spans in read_intervals(ctx.root / bed).items() for s, e in spans]


def split_targets(spans: list[tuple[str, int, int]], places: dict[str, tuple[str, int, int] | None]
                  ) -> tuple[dict[str, list[tuple[str, int, int]]], list[tuple[str, int, int]]]:
    """Panel targets that overlap each gene, and the remaining targets, which act as controls."""
    gene_targets: dict[str, list[tuple[str, int, int]]] = {locus: [] for locus in places}
    controls: list[tuple[str, int, int]] = []
    for span in spans:
        chrom, start, end = span
        hit = [locus for locus, loc in places.items()
               if loc is not None and loc[0] == chrom and start < loc[2] and loc[1] < end]
        if hit:
            for locus in hit:
                gene_targets[locus].append(span)
        else:
            controls.append(span)
    return gene_targets, controls


def median_target_depth(bam: pysam.AlignmentFile, spans: list[tuple[str, int, int]], min_mapq: int) -> float | None:
    """Median per-base depth over the given intervals (None when none can be read)."""
    depths: list[int] = []
    for chrom, start, end in spans:
        if chrom in bam.references:
            depths += depth_profile(bam, chrom, start, end, min_mapq)
    return float(statistics.median(depths)) if depths else None


def dropout_call(gene_depth: float | None, control_depth: float | None, controls: int,
                 p: dict[str, Any]) -> tuple[str, str, float | None]:
    """Presence/absence of a gene's amplicons against the panel's control amplicons.

    Returns (state, reason, ratio). A low control depth means the assay did not work for that
    sample, which is UNKNOWN, never a deletion.
    """
    min_controls = int(frozen_param(p, "da009_min_control_amplicons"))
    if controls < min_controls:
        return UNKNOWN, "TOO_FEW_CONTROL_AMPLICONS", None
    if control_depth is None or control_depth < float(frozen_param(p, "da009_min_control_depth")):
        return UNKNOWN, "CONTROL_AMPLICON_DEPTH_LOW", None
    if gene_depth is None:
        return UNKNOWN, "GENE_TARGET_NOT_IN_PANEL", None
    ratio = gene_depth / control_depth
    if ratio <= float(frozen_param(p, "da009_dropout_max_ratio")):
        return DELETED, "", ratio
    if ratio >= float(frozen_param(p, "da009_present_min_ratio")):
        return PRESENT, "", ratio
    return UNKNOWN, "AMBIGUOUS_COVERAGE", ratio


def joint_state(first: str, second: str) -> str:
    """Joint pfhrp2/pfhrp3 state; either gene UNKNOWN gives INDETERMINATE."""
    if UNKNOWN in (first, second):
        return JointDeletion.INDETERMINATE.value
    if first == DELETED and second == DELETED:
        return JointDeletion.DOUBLE_DELETION.value
    if first == DELETED:
        return JointDeletion.PFHRP2_ONLY.value
    if second == DELETED:
        return JointDeletion.PFHRP3_ONLY.value
    return JointDeletion.NO_DELETION.value


def carries(state: str, measure: str) -> bool | None:
    """Whether a joint state counts towards *measure*; None when it cannot be judged."""
    if state == JointDeletion.INDETERMINATE.value:
        return None
    deleted2 = state in (JointDeletion.DOUBLE_DELETION.value, JointDeletion.PFHRP2_ONLY.value)
    deleted3 = state in (JointDeletion.DOUBLE_DELETION.value, JointDeletion.PFHRP3_ONLY.value)
    return {"pfhrp2_deleted": deleted2, "pfhrp3_deleted": deleted3,
            "double_deletion": deleted2 and deleted3, "any_deletion": deleted2 or deleted3}[measure]


def prevalence(rows: list[dict[str, Any]], level: float, by_evidence: bool) -> list[dict[str, Any]]:
    """Joint prevalence per group, over samples whose joint state could be judged.

    With *by_evidence* each group is also split by evidence class; otherwise the evidence
    classes present are pooled and listed in the row, because an amplicon-dropout call and a
    depth-ratio call are not the same measurement.
    """
    z = z_score(level)
    groups = sorted({r["population"] for r in rows if r["population"]}) + ["ALL"]
    classes = sorted({r["evidence_class"] for r in rows}) if by_evidence else [None]
    out: list[dict[str, Any]] = []
    for group in groups:
        for evidence in classes:
            sel = [r for r in rows if (group == "ALL" or r["population"] == group)
                   and (evidence is None or r["evidence_class"] == evidence)]
            judged = [r for r in sel if r["joint_state"] != JointDeletion.INDETERMINATE.value]
            present = sorted({r["evidence_class"] for r in judged})
            for measure in MEASURES:
                k = sum(1 for r in judged if carries(r["joint_state"], measure))
                lo, hi = wilson(k, len(judged), z)
                row = {"group": group, "measure": measure, "judged": len(judged),
                       "carriers": k, "indeterminate": len(sel) - len(judged),
                       "prevalence": None if not judged else k / len(judged),
                       "ci_lower": lo, "ci_upper": hi, "interval_level": level,
                       "status": "COMPUTED" if judged else "NO_JUDGED_SAMPLES"}
                if evidence is None:
                    row["evidence_classes"] = ",".join(present)
                    row["pooled_evidence_classes"] = len(present)
                else:
                    row["evidence_class"] = evidence
                out.append(row)
    return out


def cooccurrence(rows: list[dict[str, Any]], level: float) -> list[dict[str, Any]]:
    """Do the two deletions occur together more often than independence predicts?

    Fisher's exact test on the 2x2 table of pfhrp2 by pfhrp3 deletion, with the double-deletion
    count expected if the two were independent. This describes the cohort's joint pattern; it is
    not evidence of selection, and a cohort is rarely a random sample of infections.
    """
    out: list[dict[str, Any]] = []
    for group in sorted({r["population"] for r in rows if r["population"]}) + ["ALL"]:
        judged = [r for r in rows if (group == "ALL" or r["population"] == group)
                  and r["joint_state"] != JointDeletion.INDETERMINATE.value]
        n = len(judged)
        both = sum(1 for r in judged if r["joint_state"] == JointDeletion.DOUBLE_DELETION.value)
        only2 = sum(1 for r in judged if r["joint_state"] == JointDeletion.PFHRP2_ONLY.value)
        only3 = sum(1 for r in judged if r["joint_state"] == JointDeletion.PFHRP3_ONLY.value)
        neither = sum(1 for r in judged if r["joint_state"] == JointDeletion.NO_DELETION.value)
        row: dict[str, Any] = {"group": group, "judged": n, "both_deleted": both, "pfhrp2_only": only2,
                               "pfhrp3_only": only3, "neither_deleted": neither,
                               "expected_both_if_independent": None, "odds_ratio": None, "p_value": None,
                               "interval_level": level, "status": "NO_JUDGED_SAMPLES" if not n else "COMPUTED"}
        if n:
            row["expected_both_if_independent"] = (both + only2) * (both + only3) / n
        if not n or (both + only2) in (0, n) or (both + only3) in (0, n):
            # One margin is fixed: no association can be estimated from this table.
            row["status"] = "NO_VARIATION" if n else "NO_JUDGED_SAMPLES"
            out.append(row)
            continue
        from scipy import stats

        # rows: pfhrp2 deleted or not; columns: pfhrp3 deleted or not
        result = stats.fisher_exact([[both, only2], [only3, neither]])
        odds = float(result.statistic)
        row["odds_ratio"] = None if math.isnan(odds) else odds
        row["p_value"] = float(result.pvalue)
        out.append(row)
    return out


def trend(rows: list[dict[str, Any]], min_per_unit: int, level: float) -> list[dict[str, Any]]:
    """Joint prevalence by country and month: how the joint state stands over place and time.

    Units with fewer than *min_per_unit* judged samples report their counts with the prevalence
    withheld, so a single infection never becomes a percentage.
    """
    z = z_score(level)
    out: list[dict[str, Any]] = []
    units = sorted({(r["country"], r["month"]) for r in rows if r["country"] and r["month"]})
    for country, month in units:
        judged = [r for r in rows if r["country"] == country and r["month"] == month
                  and r["joint_state"] != JointDeletion.INDETERMINATE.value]
        for measure in ("double_deletion", "any_deletion"):
            k = sum(1 for r in judged if carries(r["joint_state"], measure))
            enough = len(judged) >= min_per_unit
            lo, hi = wilson(k, len(judged), z) if enough else (None, None)
            out.append({"country": country, "month": month, "measure": measure, "judged": len(judged),
                        "carriers": k, "prevalence": (k / len(judged)) if enough and judged else None,
                        "ci_lower": lo, "ci_upper": hi, "interval_level": level,
                        "evidence_classes": ",".join(sorted({r["evidence_class"] for r in judged})),
                        "status": "COMPUTED" if enough else "INSUFFICIENT_SAMPLES"})
    return out


def _month(value: str) -> str:
    """Collection month (YYYY-MM) of an ISO collection date, else empty."""
    text = (value or "").strip()
    return text[:7] if len(text) >= 7 and text[:4].isdigit() and text[5:7].isdigit() else ""


def fmt(x: float | None) -> str:
    return "" if x is None else f"{x:.4f}"


def gene_calls(ctx: Context, eligible: list[str], genes: list[dict[str, str]],
               places: dict[str, tuple[str, int, int] | None], evidence: str, threads: int = 1
               ) -> list[dict[str, Any]]:
    """Per-sample, per-gene deletion state with the evidence the model used."""
    from pfvarx.quality.metrics import mosdepth_region_dist

    p = ctx.params
    min_mapq = int(ctx.eff["parameters"]["layer2.min_mapping_quality"]["value"])
    min_depth = int(ctx.eff["parameters"]["layer2.min_usable_depth"]["value"])
    seed = int(ctx.eff["config"].get("random_seed") or 1)
    pop_col = ctx.eff["config"]["analyses"].get("population_column") or "population"
    fasta_path = ctx.root / ctx.plan["reference"]["parasite_fasta"]
    targets = panel_target_spans(ctx)
    gene_targets, controls = split_targets(targets, places)
    core_bed = ctx.root / ctx.eff["config"]["reference"]["core_mask_bed"]
    core_spans = [(c[0], int(c[1]), int(c[2])) for c in
                  (ln.split("\t") for ln in core_bed.read_text().splitlines()
                   if ln.strip() and not ln.startswith(("#", "track")))]
    min_core = float(p["da008_min_core_depth"])
    calls: list[dict[str, Any]] = []
    for sid in eligible:
        meta = ctx.metadata.get(sid, {})
        dist = ctx.results / "qc" / "samples" / sid / "core.mosdepth.region.dist.txt"
        core = float(mosdepth_region_dist(dist, min_depth).get("median_depth_core", 0)) if dist.is_file() else 0.0
        with pysam.AlignmentFile(str(ctx.results / "alignments" / f"{sid}.bam")) as bam:
            control_depth = (median_target_depth(bam, controls, min_mapq)
                             if evidence == HrpEvidence.AMPLICON_DROPOUT else None)
            model = None
            fasta = None
            if evidence == HrpEvidence.DEPTH_RATIO and core >= min_core:
                fasta = pysam.FastaFile(str(fasta_path))
                model = GcModel(bam, fasta, core_spans, int(p["da008_window_bp"]), int(p["da008_gc_windows"]),
                                min_mapq, seed, int(frozen_param(p, "da008_gc_min_windows_per_bin")))
            try:
                for i, gene in enumerate(genes):
                    locus, loc = gene["locus_id"], places[gene["locus_id"]]
                    row: dict[str, Any] = {
                        "sample_id": sid, "locus_id": locus, "gene_id": gene["gene_id"],
                        "gene": gene["gene_name"], "evidence_class": evidence,
                        "design": ctx.design.value,
                        "location": "" if loc is None else f"{loc[0]}:{loc[1] + 1}-{loc[2]}",
                        "population": meta.get(pop_col, ""), "country": meta.get("country", ""),
                        "month": _month(meta.get("collection_date", "")),
                    }
                    if loc is None:
                        state, reason, extra = UNKNOWN, "GENE_NOT_IN_GFF", {}
                    elif evidence == HrpEvidence.AMPLICON_DROPOUT:
                        gene_depth = median_target_depth(bam, gene_targets.get(locus, []), min_mapq)
                        state, reason, ratio = dropout_call(gene_depth, control_depth, len(controls), p)
                        extra = {"gene_median_depth": fmt(gene_depth), "control_median_depth": fmt(control_depth),
                                 "control_amplicons": len(controls), "coverage_ratio": fmt(ratio),
                                 "gene_targets": len(gene_targets.get(locus, []))}
                    elif core < min_core:
                        state, reason, extra = UNKNOWN, "LOW_CORE_DEPTH", {"core_median_depth": fmt(core)}
                    else:
                        m = depth_measure(bam, loc, core, p, min_mapq, seed + i, fasta, model)
                        state, reason = call_locus("DELETION", m["depth_ratio"], m["ratio_lo"],
                                                   m["flank_ratio"], p)
                        extra = {"gene_median_depth": fmt(m["gene_depth"]), "core_median_depth": fmt(core),
                                 "depth_ratio": fmt(m["depth_ratio"]), "ratio_lo": fmt(m["ratio_lo"]),
                                 "ratio_hi": fmt(m["ratio_hi"]), "flank_ratio": fmt(m["flank_ratio"]),
                                 "spanning_pairs": m["spanning_pairs"]}
                    calls.append({**row, "state": state, "reason": reason, **extra})
            finally:
                if fasta is not None:
                    fasta.close()
    return calls


def joint_rows(calls: list[dict[str, Any]], genes: list[dict[str, str]]) -> list[dict[str, Any]]:
    """One row per sample: both gene states and the joint state."""
    first, second = genes[0]["locus_id"], genes[1]["locus_id"]
    by_sample: dict[str, dict[str, dict[str, Any]]] = {}
    for call in calls:
        by_sample.setdefault(call["sample_id"], {})[call["locus_id"]] = call
    rows = []
    for sid, per_locus in sorted(by_sample.items()):
        a, b = per_locus.get(first, {}), per_locus.get(second, {})
        state = joint_state(a.get("state", UNKNOWN), b.get("state", UNKNOWN))
        reasons = [r for r in (a.get("reason"), b.get("reason")) if r]
        rows.append({"sample_id": sid, "joint_state": state,
                     f"{genes[0]['gene_name']}_state": a.get("state", UNKNOWN),
                     f"{genes[1]['gene_name']}_state": b.get("state", UNKNOWN),
                     "evidence_class": a.get("evidence_class", ""), "design": a.get("design", ""),
                     "indeterminate_reason": ",".join(dict.fromkeys(reasons)),
                     "population": a.get("population", ""), "country": a.get("country", ""),
                     "month": a.get("month", "")})
    return rows


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA009")
    registry_contract = _contract(ctx)
    dec = decisions(ctx)["analyses"]["DA009"]
    eligible = sorted({s for c in dec["scales"].values() for s in c["eligible"]})
    excluded = {"INELIGIBLE": sorted(set(ctx.samples) - set(eligible))}
    cat_path = catalog_path(ctx)
    meta, catalogue = load_catalog(cat_path)
    status = meta.get("status", "UNKNOWN")
    p = ctx.params
    level = float(frozen_param(p, "da009_interval_level"))

    def stop(terminal: str, reasons: list[str], caveats: list[str]) -> None:
        write_manifest(ctx, "DA009", terminal, reason_codes=reasons, requested=ctx.samples, eligible=eligible,
                       excluded=excluded, analysed=[], outputs=[], caveats=caveats,
                       method={"catalog": str(cat_path), "catalog_meta": meta,
                               "evidence_class": _evidence(ctx).value})
    if status not in ("MAINTAINER_VERIFIED", "TEST_FIXTURE"):
        stop(Terminal.BLOCKED_MISSING_PREREQUISITE.value, [Reason.CATALOG_NOT_VERIFIED.value],
             [f"catalogue {cat_path.name} has status {status}; verify the pfhrp2 and pfhrp3 gene "
              "identifiers against the literature before relying on DA009"])
        return
    try:
        genes = gene_rows(catalogue, contract_loci(registry_contract))
    except ValueError as exc:
        stop(Terminal.BLOCKED_MISSING_PREREQUISITE.value, [Reason.CATALOG_MALFORMED.value],
             [f"{exc}. DA009 reports the two loci named in the analysis registry; add them to the "
              "copy-number catalogue with event DELETION, or do not request DA009."])
        return
    gff = ctx.root / ctx.eff["config"]["reference"]["gff"] if ctx.eff["config"]["reference"].get("gff") else None
    places = {g["locus_id"]: (gene_location(gff, g["gene_id"]) if gff and gff.is_file() else None) for g in genes}
    evidence = _evidence(ctx)
    if evidence == HrpEvidence.AMPLICON_DROPOUT:
        gene_targets, controls = split_targets(panel_target_spans(ctx), places)
        covered = [locus for locus, spans in gene_targets.items() if spans]
        min_controls = int(frozen_param(p, "da009_min_control_amplicons"))
        if not covered:
            stop(Terminal.SKIPPED_INELIGIBLE.value, [Reason.PANEL_LACKS_HRP_TARGETS.value],
                 ["The panel has no target region overlapping pfhrp2 or pfhrp3, so amplicon dropout cannot "
                  "be measured. Use a panel that carries both genes, or whole-genome sequencing."])
            return
        if len(controls) < min_controls:
            stop(Terminal.SKIPPED_INELIGIBLE.value, [Reason.CONTROL_AMPLICONS_INSUFFICIENT.value],
                 [f"The panel has {len(controls)} control target(s) outside pfhrp2/pfhrp3 and at least "
                  f"{min_controls} are needed: without controls a missing amplicon cannot be told from a "
                  "failed PCR."])
            return
    calls = gene_calls(ctx, eligible, genes, places, evidence.value, threads)
    joint = joint_rows(calls, genes)
    pooled = prevalence(joint, level, by_evidence=False)
    by_evidence = prevalence(joint, level, by_evidence=True)
    together = cooccurrence(joint, level)
    over_time = trend(joint, int(frozen_param(p, "da009_min_samples_per_unit")), level)
    write_tsv(out / "gene_calls.tsv", calls, _columns(calls, ["sample_id", "locus_id", "state"]))
    write_tsv(out / "joint_calls.tsv", joint, _columns(joint, ["sample_id", "joint_state"]))
    write_tsv(out / "joint_prevalence.tsv", pooled, _columns(pooled, ["group", "measure"]))
    write_tsv(out / "joint_prevalence_by_evidence.tsv", by_evidence, _columns(by_evidence, ["group", "measure"]))
    write_tsv(out / "cooccurrence.tsv", together, _columns(together, ["group"]))
    write_tsv(out / "trend.tsv", over_time, _columns(over_time, ["country", "month", "measure"]))
    write_json(out / "catalog_provenance.json", {"path": str(cat_path), "meta": meta,
                                                 "loci": [g["locus_id"] for g in genes]})
    counts = {s.value: sum(1 for r in joint if r["joint_state"] == s.value) for s in JointDeletion}
    caveats = _caveats(ctx, evidence, status, meta)
    outputs = [out / "gene_calls.tsv", out / "joint_calls.tsv", out / "joint_prevalence.tsv",
               out / "joint_prevalence_by_evidence.tsv", out / "cooccurrence.tsv", out / "trend.tsv",
               out / "catalog_provenance.json"]
    indeterminate = counts[JointDeletion.INDETERMINATE.value]
    write_manifest(ctx, "DA009",
                   (Terminal.COMPLETED_WITH_EXCLUSIONS if excluded["INELIGIBLE"] or indeterminate
                    else Terminal.COMPLETED).value,
                   reason_codes=[], requested=ctx.samples, eligible=eligible, excluded=excluded,
                   analysed=eligible, outputs=outputs, caveats=caveats,
                   method={"catalog": meta, "evidence_class": evidence.value,
                           "genes": {g["gene_name"]: g["gene_id"] for g in genes},
                           "model": ("gene depth against GC-matched core-genome windows (as DA008)"
                                     if evidence == HrpEvidence.DEPTH_RATIO else
                                     "gene amplicon depth against the panel's control amplicons"),
                           "interval": f"Wilson score interval at {level}",
                           "association": "Fisher exact test of pfhrp2 by pfhrp3 deletion",
                           **({"dropout_max_ratio": float(frozen_param(p, "da009_dropout_max_ratio")),
                               "present_min_ratio": float(frozen_param(p, "da009_present_min_ratio")),
                               "min_control_depth": float(frozen_param(p, "da009_min_control_depth")),
                               "min_control_amplicons": int(frozen_param(p, "da009_min_control_amplicons"))}
                              if evidence == HrpEvidence.AMPLICON_DROPOUT else
                              {"deletion_ratio": float(p["da008_deletion_ratio"]),
                               "present_ratio": float(p["da008_present_ratio"]),
                               "min_core_depth": float(p["da008_min_core_depth"])})},
                   results_summary={"gene_calls": len(calls), "samples": len(joint), **counts})


def _contract(ctx: Context) -> dict[str, Any]:
    from pfvarx.eligibility.engine import load_registry

    return dict(load_registry(ctx.root)["analyses"]["DA009"])


def _evidence(ctx: Context) -> HrpEvidence:
    return HrpEvidence.AMPLICON_DROPOUT if ctx.design.targeted else HrpEvidence.DEPTH_RATIO


def _columns(rows: list[dict[str, Any]], fallback: list[str]) -> list[str]:
    seen: dict[str, None] = {}
    for row in rows:
        for key in row:
            seen.setdefault(key)
    return list(seen) or fallback


def _caveats(ctx: Context, evidence: HrpEvidence, status: str, meta: dict[str, str]) -> list[str]:
    caveats = [
        "pfhrp2/pfhrp3 deletion status supports rapid-diagnostic-test surveillance at cohort level; it is "
        "not a diagnosis of an individual infection and does not establish test failure in a patient.",
        "A gene that could not be judged is UNKNOWN and the sample is INDETERMINATE; a gene is never "
        "assumed present.",
        "Prevalence is reported pooled and stratified by evidence class. A depth-ratio call and an "
        "amplicon-dropout call are different measurements, so pooled figures from several runs should be "
        "read with the evidence composition in each row.",
        "The co-occurrence test describes this cohort's joint pattern against independence. It is not "
        "evidence of selection, and a convenience sample is not a random sample of infections.",
        "Thresholds are UNCALIBRATED_PROVISIONAL.",
    ]
    if evidence == HrpEvidence.AMPLICON_DROPOUT:
        caveats[1:1] = [
            "Targeted panel: amplicon depth does not track copy number, so this is a presence/absence "
            "statement about the gene's amplicons, not a copy number. A mutated primer site or a failed "
            "PCR produces the same dropout as a deletion, which the control amplicons reduce but do not "
            "remove.",
            PANEL_CAVEAT,
        ]
    else:
        caveats[1:1] = [
            "Whole-genome sequencing: copy number is estimated from read depth against core-genome windows "
            "of the same sample with similar GC content; pfhrp2 and pfhrp3 lie in subtelomeric, repetitive "
            "regions where depth is less reliable.",
        ]
    if status == "TEST_FIXTURE":
        caveats.insert(0, "TEST FIXTURE CATALOGUE: results are for software testing, not for surveillance.")
    elif meta.get("version"):
        caveats.append(f"Catalogue {meta.get('catalog_id', '')} version {meta['version']}.")
    return caveats
