"""DA005 population structure: haploid PCA and Hudson FST (scikit-allel).

Biallelic SNV sites observed in every analysed sample, dominant-allele encoding with
ploidy set to 1 explicitly; linkage pruning with ``locate_unlinked`` (scikit-allel
defaults window 100, step 20, r^2 0.1, iterated once); FST only between groups defined
in the metadata before analysis, with block-jackknife standard errors.
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from pfvarx.configuration.loader import frozen_param
from pfvarx.constants import Reason, Terminal
from pfvarx.downstream.common import PANEL_CAVEAT, Context, decisions, write_manifest
from pfvarx.downstream.diversity import dominant
from pfvarx.io import write_json, write_tsv

# Layer 5 memory: the cohort evidence matrix is loaded for this analysis only when it reads it
# (dominant alleles (through diversity) at every cohort site). tests/unit/test_matrix_declarations.py checks the declaration against the code.
USES_COHORT_MATRIX = True

def run(ctx: Context, threads: int = 1) -> None:
    import allel

    out = ctx.out("DA005")
    dec = decisions(ctx)["analyses"]["DA005"]
    eligible = sorted({s for c in dec["scales"].values() for s in c["eligible"]})
    excluded = {"INELIGIBLE": sorted(set(ctx.samples) - set(eligible))}
    col = ctx.eff["config"]["analyses"].get("population_column") or "population"
    rows = []
    for s in ctx.sites:
        if s["type"] != "SNV":
            continue
        g = [dominant(ctx, s["site_id"], x) for x in eligible]
        if any(v is None for v in g) or len(set(g)) < 2:
            continue
        rows.append(g)
    p = ctx.params
    min_sites = int(frozen_param(p, "da005_min_sites"))
    if len(rows) < min_sites or len(eligible) < int(frozen_param(p, "da005_min_samples")):
        write_manifest(ctx, "DA005", Terminal.SKIPPED_INELIGIBLE.value,
                       reason_codes=[Reason.INSUFFICIENT_INFORMATIVE_SITES.value], requested=ctx.samples,
                       eligible=eligible, excluded=excluded, analysed=[], outputs=[], caveats=[],
                       method={"tool": "scikit-allel"})
        return
    gn = np.array(rows, dtype="i1")
    loc = allel.locate_unlinked(gn, size=int(p["da005_ld_window"]), step=int(p["da005_ld_step"]),
                                threshold=float(p["da005_ld_threshold"]))
    gpr = gn[loc] if loc.sum() >= min_sites else gn
    ncomp = min(int(frozen_param(p, "da005_pca_components")), len(eligible) - 1, gpr.shape[0])
    coords, model = allel.pca(gpr, n_components=ncomp, scaler="patterson", ploidy=1)
    pcs = [{"sample_id": s, "group": ctx.metadata.get(s, {}).get(col, ""),
            **{f"PC{i + 1}": float(coords[j, i]) for i in range(ncomp)}} for j, s in enumerate(eligible)]
    write_tsv(out / "pca.tsv", pcs, list(pcs[0]))
    groups: dict[str, list[int]] = {}
    for j, s in enumerate(eligible):
        group_name = ctx.metadata.get(s, {}).get(col, "")
        if group_name:
            groups.setdefault(group_name, []).append(j)
    fst_rows: list[dict[str, Any]] = []
    names = sorted(groups)
    minimum = int(p["da005_min_group_size"])
    for i in range(len(names)):
        for k in range(i + 1, len(names)):
            a, b = groups[names[i]], groups[names[k]]
            if len(a) < minimum or len(b) < minimum:
                fst_rows.append({"group_1": names[i], "group_2": names[k], "status": "GROUP_TOO_SMALL"})
                continue
            ac1 = np.stack([len(a) - gn[:, a].sum(1), gn[:, a].sum(1)], axis=1)
            ac2 = np.stack([len(b) - gn[:, b].sum(1), gn[:, b].sum(1)], axis=1)
            blen = max(1, gn.shape[0] // 10)
            fst, se, _, _ = allel.average_hudson_fst(ac1, ac2, blen=blen)
            fst_rows.append({"group_1": names[i], "group_2": names[k], "status": "COMPUTED", "hudson_fst": float(fst),
                             "jackknife_se": float(se), "block_length_sites": blen, "sites": int(gn.shape[0])})
    write_tsv(out / "fst.tsv", fst_rows, ["group_1", "group_2", "status", "hudson_fst", "jackknife_se",
                                          "block_length_sites", "sites"])
    write_json(out / "pca_model.json", {"explained_variance_ratio": [float(x) for x in model.explained_variance_ratio_],
                                        "sites_before_pruning": int(gn.shape[0]), "sites_after_pruning": int(gpr.shape[0])})
    write_manifest(ctx, "DA005", (Terminal.COMPLETED_WITH_EXCLUSIONS if excluded["INELIGIBLE"] else Terminal.COMPLETED).value,
                   reason_codes=[], requested=ctx.samples, eligible=eligible, excluded=excluded, analysed=eligible,
                   outputs=[out / "pca.tsv", out / "fst.tsv", out / "pca_model.json"],
                   caveats=["Polyclonal samples are represented by their dominant allele.",
                            "PCA axes describe the analysed samples only; filtering in the dashboard does not refit."]
                   + ([PANEL_CAVEAT] if ctx.design.targeted else []),
                   method={"tool": "scikit-allel", "ploidy": 1, "scaler": "patterson",
                           "ld": {"size": p["da005_ld_window"], "step": p["da005_ld_step"], "r2": p["da005_ld_threshold"]}},
                   results_summary=json.loads(json.dumps({"fst_pairs": len(fst_rows)})))
