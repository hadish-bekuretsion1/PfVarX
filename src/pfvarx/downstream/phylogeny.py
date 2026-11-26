"""DA006 conventional single-lineage phylogenetics (IQ-TREE 2).

Phase 1 (eligibility checkpoint) admitted only MONOCLONAL samples passing QC. Here the
remaining runtime gates are re-checked on the actual alignment: per-sample and per-site
missingness, parsimony-informative sites and the four-sample floor. A tree is built
only if every gate passes; otherwise the terminal state is SKIPPED_INELIGIBLE and no
tree file exists. Variable-site alignments use ascertainment-corrected (+ASC) models.
"""

from __future__ import annotations

import shutil
import subprocess

from pfvarx.constants import Reason, Terminal
from pfvarx.downstream.common import PANEL_CAVEAT, Context, decisions, write_manifest
from pfvarx.io import write_json

# Layer 5 memory: the cohort evidence matrix is loaded for this analysis only when it reads it
# (dominant alleles at every informative site). tests/unit/test_matrix_declarations.py checks the declaration against the code.
USES_COHORT_MATRIX = True

def build_alignment(ctx: Context, samples: list[str]) -> tuple[dict[str, list[str]], list[str]]:
    seqs: dict[str, list[str]] = {s: [] for s in samples}
    used = []
    for s in ctx.sites:
        if s["type"] != "SNV":
            continue
        col = []
        for x in samples:
            c = ctx.matrix[s["site_id"]].get(x)
            if not c or c["wsaf"] is None:
                col.append("-")
            else:
                col.append(s["alt"] if c["wsaf"] >= 0.5 else s["ref"])
        for x, b in zip(samples, col):
            seqs[x].append(b)
        used.append(s["site_id"])
    return seqs, used


def filter_alignment(seqs: dict[str, list[str]], max_sample_missing: float, max_site_missing: float
                     ) -> tuple[dict[str, list[str]], list[str], list[str]]:
    excl = [s for s, v in seqs.items() if v and v.count("-") / len(v) > max_sample_missing]
    keep = {s: v for s, v in seqs.items() if s not in excl}
    if not keep:
        return {}, excl, []
    n = len(next(iter(keep.values())))
    cols = []
    for i in range(n):
        column = [v[i] for v in keep.values()]
        if column.count("-") / len(column) > max_site_missing:
            continue
        if len({b for b in column if b != "-"}) < 2:
            continue  # invariant after filtering: not allowed under +ASC
        cols.append(i)
    return {s: [v[i] for i in cols] for s, v in keep.items()}, excl, [str(i) for i in cols]


def parsimony_informative(seqs: dict[str, list[str]]) -> int:
    if not seqs:
        return 0
    n = len(next(iter(seqs.values())))
    count = 0
    for i in range(n):
        tally: dict[str, int] = {}
        for v in seqs.values():
            if v[i] != "-":
                tally[v[i]] = tally.get(v[i], 0) + 1
        if sum(1 for k in tally.values() if k >= 2) >= 2:
            count += 1
    return count


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA006")
    for stale in out.glob("tree*"):
        stale.unlink()
    dec = decisions(ctx)["analyses"]["DA006"]
    eligible = sorted({s for c in dec["scales"].values() for s in c["eligible"]})
    excluded: dict[str, list[str]] = {}
    for c in dec["scales"].values():
        for s in c["excluded"]:
            excluded.setdefault("PHASE1_EXCLUSION", []).append(s)
    p = ctx.params
    seqs, _ = build_alignment(ctx, eligible)
    kept, miss_excl, cols = filter_alignment(seqs, float(p["da006_max_sample_missing"]), float(p["da006_max_site_missing"]))
    if miss_excl:
        excluded[Reason.MISSINGNESS_ABOVE_MAXIMUM.value] = miss_excl
    pis = parsimony_informative(kept)
    gates = {"monoclonal_after_missingness": len(kept), "minimum_samples": int(p["da006_min_monoclonal"]),
             "parsimony_informative_sites": pis, "minimum_informative_sites": int(p["da006_min_informative_sites"]),
             "alignment_columns": len(cols)}
    write_json(out / "runtime_gates.json", gates)
    failing = []
    if len(kept) < gates["minimum_samples"]:
        failing.append(Reason.MONOCLONAL_COUNT_BELOW_MINIMUM.value)
    if pis < gates["minimum_informative_sites"]:
        failing.append(Reason.INSUFFICIENT_INFORMATIVE_SITES.value)
    if failing:
        write_manifest(ctx, "DA006", Terminal.SKIPPED_INELIGIBLE.value, reason_codes=failing, requested=ctx.samples,
                       eligible=eligible, excluded=excluded, analysed=[], outputs=[out / "runtime_gates.json"],
                       caveats=["No tree was produced; no placeholder tree exists."],
                       method={"tool": "IQ-TREE 2"}, results_summary=gates)
        return
    aln = out / "alignment.fasta"
    aln.write_text("".join(f">{s}\n{''.join(v)}\n" for s, v in kept.items()))
    exe = shutil.which("iqtree2") or shutil.which("iqtree")
    if exe is None:
        raise SystemExit("METHOD_NOT_INSTALLED: iqtree2")
    cmd = [exe, "-s", str(aln), "-m", "MFP+ASC", "-B", str(p["da006_ufboot"]), "--alrt", str(p["da006_alrt"]),
           "--seed", str(ctx.eff["config"]["random_seed"]), "-T", str(max(1, threads)), "--prefix",
           str(out / "tree"), "-redo", "--quiet"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    (out / "iqtree_command.txt").write_text(" ".join(cmd) + "\n" + r.stdout[-3000:] + r.stderr[-3000:])
    if r.returncode != 0 or not (out / "tree.treefile").is_file():
        raise SystemExit(f"IQ-TREE failed with exit code {r.returncode}")
    # IQ-TREE working files no step reads and that add nothing to the kept tree, report and log:
    # the resume checkpoint, the starting BIONJ tree, the model-selection cache and the bootstrap splits.
    for name in ("tree.ckp.gz", "tree.bionj", "tree.model.gz", "tree.splits.nex", "tree.uniqueseq.phy"):
        (out / name).unlink(missing_ok=True)
    model = ""
    iq = out / "tree.iqtree"
    if iq.is_file():
        for line in iq.read_text().splitlines():
            if line.startswith("Best-fit model according to"):
                model = line.split(":", 1)[1].strip()
    terminal = Terminal.COMPLETED_WITH_EXCLUSIONS if any(excluded.values()) else Terminal.COMPLETED
    write_manifest(ctx, "DA006", terminal.value, reason_codes=[], requested=ctx.samples, eligible=eligible,
                   excluded=excluded, analysed=sorted(kept), outputs=[aln, out / "tree.treefile", out / "runtime_gates.json"],
                   caveats=[f"Branch support: UFBoot >= {p['da006_ufboot_support']} and SH-aLRT >= "
                            f"{p['da006_alrt_support']} are the developers' interpretation guides.",
                            "Only MONOCLONAL samples are included; polyclonal and indeterminate samples are excluded."]
                   + ([PANEL_CAVEAT] if ctx.design.targeted else []),
                   method={"tool": "IQ-TREE 2", "model_selection": "ModelFinder restricted to +ASC", "best_model": model,
                           "ufboot": p["da006_ufboot"], "alrt": p["da006_alrt"],
                           "ufboot_support": p["da006_ufboot_support"], "alrt_support": p["da006_alrt_support"]}, results_summary=gates)
