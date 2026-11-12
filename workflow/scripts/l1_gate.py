"""Collect Layer 1 metrics for one sample and apply the frozen QC policy (checkpoint)."""

import gzip
import json
from pathlib import Path

from pfvarx.quality import metrics as m
from pfvarx.quality.policy import evaluate_sample, load_policy, stage_permitted

sm = snakemake  # noqa: F821  (injected by Snakemake)
eff = json.loads(Path(sm.params.frozen).read_text())
design = eff["config"]["assay"]["design"]
min_depth = int(eff["parameters"]["layer2.min_usable_depth"]["value"])
inp = sm.input
values: dict = {}
values.update(m.flagstat(Path(inp.flagstat)))
if "fastp" in inp.keys():
    before = after = 0
    lengths = []
    q30 = []
    for f in inp.fastp:
        d = m.fastp(Path(f))
        before += d["reads_before"]
        after += d["reads_after"]
        lengths.append(d["mean_read_length_after"])
        q30.append(d["q30_fraction_after"])
    values.update({"reads_before": before, "reads_after": after,
                   "reads_retained_fraction": after / before if before else None,
                   "mean_read_length_after": sum(lengths) / len(lengths), "q30_fraction_after": min(q30)})
if "markdup" in inp.keys():
    values.update(m.picard_duplicates(Path(inp.markdup)))
if "mosdepth" in inp.keys():
    values.update(m.mosdepth_region_dist(Path(inp.mosdepth), min_depth))
if "targets" in inp.keys():
    text = gzip.open(inp.targets, "rt").read()
    on = m.samtools_count(Path(inp.ontarget).read_text())
    mapped = int(values.get("mapped_fraction", 0) * values.get("primary_reads", 0))
    values.update(m.target_metrics(text, min_depth, on, mapped))
if "umi" in inp.keys():
    values.update(m.umi_family_sizes(Path(inp.umi)))
policy = load_policy(Path(sm.params.root), eff["quality"]["policy_id"])
qc = evaluate_sample(policy, values, design, stages={"G1", "G2", "G3"})
stages = ["L1_OUTPUT", "L2", "L3", "L4", "DA002", "DA003", "DA001", "DA004", "DA005", "DA006", "DA007"]
qc["permitted"] = {st: (True if st == "L1_OUTPUT" and qc["verdict"] != "FAIL" else stage_permitted(qc, st))
                   for st in stages}
qc["sample_id"] = sm.wildcards.sample
qc["metrics"] = values
Path(sm.output[0]).write_text(json.dumps(qc, indent=2, default=str))
