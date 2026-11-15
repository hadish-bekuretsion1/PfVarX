"""The preflight resource plan must report what the workflow rules actually request.

``planning/plan.py`` summarises each job class for the readiness and dry-run reports; the
rules in ``workflow/rules`` make the real requests. They are kept in two places, so this test
fails when one changes without the other.
"""

import re
from pathlib import Path

from pfvarx.planning.plan import JOB_CLASSES

RULES = Path(__file__).resolve().parents[2] / "workflow" / "rules"

# Job class -> the rules it summarises (the plan reports the largest request among them).
CLASS_RULES = {
    "fastp": ["fastp"],
    "bwa_mem2_align": ["bwa_mem2_lane", "merge_lanes"],
    "markdup": ["mark_duplicates"],
    "fgbio_umi": ["umi_fastq_to_bam", "umi_align_zipper", "umi_group", "umi_consensus", "umi_realign",
                  "umi_preclip_and_trim"],
    "alignment_entry": ["alignment_entry"],
    "fgbio_clip": ["trim_primers"],
    "metrics": ["mosdepth_core", "mosdepth_targets"],
    "pfvarx_call": ["pfvarx_call"],
    "cohort_merge": ["cohort_merge", "cohort_evidence"],
    "snpeff": ["snpeff_build", "snpeff_annotate"],
    "downstream": ["da001_clonality", "da_execute"],
}


def rule_requests() -> dict[str, tuple[int, int]]:
    """(threads, mem_mb) of every rule; Snakemake's defaults are 1 thread and no memory."""
    out: dict[str, tuple[int, int]] = {}
    for path in RULES.glob("*.smk"):
        blocks = re.split(r"^\s*(?:rule|checkpoint) (\w+):", path.read_text(), flags=re.M)
        for name, body in zip(blocks[1::2], blocks[2::2], strict=True):
            t = re.search(r"threads:\s*(?:threads\()?(\d+)", body)
            m = re.search(r"mem_mb=mem\((\d+)\)", body)
            out[name] = (int(t.group(1)) if t else 1, int(m.group(1)) if m else 0)
    return out


def test_every_planned_job_class_matches_its_rules() -> None:
    requests = rule_requests()
    for job, rules in CLASS_RULES.items():
        assert all(r in requests for r in rules), f"{job}: rule missing from workflow/rules"
        threads = max(requests[r][0] for r in rules)
        mem = max(requests[r][1] for r in rules)
        assert JOB_CLASSES[job]["threads"] == threads, f"{job}: plan {JOB_CLASSES[job]['threads']} threads, rules {threads}"
        assert JOB_CLASSES[job]["mem_mb"] == mem, f"{job}: plan {JOB_CLASSES[job]['mem_mb']} MB, rules {mem} MB"


def test_every_resource_requesting_rule_is_planned() -> None:
    planned = {r for rules in CLASS_RULES.values() for r in rules}
    heavy = {name for name, (t, m) in rule_requests().items() if t > 1 or m > 0}
    assert heavy <= planned, f"rules with resource requests missing from the plan: {sorted(heavy - planned)}"
