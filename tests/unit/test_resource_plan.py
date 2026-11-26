"""Thread and memory requests are written once, in ``config/registry/job_resources.yaml``.

The workflow rules read their requests from that file by rule name (``threads("rule")`` and
``mem("rule")`` in ``workflow/rules/common.smk``) and the preflight resource plan reads the
same file, so the plan always reports what the rules request. These tests fail when a rule
writes a number instead, or when a rule with requests is missing from the file.
"""

import re
from pathlib import Path

from pfvarx.planning import resources

ROOT = Path(__file__).resolve().parents[2]
RULES = ROOT / "workflow" / "rules"


def rule_bodies() -> dict[str, str]:
    out: dict[str, str] = {}
    for path in RULES.glob("*.smk"):
        blocks = re.split(r"^\s*(?:rule|checkpoint) (\w+):", path.read_text(), flags=re.M)
        out.update(zip(blocks[1::2], blocks[2::2], strict=True))
    return out


def test_rules_request_resources_by_name_only() -> None:
    for name, body in rule_bodies().items():
        threads = re.search(r"threads:\s*(.+)", body)
        if threads:
            assert threads.group(1).strip() == f'threads("{name}")', f"{name}: threads not read from job_resources.yaml"
        memory = re.search(r"mem_mb=(.+?),?\n", body)
        if memory:
            assert re.fullmatch(rf'mem\("{name}"(, host=True)?(, job="\w+")?\)', memory.group(1).strip()), \
                f"{name}: memory not read from job_resources.yaml"


def test_every_resource_requesting_rule_is_in_the_file() -> None:
    listed = resources.rule_requests(ROOT)
    requesting = {name for name, body in rule_bodies().items() if "threads(" in body or "mem(" in body}
    assert requesting <= set(listed), f"rules missing from job_resources.yaml: {sorted(requesting - set(listed))}"
    assert set(listed) <= set(rule_bodies()), f"job_resources.yaml lists unknown rules: {sorted(set(listed) - set(rule_bodies()))}"


def test_job_class_requests_are_the_largest_rule_requests() -> None:
    rules = resources.rule_requests(ROOT)
    for name, jc in resources.job_classes(ROOT).items():
        if jc["rules"]:
            assert jc["threads"] == max(rules[r]["threads"] for r in jc["rules"]), name
            assert jc["mem_mb"] == max(rules[r]["mem_mb"] for r in jc["rules"]), name


def test_planning_factors_are_present() -> None:
    factors = resources.planning(ROOT)
    for key in ("default_rule_mem_mb", "measured_memory_headroom", "tmp_disk_per_input_factor",
                "slow_step_factor", "progress_refresh_seconds"):
        assert factors[key] > 0, key
