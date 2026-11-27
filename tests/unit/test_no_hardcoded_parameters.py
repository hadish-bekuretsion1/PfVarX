"""Maintainer parameters are written once, in config/parameters/layer*.yaml.

Code reads them from the run's frozen configuration. It must not repeat a value as an inline
fallback (``params.get("key", 200)``): a run frozen before a parameter existed takes the value
recorded in ``configuration/legacy_parameters.yaml`` through ``frozen_param``.
"""

import re
from pathlib import Path

import pytest
import yaml

from pfvarx.configuration.loader import LEGACY_PARAMETERS, frozen_param, load_parameters

ROOT = Path(__file__).resolve().parents[2]
PARAMS = load_parameters(ROOT)
BARE = {key.split(".", 1)[1] for key in PARAMS}
CODE = [*(ROOT / "src" / "pfvarx").rglob("*.py"), *(ROOT / "workflow").rglob("*.smk"),
        ROOT / "workflow" / "Snakefile", *(ROOT / "workflow" / "scripts").glob("*.py")]


def test_parameter_names_are_unique_across_layers() -> None:
    assert len(BARE) == len(PARAMS), "the same parameter name is used in two layer files"


def test_numbers_are_not_read_as_text() -> None:
    """PyYAML reads 1.0e7 (no exponent sign) as a string; write 1.0e+7."""
    number = re.compile(r"^[-+]?\d*\.?\d+([eE][-+]?\d+)?$")
    for key, entry in PARAMS.items():
        values = entry["value"] if isinstance(entry["value"], list) else [entry["value"]]
        for value in values:
            assert not (isinstance(value, str) and number.match(value)), f"{key}: {value!r} is text, not a number"


def test_no_inline_fallback_for_a_parameter() -> None:
    fallback = re.compile(r"""\.get\(\s*["'](?:layer\d\.)?(\w+)["']\s*,\s*[^)\s]""")
    found = []
    for path in CODE:
        for number, line in enumerate(path.read_text().splitlines(), start=1):
            for match in fallback.finditer(line):
                if match.group(1) in BARE:
                    found.append(f"{path.relative_to(ROOT)}:{number}: {match.group(0)}")
    assert not found, "parameters repeated as inline fallbacks (use frozen_param):\n" + "\n".join(found)


def test_legacy_values_name_current_parameters() -> None:
    legacy = yaml.safe_load(LEGACY_PARAMETERS.read_text())["parameters"]
    assert set(legacy) <= BARE, f"legacy values for unknown parameters: {sorted(set(legacy) - BARE)}"


def test_frozen_param_prefers_the_frozen_value_and_never_invents_one() -> None:
    assert frozen_param({"pairhmm_max_reads": 7}, "pairhmm_max_reads") == 7
    assert frozen_param({}, "pairhmm_max_reads") == 200  # the value runs frozen before it existed used
    with pytest.raises(KeyError):
        frozen_param({}, "no_such_parameter")
