"""The CLI and the GUI configure and run the same workspace, so they must change together.

Every setting in the configuration model must be editable in the YAML templates (CLI),
explained in the GUI field guide, and listed in the configuration guide with exactly the GUI
text; every GUI tab must name a terminal command that exists, and every command the GUI runs
must be a CLI command. A failure here means one side changed without the other.
"""

from __future__ import annotations

import json
import re
import runpy
from pathlib import Path
from typing import Any

import yaml

from pfvarx.cli import app
from pfvarx.configuration import web
from pfvarx.configuration.models import Bundle

REPO = Path(__file__).resolve().parents[2]
GUI = REPO / "src/pfvarx/configuration/gui"
GUIDE_TOOL = runpy.run_path(str(REPO / "scripts/configuration_guide.py"))


def _model_settings() -> set[str]:
    schema = Bundle.model_json_schema()
    defs = schema.get("$defs", {})

    def resolve(node: dict[str, Any]) -> dict[str, Any]:
        if "$ref" in node:
            return resolve(defs[node["$ref"].split("/")[-1]])
        for key in ("allOf", "anyOf"):
            refs = [x for x in node.get(key, []) if "$ref" in x]
            if refs and "properties" not in node:
                return resolve(refs[0])
        return node

    def leaves(node: dict[str, Any], prefix: str) -> set[str]:
        node = resolve(node)
        if "properties" not in node:
            return {prefix}
        return set().union(*(leaves(v, f"{prefix}.{k}") for k, v in node["properties"].items()))

    return set().union(*(leaves(schema["properties"][top], top) for top in ("config", "quality", "compute")))


def _template_settings(name: str, prefix: str) -> set[str]:
    def flat(data: dict[str, Any], path: str) -> set[str]:
        out: set[str] = set()
        for key, value in data.items():
            here = f"{path}.{key}"
            # scales is a free-form mapping (analysis -> scales), not nested settings
            out |= flat(value, here) if isinstance(value, dict) and value and key != "scales" else {here}
        return out

    return flat(yaml.safe_load((REPO / f"config/templates/{name}.yaml").read_text()), prefix)


def _cli_commands() -> set[str]:
    return {c.name or c.callback.__name__.replace("_", "-") for c in app.registered_commands}


SETTINGS = _model_settings()


def test_every_setting_is_explained_in_the_gui() -> None:
    guide = set(json.loads((GUI / "field-guide.json").read_text()))
    assert SETTINGS - guide == set(), "settings without GUI help"
    assert guide - SETTINGS == set(), "GUI help for settings that no longer exist"


def test_every_setting_is_editable_in_the_cli_templates() -> None:
    quality_compute = _template_settings("quality", "quality") | _template_settings("compute", "compute")
    for purpose in ("analysis", "benchmarking"):
        templates = _template_settings(purpose, "config") | quality_compute
        assert SETTINGS - templates == set(), f"settings missing from the {purpose} templates"
        assert templates - SETTINGS == set(), f"{purpose} template keys the model does not accept"


def test_configuration_guide_lists_every_setting_with_the_gui_text() -> None:
    documented = GUIDE_TOOL["documented_keys"](REPO)
    for kind in ("md", "html"):
        assert documented[kind] == SETTINGS, f"configuration guide .{kind} settings differ from the model"
    stale = [path.name for path, text in GUIDE_TOOL["render"](REPO).items() if path.read_text() != text]
    assert not stale, f"{stale} differ from field-guide.json; run python scripts/configuration_guide.py"


def test_every_gui_tab_names_an_existing_terminal_command() -> None:
    page = (GUI / "index.html").read_text()
    tabs = re.findall(r'<button type="button" id="([a-z-]+)-tab" aria-controls="([a-z-]+)-panel"', page)
    commands = _cli_commands()
    assert tabs
    for _, panel in tabs:
        section = page.split(f'<section id="{panel}-panel"', 1)[1].split("</section>", 1)[0]
        named = re.findall(r"\./pfvarx ([a-z-]+)", section)
        assert "Terminal equivalent" in section and named, f"{panel} tab names no terminal command"
        assert set(named) <= commands, f"{panel} tab names unknown commands {set(named) - commands}"


def test_every_command_the_gui_runs_is_a_cli_command() -> None:
    source = (REPO / "src/pfvarx/configuration/web.py").read_text()
    checks = source.split("checks = {", 1)[1].split("}", 1)[0]  # GUI button -> CLI stages it runs
    stage_lists = re.findall(r'"[a-z-]+": \(([^)]*)\)', checks)
    stages = set(web.PREFLIGHT_STAGES) | {s for group in stage_lists for s in re.findall(r'"([a-z-]+)"', group)}
    assert stage_lists and stages <= _cli_commands(), f"GUI runs commands the CLI lacks: {stages - _cli_commands()}"


def test_orchestration_guide_documents_every_cli_command() -> None:
    """The command reference must list every CLI command, including the GUI-matching ones."""
    guide = (REPO / "docs/orchestration.md").read_text()
    missing = sorted(c for c in _cli_commands() if f"`{c}" not in guide)
    assert not missing, f"docs/orchestration.md does not document {missing}"
