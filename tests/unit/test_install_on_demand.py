"""On-demand installation: a control environment first, then the tools the saved configuration needs."""

from __future__ import annotations

import runpy
import subprocess
from pathlib import Path
from typing import Any

import pytest
import yaml

from pfvarx import installation
from pfvarx.configuration.loader import load_bundle
from pfvarx.initialization.workspace import initialize
from pfvarx.io import read_json

ROOT = Path(__file__).resolve().parents[2]
BOOTSTRAP = runpy.run_path(str(ROOT / "bootstrap.py"))


def on_demand_workspace(tmp_path: Path) -> Path:
    ws = tmp_path / "ws"
    meta = BOOTSTRAP["export_workspace"](ROOT, ws, "on-demand", "wgs", [], False, "FASTQ", None)
    assert meta["environment"] == "control" and meta["profile"] == "on-demand"
    initialize(ws, config_only=True)
    config = ws / "config/config.yaml"
    data = yaml.safe_load(config.read_text())
    data["reference"]["release_label"] = "synthetic"
    config.write_text(yaml.safe_dump(data))
    return ws


def set_config(ws: Path, **values: Any) -> None:
    path = ws / "config/config.yaml"
    data = yaml.safe_load(path.read_text())
    data.update(values)
    path.write_text(yaml.safe_dump(data))


def fake_pixi(monkeypatch: pytest.MonkeyPatch, ws: Path, calls: list[list[str]]) -> None:
    def run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        python = installation.environment_python(ws, command[-1])
        python.parent.mkdir(parents=True, exist_ok=True)
        python.write_text("")
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(installation.shutil, "which", lambda name: "/usr/bin/pixi")
    monkeypatch.setattr(installation.subprocess, "run", run)


def test_on_demand_workspace_starts_with_the_control_environment_and_all_code(tmp_path: Path) -> None:
    ws = on_demand_workspace(tmp_path)
    for code in ("caller", "cohort", "annotation", "downstream"):
        assert (ws / "src/pfvarx" / code).is_dir(), code  # every step is possible; tools come later
    assert load_bundle(ws).config.final_output.value == "FO002"  # the template decides, not the installation
    set_config(ws, final_output="FO003")
    assert load_bundle(ws).config.final_output.value == "FO003"  # no installed limit on demand


def test_steps_that_run_tools_stop_until_the_tools_are_installed(tmp_path: Path) -> None:
    ws = on_demand_workspace(tmp_path)
    state = installation.status(ws)
    assert state["current"] is False and state["required_environment"] == "core"
    assert "Install required tools" in state["detail"]
    from pfvarx.execution.preflight import prepare_inputs
    from pfvarx.execution.runner import build_plan
    for step in (lambda: prepare_inputs(ws), lambda: build_plan(ws)):
        with pytest.raises(ValueError, match="TOOLS_NOT_INSTALLED"):
            step()


def test_install_tools_installs_only_what_the_configuration_needs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ws = on_demand_workspace(tmp_path)
    calls: list[list[str]] = []
    fake_pixi(monkeypatch, ws, calls)
    set_config(ws, final_output="FO001")
    result = installation.install(ws)
    assert result["status"] == "INSTALLED" and result["environment"] == "alignment-qc"
    assert calls == [["/usr/bin/pixi", "install", "--locked", "-e", "alignment-qc"]]
    assert read_json(ws / "config/installation.json")["environment"] == "alignment-qc"
    assert installation.status(ws)["current"] is True
    assert installation.stage_python(ws, "fallback").endswith(".pixi/envs/alignment-qc/bin/python")
    assert installation.install(ws)["status"] == "ALREADY_INSTALLED" and len(calls) == 1
    set_config(ws, final_output="FO002")  # a later final output needs the caller
    assert installation.status(ws)["current"] is False
    assert installation.install(ws)["environment"] == "core" and len(calls) == 2


def test_fo004_installs_the_requested_analyses_with_their_prerequisites(tmp_path: Path) -> None:
    ws = on_demand_workspace(tmp_path)
    set_config(ws, final_output="FO004")
    data = yaml.safe_load((ws / "config/config.yaml").read_text())
    data["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    data["analyses"]["requested"] = ["DA006"]
    (ws / "config/config.yaml").write_text(yaml.safe_dump(data))
    need = installation.required(ws)
    assert "DA001" in need["analyses"]  # phylogeny needs the clonality calls
    assert need["environment"] == "selected-annotation-coi-phylogeny-reads"


def test_up_front_installations_are_unaffected(tmp_path: Path) -> None:
    ws = tmp_path / "fo002"
    BOOTSTRAP["export_workspace"](ROOT, ws, None, "wgs", [], False, "FASTQ", "FO002")
    assert installation.status(ws)["mode"] == "up-front"
    installation.require_tools(ws)  # nothing to install
    assert installation.install(ws)["status"] == "NOT_NEEDED"


def test_gui_steps_find_the_installed_tools_first_on_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ws = on_demand_workspace(tmp_path)
    fake_pixi(monkeypatch, ws, [])
    installation.install(ws)
    env = installation.stage_environment(ws, {"PATH": "/control/bin:/usr/bin"})
    tools_bin = str(ws / ".pixi/envs/core/bin")
    assert env["PATH"].split(":")[0] == tools_bin  # samtools, bwa-mem2 ... found before the control env
    assert env["CONDA_PREFIX"] == str(ws / ".pixi/envs/core")
