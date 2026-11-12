"""Integration tests on the synthetic fixture (need the pixi environment's tools)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.integration
REPO = Path(__file__).resolve().parents[2]
NEED = ("samtools", "bwa-mem2", "fastp", "snakemake", "picard")


def _ws(tmp: Path, design: str, fo: str) -> Path:
    if any(shutil.which(t) is None for t in NEED):
        pytest.skip("bioinformatics tools not on PATH (run inside 'pixi run')")
    ws = tmp / "ws"
    # Exercise the current source, including a source archive with no Git commits.
    ws.mkdir()
    for name in ("config", "workflow", "src", "pixi.toml", "pixi.lock", "pyproject.toml"):
        source = REPO / name
        if source.is_dir():
            shutil.copytree(source, ws / name)
        else:
            shutil.copyfile(source, ws / name)
    for name in ("config.yaml", "quality.yaml", "compute.yaml", "profiles/resources.yaml"):
        (ws / "config" / name).unlink(missing_ok=True)
    run = lambda *a: subprocess.run(["pfvarx", *a], cwd=ws, capture_output=True, text=True)  # noqa: E731
    assert run("config-init", "--purpose", "analysis").returncode == 0
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["final_output"] = fo
    cfg["assay"]["design"] = design
    if design != "WGS":
        cfg["assay"].update(amplification_route="AMPLICON", targets_bed="inputs/targets/targets.bed",
                            primers_tsv="inputs/targets/primers.tsv")
    if design == "TARGETED_UMI":
        cfg["assay"]["read_structure"] = "8M+T 8M+T"
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    assert run("init").returncode == 0
    assert run("fixture", "--design", design, "--depth", "15").returncode == 0
    assert run("validate").returncode == 0
    assert run("prepare").returncode == 0
    return ws


@pytest.mark.parametrize("design", ["WGS", "TARGETED_NON_UMI", "TARGETED_UMI"])
def test_dry_run_creates_no_outputs(tmp_path: Path, design: str) -> None:
    ws = _ws(tmp_path, design, "FO001")
    r = subprocess.run(["pfvarx", "dry-run"], cwd=ws, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert not (ws / "runs").exists() and not (ws / "results").exists()
