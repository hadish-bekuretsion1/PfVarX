"""Installation follows the final output: tools a final output does not use are not installed."""

from __future__ import annotations

import itertools
import runpy
import tomllib
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
BOOTSTRAP = runpy.run_path(str(ROOT / "bootstrap.py"))
export, environment_for = BOOTSTRAP["export_workspace"], BOOTSTRAP["environment_for"]
PIXI = tomllib.loads((ROOT / "pixi.toml").read_text())


def features(environment: str) -> set[str]:
    return set(PIXI["environments"][environment].get("features", []))


def test_fo001_installs_alignment_and_qc_only(tmp_path: Path) -> None:
    result = export(ROOT, tmp_path / "fo001", None, "wgs", [], False, "FASTQ", "FO001")
    ws = tmp_path / "fo001"
    assert result["environment"] == "alignment-qc" and result["max_final_output"] == "FO001"
    assert features("alignment-qc") == {"pipeline", "reads"}  # no "calling": no caller, bcftools, numpy, scipy or numba
    base = set(PIXI["dependencies"])
    assert not base & {"bcftools", "numpy", "scipy", "numba"}
    assert set(PIXI["feature"]["calling"]["dependencies"]) == {"bcftools", "numpy", "scipy", "numba"}
    assert not (ws / "src/pfvarx/caller").exists() and not (ws / "workflow/rules/l2_caller.smk").exists()
    assert (ws / "src/pfvarx/caller_adapter.py").is_file()  # imported by the runner; checks for the caller
    for absent in ("src/pfvarx/cohort", "src/pfvarx/annotation", "src/pfvarx/downstream",
                   "workflow/rules/l3_cohort.smk", "workflow/rules/l5_downstream.smk"):
        assert not (ws / absent).exists(), absent
    assert (ws / "workflow/rules/l1_align.smk").is_file()
    targeted = export(ROOT, tmp_path / "fo001t", None, "targeted_umi", [], False, "FASTQ", "FO001")
    assert targeted["environment"] == "alignment-qc-targeted" and features("alignment-qc-targeted") == {"pipeline", "reads", "targeted"}


@pytest.mark.parametrize("final_output, profile, environment, has_cohort", [
    ("FO002", "sample", "core", False),
    ("FO003", "cohort", "core", True),
])
def test_fo002_and_fo003_add_the_caller(tmp_path: Path, final_output: str, profile: str, environment: str,
                                         has_cohort: bool) -> None:
    result = export(ROOT, tmp_path / "ws", None, "wgs", [], False, "FASTQ", final_output)
    assert (result["profile"], result["environment"], result["max_final_output"]) == (profile, environment, final_output)
    assert "calling" in features(environment)
    assert (tmp_path / "ws/src/pfvarx/caller/call.py").is_file()
    assert (tmp_path / "ws/src/pfvarx/cohort").exists() is has_cohort
    assert not (tmp_path / "ws/src/pfvarx/annotation").exists()


def test_fo004_installs_only_the_requested_analyses(tmp_path: Path) -> None:
    result = export(ROOT, tmp_path / "ws", None, "wgs", ["DA002"], False, "FASTQ", "FO004")
    assert result["environment"] == "selected-annotation-reads" and result["max_final_output"] == "FO004"
    assert features(result["environment"]) == {"pipeline", "calling", "annotation", "reads"}  # no R, IQ-TREE, scikit-allel
    with pytest.raises(ValueError, match="--analyses or --auto-analyses"):
        export(ROOT, tmp_path / "none", None, "wgs", [], False, "FASTQ", "FO004")


def test_final_output_and_profile_must_agree(tmp_path: Path) -> None:
    profile_for = BOOTSTRAP["profile_for"]
    assert profile_for("FO003", None) == "cohort" and profile_for(None, None) == "sample"
    with pytest.raises(ValueError, match="installs the sample profile"):
        profile_for("FO001", "cohort")
    with pytest.raises(ValueError, match="FASTQ-only"):
        export(ROOT, tmp_path / "bam", None, "wgs", [], False, "BAM", "FO001")


def test_a_configuration_above_the_installed_final_output_is_refused(tmp_path: Path) -> None:
    from pfvarx.configuration.loader import ConfigError, load_bundle
    from pfvarx.initialization.layout import config_init

    ws = tmp_path / "fo001"
    export(ROOT, ws, None, "wgs", [], False, "FASTQ", "FO001")
    config_init(ws, "analysis")
    path = ws / "config/config.yaml"
    config = yaml.safe_load(path.read_text())
    config["final_output"] = "FO002"
    path.write_text(yaml.safe_dump(config))
    with pytest.raises(ConfigError):
        load_bundle(ws)


def test_every_environment_bootstrap_can_choose_is_defined_and_locked() -> None:
    locked = set(yaml.safe_load((ROOT / "pixi.lock").read_text())["environments"])
    chosen = set()
    analyses = ["DA001", "DA002", "DA003", "DA004", "DA005", "DA006", "DA007"]
    subsets = [list(c) for n in range(len(analyses) + 1) for c in itertools.combinations(analyses, n)]
    for design, input_type, rtg in itertools.product(("wgs", "targeted_umi", "targeted_no_umi"),
                                                     ("FASTQ", "BAM"), (False, True)):
        for selection in subsets:
            chosen.add(environment_for(design, selection, rtg, input_type))
        if input_type == "FASTQ":
            chosen.add(environment_for(design, [], rtg, input_type, "FO001"))
    missing = chosen - set(PIXI["environments"])
    assert not missing, f"bootstrap chooses environments pixi.toml does not define: {sorted(missing)}"
    assert chosen <= locked, f"not in pixi.lock (run pixi lock): {sorted(chosen - locked)}"
    for environment in chosen - {"alignment-qc", "alignment-qc-targeted"}:
        assert "calling" in features(environment), environment  # FO002 and later include the caller


def test_every_installation_has_what_it_runs_and_no_r_without_an_r_analysis(capsys: pytest.CaptureFixture[str]) -> None:
    check = runpy.run_path(str(ROOT / "scripts/check_installs.py"))
    assert check["main"]() == 0, capsys.readouterr().out
    out = capsys.readouterr().out
    assert "MISSING" not in out and "full build with R" not in out
    for line in out.splitlines():  # R only where an analysis runs R scripts
        if "→" in line and not any(da in line for da in ("DA001", "DA004", "DA005", "DA006", "DA007")):
            assert "EXTRA r-base" not in line, line
