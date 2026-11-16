from pathlib import Path

import pytest
import yaml

from pfvarx.configuration.loader import ConfigError, effective_config, load_bundle, load_parameters, param
from pfvarx.initialization.layout import check_layout, config_init, init_inputs


def _init(ws: Path, purpose: str = "analysis") -> None:
    assert len(config_init(ws, purpose)) == 3


def test_config_init_never_overwrites(workspace: Path) -> None:
    _init(workspace)
    (workspace / "config" / "config.yaml").write_text("edited: true\n")
    assert config_init(workspace, "analysis") == []
    assert (workspace / "config" / "config.yaml").read_text() == "edited: true\n"


def test_templates_validate_for_both_purposes(workspace: Path, tmp_path: Path) -> None:
    for purpose in ("analysis", "benchmarking"):
        for f in ("config", "quality", "compute"):
            (workspace / "config" / f"{f}.yaml").unlink(missing_ok=True)
        _init(workspace, purpose)
        b = load_bundle(workspace)
        assert b.config.purpose.value == purpose


def test_unknown_key_rejected_with_location(workspace: Path) -> None:
    _init(workspace)
    p = workspace / "config" / "config.yaml"
    d = yaml.safe_load(p.read_text())
    d["disable_qc"] = True
    p.write_text(yaml.safe_dump(d))
    with pytest.raises(ConfigError) as exc:
        load_bundle(workspace)
    assert any("disable_qc" in m for m in exc.value.problems)


def test_targeted_design_requires_the_library_catalog(workspace: Path) -> None:
    _init(workspace)
    p = workspace / "config" / "config.yaml"
    d = yaml.safe_load(p.read_text())
    d["assay_design"] = "targeted_umi"
    p.write_text(yaml.safe_dump(d))
    with pytest.raises(ConfigError) as exc:
        load_bundle(workspace)
    assert any("library_catalog" in m for m in exc.value.problems)
    d["assay_design"] = "TARGETED_NON_UMI"  # schema-3 spelling is read as targeted_no_umi
    d["assay"]["library_catalog"] = "inputs/metadata/library_catalog.yaml"
    p.write_text(yaml.safe_dump(d))
    assert load_bundle(workspace).config.assay_design.value == "targeted_no_umi"
    d["assay_design"] = "amplicon"
    p.write_text(yaml.safe_dump(d))
    with pytest.raises(ConfigError, match="wgs, targeted_umi or targeted_no_umi"):
        load_bundle(workspace)


def test_fo004_requires_analyses_and_da007_source(workspace: Path) -> None:
    _init(workspace)
    p = workspace / "config" / "config.yaml"
    d = yaml.safe_load(p.read_text())
    d["final_output"] = "FO004"
    p.write_text(yaml.safe_dump(d))
    with pytest.raises(ConfigError):
        load_bundle(workspace)
    d["analyses"]["requested"] = ["DA007"]
    p.write_text(yaml.safe_dump(d))
    with pytest.raises(ConfigError) as exc:
        load_bundle(workspace)
    assert any("spatial_outcome_source" in m for m in exc.value.problems)


def test_compute_envelope_checked(workspace: Path) -> None:
    _init(workspace)
    p = workspace / "config" / "compute.yaml"
    d = yaml.safe_load(p.read_text())
    d["per_job"]["max_mem_mb"] = d["max_mem_mb"] + 1
    p.write_text(yaml.safe_dump(d))
    with pytest.raises(ConfigError):
        load_bundle(workspace)


def test_every_parameter_records_source_and_unset_blocks(workspace: Path) -> None:
    params = load_parameters(workspace)
    assert params and all("source" in v for v in params.values())
    params["layer2.k_max"] = {"value": None, "source": "UNSET_BLOCKING"}
    with pytest.raises(ConfigError):
        param(params, "layer2.k_max")


def test_effective_config_hash_is_deterministic(workspace: Path) -> None:
    _init(workspace)
    assert effective_config(workspace)["hash"] == effective_config(workspace)["hash"]


def test_init_creates_no_result_dirs_and_layout_is_read_only(workspace: Path) -> None:
    _init(workspace)
    cfg = load_bundle(workspace).config
    out = init_inputs(workspace, cfg)
    assert out["result_dirs_created"] == []
    assert not (workspace / "results").exists() and not (workspace / "runs").exists()
    before = sorted(p.as_posix() for p in workspace.rglob("*"))
    rep = check_layout(workspace, cfg)
    assert rep["status"] == "FAIL" and "reference.fasta" in rep["missing_files"]
    assert sorted(p.as_posix() for p in workspace.rglob("*")) == before


def test_json_out_mirrors_stdout(tmp_path: Path, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    """PFVARX_JSON_OUT receives exactly the JSON printed on stdout, written atomically."""
    import json as _json

    import typer
    from typer.testing import CliRunner

    from pfvarx.cli import emit

    out = tmp_path / "sub" / "result.json"
    monkeypatch.setenv("PFVARX_JSON_OUT", str(out))
    app = typer.Typer()

    @app.command()
    def show() -> None:
        emit({"status": "OK", "n": 3})

    res = CliRunner().invoke(app, [])
    assert res.exit_code == 0
    assert _json.loads(res.stdout) == {"status": "OK", "n": 3}
    assert _json.loads(out.read_text()) == {"status": "OK", "n": 3}
    assert not (tmp_path / "sub" / "result.json.tmp").exists()


def test_fo004_requires_annotation_path(workspace: Path) -> None:
    _init(workspace)
    path = workspace / "config/config.yaml"
    data = yaml.safe_load(path.read_text())
    data["final_output"] = "FO004"
    data["analyses"]["requested"] = ["DA001"]
    data["reference"]["gff"] = None
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ConfigError, match="reference.gff"):
        load_bundle(workspace)
    data["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    path.write_text(yaml.safe_dump(data))
    assert load_bundle(workspace).config.reference.gff == Path("inputs/reference/Pf3D7.gff")


def test_standard_input_folders_resolve_to_template_files(workspace: Path) -> None:
    _init(workspace)
    path = workspace / "config/config.yaml"
    config = yaml.safe_load(path.read_text())
    config["reference"].update(fasta="inputs/reference/", core_mask_bed="inputs/reference/")
    config["samples_tsv"] = "inputs/metadata/"
    path.write_text(yaml.safe_dump(config))
    resolved = load_bundle(workspace).config
    assert resolved.reference.fasta == Path("inputs/reference/Pf3D7.fasta")
    assert resolved.reference.core_mask_bed == Path("inputs/reference/core_genome.bed")
    assert resolved.samples_tsv == Path("inputs/metadata/samples.tsv")


def test_snpeff_genome_id_is_safe_for_any_release_label() -> None:
    """A release label with spaces ('Pf3D7 - release unknown') broke the SnpEff database build."""
    import pytest
    from pydantic import ValidationError

    from pfvarx.annotation.snpeff_db import GENOME_ID_PATTERN, genome_id
    from pfvarx.configuration.models import ReferenceConfig

    assert genome_id(None, "Pf3D7 - release unknown") == "pfvarx_Pf3D7_release_unknown"
    assert genome_id(None, "PlasmoDB 68") == "pfvarx_PlasmoDB_68"
    assert genome_id(None, "synthetic") == "pfvarx_synthetic"
    assert genome_id(None, " ¿? ") == "pfvarx_reference"
    assert genome_id("Pf3D7_study", "anything at all") == "Pf3D7_study"
    import re
    assert all(re.match(GENOME_ID_PATTERN, genome_id(None, label)) for label in ("a b", "x/y", "release: 2020-09"))
    base = {"fasta": "f.fasta", "core_mask_bed": "c.bed", "release_label": "r"}
    assert ReferenceConfig(**base, snpeff_genome_id="Pf3D7_v3.1").snpeff_genome_id == "Pf3D7_v3.1"
    with pytest.raises(ValidationError):
        ReferenceConfig(**base, snpeff_genome_id="Pf3D7 v3")
