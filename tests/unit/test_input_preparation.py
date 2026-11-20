"""Inputs Preparation: sample metadata from the input files; Study Metadata (Yes/No) is optional."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from pfvarx.configuration.loader import ConfigError, load_bundle
from pfvarx.configuration.web import editor_state, save
from pfvarx.initialization.layout import required_dirs
from pfvarx.inputs.preparation import prepare_samples
from pfvarx.inputs.transfer import destination_for


def _workspace(workspace: Path, **config: object) -> Path:
    state = editor_state(workspace, "analysis")
    state["values"]["config"]["reference"]["release_label"] = "fixture"
    state["values"]["config"].update(config)
    save(workspace, state["values"], state["revision"], False)
    for sample in ("S1", "S2"):
        for mate in ("R1", "R2"):
            (workspace / f"inputs/fastq/{sample}_L001_{mate}.fastq.gz").touch()
    return workspace


def _rows(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    return [dict(zip(lines[0].split("\t"), line.split("\t"), strict=False)) for line in lines[1:]]


def test_input_directories_create_only_folders_for_supplied_files(workspace: Path) -> None:
    _workspace(workspace)
    cfg = load_bundle(workspace).config
    assert cfg.study_metadata is False  # Study Metadata: No by default
    assert sorted(p.name for p in (workspace / "inputs").iterdir()) == ["fastq", "reference"]
    for later in ("prepared", "metadata", "runs", "results"):  # created on demand by later steps
        assert not (workspace / later).exists()
    assert "inputs/metadata" in required_dirs(cfg.model_copy(update={"study_metadata": True}))
    targeted = cfg.model_copy(update={"assay_design": cfg.assay_design.TARGETED_NO_UMI})
    assert "inputs/metadata" in required_dirs(targeted)  # the library catalog is supplied there
    with pytest.raises(ValueError, match="unavailable"):  # no metadata uploads with Study Metadata: No
        destination_for(workspace, "study.tsv", "metadata")


def test_study_metadata_yes_creates_its_folder_for_uploads(workspace: Path) -> None:
    _workspace(workspace, study_metadata=True)
    assert (workspace / "inputs/metadata").is_dir()
    assert destination_for(workspace, "study.tsv", "metadata") == workspace / "inputs/metadata/study.tsv"


def test_study_metadata_no_continues_with_the_sample_metadata(workspace: Path) -> None:
    _workspace(workspace)
    result = prepare_samples(workspace)
    assert result["sample_sheet"]["status"] == "CREATED"
    assert result["sample_sheet"]["sample_sheet"] == "prepared/samples.tsv"  # folder created on demand
    assert result["study_metadata"]["merge"] == "NOT_USED"
    rows = _rows(workspace / "prepared/samples.tsv")
    assert [(r["sample_id"], r["lane"]) for r in rows] == [("S1", "L001"), ("S2", "L001")]
    assert json.loads((workspace / "metadata/inputs_preparation.json").read_text()) == result
    assert prepare_samples(workspace)["sample_sheet"]["status"] == "EXISTING"  # never overwritten


def test_study_metadata_yes_merges_without_overwriting(workspace: Path) -> None:
    _workspace(workspace, study_metadata=True)
    (workspace / "inputs/metadata/study.tsv").write_text(
        "sample_id\tcountry\tcollection_date\nS1\tEthiopia\t2025-01-15\nS3\tKenya\t2025-02-01\n")
    result = prepare_samples(workspace)["study_metadata"]
    assert result["merge"] == "MERGED" and result["cells_filled"] == 2
    assert result["samples_without_study_row"] == ["S2"] and result["study_rows_without_sample"] == ["S3"]
    rows = {r["sample_id"]: r for r in _rows(workspace / "prepared/samples.tsv")}
    assert rows["S1"]["country"] == "Ethiopia" and rows["S2"]["country"] == ""
    assert prepare_samples(workspace)["study_metadata"]["status"] == "NOTHING_TO_WRITE"  # safe to repeat


def test_study_metadata_yes_without_a_table_continues_and_says_so(workspace: Path) -> None:
    _workspace(workspace, study_metadata=True)
    result = prepare_samples(workspace)["study_metadata"]
    assert result["merge"] == "NO_STUDY_TABLE" and "Continued with the sample metadata only" in result["warning"]


def test_study_metadata_no_ignores_a_table_in_the_folder(workspace: Path) -> None:
    _workspace(workspace)
    (workspace / "inputs/metadata").mkdir()
    (workspace / "inputs/metadata/study.tsv").write_text("sample_id\tcountry\nS1\tEthiopia\n")
    assert prepare_samples(workspace)["study_metadata"]["merge"] == "NOT_USED"
    assert _rows(workspace / "prepared/samples.tsv")[0].get("country", "") == ""  # nothing merged


def test_several_study_tables_stop_and_ask_which_one(workspace: Path) -> None:
    _workspace(workspace, study_metadata=True)
    (workspace / "inputs/metadata/a.tsv").write_text("sample_id\tcountry\nS1\tEthiopia\n")
    (workspace / "inputs/metadata/b.tsv").write_text("sample_id\tcountry\nS1\tSudan\n")
    with pytest.raises(ValueError, match="Several study metadata tables.*--metadata-table"):
        prepare_samples(workspace)
    assert prepare_samples(workspace, "inputs/metadata/b.tsv")["study_metadata"]["merge"] == "MERGED_CHOSEN_TABLE"
    assert _rows(workspace / "prepared/samples.tsv")[0]["country"] == "Sudan"


def test_catalogues_are_not_taken_for_study_tables(workspace: Path) -> None:
    _workspace(workspace, study_metadata=True)
    (workspace / "inputs/metadata/markers.tsv").write_text("#catalog_id=x\nmarker_id\tgene_id\nm1\tg1\n")
    (workspace / "inputs/metadata/library_catalog.yaml").write_text("schema_version: pfvarx.library_catalog/1\n")
    assert prepare_samples(workspace)["study_metadata"]["merge"] == "NO_STUDY_TABLE"


def test_bam_input_gets_a_blank_sample_metadata_template(workspace: Path) -> None:
    state = editor_state(workspace, "analysis")
    state["values"]["config"]["reference"]["release_label"] = "fixture"
    state["values"]["config"].update(input_type="BAM", final_output="FO002")
    save(workspace, state["values"], state["revision"], False)
    with pytest.raises(ValueError, match="blank sample metadata sheet at prepared/samples.tsv"):
        prepare_samples(workspace)
    assert (workspace / "prepared/samples.tsv").read_text().startswith("sample_id\t")


def test_templates_default_to_generated_sample_metadata_and_no_study_metadata(workspace: Path) -> None:
    template = yaml.safe_load((Path(__file__).resolve().parents[2] / "config/templates/analysis.yaml").read_text())
    assert template["samples_tsv"] == "prepared/samples.tsv"
    assert template["study_metadata"] is False and template["study_metadata_dir"] == "inputs/metadata"
    with pytest.raises(ConfigError):  # the folder is a path, not free text
        workspace.joinpath("config").mkdir(exist_ok=True)
        (workspace / "config/config.yaml").write_text(yaml.safe_dump({**template, "study_metadata": "maybe"}))
        load_bundle(workspace)
