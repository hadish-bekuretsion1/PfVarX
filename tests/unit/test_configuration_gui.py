from pathlib import Path

import pytest

from pfvarx.configuration.loader import ConfigError, load_bundle, migrate_project
from pfvarx.configuration.web import editor_state, save
from pfvarx.execution.preflight import prepare_inputs


def test_save_creates_only_input_directories_and_preserves_yaml(workspace: Path) -> None:
    state = editor_state(workspace, "analysis")
    values = state["values"]
    with pytest.raises(ValueError, match="release label"):
        save(workspace, values, state["revision"], False)
    values["config"]["reference"]["release_label"] = "fixture"
    result = save(workspace, values, state["revision"], False)
    assert set(result["directories"]) == {"inputs/fastq", "inputs/reference", "inputs/metadata"}
    assert not list((workspace / "inputs").rglob("*.*"))
    assert not (workspace / "prepared").exists()
    assert save(workspace, values, result["revision"], False)["saved"] == []
    values["config"]["project_id"] = "ANOTHER_PROJECT"
    with pytest.raises(ValueError, match="Save with backup"):
        save(workspace, values, result["revision"], False)
    changed = save(workspace, values, result["revision"], True)
    assert len(changed["backups"]) == 1
    assert "MY_PROJECT" in (workspace / changed["backups"][0]).read_text()
    with pytest.raises(ValueError, match="changed"):
        save(workspace, values, result["revision"], True)
    with pytest.raises(ValueError, match="Created sample sheet template"):
        prepare_inputs(workspace)
    assert not (workspace / "prepared").exists()


def test_auto_selection_uses_declared_targeted_branch(workspace: Path) -> None:
    state = editor_state(workspace, "analysis")
    values = state["values"]
    config = values["config"]
    config.update(final_output="FO004")
    config["reference"].update(release_label="fixture", gff="inputs/reference/test.gff")
    config["assay"].update(design="TARGETED_NON_UMI", sequencing_type="TBS_UMI",
                           amplification_route="AMPLICON", targets_bed="inputs/targets/targets.bed",
                           primers_tsv="inputs/targets/primers.tsv")
    config["analyses"]["selection"] = "auto"
    with pytest.raises(ConfigError, match="sequencing_type is incompatible"):
        save(workspace, values, state["revision"], False)
    config["assay"]["sequencing_type"] = "TBS_NON_UMI"
    saved = save(workspace, values, state["revision"], False)
    assert "inputs/targets" in saved["directories"]
    resolved = load_bundle(workspace).config.analyses
    assert resolved.requested == ["DA001", "DA002", "DA003", "DA004", "DA007"]
    assert resolved.spatial_outcome_source == "DA002"
    assert config["analyses"]["requested"] == []


def test_legacy_config_translates_ids_without_rewriting() -> None:
    legacy = {"schema_version": "pfvarx.config/2", "analyses": {
        "requested": ["DA001", "DA002", "DA003", "DA007"],
        "scales": {"DA001": ["SAMPLE"]}, "spatial_outcome_source": "DA001",
    }}
    new = migrate_project(legacy)
    assert new["analyses"] == {"requested": ["DA002", "DA003", "DA001", "DA007"],
                               "scales": {"DA002": ["SAMPLE"]}, "spatial_outcome_source": "DA002"}
    assert legacy["schema_version"] == "pfvarx.config/2"
    assert legacy["analyses"]["spatial_outcome_source"] == "DA001"


def test_prepare_refuses_when_all_registered_samples_fail(workspace: Path) -> None:
    state = editor_state(workspace, "analysis")
    values = state["values"]
    values["config"]["reference"]["release_label"] = "fixture"
    save(workspace, values, state["revision"], False)
    (workspace / "inputs/reference/Pf3D7.fasta").write_text(">chr1\nACGTACGT\n")
    (workspace / "inputs/reference/core_genome.bed").write_text("chr1\t0\t8\n")
    reads = workspace / "inputs/fastq"
    for mate in ("R1", "R2"):
        (reads / f"invalid_{mate}.fastq").write_text("invalid FASTQ\n")
    (workspace / "inputs/metadata/samples.tsv").write_text(
        "sample_id\tlane\tlibrary_id\tfastq_1\tfastq_2\n"
        "test\tL001\tlib\tinputs/fastq/invalid_R1.fastq\tinputs/fastq/invalid_R2.fastq\n"
    )
    with pytest.raises(ValueError, match="validation failed"):
        prepare_inputs(workspace)
    assert not (workspace / "prepared").exists()


def test_bam_save_creates_alignment_folders(workspace: Path) -> None:
    (workspace / "config/installation.json").write_text(
        '{"schema":"pfvarx.installation/2","profile":"sample","purpose":"analysis",'
        '"design":"WGS","input_type":"FASTQ","max_final_output":"FO002",'
        '"analyses_requested":[],"analyses_available":[]}'
    )
    state = editor_state(workspace, "analysis")
    values = state["values"]
    values["config"]["input_type"] = "BAM"
    values["config"]["reference"]["release_label"] = "fixture"
    result = save(workspace, values, state["revision"], False)
    assert set(result["directories"]) == {"inputs/alignments", "inputs/reference", "inputs/metadata"}
    assert all((workspace / folder).is_dir() for folder in result["directories"])
    assert not (workspace / "inputs/fastq").exists()
