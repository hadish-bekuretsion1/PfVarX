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
    assert set(result["directories"]) == {"inputs/fastq", "inputs/reference"}
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
    with pytest.raises(ValueError, match="No supported paired FASTQ names"):
        prepare_inputs(workspace)
    assert not (workspace / "prepared").exists()


def test_auto_selection_uses_declared_targeted_branch(workspace: Path) -> None:
    state = editor_state(workspace, "analysis")
    values = state["values"]
    config = values["config"]
    config.update(final_output="FO004")
    config["reference"].update(release_label="fixture", gff="inputs/reference/test.gff")
    config["assay_design"] = "targeted_no_umi"
    config["analyses"]["selection"] = "auto"
    with pytest.raises(ConfigError, match="library_catalog"):
        save(workspace, values, state["revision"], False)
    config["assay"]["library_catalog"] = "inputs/metadata/library_catalog.yaml"
    saved = save(workspace, values, state["revision"], False)
    assert "inputs/targets" in saved["directories"]
    resolved = load_bundle(workspace).config.analyses
    # Every analysis applies to every design; the panel and callable evidence decide at run time.
    assert resolved.requested == ["DA001", "DA002", "DA003", "DA004", "DA005", "DA006", "DA007"]
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
    (workspace / "prepared").mkdir()
    (workspace / "prepared/samples.tsv").write_text(
        "sample_id\tlane\tlibrary_id\tfastq_1\tfastq_2\n"
        "test\tL001\tlib\tinputs/fastq/invalid_R1.fastq\tinputs/fastq/invalid_R2.fastq\n"
    )
    with pytest.raises(ValueError, match="validation failed"):
        prepare_inputs(workspace)
    assert not (workspace / "prepared/reference").exists()  # no reference derivatives were built


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
    assert set(result["directories"]) == {"inputs/bam", "inputs/reference"}
    assert all((workspace / folder).is_dir() for folder in result["directories"])
    assert not (workspace / "inputs/fastq").exists()


def test_validation_discovers_fastq_sample_sheet(workspace: Path) -> None:
    from pfvarx.execution.preflight import validate_inputs
    from pfvarx.io import read_tsv

    state = editor_state(workspace, "analysis")
    state["values"]["config"]["reference"]["release_label"] = "fixture"
    save(workspace, state["values"], state["revision"], False)
    for mate in (1, 2):
        (workspace / f"inputs/fastq/sample_{mate}.fastq").write_text("@read\nACGT\n+\nIIII\n")
    validate_inputs(workspace)
    sheet = workspace / "prepared/samples.tsv"
    row, = read_tsv(sheet)
    assert row["sample_id"] == "sample"
    assert row["fastq_1"].endswith("/sample_1.fastq")
    assert row["fastq_2"].endswith("/sample_2.fastq")
    original = sheet.read_bytes()
    validate_inputs(workspace)
    assert sheet.read_bytes() == original


def test_validation_continues_after_missing_sample_read(workspace: Path) -> None:
    from pfvarx.execution.preflight import validate_inputs

    state = editor_state(workspace, "analysis")
    state["values"]["config"]["reference"]["release_label"] = "fixture"
    save(workspace, state["values"], state["revision"], False)
    reference = workspace / "inputs/reference"
    (reference / "Pf3D7.fasta").write_text(">chr1\nACGT\n")
    (reference / "core_genome.bed").write_text("chr1\t0\t4\n")
    for mate in (1, 2):
        (workspace / f"inputs/fastq/good_{mate}.fastq").write_text("@read\nACGT\n+\nIIII\n")
    validate_inputs(workspace)
    sheet = workspace / "prepared/samples.tsv"
    text = sheet.read_text()
    sheet.write_text(text + text.splitlines()[1].replace("good", "missing") + "\n")

    report = validate_inputs(workspace)

    assert report["status"] == "PARTIAL"
    assert report["samples"] == {"good": "PASS", "missing": "FAIL"}
    assert [row["sample_id"] for row in report["registry"]] == ["good"]


def _saved_workspace(workspace: Path, input_type: str = "FASTQ") -> None:
    state = editor_state(workspace, "analysis")
    values = state["values"]
    values["config"]["reference"]["release_label"] = "fixture"
    values["config"]["input_type"] = input_type
    if input_type != "FASTQ":
        values["config"]["final_output"] = "FO002"
    save(workspace, values, state["revision"], False)


def test_gui_and_cli_share_transfer_rules(workspace: Path, tmp_path_factory: pytest.TempPathFactory) -> None:
    """Upload/download in the GUI and add-input/fetch-input on the CLI use one rule set."""
    from pfvarx.configuration import web
    from pfvarx.inputs import transfer

    assert web.raw_destination is transfer.destination_for
    _saved_workspace(workspace)
    import gzip

    outside = tmp_path_factory.mktemp("disk")
    reads = gzip.compress(b"@r1\nACGT\n+\nIIII\n")
    source = outside / "S1_R1.fastq.gz"
    source.write_bytes(reads)
    saved = transfer.import_file(workspace, source)  # default role: the reads folder
    assert saved["path"] == "inputs/fastq/S1_R1.fastq.gz" and saved["bytes"] == len(reads)
    assert (workspace / saved["path"]).read_bytes() == reads and source.exists()
    with pytest.raises(ValueError, match="refusing overwrite"):
        transfer.import_file(workspace, source)
    # Content decides, not the name: reads under any name are accepted, notes named .fastq.gz are not.
    renamed = outside / "anything.reads.gz"
    renamed.write_bytes(reads)
    assert transfer.import_file(workspace, renamed)["path"] == "inputs/fastq/anything.reads.gz"
    wrong = outside / "notes.fastq.gz"
    wrong.write_bytes(gzip.compress(b"meeting notes"))
    with pytest.raises(ValueError, match="this input needs FASTQ reads"):
        transfer.import_file(workspace, wrong)
    unnamed_gzip = outside / "S2_R1.fastq"
    unnamed_gzip.write_bytes(reads)
    with pytest.raises(ValueError, match="does not end in .gz"):
        transfer.import_file(workspace, unnamed_gzip)
    with pytest.raises(ValueError, match="unavailable for the saved configuration"):
        transfer.import_file(workspace, source, role="bam")  # FASTQ study has no such folder
    with pytest.raises(ValueError, match="unavailable for the saved configuration"):
        transfer.download(workspace, "https://example.org/x.bam", None, "bam", "x.bam")
    reference = outside / "Pf3D7.fasta"
    reference.write_text(">c\nACGT\n")
    assert transfer.import_file(workspace, reference, role="reference")["path"] == "inputs/reference/Pf3D7.fasta"


def test_alignment_content_follows_the_saved_input_type(workspace: Path, tmp_path_factory: pytest.TempPathFactory) -> None:
    from pfvarx.inputs import transfer

    _saved_workspace(workspace, "BAM")
    assert transfer.folder_kinds(workspace)["bam"] == ("BAM", "ALIGNMENT_INDEX")
    assert transfer.folder_contents(workspace)["bam"][0] == "BAM alignment"
    cram = tmp_path_factory.mktemp("disk") / "S1.bam"  # CRAM content under a .bam name
    cram.write_bytes(b"CRAM\x03\x00" + b"\x00" * 20)
    with pytest.raises(ValueError, match="contains CRAM alignment; this input needs BAM alignment"):
        transfer.import_file(workspace, cram)


def test_cli_stage_commands_match_the_gui_stages(workspace: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Every GUI tab has a CLI command, and both refuse the placeholder release label."""
    from typer.testing import CliRunner

    from pfvarx.cli import app
    from pfvarx.configuration import web

    commands = {c.name or c.callback.__name__.replace("_", "-") for c in app.registered_commands}
    assert {"input-directories", "input-files", "add-input", "fetch-input", "validate-input-files",
            "prepare-reference", "inputs-preparation", "prepare-samples", "workflow-readiness-checks", "workflow-check", "workflow-run"} <= commands
    monkeypatch.chdir(workspace)
    runner = CliRunner()
    assert runner.invoke(app, ["config-init", "--purpose", "analysis"]).exit_code == 0
    refused = runner.invoke(app, ["input-directories"])
    assert refused.exit_code == 2 and "release label" in refused.stdout
    assert not (workspace / "inputs").exists()
    state = editor_state(workspace, "analysis")
    with pytest.raises(ValueError, match="release label"):
        web.check(workspace, state["values"])
