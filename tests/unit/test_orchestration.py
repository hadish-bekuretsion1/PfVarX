"""Public contracts for selective exports, registration, verification and safe retention."""

import json
import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_sample_export_omits_population_implementation(tmp_path):
    bootstrap = runpy.run_path(str(ROOT / "bootstrap.py"))
    result = bootstrap["export_workspace"](ROOT, tmp_path / "sample", "sample", "WGS", [], False)
    ws = tmp_path / "sample"
    assert result["environment"] == "core"
    assert (ws / "src/pfvarx/caller/call.py").is_file()
    assert (ws / "src/pfvarx/inputs/registry.py").is_file()
    assert not (ws / "src/pfvarx/downstream").exists()
    assert not (ws / "workflow/rules/l5_downstream.smk").exists()
    assert not (ws / "workflow/scripts/da001_moire.R").exists()
    assert not (ws / "inputs").exists()
    with pytest.raises(ValueError, match="exists"):
        bootstrap["export_workspace"](ROOT, ws, "sample", "WGS", [], False)


def test_selected_analysis_export_includes_prerequisites_only(tmp_path):
    bootstrap = runpy.run_path(str(ROOT / "bootstrap.py"))
    result = bootstrap["export_workspace"](ROOT, tmp_path / "selected", "downstream", "WGS", ["DA005"], False)
    assert set(result["analyses_available"]) == {"DA001", "DA005"}
    assert (tmp_path / "selected/src/pfvarx/downstream/structure.py").is_file()
    assert not (tmp_path / "selected/src/pfvarx/downstream/phylogeny.py").exists()


def test_registration_links_without_copy_and_refuses_overwrite(workspace, tmp_path):
    from pfvarx.initialization.layout import config_init
    from pfvarx.inputs.organize import register

    config_init(workspace, "analysis")
    r1 = tmp_path / "S_R1.fastq.gz"
    r2 = tmp_path / "S_R2.fastq.gz"
    r1.write_bytes(b"original1")
    r2.write_bytes(b"original2")
    sheet = tmp_path / "register.tsv"
    sheet.write_text(f"sample_id\tlane\tlibrary_id\tfastq_1\tfastq_2\nS\tL1\tLIB\t{r1}\t{r2}\n")
    result = register(workspace, sheet, "link")
    assert result["samples"] == 1
    assert len(list((workspace / "inputs/fastq").glob("*"))) == 2
    assert all(p.is_symlink() for p in (workspace / "inputs/fastq").iterdir())
    with pytest.raises(ValueError, match="exists"):
        register(workspace, sheet, "link")
    assert r1.read_bytes() == b"original1"


def test_registration_rejects_path_traversal_without_partial_sheet(workspace, tmp_path):
    from pfvarx.initialization.layout import config_init
    from pfvarx.inputs.organize import register

    config_init(workspace, "analysis")
    sheet = tmp_path / "bad.tsv"
    sheet.write_text("sample_id\tlane\tlibrary_id\tfastq_1\tfastq_2\n../escape\tL1\tLIB\tmissing\tmissing\n")
    with pytest.raises(ValueError):
        register(workspace, sheet, "link")
    assert not (workspace / "inputs/metadata/samples.tsv").exists()


def test_cleanup_preview_and_apply_preserve_scientific_and_original_files(tmp_path):
    from pfvarx.execution.retention import cleanup_run

    root = tmp_path
    run = root / "runs/test"
    work = run / "results/work"
    work.mkdir(parents=True)
    (run / "status.json").write_text(json.dumps({"status": "COMPLETED"}))
    original = root / "original.fastq"
    original.write_text("original")
    (work / "temp.bam").write_text("temporary")
    (work / "external-link").symlink_to(original)
    final = run / "results/final.vcf"
    final.write_text("science")
    (run / "frozen_config.json").write_text(json.dumps({"config": {"reference": {}}}))
    (run / "frozen_plan.json").write_text(json.dumps({"samples": []}))
    preview = cleanup_run(root, run, False)
    assert (work / "temp.bam").exists() and preview["bytes_eligible"] == 9
    cleanup_run(root, run, True)
    assert not (work / "temp.bam").exists()
    assert original.read_text() == "original" and final.read_text() == "science"
    assert (work / "external-link").is_symlink()


def test_cleanup_refuses_failed_run_and_active_execution(tmp_path):
    from pfvarx.execution.retention import cleanup_run, workspace_lock

    run = tmp_path / "runs/test"
    run.mkdir(parents=True)
    (run / "status.json").write_text(json.dumps({"status": "FAILED_EXECUTION"}))
    with pytest.raises(ValueError, match="COMPLETED"):
        cleanup_run(tmp_path, run, True)
    with workspace_lock(tmp_path):
        with pytest.raises(ValueError, match="active"):
            with workspace_lock(tmp_path):
                pass


def test_discovery_requires_complete_pairs_and_explicit_registration(tmp_path):
    from pfvarx.inputs.organize import discover

    source = tmp_path / "reads"
    source.mkdir()
    root = tmp_path / "workspace"
    root.mkdir()
    (source / "S_L001_R1_001.fastq.gz").touch()
    with pytest.raises(ValueError, match="Unmatched"):
        discover(root, source)
    assert not (root / "metadata/discovered_samples.tsv").exists()
    (source / "S_L001_R2_001.fastq.gz").touch()
    result = discover(root, source)
    assert result["status"] == "REVIEW_REQUIRED" and result["samples"] == 1
    assert not (root / "inputs/metadata/samples.tsv").exists()


def test_fetch_checks_hash_and_never_overwrites(tmp_path, monkeypatch):
    import hashlib
    import io

    from pfvarx.inputs.organize import fetch

    class Response(io.BytesIO):
        def geturl(self):
            return "https://example.org/file"

    class Opener:
        def open(self, *args, **kwargs):
            return Response(b"verified bytes")

    monkeypatch.setattr("urllib.request.build_opener", lambda *args: Opener())
    digest = hashlib.sha256(b"verified bytes").hexdigest()
    with pytest.raises(ValueError, match="checksum"):
        fetch(tmp_path, "https://example.org/file", "0" * 64, "reference", "ref.fa")
    assert not (tmp_path / "inputs/reference/ref.fa").exists()
    assert not list((tmp_path / "inputs/reference").glob(".download-*"))
    fetch(tmp_path, "https://example.org/file", digest, "reference", "ref.fa")
    unverified = fetch(tmp_path, "https://example.org/file", None, "reference", "unverified.fa")
    assert unverified["status"] == "DOWNLOADED_UNVERIFIED"
    assert unverified["sha256"] == digest
    with pytest.raises(ValueError, match="exists"):
        fetch(tmp_path, "https://example.org/file", digest, "reference", "ref.fa")
    with pytest.raises(ValueError, match="HTTPS"):
        fetch(tmp_path, "http://example.org/file", digest, "reference", "other.fa")


def test_prepared_verification_rejects_corruption_and_dictionary_mismatch(workspace):
    import pysam
    import yaml

    from pfvarx.configuration.loader import load_bundle
    from pfvarx.initialization.layout import config_init
    from pfvarx.io import sha256
    from pfvarx.references.prepare import plan_preparation, required_artifacts, verify_prepared

    config_init(workspace, "analysis")
    config = workspace / "config/config.yaml"
    cfg = yaml.safe_load(config.read_text())
    cfg["input_type"] = "BAM"
    config.write_text(yaml.safe_dump(cfg))
    bundle = load_bundle(workspace)
    source = workspace / bundle.config.reference.fasta
    source.parent.mkdir(parents=True)
    source.write_text(">chr1\nACGTACGT\n")
    prep = workspace / "prepared/reference"
    prep.mkdir(parents=True)
    for name in ["combined", "parasite"]:
        p = prep / (name + ".fasta")
        p.write_text(source.read_text())
        pysam.faidx(str(p))
    dictionary = prep / "combined.dict"
    dictionary.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n")
    artifacts = required_artifacts(workspace, bundle)
    assert len(artifacts) == 5  # aligned entry requires no BWA index
    manifest = {
        "schema": "pfvarx.prepared/2",
        "source_fingerprints": {str(source.resolve()): sha256(source)},
        "artifacts": {str(p.relative_to(workspace)): sha256(p) for p in artifacts},
    }
    mp = workspace / "prepared/manifest.json"
    mp.write_text(json.dumps(manifest))
    assert verify_prepared(workspace, bundle)["status"] == "PASS"
    dictionary.write_text("@SQ\tSN:chr1\tLN:99\n")
    assert verify_prepared(workspace, bundle)["status"] == "FAIL"
    assert any(a["action"] == "picard_dict" for a in plan_preparation(workspace, bundle)["actions"])
    manifest["artifacts"]["prepared/reference/combined.dict"] = sha256(dictionary)
    mp.write_text(json.dumps(manifest))
    result = verify_prepared(workspace, bundle)
    assert result["status"] == "FAIL" and any("disagree" in s for s in result["issues"])


def test_profile_configuration_rejects_unavailable_outputs(tmp_path):
    bootstrap = runpy.run_path(str(ROOT / "bootstrap.py"))
    bootstrap["export_workspace"](ROOT, tmp_path / "sample", "sample", "WGS", [], False)
    import yaml

    from pfvarx.configuration.loader import ConfigError, load_bundle
    from pfvarx.initialization.workspace import initialize

    root = tmp_path / "sample"
    initialize(root)
    path = root / "config/config.yaml"
    cfg = yaml.safe_load(path.read_text())
    cfg["final_output"] = "FO003"
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ConfigError, match="exceeds"):
        load_bundle(root)


def test_keep_policy_changes_snakemake_retention_only(workspace):
    from pfvarx.execution.runner import snakemake_cmd

    eff = {
        "compute": {
            "max_cores": 1,
            "max_jobs": 1,
            "max_mem_mb": 8000,
            "per_job": {"max_mem_mb": 8000, "max_runtime_min": 720},
        },
        "config": {"retention": {"intermediates": "keep"}},
    }
    assert "--notemp" in snakemake_cmd(workspace, workspace / "run", eff, False)
    eff["config"]["retention"]["intermediates"] = "successful"
    assert "--notemp" not in snakemake_cmd(workspace, workspace / "run", eff, False)


def test_diversity_export_does_not_install_moire(tmp_path):
    bootstrap = runpy.run_path(str(ROOT / "bootstrap.py"))
    result = bootstrap["export_workspace"](ROOT, tmp_path / "diversity", "downstream", "WGS", ["DA003"], False)
    assert (tmp_path / "diversity/src/pfvarx/downstream/coi.py").is_file()
    assert result["analyses_available"] == ["DA003"]
    assert "coi" not in result["environment"]


def test_aligned_profile_skips_read_tool_feature_and_creates_no_fastq_directory(tmp_path):
    bootstrap = runpy.run_path(str(ROOT / "bootstrap.py"))
    result = bootstrap["export_workspace"](ROOT, tmp_path / "aligned", "sample", "WGS", [], False, "BAM")
    from pfvarx.initialization.workspace import initialize

    initialize(tmp_path / "aligned")
    assert result["environment"] == "aligned"
    assert (tmp_path / "aligned/inputs/alignments").is_dir()
    assert not (tmp_path / "aligned/inputs/fastq").exists()


def test_cram_validation_never_indexes_original_reference(tmp_path):
    import hashlib

    import pysam

    from pfvarx.inputs.validate import validate_alignment

    reference = tmp_path / "original.fasta"
    reference.write_text(">chr1\nACGTACGTACGT\n")
    pysam.faidx(str(reference))
    header = {
        "HD": {"SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 12}],
        "RG": [{"ID": "rg", "SM": "S", "LB": "L", "PL": "ILLUMINA"}],
    }
    cram = tmp_path / "input.cram"
    with pysam.AlignmentFile(str(cram), "wc", header=header, reference_filename=str(reference)) as handle:
        r = pysam.AlignedSegment()
        r.query_name = "read"
        r.query_sequence = "ACGT"
        r.flag = 0
        r.reference_id = 0
        r.reference_start = 0
        r.mapping_quality = 60
        r.cigarstring = "4M"
        r.query_qualities = pysam.qualitystring_to_array("IIII")
        r.set_tag("RG", "rg")
        handle.write(r)
    pysam.index(str(cram))
    Path(str(reference) + ".fai").unlink()
    result = validate_alignment(
        cram, {"chr1": {"length": 12, "md5": hashlib.md5(b"ACGTACGTACGT").hexdigest()}}, reference
    )
    assert result.status != "FAIL"
    assert not Path(str(reference) + ".fai").exists()


@pytest.mark.parametrize(
    "analysis,module",
    [
        ("DA002", "resistance"),
        ("DA003", "diversity"),
        ("DA001", "coi"),
        ("DA004", "relatedness"),
        ("DA005", "structure"),
        ("DA006", "phylogeny"),
        ("DA007", "spatial"),
    ],
)
def test_single_analysis_export_has_importable_transitive_python_closure(tmp_path, analysis, module):
    import os
    import subprocess
    import sys

    bootstrap = runpy.run_path(str(ROOT / "bootstrap.py"))
    ws = tmp_path / analysis
    bootstrap["export_workspace"](ROOT, ws, "downstream", "WGS", [analysis], False)
    env = dict(os.environ, PYTHONPATH=str(ws / "src"), PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run(
        [sys.executable, "-c", f"import pfvarx.downstream.{module}"], cwd=ws, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_sparse_git_source_export_uses_selected_tracked_blobs(tmp_path):
    import shutil
    import subprocess

    source = tmp_path / "source"
    shutil.copytree(ROOT, source, ignore=shutil.ignore_patterns(".git", ".pixi", "__pycache__", ".ruff_cache"))
    subprocess.run(["git", "init", "-q", str(source)], check=True)
    subprocess.run(["git", "-C", str(source), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(source),
            "-c",
            "user.name=Test",
            "-c",
            "user.email=test@example.invalid",
            "commit",
            "-qm",
            "fixture",
        ],
        check=True,
    )
    subprocess.run(["git", "-C", str(source), "sparse-checkout", "init", "--cone"], check=True)
    subprocess.run(["git", "-C", str(source), "sparse-checkout", "set", "--cone"], check=True)
    assert not (source / "src/pfvarx").exists()
    bootstrap = runpy.run_path(str(source / "bootstrap.py"))
    bootstrap["export_workspace"](source, tmp_path / "selected", "sample", "WGS", [], False)
    assert (tmp_path / "selected/src/pfvarx/inputs/registry.py").is_file()
    assert not (tmp_path / "selected/src/pfvarx/downstream").exists()


def test_aligned_resource_plan_contains_actual_entry_and_qc_not_fastq_jobs():
    from pfvarx.constants import Design, InputType
    from pfvarx.planning.plan import jobs_for

    jobs = jobs_for([0, 2, 6], Design.WGS, InputType.BAM)
    assert {"alignment_entry", "metrics", "pfvarx_call"} <= set(jobs)
    assert not {"fastp", "bwa_mem2_align", "fgbio_umi"} & set(jobs)


def test_preparation_failure_retains_real_tool_diagnostics(tmp_path):
    import sys

    from pfvarx.references.prepare import _run_preparation

    with pytest.raises(ValueError, match="preparation_logs"):
        _run_preparation(
            tmp_path, "index", [sys.executable, "-c", "import sys; print('index failure detail'); sys.exit(7)"]
        )
    log = tmp_path / "metadata/preparation_logs/index.log"
    assert "index failure detail" in log.read_text() and "Exit status: 7" in log.read_text()
    _run_preparation(tmp_path, "index", [sys.executable, "-c", "print('repaired')"])
    assert "index failure detail" in log.read_text() and "repaired" in log.read_text()


def test_cleanup_protects_registered_metadata_even_inside_work(tmp_path):
    from pfvarx.execution.retention import cleanup_run

    run = tmp_path / "runs/test"
    work = run / "results/work"
    work.mkdir(parents=True)
    original = work / "samples.tsv"
    original.write_text("original metadata")
    (run / "status.json").write_text(json.dumps({"status": "COMPLETED"}))
    (run / "frozen_config.json").write_text(json.dumps({"config": {"reference": {}, "samples_tsv": str(original)}}))
    (run / "frozen_plan.json").write_text(json.dumps({"samples": []}))
    result = cleanup_run(tmp_path, run, True)
    assert original.read_text() == "original metadata" and not result["deleted"]
