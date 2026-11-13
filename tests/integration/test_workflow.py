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
    cfg["reference"]["release_label"] = "synthetic"
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


def test_end_to_end_annotated_run_writes_indexed_bgzf(tmp_path: Path) -> None:
    """A real (small) run through alignment, calling, cohort merge, annotation and one analysis."""
    import json

    if shutil.which("snpEff") is None:
        pytest.skip("snpEff not on PATH (annotation feature not installed)")
    ws = _ws(tmp_path, "WGS", "FO002")
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["final_output"] = "FO004"
    cfg["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    cfg["analyses"].update(requested=["DA002"], resistance_catalog="inputs/metadata/test_resistance_catalog.tsv")
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    # Keep three samples: enough for a cohort, quick enough for continuous integration.
    sheet = ws / "inputs/metadata/samples.tsv"
    sheet.write_text("\n".join(sheet.read_text().splitlines()[:4]) + "\n")
    r = subprocess.run(["pfvarx", "pipeline"], cwd=ws, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    result = json.loads(r.stdout)
    assert result["status"] == "COMPLETED", result
    res = ws / "runs" / result["run_id"] / "results"
    annotated = res / "annotation/cohort.ann.vcf.gz"
    assert annotated.read_bytes()[12:16] == b"BC\x02\x00"  # BGZF extra field, not plain gzip
    assert (res / "annotation/cohort.ann.vcf.gz.tbi").is_file()
    assert (res / "cohort/cohort.vcf.gz.tbi").is_file()
    config = (ws / "prepared/snpeff/snpEff.config").read_text()
    assert str(ws) not in config  # the annotation database must survive moving the workspace
    manifest = json.loads((res / "downstream/DA002/status_manifest.json").read_text())
    assert manifest["terminal_state"].startswith("COMPLETED"), manifest


def test_umi_route_survives_fixed_read_structure_and_short_reads(tmp_path: Path) -> None:
    """Adapter-dimer pairs and a fixed '8M100T' structure used to stop fgbio FastqToBam."""
    import gzip
    import json

    if shutil.which("fgbio") is None:
        pytest.skip("fgbio not on PATH (targeted feature not installed)")
    ws = _ws(tmp_path, "TARGETED_UMI", "FO001")
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["assay"]["read_structure"] = "8M100T 8M100T"
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    sheet = ws / "inputs/metadata/samples.tsv"
    sheet.write_text("\n".join(sheet.read_text().splitlines()[:2]) + "\n")
    for mate in ("R1", "R2"):
        with gzip.open(ws / f"inputs/fastq/PF001_L001_{mate}.fastq.gz", "at") as handle:
            for i, length in enumerate((0, 5, 8, 30)):
                handle.write(f"@short:{i}\n{('ACGT' * 27)[:length]}\n+\n{'I' * length}\n")
    r = subprocess.run(["pfvarx", "pipeline"], cwd=ws, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    result = json.loads(r.stdout)
    assert result["status"] == "COMPLETED", result
    assert (ws / "runs" / result["run_id"] / "results/alignments/PF001.bam").is_file()
