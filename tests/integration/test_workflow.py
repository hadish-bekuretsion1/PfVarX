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


def _ws(tmp: Path, design: str, fo: str, depth: int = 15, fixture_args: tuple[str, ...] = ()) -> Path:
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
    assert run("fixture", "--design", design, "--depth", str(depth), *fixture_args).returncode == 0
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
    # A release label with spaces once broke the SnpEff database ID; keep one here.
    cfg["reference"]["release_label"] = "Pf3D7 - release unknown"
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


def _umi_to_resistance(tmp_path: Path, depth: int) -> tuple[Path, dict]:
    import json

    if shutil.which("fgbio") is None or shutil.which("snpEff") is None:
        pytest.skip("fgbio or snpEff not on PATH (targeted/annotation features not installed)")
    ws = _ws(tmp_path, "TARGETED_UMI", "FO002", depth=depth)
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["final_output"] = "FO004"
    cfg["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    cfg["analyses"].update(requested=["DA002"], resistance_catalog="inputs/metadata/test_resistance_catalog.tsv")
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    sheet = ws / "inputs/metadata/samples.tsv"
    sheet.write_text("\n".join(sheet.read_text().splitlines()[:4]) + "\n")
    r = subprocess.run(["pfvarx", "pipeline"], cwd=ws, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    result = json.loads(r.stdout)
    return ws / "runs" / result["run_id"] / "results", result


def test_no_sample_passing_qc_skips_cohort_instead_of_crashing(tmp_path: Path) -> None:
    """Shallow UMI families fail the amplicon-dropout gate for every sample: the run must
    finish with an explicit empty-cohort outcome, not a bcftools merge on no files."""
    import json

    res, result = _umi_to_resistance(tmp_path, depth=15)
    assert result["status"] == "COMPLETED", result
    for gate in (res / "qc/samples").glob("*/qc_decision.json"):
        assert json.loads(gate.read_text())["verdict"] == "FAIL"
    status = json.loads((res / "cohort/cohort_status.json").read_text())
    assert status["status"] == "SKIPPED_NO_ELIGIBLE_SAMPLES" and status["cohort_samples"] == []
    assert not (res / "cohort/cohort.vcf.gz").exists()
    manifest = json.loads((res / "downstream/DA002/status_manifest.json").read_text())
    assert manifest["terminal_state"] == "SKIPPED_INELIGIBLE"
    assert "NO_ELIGIBLE_SAMPLES" in manifest["reason_codes"]
    assert (res / "report/pfvarx_dashboard.html").is_file()


def test_umi_route_reaches_resistance_analysis(tmp_path: Path) -> None:
    """With deep enough UMI families the targeted UMI route runs through calling, cohort,
    annotation and DA002 (the shallow fixture above never gets past QC)."""
    import json

    res, result = _umi_to_resistance(tmp_path, depth=60)
    assert result["status"] == "COMPLETED", result
    for gate in (res / "qc/samples").glob("*/qc_decision.json"):
        assert json.loads(gate.read_text())["verdict"] != "FAIL"
    assert json.loads((res / "cohort/cohort_status.json").read_text())["status"] == "COMPLETED"
    assert (res / "cohort/cohort.vcf.gz.tbi").is_file() and (res / "annotation/cohort.ann.vcf.gz.tbi").is_file()
    manifest = json.loads((res / "downstream/DA002/status_manifest.json").read_text())
    assert manifest["terminal_state"].startswith("COMPLETED"), manifest


def test_every_downstream_analysis_executes(tmp_path: Path) -> None:
    """A WGS cohort that reaches all seven analyses (DA001-DA007), not just their skip paths.

    The synthetic genome is 48 kb, so the fixture adds variant sites (--variant-scale 4) and
    the test thins MOIRE loci at 100 bp instead of 500 bp: at 500 bp only ~90 loci exist and
    MOIRE cannot separate monoclonal from polyclonal infections. A real 23 Mb genome keeps
    tens of thousands of loci at the default spacing, so only this test copy is changed.
    """
    import json
    import re

    tools = ("snpEff", "Rscript", "hmmIBD", "iqtree2")
    if any(shutil.which(t) is None for t in tools):
        pytest.skip("downstream method tools not on PATH (default environment needed)")
    ws = _ws(tmp_path, "WGS", "FO002", depth=40, fixture_args=("--n-mono", "8", "--variant-scale", "4"))
    layer5 = ws / "config/parameters/layer5.yaml"
    text, n = re.subn(r"^(  da001_thin_bp: \{value: )\d+", r"\g<1>100", layer5.read_text(), flags=re.M)
    assert n == 1
    layer5.write_text(text)
    analyses = ["DA001", "DA002", "DA003", "DA004", "DA005", "DA006", "DA007"]
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["final_output"] = "FO004"
    cfg["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    cfg["analyses"].update(requested=analyses, resistance_catalog="inputs/metadata/test_resistance_catalog.tsv",
                           spatial_outcome_source="DA002", population_column="population")
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    r = subprocess.run(["pfvarx", "pipeline"], cwd=ws, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    result = json.loads(r.stdout)
    assert result["status"] == "COMPLETED", result
    res = ws / "runs" / result["run_id"] / "results"
    states = {da: json.loads((res / f"downstream/{da}/status_manifest.json").read_text())["terminal_state"]
              for da in analyses}
    assert all(s.startswith("COMPLETED") for s in states.values()), states
    # Clonality must not contradict the simulated truth: PF001-PF008 are single clones and
    # PF009-PF011 are mixtures. INDETERMINATE is allowed; a wrong definite call is not.
    calls = {row.split("\t")[0]: row.split("\t")[1]
             for row in (res / "downstream/DA001/clonality.tsv").read_text().splitlines()[1:]}
    truth = {f"PF{i:03d}": "MONOCLONAL" if i <= 8 else "POLYCLONAL" for i in range(1, 12)}
    assert all(calls[s] in (truth[s], "INDETERMINATE") for s in truth), calls
    assert sum(c == "MONOCLONAL" for c in calls.values()) >= 4, calls
    summary = json.loads((res / "downstream/DA004/status_manifest.json").read_text())["results_summary"]
    assert summary["pairs_estimated"] > 0, summary
    assert (res / "downstream/DA006/tree.treefile").is_file()
