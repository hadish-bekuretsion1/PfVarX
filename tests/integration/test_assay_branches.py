"""Every assay design through FASTQ and BAM entry, driven by ``assay_design`` and the manifest.

A FASTQ run per design (intermediates kept) supplies real alignments at each processing
state; BAM-entry runs then start from them with the documented ``processing_state`` and must
run only the steps that state still needs. The rules that actually ran are read from the
Snakemake log. Conflicting technique labels must be rejected at validation, naming the
sample and library.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest
import yaml

pytestmark = pytest.mark.integration
REPO = Path(__file__).resolve().parents[2]
NEED = ("samtools", "bwa-mem2", "fastp", "snakemake", "picard", "fgbio")
DEPTH = {"wgs": 15, "targeted_no_umi": 60, "targeted_umi": 60}
UMI_RULES = {"umi_fastq_to_bam", "umi_group", "umi_consensus", "umi_realign"}


def pfvarx(ws: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["pfvarx", *args], cwd=ws, capture_output=True, text=True)


def workspace(path: Path, design: str, input_type: str = "FASTQ") -> Path:
    if any(shutil.which(t) is None for t in NEED):
        pytest.skip("bioinformatics tools not on PATH (run inside 'pixi run')")
    path.mkdir(parents=True)
    for name in ("config", "workflow", "src", "pixi.toml", "pixi.lock", "pyproject.toml"):
        source = REPO / name
        (shutil.copytree if source.is_dir() else shutil.copyfile)(source, path / name)
    for name in ("config.yaml", "quality.yaml", "compute.yaml", "profiles/resources.yaml"):
        (path / "config" / name).unlink(missing_ok=True)
    assert pfvarx(path, "config-init", "--purpose", "analysis").returncode == 0
    cfg = yaml.safe_load((path / "config/config.yaml").read_text())
    cfg.update(final_output="FO002", input_type=input_type.lower(), assay_design=design)  # lower case accepted
    cfg["reference"]["release_label"] = "synthetic"
    cfg["retention"] = {"intermediates": "keep"}
    if design != "wgs":
        cfg["assay"]["library_catalog"] = "inputs/metadata/library_catalog.yaml"
    (path / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    assert pfvarx(path, "init").returncode == 0
    r = pfvarx(path, "fixture", "--design", design, "--depth", str(DEPTH[design]), "--n-mono", "2", "--n-poly", "0")
    assert r.returncode == 0, r.stdout + r.stderr
    return path


def run_pipeline(ws: Path) -> tuple[Path, set[str]]:
    for step in ("validate", "prepare"):
        r = pfvarx(ws, step)
        assert r.returncode == 0, f"{step}: {r.stdout[-3000:]}{r.stderr[-2000:]}"
    r = pfvarx(ws, "pipeline")
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    result = json.loads(r.stdout)
    assert result["status"] == "COMPLETED", result
    run = ws / "runs" / result["run_id"]
    rules = set()
    for log in run.glob("logs/snakemake.pass*.log"):
        rules |= set(re.findall(r"^(?:local)?(?:rule|checkpoint) (\w+):", log.read_text(), re.M))
    return run / "results", rules


def assert_calls(res: Path, sample: str, unit: str) -> None:
    vcf = res / f"variants/samples/{sample}.vcf.gz"
    assert vcf.is_file() and (res / f"variants/samples/{sample}.vcf.gz.tbi").is_file()
    header = subprocess.run(["bcftools", "view", "-h", str(vcf)], capture_output=True, text=True, check=True).stdout
    assert f"##PfVarX_count_unit={unit}" in header
    qc = json.loads((res / f"qc/samples/{sample}/{sample}.caller_qc.json").read_text())
    assert qc["count_unit"] == unit
    gate = json.loads((res / f"qc/samples/{sample}/qc_decision.json").read_text())
    assert gate["count_unit"] == unit
    evidence = (res / f"variants/samples/{sample}.evidence.tsv").read_text().splitlines()
    assert len(evidence) > 1 and (res / f"variants/samples/{sample}.callable.bed").stat().st_size > 0


@pytest.fixture(scope="module")
def fastq_runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, tuple[Path, Path, set[str]]]:
    out = {}
    for design in DEPTH:
        ws = workspace(tmp_path_factory.mktemp(design) / "ws", design)
        res, rules = run_pipeline(ws)
        out[design] = (ws, res, rules)
    return out


@pytest.mark.parametrize("design", list(DEPTH))
def test_fastq_entry_runs_the_design_branch(fastq_runs: dict[str, Any], design: str) -> None:
    ws, res, rules = fastq_runs[design]
    unit = "consensus_molecules" if design == "targeted_umi" else "reads"
    for sample in ("PF001", "PF002"):
        assert_calls(res, sample, unit)
    assert ("mark_duplicates" in rules) == (design == "wgs")
    assert ("trim_primers" in rules) == (design == "targeted_no_umi")
    assert (UMI_RULES <= rules) == (design == "targeted_umi")
    library = json.loads((res / "qc/samples/PF001/qc_decision.json").read_text())["library"]
    assert library["amplification_route"] == ("DIRECT" if design == "wgs" else "AMPLICON")
    if design != "wgs":
        assert (library["panel_id"], library["panel_version"]) == ("PFVARX_SYNTHETIC_PANEL", "1")


def source_bam(res: Path, sample: str, kind: str, dest: Path) -> None:
    """Copy (and coordinate-sort and index) one alignment from a FASTQ run."""
    path = {"lane": res / f"work/lanes/{sample}.L001.bam", "merged": res / f"work/merged/{sample}.bam",
            "umi_raw": res / f"work/umi/lanes/{sample}.L001.mapped.bam",
            "consensus": res / f"work/umi/{sample}.consensus.mapped.bam",
            "final": res / f"alignments/{sample}.bam"}[kind]
    assert path.is_file(), path
    subprocess.run(["samtools", "sort", "-o", str(dest), str(path)], check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(dest)], check=True)


CASES = [
    # design, source alignment, documented state, rules that must run, rules that must not
    ("wgs", "lane", "aligned", set(), {"entry_trim_primers"} | UMI_RULES),
    ("wgs", "final", "duplicates_marked", set(), {"entry_trim_primers"} | UMI_RULES),
    ("targeted_no_umi", "merged", "aligned", {"entry_trim_primers"}, UMI_RULES),
    ("targeted_no_umi", "final", "primers_trimmed", set(), {"entry_trim_primers"} | UMI_RULES),
    ("targeted_umi", "umi_raw", "aligned", {"umi_group", "umi_consensus", "umi_realign", "entry_trim_primers"},
     {"umi_fastq_to_bam"}),
    ("targeted_umi", "consensus", "consensus", {"entry_trim_primers"}, UMI_RULES),
    ("targeted_umi", "final", "consensus_primers_trimmed", set(), {"entry_trim_primers"} | UMI_RULES),
]


@pytest.mark.parametrize(("design", "kind", "state", "ran", "skipped"), CASES,
                         ids=[f"{c[0]}-{c[2]}" for c in CASES])
def test_bam_entry_runs_only_the_steps_its_state_needs(fastq_runs: dict[str, Any], tmp_path: Path, design: str,
                                                       kind: str, state: str, ran: set[str], skipped: set[str]) -> None:
    _, source, _ = fastq_runs[design]
    ws = workspace(tmp_path / "ws", design, "BAM")
    sheet = ws / "prepared/samples.tsv"
    rows = [line.split("\t") for line in sheet.read_text().splitlines()]
    header = rows[0]
    (ws / "inputs/bam").mkdir(exist_ok=True)
    for row in rows[1:]:
        sample = row[header.index("sample_id")]
        source_bam(source, sample, kind, ws / f"inputs/bam/{sample}.bam")
        row[header.index("bam")] = f"inputs/bam/{sample}.bam"
        row[header.index("fastq_1")] = row[header.index("fastq_2")] = ""
        row[header.index("processing_state")] = state
    sheet.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    res, rules = run_pipeline(ws)
    assert {"alignment_entry", "entry_finalize"} <= rules
    assert ran <= rules, f"expected {ran - rules} to run"
    assert not (skipped & rules), f"{skipped & rules} should have been skipped"
    unit = "consensus_molecules" if design == "targeted_umi" else "reads"
    for sample in ("PF001", "PF002"):
        assert_calls(res, sample, unit)
        contract = json.loads((res / f"qc/samples/{sample}/{sample}.entry_contract.json").read_text())
        assert contract["processing_state"] == state
        if design == "wgs":
            assert ("duplicates_marked" in contract["actions"]) == (state == "aligned")


def validation_failures(ws: Path) -> list[str]:
    """FAIL details of './pfvarx validate'. A sample-level failure leaves the others admitted
    (report status PARTIAL), so the command's exit code is not the test."""
    pfvarx(ws, "validate")
    report = json.loads((ws / "metadata/validation_report.json").read_text())
    return [f["detail"] for rep in report["reports"] for f in rep["findings"] if f["level"] == "FAIL"]


def set_column(ws: Path, column: str, values: dict[str, str]) -> None:
    sheet = ws / "prepared/samples.tsv"
    rows = [line.split("\t") for line in sheet.read_text().splitlines()]
    for row in rows[1:]:
        sid = row[rows[0].index("sample_id")]
        if sid in values:
            row[rows[0].index(column)] = values[sid]
    sheet.write_text("\n".join("\t".join(r) for r in rows) + "\n")


def test_conflicting_technique_labels_are_rejected(tmp_path: Path) -> None:
    ws = workspace(tmp_path / "no_umi", "targeted_no_umi")
    set_column(ws, "technique", {"PF001": "TBS_UMI", "PF002": "null"})
    failures = validation_failures(ws)
    assert any("sample PF001 (library PF001_lib1)" in f and "TBS_UMI" in f and "targeted_umi" in f for f in failures)
    assert not any("PF002" in f for f in failures)  # null means unspecified
    report = json.loads((ws / "metadata/validation_report.json").read_text())
    assert report["samples"] == {"PF001": "FAIL", "PF002": "PASS"} and report["status"] == "PARTIAL"

    wgs = workspace(tmp_path / "wgs", "wgs")
    set_column(wgs, "technique", {"PF001": "AMPLICON", "PF002": "SELECTIVE_WGS"})
    set_column(wgs, "amplification_route", {"PF001": "", "PF002": ""})
    failures = validation_failures(wgs)
    assert any("sample PF001 (library PF001_lib1)" in f and "AMPLICON" in f for f in failures)
    report = json.loads((wgs / "metadata/validation_report.json").read_text())
    assert report["samples"]["PF002"] == "PASS"
    assert report["libraries"]["PF002"]["amplification_route"] == "SWGA"  # SELECTIVE_WGS records sWGA
    set_column(wgs, "amplification_route", {"PF002": "DIRECT"})
    failures = validation_failures(wgs)
    assert any("sample PF002" in f and "disagree" in f for f in failures)

    umi = workspace(tmp_path / "umi", "targeted_umi")
    set_column(umi, "technique", {"PF001": "TBS_NON_UMI"})
    failures = validation_failures(umi)
    assert any("sample PF001 (library PF001_lib1)" in f and "TBS_NON_UMI" in f for f in failures)


def test_consensus_state_on_raw_umi_reads_is_rejected(fastq_runs: dict[str, Any], tmp_path: Path) -> None:
    """A documented state the reads contradict is refused before anything runs."""
    _, source, _ = fastq_runs["targeted_umi"]
    ws = workspace(tmp_path / "ws", "targeted_umi", "BAM")
    sheet = ws / "prepared/samples.tsv"
    rows = [line.split("\t") for line in sheet.read_text().splitlines()]
    header = rows[0]
    (ws / "inputs/bam").mkdir(exist_ok=True)
    for row in rows[1:]:
        sample = row[header.index("sample_id")]
        source_bam(source, sample, "umi_raw", ws / f"inputs/bam/{sample}.bam")
        row[header.index("bam")] = f"inputs/bam/{sample}.bam"
        row[header.index("processing_state")] = "consensus"
    sheet.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    failures = validation_failures(ws)
    assert any("sample PF001" in f and "consensus" in f and "cD" in f for f in failures)
    report = json.loads((ws / "metadata/validation_report.json").read_text())
    assert report["status"] == "FAIL" and set(report["samples"].values()) == {"FAIL"}
