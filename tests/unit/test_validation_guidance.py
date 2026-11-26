"""Failed validation must say which input, why it failed and how to fix it."""

from pathlib import Path

from pfvarx.configuration.loader import load_parameters
from pfvarx.execution.preflight import summarise
from pfvarx.inputs import validate as v
from pfvarx.inputs.guidance import explain

CONTIGS = {"Pf3D7_01_v3": {"length": 10_000}, "Pf3D7_02_v3": {"length": 10_000}}
LIMITS = v.CoreMaskLimits.from_parameters(load_parameters(Path(__file__).resolve().parents[2]))


def test_missing_core_mask_explains_folder_shortcut_and_lists_candidates(tmp_path: Path) -> None:
    ref = tmp_path / "inputs/reference"
    ref.mkdir(parents=True)
    (ref / "pf3D7_genes.bed").write_text("Pf3D7_01_v3\t0\t10\tPF3D7_0100100\n")
    report = v.validate_bed(ref / "core_genome.bed", CONTIGS, "CORE_MASK", LIMITS).as_dict()
    finding = report["findings"][0]
    advice = explain(finding, "CORE_MASK", tmp_path, {"reference": {"core_mask_bed": "inputs/reference/"}})
    assert advice["file"] == "Core-genome mask BED" and "accessible core" in advice["used_for"]
    assert "names only the folder inputs/reference/" in advice["why"]
    text = " ".join(advice["actions"])
    assert "pf3D7_genes.bed" in text and "Miles et al. 2016" in text and "$4==\"Core\"" in text
    assert advice["actions"][-1].startswith("Then run ./pfvarx validate again")


def test_gene_bed_used_as_core_mask_is_flagged(tmp_path: Path) -> None:
    genes = tmp_path / "core_genome.bed"
    rows = [f"Pf3D7_01_v3\t{i * 4}\t{i * 4 + 6}\tPF3D7_01{i:05d}\t0\t+" for i in range(1200)]
    genes.write_text("\n".join(rows) + "\n")
    rep = v.validate_bed(genes, {"Pf3D7_01_v3": {"length": 10_000}}, "CORE_MASK", LIMITS)
    assert rep.status == "WARNING"
    (warning,) = rep.findings
    assert warning.reason == "CORE_MASK_SUSPECT"
    assert "gene BED" in warning.detail and "1200 intervals" in warning.detail and "overlap" in warning.detail


def test_genuine_core_mask_passes_without_warning(tmp_path: Path) -> None:
    core = tmp_path / "core_genome.bed"
    core.write_text("Pf3D7_01_v3\t500\t9500\nPf3D7_02_v3\t400\t9300\n")
    rep = v.validate_bed(core, CONTIGS, "CORE_MASK", LIMITS)
    assert rep.status == "PASS" and not rep.findings
    # The same content as a targets BED is never judged as a mask.
    assert v.validate_bed(core, CONTIGS, "TARGETS").status == "PASS"


def test_summary_reads_as_problem_reason_and_numbered_steps() -> None:
    finding = {"level": "FAIL", "reason": "FASTQ_TRUNCATED", "path": "a.fq.gz | b.fq.gz",
               "detail": "file ends inside a record", "record": None}
    finding.update(explain(finding, "FASTQ_PAIR", Path("."), {}))
    lines = summarise({"status": "FAIL", "reports": [{"findings": [finding]}]})
    assert lines[0] == "Validation FAIL: 1 problem(s) must be fixed before the pipeline can run."
    assert lines[1].startswith("Problem 1: Paired FASTQ files.") and "file ends inside a record" in lines[1]
    assert lines[3].startswith("  1. The file ends inside a read")
    assert summarise({"status": "PASS", "reports": []}) == ["Validation PASS: inputs are ready for ./pfvarx prepare."]
