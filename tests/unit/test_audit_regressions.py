import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from pfvarx.caller_adapter import EVIDENCE_COLUMNS, validate_outputs
from pfvarx.downstream import run
from pfvarx.downstream.common import Context
from pfvarx.io import read_tsv


def test_independent_method_failure_records_terminal_and_releases_coi_checkpoint(
    workspace: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = workspace / "runs/test"
    ctx = Context(root=workspace, run=run_dir, eff={"config": {
        "run_policy": "CONTINUE_INDEPENDENT", "project_id": "TEST",
    }}, plan={"samples": [{"sample_id": "S1"}]}, params={}, samples=["S1"], metadata={})

    def fail(*args: object, **kwargs: object) -> None:
        raise SystemExit("MOIRE failed in test")

    import importlib
    monkeypatch.setattr(importlib, "import_module", lambda name: SimpleNamespace(run=fail))
    run.run_da(ctx, "DA001", 1)
    manifest = json.loads((ctx.results / "downstream/DA001/status_manifest.json").read_text())
    assert manifest["terminal_state"] == "FAILED_EXECUTION"
    assert manifest["samples"]["analysed"] == []
    assert manifest["outputs"] == []
    assert (ctx.results / "downstream/DA001/clonality.tsv").read_text() == "sample_id\tclonality\n"
    ctx.eff["config"]["run_policy"] = "STRICT"
    with pytest.raises(SystemExit, match="MOIRE failed"):
        run.run_da(ctx, "DA001", 1)


def test_external_vcf_rejects_unusable_evidence_values(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    pysam.faidx(str(ref))
    paths = {name: tmp_path / name for name in ("vcf", "evidence", "callable", "focal", "qc")}
    headers = {"GT": ("1", "String"), "AD": ("R", "Integer"), "DP": ("1", "Integer"),
               "PST": ("1", "String"), "WSAF": ("1", "Float"), "WSAFLO": ("1", "Float"),
               "WSAFHI": ("1", "Float")}
    prefix = '##fileformat=VCFv4.3\n##contig=<ID=chr1,length=4>\n'
    prefix += "".join(f'##FORMAT=<ID={key},Number={n},Type={kind},Description="{key}">\n'
                      for key, (n, kind) in headers.items())
    prefix += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
    record = 'chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:AD:DP:PST:WSAF:WSAFLO:WSAFHI\t'
    paths["evidence"].write_text("\t".join(EVIDENCE_COLUMNS) + "\n")
    paths["callable"].write_text("chr1\t0\t4\n")
    paths["focal"].write_text("")
    paths["qc"].write_text(json.dumps({"schema": "pfvarx.caller_qc/1", "sample": "S1",
                                      "mixture": {"status": "UNAVAILABLE"}, "status_counts": {}}))
    paths["vcf"].write_text(prefix + record + '0:9,1:10:PASS:0.1:0.0:0.4\n')
    validate_outputs(paths, "S1", ref)
    paths["vcf"].write_text(prefix + record + '0:9,1:10:PASS:.:0.0:0.4\n')
    with pytest.raises(ValueError, match="allele-fraction"):
        validate_outputs(paths, "S1", ref)
    paths["vcf"].write_text(prefix + record + '0:9,1:8:PASS:0.1:0.0:0.4\n')
    with pytest.raises(ValueError, match="AD and DP"):
        validate_outputs(paths, "S1", ref)


def test_f1_is_zero_when_truth_has_variants_but_no_calls(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pfvarx.benchmarking import evaluate as benchmark

    truth = tmp_path / "inputs/truth"
    truth.mkdir(parents=True)
    (truth / "truth.json").write_text(json.dumps({
        "clones": {}, "samples": {"S1": {"clonality": "MONOCLONAL"}},
    }))
    (truth / "confident_regions.bed").write_text("chr1\t0\t4\n")
    run_dir = tmp_path / "run"
    variants = run_dir / "results/variants/samples"
    variants.mkdir(parents=True)
    (run_dir / "frozen_plan.json").write_text(json.dumps({"reference": {"parasite_fasta": "ref.fa"}}))
    header = '##fileformat=VCFv4.3\n##contig=<ID=chr1,length=4>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    (truth / "S1.truth.vcf").write_text(header + 'chr1\t1\t.\tA\tC\t.\tPASS\t.\n')
    (variants / "S1.vcf.gz").write_bytes(gzip.compress(header.encode()))
    (variants / "S1.callable.bed").write_text("chr1\t0\t4\n")
    monkeypatch.setattr(benchmark, "load_partition", lambda root: {"benchmarking": ["S1"], "sha256": "test"})
    monkeypatch.setattr(benchmark, "_norm", lambda vcf, ref, output: vcf)
    result = benchmark.evaluate(tmp_path, run_dir, tmp_path / "metrics")
    assert result["accuracy"]["SNV"]["f1"] == 0.0
    assert result["accuracy"]["INDEL"]["f1"] is None


@pytest.mark.parametrize(("fit_text", "expected"), [
    ("marker_id\tt\tfit\n", "INSUFFICIENT_MODEL_SUPPORT"),
    ("marker_id\tt\tfit\nM1\t0\t0.25\n", "FITTED"),
])
def test_spatial_manifest_only_claims_a_fit_with_rows(
    workspace: Path, monkeypatch: pytest.MonkeyPatch, fit_text: str, expected: str
) -> None:
    from pfvarx.downstream import spatial

    samples = [f"S{i}" for i in range(10)]
    ctx = Context(
        root=workspace, run=workspace / "runs/test",
        eff={"config": {"project_id": "TEST", "analyses": {"spatial_outcome_source": "DA002"}}},
        plan={"samples": [{"sample_id": sid} for sid in samples]},
        params={"da007_min_samples_per_unit": 1}, samples=samples,
        metadata={sid: {"country": "TEST", "latitude": "0", "longitude": "0",
                        "collection_date": "2020-01-01"} for sid in samples},
    )
    markers = ctx.out("DA002")
    (markers / "status_manifest.json").write_text('{"terminal_state":"COMPLETED"}')
    (markers / "marker_calls.tsv").write_text(
        "sample_id\tmarker_id\tstate\n" + "".join(f"{sid}\tM1\tMUTANT\n" for sid in samples)
    )
    monkeypatch.setattr(spatial, "in_country", lambda *args: True)
    monkeypatch.setattr(spatial.shutil, "which", lambda name: "/test/Rscript")

    def fake_r(command: list[str], **kwargs: object) -> SimpleNamespace:
        Path(command[3]).write_text(fit_text)  # Rscript script input OUTPUT [settings...]
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(spatial.subprocess, "run", fake_r)
    spatial.run(ctx)
    manifest = json.loads((ctx.results / "downstream/DA007/status_manifest.json").read_text())
    assert manifest["method"]["model"] == expected
    assert any(item["path"].endswith("/gam_fit.tsv") for item in manifest["outputs"]) == (expected == "FITTED")


def test_coi_withholds_summary_below_informative_locus_gate(
    workspace: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pfvarx.downstream import coi

    ctx = Context(
        root=workspace, run=workspace / "runs/test", eff={"config": {"project_id": "TEST"}},
        plan={"samples": [{"sample_id": "S1"}]}, samples=["S1"], metadata={},
        params={"da001_thin_bp": 500, "da001_allele_presence_floor": 0.05,
                "da001_credibility": 0.95, "da001_min_informative_loci": 2, "population_frequency_min_samples": 10},
        sites=[{"type": "SNV", "n_pass": "1", "site_id": "SITE", "chrom": "chr1", "pos": "1"}],
        matrix={"SITE": {"S1": {"wsaf": 0.4}}},
    )
    gate = ctx.results / "qc/samples/S1/qc_decision.json"
    gate.parent.mkdir(parents=True)
    gate.write_text('{"permitted":{"DA001":true}}')

    def unexpected_fit(*args: object, **kwargs: object) -> None:
        pytest.fail("insufficient loci must stop before fitting")

    monkeypatch.setattr(coi.subprocess, "run", unexpected_fit)
    coi.run(ctx)
    out = ctx.results / "downstream/DA001"
    manifest = json.loads((out / "status_manifest.json").read_text())
    assert manifest["terminal_state"] == "SKIPPED_INELIGIBLE"
    assert manifest["samples"]["excluded"]["INSUFFICIENT_INFORMATIVE_SITES"] == ["S1"]
    assert manifest["denominators"]["analysed"] == 0
    assert not (out / "moire_coi.tsv").exists()


def test_reference_fingerprint_survives_moving_the_workspace(tmp_path):
    from pfvarx.references.prepare import _source_key

    first, second = tmp_path / "a", tmp_path / "b"
    for root in (first, second):
        (root / "inputs/reference").mkdir(parents=True)
        (root / "inputs/reference/ref.fasta").write_text(">c\nACGT\n")
    relative = Path("inputs/reference/ref.fasta")
    assert _source_key(first, relative) == _source_key(second, relative) == "inputs/reference/ref.fasta"
    # Manifests written by earlier releases recorded the absolute path of the same file.
    assert _source_key(first, (first / relative).resolve()) == "inputs/reference/ref.fasta"
    outside = tmp_path / "shared.fasta"
    outside.write_text(">c\nACGT\n")
    assert _source_key(first, outside) == str(outside.resolve())


def test_snpeff_configuration_is_relocatable_and_uses_the_mitochondrial_code(tmp_path):
    from pfvarx.annotation import snpeff_db

    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">Pf3D7_01_v3 | organism=Plasmodium\nACGT\n>Pf3D7_MIT_v3\nACGT\n>Pf3D7_API_v3\nACGT\n")
    assert snpeff_db.mitochondrial_contigs(fasta) == ["Pf3D7_MIT_v3"]
    text = snpeff_db.write_config(tmp_path, "pf", ["Pf3D7_MIT_v3"]).read_text()
    assert "data.dir = ./data/" in text and str(tmp_path) not in text
    assert "pf.Pf3D7_MIT_v3.codonTable : Protozoan_Mitochondrial" in text
    codons = dict(c.strip().split("/") for c in snpeff_db.MITOCHONDRIAL_CODONS.split(","))
    assert len(codons) == 64 and codons["TGA"] == "W" and codons["TAA"] == "*"
    assert "codonTable" not in snpeff_db.write_config(tmp_path, "pf", []).read_text()


def test_spatial_boundaries_are_read_from_the_workspace():
    from pfvarx.downstream import spatial

    assert not spatial.GEO.is_absolute()


def test_who_k13_flag_needs_enough_callable_infections() -> None:
    from pfvarx.downstream.resistance import prevalence

    marker = {"marker_id": "k13_C580Y", "gene_name": "kelch13", "gene_id": "PF3D7_1343700",
              "evidence_class": "WHO validated"}

    def calls(states: list[str], group: str) -> list[dict[str, str]]:
        return [{"marker_id": "k13_C580Y", "population": group, "state": s, "alt_fraction": ""} for s in states]

    rows = prevalence(calls(["MUTANT"], "A") + calls(["MUTANT", "WILD_TYPE", "WILD_TYPE", "UNKNOWN"], "B"),
                      [marker], "population", 0.05, 3)
    flag = {r["group"]: (r["callable"], r["prevalence_infections"], r["who_k13_alert"]) for r in rows}
    assert flag["A"] == (1, "1.0000", "INSUFFICIENT_CALLABLE")  # one sample never raises the WHO flag
    assert flag["B"] == (3, "0.3333", "ABOVE_5_PERCENT")
    assert flag["ALL"] == (4, "0.5000", "ABOVE_5_PERCENT")


def test_small_cohort_uses_the_caller_mixture_instead_of_moire(
    workspace: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Below the population minimum MOIRE is not identified (on a 9-sample truth cohort its
    P(COI>1) for single-clone samples ranged 0.00-1.00 across seeds), so it is not run: the
    call comes from the caller's mixture test and the COI count and FWS are withheld."""
    from pfvarx.downstream import coi

    sites = [{"type": "SNV", "n_pass": "1", "site_id": f"S{i}", "chrom": "chr1", "pos": str(1000 * i)} for i in range(1, 4)]
    samples = ["A", "B", "C"]
    ctx = Context(
        root=workspace, run=workspace / "runs/test", eff={"config": {"project_id": "TEST", "random_seed": 1}},
        plan={"samples": [{"sample_id": s} for s in samples]}, samples=samples, metadata={},
        params={"da001_thin_bp": 500, "da001_allele_presence_floor": 0.05, "da001_credibility": 0.95,
                "da001_min_informative_loci": 2, "population_frequency_min_samples": 10,
                "da001_moire_burnin": 10, "da001_moire_samples": 10, "da001_fws_reference_cutoff": 0.95},
        sites=sites, matrix={s["site_id"]: {sid: {"wsaf": 0.4} for sid in samples} for s in sites},
    )
    mixtures = {"A": {"status": "ESTIMATED", "k_selected": 1}, "B": {"status": "ESTIMATED", "k_selected": 3},
                "C": {"status": "INSUFFICIENT_SITES", "k_selected": None}}
    for sid in samples:
        d = ctx.results / f"qc/samples/{sid}"
        d.mkdir(parents=True)
        (d / "qc_decision.json").write_text('{"permitted":{"DA001":true}}')
        (d / f"{sid}.caller_qc.json").write_text(json.dumps({"mixture": mixtures[sid]}))

    def no_moire(*args: object, **kwargs: object) -> None:
        pytest.fail("MOIRE must not run below the population minimum")

    monkeypatch.setattr(coi.subprocess, "run", no_moire)
    coi.run(ctx)
    out = ctx.results / "downstream/DA001"
    rows = {r["sample_id"]: r for r in read_tsv(out / "clonality.tsv")}
    assert (rows["A"]["clonality"], rows["A"]["basis"]) == ("MONOCLONAL", "CALLER_MIXTURE_LRT")
    assert (rows["B"]["clonality"], rows["B"]["basis"]) == ("POLYCLONAL", "CALLER_MIXTURE_LRT")
    assert (rows["C"]["clonality"], rows["C"]["basis"]) == ("INDETERMINATE", "")
    assert all(r["prob_polyclonal"] == r["coi_median"] == r["fws"] == "" for r in rows.values())
    manifest = json.loads((out / "status_manifest.json").read_text())
    assert manifest["method"]["moire_run"] is False and manifest["samples"]["analysed"] == ["A", "B"]
    assert not (out / "moire_coi.tsv").exists()
