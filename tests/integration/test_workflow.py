"""Integration tests on the synthetic fixture (need the pixi environment's tools)."""

from __future__ import annotations

import gzip
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
    cfg["assay_design"] = design
    cfg["reference"]["release_label"] = "synthetic"
    if design != "wgs":
        # The fixture writes the panel and UMI protocol into the library catalog.
        cfg["assay"]["library_catalog"] = "inputs/metadata/library_catalog.yaml"
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    assert run("init").returncode == 0
    assert run("fixture", "--design", design, "--depth", str(depth), *fixture_args).returncode == 0
    assert run("validate").returncode == 0
    assert run("prepare").returncode == 0
    return ws


@pytest.mark.parametrize("design", ["wgs", "targeted_no_umi", "targeted_umi"])
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
    ws = _ws(tmp_path, "wgs", "FO002")
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["final_output"] = "FO004"
    cfg["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    # A release label with spaces once broke the SnpEff database ID; keep one here.
    cfg["reference"]["release_label"] = "Pf3D7 - release unknown"
    cfg["analyses"].update(requested=["DA002"], resistance_catalog="inputs/metadata/test_resistance_catalog.tsv")
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    # Keep three samples: enough for a cohort, quick enough for continuous integration.
    sheet = ws / "prepared/samples.tsv"
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
    ws = _ws(tmp_path, "targeted_umi", "FO001")
    catalog_path = ws / "inputs/metadata/library_catalog.yaml"
    catalog = yaml.safe_load(catalog_path.read_text())
    catalog["library_protocols"][0]["umi"]["read_structure"] = "8M100T 8M100T"
    catalog_path.write_text(yaml.safe_dump(catalog))
    sheet = ws / "prepared/samples.tsv"
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
    ws = _ws(tmp_path, "targeted_umi", "FO002", depth=depth)
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["final_output"] = "FO004"
    cfg["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    cfg["analyses"].update(requested=["DA002"], resistance_catalog="inputs/metadata/test_resistance_catalog.tsv")
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    sheet = ws / "prepared/samples.tsv"
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
    """A WGS cohort that reaches all nine analyses (DA001-DA009), not just their skip paths.

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
    ws = _ws(tmp_path, "wgs", "FO002", depth=40, fixture_args=("--n-mono", "8", "--variant-scale", "4"))
    layer5 = ws / "config/parameters/layer5.yaml"
    text, n = re.subn(r"^(  da001_thin_bp: \{value: )\d+", r"\g<1>100", layer5.read_text(), flags=re.M)
    assert n == 1
    layer5.write_text(text)
    analyses = ["DA001", "DA002", "DA003", "DA004", "DA005", "DA006", "DA007", "DA008", "DA009"]
    # DA008 on the synthetic genes: no copy-number change was simulated, and one locus is absent.
    (ws / "config/registry/cnv_catalog_pfvarx_v1.tsv").write_text(
        "#catalog_id=test_cnv\n#version=test\n#status=TEST_FIXTURE\n"
        "locus_id\tgene_id\tgene_name\tevent\tmeaning\tevidence_class\tsource\n"
        "mdr1_amp\tPFX_MDR1\tmdr1\tAMPLIFICATION\ttest\tTEST\tfixture\n"
        "crt_del\tPFX_CRT\tcrt\tDELETION\ttest\tTEST\tfixture\n"
        "absent_amp\tPFX_NONE\tnone\tAMPLIFICATION\ttest\tTEST\tfixture\n"
        # DA009 reads these two loci: both genes exist in the synthetic GFF, so no deletion is simulated.
        "pfhrp2_deletion\tPFX_CRT\tpfhrp2\tDELETION\ttest\tTEST\tfixture\n"
        "pfhrp3_deletion\tPFX_K13\tpfhrp3\tDELETION\ttest\tTEST\tfixture\n")
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
    cnv_rows = [dict(zip(lines[0].split("\t"), row.split("\t")))
                for lines in [(res / "downstream/DA008/cnv_calls.tsv").read_text().splitlines()] for row in lines[1:]]
    by_locus: dict[str, set[str]] = {}
    for row in cnv_rows:
        by_locus.setdefault(row["locus_id"], set()).add(row["state"])
    assert by_locus["mdr1_amp"] == {"SINGLE_COPY"}, cnv_rows
    assert by_locus["crt_del"] == {"PRESENT"}, cnv_rows
    assert {r["reason"] for r in cnv_rows if r["locus_id"] == "absent_amp"} == {"GENE_NOT_IN_GFF"}
    # Short synthetic genes (2-3 windows at 40x) give ratios of about 0.7-1.4; the interval check keeps
    # them from being called AMPLIFIED, and nothing is called DELETED.
    assert all(0.5 <= float(r["depth_ratio"]) <= 1.6 for r in cnv_rows if r["locus_id"] != "absent_amp"), cnv_rows
    assert not {r["state"] for r in cnv_rows} & {"AMPLIFIED", "DELETED"}, cnv_rows
    # DA009: both genes are present in the simulated data, so every judged sample is NO_DELETION,
    # the WGS evidence model is recorded, and the pooled prevalence is zero with its interval.
    joint = [dict(zip(lines[0].split("\t"), row.split("\t")))
             for lines in [(res / "downstream/DA009/joint_calls.tsv").read_text().splitlines()] for row in lines[1:]]
    assert joint and {r["joint_state"] for r in joint} <= {"NO_DELETION", "INDETERMINATE"}, joint
    assert {r["evidence_class"] for r in joint} == {"DEPTH_RATIO"}, joint
    hrp_prev = [dict(zip(lines[0].split("\t"), row.split("\t")))
                for lines in [(res / "downstream/DA009/joint_prevalence.tsv").read_text().splitlines()]
                for row in lines[1:]]
    pooled = [r for r in hrp_prev if r["group"] == "ALL" and r["measure"] == "double_deletion"]
    assert len(pooled) == 1 and pooled[0]["carriers"] == "0" and pooled[0]["evidence_classes"] == "DEPTH_RATIO"
    assert (res / "downstream/DA009/joint_prevalence_by_evidence.tsv").is_file()
    assert (res / "downstream/DA009/cooccurrence.tsv").is_file()
    # Microhaplotypes: the default is off, so the artifact is empty with its reason recorded and
    # DA001 says which locus basis it used.
    mh = json.loads((res / "cohort/microhaplotypes.json").read_text())
    assert mh["loci_used"] == 0 and mh["reason"] == "MICROHAPLOTYPES_OFF" and mh["mode"] == "off"
    assert (res / "cohort/microhaplotype_loci.tsv").is_file()
    da001_method = json.loads((res / "downstream/DA001/status_manifest.json").read_text())["method"]
    assert da001_method["locus_basis"] == "BIALLELIC_SNV", da001_method
    assert "no microhaplotype locus" in da001_method["basis_reason"], da001_method

    # Both other modes, run against this run's own alignments and sites (no second pipeline run).
    run_dir = ws / "runs" / result["run_id"]
    plan = json.loads((run_dir / "frozen_plan.json").read_text())

    def microhaplotypes(mode: str) -> dict:
        work = tmp_path / f"mh_{mode}"
        work.mkdir()
        frozen = json.loads((run_dir / "frozen_config.json").read_text())
        frozen["config"]["analyses"]["microhap_loci"] = mode
        (work / "frozen_config.json").write_text(json.dumps(frozen))
        shutil.copy(run_dir / "frozen_plan.json", work / "frozen_plan.json")
        done = subprocess.run(["python", "-m", "pfvarx.cohort.microhap",
                               "--sites", str(res / "cohort/sites.tsv"),
                               "--samples", str(res / "cohort/cohort_samples.txt"),
                               "--bams", *[str(p) for p in sorted((res / "alignments").glob("*.bam"))],
                               "--callable", *[str(p) for p in sorted((res / "variants/samples").glob("*.callable.bed"))],
                               "--reference", str(ws / plan["reference"]["parasite_fasta"]),
                               "--root", str(ws),
                               "--frozen-config", str(work / "frozen_config.json"),
                               "--out-haplotypes", str(work / "mh.tsv.gz"),
                               "--out-loci", str(work / "loci.tsv"),
                               "--out-summary", str(work / "summary.json")], cwd=ws, capture_output=True, text=True)
        assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
        return {"summary": json.loads((work / "summary.json").read_text()), "dir": work}

    # Derived: loci come from this run's own core-genome windows, with no panel file at all.
    derived = microhaplotypes("derived")
    assert derived["summary"]["locus_source"] == "DERIVED_CORE_WINDOWS", derived["summary"]
    assert derived["summary"]["loci_used"] >= 1 and derived["summary"]["haplotypes"] >= 1, derived["summary"]
    assert "comparable only with cohorts" in derived["summary"]["comparability"]

    # Catalogue: a fixture panel over the same cohort's SNV-dense windows.
    sites = [row.split("\t") for row in (res / "cohort/sites.tsv").read_text().splitlines()[1:]]
    windows: dict[tuple[str, int], int] = {}
    for chrom, pos, _ref, _alt, kind, n_pass in ((r[1], int(r[2]), r[3], r[4], r[5], int(r[6])) for r in sites):
        if kind == "SNV" and n_pass > 0:
            windows[(chrom, pos // 300)] = windows.get((chrom, pos // 300), 0) + 1
    chosen = [(chrom, block) for (chrom, block), n in sorted(windows.items()) if n >= 2][:12]
    if not chosen:  # a locus with one variable position is still a valid (if uninformative) locus
        chosen = [(chrom, block) for (chrom, block), n in sorted(windows.items()) if n >= 1][:12]
    assert chosen, "the synthetic cohort should have PASS SNVs to build loci from"
    (ws / "config/registry/microhap_panel_pfvarx_v1.tsv").write_text(
        "#catalog_id=test_microhap\n#version=test\n#status=TEST_FIXTURE\n"
        "locus_id\tchrom\tstart\tend\tgene_id\tname\tsource\n"
        + "".join(f"mh{i}\t{chrom}\t{block * 300 + 1}\t{block * 300 + 300}\t\twindow\tfixture\n"
                  for i, (chrom, block) in enumerate(chosen)))
    catalogue = microhaplotypes("catalogue")
    summary, work = catalogue["summary"], catalogue["dir"]
    assert summary["panel_verified"] and summary["loci_used"] >= 1, summary
    assert summary["haplotypes"] >= 1 and summary["count_unit"] == "reads", summary
    called = [row.split("\t") for row in
              gzip.decompress((work / "mh.tsv.gz").read_bytes()).decode().splitlines()[1:]]
    assert any(r[5] == "CALLED" and r[2] for r in called), called[:5]
    # Every haplotype is as long as its locus has variable positions, and holds no unresolved base.
    positions = {r.split("\t")[0]: int(r.split("\t")[7] or 0)
                 for r in (work / "loci.tsv").read_text().splitlines()[1:]}
    assert all(len(r[2]) == positions[r[0]] and "N" not in r[2] for r in called if r[5] == "CALLED")
    # Cohort QC (G5) ran and gave each cohort sample a verdict.
    cohort_qc = json.loads((res / "cohort/cohort_qc.json").read_text())
    assert set(cohort_qc["samples"]) == set((res / "cohort/cohort_samples.txt").read_text().split())
    assert all(q["verdict"] in ("PASS", "WARNING", "FAIL") for q in cohort_qc["samples"].values())
    # Tool byproducts that no step reads are not kept.
    unused = [p for pattern in ("qc/samples/*/core.mosdepth.global.dist.txt", "qc/samples/*/core.regions.bed.gz*",
                                "downstream/DA006/tree.ckp.gz", "downstream/DA006/tree.bionj",
                                "downstream/DA006/tree.model.gz", "downstream/DA006/tree.splits.nex")
              for p in res.glob(pattern)]
    assert not unused, unused
    assert list(res.glob("qc/samples/*/core.mosdepth.region.dist.txt"))  # the file the QC gate reads stays
    report = (res / "report/clonality_report.txt").read_text()
    assert all(f"Sample PF{i:03d} is" in report for i in range(1, 12)), report


def test_human_dna_removal_drops_host_pairs_and_keeps_a_parasite_header(tmp_path: Path) -> None:
    """host_removal: true aligns to parasite+host, removes host pairs and counts them in QC."""
    import gzip
    import json
    import random

    ws = _ws(tmp_path, "wgs", "FO001")
    rng = random.Random(3)
    genome = "".join(rng.choices("ACGT", weights=[29.5, 20.5, 20.5, 29.5], k=40_000))
    (ws / "inputs/reference/host.fasta").write_text(">chr1\n" + genome + "\n")
    comp = str.maketrans("ACGT", "TGCA")
    sheet = ws / "prepared/samples.tsv"
    sheet.write_text("\n".join(sheet.read_text().splitlines()[:3]) + "\n")  # two samples
    for r1 in sorted((ws / "inputs/fastq").glob("PF00[12]_*_R1.fastq.gz")):
        r2 = Path(str(r1).replace("_R1.", "_R2."))
        old1, old2 = gzip.open(r1, "rt").read(), gzip.open(r2, "rt").read()
        n = len(old1.splitlines()[1])
        with gzip.open(r1, "wt") as h1, gzip.open(r2, "wt") as h2:
            h1.write(old1)
            h2.write(old2)
            for i in range(300):
                start = rng.randrange(0, len(genome) - 400)
                frag = genome[start:start + 300]
                h1.write(f"@host{i} 1:N:0:1\n{frag[:n]}\n+\n{'I' * n}\n")
                h2.write(f"@host{i} 2:N:0:1\n{frag[::-1].translate(comp)[:n]}\n+\n{'I' * n}\n")
    cfg = yaml.safe_load((ws / "config/config.yaml").read_text())
    cfg["reference"].update(host_removal=True, host_fasta="inputs/reference/host.fasta")
    (ws / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    # The tiny synthetic host index needs no extra memory (a human index does).
    params = ws / "config/parameters/layer1.yaml"
    params.write_text(params.read_text().replace("bwa_mem2_host_extra_mem_mb: {value: 12000",
                                                 "bwa_mem2_host_extra_mem_mb: {value: 0"))
    run = lambda *a: subprocess.run(["pfvarx", *a], cwd=ws, capture_output=True, text=True)  # noqa: E731
    assert run("validate").returncode == 0
    first = run("prepare")  # PfVarX never builds the host index itself
    assert first.returncode != 0 and "HOST_INDEX_REQUIRED" in first.stdout
    subprocess.run(["bwa-mem2", "index", "combined.fasta"], cwd=ws / "prepared/reference",
                   capture_output=True, check=True)
    assert run("prepare").returncode == 0
    assert json.loads(run("verify-prepared").stdout)["status"] == "PASS"
    r = run("pipeline")
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    res = ws / "runs" / json.loads(r.stdout)["run_id"] / "results"
    for sample in ("PF001", "PF002"):
        metrics = json.loads((res / f"qc/samples/{sample}/qc_decision.json").read_text())["metrics"]
        assert metrics["host_templates_removed"] == 300
        header = subprocess.run(["samtools", "view", "-H", str(res / f"alignments/{sample}.bam")],
                                capture_output=True, text=True, check=True).stdout
        assert [line.split("\t")[1] for line in header.splitlines() if line.startswith("@SQ")] == ["SN:ctg1", "SN:ctg2"]
