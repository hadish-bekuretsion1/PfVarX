"""Layer 3 at cohort scale: the merge under open-file limits, compact caller statuses, and
cohort evidence and microhaplotypes in worker processes, each with the serial result."""

from __future__ import annotations

import gzip
import json
import random
import resource
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from pfvarx.cohort import evidence, microhap
from pfvarx.configuration.loader import load_parameters
from pfvarx.constants import Design

ROOT = Path(__file__).resolve().parents[2]
LENGTH = 6000
pytestmark = pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is not installed")


def _sample_vcfs(directory: Path, n: int, seed: int = 3) -> list[Path]:
    """Per-sample VCFs with several ALTs and indels sharing positions across samples."""
    rng = random.Random(seed)
    positions = sorted(rng.sample(range(50, LENGTH - 50), 60))
    choices = [("A", "C"), ("A", "G"), ("ATT", "A"), ("AT", "A"), ("A", "AT")]
    out = []
    for i in range(n):
        header = pysam.VariantHeader()
        header.add_line(f"##contig=<ID=ctg,length={LENGTH}>")
        header.add_line('##INFO=<ID=TYPE,Number=1,Type=String,Description="t">')
        for key, number, kind in (("GT", "1", "String"), ("AD", "R", "Integer"), ("DP", "1", "Integer"),
                                  ("PST", "1", "String")):
            header.add_line(f'##FORMAT=<ID={key},Number={number},Type={kind},Description="x">')
        header.add_sample(f"S{i:03d}")
        path = directory / f"S{i:03d}.vcf.gz"
        records = [(p, ref, alt) for p in positions if rng.random() < .4
                   for ref, alt in sorted(rng.sample(choices, rng.choice([1, 1, 2])), key=lambda x: (len(x[0]), x))]
        with pysam.VariantFile(str(path), "wz", header=header) as vcf:
            for p, ref, alt in sorted(records, key=lambda r: r[0]):
                rec = vcf.new_record(contig="ctg", start=p - 1, alleles=(ref, alt))
                rec.info["TYPE"] = "SNV" if len(ref) == len(alt) else "INDEL"
                s = rec.samples[0]
                s["GT"], s["AD"], s["DP"] = (1,), (2, 5), 7
                s["PST"] = rng.choice(["PASS", "PASS", "LOW_CONFIDENCE", "FILTERED"])
                vcf.write(rec)
        pysam.tabix_index(str(path), preset="vcf", force=True)
        out.append(path)
    return out


def _records(path: Path) -> list[str]:
    return subprocess.run(["bcftools", "view", "-H", str(path)], check=True, capture_output=True,
                          text=True).stdout.splitlines()


def _samples(path: Path) -> list[str]:
    return subprocess.run(["bcftools", "query", "-l", str(path)], check=True, capture_output=True,
                          text=True).stdout.split()


def _merge(tmp_path: Path, vcfs: list[Path], name: str, limit: tuple[int, int] | None,
           headroom: int = 8) -> subprocess.CompletedProcess[str]:
    listing = tmp_path / f"{name}.txt"
    listing.write_text("".join(f"{v}\n" for v in vcfs))

    def restrict() -> None:
        if limit is not None:
            resource.setrlimit(resource.RLIMIT_NOFILE, limit)

    return subprocess.run([sys.executable, "-m", "pfvarx.cohort.merge", "--list", str(listing),
                           "--out", str(tmp_path / f"{name}.vcf.gz"), "--headroom", str(headroom)],
                          capture_output=True, text=True, preexec_fn=restrict)


@pytest.fixture(scope="module")
def vcfs(tmp_path_factory: pytest.TempPathFactory) -> list[Path]:
    return _sample_vcfs(tmp_path_factory.mktemp("vcfs"), 30)


def test_merge_raises_a_low_soft_limit_and_gives_the_direct_merge(tmp_path: Path, vcfs: list[Path]) -> None:
    listing = tmp_path / "all.txt"
    listing.write_text("".join(f"{v}\n" for v in vcfs))
    subprocess.run(["bcftools", "merge", "-m", "none", "-i", "-", "-Oz", "-o", str(tmp_path / "plain.vcf.gz"),
                    "-l", str(listing)], check=True)
    hard = resource.getrlimit(resource.RLIMIT_NOFILE)[1]
    result = _merge(tmp_path, vcfs, "raised", (20, hard))  # 30 inputs under a soft limit of 20
    assert result.returncode == 0, result.stderr
    assert "strategy direct" in result.stderr
    assert _records(tmp_path / "raised.vcf.gz") == _records(tmp_path / "plain.vcf.gz")  # same file content, same order
    assert (tmp_path / "raised.vcf.gz.tbi").is_file()


def test_merge_beyond_the_hard_limit_batches_with_the_same_records(tmp_path: Path, vcfs: list[Path]) -> None:
    direct = _merge(tmp_path, vcfs, "direct", None)
    batched = _merge(tmp_path, vcfs, "batched", (24, 24))  # room for 16 inputs at a time
    assert direct.returncode == 0 and batched.returncode == 0, batched.stderr
    assert "strategy batched:16" in batched.stderr
    assert _samples(tmp_path / "batched.vcf.gz") == _samples(tmp_path / "direct.vcf.gz") == [
        f"S{i:03d}" for i in range(30)]
    assert sorted(_records(tmp_path / "batched.vcf.gz")) == sorted(_records(tmp_path / "direct.vcf.gz"))


def test_merge_stops_with_a_clear_message_when_the_limit_is_too_small(tmp_path: Path, vcfs: list[Path]) -> None:
    result = _merge(tmp_path, vcfs, "tiny", (12, 12))  # 4 inputs at a time: 8 batches > 4
    assert result.returncode != 0 and "ulimit -n" in result.stderr


def test_merge_copies_a_single_sample(tmp_path: Path, vcfs: list[Path]) -> None:
    assert _merge(tmp_path, vcfs[:1], "one", None).returncode == 0
    assert (tmp_path / "one.vcf.gz").read_bytes() == vcfs[0].read_bytes()


def test_site_table_holds_what_site_list_holds(tmp_path: Path, vcfs: list[Path]) -> None:
    assert _merge(tmp_path, vcfs, "cohort", None).returncode == 0
    cohort = tmp_path / "cohort.vcf.gz"
    old, new = evidence.site_list(cohort), evidence.site_table(cohort)
    assert [(s["chrom"], s["pos"], s["ref"], s["alt"], s["type"], s["n_pass"]) for s in old] == new.sites
    for i, site in enumerate(old):
        for name, column in new.columns.items():
            code = int(new.status[i, column])
            assert (site["statuses"][name] or ".") == (new.vocabulary[code - 1] if code else ".")
    assert new.status.dtype.itemsize == 1


def _bam(path: Path, seed: int) -> Path:
    rng = random.Random(seed)
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"},
                                              "SQ": [{"SN": "ctg", "LN": LENGTH}]})
    raw = path.with_suffix(".raw.bam")
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        for i in range(900):
            read = pysam.AlignedSegment(header)
            start = rng.randrange(0, LENGTH - 100)
            read.query_name, read.reference_id, read.reference_start = f"r{i}", 0, start
            read.query_sequence = "".join(rng.choice("ACGT") for _ in range(100))
            read.cigarstring, read.mapping_quality = "100M", 60
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.is_reverse = rng.random() < .5
            out.write(read)
    pysam.sort("-o", str(path), str(raw))
    pysam.index(str(path))
    return path


def test_cohort_evidence_is_the_same_with_worker_processes(tmp_path: Path, vcfs: list[Path]) -> None:
    group = vcfs[:6]
    assert _merge(tmp_path, group, "cohort", None).returncode == 0
    samples = [v.name.split(".")[0] for v in group]
    (tmp_path / "samples.txt").write_text("\n".join(samples) + "\n")
    bams = [_bam(tmp_path / f"{s}.bam", i) for i, s in enumerate(samples)]
    beds = []
    for i, s in enumerate(samples):  # callable regions differ between samples
        bed = tmp_path / f"{s}.callable.bed"
        bed.write_text(f"ctg\t{100 * i}\t{LENGTH - 300 * i}\n")
        beds.append(bed)
    (tmp_path / "bams.txt").write_text("".join(f"{b}\n" for b in bams))
    (tmp_path / "beds.txt").write_text("".join(f"{b}\n" for b in beds))
    frozen = tmp_path / "frozen.json"
    frozen.write_text(json.dumps({"parameters": load_parameters(ROOT)}))
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">ctg\n" + "A" * LENGTH + "\n")
    outputs = {}
    for label, extra in (("args", ["--bams", *map(str, bams), "--callable", *map(str, beds)]),
                         ("lists1", ["--bam-list", str(tmp_path / "bams.txt"), "--callable-list", str(tmp_path / "beds.txt")]),
                         ("lists3", ["--bam-list", str(tmp_path / "bams.txt"), "--callable-list",
                                     str(tmp_path / "beds.txt"), "--threads", "3"])):
        out = tmp_path / label
        out.mkdir()
        evidence.main(["--vcf", str(tmp_path / "cohort.vcf.gz"), "--samples", str(tmp_path / "samples.txt"),
                       "--reference", str(fasta), "--frozen-config", str(frozen), *extra,
                       "--out-matrix", str(out / "m.tsv.gz"), "--out-missing", str(out / "miss.tsv"),
                       "--out-sites", str(out / "sites.tsv")])
        outputs[label] = (gzip.decompress((out / "m.tsv.gz").read_bytes()), (out / "miss.tsv").read_bytes(),
                          (out / "sites.tsv").read_bytes())
    assert outputs["args"] == outputs["lists1"] == outputs["lists3"]
    states = {line.split(b"\t")[2] for line in outputs["args"][0].splitlines()[1:]}
    assert b"UNCALLABLE" in states and len(states) >= 3, states  # callable and uncallable cells both tested


def test_microhaplotypes_are_the_same_with_worker_processes(tmp_path: Path) -> None:
    from tests.unit.test_microhap import PARAMS, POSITIONS, _spanning
    from tests.unit.test_microhap import _bam as hap_bam

    sample_bams = {}
    for i, (name, mix) in enumerate((("RICH", 12), ("MID", 7), ("THIN", 3), ("EMPTY", 0))):
        (tmp_path / name).mkdir()
        reads = [_spanning(f"{name}a{k}", ("A", "C", "G")) for k in range(mix)]
        reads += [_spanning(f"{name}b{k}", ("A", "T", "G")) for k in range(mix // 2 + i)]
        sample_bams[name] = hap_bam(tmp_path / name, reads or [_spanning("x", ("A", "A", "A"), start=2000)])
    bed = tmp_path / "mid.bed"
    bed.write_text("ctg\t0\t1020\n")  # MID: one position outside its callable region
    sites = [{"chrom": "ctg", "pos": str(p + 1), "ref": "A", "alt": "C", "n_pass": "2"} for p in POSITIONS]
    panel = [{"locus_id": f"m{k}", "chrom": "ctg", "start": str(950 + k), "end": "1150", "gene_id": "",
              "name": "", "source": "fixture"} for k in range(3)]
    panel.append({"locus_id": "out", "chrom": "ctg", "start": "3000", "end": "3200", "gene_id": "",
                  "name": "", "source": "fixture"})
    args = dict(sites=sites, panel=panel, samples=list(sample_bams), bams=sample_bams, callable_beds={"MID": bed},
                lengths={"ctg": LENGTH}, territory={"ctg": [(900, 1200)]}, design=Design.TARGETED_UMI, gff=None,
                params=PARAMS, min_mapq=20, min_bq=13)
    serial = microhap.build(**args, threads=1)  # type: ignore[arg-type]
    parallel = microhap.build(**args, threads=3)  # type: ignore[arg-type]
    assert serial == parallel
    reasons = {r["reason"] for r in serial[0]}
    assert {"", "LOW_FRAGMENT_DEPTH", "POSITION_NOT_CALLABLE"} <= reasons
