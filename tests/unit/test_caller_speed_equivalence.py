"""The faster caller paths give exactly the results of the reference paths.

* ``extract`` (pysam column accessors) against ``extract_reference`` (one read object per
  read and column), on reads with every case the evidence rules treat specially;
* indel rescoring and linkage blocks in worker processes against the serial functions.
"""

import dataclasses
import random
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pysam
import pytest

from pfvarx.caller import call
from pfvarx.caller.call import Candidate, Params, extract, extract_reference
from pfvarx.configuration.loader import load_parameters

ROOT = Path(__file__).resolve().parents[2]
CONTIGS = {"chrA": 3000, "chrB": 1800}


@pytest.fixture(scope="module")
def params() -> Params:
    return Params.from_frozen({"parameters": load_parameters(ROOT), "config": {"random_seed": 20261106}})


def _reference(rng: random.Random) -> dict[str, str]:
    """AT-rich sequence with long homopolymer tracts, as in P. falciparum."""
    out = {}
    for name, length in CONTIGS.items():
        parts: list[str] = []
        while sum(map(len, parts)) < length:
            parts.append(rng.choice("AT") * rng.randint(1, 14) if rng.random() < .3 else rng.choice("AAATTTCG"))
        out[name] = "".join(parts)[:length]
    return out


def _clone(ref: dict[str, str]) -> dict[str, str]:
    """A second clone: substitutions every 37 bp, so reads from it carry linked alleles."""
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    return {name: "".join(swap[b] if i % 37 == 5 else b for i, b in enumerate(seq)) for name, seq in ref.items()}


def _read(rng: random.Random, header: pysam.AlignmentHeader, ref: dict[str, str], i: int,
          clone: dict[str, str]) -> pysam.AlignedSegment:
    chrom = rng.choice(list(ref))
    seq = clone[chrom] if rng.random() < .35 else ref[chrom]
    start = rng.randint(0, len(seq) - 200)
    ops: list[tuple[int, int]] = []
    query: list[str] = []
    if rng.random() < .15:  # leading soft clip
        k = rng.randint(1, 8)
        ops.append((4, k))
        query += rng.choices("ACGTN", k=k)
    rpos = start
    while len(query) < 100:
        r = rng.random()
        if r < .08 and ops and ops[-1][0] == 0:  # insertion, up to 12 bases (two-digit lengths)
            k = rng.randint(1, 12)
            ops.append((1, k))
            query += rng.choices("AT" if rng.random() < .7 else "ACGTN", k=k)
        elif r < .16 and ops and ops[-1][0] == 0:  # deletion, up to 12 bases
            k = rng.randint(1, 12)
            ops.append((2, k))
            rpos += k
        elif r < .18 and ops and ops[-1][0] == 0:  # reference skip
            k = rng.randint(5, 30)
            ops.append((3, k))
            rpos += k
        else:
            k = rng.randint(5, 30)
            for _ in range(k):
                base = seq[rpos]
                if rng.random() < .04:
                    base = rng.choice("ACGTN")
                query.append(base)
                rpos += 1
            ops.append((0, k))
    if rng.random() < .15:  # trailing soft clip
        k = rng.randint(1, 8)
        ops.append((4, k))
        query += rng.choices("ACGT", k=k)
    read = pysam.AlignedSegment(header)
    read.query_name = f"r{i // 2 if rng.random() < .3 else i}"  # some names shared: mates and overlaps
    read.reference_id = header.get_tid(chrom)
    read.reference_start = start
    read.mapping_quality = rng.choice([0, 5, 20, 40, 60, 60, 60])
    read.cigartuples = ops
    flag = 0x10 if rng.random() < .5 else 0
    flag |= rng.choice([0, 0, 0, 0, 0x400, 0x100, 0x200, 0x800])
    if rng.random() < .4:
        flag |= 0x1 | 0x2 | (0x40 if rng.random() < .5 else 0x80)
        read.next_reference_id = read.reference_id
        read.next_reference_start = max(0, start + rng.randint(-60, 60))
    read.flag = flag
    kind = rng.random()
    if kind < .03:  # no sequence stored
        read.query_sequence = None
        read.cigartuples = [(0, rpos - start)]
    else:
        read.query_sequence = "".join(query)
        if kind > .06:  # otherwise: no base qualities stored
            read.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(33 + rng.choice([2, 8, 12, 13, 14, 20, 30, 37, 40, 41])) for _ in query))
    return read


@pytest.fixture(scope="module")
def data(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, dict[str, str]]:
    rng = random.Random(20261127)
    directory = tmp_path_factory.mktemp("equivalence")
    ref = _reference(rng)
    fasta = directory / "ref.fa"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in ref.items()))
    pysam.faidx(str(fasta))
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": k, "LN": v} for k, v in CONTIGS.items()]})
    clone = _clone(ref)
    reads = sorted((_read(rng, header, ref, i, clone) for i in range(2500)),
                   key=lambda r: (r.reference_id, r.reference_start))
    unsorted = directory / "reads.unsorted.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for read in reads:
            out.write(read)
    bam = directory / "reads.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    return bam, fasta, ref


def _assert_same(a: call.Evidence, b: call.Evidence) -> None:
    for name in ("n", "nonref", "indel", "in_core"):
        assert np.array_equal(getattr(a, name), getattr(b, name)), name
    assert list(a.alleles.items()) == list(b.alleles.items())  # insertion order too
    assert a.refstrand == b.refstrand


@pytest.mark.parametrize("min_bq", [1, 10, 13, 20])
def test_column_evidence_matches_reference(data: tuple[Path, Path, dict[str, str]], params: Params,
                                           min_bq: int) -> None:
    bam_path, fasta_path, _ = data
    p = dataclasses.replace(params, min_bq=min_bq, candidate_min_alt_reads=1)
    with pysam.AlignmentFile(str(bam_path)) as bam, pysam.FastaFile(str(fasta_path)) as fasta:
        for chrom, length in CONTIGS.items():
            for intervals in ([(0, length)], [(0, 700), (701, 1500), (1700, length)]):
                fast, seq = extract(bam, fasta, chrom, intervals, p)
                ref, ref_seq = extract_reference(bam, fasta, chrom, intervals, p)
                assert seq == ref_seq
                _assert_same(fast, ref)
                assert fast.alleles, "the test reads reach no candidate alleles"
                kinds = {kind for alleles in fast.alleles.values() for kind, _ in alleles}
                assert kinds == {"SNV", "INS", "DEL"}


def test_zero_quality_floor_uses_reference(data: tuple[Path, Path, dict[str, str]], params: Params) -> None:
    bam_path, fasta_path, _ = data
    p = dataclasses.replace(params, min_bq=0, candidate_min_alt_reads=1)
    with pysam.AlignmentFile(str(bam_path)) as bam, pysam.FastaFile(str(fasta_path)) as fasta:
        _assert_same(extract(bam, fasta, "chrB", [(0, CONTIGS["chrB"])], p)[0],
                     extract_reference(bam, fasta, "chrB", [(0, CONTIGS["chrB"])], p)[0])


def _candidates(bam_path: Path, fasta_path: Path, p: Params) -> list[Candidate]:
    out = []
    with pysam.AlignmentFile(str(bam_path)) as bam, pysam.FastaFile(str(fasta_path)) as fasta:
        for chrom, length in CONTIGS.items():
            ev, seq = extract(bam, fasta, chrom, [(0, length)], p)
            for pos, alleles in ev.alleles.items():
                for (kind, allele), (y, af, ar, _) in alleles.items():
                    ref, alt = ((seq[pos], allele) if kind == "SNV" else (seq[pos], seq[pos] + allele)
                                if kind == "INS" else (seq[pos] + allele, seq[pos]))
                    status = ("INSUFFICIENT_DEPTH" if ev.n[pos] < 3 else "PASS" if pos % 3 else "LOW_CONFIDENCE")
                    out.append(Candidate(chrom, pos, ref, alt, kind, int(y), int(ev.n[pos]), int(af), int(ar),
                                         status=status))
    return out


def test_parallel_rescoring_and_linkage_match_serial(data: tuple[Path, Path, dict[str, str]],
                                                     params: Params) -> None:
    bam_path, fasta_path, ref = data
    p = dataclasses.replace(params, candidate_min_alt_reads=1, window=150)
    serial = _candidates(bam_path, fasta_path, p)
    parallel = [dataclasses.replace(c) for c in serial]
    assert sum(c.kind != "SNV" for c in serial) > 20
    with pysam.AlignmentFile(str(bam_path)) as bam:
        call.rescore_indels(bam, ref, serial, p)
        call.linkage_blocks(bam, serial, p)
    with ProcessPoolExecutor(2, initializer=call._init_worker, initargs=(str(bam_path), str(fasta_path), p)) as pool, \
            pysam.AlignmentFile(str(bam_path)) as bam:
        parallel = call.rescore_parallel(pool, ref, parallel)
        call.linkage_blocks(bam, parallel, p, pool)
    assert [dataclasses.asdict(c) for c in parallel] == [dataclasses.asdict(c) for c in serial]
    assert any(c.pairhmm_reads for c in serial)
    linked = [c.lhid for c in serial if c.lhid != "."]
    assert len(set(linked)) < len(linked), "no two alleles were linked"


def test_parallel_error_model_matches_serial(data: tuple[Path, Path, dict[str, str]], params: Params) -> None:
    bam_path, fasta_path, ref = data
    p = dataclasses.replace(params, error_fit_min_sites=10)
    with pysam.AlignmentFile(str(bam_path)) as bam, pysam.FastaFile(str(fasta_path)) as fasta:
        evs = {chrom: extract(bam, fasta, chrom, [(0, length)], p)[0] for chrom, length in CONTIGS.items()}
    serial = call.fit_error_model(evs, ref, p)
    exclude = {(c.chrom, c.pos) for c in call.screen(evs, ref, serial, p)}
    with ProcessPoolExecutor(2, initializer=call._init_worker, initargs=(str(bam_path), str(fasta_path), p)) as pool:
        assert call.fit_error_model(evs, ref, p, pool=pool) == serial
        assert call.fit_error_model(evs, ref, p, exclude=exclude, pool=pool) == call.fit_error_model(
            evs, ref, p, exclude=exclude)
    assert serial["snv_any"]["fitted"]
