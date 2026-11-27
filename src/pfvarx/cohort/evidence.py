"""Layer 3 cohort evidence matrix and missingness semantics.

The merged VCF is variant-only: an absent record is NOT a reference call. For
downstream estimators that need a genotype at every cohort site, each sample's own
reads are re-queried at the union of cohort sites (same admission filters as the
caller). The resulting per-sample allele counts are *site-forced per-sample evidence*,
not joint genotyping; positions outside a sample's callable region are UNCALLABLE and
positions below the minimum usable depth are MISSING.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
from collections import deque
from collections.abc import Iterator
from concurrent.futures import Future, ProcessPoolExecutor
from pathlib import Path
from typing import Any, NamedTuple, TypedDict

import numpy as np
import pysam

from pfvarx.configuration.loader import frozen_param
from pfvarx.genomic_evidence import FLAG_FILTER, in_intervals, read_bed

STATES = ("ALT_DOMINANT", "MIXED", "REF_DOMINANT", "MISSING_LOW_DEPTH", "UNCALLABLE")


class CohortSite(TypedDict):
    chrom: str
    pos: int
    ref: str
    alt: str
    type: object
    n_pass: int
    statuses: dict[str, str | None]


def site_list(vcf: Path) -> list[CohortSite]:
    sites: list[CohortSite] = []
    with pysam.VariantFile(str(vcf)) as f:
        for rec in f:
            if rec.alts is None or rec.ref is None:
                continue
            statuses = {str(s): rec.samples[s].get("PST") for s in rec.samples}
            sites.append({"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0],
                          "type": rec.info.get("TYPE"), "n_pass": sum(1 for v in statuses.values() if v == "PASS"),
                          "statuses": statuses})
    return sites


def allele_counts(bam: pysam.AlignmentFile, site: CohortSite, min_mapq: int, min_bq: int,
                  max_depth: int) -> tuple[int, int, int]:
    """(ref, alt, depth) at one site; indels counted from CIGAR at the anchor; *max_depth* is
    layer2.pileup_max_depth."""
    chrom, pos = str(site["chrom"]), int(site["pos"]) - 1
    ref, alt = str(site["ref"]), str(site["alt"])
    kind = "SNV" if len(ref) == len(alt) == 1 else ("INS" if len(alt) > len(ref) else "DEL")
    r = a = d = 0
    for col in bam.pileup(chrom, pos, pos + 1, truncate=True, stepper="nofilter", flag_filter=FLAG_FILTER,
                          min_mapping_quality=min_mapq, min_base_quality=0, ignore_overlaps=True,
                          ignore_orphans=False, max_depth=max_depth, compute_baq=False):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            sequence = pr.alignment.query_sequence
            qualities = pr.alignment.query_qualities
            if sequence is None or qualities is None:
                continue
            q = qualities[pr.query_position]
            if q < min_bq:
                continue
            d += 1
            base = sequence[pr.query_position].upper()
            if kind == "SNV":
                if base == alt:
                    a += 1
                elif base == ref:
                    r += 1
            else:
                ind = pr.indel
                want = len(alt) - len(ref)
                if ind == want and (kind == "DEL" or sequence[
                        pr.query_position + 1: pr.query_position + 1 + ind].upper() == alt[1:]):
                    a += 1
                elif ind == 0:
                    r += 1
    return r, a, d


class SiteTable(NamedTuple):
    """The cohort sites and every sample's caller status, held compactly.

    ``status`` has one byte per site and sample: 0 for no call, otherwise 1 + the index of
    the value in ``vocabulary`` (PASS, LOW_CONFIDENCE, ...). One byte per cell instead of a
    Python dictionary entry (about 87 bytes) keeps large cohorts within memory.
    """

    sites: list[tuple[str, int, str, str, object, int]]  # chrom, pos (1-based), ref, alt, TYPE, n_pass
    status: np.ndarray  # (sites, samples) uint8
    vocabulary: list[str]
    columns: dict[str, int]  # sample ID -> column of ``status``


def site_table(vcf: Path) -> SiteTable:
    """Cohort sites with each sample's caller status (PST), as ``site_list`` but compact."""
    codes: dict[str, int] = {}
    rows: list[bytes] = []
    sites: list[tuple[str, int, str, str, object, int]] = []
    with pysam.VariantFile(str(vcf)) as f:
        names = [str(x) for x in f.header.samples]
        for rec in f:
            if rec.alts is None or rec.ref is None:
                continue
            row = bytearray(len(names))
            n_pass = 0
            for j, sample in enumerate(rec.samples.values()):
                value = sample.get("PST")
                if not value:
                    continue
                if value not in codes:
                    if len(codes) == 255:
                        raise ValueError("more than 255 distinct caller statuses in the cohort VCF")
                    codes[value] = len(codes) + 1
                row[j] = codes[value]
                n_pass += value == "PASS"
            rows.append(bytes(row))
            sites.append((rec.chrom, rec.pos, rec.ref, rec.alts[0], rec.info.get("TYPE"), n_pass))
    status = (np.frombuffer(b"".join(rows), dtype=np.uint8).reshape(len(rows), len(names)) if rows
              else np.zeros((0, len(names)), dtype=np.uint8))
    vocabulary = [value for value, _ in sorted(codes.items(), key=lambda item: item[1])]
    return SiteTable(sites, status, vocabulary, {name: j for j, name in enumerate(names)})


class Settings(NamedTuple):
    min_mapq: int
    min_bq: int
    max_depth: int
    min_depth: int
    dominant: float


_WORKER: dict[str, Any] = {}


def _init_worker(sites: list[tuple[str, int, str, str, object, int]], vocabulary: list[str],
                 settings: Settings) -> None:
    _WORKER.update(sites=sites, vocabulary=vocabulary, settings=settings)


def sample_rows(job: tuple[str, str, str, bytes]) -> tuple[str, dict[str, int]]:
    """One sample's matrix rows (in site order) and state counts.

    Each sample reads only its own alignment and callable regions, so samples are independent
    and can run in any process; the caller writes them in cohort order.
    """
    sid, bam_path, callable_path, status_column = job
    sites, vocabulary, cfg = _WORKER["sites"], _WORKER["vocabulary"], _WORKER["settings"]
    callable_iv = read_bed(Path(callable_path))
    counts = {k: 0 for k in STATES}
    lines = []
    with pysam.AlignmentFile(bam_path) as bam:
        for i, (chrom, pos, ref, alt, _type, _n_pass) in enumerate(sites):
            code = status_column[i]
            status = vocabulary[code - 1] if code else "."
            if not in_intervals(callable_iv.get(str(chrom), []), int(pos) - 1):
                state, r, al, d, w = "UNCALLABLE", 0, 0, 0, ""
            else:
                site = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
                r, al, d = allele_counts(bam, site, cfg.min_mapq, cfg.min_bq, cfg.max_depth)  # type: ignore[arg-type]
                tot = r + al
                if tot < cfg.min_depth:
                    state, w = "MISSING_LOW_DEPTH", ""
                else:
                    f = al / tot
                    w = f"{f:.4f}"
                    state = ("ALT_DOMINANT" if f >= cfg.dominant
                             else ("REF_DOMINANT" if f <= 1 - cfg.dominant else "MIXED"))
            counts[state] += 1
            lines.append(f"S{i}\t{sid}\t{state}\t{r}\t{al}\t{d}\t{w}\t{status}\n")
    return "".join(lines), counts


def _paths(listing: str | None, values: list[str] | None) -> list[str]:
    """File paths given on the command line or, for large cohorts, one per line in a file."""
    if listing:
        return [line for line in Path(listing).read_text().splitlines() if line.strip()]
    return list(values or [])


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--bams", nargs="+")
    ap.add_argument("--bam-list", help="file with one alignment per line, in cohort order (instead of --bams)")
    ap.add_argument("--callable", nargs="+")
    ap.add_argument("--callable-list", help="file with one callable BED per line, in cohort order")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--frozen-config", required=True)
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-missing", required=True)
    ap.add_argument("--out-sites", required=True)
    ap.add_argument("--threads", type=int, default=1, help="worker processes; the output does not depend on it")
    a = ap.parse_args(argv)
    eff = json.loads(Path(a.frozen_config).read_text())
    prm = {k: v["value"] for k, v in eff["parameters"].items()}
    settings = Settings(
        min_mapq=int(prm["layer2.min_mapping_quality"]), min_bq=int(prm["layer2.min_base_quality"]),
        max_depth=int(frozen_param({k.split(".", 1)[1]: v for k, v in prm.items() if k.startswith("layer2.")},
                                   "pileup_max_depth")),
        min_depth=int(prm["layer2.min_usable_depth"]), dominant=float(prm["layer3.dominant_allele_fraction"]))
    samples = Path(a.samples).read_text().split()
    bams, callables = _paths(a.bam_list, a.bams), _paths(a.callable_list, a.callable)
    if not (len(samples) == len(bams) == len(callables)):
        raise SystemExit(f"{len(samples)} samples, {len(bams)} alignments and {len(callables)} callable BEDs")
    table = site_table(Path(a.vcf))
    with Path(a.out_sites).open("w") as h:
        h.write("site_id\tchrom\tpos\tref\talt\ttype\tn_pass\n")
        for i, (chrom, pos, ref, alt, kind, n_pass) in enumerate(table.sites):
            h.write(f"S{i}\t{chrom}\t{pos}\t{ref}\t{alt}\t{kind}\t{n_pass}\n")
    empty = bytes(len(table.sites))
    jobs = [(sid, bam, cal, table.status[:, table.columns[sid]].tobytes() if sid in table.columns else empty)
            for sid, bam, cal in zip(samples, bams, callables, strict=True)]
    threads = max(1, a.threads)
    summary = []
    with gzip.open(a.out_matrix, "wt") as out:
        out.write("site_id\tsample_id\tstate\tref_count\talt_count\tdepth\twsaf\tcaller_status\n")
        for (sid, *_), (text, counts) in zip(jobs, _run(jobs, table, settings, threads), strict=True):
            out.write(text)
            n = len(table.sites) or 1
            summary.append({"sample_id": sid, **counts,
                            "missing_fraction": (counts["MISSING_LOW_DEPTH"] + counts["UNCALLABLE"]) / n})
    cols = ["sample_id", *STATES, "missing_fraction"]
    with Path(a.out_missing).open("w") as h:
        h.write("\t".join(cols) + "\n")
        for row in summary:
            h.write("\t".join(f"{row[c]:.4f}" if c == "missing_fraction" else str(row[c]) for c in cols) + "\n")


def _run(jobs: list[tuple[str, str, str, bytes]], table: SiteTable, settings: Settings,
         threads: int) -> Iterator[tuple[str, dict[str, int]]]:
    """Each sample's rows, in cohort order. With several threads, samples run in worker
    processes, and at most two per worker are finished ahead of the one being written, so
    memory stays bounded however large the cohort."""
    if threads == 1:
        _init_worker(table.sites, table.vocabulary, settings)
        yield from (sample_rows(job) for job in jobs)
        return
    with ProcessPoolExecutor(threads, initializer=_init_worker,
                             initargs=(table.sites, table.vocabulary, settings)) as pool:
        pending: deque[Future[tuple[str, dict[str, int]]]] = deque()
        queued = iter(jobs)
        for job in itertools.islice(queued, 2 * threads):
            pending.append(pool.submit(sample_rows, job))
        while pending:
            result = pending.popleft().result()
            for job in itertools.islice(queued, 1):
                pending.append(pool.submit(sample_rows, job))
            yield result


if __name__ == "__main__":
    main()
