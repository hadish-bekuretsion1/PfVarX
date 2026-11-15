"""PfVarX-Call: per-sample SNP and small-indel caller (Layer 2 orchestration).

Implemented in this release (see docs/status.md for statuses):
  * L2.2 evidence extraction with pysam pileup (flag filter 0xF04, MAPQ/BQ admission,
    overlapping mates counted once via htslib overlap detection);
  * L2.4 per-sample beta-binomial error model fitted on known-reference positions
    (substitution: one class; indels: four homopolymer tract-length classes), with
    tail trimming and one re-estimation after removing screened candidates;
  * L2.5 pair-HMM re-scoring of indel candidates against local reference and
    alternative haplotypes, and read-backed linkage blocks within windows;
  * L2.6 pass-one clone mixture with parametric-bootstrap selection of K;
  * L2.7 WSAF profile-likelihood intervals; six decision statuses;
  * L2.9 VCFv4.3 output with declared FORMAT fields, callability BED, evidence TSV, QC JSON.
Not implemented (SPECIFIED only): the full covariate logit error model, combinatorial
window haplotype enumeration, the pass-two window-level mixture factor, representation-
bias correction of WSAF, and the optional classifier.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pysam

from pfvarx import __version__
from pfvarx.caller import mixture, pairhmm
from pfvarx.caller.stats import (
    betabinom_sf_log10,
    fisher_strand_neglog10,
    fit_betabinom,
    profile_interval,
)
from pfvarx.constants import CallStatus
from pfvarx.genomic_evidence import FLAG_FILTER, in_intervals, read_bed

BASES = "ACGT"
TRACT_CLASSES = ((1, 1), (2, 4), (5, 9), (10, 10_000))


def tract_class(length: int) -> int:
    for i, (lo, hi) in enumerate(TRACT_CLASSES):
        if lo <= length <= hi:
            return i
    return len(TRACT_CLASSES) - 1


@dataclass
class Candidate:
    chrom: str
    pos: int  # 0-based anchor position
    ref: str
    alt: str
    kind: str  # SNV | INS | DEL
    y: int
    n: int
    alt_f: int = 0
    alt_r: int = 0
    ref_f: int = 0
    ref_r: int = 0
    mean_alt_bq: float = 0.0
    tract: int = 1
    screen: float = 0.0
    status: str = ""
    filter_reasons: list[str] = field(default_factory=list)
    wsaf: float = 0.0
    wsaf_lo: float = 0.0
    wsaf_hi: float = 0.0
    lhid: str = "."
    lhwin: str = "."
    pairhmm_reads: int = 0


@dataclass
class Params:
    min_mapq: int
    min_bq: int
    window: int
    min_depth: int
    lam: float
    pass_min: float
    strand_max: float
    min_fragments: int
    eps0: float
    phi0: float
    trim: float
    gap_open: float
    gap_extend: float
    k_max: int
    alpha: float
    replicates: int
    restarts: int
    agreement: float
    em_tol: float
    em_iter: int
    min_mix_sites: int
    level: float
    seed: int
    min_weight: float
    candidate_min_alt_reads: int = 2
    pairhmm_max_reads: int = 200
    pairhmm_flank_bp: int = 60
    focal_min_length_bp: int = 50

    @classmethod
    def from_frozen(cls, eff: dict[str, Any]) -> Params:
        p = {k.split(".", 1)[1]: v["value"] for k, v in eff["parameters"].items() if k.startswith("layer2.")}
        return cls(min_mapq=int(p["min_mapping_quality"]), min_bq=int(p["min_base_quality"]),
                   window=int(p["window_size_bp"]), min_depth=int(p["min_usable_depth"]),
                   lam=float(p["screening_threshold_neglog10"]), pass_min=float(p["pass_min_neglog10"]),
                   strand_max=float(p["strand_bias_max_neglog10"]), min_fragments=int(p["pass_min_fragments"]),
                   eps0=float(p["error_mean_default"]), phi0=float(p["overdispersion_phi_default"]),
                   trim=float(p["known_ref_trim_fraction"]), gap_open=float(p["pairhmm_gap_open"]),
                   gap_extend=float(p["pairhmm_gap_extend"]), k_max=int(p["k_max"]),
                   alpha=float(p["k_test_alpha"]), replicates=int(p["k_bootstrap_replicates_analysis"]),
                   restarts=int(p["em_restarts"]), agreement=float(p["em_restart_agreement"]),
                   em_tol=float(p["em_tolerance"]), em_iter=int(p["em_max_iterations"]),
                   min_mix_sites=int(p["mixture_informative_min_sites"]),
                   level=float(p["wsaf_interval_level"]), seed=int(eff["config"]["random_seed"]),
                   min_weight=float(p["mixture_min_component_weight"]),
                   # .get keeps runs frozen before these became parameters reproducible.
                   candidate_min_alt_reads=int(p.get("candidate_min_alt_reads", 2)),
                   pairhmm_max_reads=int(p.get("pairhmm_max_reads", 200)),
                   pairhmm_flank_bp=int(p.get("pairhmm_flank_bp", 60)),
                   focal_min_length_bp=int(p.get("focal_loss_min_length_bp", 50)))


def homopolymer_after(seq: str, i: int) -> int:
    """Length of the homopolymer run starting at i+1 (the context an indel at anchor i extends)."""
    if i + 1 >= len(seq):
        return 1
    b = seq[i + 1]
    j = i + 1
    while j < len(seq) and seq[j] == b:
        j += 1
    return j - (i + 1)


def run_remaining(seq: str) -> np.ndarray:
    """rem[j] = length of the run of identical bases starting at j (vectorised)."""
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    if a.size == 0:
        return np.zeros(0, dtype=np.int32)
    change = np.concatenate([a[1:] != a[:-1], [True]])  # True at last index of each run
    ends = np.where(change)[0]
    run_end = np.repeat(ends, np.diff(np.concatenate([[-1], ends])))
    return (run_end - np.arange(a.size) + 1).astype(np.int32)


def tract_classes_array(seq: str) -> np.ndarray:
    rem = run_remaining(seq)
    after = np.concatenate([rem[1:], [1]])
    cls = np.full(after.shape, len(TRACT_CLASSES) - 1, dtype=np.int8)
    for i, (lo, hi) in reversed(list(enumerate(TRACT_CLASSES))):
        cls[(after >= lo) & (after <= hi)] = i
    return cls


class Evidence:
    """Per-position arrays over a contig plus sparse candidate allele tallies."""

    def __init__(self, length: int):
        self.n = np.zeros(length, dtype=np.int32)
        self.nonref = np.zeros(length, dtype=np.int32)
        self.indel = np.zeros(length, dtype=np.int32)
        self.in_core = np.zeros(length, dtype=bool)
        self.alleles: dict[int, dict[tuple[str, str], list[float]]] = defaultdict(dict)
        self.refstrand: dict[int, list[int]] = {}


def extract(bam: pysam.AlignmentFile, fasta: pysam.FastaFile, chrom: str,
            intervals: list[tuple[int, int]], p: Params) -> tuple[Evidence, str]:
    seq = fasta.fetch(chrom).upper()
    ev = Evidence(len(seq))
    for start, end in intervals:
        ev.in_core[start:end] = True
        for col in bam.pileup(chrom, start, end, truncate=True, stepper="nofilter",
                              flag_filter=FLAG_FILTER, min_mapping_quality=p.min_mapq,
                              min_base_quality=0, ignore_overlaps=True, ignore_orphans=False,
                              max_depth=1_000_000, compute_baq=False):
            pos = col.reference_pos
            ref = seq[pos]
            n = 0
            nonref = 0
            rf = rr = 0
            tallies = ev.alleles[pos] if pos in ev.alleles else {}
            local: dict[tuple[str, str], list[float]] = {}
            indel_n = 0
            for pr in col.pileups:
                if pr.is_refskip or pr.is_del:
                    continue
                aln = pr.alignment
                qpos = pr.query_position
                if qpos is None or aln.query_sequence is None:
                    continue
                bq = aln.query_qualities[qpos] if aln.query_qualities is not None else 0
                if bq < p.min_bq:
                    continue
                n += 1
                base = aln.query_sequence[qpos].upper()
                rev = aln.is_reverse
                if base == ref:
                    rr += int(rev)
                    rf += int(not rev)
                elif base in BASES:
                    nonref += 1
                    key = ("SNV", base)
                    t = local.setdefault(key, [0, 0, 0, 0.0])
                    t[0] += 1
                    t[1 if not rev else 2] += 1
                    t[3] += bq
                if pr.indel > 0:
                    ins = aln.query_sequence[qpos + 1: qpos + 1 + pr.indel].upper()
                    key = ("INS", ins)
                    indel_n += 1
                    t = local.setdefault(key, [0, 0, 0, 0.0])
                    t[0] += 1
                    t[1 if not rev else 2] += 1
                    t[3] += bq
                elif pr.indel < 0:
                    dele = seq[pos + 1: pos + 1 - pr.indel]
                    key = ("DEL", dele)
                    indel_n += 1
                    t = local.setdefault(key, [0, 0, 0, 0.0])
                    t[0] += 1
                    t[1 if not rev else 2] += 1
                    t[3] += bq
            ev.n[pos] = n
            ev.nonref[pos] = nonref
            ev.indel[pos] = indel_n
            keep = {k: v for k, v in local.items() if v[0] >= p.candidate_min_alt_reads}
            if keep:
                tallies.update(keep)
                ev.alleles[pos] = tallies
                ev.refstrand[pos] = [rf, rr]
    return ev, seq


def fit_error_model(evs: dict[str, Evidence], seqs: dict[str, str], p: Params,
                    exclude: set[tuple[str, int]] | None = None) -> dict[str, Any]:
    ys, ns = [], []
    iy: dict[int, list[np.ndarray]] = defaultdict(list)
    inn: dict[int, list[np.ndarray]] = defaultdict(list)
    for chrom, ev in evs.items():
        mask = ev.in_core & (ev.n > 0)
        if exclude:
            for excluded_chrom, pos in exclude:
                if excluded_chrom == chrom:
                    mask[pos] = False
        ys.append(ev.nonref[mask])
        ns.append(ev.n[mask])
        idx = np.where(mask)[0]
        cls = tract_classes_array(seqs[chrom])[idx]
        for c in range(len(TRACT_CLASSES)):
            sel = idx[cls == c]
            iy[c].append(ev.indel[sel])
            inn[c].append(ev.n[sel])
    y, n = np.concatenate(ys) if ys else np.array([]), np.concatenate(ns) if ns else np.array([])
    snv = fit_betabinom(y, n, p.trim, p.eps0, p.phi0)
    indel = {}
    for c in range(len(TRACT_CLASSES)):
        yy = np.concatenate(iy[c]) if iy[c] else np.array([])
        nn = np.concatenate(inn[c]) if inn[c] else np.array([])
        f = fit_betabinom(yy, nn, p.trim, p.eps0 / 5, p.phi0)
        f["mu"] = max(f["mu"], 1e-6)
        indel[str(c)] = f
    return {"snv_any": snv, "indel_by_tract_class": indel,
            "tract_classes": [list(t) for t in TRACT_CLASSES],
            "covariates_used": ["substitution_any", "homopolymer_tract_class"],
            "calibration": "PER_SAMPLE_KNOWN_REFERENCE_FIT"}


def screen(evs: dict[str, Evidence], seqs: dict[str, str], model: dict[str, Any], p: Params) -> list[Candidate]:
    out: list[Candidate] = []
    mu_snv = model["snv_any"]["mu"] / 3.0
    phi = model["snv_any"]["phi"]
    for chrom, ev in evs.items():
        seq = seqs[chrom]
        for pos, alleles in ev.alleles.items():
            n = int(ev.n[pos])
            rf, rr = ev.refstrand.get(pos, [0, 0])
            for (kind, allele), (y, af, ar, bqs) in alleles.items():
                tract = homopolymer_after(seq, pos)
                if kind == "SNV":
                    ref, alt = seq[pos], allele
                    mu, ph = mu_snv, phi
                elif kind == "INS":
                    ref, alt = seq[pos], seq[pos] + allele
                    m = model["indel_by_tract_class"][str(tract_class(tract))]
                    mu, ph = m["mu"], m["phi"]
                else:
                    ref, alt = seq[pos] + allele, seq[pos]
                    m = model["indel_by_tract_class"][str(tract_class(tract))]
                    mu, ph = m["mu"], m["phi"]
                c = Candidate(chrom, pos, ref, alt, kind, int(y), n, int(af), int(ar), rf, rr,
                              bqs / max(y, 1), tract)
                if n < p.min_depth:
                    c.status = CallStatus.INSUFFICIENT_DEPTH.value
                    out.append(c)
                    continue
                c.screen = betabinom_sf_log10(int(y), n, mu, ph)
                if c.screen >= p.lam:
                    out.append(c)
    out.sort(key=lambda c: (c.chrom, c.pos, c.kind, c.alt))
    return out


def rescore_indels(bam: pysam.AlignmentFile, seqs: dict[str, str], cands: list[Candidate],
                   p: Params, max_reads: int | None = None, flank: int | None = None) -> None:
    """Replace CIGAR counts of indel candidates with pair-HMM read support counts."""
    max_reads = p.pairhmm_max_reads if max_reads is None else max_reads
    flank = p.pairhmm_flank_bp if flank is None else flank
    for c in cands:
        if c.kind == "SNV" or c.status == CallStatus.INSUFFICIENT_DEPTH:
            continue
        seq = seqs[c.chrom]
        lo, hi = max(0, c.pos - flank), min(len(seq), c.pos + len(c.ref) + flank)
        ref_hap = seq[lo:hi]
        alt_hap = seq[lo:c.pos] + c.alt + seq[c.pos + len(c.ref):hi]
        alt = refc = 0
        af = ar = 0
        used = 0
        for read in bam.fetch(c.chrom, c.pos, c.pos + 1):
            if read.flag & FLAG_FILTER or read.mapping_quality < p.min_mapq or read.query_sequence is None or read.query_qualities is None:
                continue
            if read.reference_start > c.pos - 5 or read.reference_end is None or read.reference_end < c.pos + len(c.ref) + 5:
                continue  # read must span the event with anchoring bases
            used += 1
            if used > max_reads:
                break
            s = pairhmm.support(read.query_sequence, list(read.query_qualities), ref_hap, alt_hap,
                                p.gap_open, p.gap_extend)
            if s > 0:
                alt += 1
                af += int(not read.is_reverse)
                ar += int(read.is_reverse)
            elif s < 0:
                refc += 1
        c.pairhmm_reads = min(used, max_reads)
        if alt + refc > 0:
            c.y, c.n, c.alt_f, c.alt_r = alt, alt + refc, af, ar


def linkage_blocks(bam: pysam.AlignmentFile, cands: list[Candidate], p: Params) -> None:
    """Connected components of alternative alleles co-observed on >= 2 fragments per window."""
    by_window: dict[tuple[str, int], list[Candidate]] = defaultdict(list)
    for c in cands:
        if c.kind == "SNV" and c.status in (CallStatus.PASS, CallStatus.LOW_CONFIDENCE):
            by_window[(c.chrom, c.pos // p.window)].append(c)
    for (chrom, w), group in by_window.items():
        start, end = w * p.window, (w + 1) * p.window
        for c in group:
            c.lhwin = f"{chrom}_{start + 1}-{end}"
        if len(group) < 2:
            for c in group:
                c.lhid = f"{chrom}_{w}_{c.pos + 1}"
            continue
        positions = {c.pos: c for c in group}
        pair_counts: dict[tuple[int, int], int] = defaultdict(int)
        seen: dict[str, set[int]] = defaultdict(set)
        for read in bam.fetch(chrom, start, end):
            if read.flag & FLAG_FILTER or read.mapping_quality < p.min_mapq:
                continue
            if read.query_sequence is None or read.query_name is None:
                continue
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                cc = positions.get(rpos)
                if cc is not None and read.query_sequence[qpos].upper() == cc.alt:
                    seen[read.query_name].add(rpos)
        for alts in seen.values():
            alts_s = sorted(alts)
            for i in range(len(alts_s)):
                for j in range(i + 1, len(alts_s)):
                    pair_counts[(alts_s[i], alts_s[j])] += 1
        parent = {pos: pos for pos in positions}

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for (a, b), k in pair_counts.items():
            if k >= 2:
                parent[find(a)] = find(b)
        for pos, c in positions.items():
            c.lhid = f"{chrom}_{w}_{find(pos) + 1}"


def decide(cands: list[Candidate], model: dict[str, Any], p: Params, unreliable: dict[str, list[tuple[int, int]]],
           focal: dict[str, list[tuple[int, int]]]) -> None:
    for c in cands:
        if c.kind == "SNV":
            eps, phi = model["snv_any"]["mu"] / 3.0, model["snv_any"]["phi"]
        else:
            m = model["indel_by_tract_class"][str(tract_class(c.tract))]
            eps, phi = m["mu"], m["phi"]
        if c.n > 0:
            c.wsaf, c.wsaf_lo, c.wsaf_hi = profile_interval(c.y, c.n, min(eps, 0.25), phi, p.level)
        if c.status == CallStatus.INSUFFICIENT_DEPTH:
            if in_intervals(focal.get(c.chrom, []), c.pos):
                c.status = CallStatus.FOCAL_CALLABILITY_LOSS.value
            continue
        if in_intervals(unreliable.get(c.chrom, []), c.pos):
            c.status = CallStatus.CONTEXT_UNRELIABLE.value
            continue
        c.screen = betabinom_sf_log10(c.y, c.n, eps, phi) if c.n else 0.0
        sb = fisher_strand_neglog10(c.ref_f, c.ref_r, c.alt_f, c.alt_r)
        if sb > p.strand_max:
            c.status = CallStatus.FILTERED.value
            c.filter_reasons.append("STRAND_BIAS")
        elif c.screen >= p.pass_min and c.y >= p.min_fragments:
            c.status = CallStatus.PASS.value
        else:
            c.status = CallStatus.LOW_CONFIDENCE.value


def callable_intervals(evs: dict[str, Evidence], min_depth: int) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for chrom, ev in evs.items():
        ok = ev.in_core & (ev.n >= min_depth)
        out[chrom] = _runs(ok)
    return out


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    if not mask.any():
        return []
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    return list(zip(starts.tolist(), ends.tolist()))


def focal_losses(evs: dict[str, Evidence], min_depth: int, median: float, min_len: int = 50) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    if median < min_depth:
        return out
    for chrom, ev in evs.items():
        low = ev.in_core & (ev.n < min_depth)
        out[chrom] = [(s, e) for s, e in _runs(low) if e - s >= min_len]
    return out


def write_bed(path: Path, iv: dict[str, list[tuple[int, int]]]) -> None:
    with path.open("w") as h:
        for chrom in sorted(iv):
            for s, e in iv[chrom]:
                h.write(f"{chrom}\t{s}\t{e}\n")


def write_vcf(path: Path, sample: str, cands: list[Candidate], fasta: pysam.FastaFile,
              contig_md5: dict[str, str], model: dict[str, Any], calibration: str) -> None:
    lines = [
        "##fileformat=VCFv4.3",
        f"##source=PfVarX-Call {__version__}",
        f"##pfvarx_calibration={calibration}",
        "##pfvarx_record_set=variant-only; absent records are NOT reference calls",
        f"##pfvarx_error_model={json.dumps(model['snv_any'], sort_keys=True)}",
        "##pfvarx_feature_version=1",
    ]
    for name in fasta.references:
        lines.append(f"##contig=<ID={name},length={fasta.get_reference_length(name)},md5={contig_md5.get(name, '.')}>")
    for st in CallStatus:
        if st != CallStatus.PASS:
            lines.append(f'##FILTER=<ID={st.value},Description="PfVarX-Call decision status {st.value}">')
    lines += [
        '##FILTER=<ID=STRAND_BIAS,Description="Strand imbalance above the operating point">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="SNV, INS or DEL">',
        '##INFO=<ID=SCR,Number=1,Type=Float,Description="-log10 beta-binomial tail probability (Eq. 3.10)">',
        '##INFO=<ID=HP,Number=1,Type=Integer,Description="Homopolymer tract length following the anchor">',
        '##INFO=<ID=PHMM,Number=1,Type=Integer,Description="Reads re-scored by the pair-HMM (indels)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype of the dominant allele">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref, alt)">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Usable depth">',
        '##FORMAT=<ID=PST,Number=1,Type=String,Description="PfVarX-Call status">',
        '##FORMAT=<ID=WSAF,Number=1,Type=Float,Description="Within-sample alternative-allele DNA fraction (MLE)">',
        '##FORMAT=<ID=WSAFLO,Number=1,Type=Float,Description="Lower profile-likelihood bound">',
        '##FORMAT=<ID=WSAFHI,Number=1,Type=Float,Description="Upper profile-likelihood bound">',
        '##FORMAT=<ID=WSAFCOR,Number=1,Type=String,Description="Representation-bias correction state">',
        '##FORMAT=<ID=LHID,Number=1,Type=String,Description="Local haplotype (read-backed linkage block) id">',
        '##FORMAT=<ID=LHWIN,Number=1,Type=String,Description="Window (contig_start-end, 1-based) within which the haplotype claim holds">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample,
    ]
    for c in sorted(cands, key=lambda c: (fasta.references.index(c.chrom), c.pos, c.alt)):
        filt = "PASS" if c.status == CallStatus.PASS else ";".join([c.status] + c.filter_reasons)
        gt = "1" if c.wsaf >= 0.5 else "0"
        qual = "." if not math.isfinite(c.screen) else f"{min(c.screen, 9999):.2f}"
        info = f"TYPE={c.kind};SCR={qual if qual != '.' else 0};HP={c.tract}"
        if c.kind != "SNV":
            info += f";PHMM={c.pairhmm_reads}"
        fmt = "GT:AD:DP:PST:WSAF:WSAFLO:WSAFHI:WSAFCOR:LHID:LHWIN"
        val = (f"{gt}:{c.n - c.y},{c.y}:{c.n}:{c.status}:{c.wsaf:.4f}:{c.wsaf_lo:.4f}:{c.wsaf_hi:.4f}"
               f":UNCORRECTED:{c.lhid}:{c.lhwin}")
        lines.append(f"{c.chrom}\t{c.pos + 1}\t.\t{c.ref}\t{c.alt}\t{qual}\t{filt}\t{info}\t{fmt}\t{val}")
    path.write_text("\n".join(lines) + "\n")


_WORKER: dict[str, Any] = {}
_T0 = [0.0]


def _progress(message: str) -> None:
    """Stage progress on stderr (the job log); never written to the result files."""
    print(f"[pfvarx-call {time.monotonic() - _T0[0]:8.1f} s] {message}", file=sys.stderr, flush=True)


def _init_worker(bam_path: str, ref_path: str, p: Params) -> None:
    _WORKER.update(bam=pysam.AlignmentFile(bam_path), fasta=pysam.FastaFile(ref_path), p=p)


EXTRACT_CHUNK_BP = 250_000  # work unit for parallel evidence collection; results do not depend on it


def _extract_task(job: tuple[str, int, int]) -> tuple[str, int, int, dict[str, Any]]:
    """Evidence for one core interval piece: each position's evidence depends only on the
    reads covering it, so pieces can be collected independently and joined exactly."""
    chrom, start, end = job
    ev, _ = extract(_WORKER["bam"], _WORKER["fasta"], chrom, [(start, end)], _WORKER["p"])
    part = {"n": ev.n[start:end].copy(), "nonref": ev.nonref[start:end].copy(), "indel": ev.indel[start:end].copy(),
            "alleles": {pos: a for pos, a in ev.alleles.items() if a},
            "refstrand": dict(ev.refstrand)}
    return chrom, start, end, part


def _pieces(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out = []
    for start, end in intervals:
        for s in range(start, end, EXTRACT_CHUNK_BP):
            out.append((s, min(end, s + EXTRACT_CHUNK_BP)))
    return out


def _rescore_task(job: tuple[str, str, list[Candidate]]) -> list[Candidate]:
    chrom, seq, cands = job
    rescore_indels(_WORKER["bam"], {chrom: seq}, cands, _WORKER["p"])
    return cands


def rescore_parallel(pool: ProcessPoolExecutor, seqs: dict[str, str], cands: list[Candidate]) -> list[Candidate]:
    """Indel rescoring by contig in worker processes. Candidates are matched back by their
    index in *cands*, never by allele content: two deletions sharing an anchor base (TTA>T
    and TTATA>T) have the same ALT and would otherwise collide."""
    groups: dict[str, list[int]] = defaultdict(list)
    for i, c in enumerate(cands):
        groups[c.chrom].append(i)
    out = list(cands)
    jobs = [(contig, seqs[contig], [cands[i] for i in idx]) for contig, idx in groups.items()]
    for (contig, _, _), rescored in zip(jobs, pool.map(_rescore_task, jobs), strict=True):
        for i, c in zip(groups[contig], rescored, strict=True):
            out[i] = c
    return out


def run(args: argparse.Namespace) -> dict[str, Any]:
    # Fault injection for failure-isolation tests only (documented in docs/testing.md).
    if os.environ.get("PFVARX_FAULT_INJECT") == f"pfvarx_call:{args.sample}":
        print(f"PFVARX_FAULT_INJECT: deliberate failure for {args.sample} (test only)", file=sys.stderr)
        raise SystemExit(70)
    eff = json.loads(Path(args.frozen_config).read_text())
    p = Params.from_frozen(eff)
    core = read_bed(Path(args.core_mask))
    unreliable = read_bed(Path(args.context_unreliable)) if args.context_unreliable else {}
    fasta = pysam.FastaFile(args.reference)
    bam = pysam.AlignmentFile(args.bam)
    from pfvarx.inputs.validate import fasta_contigs
    md5 = {k: v["md5"] for k, v in fasta_contigs(Path(args.reference)).items()}
    evs: dict[str, Evidence] = {}
    seqs: dict[str, str] = {}
    threads = max(1, int(getattr(args, "threads", 1) or 1))
    _T0[0] = time.monotonic()
    _progress(f"sample {args.sample}: {threads} worker(s)")
    chroms = [c for c in fasta.references if c in core]
    # Chromosomes, indel candidates and bootstrap replicates are independent: with several
    # threads they run in worker processes, with results identical to the serial run.
    pool = ProcessPoolExecutor(threads, initializer=_init_worker, initargs=(args.bam, args.reference, p)) \
        if threads > 1 else None
    if pool is not None:
        for contig in chroms:
            seqs[contig] = fasta.fetch(contig).upper()
            evs[contig] = Evidence(len(seqs[contig]))
        jobs = [(contig, s, e) for contig in chroms for s, e in _pieces(core[contig])]
        parts: dict[tuple[str, int], dict[str, Any]] = {}
        for chrom, start, end, part in pool.map(_extract_task, jobs, chunksize=1):
            ev = evs[chrom]
            ev.n[start:end], ev.nonref[start:end], ev.indel[start:end] = part["n"], part["nonref"], part["indel"]
            ev.in_core[start:end] = True
            parts[(chrom, start)] = part
        for chrom, start, _ in jobs:  # same insertion order as a serial pass over the intervals
            part = parts[(chrom, start)]
            for pos, alleles in sorted(part["alleles"].items()):
                evs[chrom].alleles[pos] = alleles
                evs[chrom].refstrand[pos] = part["refstrand"][pos]
    else:
        for chrom in chroms:
            evs[chrom], seqs[chrom] = extract(bam, fasta, chrom, core[chrom], p)
    _progress(f"read evidence collected for {len(evs)} contig(s)")
    model = fit_error_model(evs, seqs, p)
    cands = screen(evs, seqs, model, p)
    # one re-estimation after removing screened candidate positions (field samples)
    excl = {(c.chrom, c.pos) for c in cands}
    model = fit_error_model(evs, seqs, p, exclude=excl)
    model["iterations"] = 2
    cands = screen(evs, seqs, model, p)
    _progress(f"error model fitted twice; {len(cands)} candidates screened")
    if pool is not None:
        cands = rescore_parallel(pool, seqs, cands)
    else:
        rescore_indels(bam, seqs, cands, p)
    _progress("indel candidates rescored with the pair-HMM")
    depth_all = np.concatenate([ev.n[ev.in_core] for ev in evs.values()]) if evs else np.array([0])
    median = float(np.median(depth_all)) if depth_all.size else 0.0
    focal = focal_losses(evs, p.min_depth, median, p.focal_min_length_bp)
    decide(cands, model, p, unreliable, focal)
    linkage_blocks(bam, cands, p)
    callable_iv = callable_intervals(evs, p.min_depth)
    core_bases = int(sum(ev.in_core.sum() for ev in evs.values()))
    callable_bases = sum(e - s for iv in callable_iv.values() for s, e in iv)
    # pass one mixture on PASS biallelic SNVs, thinned to >= 500 bp apart
    snv = [c for c in cands if c.kind == "SNV" and c.status == CallStatus.PASS]
    by_pos: dict[tuple[str, int], int] = defaultdict(int)
    for c in snv:
        by_pos[(c.chrom, c.pos)] += 1
    snv = [c for c in snv if by_pos[(c.chrom, c.pos)] == 1]
    thinned: list[Candidate] = []
    for c in snv:
        if not thinned or thinned[-1].chrom != c.chrom or c.pos - thinned[-1].pos >= 500:
            thinned.append(c)
    mix: dict[str, Any]
    _progress(f"decisions and linkage done; mixture model on {len(thinned)} sites")
    if len(thinned) >= p.min_mix_sites:
        y = np.array([c.y for c in thinned])
        n = np.array([c.n for c in thinned])
        mix = mixture.select_k(y, n, model["snv_any"]["mu"] / 3.0, model["snv_any"]["phi"],
                               p.k_max, p.alpha, p.replicates, p.seed, p.restarts, p.agreement,
                               p.em_tol, p.em_iter, p.min_weight, executor=pool, batch=threads)
        mix["status"] = "ESTIMATED"
        mix["informative_sites"] = len(thinned)
        mix["bootstrap_basis"] = (f"{p.replicates} replicates (UNCALIBRATED_PROVISIONAL; "
                                  "release design uses 1999 and a precomputed surface)")
    else:
        mix = {"status": "INSUFFICIENT_INFORMATIVE_SITES", "informative_sites": len(thinned),
               "required": p.min_mix_sites, "label": "technical mixture components; NOT biological COI"}
    _progress("mixture model selected")
    if pool is not None:
        pool.shutdown()
    calibration = "UNCALIBRATED"
    out_vcf = Path(args.out_vcf)
    write_vcf(out_vcf, args.sample, cands, fasta, md5, model, calibration)
    write_bed(Path(args.out_callable), callable_iv)
    write_bed(Path(args.out_focal), focal)
    cols = ["chrom", "pos", "ref", "alt", "kind", "y", "n", "alt_f", "alt_r", "ref_f", "ref_r",
            "mean_alt_bq", "tract", "screen", "status", "wsaf", "wsaf_lo", "wsaf_hi", "lhid", "pairhmm_reads"]
    with Path(args.out_evidence).open("w") as h:
        h.write("\t".join(cols) + "\n")
        for c in cands:
            h.write("\t".join(str(getattr(c, k) if k != "pos" else c.pos + 1) for k in cols) + "\n")
    counts: dict[str, int] = defaultdict(int)
    for c in cands:
        counts[f"{c.kind}:{c.status}"] += 1
    qc = {"schema": "pfvarx.caller_qc/1", "sample": args.sample, "caller_version": __version__,
          "calibration": calibration, "core_bases": core_bases, "callable_bases": callable_bases,
          "callable_fraction_core_caller": callable_bases / core_bases if core_bases else 0.0,
          "median_usable_depth_core": median, "records": len(cands), "status_counts": dict(counts),
          "pass_snv_per_callable_mb": (sum(1 for c in cands if c.kind == 'SNV' and c.status == 'PASS')
                                       / (callable_bases / 1e6) if callable_bases else None),
          "error_model": model, "mixture": mix,
          "mixture_converged": int(bool(mix.get("fit", {}).get("converged", False))) if mix["status"] == "ESTIMATED" else None,
          "parameters": p.__dict__}
    Path(args.out_qc).write_text(json.dumps(qc, indent=2, default=float))
    return qc


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="pfvarx-call", description=__doc__.splitlines()[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--core-mask", required=True)
    ap.add_argument("--frozen-config", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--context-unreliable")
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--out-evidence", required=True)
    ap.add_argument("--out-callable", required=True)
    ap.add_argument("--out-focal", required=True)
    ap.add_argument("--out-qc", required=True)
    ap.add_argument("--threads", type=int, default=1, help="worker processes; results do not depend on it")
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
