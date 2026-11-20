"""Deterministic synthetic fixtures for software testing (NOT biological validation).

Generates an AT-rich reference with a GFF annotation, haploid clones carrying known
SNPs and indels, monoclonal and polyclonal samples with known clone proportions,
paired-end FASTQ for WGS, targeted (amplicon) and targeted-UMI designs, metadata with
populations, coordinates and dates, a core-genome mask, a test resistance catalogue
whose markers lie in the synthetic genes, and truth files for the benchmarking harness.
"""

from __future__ import annotations

import gzip
import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from pfvarx.configuration.library import CATALOG_SCHEMA
from pfvarx.constants import Design

CODON = {"Cys": "TGT", "Tyr": "TAT", "Lys": "AAA", "Thr": "ACA", "Asn": "AAT", "Ile": "ATA",
         "Ser": "TCA", "Ala": "GCA", "Arg": "AGA", "Gly": "GGA", "Leu": "TTA", "Phe": "TTT",
         "Met": "ATG", "Stop": "TAA"}
COMP = str.maketrans("ACGTN", "TGCAN")

# Synthetic genes (forward strand) and test markers (gene, codon, ref aa, alt aa).
GENES = [("PFX_K13", "ctg1", 6001, 300), ("PFX_CRT", "ctg1", 14001, 150), ("PFX_MDR1", "ctg2", 8001, 200)]
MARKERS = [("PFX_K13", 120, "Cys", "Tyr"), ("PFX_K13", 200, "Arg", "Lys"), ("PFX_CRT", 76, "Lys", "Thr"),
           ("PFX_MDR1", 86, "Asn", "Tyr")]


@dataclass
class Variant:
    chrom: str
    pos: int  # 1-based VCF position
    ref: str
    alt: str


def make_reference(rng: random.Random, lengths: dict[str, int]) -> dict[str, list[str]]:
    seqs: dict[str, list[str]] = {}
    for name, length in lengths.items():
        s: list[str] = []
        while len(s) < length:
            if rng.random() < 0.03:
                s.extend(rng.choice("AT") * rng.randint(6, 14))  # homopolymer tracts
            else:
                s.append(rng.choices("ATCG", weights=[0.4, 0.4, 0.1, 0.1])[0])
        seqs[name] = s[:length]
    # Write genes: ATG + codons + stop, avoiding internal stops; marker codons fixed.
    for gene, chrom, start, ncodons in GENES:
        codons = ["ATG"]
        for i in range(2, ncodons):
            mk = [m for m in MARKERS if m[0] == gene and m[1] == i]
            if mk:
                codons.append(CODON[mk[0][2]])
                continue
            while True:
                c = "".join(rng.choices("ATCG", weights=[0.35, 0.35, 0.15, 0.15], k=3))
                if c not in ("TAA", "TAG", "TGA"):
                    break
            codons.append(c)
        codons.append("TAA")
        cds = "".join(codons)
        seqs[chrom][start - 1:start - 1 + len(cds)] = list(cds)
    return seqs


def marker_positions() -> list[dict[str, Any]]:
    out = []
    for gene, codon, ref, alt in MARKERS:
        g = next(x for x in GENES if x[0] == gene)
        first = g[2] + (codon - 1) * 3  # 1-based first base of codon
        out.append({"gene": gene, "chrom": g[1], "codon": codon, "ref_aa": ref, "alt_aa": alt,
                    "codon_start": first, "ref_codon": CODON[ref], "alt_codon": CODON[alt]})
    return out


def make_clones(rng: random.Random, ref: dict[str, list[str]], n_clones: int, core: dict[str, tuple[int, int]],
                variant_scale: int = 1) -> tuple[list[dict[str, list[str]]], list[list[Variant]]]:
    """Clones from two populations: population-specific and shared SNPs, plus indels.

    ``variant_scale`` multiplies the number of simulated sites per contig (1 keeps the
    original fixture exactly); larger values give enough informative markers for the
    relatedness, structure and phylogeny analyses."""
    pop_sites: dict[int, list[Variant]] = {0: [], 1: []}
    shared: list[Variant] = []
    taken: set[tuple[str, int]] = set()
    marker_cod = {(m["chrom"], m["codon_start"] + k) for m in marker_positions() for k in range(3)}

    def pick_site(chrom: str) -> int:
        lo, hi = core[chrom]
        while True:
            p = rng.randint(lo + 200, hi - 200)
            if all((chrom, p + d) not in taken for d in range(-12, 13)) and (chrom, p) not in marker_cod:
                return p

    def snv(chrom: str) -> Variant:
        p = pick_site(chrom)
        taken.add((chrom, p))
        r = ref[chrom][p - 1]
        a = rng.choice([b for b in "ACGT" if b != r])
        return Variant(chrom, p, r, a)

    def indel(chrom: str) -> Variant:
        p = pick_site(chrom)
        taken.add((chrom, p))
        anchor = ref[chrom][p - 1]
        if rng.random() < 0.5:
            k = rng.randint(1, 4)
            return Variant(chrom, p, anchor + "".join(ref[chrom][p:p + k]), anchor)
        ins = "".join(rng.choices("AT", k=rng.randint(1, 4)))
        return Variant(chrom, p, anchor, anchor + ins)

    for chrom in ref:
        for _ in range(25 * variant_scale):
            shared.append(snv(chrom))
        for pop in (0, 1):
            for _ in range(35 * variant_scale):
                pop_sites[pop].append(snv(chrom))
            for _ in range(4 * variant_scale):
                pop_sites[pop].append(indel(chrom))
    clones_seq, clones_var = [], []
    markers = marker_positions()
    for i in range(n_clones):
        pop = i % 2
        vs: list[Variant] = []
        vs += [v for v in shared if rng.random() < 0.5]
        vs += [v for v in pop_sites[pop] if rng.random() < 0.7]
        vs += [v for v in pop_sites[1 - pop] if rng.random() < 0.1]
        for m in markers:  # resistance markers carried by some clones
            if (i + m["codon"]) % 3 == 0:
                cs = m["codon_start"]
                for k in range(3):
                    if m["ref_codon"][k] != m["alt_codon"][k]:
                        vs.append(Variant(m["chrom"], cs + k, m["ref_codon"][k], m["alt_codon"][k]))
        vs.sort(key=lambda v: (v.chrom, v.pos))
        seq = {c: list(s) for c, s in ref.items()}
        for v in sorted(vs, key=lambda v: (v.chrom, -v.pos)):
            s = seq[v.chrom]
            s[v.pos - 1:v.pos - 1 + len(v.ref)] = list(v.alt)
        clones_seq.append(seq)
        clones_var.append(vs)
    return clones_seq, clones_var


def phred_char(q: int) -> str:
    return chr(q + 33)


def simulate_pairs(rng: random.Random, genomes: list[tuple[float, dict[str, list[str]]]], n_pairs: int,
                   read_len: int, insert_mean: int, insert_sd: int, regions: list[tuple[str, int, int]] | None,
                   name_prefix: str, umi_len: int = 0, family_mean: float = 0.0
                   ) -> tuple[list[str], list[str]]:
    r1: list[str] = []
    r2: list[str] = []
    weights = [w for w, _ in genomes]
    strings = [{c: "".join(s) for c, s in g.items()} for _, g in genomes]
    pid = 0
    while pid < n_pairs:
        gi = rng.choices(range(len(genomes)), weights=weights)[0]
        g = strings[gi]
        if regions:
            chrom, region_start, region_end = rng.choice(regions)
            frag = g[chrom][region_start:region_end]
        else:
            chrom = rng.choice(list(g))
            ins = max(read_len + 20, int(rng.gauss(insert_mean, insert_sd)))
            start = rng.randint(0, len(g[chrom]) - ins - 1)
            frag = g[chrom][start:start + ins]
        if len(frag) < read_len:
            continue
        if rng.random() < 0.5:
            frag = frag.translate(COMP)[::-1]
        umi = "".join(rng.choices("ACGT", k=umi_len)) if umi_len else ""
        copies = 1 + (int(rng.expovariate(1 / family_mean)) if family_mean else 0)
        for _ in range(copies):
            if pid >= n_pairs:
                break
            pid += 1
            for out, s in ((r1, frag[:read_len]), (r2, frag.translate(COMP)[::-1][:read_len])):
                if umi_len:
                    s = umi + s[: read_len - umi_len]
                seq, qual = [], []
                for i, b in enumerate(s):
                    q = max(2, min(40, int(rng.gauss(36 - 8 * i / read_len, 3))))
                    if rng.random() < 10 ** (-q / 10):
                        b = rng.choice([x for x in "ACGT" if x != b])
                    seq.append(b)
                    qual.append(phred_char(q))
                out.append(f"@{name_prefix}:{pid}\n{''.join(seq)}\n+\n{''.join(qual)}\n")
    return r1, r2


def write_fastq(path: Path, records: list[str]) -> None:
    with gzip.open(path, "wt", compresslevel=1) as h:
        h.writelines(records)


def write_fasta(path: Path, seqs: dict[str, list[str]]) -> None:
    with path.open("w") as h:
        for name, s in seqs.items():
            h.write(f">{name}\n")
            seq = "".join(s)
            for i in range(0, len(seq), 60):
                h.write(seq[i:i + 60] + "\n")


FIXTURE_PANEL = ("PFVARX_SYNTHETIC_PANEL", "1")
FIXTURE_PROTOCOL = "SYNTHETIC_UMI8"
FIXTURE_READ_STRUCTURE = "8M+T 8M+T"
FIXTURE_TECHNIQUE = {Design.WGS: "WGS", Design.TARGETED_NO_UMI: "AMPLICON", Design.TARGETED_UMI: "TBS_UMI"}


def build(dest: Path, design: str = "wgs", seed: int = 7, depth: int = 40, n_mono: int = 6,
          n_poly: int = 3, purpose: str = "analysis", low_depth_sample: bool = False,
          variant_scale: int = 1) -> dict[str, Any]:
    if variant_scale < 1:
        raise ValueError("variant_scale must be at least 1")
    design = Design.parse(design)
    targeted = design != Design.WGS
    rng = random.Random(seed)
    lengths = {"ctg1": 30000, "ctg2": 24000}
    ref = make_reference(rng, lengths)
    core = {c: (1500, L - 1500) for c, L in lengths.items()}  # 'subtelomeric' ends masked
    clones, cvars = make_clones(rng, ref, n_mono + 2, core, variant_scale)
    (dest / "inputs/reference").mkdir(parents=True, exist_ok=True)
    (dest / "inputs/metadata").mkdir(parents=True, exist_ok=True)
    (dest / "inputs/fastq").mkdir(parents=True, exist_ok=True)
    write_fasta(dest / "inputs/reference/Pf3D7.fasta", ref)
    with (dest / "inputs/reference/core_genome.bed").open("w") as h:
        for c, (a, b) in core.items():
            h.write(f"{c}\t{a}\t{b}\n")
    with (dest / "inputs/reference/Pf3D7.gff").open("w") as h:
        h.write("##gff-version 3\n")
        for gene, chrom, start, ncod in GENES:
            end = start + ncod * 3 - 1
            h.write(f"{chrom}\tsynthetic\tgene\t{start}\t{end}\t.\t+\t.\tID={gene};Name={gene}\n")
            h.write(f"{chrom}\tsynthetic\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.1;Parent={gene}\n")
            h.write(f"{chrom}\tsynthetic\texon\t{start}\t{end}\t.\t+\t.\tID={gene}.1-E1;Parent={gene}.1\n")
            h.write(f"{chrom}\tsynthetic\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene}.1-CDS;Parent={gene}.1\n")
    regions = None
    if targeted:
        (dest / "inputs/targets").mkdir(parents=True, exist_ok=True)
        amps = []
        for chrom, (a, b) in core.items():
            for k in range(a + 500, b - 500, 1100):
                amps.append((chrom, k, k + 220))
        for m in marker_positions():
            s = m["codon_start"] - 100
            amps.append((m["chrom"], s, s + 220))
        amps = sorted(set(amps))
        regions = [(c, a, b) for c, a, b in amps]
        with (dest / "inputs/targets/targets.bed").open("w") as h:
            for c, a, b in amps:
                h.write(f"{c}\t{a}\t{b}\n")
        with (dest / "inputs/targets/primers.tsv").open("w") as h:
            h.write("chrom\tleft_start\tleft_end\tright_start\tright_end\n")
            for c, a, b in amps:
                h.write(f"{c}\t{a + 1}\t{a + 20}\t{b - 19}\t{b}\n")
    samples = []
    truth: dict[str, Any] = {"clones": {}, "samples": {}}
    for i, vs in enumerate(cvars):
        truth["clones"][f"clone{i}"] = [v.__dict__ for v in vs]
    countries = [("Ethiopia", 9.0, 38.7), ("Ethiopia", 7.0, 39.5), ("Kenya", -0.1, 34.7), ("Kenya", 0.5, 35.3)]
    plan = [(f"PF{i + 1:03d}", [(1.0, i)]) for i in range(n_mono)]
    mixes = [[(0.7, 0), (0.3, 1)], [(0.5, 2), (0.3, 3), (0.2, 4)], [(0.85, n_mono), (0.15, n_mono + 1)]]
    for j in range(n_poly):
        plan.append((f"PF{n_mono + j + 1:03d}", mixes[j % len(mixes)]))
    if low_depth_sample:  # a deliberately under-sequenced sample to exercise QC failure isolation
        plan.append((f"PF{n_mono + n_poly + 1:03d}", [(1.0, 0)]))
    read_len = 100
    for sid, comp in plan:
        genomes = [(w, clones[c]) for w, c in comp]
        d = 3 if low_depth_sample and sid == plan[-1][0] else depth
        if not targeted:
            n_pairs = int(d * sum(lengths.values()) / (2 * read_len))
        else:
            n_pairs = int(d * 2 * len(regions or []) * 220 / (2 * read_len))
        umi = 8 if design == Design.TARGETED_UMI else 0
        fam = 3.0 if design == Design.TARGETED_UMI else 0.0
        r1, r2 = simulate_pairs(rng, genomes, n_pairs, read_len + umi if umi else read_len, 300, 40,
                                regions, sid, umi, fam)
        write_fastq(dest / f"inputs/fastq/{sid}_L001_R1.fastq.gz", r1)
        write_fastq(dest / f"inputs/fastq/{sid}_L001_R2.fastq.gz", r2)
        idx = int(sid[2:]) - 1
        pop = "POP_A" if comp[0][1] % 2 == 0 else "POP_B"
        country, lat, lon = countries[idx % len(countries)]
        samples.append({"sample_id": sid, "lane": "L001", "fastq_1": f"inputs/fastq/{sid}_L001_R1.fastq.gz",
                        "fastq_2": f"inputs/fastq/{sid}_L001_R2.fastq.gz", "bam": "", "vcf": "",
                        "library_id": f"{sid}_lib1", "flowcell_type": "unpatterned", "population": pop,
                        "technique": FIXTURE_TECHNIQUE[design], "amplification_route": "AMPLICON" if targeted else "DIRECT",
                        "panel_id": FIXTURE_PANEL[0] if targeted else "", "panel_version": FIXTURE_PANEL[1] if targeted else "",
                        "library_protocol_id": FIXTURE_PROTOCOL if design == Design.TARGETED_UMI else "",
                        "processing_state": "", "adapter_kit": "",
                        "country": country, "latitude": f"{lat + 0.01 * idx:.3f}",
                        "longitude": f"{lon + 0.01 * idx:.3f}",
                        "collection_date": f"2025-{(idx % 12) + 1:02d}-15"})
        truth["samples"][sid] = {"components": [{"clone": f"clone{c}", "proportion": w} for w, c in comp],
                                 "true_k": len(comp),
                                 "clonality": "MONOCLONAL" if len(comp) == 1 else "POLYCLONAL"}
    cols = ["sample_id", "lane", "fastq_1", "fastq_2", "bam", "vcf", "library_id", "flowcell_type",
            "technique", "amplification_route", "panel_id", "panel_version", "library_protocol_id",
            "processing_state", "adapter_kit", "population", "country", "latitude", "longitude", "collection_date"]
    if targeted:
        # Library catalog: the synthetic panel and, for UMI libraries, the UMI protocol.
        catalog_doc = {"schema_version": CATALOG_SCHEMA,
                       "panels": [{"panel_id": FIXTURE_PANEL[0], "panel_version": FIXTURE_PANEL[1],
                                   "targets_bed": "inputs/targets/targets.bed",
                                   "primers_tsv": "inputs/targets/primers.tsv",
                                   "description": "synthetic amplicon panel (software test)"}],
                       "library_protocols": ([{"library_protocol_id": FIXTURE_PROTOCOL, "technique": "TBS_UMI",
                                               "amplification_route": "AMPLICON",
                                               "umi": {"read_structure": FIXTURE_READ_STRUCTURE},
                                               "description": "8-base inline UMI on both reads"}]
                                             if design == Design.TARGETED_UMI else [])}
        (dest / "inputs/metadata/library_catalog.yaml").write_text(yaml.safe_dump(catalog_doc, sort_keys=False))
    # The sample sheet goes where the workspace configuration expects it (prepared/samples.tsv by default).
    config = dest / "config/config.yaml"
    sheet_path = dest / ((yaml.safe_load(config.read_text()) or {}).get("samples_tsv", "prepared/samples.tsv")
                         if config.is_file() else "prepared/samples.tsv")
    sheet_path.parent.mkdir(parents=True, exist_ok=True)
    with sheet_path.open("w") as h:
        h.write("\t".join(cols) + "\n")
        for s in samples:
            h.write("\t".join(s[c] for c in cols) + "\n")
    catalog = dest / "inputs/metadata/test_resistance_catalog.tsv"
    with catalog.open("w") as h:
        h.write("#catalog_id=pfvarx_synthetic_test_catalog\n#version=1.0.0\n#status=TEST_FIXTURE\n")
        h.write("marker_id\tgene_id\tgene_name\tchrom\tcodon\tcodon_start\tstrand\tref_aa\talt_aa\tevidence_class\tsource\n")
        for m in marker_positions():
            h.write(f"{m['gene']}_{m['ref_aa']}{m['codon']}{m['alt_aa']}\t{m['gene']}\t{m['gene']}\t{m['chrom']}\t"
                    f"{m['codon']}\t{m['codon_start']}\t+\t{m['ref_aa']}\t{m['alt_aa']}\tTEST_FIXTURE\tsynthetic\n")
    if purpose == "benchmarking":
        (dest / "inputs/truth").mkdir(parents=True, exist_ok=True)
        (dest / "inputs/truth/truth.json").write_text(json.dumps(truth, indent=1))
        with (dest / "inputs/truth/confident_regions.bed").open("w") as h:
            for c, (a, b) in core.items():
                h.write(f"{c}\t{a}\t{b}\n")
        for sid, info in truth["samples"].items():
            if info["clonality"] == "MONOCLONAL":
                clone = int(info["components"][0]["clone"][5:])
                lines = (["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">',
                          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
                         + [f"##contig=<ID={c},length={L}>" for c, L in lengths.items()])
                lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sid)
                for v in cvars[clone]:
                    lines.append(f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\t.\tGT\t1")
                (dest / f"inputs/truth/{sid}.truth.vcf").write_text("\n".join(lines) + "\n")
    return {"samples": [s["sample_id"] for s in samples], "design": design.value, "truth": truth["samples"]}
