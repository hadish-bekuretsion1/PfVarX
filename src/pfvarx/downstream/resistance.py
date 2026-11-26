"""DA002 drug-resistance markers with callable denominators.

Marker genotypes are read directly from each sample's alignment: reads that span all
three codon bases are translated, and the fraction encoding the catalogue's alternative
amino acid is the marker's within-sample fraction. An uncallable codon is UNKNOWN,
never WILD_TYPE. The catalogue must declare status MAINTAINER_VERIFIED, or TEST_FIXTURE
(outputs then carry a test-fixture warning); any other status blocks DA002.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pysam

from pfvarx.constants import MarkerState, Reason, Terminal
from pfvarx.downstream.common import Context, decisions, wilson, write_manifest
from pfvarx.genomic_evidence import FLAG_FILTER
from pfvarx.io import write_json, write_tsv

# Layer 5 memory: the cohort evidence matrix is loaded for this analysis only when it reads it
# (reads the annotated VCF and each sample's reads, not the matrix). tests/unit/test_matrix_declarations.py checks the declaration against the code.
USES_COHORT_MATRIX = False
CODE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met", "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop", "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys", "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp", "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg", "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}
COMP = str.maketrans("ACGT", "TGCA")


def load_catalog(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    meta: dict[str, str] = {}
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    for line in path.read_text().splitlines():
        if line.startswith("#") and "=" in line:
            k, v = line[1:].split("=", 1)
            meta[k.strip()] = v.strip()
        elif line.strip() and header is None:
            header = line.split("\t")
        elif line.strip():
            rows.append(dict(zip(header or [], line.split("\t"))))
    return meta, rows


def catalog_path(ctx: Context) -> Path:
    return ctx.root / (ctx.eff["config"]["analyses"].get("resistance_catalog")
                       or "config/registry/resistance_catalog_pfvarx_v1.tsv")


def codon_positions(ctx: Context, markers: list[dict[str, str]]) -> list[list[int]]:
    """0-based genomic positions of each marker codon (empty when it cannot be placed).
    Sets the marker's chrom and strand when they come from the GFF."""
    gff = ctx.root / ctx.eff["config"]["reference"]["gff"] if ctx.eff["config"]["reference"].get("gff") else None
    out: list[list[int]] = []
    for m in markers:
        if m.get("codon_start"):
            st = int(m["codon_start"]) - 1
            out.append([st, st + 1, st + 2] if m.get("strand", "+") == "+" else [st + 2, st + 1, st])
        elif gff is not None and (res := resolve_codon(gff, m["gene_id"], int(m["codon"]))):
            m["chrom"], positions, m["strand"] = res
            out.append(positions)
        else:
            out.append([])
    return out


def resolve_codon(gff: Path, gene_id: str, codon: int) -> tuple[str, list[int], str] | None:
    """Genomic 0-based positions of a codon from GFF3 CDS features (handles strand and exons)."""
    cds: list[tuple[str, int, int, str]] = []
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        parent = attrs.get("Parent", "")
        if parent.split(".")[0] == gene_id or parent == f"{gene_id}.1":
            cds.append((c[0], int(c[3]) - 1, int(c[4]), c[6]))
    if not cds:
        return None
    strand = cds[0][3]
    bases: list[int] = []
    for chrom, s0, e, _ in sorted(cds, key=lambda x: x[1]):
        bases.extend(range(s0, e))
    if strand == "-":
        bases = bases[::-1]
    i = (codon - 1) * 3
    if i + 3 > len(bases):
        return None
    return cds[0][0], bases[i:i + 3], strand


def codon_reads(bam: pysam.AlignmentFile, chrom: str, pos: list[int], strand: str, min_mapq: int,
                min_bq: int) -> list[str]:
    """Codons (in transcript orientation) observed on fragments covering all three codon bases."""
    out = []
    seen: set[str] = set()
    for read in bam.fetch(chrom, min(pos), max(pos) + 1):
        if read.flag & FLAG_FILTER or read.mapping_quality < min_mapq or read.query_name in seen:
            continue
        if read.query_name is None or read.query_sequence is None or read.query_qualities is None:
            continue
        pairs = {r: q for q, r in read.get_aligned_pairs(matches_only=True) if r is not None and q is not None}
        if not all(p in pairs for p in pos):
            continue
        quals = read.query_qualities
        if any(quals[pairs[p]] < min_bq for p in pos):
            continue
        seen.add(read.query_name)  # one observation per fragment
        codon = "".join(read.query_sequence[pairs[p]] for p in pos).upper()
        out.append(codon if strand == "+" else codon.translate(COMP))
    return out


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA002")
    dec = decisions(ctx)["analyses"]["DA002"]
    eligible = sorted({s for c in dec["scales"].values() for s in c["eligible"]})
    excluded = {"INELIGIBLE": sorted(set(ctx.samples) - set(eligible))}
    cat_path = catalog_path(ctx)
    meta, markers = load_catalog(cat_path)
    status = meta.get("status", "UNKNOWN")
    if status not in ("MAINTAINER_VERIFIED", "TEST_FIXTURE"):
        write_manifest(ctx, "DA002", Terminal.BLOCKED_MISSING_PREREQUISITE.value,
                       reason_codes=[Reason.CATALOG_NOT_VERIFIED.value], requested=ctx.samples, eligible=eligible,
                       excluded=excluded, analysed=[], outputs=[],
                       caveats=[f"catalogue {cat_path.name} has status {status}; verify it against the WHO compendium"],
                       method={"catalog": str(cat_path), "catalog_meta": meta})
        return
    p = ctx.params
    min_depth = int(ctx.eff["parameters"]["layer2.min_usable_depth"]["value"])
    min_mapq = int(ctx.eff["parameters"]["layer2.min_mapping_quality"]["value"])
    min_bq = int(ctx.eff["parameters"]["layer2.min_base_quality"]["value"])
    floor, mixed_thr = float(p["da002_reporting_floor"]), float(p["da002_mixed_threshold"])
    marker_positions = codon_positions(ctx, markers)
    calls: list[dict[str, Any]] = []
    for sid in eligible:
        with pysam.AlignmentFile(str(ctx.results / "alignments" / f"{sid}.bam")) as bam:
            for m, positions in zip(markers, marker_positions):
                codons = (codon_reads(bam, m["chrom"], positions, m.get("strand", "+"), min_mapq, min_bq)
                          if positions else [])
                n = len(codons)
                aas = [CODE.get(c, "X") for c in codons]
                alt = sum(1 for a in aas if a == m["alt_aa"])
                refn = sum(1 for a in aas if a == m["ref_aa"])
                if n < min_depth:
                    state, frac = MarkerState.UNKNOWN, None
                else:
                    frac = alt / n
                    if frac >= mixed_thr:
                        state = MarkerState.MUTANT
                    elif frac >= floor:
                        state = MarkerState.MIXED
                    elif alt > 0:
                        state = MarkerState.BELOW_REPORTING_FLOOR
                    else:
                        state = MarkerState.WILD_TYPE
                lo, hi = wilson(alt, n) if n else (None, None)
                calls.append({"sample_id": sid, "marker_id": m["marker_id"], "gene": m["gene_name"],
                              "codon": m["codon"], "ref_aa": m["ref_aa"], "alt_aa": m["alt_aa"], "state": state.value,
                              "codon_reads": n, "alt_reads": alt, "ref_reads": refn, "other_reads": n - alt - refn,
                              "alt_fraction": "" if frac is None else f"{frac:.4f}",
                              "fraction_lo": "" if lo is None else f"{lo:.4f}", "fraction_hi": "" if hi is None else f"{hi:.4f}",
                              "fraction_correction": "UNCORRECTED",
                              "population": ctx.metadata.get(sid, {}).get(ctx.eff["config"]["analyses"].get("population_column") or "population", ""),
                              "country": ctx.metadata.get(sid, {}).get("country", ""),
                              "evidence_class": m.get("evidence_class", "")})
    cols = list(calls[0]) if calls else ["sample_id", "marker_id", "state"]
    write_tsv(out / "marker_calls.tsv", calls, cols)
    prev = prevalence(calls, markers, "population", float(p["da002_k13_prevalence_alert"]),
                      int(p["da002_min_callable_per_group"]))
    write_tsv(out / "prevalence.tsv", prev, list(prev[0]) if prev else ["marker_id"])
    caveats = ["Molecular markers are not evidence of an individual patient's treatment outcome.",
               "An UNKNOWN marker is uncallable; it is never counted as wild type.",
               "Within-sample fractions are DNA fractions and UNCORRECTED for representation bias.",
               "Reporting floor (5%) and MIXED threshold (95%) follow the MalariaGEN Pf7 read-fraction rule; "
               "they are not calibrated for this caller's detection limit.",
               "who_k13_alert marks validated PfKelch13 markers above the WHO 5% level for follow-up; "
               "confirmation needs delayed-clearance evidence from a therapeutic efficacy study. "
               f"Groups with fewer than {int(p['da002_min_callable_per_group'])} callable infections are "
               "INSUFFICIENT_CALLABLE."]
    if status == "TEST_FIXTURE":
        caveats.insert(0, "TEST FIXTURE CATALOGUE: results are for software testing, not for surveillance.")
    write_json(out / "catalog_provenance.json", {"path": str(cat_path), "meta": meta, "markers": len(markers)})
    unknown = sum(1 for c in calls if c["state"] == "UNKNOWN")
    write_manifest(ctx, "DA002", (Terminal.COMPLETED_WITH_EXCLUSIONS if excluded["INELIGIBLE"] or unknown
                                  else Terminal.COMPLETED).value,
                   reason_codes=[], requested=ctx.samples, eligible=eligible, excluded=excluded, analysed=eligible,
                   outputs=[out / "marker_calls.tsv", out / "prevalence.tsv", out / "catalog_provenance.json"],
                   caveats=caveats, method={"catalog": meta, "min_codon_reads": min_depth,
                           "count_unit": sorted({(s.get("assay") or {}).get("count_unit", "reads") for s in ctx.plan["samples"]}), "reporting_floor": floor,
                                            "mixed_threshold": mixed_thr, "interval": "Wilson (simple random sampling)"},
                   results_summary={"marker_calls": len(calls), "unknown_calls": unknown})


def is_validated_k13(marker: dict[str, str]) -> bool:
    gene = (marker.get("gene_name", "") + " " + marker.get("gene_id", "")).lower()
    return ("k13" in gene or "kelch13" in gene or "pf3d7_1343700" in gene) and "valid" in marker.get("evidence_class", "").lower()


def prevalence(calls: list[dict[str, Any]], markers: list[dict[str, str]], group_key: str,
               k13_alert: float | None = None, min_callable: int = 1) -> list[dict[str, Any]]:
    """Marker prevalence per group. ``k13_alert`` flags validated PfKelch13 markers above the WHO level.

    The flag needs at least *min_callable* callable infections in the group: one mutant among
    one or two callable samples is a prevalence of 50-100% and would otherwise raise the WHO
    flag for a whole population. Below the minimum the flag is INSUFFICIENT_CALLABLE; the
    prevalence and its Wilson interval are still reported."""
    rows = []
    groups = sorted({c[group_key] for c in calls}) + ["ALL"]
    for m in markers:
        for g in groups:
            sel = [c for c in calls if c["marker_id"] == m["marker_id"] and (g == "ALL" or c[group_key] == g)]
            callable_ = [c for c in sel if c["state"] != "UNKNOWN"]
            k = sum(1 for c in callable_ if c["state"] in ("MUTANT", "MIXED"))
            n = len(callable_)
            lo, hi = wilson(k, n)
            fr = [float(c["alt_fraction"]) for c in callable_ if c["alt_fraction"]]
            rows.append({"marker_id": m["marker_id"], "group": g, "samples": len(sel), "callable": n,
                         "unknown": len(sel) - n, "mutant_or_mixed": k,
                         "prevalence_infections": "" if not n else f"{k / n:.4f}",
                         "wilson_lo": "" if lo is None else f"{lo:.4f}", "wilson_hi": "" if hi is None else f"{hi:.4f}",
                         "prevalence_clones_dna_weighted": "" if not fr else f"{sum(fr) / len(fr):.4f}",
                         # WHO molecular criterion for suspected artemisinin partial resistance; a flag
                         # for follow-up (therapeutic efficacy study), never a resistance decision.
                         "who_k13_alert": ("" if k13_alert is None or not n or not is_validated_k13(m)
                                           else "INSUFFICIENT_CALLABLE" if n < min_callable
                                           else "ABOVE_5_PERCENT" if k / n > k13_alert else "NOT_ABOVE")})
    return rows
