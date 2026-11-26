"""Layer 3 multiplexed microhaplotypes: read-backed haplotypes at catalogue loci.

A microhaplotype is the combination of alleles a single DNA fragment carries across the
variable positions of one short locus (an amplicon, or a window of the core genome). Read
together they form one multiallelic marker instead of several biallelic ones, which is what
MOIRE (DA001) and Dcifer (DA004) are built for: a locus with many haplotypes separates clones
far better than the same sites treated one at a time.

How a locus is called, for every assay design:

* **Loci** come from a maintainer catalogue (``analyses.microhap_panel``, else the shipped
  registry panel), given as coordinates or as a gene identifier placed with the GFF. A locus is
  used only where the run can speak about it: inside the panel's targets for a targeted design,
  inside the core mask for whole-genome sequencing. Anything else is recorded, not silently
  dropped.
* **Variable positions** are this cohort's PASS SNVs inside the locus (``cohort/sites.tsv``).
  No variant is discovered here: the caller's decisions are reused, so a haplotype cannot
  contain a site the caller did not pass. Indels are left out, because they move the positions
  of everything after them within a read.
* **Fragments, not reads.** Observations are grouped by read name, so a pair whose mates overlap
  counts once; mates that disagree at a position make that position unresolved. A fragment that
  does not cover every variable position with a good enough base is ``incomplete`` and counted
  as such, never guessed.
* **Floors.** A sample needs ``microhap_min_depth`` complete fragments, and a haplotype needs
  ``microhap_min_haplotype_reads`` fragments and ``microhap_min_haplotype_fraction`` of them;
  below that it is treated as sequencing error. A locus that cannot be judged is UNCALLABLE for
  that sample and is never read as "no haplotype".

With UMI designs the fragments are consensus molecules, so the counts are molecules; the unit is
recorded in the summary. The shipped panel contains no loci: published panel coordinates are not
invented here, so a site supplies its own verified panel and the artifact is empty until then
(DA001 and DA004 then use biallelic SNVs, as before, and record that they did).
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pysam

from pfvarx.configuration.loader import frozen_param
from pfvarx.constants import CountUnit, Design
from pfvarx.genomic_evidence import FLAG_FILTER, in_intervals, read_bed

SCHEMA = "pfvarx.microhaplotypes/1"
PANEL = "config/registry/microhap_panel_pfvarx_v1.tsv"
PANEL_COLUMNS = ("locus_id", "chrom", "start", "end", "gene_id", "name", "source")
HAPLOTYPE_COLUMNS = ("locus_id", "sample_id", "haplotype", "fragments", "fraction", "state")
LOCUS_COLUMNS = ("locus_id", "chrom", "start", "end", "name", "status", "reason", "variable_positions",
                 "positions", "samples_called", "samples_uncallable", "haplotypes", "source")
UNRESOLVED = "N"


def panel_path(root: Path, configured: str | None) -> Path:
    return root / (configured or PANEL)


def resolve_locus(row: dict[str, str], gff: Path | None, lengths: dict[str, int]
                  ) -> tuple[tuple[str, int, int] | None, str]:
    """(contig, 0-based start, end) of a catalogue locus, with the reason when it cannot be placed.

    Coordinates are 1-based inclusive in the catalogue, as in a BED-less genomic range; a locus
    may instead name a gene, placed with the GFF like the resistance and copy-number catalogues.
    """
    chrom, start, end = row.get("chrom", ""), row.get("start", ""), row.get("end", "")
    if chrom and start and end:
        try:
            first, last = int(start), int(end)
        except ValueError:
            return None, "COORDINATES_MALFORMED"
        if chrom not in lengths:
            return None, "CONTIG_NOT_IN_REFERENCE"
        if not 0 < first <= last <= lengths[chrom]:
            return None, "COORDINATES_OUT_OF_RANGE"
        return (chrom, first - 1, last), ""
    gene = row.get("gene_id", "")
    if not gene:
        return None, "NO_COORDINATES_OR_GENE"
    if gff is None or not gff.is_file():
        return None, "GENE_LOCUS_NEEDS_GFF"
    from pfvarx.downstream.cnv import gene_location

    place = gene_location(gff, gene)
    return (place, "") if place else (None, "GENE_NOT_IN_GFF")


def territory_status(locus: tuple[str, int, int], territory: dict[str, list[tuple[int, int]]],
                     design: Design) -> tuple[bool, str]:
    """Whether the whole locus lies in what the run sequenced, and why not when it does not."""
    chrom, start, end = locus
    spans = territory.get(chrom, [])
    inside = any(s <= start and end <= e for s, e in spans)
    if inside:
        return True, ""
    return False, "OUTSIDE_PANEL_TARGETS" if design.targeted else "OUTSIDE_CORE_MASK"


def variable_positions(sites: list[dict[str, str]], locus: tuple[str, int, int], cap: int
                       ) -> tuple[list[int], str]:
    """This cohort's PASS SNV positions (0-based) inside the locus, and why they were limited."""
    chrom, start, end = locus
    found = sorted({int(s["pos"]) - 1 for s in sites
                    if s["chrom"] == chrom and start <= int(s["pos"]) - 1 < end
                    and len(s["ref"]) == 1 and len(s["alt"]) == 1 and int(s["n_pass"]) > 0})
    if len(found) > cap:
        return found[:cap], "VARIABLE_POSITIONS_CAPPED"
    return found, ""


def derive_loci(sites: list[dict[str, str]], territory: dict[str, list[tuple[int, int]]], design: Design,
                window_bp: int, min_sites: int, max_loci: int) -> tuple[list[dict[str, str]], str, int]:
    """Microhaplotype loci from the run's own inputs, when no catalogue is supplied.

    A targeted panel already defines them: one amplicon is one locus, so each target region
    holding at least *min_sites* of this cohort's PASS SNVs becomes a locus. Whole-genome runs
    have no such definition, so the core mask is tiled into windows of *window_bp* and the
    windows with enough PASS SNVs are used.

    Nothing is invented: every coordinate comes from the panel the run was given or from the core
    mask it declares. The loci are specific to this run's panel or window length, so haplotype
    frequencies are comparable with another cohort only if that cohort used the same definition -
    which is why a published panel's own catalogue remains the better input when one exists.
    """
    by_contig: dict[str, list[int]] = {}
    for site in sites:
        if len(site["ref"]) == 1 and len(site["alt"]) == 1 and int(site["n_pass"]) > 0:
            by_contig.setdefault(site["chrom"], []).append(int(site["pos"]) - 1)
    for positions in by_contig.values():
        positions.sort()
    candidates: list[tuple[str, int, int]] = []
    if design.targeted:
        source = "DERIVED_PANEL_TARGETS"
        for chrom, spans in sorted(territory.items()):
            for start, end in spans:
                candidates.append((chrom, start, end))
    else:
        source = "DERIVED_CORE_WINDOWS"
        for chrom, spans in sorted(territory.items()):
            for start, end in spans:
                for window_start in range(start, end, window_bp):
                    candidates.append((chrom, window_start, min(end, window_start + window_bp)))
    loci: list[dict[str, str]] = []
    for chrom, start, end in candidates:
        # The positions are sorted, so each window's count is two binary searches. A whole genome
        # tiles into of the order of 100 000 windows against as many sites, and a scan per window
        # would not finish.
        positions = by_contig.get(chrom, [])
        inside = bisect.bisect_left(positions, end) - bisect.bisect_left(positions, start)
        if inside >= min_sites:
            loci.append({"locus_id": f"derived_{chrom}_{start + 1}_{end}", "chrom": chrom,
                         "start": str(start + 1), "end": str(end), "gene_id": "",
                         "name": f"{inside} PASS SNVs", "source": source.lower()})
    capped = max(0, len(loci) - max_loci)
    return loci[:max_loci], source, capped


def fragment_haplotypes(bam: pysam.AlignmentFile, chrom: str, positions: list[int],
                        min_mapq: int, min_bq: int) -> tuple[Counter[str], int]:
    """Haplotype counts over *positions* by fragment, and the number of incomplete fragments.

    Both mates of a pair contribute to one fragment: an overlap counts once and a disagreement
    leaves that position unresolved, which makes the fragment incomplete.
    """
    observed: dict[str, dict[int, str]] = {}
    wanted = set(positions)
    for read in bam.fetch(chrom, min(positions), max(positions) + 1):
        if read.flag & FLAG_FILTER or read.mapping_quality < min_mapq:
            continue
        if read.query_sequence is None or read.query_qualities is None:
            continue
        if read.query_name is None:
            continue
        fragment = observed.setdefault(read.query_name, {})
        for query_index, reference_position in read.get_aligned_pairs(matches_only=True):
            if reference_position not in wanted or query_index is None:
                continue
            if read.query_qualities[query_index] < min_bq:
                continue
            base = read.query_sequence[query_index].upper()
            seen = fragment.get(reference_position)
            fragment[reference_position] = base if seen is None or seen == base else UNRESOLVED
    counts: Counter[str] = Counter()
    incomplete = 0
    for fragment in observed.values():
        haplotype = "".join(fragment.get(position, UNRESOLVED) for position in positions)
        if UNRESOLVED in haplotype:
            incomplete += 1
        else:
            counts[haplotype] += 1
    return counts, incomplete


def call_sample(counts: Counter[str], params: dict[str, Any]) -> tuple[str, str, list[tuple[str, int, float]]]:
    """(state, reason, kept haplotypes) for one sample at one locus."""
    total = sum(counts.values())
    if total < int(frozen_param(params, "microhap_min_depth")):
        return "UNCALLABLE", "LOW_FRAGMENT_DEPTH", []
    min_reads = int(frozen_param(params, "microhap_min_haplotype_reads"))
    min_fraction = float(frozen_param(params, "microhap_min_haplotype_fraction"))
    kept = [(haplotype, n, n / total) for haplotype, n in counts.most_common()
            if n >= min_reads and n / total >= min_fraction]
    if not kept:
        return "UNCALLABLE", "NO_HAPLOTYPE_ABOVE_FLOOR", []
    return "CALLED", "", kept


def build(sites: list[dict[str, str]], panel: list[dict[str, str]], samples: list[str], bams: dict[str, Path],
          callable_beds: dict[str, Path], lengths: dict[str, int], territory: dict[str, list[tuple[int, int]]],
          design: Design, gff: Path | None, params: dict[str, Any], min_mapq: int, min_bq: int
          ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Haplotype rows and locus rows for the whole cohort."""
    cap = int(frozen_param(params, "microhap_max_variable_positions"))
    haplotype_rows: list[dict[str, Any]] = []
    locus_rows: list[dict[str, Any]] = []
    for row in panel:
        locus_id = row["locus_id"]
        place, reason = resolve_locus(row, gff, lengths)
        record: dict[str, Any] = {"locus_id": locus_id, "name": row.get("name", ""),
                                  "source": row.get("source", ""), "chrom": "", "start": "", "end": "",
                                  "status": "EXCLUDED", "reason": reason, "variable_positions": 0,
                                  "positions": "", "samples_called": 0, "samples_uncallable": 0, "haplotypes": 0}
        if place is None:
            locus_rows.append(record)
            continue
        chrom, start, end = place
        record.update(chrom=chrom, start=start + 1, end=end)
        inside, why = territory_status(place, territory, design)
        if not inside:
            record["reason"] = why
            locus_rows.append(record)
            continue
        positions, capped = variable_positions(sites, place, cap)
        record.update(variable_positions=len(positions),
                      positions=",".join(str(p + 1) for p in positions))
        if not positions:
            record["reason"] = "NO_VARIABLE_SITES_IN_COHORT"
            locus_rows.append(record)
            continue
        record["status"] = "USED"
        record["reason"] = capped
        haplotypes: set[str] = set()
        for sample in samples:
            with pysam.AlignmentFile(str(bams[sample])) as bam:
                if chrom not in bam.references:
                    haplotype_rows.append({"locus_id": locus_id, "sample_id": sample, "haplotype": "",
                                           "fragments": 0, "fraction": "", "state": "UNCALLABLE",
                                           "reason": "CONTIG_NOT_IN_ALIGNMENT"})
                    record["samples_uncallable"] += 1
                    continue
                intervals = read_bed(callable_beds[sample]) if sample in callable_beds else {}
                if intervals and not all(in_intervals(intervals.get(chrom, []), p) for p in positions):
                    haplotype_rows.append({"locus_id": locus_id, "sample_id": sample, "haplotype": "",
                                           "fragments": 0, "fraction": "", "state": "UNCALLABLE",
                                           "reason": "POSITION_NOT_CALLABLE"})
                    record["samples_uncallable"] += 1
                    continue
                counts, incomplete = fragment_haplotypes(bam, chrom, positions, min_mapq, min_bq)
            state, reason, kept = call_sample(counts, params)
            if state != "CALLED":
                haplotype_rows.append({"locus_id": locus_id, "sample_id": sample, "haplotype": "",
                                       "fragments": sum(counts.values()), "fraction": "", "state": state,
                                       "reason": reason, "incomplete_fragments": incomplete})
                record["samples_uncallable"] += 1
                continue
            record["samples_called"] += 1
            for haplotype, n, fraction in kept:
                haplotypes.add(haplotype)
                haplotype_rows.append({"locus_id": locus_id, "sample_id": sample, "haplotype": haplotype,
                                       "fragments": n, "fraction": round(fraction, 4), "state": "CALLED",
                                       "reason": "", "incomplete_fragments": incomplete})
        record["haplotypes"] = len(haplotypes)
        locus_rows.append(record)
    return haplotype_rows, locus_rows


def write_haplotypes(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = (*HAPLOTYPE_COLUMNS, "reason", "incomplete_fragments")
    with gzip.open(path, "wt") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join("" if row.get(c) is None else str(row.get(c, "")) for c in columns) + "\n")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sites", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--bams", nargs="+", required=True)
    ap.add_argument("--callable", nargs="*", default=[])
    ap.add_argument("--reference", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--frozen-config", required=True)
    ap.add_argument("--out-haplotypes", required=True)
    ap.add_argument("--out-loci", required=True)
    ap.add_argument("--out-summary", required=True)
    a = ap.parse_args(argv)
    from pfvarx.io import read_tsv, write_json, write_tsv

    root = Path(a.root)
    eff = json.loads(Path(a.frozen_config).read_text())
    plan = json.loads((Path(a.frozen_config).parent / "frozen_plan.json").read_text())
    params = {k.split(".", 1)[1]: v["value"] for k, v in eff["parameters"].items() if k.startswith("layer3.")}
    layer2 = {k.split(".", 1)[1]: v["value"] for k, v in eff["parameters"].items() if k.startswith("layer2.")}
    min_mapq, min_bq = int(layer2["min_mapping_quality"]), int(layer2["min_base_quality"])
    design = Design.parse(plan.get("design") or eff["config"].get("assay_design"))
    samples = Path(a.samples).read_text().split()
    bams = {Path(b).name.split(".")[0]: Path(b) for b in a.bams}
    callable_beds = {Path(c).name.split(".")[0]: Path(c) for c in a.callable}
    mode = str(eff["config"]["analyses"].get("microhap_loci") or "off")
    path = panel_path(root, eff["config"]["analyses"].get("microhap_panel"))
    from pfvarx.downstream.resistance import load_catalog

    meta: dict[str, str] = {}
    panel: list[dict[str, str]] = []
    if mode == "catalogue":
        meta, panel = load_catalog(path) if path.is_file() else ({}, [])
    with pysam.FastaFile(a.reference) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
    gff = root / eff["config"]["reference"]["gff"] if eff["config"]["reference"].get("gff") else None
    if design.targeted:
        panel_targets = (plan.get("panel") or {}).get("targets_bed")
        territory = read_bed(root / panel_targets) if panel_targets else {}
    else:
        territory = read_bed(root / eff["config"]["reference"]["core_mask_bed"])
    sites = read_tsv(Path(a.sites))
    source, capped = ("CATALOGUE" if mode == "catalogue" else "OFF"), 0
    if mode == "derived":
        panel, source, capped = derive_loci(sites, territory, design,
                                            int(frozen_param(params, "microhap_derived_window_bp")),
                                            int(frozen_param(params, "microhap_derived_min_sites")),
                                            int(frozen_param(params, "microhap_derived_max_loci")))
        meta = {"catalog_id": "derived", "version": source.lower(), "status": "DERIVED_FROM_RUN_INPUTS"}
    status = meta.get("status", "NO_PANEL" if not panel else "UNKNOWN")
    # Derived loci need no maintainer verification: their coordinates are the run's own panel
    # targets or core mask, not a claim about a published panel.
    verified = status in ("MAINTAINER_VERIFIED", "TEST_FIXTURE", "DERIVED_FROM_RUN_INPUTS")
    haplotype_rows: list[dict[str, Any]] = []
    locus_rows: list[dict[str, Any]] = []
    if panel and verified:
        haplotype_rows, locus_rows = build(sites, panel, samples, bams, callable_beds, lengths, territory,
                                          design, gff, params, min_mapq, min_bq)
    write_haplotypes(Path(a.out_haplotypes), haplotype_rows)
    write_tsv(Path(a.out_loci), locus_rows, list(LOCUS_COLUMNS))
    used = [r for r in locus_rows if r["status"] == "USED"]
    reason = ("" if panel and verified else
              "MICROHAPLOTYPES_OFF" if mode == "off" else
              ("NO_DERIVED_LOCUS_HAS_ENOUGH_SITES" if mode == "derived" else "NO_LOCI_IN_PANEL")
              if not panel else "PANEL_NOT_VERIFIED")
    write_json(Path(a.out_summary), {
        "schema": SCHEMA, "panel": str(path) if mode == "catalogue" else None, "panel_meta": meta,
        "panel_status": status, "panel_verified": verified, "reason": reason,
        "locus_source": source, "mode": mode, "loci_over_cap": capped,
        "comparability": ("Loci come from this run's own panel targets or core mask, so haplotype "
                          "frequencies are comparable only with cohorts using the same definition."
                          if source.startswith("DERIVED") else
                          "Loci come from the supplied catalogue." if source == "CATALOGUE" else ""),
        "count_unit": (CountUnit.CONSENSUS_MOLECULES if design == Design.TARGETED_UMI else CountUnit.READS).value,
        "design": design.value, "loci_in_panel": len(panel), "loci_used": len(used),
        "loci_excluded": {r["locus_id"]: r["reason"] for r in locus_rows if r["status"] != "USED"},
        "samples": len(samples),
        "haplotypes": len({(r["locus_id"], r["haplotype"]) for r in haplotype_rows if r["state"] == "CALLED"}),
        "mean_haplotypes_per_locus": (round(sum(r["haplotypes"] for r in used) / len(used), 2) if used else None),
        "parameters": {k: v for k, v in params.items() if k.startswith("microhap_")},
        "admission": {"min_mapping_quality": min_mapq, "min_base_quality": min_bq},
    })


if __name__ == "__main__":
    main()
