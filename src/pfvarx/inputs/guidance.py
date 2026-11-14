"""Plain-language explanations and fixes for input-validation failures.

Every failed check in ``metadata/validation_report.json`` (and in the ``validate`` command
and GUI output) gets: which input it is, what PfVarX uses it for, why it failed and the
steps that fix it. The checks themselves live in ``validate.py``; this module only explains.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

# kind -> (input name, configuration key, what the pipeline uses it for, accepted content kinds)
ROLES: dict[str, tuple[str, str, str, tuple[str, ...]]] = {
    "FASTA": ("Reference genome FASTA", "reference.fasta",
              "Reads are aligned to it and every variant is described against it.", ("FASTA",)),
    "CORE_MASK": ("Core-genome mask BED", "reference.core_mask_bed",
                  "Variants are called, and depth and callability QC are measured, only inside these regions: the "
                  "accessible core of the genome, without the subtelomeric and internal hypervariable regions "
                  "(var, rifin, stevor) where reads map unreliably.", ("BED",)),
    "GFF": ("Gene annotation GFF", "reference.gff",
            "SnpEff builds its annotation database from it, and drug-resistance codons (DA002) are located with it.",
            ("GFF",)),
    "TARGETS": ("Target regions BED", "assay.targets_bed",
                "Defines the amplicons of a targeted panel for on-target and dropout QC.", ("BED",)),
    "PRIMERS": ("Primer table", "assay.primers_tsv",
                "Primer positions are clipped from targeted reads before calling.", ("PRIMERS",)),
    "SAMPLE_METADATA": ("Sample sheet", "samples_tsv",
                        "Lists every sample and lane with its input files and study metadata.", ("SAMPLE_SHEET",)),
    "CATALOGUE": ("Resistance marker catalogue", "analyses.resistance_catalog",
                  "DA002 genotypes these markers; its status line decides whether real marker analysis may run.",
                  ("CATALOGUE",)),
    "FASTQ_PAIR": ("Paired FASTQ files", "samples.tsv fastq_1/fastq_2", "Raw reads for one sample and lane.",
                   ("FASTQ",)),
    "ALIGNMENT": ("Alignment (BAM/CRAM)", "samples.tsv bam", "Existing alignment used instead of FASTQ.",
                  ("BAM", "CRAM")),
    "UMI_READ_STRUCTURE": ("UMI read structure", "assay.read_structure",
                           "Tells fgbio where the UMI and template bases are in each read.", ()),
    "VCF": ("VCF input", "input_type", "VCF entry is not supported; use FASTQ, BAM or CRAM.", ("VCF",)),
}

GENERIC: dict[str, list[str]] = {
    "REFERENCE_MISMATCH": [
        "Use files made for the same reference assembly and release as the configured FASTA.",
        "Contig names must match exactly (for example Pf3D7_01_v3, not chr1 or Pf3D7_01).",
    ],
    "REFERENCE_MALFORMED": ["Check the file format against the reported detail and replace the file."],
    "BED_MALFORMED": [
        "Each line must be: contig<TAB>start<TAB>end, with whole-number coordinates.",
        "BED is zero-based and half-open: start counts from 0 and end is exclusive. Convert one-based "
        "coordinates by subtracting 1 from the start.",
        "Every interval must lie inside its contig, and the file must hold at least one interval.",
    ],
    "PRIMER_TABLE_MALFORMED": ["Correct the primer table at the reported line, using one-based inclusive coordinates."],
    "METADATA_MISSING_FIELD": ["Fill in or correct the reported field in the sample sheet."],
    "METADATA_DUPLICATE_SAMPLE": ["Each sample_id and lane pair may appear only once; remove or rename the duplicate."],
    "ALIGNMENT_UNSORTED": ["Sort by coordinate: samtools sort -o sorted.bam input.bam, then samtools index sorted.bam."],
    "ALIGNMENT_NO_INDEX": ["Create the index next to the file: samtools index FILE."],
    "ALIGNMENT_READ_GROUP": ["Every @RG needs ID, SM, LB and PL, and SM must equal the sample_id in the sample sheet."],
    "ALIGNMENT_PROCESSING_STATE": ["Supply an alignment processed as the assay requires (see docs/orchestration.md)."],
    "FASTQ_TRUNCATED": ["The file ends inside a read: copy or download it again and check its checksum."],
    "FASTQ_MATE_COUNT_MISMATCH": ["Read 1 and read 2 must come from the same sample and lane, with equal read counts."],
    "FASTQ_MATE_NAME_MISMATCH": ["Read 1 and read 2 must list the same read names in the same order."],
    "FASTQ_QUALITY_ENCODING": ["Qualities must be Phred+33 (Illumina 1.8 and later); re-export the reads in that encoding."],
    "FASTQ_MALFORMED_RECORD": ["Replace the file with a verified export; each read needs four lines (@name, bases, +, qualities)."],
    "UMI_READ_STRUCTURE": ["Set assay.read_structure from the library protocol, for example '8M+T 8M+T' "
                           "(8 UMI bases, then template)."],
    "VCF_INCOMPATIBLE": ["Use FASTQ, BAM or CRAM input instead."],
    "CONTENT_MISMATCH": [
        "PfVarX recognises files by their content, not their name or extension. Check that this is the file "
        "you meant, and replace it with one that contains what the input needs.",
        "To see what a file contains, look at its first lines (for a compressed file: zcat FILE | head).",
    ],
    "GFF_MALFORMED": [
        "Use the GFF3 from the same PlasmoDB release as the FASTA: nine tab-separated columns, one-based "
        "coordinates, with gene and CDS features.",
    ],
    "CATALOG_MALFORMED": [
        "The catalogue needs '#key=value' lines including '#status=' and the columns marker_id, gene_id, "
        "gene_name, chrom, codon, codon_start, strand, ref_aa and alt_aa (see the shipped catalogue in "
        "config/registry/).",
    ],
}

CORE_MASK_SOURCE = [
    "Use the standard core-genome definition (Miles et al. 2016, Genome Research), which MalariaGEN Pf3k, Pf7 "
    "and Pf8 use: the 3D7 region classification file (regions-20130225.onebased.txt), whose rows are typed Core, "
    "SubtelomericRepeat, SubtelomericHypervariable, InternalHypervariable or Centromere.",
    "Keep only the Core rows and convert to zero-based BED: awk -F'\\t' '$4==\"Core\"{print $1\"\\t\"$2-1\"\\t\"$3}' "
    "regions-20130225.onebased.txt > inputs/reference/core_genome.bed",
    "Do not use a gene or feature BED (for example one made from the GFF) as the core mask: calling would be "
    "limited to genes and the hypervariable var, rifin and stevor genes would stay in.",
]


def _raw(config: dict[str, Any], key: str) -> Any:
    value: Any = config
    for part in key.split("."):
        value = value.get(part) if isinstance(value, dict) else None
    return value


def _neighbours(root: Path, path: Path, kinds: tuple[str, ...]) -> list[str]:
    """Files in the same folder whose content is the right kind, whatever they are called."""
    from pfvarx.inputs.content import find_in_folder

    folder = (path if path.is_absolute() else root / path).parent
    return [p.name for p in find_in_folder(folder, kinds) if p.name != path.name][:10] if kinds else []


def explain(finding: dict[str, Any], kind: str, root: Path, raw_config: dict[str, Any]) -> dict[str, Any]:
    """Return ``file``, ``used_for``, ``why`` and ``actions`` for one failed or warned check."""
    if kind.startswith("CONTENT:"):  # a configured file whose content is the wrong kind
        setting = kind.split(":", 1)[1]
        kind = next((k for k, role in ROLES.items() if role[1] == setting), kind)
    name, key, purpose, kinds = ROLES.get(kind, ("Input", "", "", ()))
    reason, detail, path = finding.get("reason", ""), finding.get("detail", ""), finding.get("path", "")
    actions: list[str] = []
    if reason == "FILE_MISSING":
        why = f"{name} was not found at {path}."
        in_config = bool(key) and " " not in key and "/" not in key  # a config/config.yaml setting
        configured = _raw(raw_config, key) if in_config else None
        if isinstance(configured, str) and configured.endswith("/"):
            why += (f" {key} in config/config.yaml names only the folder {configured}, so PfVarX looks for the "
                    f"default file name {Path(path).name} there.")
            actions.append(f"Either save the file as {path}, or set {key} in config/config.yaml to the full path "
                           "of your file.")
        elif in_config:
            actions.append(f"Put the file at {path}, or correct {key} in config/config.yaml to point to it.")
        else:
            actions.append(f"Put the file at {path}, or correct its path in the sample sheet.")
        found = _neighbours(root, Path(path), kinds) if path and " | " not in path else []
        if len(found) == 1:
            actions.append(f"{found[0]} in that folder has the right kind of content. Use it only if it really is "
                           f"the {name}.")
        elif found:
            actions.append(f"Several files in that folder have the right kind of content: {', '.join(found)}. "
                           f"Set {key} to the full path of the one you mean.")
        if kind == "CORE_MASK":
            actions += CORE_MASK_SOURCE
        elif kind == "FASTA":
            actions.append("Download the P. falciparum 3D7 genome FASTA from PlasmoDB and record its release in "
                           "reference.release_label.")
        elif kind == "GFF":
            actions.append("Use the GFF from the same PlasmoDB release as the FASTA, or set reference.gff to null "
                           "if you do not need annotation or drug-resistance analysis.")
        elif kind == "SAMPLE_METADATA":
            actions.append("Register samples with ./pfvarx register-inputs (or the GUI Input files tab), or write "
                           "inputs/metadata/samples.tsv with sample_id, lane and library_id columns.")
    else:
        if reason == "CORE_MASK_SUSPECT":
            why = f"{name} may be the wrong file: {detail}."
        elif reason == "CONTENT_MISMATCH":
            why = f"{name} has the wrong content: {detail}."
        else:
            why = f"{name} failed the {reason.replace('_', ' ').lower()} check: {detail}."
        if finding.get("record"):
            why += f" (line {finding['record']})"
        actions += GENERIC.get(reason, ["Correct the input as described in the detail."])
        if kind == "CORE_MASK" and reason in ("BED_MALFORMED", "REFERENCE_MISMATCH"):
            actions += CORE_MASK_SOURCE[:2]
        if reason == "CORE_MASK_SUSPECT":
            actions = ["Check that this file is the core-genome mask and not a gene or feature list. With a gene "
                       "list, calling and QC would cover genes only and keep the hypervariable gene families.",
                       *CORE_MASK_SOURCE, "If this really is your study's own core mask, record why in the study "
                       "protocol; the warning does not stop the run."]
    if finding.get("level") == "FAIL":
        actions.append("Then run ./pfvarx validate again (GUI: Input Files Validation).")
    return {"file": name, "used_for": purpose, "why": why, "actions": actions}
