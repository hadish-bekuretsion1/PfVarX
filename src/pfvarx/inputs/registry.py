"""Validate every configured input and build the sample/lane/read-group registry."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pfvarx.configuration.loader import load_parameters, param
from pfvarx.configuration.models import Bundle
from pfvarx.constants import Design, InputType
from pfvarx.inputs import validate as v
from pfvarx.inputs.content import SETTING_KINDS, content_problem


def required_metadata(bundle: Bundle) -> set[str]:
    req: set[str] = set()
    an = bundle.config.analyses
    if an.selection == "auto":
        return req  # optional downstream metadata is assessed by eligibility
    if an.population_column and any(d in an.requested for d in ("DA003", "DA005")):
        req.add(an.population_column)
    if "DA007" in an.requested:
        req |= {"latitude", "longitude", "collection_date", "country"}
    return req


def validate_all(root: Path, bundle: Bundle) -> dict[str, Any]:
    cfg = bundle.config
    params = load_parameters(root)
    reports: list[v.FileReport] = []
    # Each configured file must contain what its setting needs; the file name is not checked.
    for setting, path in (
        ("reference.fasta", cfg.reference.fasta), ("reference.core_mask_bed", cfg.reference.core_mask_bed),
        ("reference.gff", cfg.reference.gff), ("assay.targets_bed", cfg.assay.targets_bed),
        ("assay.primers_tsv", cfg.assay.primers_tsv), ("samples_tsv", cfg.samples_tsv),
        ("analyses.resistance_catalog", cfg.analyses.resistance_catalog),
    ):
        if path is not None and (root / path).is_file():
            problem = content_problem(root / path, SETTING_KINDS[setting])
            if problem:
                report = v.FileReport(path=str(path), kind="CONTENT:" + setting)
                report.fail(v.Reason.CONTENT_MISMATCH, problem)
                reports.append(report)
    fasta = v.validate_fasta(root / cfg.reference.fasta)
    reports.append(fasta)
    contigs = fasta.facts.get("contigs", {})
    if contigs:
        reports.append(v.validate_bed(root / cfg.reference.core_mask_bed, contigs, "CORE_MASK"))
        if cfg.assay.targets_bed:
            reports.append(v.validate_bed(root / cfg.assay.targets_bed, contigs, "TARGETS"))
        if cfg.assay.primers_tsv:
            reports.append(v.validate_primers(root / cfg.assay.primers_tsv, contigs))
    if cfg.reference.gff:
        if contigs:  # structure and contig checks need the reference contigs
            reports.append(v.validate_gff(root / cfg.reference.gff, contigs))
        elif not (root / cfg.reference.gff).is_file():
            gff = v.FileReport(path=str(cfg.reference.gff), kind="GFF")
            gff.fail(v.Reason.FILE_MISSING, "GFF not found")
            reports.append(gff)
    if cfg.analyses.resistance_catalog:
        reports.append(v.validate_catalog(root / cfg.analyses.resistance_catalog))
    meta, rows = v.validate_samples(root / cfg.samples_tsv, root, cfg.input_type,
                                    cfg.assay.design, required_metadata(bundle))
    reports.append(meta)
    sample_status: dict[str, str] = {}
    metadata_global_failure = False
    for finding in meta.findings:
        if finding.level != "FAIL":
            continue
        if finding.record is None:
            metadata_global_failure = True
        else:
            sid = rows[finding.record - 2].get("sample_id")
            if sid:
                sample_status[sid] = "FAIL"
    if not metadata_global_failure:
        nrec = int(param(params, "layer0.fastq_records_checked_for_encoding"))
        for row in rows:
            sid = row.get("sample_id")
            if not sid or sample_status.get(sid) == "FAIL":
                continue
            if cfg.input_type == InputType.FASTQ:
                r1, r2 = root / row["fastq_1"], root / row["fastq_2"]
                rep = v.validate_fastq_pair(r1, r2, nrec)
                rep.facts["sample_id"] = sid
                reports.append(rep)
                if rep.status != "FAIL" and cfg.assay.design == Design.TARGETED_UMI and cfg.assay.read_structure:
                    lengths = (v.first_read_length(r1), v.first_read_length(r2))
                    us = v.validate_read_structure(cfg.assay.read_structure, lengths)
                    us.facts["sample_id"] = sid
                    reports.append(us)
                    if us.status == "FAIL":
                        rep.status = "FAIL"
            elif cfg.input_type in (InputType.BAM, InputType.CRAM) and contigs:
                rep = v.validate_alignment(root / row["bam"], contigs, root / cfg.reference.fasta)
                rep.facts["sample_id"] = sid
                if rep.facts.get("samples") and rep.facts["samples"] != [sid]:
                    rep.fail(v.Reason.ALIGNMENT_READ_GROUP,
                             f"@RG SM {rep.facts['samples']} does not match sample_id {sid}")
                reports.append(rep)
            else:
                rep = v.FileReport(path=row.get("vcf", ""), kind="VCF")
                reports.append(rep)
            if cfg.input_type in (InputType.BAM, InputType.CRAM) and (root / row["bam"]).is_file():
                # FASTQ content is checked inside the pair validation; alignments are checked here.
                problem = content_problem(root / row["bam"], (cfg.input_type.value,))
                if problem:
                    rep.fail(v.Reason.CONTENT_MISMATCH, problem)
            prev = sample_status.get(sid, "PASS")
            sample_status[sid] = "FAIL" if "FAIL" in (prev, rep.status) else (
                "WARNING" if "WARNING" in (prev, rep.status) else "PASS")
    global_fail = metadata_global_failure or any(
        r.status == "FAIL" for r in reports
        if r.kind not in ("SAMPLE_METADATA", "FASTQ_PAIR", "ALIGNMENT", "UMI_READ_STRUCTURE")
    )
    admitted = sum(status != "FAIL" for status in sample_status.values())
    if global_fail or not admitted:
        status = "FAIL"
    else:
        status = "PASS" if meta.status != "FAIL" and admitted == len(sample_status) else "PARTIAL"
    return {
        "schema": "pfvarx.validation/1",
        "status": status,
        "reference_contigs": contigs,
        "samples": sample_status,
        "registry": [row for row in rows if sample_status.get(row.get("sample_id", "")) in ("PASS", "WARNING")],
        "reports": [r.as_dict() for r in reports],
    }
