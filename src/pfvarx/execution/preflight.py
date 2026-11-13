"""Config-driven preparation and auditable verification before scientific execution."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pfvarx.configuration.loader import effective_config, load_bundle
from pfvarx.inputs.organize import discover
from pfvarx.inputs.registry import validate_all
from pfvarx.initialization.layout import SAMPLE_COLUMNS
from pfvarx.io import write_json, write_tsv
from pfvarx.references.prepare import prepare, verify_prepared


VALIDATION_FIXES = {
    "FILE_MISSING": "Add the required file or correct its saved path, then validate again.",
    "REFERENCE_MISMATCH": "Use files that match the configured reference assembly and contig names.",
    "REFERENCE_MALFORMED": "Check the file format against the finding and replace malformed inputs.",
    "BED_MALFORMED": "Check tab-separated BED coordinates, contig names, and reference bounds.",
    "PRIMER_TABLE_MALFORMED": "Correct the primer TSV using the assay definitions and the reported field.",
    "METADATA_MISSING_FIELD": "Fill in the reported field in the sample sheet.",
    "METADATA_DUPLICATE_SAMPLE": "Review the sample sheet and correct duplicate sample/lane records.",
    "ALIGNMENT_UNSORTED": "Coordinate-sort the alignment and regenerate its index.",
    "ALIGNMENT_NO_INDEX": "Create the matching index and place it next to the alignment.",
    "ALIGNMENT_READ_GROUP": "Correct read group and sample identifiers to match the sample sheet.",
    "ALIGNMENT_PROCESSING_STATE": "Supply an alignment that meets the configured processing requirements.",
    "FASTQ_TRUNCATED": "Replace the incomplete FASTQ with a verified copy.",
    "FASTQ_MATE_COUNT_MISMATCH": "Supply both mates from the same sample and lane with matching record counts.",
    "FASTQ_MATE_NAME_MISMATCH": "Supply matched mates with corresponding read IDs in the same order.",
    "FASTQ_QUALITY_ENCODING": "Check the input quality encoding against the sequencing export and configuration.",
    "FASTQ_MALFORMED_RECORD": "Replace the malformed FASTQ with a verified export.",
    "UMI_READ_STRUCTURE": "Match the configured read structure to the documented library and read lengths.",
    "VCF_INCOMPATIBLE": "Supply a VCF compatible with the reference and sample configuration.",
}


def validate_inputs(root: Path) -> dict[str, Any]:
    bundle = load_bundle(root)
    cfg = bundle.config
    sheet = root / cfg.samples_tsv
    if cfg.input_type.value == "FASTQ" and not sheet.exists() and not sheet.is_symlink():
        discover(root, root / "inputs/fastq", sheet)
    report = validate_all(root, bundle)
    for entry in report["reports"]:
        for finding in entry["findings"]:
            finding["solution"] = VALIDATION_FIXES.get(
                finding["reason"], "Review the reported detail and correct the affected input or configuration."
            )
    write_json(root / "metadata/validation_report.json", report)
    matrix = []
    for entry in report["reports"]:
        for metric, value in entry.get("facts", {}).items():
            if isinstance(value, (str, int, float, bool)) or value is None:
                matrix.append(
                    {
                        "path": entry["path"],
                        "kind": entry["kind"],
                        "status": entry["status"],
                        "metric": metric,
                        "observed": value,
                        "source": "input validation; not post-alignment QC",
                    }
                )
    write_tsv(
        root / "metadata/input_quality_matrix.tsv", matrix, ["path", "kind", "status", "metric", "observed", "source"]
    )
    eff = effective_config(root)
    write_json(root / "metadata/configuration_summary.json", eff)
    return report


def prepare_inputs(root: Path, dry_run: bool = False) -> dict[str, Any]:
    bundle = load_bundle(root)
    if dry_run:
        return prepare(root, bundle, True)
    sheet = root / bundle.config.samples_tsv
    if bundle.config.input_type.value != "FASTQ" and not sheet.exists() and not sheet.is_symlink():
        sheet.parent.mkdir(parents=True, exist_ok=True)
        with sheet.open("x", encoding="utf-8") as handle:
            handle.write("\t".join(SAMPLE_COLUMNS) + "\n")
        raise ValueError(f"Created sample sheet template at {sheet}. Register or enter the samples, then rerun preparation.")
    validation = validate_inputs(root)
    if validation["status"] == "FAIL":
        failures = [
            f"{finding['path']}: {finding['detail']} Suggested fix: {finding['solution']}"
            for entry in validation["reports"]
            for finding in entry["findings"]
            if finding["level"] == "FAIL"
        ]
        details = "\n".join(failures)
        raise ValueError(f"Input validation failed. See metadata/validation_report.json.\n{details}")
    result = prepare(root, bundle, False)
    verification = verify_prepared(root, bundle)
    write_json(root / "metadata/preparation_report.json", {"preparation": result, "verification": verification})
    if verification["status"] != "PASS":
        raise ValueError("Prepared verification failed; inspect metadata/preparation_report.json")
    return {**result, "verification": verification}


def verify_inputs(root: Path) -> dict[str, Any]:
    validation = validate_inputs(root)
    verification = verify_prepared(root, load_bundle(root))
    result = {
        "status": "PASS"
        if validation["status"] in ("PASS", "PARTIAL") and verification["status"] == "PASS"
        else "FAIL",
        "input_status": validation["status"],
        "sample_status": validation["samples"],
        "reference": verification,
    }
    write_json(root / "metadata/prepared_verification.json", result)
    return result


def orchestrate(root: Path, through: str = "execute") -> dict[str, Any]:
    from pfvarx.execution.retention import workspace_lock

    with workspace_lock(root):
        return _orchestrate(root, through)


def _orchestrate(root: Path, through: str = "execute") -> dict[str, Any]:
    from pfvarx.execution.runner import _execute, dry_run
    from pfvarx.initialization.workspace import initialize

    allowed = ["initialize", "validate", "prepare", "verify", "execute"]
    if through not in allowed:
        raise ValueError("through must be " + ", ".join(allowed))
    stages = [{"stage": "initialize", "result": initialize(root)}]
    if through == "initialize":
        return {"status": "READY_FOR_CONFIGURATION", "stages": stages}
    report = validate_inputs(root)
    stages.append({"stage": "validate", "status": report["status"]})
    if report["status"] == "FAIL":
        raise ValueError("Validation failed; see metadata/validation_report.json")
    if through == "validate":
        return {"status": "VALIDATED", "stages": stages}
    stages.append({"stage": "prepare", "result": prepare_inputs(root)})
    if through == "prepare":
        return {"status": "PREPARED", "stages": stages}
    verification = verify_inputs(root)
    stages.append({"stage": "verify", "result": verification})
    if verification["status"] != "PASS":
        raise ValueError("Prepared verification failed; see metadata/prepared_verification.json")
    if through == "verify":
        return {"status": "VERIFIED", "stages": stages}
    preflight = dry_run(root)
    if preflight["status"] != "OK":
        raise ValueError("Workflow preflight failed: " + str(preflight))
    return {**_execute(root), "orchestration_stages": stages}
