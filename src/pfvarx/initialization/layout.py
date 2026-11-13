"""Create purpose- and assay-specific input folders; never create result folders.

``config_init`` copies the versioned templates for the chosen purpose without
overwriting an existing file. ``init`` creates only input directories and
metadata templates. ``check_layout`` is read-only.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from pfvarx.configuration.models import ProjectConfig
from pfvarx.constants import Design, InputType, Purpose
from pfvarx.io import atomic_write_text, within

SAMPLE_COLUMNS = ["sample_id", "lane", "fastq_1", "fastq_2", "bam", "vcf", "library_id",
                  "flowcell_type", "population", "country", "latitude", "longitude",
                  "collection_date"]


def config_init(root: Path, purpose: str) -> list[Path]:
    try:
        chosen = Purpose(purpose)
    except ValueError as exc:
        raise ValueError("--purpose must be 'analysis' or 'benchmarking'") from exc
    src = root / "config" / "templates"
    sources = {"config.yaml": src / f"{chosen.value}.yaml",
               "quality.yaml": src / "quality.yaml", "compute.yaml": src / "compute.yaml"}
    created: list[Path] = []
    for name in ("config.yaml", "quality.yaml", "compute.yaml"):
        dst = root / "config" / name
        if dst.exists():
            continue
        shutil.copyfile(sources[name], dst)
        created.append(dst)
    return created


def required_dirs(cfg: ProjectConfig) -> list[str]:
    dirs = ["inputs/reference", "inputs/metadata"]
    if cfg.input_type == InputType.FASTQ:
        dirs.append("inputs/fastq")
    elif cfg.input_type in (InputType.BAM, InputType.CRAM):
        dirs.append("inputs/alignments")
    else:
        dirs.append("inputs/vcf")
    if cfg.assay.design != Design.WGS:
        dirs.append("inputs/targets")
    if cfg.purpose == Purpose.BENCHMARKING:
        dirs += ["inputs/truth", "inputs/partitions"]
    return dirs


def create_input_dirs(root: Path, directories: list[str]) -> list[str]:
    paths = [(rel, within(root, rel)) for rel in directories]
    for rel, path in paths:
        if any(part.exists() and not part.is_dir() for part in (path, *path.parents)):
            raise ValueError(f"Cannot create input folder: {rel}")
    created = []
    for rel, path in paths:
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
            created.append(rel)
    return created


def init_inputs(root: Path, cfg: ProjectConfig) -> dict[str, Any]:
    directories = required_dirs(cfg)
    created = create_input_dirs(root, directories)
    samples = root / cfg.samples_tsv
    template_written = False
    if not samples.exists():
        atomic_write_text(samples, "\t".join(SAMPLE_COLUMNS) + "\n")
        template_written = True
    return {"directories": [rel + "/" for rel in directories], "created_dirs": created, "samples_template_written": template_written,
            "result_dirs_created": []}


def check_layout(root: Path, cfg: ProjectConfig) -> dict[str, Any]:
    missing = [rel for rel in required_dirs(cfg) if not (root / rel).is_dir()]
    files: dict[str, Path | None] = {"reference.fasta": cfg.reference.fasta, "reference.core_mask_bed": cfg.reference.core_mask_bed,
             "samples_tsv": cfg.samples_tsv}
    if cfg.final_output.value == "FO004" and (cfg.reference.gff or cfg.analyses.selection != "auto"):
        files["reference.gff"] = cfg.reference.gff
    if cfg.assay.targets_bed:
        files["assay.targets_bed"] = cfg.assay.targets_bed
    if cfg.assay.primers_tsv:
        files["assay.primers_tsv"] = cfg.assay.primers_tsv
    missing_files = {k: str(v) for k, v in files.items() if v is None or not (root / v).is_file()}
    misplaced = []
    for unexpected in ("results", "work"):
        if (root / unexpected).exists() and not any((root / "runs").glob("*")):
            misplaced.append(f"{unexpected}/ exists before any run")
    status = "PASS" if not (missing or missing_files or misplaced) else "FAIL"
    return {"status": status, "missing_dirs": missing, "missing_files": missing_files,
            "misplaced": misplaced}
