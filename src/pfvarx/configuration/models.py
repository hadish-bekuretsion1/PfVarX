"""Typed schemas for the three user-edited YAML files.

User YAML selects the scientific question and the computational envelope. It cannot
switch off a scientific gate: unknown keys are rejected (``extra="forbid"``) and the
scientific policy is loaded from maintainer-owned registries, not from user files.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from pfvarx.constants import (
    DOWNSTREAM_IDS,
    Design,
    FinalOutput,
    InputType,
    Purpose,
    Scale,
)

CONFIG_SCHEMA = "pfvarx.config/4"
QUALITY_SCHEMA = "pfvarx.quality/2"
COMPUTE_SCHEMA = "pfvarx.compute/2"


class Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class AssayConfig(Strict):
    """Library-wide settings. What each library is (panel, protocol, technique, route,
    imported-BAM state) is recorded per sample in the manifest and the library catalog."""

    library_catalog: Path | None = Field(
        default=None, description="Panels and library protocols (YAML); required for targeted designs"
    )
    two_colour_chemistry: bool = False
    platform: str = "ILLUMINA"


class ReferenceConfig(Strict):
    fasta: Path
    gff: Path | None = None
    core_mask_bed: Path
    host_removal: bool = Field(
        default=False,
        description="Human DNA removal (wgs only): true aligns to the parasite and host genomes together and removes host reads",
    )
    host_fasta: Path | None = Field(
        default=None,
        description="Host (human) genome FASTA used when host_removal is true; its bwa-mem2 index is supplied, not built",
    )
    release_label: str = Field(description="PlasmoDB release, e.g. '68'; recorded in provenance")
    snpeff_genome_id: str | None = Field(
        default=None, pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]*$",
        description="SnpEff database ID: letters, digits, dots, underscores and hyphens only",
    )

    @model_validator(mode="after")
    def _host(self) -> ReferenceConfig:
        if self.host_removal and self.host_fasta is None:
            raise ValueError("reference.host_removal is true: set reference.host_fasta to the human genome FASTA "
                             "(or set host_removal to false)")
        return self

    @property
    def host_fasta_used(self) -> Path | None:
        """The host FASTA when host-read removal is switched on, else None."""
        return self.host_fasta if self.host_removal else None


class AnalysesConfig(Strict):
    selection: Literal["explicit", "auto"] = "explicit"
    requested: list[str] = Field(default_factory=list)
    scales: dict[str, list[Scale]] = Field(default_factory=dict)
    population_column: str | None = None
    resistance_catalog: Path | None = Field(
        default=None, description="Versioned marker catalogue; default config/registry/resistance_catalog_pfvarx_v1.tsv"
    )
    spatial_outcome_source: str | None = Field(
        default=None, description="Source DA whose outcome DA007 summarises, e.g. DA002"
    )

    @field_validator("requested")
    @classmethod
    def _known(cls, value: list[str]) -> list[str]:
        unknown = sorted(set(value) - set(DOWNSTREAM_IDS))
        if unknown:
            raise ValueError(f"unknown downstream analyses: {unknown}")
        if len(set(value)) != len(value):
            raise ValueError("duplicate downstream analyses requested")
        return value

    @model_validator(mode="after")
    def _scales(self) -> AnalysesConfig:
        if self.selection == "auto" and (self.requested or self.scales):
            raise ValueError("automatic selection requires empty requested and scales")
        for key in self.scales:
            if key not in self.requested:
                raise ValueError(f"scale given for analysis not requested: {key}")
        if "DA007" in self.requested and self.spatial_outcome_source != "DA002":
            raise ValueError("DA007 requires analyses.spatial_outcome_source: DA002")
        return self


class RetentionConfig(Strict):
    intermediates: Literal["successful", "keep"] = "successful"


class CallerConfig(Strict):
    backend: Literal["native", "external"] = "native"
    executable: Path | None = None
    sha256: str | None = Field(default=None, pattern=r"^[a-f0-9]{64}$")

    @model_validator(mode="after")
    def _backend(self) -> CallerConfig:
        if self.backend == "external":
            if self.executable is None or not self.executable.is_absolute() or self.sha256 is None:
                raise ValueError("external caller requires an absolute executable path and SHA-256")
        elif self.executable is not None or self.sha256 is not None:
            raise ValueError("executable and sha256 are for the external caller only")
        return self


class ProjectConfig(Strict):
    schema_version: str
    purpose: Purpose
    project_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")
    input_type: InputType
    final_output: FinalOutput
    assay_design: Design = Field(description="Selects the processing branch: wgs | targeted_umi | targeted_no_umi")
    assay: AssayConfig = Field(default_factory=AssayConfig)
    reference: ReferenceConfig
    samples_tsv: Path
    study_metadata: bool = Field(
        default=False,
        description="Study Metadata: true merges a study table from study_metadata_dir into the sample metadata",
    )
    study_metadata_dir: Path = Field(
        default=Path("inputs/metadata"),
        description="Folder holding the study metadata table, used when study_metadata is true",
    )
    analyses: AnalysesConfig = Field(default_factory=AnalysesConfig)
    run_policy: Literal["CONTINUE_INDEPENDENT", "STRICT"] = "CONTINUE_INDEPENDENT"
    random_seed: int = Field(default=20261106, ge=0)
    caller: CallerConfig = Field(default_factory=CallerConfig)
    retention: RetentionConfig = Field(default_factory=RetentionConfig)

    @field_validator("schema_version")
    @classmethod
    def _schema(cls, value: str) -> str:
        if value != CONFIG_SCHEMA:
            raise ValueError(f"config schema_version must be {CONFIG_SCHEMA}")
        return value

    @field_validator("assay_design", mode="before")
    @classmethod
    def _design(cls, value: object) -> Design:
        try:
            return Design.parse(value)
        except ValueError:
            raise ValueError("assay_design must be wgs, targeted_umi or targeted_no_umi") from None

    @field_validator("input_type", mode="before")
    @classmethod
    def _input(cls, value: object) -> object:
        return value.upper() if isinstance(value, str) else value  # fastq, bam and cram are accepted

    @model_validator(mode="after")
    def _mode(self) -> ProjectConfig:
        if self.final_output == FinalOutput.FO004 and self.analyses.selection == "explicit" and not self.analyses.requested:
            raise ValueError("Select at least one downstream analysis for explicit FO004")
        if self.final_output == FinalOutput.FO004 and self.reference.gff is None and self.analyses.selection != "auto":
            raise ValueError("FO004 requires reference.gff for reference-matched annotation")
        if self.final_output != FinalOutput.FO004 and (self.analyses.requested or self.analyses.selection == "auto"):
            raise ValueError("downstream analyses can only be requested with FO004")
        if self.input_type in (InputType.BAM, InputType.CRAM) and self.final_output == FinalOutput.FO001:
            raise ValueError("BAM/CRAM enter at Layer 2. Choose FO002, FO003 or FO004.")
        if self.input_type == InputType.VCF:
            raise ValueError("VCF entry is not implemented; supply FASTQ, BAM or CRAM")
        if self.assay_design.targeted and self.assay.library_catalog is None:
            raise ValueError("targeted designs require assay.library_catalog (panels and library protocols)")
        if self.reference.host_removal and self.assay_design != Design.WGS:
            raise ValueError(f"Human DNA removal (reference.host_removal: true) is eligible only for assay_design "
                             f"wgs, not {self.assay_design.value}; set reference.host_removal to false")
        return self


class QualityConfig(Strict):
    schema_version: str
    policy_id: str = Field(description="Maintainer-registered QC policy identifier")
    report_not_applicable: bool = True

    @field_validator("schema_version")
    @classmethod
    def _schema(cls, value: str) -> str:
        if value != QUALITY_SCHEMA:
            raise ValueError(f"quality schema_version must be {QUALITY_SCHEMA}")
        return value


class JobLimits(Strict):
    max_threads: int = Field(ge=1)
    max_mem_mb: int = Field(ge=256)
    max_runtime_min: int = Field(ge=1)


class ComputeConfig(Strict):
    schema_version: str
    executor: Literal["local", "hpc", "server", "cloud", "slurm"] = Field(
        description=("Deployment target, recorded with the run. Every target runs the same single-host "
                     "executor: HPC uses one allocated node; server and cloud use the host.")
    )
    max_cores: int = Field(ge=1)
    max_mem_mb: int = Field(ge=512)
    max_jobs: int = Field(ge=1)
    tmp_disk_mb: int = Field(ge=256)
    tmp_dir: str | None = Field(
        default=None,
        description="Absolute directory for job temporary files (e.g. node-local scratch); null uses the system default.",
    )
    per_job: JobLimits
    max_retries: int = Field(default=2, ge=0, le=2)
    retry_memory_factor: float = Field(default=1.5, ge=1.0, le=2.0)
    log_level: str = Field(default="INFO", pattern=r"^(DEBUG|INFO|WARNING|ERROR)$")
    keep_attempt_logs: bool = True
    slurm_partition: str | None = None

    @field_validator("schema_version")
    @classmethod
    def _schema(cls, value: str) -> str:
        if value != COMPUTE_SCHEMA:
            raise ValueError(f"compute schema_version must be {COMPUTE_SCHEMA}")
        return value

    @model_validator(mode="after")
    def _envelope(self) -> ComputeConfig:
        if self.per_job.max_threads > self.max_cores:
            raise ValueError("per_job.max_threads exceeds max_cores")
        if self.per_job.max_mem_mb > self.max_mem_mb:
            raise ValueError("per_job.max_mem_mb exceeds max_mem_mb")
        if self.tmp_dir is not None and not self.tmp_dir.startswith("/"):
            raise ValueError("tmp_dir must be an absolute path, or null for the system temporary directory")
        if self.executor == "slurm" and not self.slurm_partition:
            raise ValueError("executor slurm requires slurm_partition")
        if self.executor == "slurm":
            raise ValueError(
                "executor slurm is specified but not yet tested in this release (UNSUPPORTED)"
            )
        return self


class Bundle(Strict):
    config: ProjectConfig
    quality: QualityConfig
    compute: ComputeConfig
