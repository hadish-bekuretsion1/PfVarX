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
    Route,
    Scale,
)

CONFIG_SCHEMA = "pfvarx.config/3"
QUALITY_SCHEMA = "pfvarx.quality/2"
COMPUTE_SCHEMA = "pfvarx.compute/2"


class Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class AssayConfig(Strict):
    design: Design
    sequencing_type: Literal["WGS", "SWGS", "AMPLICON", "TBS_UMI", "TBS_NON_UMI"] | None = None
    amplification_route: Route
    read_structure: str | None = Field(
        default=None, description="fgbio read structure, e.g. '8M143T 8M143T'; UMI only"
    )
    targets_bed: Path | None = None
    primers_tsv: Path | None = None
    two_colour_chemistry: bool = False
    platform: str = "ILLUMINA"

    @model_validator(mode="after")
    def _design_consistency(self) -> AssayConfig:
        compatible = {
            "WGS": (Design.WGS,), "SWGS": (Design.WGS,),
            "AMPLICON": (Design.TARGETED_NON_UMI, Design.TARGETED_UMI),
            "TBS_UMI": (Design.TARGETED_UMI,),
            "TBS_NON_UMI": (Design.TARGETED_NON_UMI,),
        }
        if self.sequencing_type and self.design not in compatible[self.sequencing_type]:
            raise ValueError("sequencing_type is incompatible with assay.design")
        targeted = self.design in (Design.TARGETED_NON_UMI, Design.TARGETED_UMI)
        if targeted and (self.targets_bed is None or self.primers_tsv is None):
            raise ValueError("targeted designs require assay.targets_bed and assay.primers_tsv")
        if targeted and self.amplification_route != Route.AMPLICON:
            raise ValueError("targeted designs require amplification_route: AMPLICON")
        if self.design == Design.WGS and self.amplification_route == Route.AMPLICON:
            raise ValueError("WGS cannot use amplification_route: AMPLICON")
        if self.design == Design.TARGETED_UMI and not self.read_structure:
            raise ValueError("TARGETED_UMI requires assay.read_structure")
        if self.design != Design.TARGETED_UMI and self.read_structure:
            raise ValueError("assay.read_structure is only valid for TARGETED_UMI")
        return self


class ReferenceConfig(Strict):
    fasta: Path
    gff: Path | None = None
    core_mask_bed: Path
    host_fasta: Path | None = None
    release_label: str = Field(description="PlasmoDB release, e.g. '68'; recorded in provenance")
    snpeff_genome_id: str | None = None


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
    intermediates: str = Field(default="successful", pattern=r"^(successful|keep)$")


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
    assay: AssayConfig
    reference: ReferenceConfig
    samples_tsv: Path
    analyses: AnalysesConfig = Field(default_factory=AnalysesConfig)
    run_policy: str = Field(default="CONTINUE_INDEPENDENT", pattern=r"^(CONTINUE_INDEPENDENT|STRICT)$")
    random_seed: int = Field(default=20261106, ge=0)
    caller: CallerConfig = Field(default_factory=CallerConfig)
    retention: RetentionConfig = Field(default_factory=RetentionConfig)

    @field_validator("schema_version")
    @classmethod
    def _schema(cls, value: str) -> str:
        if value != CONFIG_SCHEMA:
            raise ValueError(f"config schema_version must be {CONFIG_SCHEMA}")
        return value

    @model_validator(mode="after")
    def _mode(self) -> ProjectConfig:
        if self.final_output == FinalOutput.FO004 and not self.analyses.requested and self.analyses.selection != "auto":
            raise ValueError("FO004 requires at least one requested downstream analysis")
        if self.final_output == FinalOutput.FO004 and self.reference.gff is None and self.analyses.selection != "auto":
            raise ValueError("FO004 requires reference.gff for reference-matched annotation")
        if self.final_output != FinalOutput.FO004 and (self.analyses.requested or self.analyses.selection == "auto"):
            raise ValueError("downstream analyses can only be requested with FO004")
        if self.input_type in (InputType.BAM, InputType.CRAM) and self.final_output == FinalOutput.FO001:
            raise ValueError("BAM/CRAM enter at Layer 2. Choose FO002, FO003 or FO004.")
        if self.input_type == InputType.VCF:
            raise ValueError("VCF entry is not implemented; supply FASTQ, BAM or CRAM")
        if self.assay.design != Design.WGS and set(self.analyses.requested) & {"DA005", "DA006"}:
            raise ValueError("DA005/DA006 require WGS; remove them or choose the actual WGS assay")
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
        description="Deployment target. HPC uses one allocated node; server and cloud use the local executor on the host."
    )
    max_cores: int = Field(ge=1)
    max_mem_mb: int = Field(ge=512)
    max_jobs: int = Field(ge=1)
    tmp_disk_mb: int = Field(ge=256)
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
