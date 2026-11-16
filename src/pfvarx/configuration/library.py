"""Library catalog: targeted panels and library protocols, referenced from the sample manifest.

The run settings name only the assay design (``assay_design``). What a library *is* lives
here and in the sample manifest:

* a **panel** (``panel_id`` + ``panel_version``) defines the targeted intervals and the
  primer (amplicon) or probe (capture) definitions;
* a **library protocol** (``library_protocol_id``) records how the library was made: an
  optional technique label and amplification route, and for UMI libraries the read
  structure used for UMI extraction plus any consensus settings that differ from the
  maintainer defaults in ``config/parameters/layer1.yaml``.

The catalog is a YAML file (``assay.library_catalog``, conventionally
``inputs/metadata/library_catalog.yaml``). It is required for targeted designs and optional
for ``wgs``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

from pfvarx.constants import Route, Technique

CATALOG_SCHEMA = "pfvarx.library_catalog/1"
ID_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$"

# Consensus settings a protocol may set, and the layer-1 parameter each replaces.
UMI_SETTINGS = {
    "group_strategy": "fgbio_group_strategy",
    "group_edits": "fgbio_group_edits",
    "group_min_map_q": "fgbio_group_min_map_q",
    "consensus_min_reads": "fgbio_consensus_min_reads",
    "consensus_error_pre_umi": "fgbio_consensus_error_pre_umi",
    "consensus_error_post_umi": "fgbio_consensus_error_post_umi",
    "consensus_min_input_base_quality": "fgbio_consensus_min_input_base_quality",
    "filter_min_reads": "fgbio_filter_min_reads",
    "filter_min_base_quality": "fgbio_filter_min_base_quality",
    "filter_max_read_error_rate": "fgbio_filter_max_read_error_rate",
    "filter_max_base_error_rate": "fgbio_filter_max_base_error_rate",
    "filter_max_no_call_fraction": "fgbio_filter_max_no_call_fraction",
}


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class Panel(_Strict):
    panel_id: str = Field(pattern=ID_PATTERN)
    panel_version: str = Field(pattern=ID_PATTERN)
    targets_bed: Path = Field(description="Targeted intervals (BED, 0-based half-open)")
    primers_tsv: Path | None = Field(default=None, description="Amplicon primers (fgbio TrimPrimers table)")
    probes_bed: Path | None = Field(default=None, description="Capture probes (BED); recorded, not trimmed")
    description: str | None = None

    @model_validator(mode="after")
    def _definitions(self) -> Panel:
        if self.primers_tsv is None and self.probes_bed is None:
            raise ValueError("a panel needs primers_tsv (amplicon) or probes_bed (capture)")
        return self

    @property
    def key(self) -> tuple[str, str]:
        return self.panel_id, self.panel_version


class UmiSettings(_Strict):
    read_structure: str = Field(description="fgbio read structures for R1 and R2, e.g. '8M+T 8M+T'")
    group_strategy: Literal["Identity", "Edit", "Adjacency", "Paired"] | None = None
    group_edits: int | None = Field(default=None, ge=0)
    group_min_map_q: int | None = Field(default=None, ge=0)
    consensus_min_reads: int | None = Field(default=None, ge=1)
    consensus_error_pre_umi: int | None = Field(default=None, ge=1)
    consensus_error_post_umi: int | None = Field(default=None, ge=1)
    consensus_min_input_base_quality: int | None = Field(default=None, ge=0)
    filter_min_reads: int | None = Field(default=None, ge=1)
    filter_min_base_quality: int | None = Field(default=None, ge=0)
    filter_max_read_error_rate: float | None = Field(default=None, ge=0, le=1)
    filter_max_base_error_rate: float | None = Field(default=None, ge=0, le=1)
    filter_max_no_call_fraction: float | None = Field(default=None, ge=0, le=1)


class LibraryProtocol(_Strict):
    library_protocol_id: str = Field(pattern=ID_PATTERN)
    technique: Technique | None = None
    amplification_route: Route | None = None
    umi: UmiSettings | None = None
    description: str | None = None


class LibraryCatalog(_Strict):
    schema_version: str
    panels: list[Panel] = Field(default_factory=list)
    library_protocols: list[LibraryProtocol] = Field(default_factory=list)

    @field_validator("schema_version")
    @classmethod
    def _schema(cls, value: str) -> str:
        if value != CATALOG_SCHEMA:
            raise ValueError(f"library catalog schema_version must be {CATALOG_SCHEMA}")
        return value

    @model_validator(mode="after")
    def _unique(self) -> LibraryCatalog:
        keys = [p.key for p in self.panels]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate panel_id + panel_version")
        ids = [p.library_protocol_id for p in self.library_protocols]
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate library_protocol_id")
        return self

    def panel(self, panel_id: str, version: str) -> Panel | None:
        return next((p for p in self.panels if p.key == (panel_id, version)), None)

    def protocol(self, protocol_id: str) -> LibraryProtocol | None:
        return next((p for p in self.library_protocols if p.library_protocol_id == protocol_id), None)


class CatalogError(ValueError):
    pass


def load_catalog(path: Path) -> LibraryCatalog:
    """Parse and validate a library catalog; problems are raised as one CatalogError."""
    if not path.is_file():
        raise CatalogError(f"{path}: library catalog not found")
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise CatalogError(f"{path}: YAML syntax error: {exc}") from exc
    if not isinstance(data, dict):
        raise CatalogError(f"{path}: top level must be a mapping")
    try:
        return LibraryCatalog.model_validate(data)
    except ValidationError as exc:
        problems = "; ".join(f"{'.'.join(str(p) for p in e['loc']) or '<root>'}: {e['msg']}" for e in exc.errors())
        raise CatalogError(f"{path}: {problems}") from exc


def umi_settings(protocol: LibraryProtocol | None, params: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Resolved UMI extraction and consensus settings, each with where its value came from."""
    out: dict[str, Any] = {}
    for name, parameter in UMI_SETTINGS.items():
        given = getattr(protocol.umi, name) if protocol and protocol.umi else None
        out[name] = {"value": params[f"layer1.{parameter}"]["value"] if given is None else given,
                     "source": "PARAMETER_DEFAULT" if given is None else "LIBRARY_PROTOCOL"}
    out["read_structure"] = {"value": protocol.umi.read_structure if protocol and protocol.umi else None,
                             "source": "LIBRARY_PROTOCOL"}
    return out
