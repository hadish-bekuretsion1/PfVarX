"""Load, validate and freeze the three user YAML files plus maintainer parameters.

Precedence (fixed, highest first): scientific invariants in code, maintainer
registries and parameter files, the user profile (config/quality/compute YAML).
There are no command-line overrides of scientific settings.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ValidationError

from pfvarx.configuration.models import Bundle, ComputeConfig, ProjectConfig, QualityConfig
from pfvarx.constants import Design, ParameterSource


class ConfigError(ValueError):
    """Raised with a list of precise, file- and key-localised problems."""

    def __init__(self, problems: list[str]):
        super().__init__("; ".join(problems))
        self.problems = problems


USER_FILES = {"config": "config/config.yaml", "quality": "config/quality.yaml", "compute": "config/compute.yaml"}


def read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ConfigError([f"{path}: file not found (run 'pixi run config-init')"])
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError([f"{path}: YAML syntax error: {exc}"]) from exc
    if not isinstance(data, dict):
        raise ConfigError([f"{path}: top level must be a mapping"])
    return data


def _errors(name: str, exc: ValidationError) -> list[str]:
    out = []
    for err in exc.errors():
        loc = ".".join(str(p) for p in err["loc"]) or "<root>"
        out.append(f"{USER_FILES[name]}: {loc}: {err['msg']}")
    return out


LEGACY_IDS = {"DA001": "DA002", "DA002": "DA003", "DA003": "DA001"}


def migrate_project(project: dict[str, Any]) -> dict[str, Any]:
    if project.get("schema_version") != "pfvarx.config/2":
        return project
    analyses = dict(project.get("analyses") or {})
    analyses["requested"] = [LEGACY_IDS.get(key, key) for key in analyses.get("requested", [])]
    analyses["scales"] = {LEGACY_IDS.get(key, key): value for key, value in analyses.get("scales", {}).items()}
    source = analyses.get("spatial_outcome_source")
    if source:
        analyses["spatial_outcome_source"] = LEGACY_IDS.get(source, source)
    return {**project, "schema_version": "pfvarx.config/3", "analyses": analyses}


# Folder shortcuts (a setting that names only a folder) and the file name used when the
# folder holds no single file with the right content (validation then reports it).
FOLDER_SETTINGS = {
    "samples_tsv": "inputs/metadata/samples.tsv",
    "reference.fasta": "inputs/reference/Pf3D7.fasta",
    "reference.core_mask_bed": "inputs/reference/core_genome.bed",
    "reference.gff": "inputs/reference/Pf3D7.gff",
    "assay.library_catalog": "inputs/metadata/library_catalog.yaml",
    "analyses.resistance_catalog": "inputs/metadata/verified_catalog.tsv",
}


def resolve_input_defaults(project: dict[str, Any], root: Path | None = None) -> dict[str, Any]:
    """Expand a folder-only file setting to the file it means.

    The file is chosen by *content*, not by name: if the folder holds exactly one file whose
    content suits the setting (for example one BED for the core mask), that file is used.
    Otherwise the conventional file name is used, so validation can report the missing or
    ambiguous file with the candidates it found.
    """
    from pfvarx.inputs.content import SETTING_KINDS, find_in_folder

    resolved = {k: dict(v) if isinstance(v, dict) else v for k, v in project.items()}
    for setting, default in FOLDER_SETTINGS.items():
        section, _, key = setting.rpartition(".")
        group = resolved.get(section) if section else resolved
        if not isinstance(group, dict) or not isinstance(group.get(key), str):
            continue
        value = group[key]
        folder = Path(value)
        is_folder = value.endswith("/") or (root is not None and (root / folder).is_dir())
        if not is_folder:
            continue
        matches = find_in_folder(root / folder, SETTING_KINDS[setting]) if root is not None else []
        if len(matches) == 1:
            group[key] = str(matches[0].relative_to(root)) if root and matches[0].is_relative_to(root) \
                else str(matches[0])
        else:
            group[key] = str(folder / Path(default).name)
    return resolved


def validate_bundle(root: Path, parsed: dict[str, Any]) -> Bundle:
    problems: list[str] = []
    project = migrate_project(parsed["config"])
    if project.get("schema_version") == "pfvarx.config/3":
        raise ConfigError([
            f"{USER_FILES['config']}: uses configuration schema 3. Run './pfvarx config-migrate --apply' "
            "(backs up the old files) to move assay.design to assay_design and the panel, UMI and "
            "technique settings into the library catalog and sample manifest; see docs/configuration.md."
        ])
    parsed = {**parsed, "config": resolve_input_defaults(project, root)}
    models: dict[str, type[BaseModel]] = {"config": ProjectConfig, "quality": QualityConfig, "compute": ComputeConfig}
    built: dict[str, Any] = {}
    for name, model in models.items():
        try:
            built[name] = model.model_validate(parsed[name])
        except ValidationError as exc:
            problems += _errors(name, exc)
    if problems:
        raise ConfigError(problems)
    bundle = Bundle(**built)
    analyses = bundle.config.analyses
    if analyses.selection == "auto":
        registry = read_yaml(root / "config/registry/downstream_analyses.yaml")["analyses"]
        selected = sorted(da for da, contract in registry.items()
                          if bundle.config.assay_design.value in [Design.parse(d).value for d in contract["designs"]])
        if analyses.spatial_outcome_source not in (None, "DA002"):
            raise ConfigError(["Automatic selection supports only DA002 as the spatial outcome source"])
        changes = {
            "requested": selected,
            "population_column": analyses.population_column or "population",
            "spatial_outcome_source": "DA002",
        }
        bundle = bundle.model_copy(update={"config": bundle.config.model_copy(
            update={"analyses": analyses.model_copy(update=changes)}
        )})
    installation = root / "config/installation.json"
    if installation.is_file():
        from pfvarx.io import read_json

        metadata = read_json(installation)
        if metadata["profile"] not in ("development", "on-demand"):  # on demand: tools follow the config
            available = metadata["analyses_available"]
            if metadata["schema"] == "pfvarx.installation/1":
                available = [LEGACY_IDS.get(da, da) for da in available]
            from pfvarx.eligibility.engine import internal_prerequisites, load_registry

            effective = bundle.config.analyses
            prerequisites = internal_prerequisites(effective.requested, load_registry(root), effective.spatial_outcome_source)
            unavailable = set(effective.requested + prerequisites) - set(available)
            installed_input = metadata.get("input_type", "FASTQ")
            supported_inputs = {"FASTQ", "BAM", "CRAM"} if installed_input == "FASTQ" else {installed_input}
            reasons = []
            if bundle.config.final_output.value > metadata["max_final_output"]:
                reasons.append(f"final_output {bundle.config.final_output.value} needs tools this workspace did not "
                               f"install (installed for {metadata['max_final_output']}); bootstrap a new workspace with "
                               f"--final-output {bundle.config.final_output.value}")
            if unavailable:
                reasons.append(f"analyses {sorted(unavailable)} are not installed; bootstrap with --analyses for them")
            if bundle.config.assay_design != Design.parse(metadata["design"]):
                reasons.append(f"assay_design {bundle.config.assay_design.value} differs from the installed design "
                               f"{metadata['design']}; bootstrap with --design {bundle.config.assay_design.value}")
            if bundle.config.input_type.value not in supported_inputs:
                reasons.append(f"input_type {bundle.config.input_type.value} is not installed; bootstrap with "
                               f"--input-type {bundle.config.input_type.value}")
            if reasons:
                raise ConfigError(["Configuration exceeds this selected installation: " + "; ".join(reasons)
                                   + ". Do not delete installation metadata."])
    return bundle


def load_bundle(root: Path) -> Bundle:
    return validate_bundle(root, {name: read_yaml(root / rel) for name, rel in USER_FILES.items()})


def load_parameters(root: Path) -> dict[str, dict[str, Any]]:
    """Load maintainer parameter files and verify every entry records its source class."""
    params: dict[str, dict[str, Any]] = {}
    problems: list[str] = []
    allowed = {s.value for s in ParameterSource}
    for layer in range(7):
        path = root / "config" / "parameters" / f"layer{layer}.yaml"
        data = read_yaml(path)
        for key, entry in (data.get("parameters") or {}).items():
            if not isinstance(entry, dict) or "value" not in entry or "source" not in entry:
                problems.append(f"{path}: {key}: needs 'value' and 'source'")
                continue
            if entry["source"] not in allowed:
                problems.append(f"{path}: {key}: source '{entry['source']}' not allowed")
            if entry["source"] == ParameterSource.UNSET_BLOCKING and entry["value"] is not None:
                problems.append(f"{path}: {key}: UNSET_BLOCKING parameters must have value null")
            params[f"layer{layer}.{key}"] = entry
    if problems:
        raise ConfigError(problems)
    return params


def param(params: dict[str, dict[str, Any]], key: str) -> Any:
    """Return a parameter value, refusing UNSET_BLOCKING values."""
    entry = params[key]
    if entry["source"] == ParameterSource.UNSET_BLOCKING:
        raise ConfigError([f"parameter {key} is UNSET_BLOCKING; calibrate it before use"])
    return entry["value"]


def canonical_json(data: Any) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def effective_config(root: Path) -> dict[str, Any]:
    bundle = load_bundle(root)
    params = load_parameters(root)
    eff: dict[str, Any] = {
        "config": bundle.config.model_dump(mode="json"),
        "quality": bundle.quality.model_dump(mode="json"),
        "compute": bundle.compute.model_dump(mode="json"),
        "parameters": params,
    }
    eff["hash"] = sha256_text(canonical_json(eff))
    return eff
