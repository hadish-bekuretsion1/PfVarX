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
from pfvarx.constants import ParameterSource


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


def validate_bundle(root: Path, parsed: dict[str, Any]) -> Bundle:
    problems: list[str] = []
    parsed = {**parsed, "config": migrate_project(parsed["config"])}
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
        selected = sorted(da for da, contract in registry.items() if bundle.config.assay.design.value in contract["designs"])
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
    if bundle.config.final_output.value == "FO004" and analyses.selection == "explicit":
        requested = ["DA001", *(da for da in analyses.requested if da != "DA001")]
        bundle = bundle.model_copy(update={"config": bundle.config.model_copy(
            update={"analyses": analyses.model_copy(update={"requested": requested})}
        )})
    installation = root / "config/installation.json"
    if installation.is_file():
        from pfvarx.io import read_json

        metadata = read_json(installation)
        if metadata["profile"] != "development":
            available = metadata["analyses_available"]
            if metadata["schema"] == "pfvarx.installation/1":
                available = [LEGACY_IDS.get(da, da) for da in available]
            unavailable = set(bundle.config.analyses.requested) - set(available)
            if (
                unavailable
                or bundle.config.final_output.value > metadata["max_final_output"]
                or bundle.config.assay.design.value != metadata["design"]
                or bundle.config.input_type.value != metadata.get("input_type", "FASTQ")
            ):
                raise ConfigError(
                    [
                        "Configuration exceeds this selected installation. Bootstrap a new matching profile/assay with the required --analyses; do not delete installation metadata."
                    ]
                )
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
