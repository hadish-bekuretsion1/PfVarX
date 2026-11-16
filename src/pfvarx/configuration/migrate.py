"""Migrate a workspace from configuration schema 3 to 4 (``./pfvarx config-migrate``).

Schema 4 makes ``assay_design`` the setting that selects the processing branch and moves
what describes each library out of the run settings:

* ``assay.design``                  -> ``assay_design`` (wgs | targeted_umi | targeted_no_umi)
* ``assay.sequencing_type``         -> manifest column ``technique``
* ``assay.amplification_route``     -> manifest column ``amplification_route``
* ``assay.targets_bed/primers_tsv`` -> a panel in the library catalog (``LEGACY_PANEL``,
                                       version ``config3``), named in each manifest row
* ``assay.read_structure``          -> a UMI library protocol (``LEGACY_UMI_PROTOCOL``)

For imported targeted BAM/CRAM files, schema 3 required primer trimming to have been done
(and, for UMI libraries, consensus), so their rows get the matching ``processing_state``.
WGS imports keep the header-based duplicate-marking check, so their state is left empty.

Only empty manifest cells are filled; existing values are never overwritten. Every file
that is rewritten is first copied to ``<name>.before-config4``.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import yaml

from pfvarx.configuration.library import CATALOG_SCHEMA
from pfvarx.constants import Design, ProcessingState
from pfvarx.io import atomic_write_text, read_tsv, write_tsv

LEGACY_PANEL = ("LEGACY_PANEL", "config3")
LEGACY_UMI_PROTOCOL = "LEGACY_UMI_PROTOCOL"
CATALOG_PATH = "inputs/metadata/library_catalog.yaml"
MANIFEST_COLUMNS = ("technique", "amplification_route", "panel_id", "panel_version",
                    "library_protocol_id", "processing_state")
LEGACY_TECHNIQUES = {"SWGS": "SELECTIVE_WGS"}


def plan_v3(project: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any] | None, dict[str, str]]:
    """Schema-4 project, the library catalog to write (or None), and manifest defaults."""
    assay = dict(project.get("assay") or {})
    design = Design.parse(assay.pop("design"))
    technique = assay.pop("sequencing_type", None)
    route = assay.pop("amplification_route", None)
    targets, primers = assay.pop("targets_bed", None), assay.pop("primers_tsv", None)
    read_structure = assay.pop("read_structure", None)
    defaults: dict[str, str] = {}
    if technique:
        defaults["technique"] = LEGACY_TECHNIQUES.get(str(technique), str(technique))
    if route:
        defaults["amplification_route"] = route
    catalog: dict[str, Any] | None = None
    if design.targeted:
        catalog = {"schema_version": CATALOG_SCHEMA,
                   "panels": [{"panel_id": LEGACY_PANEL[0], "panel_version": LEGACY_PANEL[1],
                               "targets_bed": targets, "primers_tsv": primers,
                               "description": "Migrated from assay.targets_bed and assay.primers_tsv (config 3)"}],
                   "library_protocols": []}
        defaults.update(panel_id=LEGACY_PANEL[0], panel_version=LEGACY_PANEL[1])
        if design == Design.TARGETED_UMI:
            catalog["library_protocols"] = catalog["library_protocols"] + [
                {"library_protocol_id": LEGACY_UMI_PROTOCOL, "umi": {"read_structure": read_structure},
                 "description": "Migrated from assay.read_structure (config 3)"}]
            defaults["library_protocol_id"] = LEGACY_UMI_PROTOCOL
        assay["library_catalog"] = CATALOG_PATH
        if project.get("input_type") in ("BAM", "CRAM"):
            defaults["processing_state"] = (ProcessingState.CONSENSUS_PRIMERS_TRIMMED.value
                                            if design == Design.TARGETED_UMI else
                                            ProcessingState.PRIMERS_TRIMMED.value)
    migrated: dict[str, Any] = {"schema_version": "pfvarx.config/4"}
    for key, value in project.items():
        if key == "schema_version":
            continue
        if key == "assay":
            migrated["assay_design"] = design.value
            migrated["assay"] = assay
        else:
            migrated[key] = value
    return migrated, catalog, defaults


def migrate_workspace(root: Path, apply: bool = False) -> dict[str, Any]:
    """Report (and with *apply*, perform) the schema 3 -> 4 migration of one workspace."""
    from pfvarx.configuration.loader import USER_FILES, migrate_project, read_yaml

    config_path = root / USER_FILES["config"]
    project = migrate_project(read_yaml(config_path))  # schema 2 -> 3 first
    if project.get("schema_version") != "pfvarx.config/3":
        return {"status": "NOT_NEEDED", "schema_version": project.get("schema_version")}
    migrated, catalog, defaults = plan_v3(project)
    manifest = root / str(project.get("samples_tsv") or "inputs/metadata/samples.tsv")
    catalog_path = root / CATALOG_PATH
    if catalog and catalog_path.exists():
        return {"status": "BLOCKED", "detail": f"{CATALOG_PATH} already exists; merge the panel by hand",
                "catalog": catalog}
    rows = read_tsv(manifest) if manifest.is_file() else []
    filled = 0
    for row in rows:
        for column, value in defaults.items():
            if not row.get(column):
                row[column] = value
                filled += 1
    report = {"status": "APPLIED" if apply else "PLANNED", "config": str(config_path.relative_to(root)),
              "assay_design": migrated["assay_design"], "manifest_defaults": defaults,
              "manifest_cells_filled": filled, "catalog_written": CATALOG_PATH if catalog else None}
    if defaults and not rows:
        # Nothing to fill yet: say what the old settings recorded, so it is not lost.
        report["manifest_not_found"] = (
            f"{manifest.relative_to(root) if manifest.is_relative_to(root) else manifest} has no rows yet; when you "
            f"create it, add these values to every row: {defaults}")
    if not apply:
        return report
    backups = [config_path] + ([manifest] if rows and filled else [])
    for path in backups:
        shutil.copy2(path, path.with_name(path.name + ".before-config4"))
    if catalog:
        catalog_path.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_text(catalog_path, yaml.safe_dump(catalog, sort_keys=False))
    if rows and filled:
        columns = list(rows[0]) + [c for c in MANIFEST_COLUMNS if c in defaults and c not in rows[0]]
        write_tsv(manifest, rows, columns)
    atomic_write_text(config_path, yaml.safe_dump(migrated, sort_keys=False))
    report["backups"] = [str(p.with_name(p.name + ".before-config4").relative_to(root)) for p in backups]
    return report
