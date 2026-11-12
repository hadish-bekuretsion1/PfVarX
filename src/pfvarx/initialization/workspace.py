"""Profile-aware initialization after the locked environment has been installed."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from pfvarx.configuration.loader import load_bundle
from pfvarx.initialization.layout import config_init, init_inputs
from pfvarx.io import read_json, write_json


def initialize(root: Path) -> dict[str, Any]:
    path = root / "config/installation.json"
    metadata = read_json(path) if path.exists() else None
    purpose = metadata["purpose"] if metadata else "analysis"
    created = config_init(root, purpose)
    if metadata and root / "config/config.yaml" in created:
        config_path = root / "config/config.yaml"
        config = yaml.safe_load(config_path.read_text())
        config["final_output"] = metadata["max_final_output"]
        config["input_type"] = metadata.get("input_type", "FASTQ")
        config["assay"]["design"] = metadata["design"]
        if metadata["design"] != "WGS":
            config["assay"].update(
                amplification_route="AMPLICON",
                targets_bed="inputs/targets/targets.bed",
                primers_tsv="inputs/targets/primers.tsv",
            )
        if metadata["design"] == "TARGETED_UMI":
            config["assay"]["read_structure"] = "8M+T 8M+T"
        if metadata["analyses_requested"]:
            from pfvarx.configuration.loader import LEGACY_IDS

            requested = metadata["analyses_requested"]
            if metadata["schema"] == "pfvarx.installation/1":
                requested = [LEGACY_IDS.get(da, da) for da in requested]
            config["analyses"]["requested"] = [] if metadata.get("auto_analyses") else requested
            if metadata.get("auto_analyses"):
                config["analyses"]["selection"] = "auto"
            config["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
            if "DA007" in requested:
                config["analyses"]["spatial_outcome_source"] = "DA002"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    result = init_inputs(root, load_bundle(root).config)
    (root / "metadata").mkdir(exist_ok=True)
    report = {
        "status": "INITIALIZED",
        "configuration_created": [str(p.relative_to(root)) for p in created],
        **result,
        "next": "Review all configuration and real library/UMI metadata before registering inputs",
    }
    write_json(root / "metadata/initialization.json", report)
    return report
