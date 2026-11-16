"""Build an optional selected-profile Apptainer image, or emit its complete build context.
Context emission is not image validation. Requires the local pinned Pixi executable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import runpy
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


DEFAULT_BASE = "ubuntu:24.04"


def context(destination: Path, profile: str, design: str, analyses: list[str], rtg: bool = False,
            base: str = DEFAULT_BASE) -> dict:
    if "\n" in base or not base.strip():
        raise ValueError("Invalid base image reference")
    if destination.exists():
        raise ValueError("Build context already exists")
    pixi = shutil.which("pixi")
    if not pixi:
        raise ValueError("Install the pinned Pixi executable first")
    bootstrap = runpy.run_path(str(ROOT / "bootstrap.py"))
    destination.mkdir(parents=True)
    metadata = bootstrap["export_workspace"](ROOT, destination / "workspace", profile, design, analyses, rtg)
    shutil.copy2(pixi, destination / "pixi")
    env = metadata["environment"]
    recipe = f"""Bootstrap: docker
From: {base}

%files
    {destination / "pixi"} /usr/local/bin/pixi
    {destination / "workspace"} /opt/pfvarx

%post
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends ca-certificates
    rm -rf /var/lib/apt/lists/*
    cd /opt/pfvarx
    pixi install --locked -e {env}
    pixi run --locked -e {env} pfvarx version

%environment
    export PATH=/opt/pfvarx/.pixi/envs/{env}/bin:$PATH
    export PYTHONDONTWRITEBYTECODE=1
    export PYTHONPATH=/workspace/src

%runscript
    cd /workspace || exit 2
    exec /opt/pfvarx/.pixi/envs/{env}/bin/python -m pfvarx.cli "$@"

%test
    /opt/pfvarx/.pixi/envs/{env}/bin/pfvarx --help
"""
    (destination / "Apptainer.def").write_text(recipe)
    report = {
        "installation": metadata,
        "pixi_sha256": hashlib.sha256(Path(pixi).read_bytes()).hexdigest(),
        "lock_sha256": hashlib.sha256((ROOT / "pixi.lock").read_bytes()).hexdigest(),
        "base_image": base,
        "base_image_pinned": "@sha256:" in base,
        "status": "CONTEXT_ONLY_NOT_IMAGE_VALIDATION",
        "base_image_note": ("Pinned by digest" if "@sha256:" in base else
                            "Not pinned: pass --base ubuntu:24.04@sha256:<digest> for a reproducible image"),
    }
    (destination / "build_manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    target = p.add_mutually_exclusive_group(required=True)
    target.add_argument("--image", type=Path)
    target.add_argument("--emit-context", type=Path)
    p.add_argument("--profile", choices=["sample", "cohort", "downstream", "benchmark"], default="sample")
    p.add_argument("--design", choices=["wgs", "targeted_umi", "targeted_no_umi"], default="wgs")
    p.add_argument("--analyses", nargs="*", default=[])
    p.add_argument("--rtg", action="store_true")
    p.add_argument("--base", default=DEFAULT_BASE,
                   help="base image; pin it by digest (ubuntu:24.04@sha256:...) for a reproducible build")
    a = p.parse_args()
    if a.emit_context:
        print(json.dumps(context(a.emit_context.resolve(), a.profile, a.design, a.analyses, a.rtg, a.base), indent=2))
        return
    if a.image.exists():
        p.error("Image already exists; no overwrite")
    if not shutil.which("apptainer"):
        p.error("Apptainer is not installed; --emit-context is available for review")
    with tempfile.TemporaryDirectory(prefix="pfvarx-image-") as tmp:
        build = Path(tmp) / "context"
        report = context(build, a.profile, a.design, a.analyses, a.rtg, a.base)
        subprocess.run(["apptainer", "build", str(a.image.resolve()), str(build / "Apptainer.def")], check=True)
        report["status"] = "BUILT_REQUIRES_E2E_QUALIFICATION"
        report["image_sha256"] = hashlib.sha256(a.image.read_bytes()).hexdigest()
        a.image.with_suffix(".build.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
