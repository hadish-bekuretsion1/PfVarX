"""Check that each installation gets the tools its final output uses, and no others.

For each typical installation (FO001 to FO004, the three designs, BAM input, chosen
analyses, benchmarking) this asks bootstrap.py which environment and source files it would
install, reads that environment's packages from pixi.lock (nothing is installed or
downloaded), and compares them with the tools that final output actually runs:

    pixi run python scripts/check_installs.py

    MISSING  a tool the final output runs is not in its environment (an error)
    EXTRA    a tool is installed that this final output never runs (reported)
    ok       needed and installed, or neither

It also lists which parts of the pipeline code each installation copies. The exit code is 1
when anything is MISSING. To check a real installation instead, see docs/deployment.md
("Check what an installation contains").
"""

from __future__ import annotations

import runpy
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
BOOTSTRAP = runpy.run_path(str(ROOT / "bootstrap.py"))

# Tools every final output runs, and the tools each part of the pipeline adds.
BASE = {"python", "snakemake-minimal", "samtools", "htslib", "pysam", "picard-slim", "mosdepth"}
TOOLS = ["fastp", "bwa-mem2", "fgbio", "bcftools", "numpy", "scipy", "snpeff", "r-base", "r-moire",
         "r-dcifer", "hmmibd", "scikit-allel", "iqtree", "r-mgcv", "rtg-tools"]
# The conda package that provides each tool, where the names differ.
PACKAGE = {"fgbio": "fgbio-minimal"}
# Full builds that bring R for optional plots; PfVarX uses the -slim/-minimal builds instead.
WITH_R = ("picard", "fgbio")
CODE = {"caller": "src/pfvarx/caller", "cohort": "src/pfvarx/cohort", "annotation": "src/pfvarx/annotation",
        "downstream": "src/pfvarx/downstream", "benchmarking": "src/pfvarx/benchmarking"}

# (label, bootstrap arguments: profile, design, analyses, rtg, input type, final output)
SCENARIOS = [
    ("FO001 WGS FASTQ", (None, "wgs", [], False, "FASTQ", "FO001")),
    ("FO001 targeted UMI FASTQ", (None, "targeted_umi", [], False, "FASTQ", "FO001")),
    ("FO002 WGS FASTQ", (None, "wgs", [], False, "FASTQ", "FO002")),
    ("FO002 WGS BAM", (None, "wgs", [], False, "BAM", "FO002")),
    ("FO002 targeted no-UMI FASTQ", (None, "targeted_no_umi", [], False, "FASTQ", "FO002")),
    ("FO003 WGS FASTQ", (None, "wgs", [], False, "FASTQ", "FO003")),
    ("FO004 WGS DA002 (resistance)", (None, "wgs", ["DA002"], False, "FASTQ", "FO004")),
    ("FO004 WGS DA001 (clonality)", (None, "wgs", ["DA001"], False, "FASTQ", "FO004")),
    ("FO004 WGS DA006 (phylogeny)", (None, "wgs", ["DA006"], False, "FASTQ", "FO004")),
    ("FO004 WGS DA007 (spatial)", (None, "wgs", ["DA007"], False, "FASTQ", "FO004")),
    ("Benchmark with RTG", ("benchmark", "wgs", [], True, "FASTQ", None)),
]


def needed(design: str, analyses: set[str], rtg: bool, input_type: str, final_output: str) -> set[str]:
    """The tools this installation's final output runs (workflow rules and analysis scripts)."""
    tools = set(BASE)
    if input_type == "FASTQ":
        tools |= {"fastp", "bwa-mem2"}
    if design != "wgs":
        tools.add("fgbio")
        if design == "targeted_umi":
            tools.add("bwa-mem2")  # consensus reads are realigned
    if final_output >= "FO002":
        tools |= {"bcftools", "numpy", "scipy"}
    if final_output == "FO004":
        tools.add("snpeff")
    if analyses & {"DA001", "DA004", "DA005", "DA006"}:
        tools |= {"r-base", "r-moire"}
    if "DA004" in analyses:
        tools |= {"r-dcifer", "hmmibd"}
    if "DA005" in analyses:
        tools.add("scikit-allel")
    if "DA006" in analyses:
        tools.add("iqtree")
    if "DA007" in analyses:
        tools |= {"r-base", "r-mgcv"}
    if rtg:
        tools.add("rtg-tools")
    return tools


def locked_packages(environment: str) -> set[str]:
    lock = yaml.safe_load((ROOT / "pixi.lock").read_text())
    names = set()
    for packages in lock["environments"][environment]["packages"].values():
        for package in packages:
            url = package.get("conda") or package.get("pypi") or ""
            names.add(url.rsplit("/", 1)[-1].rsplit("-", 2)[0])
    return names


def main() -> int:
    missing_any = False
    with tempfile.TemporaryDirectory(prefix="pfvarx-installs-") as tmp:
        for n, (label, args) in enumerate(SCENARIOS):
            profile, design, analyses, rtg, input_type, final_output = args
            profile = BOOTSTRAP["profile_for"](final_output, profile)
            meta = BOOTSTRAP["export_workspace"](ROOT, Path(tmp) / f"ws{n}", profile, design, analyses, rtg,
                                                 input_type, final_output)
            final = meta["max_final_output"]
            installed = locked_packages(meta["environment"])
            want = needed(design, set(meta["analyses_available"]), rtg, input_type, final)
            rows, missing, extra = [], [], []
            for tool in sorted(BASE) + TOOLS:
                has, uses = PACKAGE.get(tool, tool) in installed, tool in want
                if uses and not has:
                    missing.append(tool)
                elif has and not uses:
                    extra.append(tool)
                rows.append((tool, uses, has))
            extra += [f"{name} (full build with R)" for name in WITH_R if name in installed]
            code = [name for name, path in CODE.items() if (Path(tmp) / f"ws{n}" / path).exists()]
            status = "MISSING " + ", ".join(missing) if missing else ("ok" if not extra else "ok, EXTRA " + ", ".join(extra))
            print(f"\n{label}  →  environment {meta['environment']}, up to {final}: {status}")
            print("  installed: " + ", ".join(t for t, _, has in rows if has))
            print("  not installed: " + ", ".join(t for t, _, has in rows if not has))
            print("  pipeline code copied: core" + "".join(", " + c for c in code))
            missing_any |= bool(missing)
    print("\nAll installations have the tools their final output runs." if not missing_any
          else "\nSome installations lack a tool they run (MISSING above).")
    return 1 if missing_any else 0


if __name__ == "__main__":
    sys.exit(main())
