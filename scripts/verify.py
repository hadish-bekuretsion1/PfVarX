"""Verify an unpacked/current checkout using isolated synthetic workspaces.

Run with `pixi run python scripts/verify.py`. No Git commit is required.
The fixture catalogue is for software testing only. Original inputs are never used.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def verify(scenario: str, destination: Path, dry_run_only: bool = False) -> dict:
    ws = destination / scenario
    if ws.exists():
        raise ValueError(f"Verification workspace already exists: {ws}; choose a new --output")
    ws.mkdir(parents=True)
    for name in ("src", "config", "workflow", "pixi.toml", "pixi.lock", "pyproject.toml"):
        source = ROOT / name
        if source.is_dir():
            shutil.copytree(source, ws / name, ignore=shutil.ignore_patterns("__pycache__"))
        else:
            shutil.copyfile(source, ws / name)
    for name in ("config.yaml", "quality.yaml", "compute.yaml", "profiles/resources.yaml"):
        (ws / "config" / name).unlink(missing_ok=True)
    env = dict(os.environ, PATH=f"{Path(sys.executable).parent}:{os.environ['PATH']}",
               PYTHONPATH=str(ws / "src"))
    commands = []

    def run(*args: str, expected: int = 0) -> dict:
        label = " ".join(args)
        print(f"[{scenario}] {label}", flush=True)
        result = subprocess.run([sys.executable, "-m", "pfvarx.cli", *args], cwd=ws,
                                env=env, capture_output=True, text=True)
        log = ws / f"{len(commands):02d}-{args[0]}.log"
        log.write_text(result.stdout + "\n" + result.stderr)
        commands.append({"command": label, "exit_code": result.returncode, "log": str(log)})
        if result.returncode != expected:
            raise RuntimeError(f"{scenario}: {label} exited {result.returncode}; see {log}")
        return json.loads(result.stdout)

    purpose = "benchmarking" if scenario == "benchmark" else "analysis"
    design = {"targeted": "TARGETED_NON_UMI", "auto": "TARGETED_NON_UMI", "umi": "TARGETED_UMI"}.get(scenario, "WGS")
    run("config-init", "--purpose", purpose)
    cfgfile = ws / "config/config.yaml"
    cfg = yaml.safe_load(cfgfile.read_text())
    final_output = {"fo001": "FO001", "fo003": "FO003", "auto": "FO004"}.get(
        scenario, "FO004" if design == "WGS" else "FO002")
    cfg.update(project_id=f"VERIFY_{scenario.upper()}", final_output=final_output)
    cfg["assay"]["design"] = design
    cfg["reference"]["release_label"] = "synthetic"
    if design != "WGS":
        cfg["assay"].update(amplification_route="AMPLICON", targets_bed="inputs/targets/targets.bed",
                             primers_tsv="inputs/targets/primers.tsv")
    if design == "TARGETED_UMI":
        cfg["assay"]["read_structure"] = "8M+T 8M+T"
    if scenario == "auto":
        cfg["assay"]["sequencing_type"] = "AMPLICON"
        cfg["analyses"]["selection"] = "auto"
    elif final_output == "FO004":
        cfg["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
        cfg["analyses"].update(requested=[f"DA{i:03d}" for i in range(1, 8)],
                               scales={"DA002": ["SAMPLE", "POPULATION"], "DA003": ["POPULATION"],
                                       "DA005": ["BETWEEN_POPULATIONS"]},
                               population_column="population", spatial_outcome_source="DA002",
                               resistance_catalog="inputs/metadata/test_resistance_catalog.tsv")
    cfgfile.write_text(yaml.safe_dump(cfg))
    compfile = ws / "config/compute.yaml"
    comp = yaml.safe_load(compfile.read_text())
    comp["max_cores"] = min(4, os.cpu_count() or 1)
    comp["per_job"]["max_threads"] = comp["max_cores"]
    compfile.write_text(yaml.safe_dump(comp))
    run("init")
    fixture = ["fixture", "--design", design, "--purpose", purpose]
    if design != "WGS" or scenario in ("fo001", "fo003"):
        fixture += ["--n-mono", "2", "--n-poly", "0"]
    if purpose == "benchmarking":
        fixture += ["--n-mono", "10", "--low-depth-sample"]
    fixture_info = run(*fixture)
    if purpose == "benchmarking":
        run("benchmark-partition", "--holdout-fraction", "0.25")
    if scenario == "auto":
        sheet = ws / cfg["samples_tsv"]
        rows = [line.split("\t") for line in sheet.read_text().splitlines()]
        optional = {"population", "country", "latitude", "longitude", "collection_date"}
        keep = [i for i, name in enumerate(rows[0]) if name not in optional]
        sheet.write_text("\n".join("\t".join(row[i] for i in keep) for row in rows) + "\n")
    run("check-layout")
    run("validate")
    run("prepare")
    preparation = run("prepare", "--dry-run")
    assert not preparation["actions"], preparation
    assert run("ready")["ready"]
    dry = run("dry-run")
    assert dry["status"] == "OK" and dry["snakemake_exit"] == 0, dry
    assert not (ws / "runs").exists() and not (ws / "results").exists()
    if dry_run_only:
        report = {"scenario": scenario, "status": "PASSED", "scope": "dry-run",
                  "snakemake_exit": dry["snakemake_exit"], "runs_or_results_created": False,
                  "commands": commands}
        (ws / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
        return report
    result = run("pipeline")
    assert result["status"] == "COMPLETED", result
    run("run-verify", "--run-id", result["run_id"])
    expected_samples = set(fixture_info["samples"])
    alignments = {p.stem for p in (ws / "results/alignments").glob("*.bam")}
    assert alignments == expected_samples, (alignments, expected_samples)
    if final_output != "FO001":
        calls = {p.name.removesuffix(".vcf.gz") for p in
                 (ws / "results/variants/samples").glob("*.vcf.gz")}
        if purpose == "benchmarking":
            # The deliberately under-sequenced sample must stop at Layer 1.
            blocked = {p.parent.name for p in (ws / "results/qc/samples").glob("*/qc_decision.json")
                       if json.loads(p.read_text())["verdict"] == "FAIL"}
            assert len(blocked) == 1, blocked
            expected_samples -= blocked
        assert calls == expected_samples, (calls, expected_samples)
    if final_output in ("FO003", "FO004"):
        assert (ws / "results/cohort/cohort.vcf.gz").is_file()
        assert (ws / "results/cohort/evidence_matrix.tsv.gz").is_file()
    report = ws / "results/report/pfvarx_dashboard.html"
    assert report.stat().st_size > 1000
    manifests = {p.parent.name: json.loads(p.read_text()) for p in
                 (ws / "results/downstream").glob("*/status_manifest.json")}
    if scenario == "auto":
        assert set(manifests) == {"DA001", "DA002", "DA003", "DA004", "DA007"}
        assert manifests["DA001"]["terminal_state"].startswith("COMPLETED")
        for da in ("DA002", "DA007"):
            assert manifests[da]["terminal_state"] == "BLOCKED_MISSING_PREREQUISITE"
        assert not (ws / "results/annotation").exists()
    elif final_output == "FO004":
        assert set(manifests) == set(cfg["analyses"]["requested"]), manifests.keys()
        assert all(m["terminal_state"] != "FAILED_EXECUTION" for m in manifests.values()), manifests
    if purpose == "benchmarking":
        run("benchmark-evaluate", "--out", "benchmarks/development")
        run("benchmark-evaluate", "--holdout", "--out", "benchmarks/heldout")
        run("benchmark-evaluate", "--holdout", "--out", "benchmarks/heldout", expected=2)
    summary = {"scenario": scenario, "status": "PASSED", "workspace": str(ws),
               "run_id": result["run_id"], "dashboard": str(report),
               "analyses": {k: v["terminal_state"] for k, v in manifests.items()}, "commands": commands}
    (ws / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    choices = ["wgs", "targeted", "umi", "benchmark", "fo001", "fo003", "auto"]
    parser.add_argument("scenarios", nargs="*", metavar="SCENARIO",
                        help="scenarios: " + ", ".join(choices) + "; default: all")
    parser.add_argument("--output", type=Path, default=ROOT / "work/verification")
    parser.add_argument("--dry-run-only", action="store_true", help="Stop after Snakemake preflight")
    args = parser.parse_args()
    scenarios = args.scenarios or choices
    if unknown := sorted(set(scenarios) - set(choices)):
        parser.error("unknown scenarios: " + ", ".join(unknown))
    for scenario in scenarios:
        verify(scenario, args.output.resolve(), args.dry_run_only)


if __name__ == "__main__":
    main()
