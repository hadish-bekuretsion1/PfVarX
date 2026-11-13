"""Conditional, idempotent derivation of reference indexes (L0.4).

Raw inputs are never modified. Derivatives are written under ``inputs/reference/prepared/`` and
built only when missing or when the source checksum changed.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

from pfvarx.configuration.models import Bundle
from pfvarx.constants import InputType
from pfvarx.io import read_json, sha256, write_json

BWA_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")


def _tool(name: str) -> str:
    found = shutil.which(name)
    if not found:
        raise FileNotFoundError(f"TOOL_NOT_FOUND: {name} is not on PATH (run inside 'pixi run')")
    return found


def _source_key(root: Path, source: Path) -> str:
    """Workspace-relative key for sources inside the workspace, so a moved workspace stays current."""
    resolved = (root / source).resolve()
    base = root.resolve()
    return resolved.relative_to(base).as_posix() if resolved.is_relative_to(base) else str(resolved)


def plan_preparation(root: Path, bundle: Bundle) -> dict[str, Any]:
    cfg = bundle.config
    src = root / cfg.reference.fasta
    prep = root / "inputs/reference/prepared"
    manifest_path = root / "prepared" / "manifest.json"
    old = read_json(manifest_path) if manifest_path.is_file() else {}
    fingerprints = {_source_key(root, src): sha256(src)}
    recorded = {_source_key(root, Path(k)): v for k, v in (old.get("source_fingerprints") or {}).items()}
    stale = recorded != fingerprints
    sealed = old.get("artifacts", {})
    def damaged(path: Path) -> bool:
        rel = str(path.relative_to(root))
        return not path.is_file() or path.stat().st_size == 0 or (rel in sealed and sha256(path) != sealed[rel])
    target = prep / "combined.fasta"
    actions = []
    if stale or damaged(target):
        actions.append({"action": "write_combined_fasta", "output": str(target)})
    if stale or damaged(Path(f"{target}.fai")):
        actions.append({"action": "samtools_faidx", "output": f"{target}.fai"})
    if stale or damaged(prep / "combined.dict"):
        actions.append({"action": "picard_dict", "output": str(prep / "combined.dict")})
    if cfg.input_type == InputType.FASTQ and (stale or any(damaged(Path(f"{target}{s}")) for s in BWA_SUFFIXES)):
        actions.append({"action": "bwa_mem2_index", "output": f"{target}.*"})
    parasite = prep / "parasite.fasta"
    if stale or damaged(parasite) or damaged(Path(f"{parasite}.fai")):
        actions.append({"action": "write_parasite_fasta", "output": str(parasite)})
    if not actions and not sealed:
        actions.append({"action": "seal_artifacts", "output": str(manifest_path)})
    return {"stale": stale, "actions": actions, "source_fingerprints": fingerprints,
            "manifest": str(manifest_path)}


def _run_preparation(root: Path, action: str, command: list[str]) -> None:
    log = root / "metadata/preparation_logs" / (action + ".log")
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a") as handle:
        handle.write("Starting reference action: " + action + "\n")
        handle.flush()
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        handle.write("Exit status: " + str(result.returncode) + "\n")
    if result.returncode:
        raise ValueError(f"Reference preparation failed ({result.returncode}); inspect {log.relative_to(root)}")


def prepare(root: Path, bundle: Bundle, dry_run: bool) -> dict[str, Any]:
    plan = plan_preparation(root, bundle)
    if dry_run or not plan["actions"]:
        return {"dry_run": dry_run, **plan}
    cfg = bundle.config
    prep = root / "inputs/reference/prepared"
    prep.mkdir(parents=True, exist_ok=True)
    combined = prep / "combined.fasta"
    parasite = prep / "parasite.fasta"
    done = []
    for act in plan["actions"]:
        a = act["action"]
        if a == "write_combined_fasta":
            with combined.open("wb") as out:
                out.write((root / cfg.reference.fasta).read_bytes().rstrip(b"\n") + b"\n")
        elif a == "write_parasite_fasta":
            shutil.copyfile(root / cfg.reference.fasta, parasite)
            _run_preparation(root, a, [_tool("samtools"), "faidx", str(parasite)])
        elif a == "samtools_faidx":
            _run_preparation(root, a, [_tool("samtools"), "faidx", str(combined)])
        elif a == "picard_dict":
            dictionary = prep / "combined.dict"
            dictionary.unlink(missing_ok=True)
            _run_preparation(root, a, [_tool("picard"), "CreateSequenceDictionary", "-R", str(combined),
                            "-O", str(dictionary)])
        elif a == "bwa_mem2_index":
            _run_preparation(root, a, [_tool("bwa-mem2"), "index", str(combined)])
        done.append(a)
    manifest = {
        "schema": "pfvarx.prepared/2",
        "source_fingerprints": plan["source_fingerprints"],
        "reference_sha256": sha256(root / cfg.reference.fasta),
        "combined_fasta": str(combined.relative_to(root)),
        "parasite_fasta": str(parasite.relative_to(root)),
        "actions_performed": done,
    }
    artifacts = required_artifacts(root, bundle)
    if any(not p.is_file() or p.stat().st_size == 0 for p in artifacts):
        raise ValueError("Reference preparation produced missing/empty artifacts; inspect tool output and rerun prepare")
    manifest["artifacts"] = {str(p.relative_to(root)): sha256(p) for p in artifacts}
    write_json(root / "prepared" / "manifest.json", manifest)
    return {"dry_run": False, **plan, "performed": done}


def required_artifacts(root: Path, bundle: Bundle) -> list[Path]:
    prep = root / "inputs/reference/prepared"
    paths = [prep / name for name in ("combined.fasta", "combined.fasta.fai", "combined.dict", "parasite.fasta", "parasite.fasta.fai")]
    if bundle.config.input_type == InputType.FASTQ:
        paths += [Path(str(prep / "combined.fasta") + suffix) for suffix in BWA_SUFFIXES]
    return paths


def verify_prepared(root: Path, bundle: Bundle) -> dict[str, Any]:
    """Verify source fingerprints, derived hashes and FASTA/FAI/dictionary compatibility."""
    import pysam

    issues = []
    path = root / "prepared/manifest.json"
    if not path.is_file():
        return {"status": "FAIL", "issues": ["Missing prepared/manifest.json; run prepare"], "artifacts_checked": 0}
    try:
        manifest = read_json(path)
        plan = plan_preparation(root, bundle)
        if plan["actions"]:
            issues.append("Missing, stale, unsealed or damaged preparation; run prepare")
        paths = required_artifacts(root, bundle)
        for artifact in paths:
            rel = str(artifact.relative_to(root))
            if not artifact.is_file() or artifact.stat().st_size == 0:
                issues.append("Missing/empty artifact: " + rel)
            elif manifest.get("artifacts", {}).get(rel) != sha256(artifact):
                issues.append("Artifact hash mismatch or not sealed: " + rel)
        if not issues:
            with pysam.FastaFile(str(root / "inputs/reference/prepared/combined.fasta")) as combined:
                expected = dict(zip(combined.references, combined.lengths))
                dictionary = {}
                for line in (root / "inputs/reference/prepared/combined.dict").read_text().splitlines():
                    if line.startswith("@SQ\t"):
                        fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
                        dictionary[fields["SN"]] = int(fields["LN"])
                if dictionary != expected:
                    issues.append("Combined FASTA index and sequence dictionary disagree")
            with pysam.FastaFile(str(root / "inputs/reference/prepared/parasite.fasta")) as parasite:
                if not parasite.references:
                    issues.append("Empty parasite FASTA index")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        issues.append("Cannot verify preparation: " + str(exc))
    return {"status": "FAIL" if issues else "PASS", "issues": issues,
            "artifacts_checked": len(required_artifacts(root, bundle))}
