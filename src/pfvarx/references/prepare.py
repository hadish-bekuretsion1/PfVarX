"""Conditional, idempotent derivation of reference indexes (L0.4).

Raw inputs are never modified, and ``inputs/`` holds only the files a user supplies.
Derivatives are written under ``prepared/reference/`` (sealed in ``prepared/manifest.json``)
and built only when missing or when the source checksum changed. Workspaces prepared before
this layout kept them in ``inputs/reference/prepared/``; ``prepare`` moves them once.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

from pfvarx.configuration.models import Bundle
from pfvarx.constants import Design, InputType
from pfvarx.io import read_json, sha256, write_json
from pfvarx.references import host as host_reference

BWA_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
PREPARED_REFERENCE = "prepared/reference"
LEGACY_PREPARED_REFERENCE = "inputs/reference/prepared"
# Files at least this large (a human-sized combined FASTA and its index) are sealed by size and
# modification time instead of SHA-256, so verification does not re-read tens of GB each time.
LARGE_ARTIFACT_BYTES = 1 << 30


def fingerprint(path: Path) -> str:
    """SHA-256 of a file, or ``size:<bytes>:mtime:<ns>`` for files of 1 GiB or more."""
    st = path.stat()
    if st.st_size >= LARGE_ARTIFACT_BYTES:
        return f"size:{st.st_size}:mtime:{st.st_mtime_ns}"
    return sha256(path)


def needs_aligner_index(cfg: Any) -> bool:
    """FASTQ entry aligns reads; targeted_umi BAM entry may realign new consensus reads."""
    return cfg.input_type == InputType.FASTQ or cfg.assay_design == Design.TARGETED_UMI


def legacy_location(root: Path) -> bool:
    """Derivatives still in the pre-0.0.0-layout folder inside inputs/, not yet moved."""
    return (root / LEGACY_PREPARED_REFERENCE).is_dir() and not (root / PREPARED_REFERENCE).exists()


def relocate_legacy(root: Path) -> bool:
    """Move derivatives from inputs/reference/prepared/ to prepared/reference/ (no rebuild)."""
    if not legacy_location(root):
        return False
    (root / PREPARED_REFERENCE).parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(root / LEGACY_PREPARED_REFERENCE), str(root / PREPARED_REFERENCE))
    manifest_path = root / "prepared" / "manifest.json"
    if manifest_path.is_file():
        manifest = read_json(manifest_path)
        old, new = LEGACY_PREPARED_REFERENCE + "/", PREPARED_REFERENCE + "/"
        manifest["artifacts"] = {k.replace(old, new, 1): v for k, v in (manifest.get("artifacts") or {}).items()}
        for key in ("combined_fasta", "parasite_fasta"):
            if isinstance(manifest.get(key), str):
                manifest[key] = manifest[key].replace(old, new, 1)
        write_json(manifest_path, manifest)
    return True


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
    prep = root / PREPARED_REFERENCE
    manifest_path = root / "prepared" / "manifest.json"
    old = read_json(manifest_path) if manifest_path.is_file() else {}
    fingerprints = {_source_key(root, src): sha256(src)}
    host = cfg.reference.host_fasta_used
    if host is not None:
        fingerprints[_source_key(root, host)] = fingerprint(root / host)
    recorded = {_source_key(root, Path(k)): v for k, v in (old.get("source_fingerprints") or {}).items()}
    stale = recorded != fingerprints
    sealed = old.get("artifacts", {})
    def damaged(path: Path) -> bool:
        rel = str(path.relative_to(root))
        return not path.is_file() or path.stat().st_size == 0 or (rel in sealed and fingerprint(path) != sealed[rel])
    target = prep / "combined.fasta"
    actions = []
    if legacy_location(root):
        # Moved first; prepare then re-plans, so only what is really missing gets rebuilt.
        return {"stale": stale, "source_fingerprints": fingerprints, "manifest": str(manifest_path),
                "actions": [{"action": "relocate_from_inputs", "output": PREPARED_REFERENCE + "/"}]}
    if stale or damaged(target):
        actions.append({"action": "write_combined_fasta", "output": str(target)})
    if stale or damaged(Path(f"{target}.fai")):
        actions.append({"action": "samtools_faidx", "output": f"{target}.fai"})
    if stale or damaged(prep / "combined.dict"):
        actions.append({"action": "picard_dict", "output": str(prep / "combined.dict")})
    index_damaged = stale or any(damaged(Path(f"{target}{s}")) for s in BWA_SUFFIXES)
    if needs_aligner_index(cfg) and host is not None:
        # A host index is supplied by the user (too large to build here); check it matches.
        unsealed = any(str(Path(f"{target}{s}").relative_to(root)) not in sealed for s in BWA_SUFFIXES)
        if index_damaged or unsealed or host_reference.index_problems(target, BWA_SUFFIXES):
            actions.append({"action": "check_host_index", "output": f"{target}.*"})
    elif needs_aligner_index(cfg) and index_damaged:
        actions.append({"action": "bwa_mem2_index", "output": f"{target}.*"})
    parasite = prep / "parasite.fasta"
    if stale or damaged(parasite) or damaged(Path(f"{parasite}.fai")):
        actions.append({"action": "write_parasite_fasta", "output": str(parasite)})
    if host is not None and (stale or damaged(prep / "parasite.dict")):
        actions.append({"action": "picard_dict_parasite", "output": str(prep / "parasite.dict")})
    if not actions and not sealed:
        actions.append({"action": "seal_artifacts", "output": str(manifest_path)})
    return {"stale": stale, "actions": actions, "source_fingerprints": fingerprints,
            "manifest": str(manifest_path)}


def old_host_facts(root: Path) -> dict[str, Any] | None:
    """Host facts recorded by the last prepare (kept when the combined FASTA is not rewritten)."""
    path = root / "prepared" / "manifest.json"
    return read_json(path).get("host") if path.is_file() else None


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
    moved = relocate_legacy(root)
    if moved:
        plan = plan_preparation(root, bundle)  # the moved derivatives are reused, not rebuilt
        if not plan["actions"]:
            return {"dry_run": False, **plan, "performed": ["relocate_from_inputs"]}
    cfg = bundle.config
    prep = root / PREPARED_REFERENCE
    prep.mkdir(parents=True, exist_ok=True)
    combined = prep / "combined.fasta"
    parasite = prep / "parasite.fasta"
    done = []
    host_facts = old_host_facts(root)
    index_missing: list[str] = []
    for act in plan["actions"]:
        a = act["action"]
        if a == "write_combined_fasta" and cfg.reference.host_fasta_used is not None:
            host_facts = host_reference.write_combined(root / cfg.reference.fasta, root / cfg.reference.host_fasta_used,
                                                       combined)
            host_facts["fasta"] = str(cfg.reference.host_fasta_used)
        elif a == "write_combined_fasta":
            host_facts = None
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
        elif a == "picard_dict_parasite":
            dictionary = prep / "parasite.dict"
            dictionary.unlink(missing_ok=True)
            _run_preparation(root, a, [_tool("picard"), "CreateSequenceDictionary", "-R", str(parasite),
                            "-O", str(dictionary)])
        elif a == "bwa_mem2_index":
            _run_preparation(root, a, [_tool("bwa-mem2"), "index", str(combined)])
        elif a == "check_host_index":
            index_missing = host_reference.index_problems(combined, BWA_SUFFIXES)
            if index_missing:
                continue
        done.append(a)
    if moved:
        done.insert(0, "relocate_from_inputs")
    manifest = {
        "schema": "pfvarx.prepared/2",
        "source_fingerprints": plan["source_fingerprints"],
        "reference_sha256": sha256(root / cfg.reference.fasta),
        "combined_fasta": str(combined.relative_to(root)),
        "parasite_fasta": str(parasite.relative_to(root)),
        "host": host_facts,
        "actions_performed": done,
    }
    artifacts = required_artifacts(root, bundle)
    if index_missing:
        # Seal what was prepared so the next run only checks the supplied index, then stop.
        present = [p for p in artifacts if p.is_file() and p.stat().st_size > 0 and not
                   any(str(p).endswith(s) for s in BWA_SUFFIXES)]
        manifest["artifacts"] = {str(p.relative_to(root)): fingerprint(p) for p in present}
        manifest["host_index"] = "REQUIRED"
        write_json(root / "prepared" / "manifest.json", manifest)
        steps = host_reference.index_instructions(prep)
        raise ValueError("HOST_INDEX_REQUIRED: " + "; ".join(index_missing) + ". Steps: "
                         + " ".join(f"{i}. {t}" for i, t in enumerate(steps, 1)))
    if any(not p.is_file() or p.stat().st_size == 0 for p in artifacts):
        raise ValueError("Reference preparation produced missing/empty artifacts; inspect tool output and rerun prepare")
    manifest["artifacts"] = {str(p.relative_to(root)): fingerprint(p) for p in artifacts}
    write_json(root / "prepared" / "manifest.json", manifest)
    return {"dry_run": False, **plan, "performed": done}


def required_artifacts(root: Path, bundle: Bundle) -> list[Path]:
    prep = root / PREPARED_REFERENCE
    paths = [prep / name for name in ("combined.fasta", "combined.fasta.fai", "combined.dict", "parasite.fasta", "parasite.fasta.fai")]
    if bundle.config.reference.host_fasta_used is not None:
        paths.append(prep / "parasite.dict")
    if needs_aligner_index(bundle.config):
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
        if bundle.config.reference.host_fasta_used is not None and needs_aligner_index(bundle.config):
            combined_path = root / PREPARED_REFERENCE / "combined.fasta"
            if combined_path.is_file():
                for problem in host_reference.index_problems(combined_path, BWA_SUFFIXES):
                    issues.append("HOST_INDEX_REQUIRED: " + problem + "; see ./pfvarx prepare for the steps")
        paths = required_artifacts(root, bundle)
        for artifact in paths:
            rel = str(artifact.relative_to(root))
            if not artifact.is_file() or artifact.stat().st_size == 0:
                issues.append("Missing/empty artifact: " + rel)
            elif manifest.get("artifacts", {}).get(rel) != fingerprint(artifact):
                issues.append("Artifact hash mismatch or not sealed: " + rel)
        if not issues:
            with pysam.FastaFile(str(root / PREPARED_REFERENCE / "combined.fasta")) as combined:
                expected = dict(zip(combined.references, combined.lengths))
                dictionary = {}
                for line in (root / PREPARED_REFERENCE / "combined.dict").read_text().splitlines():
                    if line.startswith("@SQ\t"):
                        fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
                        dictionary[fields["SN"]] = int(fields["LN"])
                if dictionary != expected:
                    issues.append("Combined FASTA index and sequence dictionary disagree")
            with pysam.FastaFile(str(root / PREPARED_REFERENCE / "parasite.fasta")) as parasite:
                if not parasite.references:
                    issues.append("Empty parasite FASTA index")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        issues.append("Cannot verify preparation: " + str(exc))
    return {"status": "FAIL" if issues else "PASS", "issues": issues,
            "artifacts_checked": len(required_artifacts(root, bundle))}
