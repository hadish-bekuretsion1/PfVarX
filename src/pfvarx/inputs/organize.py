"""Explicit registration/discovery and checksum-verified HTTPS acquisition. Originals are read-only."""

from __future__ import annotations

import hashlib
import re
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import Any

from pfvarx.configuration.loader import load_bundle
from pfvarx.initialization.layout import SAMPLE_COLUMNS
from pfvarx.inputs.content import FOLDER_KINDS, require
from pfvarx.io import read_tsv, write_json, write_tsv

ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")
PAIR = re.compile(
    r"^(?P<sample>[A-Za-z0-9][A-Za-z0-9_-]*?)(?:_(?P<lane>L\d{3}))?_R?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$"
)


def discover(root: Path, directory: Path, target: Path | None = None) -> dict[str, Any]:
    directory = directory.resolve()
    if not directory.is_dir():
        raise ValueError("FASTQ directory does not exist")
    pairs: dict[tuple[str, str], dict[str, str]] = {}
    ignored = []
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        match = PAIR.fullmatch(path.name)
        if not match:
            ignored.append(path.name)
            continue
        sample, lane = match["sample"], match["lane"] or "L001"
        if not ID.fullmatch(sample):
            raise ValueError("Invalid discovered sample ID: " + sample)
        row = pairs.setdefault((sample, lane), {"sample_id": sample, "lane": lane, "library_id": sample})
        field = "fastq_" + match["mate"]
        if field in row:
            raise ValueError("Ambiguous duplicate mate: " + sample + " " + lane)
        # Files inside the workspace are recorded relative to it, so the workspace can be moved.
        base = root.resolve()
        row[field] = str(path.relative_to(base)) if path.is_relative_to(base) else str(path.resolve())
    if not pairs:
        raise ValueError("No supported paired FASTQ names (expected _R1/_R2 or _1/_2); provide an explicit registration TSV")
    if any(not {"fastq_1", "fastq_2"} <= set(row) for row in pairs.values()):
        raise ValueError("Unmatched FASTQ mate; no sheet written")
    target = target if target is not None else root / "metadata/discovered_samples.tsv"
    if target.exists():
        raise ValueError("Discovery sheet exists; review it or choose a new workspace")
    write_tsv(target, list(pairs.values()), SAMPLE_COLUMNS)
    return {
        "status": "REVIEW_REQUIRED",
        "candidate_sheet": str(target.relative_to(root)),
        "samples": len({k[0] for k in pairs}),
        "lanes": len(pairs),
        "ignored_files": ignored,
        "warning": "Sample/lane IDs come from filenames; library_id defaults to sample_id. Verify wet-lab identities and replace generated IDs before register-inputs.",
    }


def register(root: Path, manifest: Path, mode: str = "reference") -> dict[str, Any]:
    if mode not in {"reference", "link"}:
        raise ValueError("mode must be reference or link")
    cfg = load_bundle(root).config
    target = root / cfg.samples_tsv
    if target.exists() and read_tsv(target):
        raise ValueError("Registered sample sheet exists; edit/review it explicitly rather than overwriting")
    manifest = manifest.resolve()
    rows = read_tsv(manifest)
    if not rows:
        raise ValueError("Registration sheet has no samples")
    fields = ("fastq_1", "fastq_2") if cfg.input_type.value == "FASTQ" else ("bam",)
    operations = []
    seen = set()
    destinations = set()
    for row in rows:
        for key in ("sample_id", "lane", "library_id"):
            if not ID.fullmatch(row.get(key, "")):
                raise ValueError("Invalid/missing registration " + key)
        sample_lane = (row["sample_id"], row["lane"])
        if sample_lane in seen:
            raise ValueError("Duplicate sample/lane registration")
        seen.add(sample_lane)
        for field in fields:
            if not row.get(field):
                raise ValueError("Missing registration file " + field)
            path = (manifest.parent / row[field]).resolve()
            if not path.is_file():
                raise ValueError("Missing input: " + str(path))
            if mode == "reference":
                row[field] = str(path)
            else:
                role = "fastq" if field.startswith("fastq") else "bam"
                suffix = "".join(path.suffixes)
                name = f"{row['sample_id']}.{row['lane']}.{field}{suffix}"
                require(path, FOLDER_KINDS[role], name)
                destination = root / "inputs" / role / name
                if destination.exists() or destination.is_symlink() or destination in destinations:
                    raise ValueError("Organized path exists or duplicated: " + str(destination))
                if not destination.parent.resolve().is_relative_to(root.resolve()):
                    raise ValueError("Input destination escapes workspace")
                destinations.add(destination)
                operations.append((path, destination))
                row[field] = str(destination.relative_to(root))
                # Input alignments require an existing sibling index; link that index too.
                if field == "bam":
                    indices = [Path(str(path) + ext) for ext in (".bai", ".csi", ".crai")]
                    indices += [path.with_suffix(".bai"), path.with_suffix(".crai")]
                    index = next((p for p in indices if p.is_file()), None)
                    if index is not None:
                        ext = ".crai" if index.suffix == ".crai" else index.suffix
                        linked = Path(str(destination) + ext)
                        if linked.exists() or linked.is_symlink():
                            raise ValueError("Organized index exists")
                        operations.append((index, linked))
    columns = SAMPLE_COLUMNS + [k for row in rows for k in row if k not in SAMPLE_COLUMNS]
    columns = list(dict.fromkeys(columns))
    created = []
    try:
        for source, destination in operations:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.symlink_to(source)
            created.append(destination)
        write_tsv(target, rows, columns)
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise
    report = {
        "status": "REGISTERED",
        "mode": mode,
        "samples": len({r["sample_id"] for r in rows}),
        "lanes": len(rows),
        "sample_sheet": str(target.relative_to(root)),
        "source_manifest": str(manifest),
        "originals_modified": False,
    }
    write_json(root / "metadata/registration.json", report)
    return report


class HTTPSRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != "https":
            raise ValueError("HTTPS downgrade refused")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def fetch(root: Path, url: str, expected: str | None, role: str, name: str,
          check: Callable[[Path], None] | None = None) -> dict[str, Any]:
    """Download one HTTPS resource; *check* inspects the downloaded content before it is kept."""
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https" or parsed.username or parsed.password:
        raise ValueError("Use an HTTPS URL without embedded credentials")
    if expected is not None and not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
        raise ValueError("A trusted SHA-256 checksum is required")
    if role not in FOLDER_KINDS:
        raise ValueError("Unknown resource role")
    if Path(name).name != name or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
        raise ValueError("Use a plain safe filename")
    destination = root / "inputs" / role / name
    if not destination.parent.resolve().is_relative_to(root.resolve()):
        raise ValueError("Download path escapes workspace")
    if destination.exists() or destination.is_symlink():
        raise ValueError("Destination exists; refusing overwrite")
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".download-", delete=False)
    temp = Path(handle.name)
    digest = hashlib.sha256()
    size = 0
    try:
        opener = urllib.request.build_opener(HTTPSRedirect())
        with handle, opener.open(url, timeout=60) as response:
            if urllib.parse.urlsplit(response.geturl()).scheme != "https":
                raise ValueError("HTTPS required")
            while block := response.read(1024 * 1024):
                handle.write(block)
                digest.update(block)
                size += len(block)
        if expected is not None and digest.hexdigest() != expected.lower():
            raise ValueError("Download checksum mismatch; partial resource removed")
        (check or (lambda staged: require(staged, FOLDER_KINDS[role], name)))(temp)
        # Atomic no-clobber publication: hard-link refuses a racing existing destination.
        import os

        os.link(temp, destination)
    except (urllib.error.URLError, TimeoutError) as exc:
        raise ValueError("Download failed (" + type(exc).__name__ + "); check HTTPS source and connectivity") from exc
    finally:
        temp.unlink(missing_ok=True)
    report = {
        "status": "VERIFIED" if expected is not None else "DOWNLOADED_UNVERIFIED",
        "path": str(destination.relative_to(root)),
        "bytes": size,
        "sha256": digest.hexdigest(),
        "source": urllib.parse.urlunsplit((parsed.scheme, parsed.netloc, parsed.path, "", "")),
    }
    write_json(root / "metadata/downloads" / (name + ".json"), report)
    return report
