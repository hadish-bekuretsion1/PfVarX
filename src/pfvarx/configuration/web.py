"""Loopback-only configuration editor and CLI preflight launcher."""

from __future__ import annotations

import json
import os
import secrets
import subprocess
import sys
import webbrowser
from functools import partial
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlsplit

import yaml

from pfvarx.configuration.loader import (
    USER_FILES,
    canonical_json,
    migrate_project,
    read_yaml,
    sha256_text,
    validate_bundle,
)
from pfvarx.configuration.models import Bundle
from pfvarx.execution.retention import workspace_lock
from pfvarx.initialization.layout import create_input_dirs, required_dirs
from pfvarx.io import atomic_write_text, read_json, within


def revision(root: Path) -> str:
    files = {name: (root / rel).read_text() if (root / rel).is_file() else None
             for name, rel in USER_FILES.items()}
    return sha256_text(canonical_json(files))


def editor_state(root: Path, purpose: str) -> dict[str, Any]:
    if purpose not in {"analysis", "benchmarking"}:
        raise ValueError("Choose analysis or benchmarking")
    templates = {"config": purpose, "quality": "quality", "compute": "compute"}
    parsed = {}
    for name, rel in USER_FILES.items():
        path = root / rel
        parsed[name] = read_yaml(path if path.is_file() else root / f"config/templates/{templates[name]}.yaml")
    project = migrate_project(parsed["config"])
    installation = root / "config/installation.json"
    if installation.is_file() and not (root / USER_FILES["config"]).is_file():
        meta = read_json(installation)
        project.update(final_output=meta["max_final_output"], input_type=meta.get("input_type", "FASTQ"))
        project["assay"]["design"] = meta["design"]
        if meta["design"] != "WGS":
            project["assay"].update(amplification_route="AMPLICON", targets_bed="inputs/targets/targets.bed",
                                    primers_tsv="inputs/targets/primers.tsv")
        if meta.get("auto_analyses"):
            project["analyses"]["selection"] = "auto"
        else:
            ids = meta["analyses_requested"]
            if meta["schema"] == "pfvarx.installation/1":
                from pfvarx.configuration.loader import LEGACY_IDS

                ids = [LEGACY_IDS.get(da, da) for da in ids]
            project["analyses"]["requested"] = ids
        if project["final_output"] == "FO004":
            project["reference"]["gff"] = "inputs/reference/Pf3D7.gff"
    parsed["config"] = project
    return {
        "values": parsed, "schema": Bundle.model_json_schema(), "revision": revision(root),
        "installation": read_json(installation) if installation.is_file() else None,
        "input_folders_created": (root / "inputs").is_dir(),
        "field_guide": read_json(Path(__file__).with_name("gui") / "field-guide.json"),
        "analyses": read_yaml(root / "config/registry/downstream_analyses.yaml")["analyses"],
    }


def check(root: Path, values: dict[str, Any]) -> list[str]:
    bundle = validate_bundle(root, values)
    if bundle.config.reference.release_label.strip() in {"", "SET_REFERENCE_RELEASE"}:
        raise ValueError("Enter the exact reference release label")
    return required_dirs(bundle.config)


def save(root: Path, values: dict[str, Any], expected: str, replace: bool) -> dict[str, Any]:
    dirs = check(root, values)
    with workspace_lock(root):
        current = revision(root)
        if current != expected:
            raise ValueError("Configuration changed since this form loaded. Reload the page.")
        # Check every destination before creating any files or directories.
        destinations = {name: within(root, rel) for name, rel in USER_FILES.items()}
        input_dirs = [within(root, rel) for rel in dirs]
        for path in input_dirs:
            if path.exists() and not path.is_dir():
                raise ValueError(f"Required input directory is a file: {path}")
            if any(parent.exists() and not parent.is_dir() for parent in path.parents):
                raise ValueError(f"Input directory has a non-directory parent: {path}")
        changes = {}
        originals: dict[str, str] = {}
        for name, path in destinations.items():
            if path.exists():
                old = read_yaml(path)
                if old == values[name]:
                    continue
                if not replace:
                    raise ValueError("Existing configuration differs. Select Save with backup to update it.")
                originals[name] = path.read_text()
            changes[name] = yaml.safe_dump(values[name], sort_keys=False)
        backups = []
        for name, content in originals.items():
            backup = within(root, f"config/backups/{current}/{name}.yaml")
            if not backup.exists():
                atomic_write_text(backup, content)
            backups.append(str(backup.relative_to(root)))
        create_input_dirs(root, dirs)
        written = []
        try:
            for name, content in changes.items():
                atomic_write_text(destinations[name], content)
                written.append(name)
        except OSError:
            for name in written:
                if name in originals:
                    atomic_write_text(destinations[name], originals[name])
                else:
                    destinations[name].unlink()
            raise
        return {"directories": dirs, "saved": [USER_FILES[name] for name in changes],
                "backups": backups, "revision": revision(root)}


PREFLIGHT_STAGES = ("validate", "prepare", "verify-prepared", "ready")


def run_preflight(root: Path) -> dict[str, Any]:
    stages = []
    env = os.environ.copy()
    env.pop("PFVARX_JSON_OUT", None)
    for stage in PREFLIGHT_STAGES:
        result = subprocess.run(
            [sys.executable, "-m", "pfvarx.cli", stage],
            cwd=root, env=env, capture_output=True, text=True, check=False,
        )
        stages.append({
            "command": stage, "exit_code": result.returncode,
            "stdout": result.stdout, "stderr": result.stderr,
        })
        if result.returncode:
            break
    return {"stages": stages, "passed": len(stages) == len(PREFLIGHT_STAGES) and not stages[-1]["exit_code"]}


class Handler(BaseHTTPRequestHandler):
    def __init__(self, *args: Any, root: Path, token: str, **kwargs: Any):
        self.root = root
        self.prefix = "/" + token + "/"
        super().__init__(*args, **kwargs)

    def respond(self, data: Any, status: int = 200, media: str = "application/json") -> None:
        body = data if isinstance(data, bytes) else json.dumps(data).encode()
        self.send_response(status)
        self.send_header("Content-Type", media)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'self'; frame-ancestors 'none'; base-uri 'none'")
        self.end_headers()
        self.wfile.write(body)

    def route(self) -> tuple[str, dict[str, list[str]]]:
        url = urlsplit(self.path)
        if not url.path.startswith(self.prefix):
            raise ValueError("Invalid configuration session")
        return url.path.removeprefix(self.prefix), parse_qs(url.query)

    def log_message(self, format: str, *args: Any) -> None:
        pass  # Keep the session token out of request logs.

    def do_GET(self) -> None:
        try:
            name, query = self.route()
            assets = {"": ("index.html", "text/html"), "app.js": ("app.js", "text/javascript"),
                      "style.css": ("style.css", "text/css")}
            if name in assets:
                filename, media = assets[name]
                self.respond(Path(__file__).with_name("gui").joinpath(filename).read_bytes(), media=media)
            elif name == "state":
                self.respond(editor_state(self.root, query.get("purpose", ["analysis"])[0]))
            elif name == "browse":
                root = within(self.root, "inputs")
                path = within(root, query.get("path", ["."])[0])
                entries = sorted(p for p in path.iterdir() if not p.name.startswith(".") and not p.is_symlink()) if path.is_dir() else []
                self.respond([{"path": str(p.relative_to(root)), "directory": p.is_dir()} for p in entries[:500]])
            else:
                self.respond({"error": "Unknown resource"}, 404)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.respond({"error": str(exc)}, 400)

    def do_POST(self) -> None:
        try:
            name, _ = self.route()
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= 262144 or self.headers.get_content_type() != "application/json":
                raise ValueError("Expected a configuration JSON request")
            payload = json.loads(self.rfile.read(size))
            if name == "preflight":
                self.respond(run_preflight(self.root))
                return
            values = payload["values"]
            if name == "check":
                self.respond({"directories": check(self.root, values)})
            elif name == "save":
                self.respond(save(self.root, values, payload["revision"], payload.get("replace") is True))
            else:
                self.respond({"error": "Unknown resource"}, 404)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.respond({"error": str(exc)}, 400)


def serve(root: Path, port: int, open_browser: bool) -> None:
    token = secrets.token_urlsafe(24)
    with HTTPServer(("127.0.0.1", port), partial(Handler, root=root, token=token)) as server:
        url = f"http://127.0.0.1:{server.server_port}/{token}/"
        print(f"Configuration editor: {url}\nPress Ctrl+C to stop.", flush=True)
        if open_browser:
            webbrowser.open(url)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
