"""Loopback-only configuration editor and CLI preflight launcher."""

from __future__ import annotations

import json
import os
import secrets
import subprocess
import sys
import tempfile
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
    load_bundle,
    migrate_project,
    read_yaml,
    sha256_text,
    validate_bundle,
)
from pfvarx.configuration.migrate import migrate_workspace, plan_v3
from pfvarx.configuration.models import Bundle
from pfvarx.constants import Design
from pfvarx.execution.retention import workspace_lock
from pfvarx.initialization.layout import check_release_label, create_input_dirs, required_dirs
from pfvarx.inputs.transfer import check_content, destination_for, download, folder_contents
from pfvarx.installation import stage_environment, stage_python
from pfvarx.installation import status as tools_status
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
    migration_needed = project.get("schema_version") == "pfvarx.config/3"
    if migration_needed:
        # Show the schema-4 form; saving migrates the catalog and manifest first (config-migrate).
        project = plan_v3(project)[0]
    installation = root / "config/installation.json"
    if installation.is_file() and not (root / USER_FILES["config"]).is_file() \
            and read_json(installation).get("profile") != "on-demand":
        meta = read_json(installation)
        project.update(final_output=meta["max_final_output"], input_type=meta.get("input_type", "FASTQ"))
        project["assay_design"] = Design.parse(meta["design"]).value
        if Design.parse(meta["design"]).targeted:
            project["assay"]["library_catalog"] = "inputs/metadata/library_catalog.yaml"
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
        "migration_needed": migration_needed,
        "installation": read_json(installation) if installation.is_file() else None,
        "input_folders_created": (root / "inputs").is_dir(),
        "field_guide": read_json(Path(__file__).with_name("gui") / "field-guide.json"),
        "analyses": read_yaml(root / "config/registry/downstream_analyses.yaml")["analyses"],
    }


def check(root: Path, values: dict[str, Any]) -> list[str]:
    bundle = validate_bundle(root, values)
    check_release_label(bundle.config)
    return required_dirs(bundle.config)


def save(root: Path, values: dict[str, Any], expected: str, replace: bool) -> dict[str, Any]:
    dirs = check(root, values)
    with workspace_lock(root):
        current = revision(root)
        if current != expected:
            raise ValueError("Configuration changed since this form loaded. Reload the page.")
        migration = None
        on_disk = within(root, USER_FILES["config"])
        if on_disk.is_file() and migrate_project(read_yaml(on_disk)).get("schema_version") == "pfvarx.config/3":
            # Same step as './pfvarx config-migrate --apply': library catalog and manifest columns.
            migration = migrate_workspace(root, apply=True)
            if migration["status"] == "BLOCKED":
                raise ValueError(migration["detail"])
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
                "backups": backups, "revision": revision(root), "migration": migration}


PREFLIGHT_STAGES = ("prepare-samples", "validate", "prepare", "verify-prepared", "ready")
# Inputs Preparation: sample sheet and optional metadata merge, then references from the input
# files and final output, verification and readiness (the same as ./pfvarx inputs-preparation).
INPUTS_PREPARATION_STAGES = ("prepare-samples", "prepare", "verify-prepared", "ready")


def run_preflight(root: Path, stages_to_run: tuple[str | tuple[str, ...], ...] = PREFLIGHT_STAGES) -> dict[str, Any]:
    """Run CLI stages in order and stop at the first failure; a stage may carry its own options."""
    stages = []
    env = os.environ.copy()
    env.pop("PFVARX_JSON_OUT", None)
    tools_env = stage_environment(root, env)  # after Install required tools: its bin first on PATH
    for stage in stages_to_run:
        args = [stage] if isinstance(stage, str) else list(stage)
        # Steps run in the installed tools environment (on demand: after Install required tools).
        python = sys.executable if args[0] == "install-tools" else stage_python(root, sys.executable)
        result = subprocess.run(
            [python, "-m", "pfvarx.cli", *args],
            cwd=root, env=env if python == sys.executable else tools_env, capture_output=True, text=True, check=False,
        )
        stages.append({
            "command": " ".join(args), "exit_code": result.returncode,
            "stdout": result.stdout, "stderr": result.stderr,
        })
        if result.returncode:
            break
    return {"stages": stages, "passed": len(stages) == len(stages_to_run) and not stages[-1]["exit_code"]}


raw_destination = destination_for  # one rule set for GUI uploads/downloads and the CLI


WORKFLOW_RESULT = "metadata/workflow_run.json"
WORKFLOW_LOG = "metadata/workflow_run.log"
_workflow: dict[str, Any] = {"process": None}


def start_workflow(root: Path) -> dict[str, Any]:
    """Start './pfvarx workflow-run' in the background; the page polls workflow_status."""
    running = _workflow["process"]
    if running is not None and running.poll() is None:
        raise ValueError("A workflow run started from this page is still in progress")
    result = within(root, WORKFLOW_RESULT)
    result.parent.mkdir(parents=True, exist_ok=True)
    result.unlink(missing_ok=True)
    env = os.environ.copy()
    env["PFVARX_JSON_OUT"] = str(result)
    env = stage_environment(root, env)
    with within(root, WORKFLOW_LOG).open("w") as log:
        _workflow["process"] = subprocess.Popen(
            [stage_python(root, sys.executable), "-m", "pfvarx.cli", "workflow-run"], cwd=root, env=env,
            stdout=subprocess.DEVNULL, stderr=log, start_new_session=True,
        )
    return {"started": True, "command": "./pfvarx workflow-run"}


def workflow_status(root: Path) -> dict[str, Any]:
    from pfvarx.runs.registry import latest

    process = _workflow["process"]
    running = process is not None and process.poll() is None
    result = root / WORKFLOW_RESULT
    run = latest(root) if (root / "runs").is_dir() else None
    status = read_json(run / "status.json") if run and (run / "status.json").is_file() else None
    return {
        "running": running,
        "exit_code": None if process is None or running else process.returncode,
        "result": read_json(result) if result.is_file() and not running else None,
        "latest_run": run.name if run else None,
        "latest_status": status.get("status") if status else None,
        "log": WORKFLOW_LOG,
    }


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
            elif name == "input-folders":
                cfg = load_bundle(self.root).config
                folders = [folder for folder in required_dirs(cfg) if (self.root / folder).is_dir()]
                self.respond({
                    "folders": folders,
                    "contents": folder_contents(self.root),
                    "annotation": cfg.reference.gff is not None,
                })
            elif name == "tools-status":
                self.respond(tools_status(self.root))
            elif name == "workflow-status":
                self.respond(workflow_status(self.root))
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
            name, query = self.route()
            if name == "upload-raw":
                size = int(self.headers.get("Content-Length", "0"))
                if size <= 0 or self.headers.get_content_type() != "application/octet-stream":
                    raise ValueError("Choose a nonempty raw data file")
                with workspace_lock(self.root):
                    destination = raw_destination(self.root, query.get("filename", [""])[0], query.get("role", ["raw"])[0])
                    with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".upload-") as staged:
                        remaining = size
                        while remaining:
                            block = self.rfile.read(min(remaining, 1024 * 1024))
                            if not block:
                                raise ValueError("Upload ended before the complete file arrived")
                            staged.write(block)
                            remaining -= len(block)
                        staged.flush()
                        check_content(self.root, Path(staged.name), destination)
                        os.link(staged.name, destination)
                self.respond({"path": str(destination.relative_to(self.root)), "bytes": size})
                return
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= 262144 or self.headers.get_content_type() != "application/json":
                raise ValueError("Expected a configuration JSON request")
            payload = json.loads(self.rfile.read(size))
            if name == "download-raw":
                with workspace_lock(self.root):
                    result = download(self.root, payload["url"], (payload.get("sha256") or "").strip() or None,
                                      payload.get("role", "raw"), payload["filename"])
                self.respond(result)
                return
            if name == "workflow-run":
                self.respond(start_workflow(self.root))
                return
            checks = {
                "preflight": PREFLIGHT_STAGES,
                "validate-inputs": ("validate",),
                "prepare-inputs": INPUTS_PREPARATION_STAGES,
                "workflow-check": ("dry-run",),
                "install-tools": ("install-tools",),
            }
            if name in checks:
                self.respond(run_preflight(self.root, checks[name]))
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
