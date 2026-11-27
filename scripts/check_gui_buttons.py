"""Check that every button of the configuration GUI does what it says (software check).

It starts the real GUI server (``./pfvarx configure --gui``) on a new temporary workspace,
sends each button's request exactly as the page does, and then checks the result on disk:
which folders and files exist, and which do not exist yet. Synthetic test reads come from
``pfvarx fixture``; nothing is downloaded and no browser is needed.

    pixi run python scripts/check_gui_buttons.py            # keeps nothing
    pixi run python scripts/check_gui_buttons.py --keep     # keeps the workspaces to inspect

Each line reads PASS or FAIL with the button, what was checked and what was found. The exit
code is 0 only when every check passes. It covers the layout and wiring of the buttons, not
the visual design; tests/dashboard/check_configure.py clicks them in a real browser.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
RESULTS: list[tuple[bool, str, str]] = []


def record(ok: bool, button: str, detail: str) -> bool:
    RESULTS.append((ok, button, detail))
    print(f"{'PASS' if ok else 'FAIL'}  {button:<46} {detail}", flush=True)
    return ok


def folders(ws: Path) -> list[str]:
    """Top-level workspace folders that steps create (config/ and the source are excluded)."""
    return sorted(p.name for p in ws.iterdir() if p.is_dir() and p.name not in {"config", "workflow", "src"})


class Gui:
    """The GUI server for one workspace, called the way the page calls it."""

    def __init__(self, ws: Path) -> None:
        env = dict(os.environ, PYTHONPATH=str(REPO / "src"), PYTHONDONTWRITEBYTECODE="1")
        self.process = subprocess.Popen([sys.executable, "-m", "pfvarx.cli", "configure", "--gui", "--no-browser"],
                                        cwd=ws, env=env, stdout=subprocess.PIPE, text=True)
        assert self.process.stdout is not None
        self.url = self.process.stdout.readline().strip().split()[-1]

    def get(self, path: str) -> Any:
        with urllib.request.urlopen(self.url + path, timeout=600) as response:
            return json.load(response)

    def post(self, path: str, payload: Any = None, data: bytes | None = None) -> Any:
        body = data if data is not None else json.dumps(payload).encode()
        kind = "application/octet-stream" if data is not None else "application/json"
        request = urllib.request.Request(self.url + path, data=body, headers={"Content-Type": kind}, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=1800) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            return {"error": json.load(error).get("error", str(error))}

    def stop(self) -> None:
        self.process.terminate()
        self.process.wait(timeout=10)


def new_workspace(base: Path, name: str) -> Path:
    ws = base / name
    ws.mkdir()
    for item in ("config", "workflow", "src"):
        shutil.copytree(REPO / item, ws / item, ignore=shutil.ignore_patterns("__pycache__"))
    for item in ("pixi.toml", "pixi.lock", "pyproject.toml"):  # a real run needs the locked environment
        shutil.copyfile(REPO / item, ws / item)
    for private in ("config.yaml", "quality.yaml", "compute.yaml", "profiles/resources.yaml"):
        (ws / "config" / private).unlink(missing_ok=True)
    return ws


def synthetic_inputs(base: Path) -> Path:
    """A donor workspace whose inputs/ holds synthetic reads and reference (software test data)."""
    donor = new_workspace(base, "donor")
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    for args in (["config-init", "--purpose", "analysis"], ["fixture", "--design", "wgs", "--depth", "20"]):
        subprocess.run([sys.executable, "-m", "pfvarx.cli", *args], cwd=donor, env=env, check=True,
                       capture_output=True)
    return donor


def stages(result: dict[str, Any]) -> str:
    return ", ".join(f"{s['command']}={'ok' if s['exit_code'] == 0 else 'failed'}" for s in result.get("stages", []))


def save_form(gui: Gui, **config: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load the form, change settings as a user would, and click Create Configuration and Input Directories."""
    state = gui.get("state?purpose=analysis")
    values = state["values"]
    values["config"]["reference"]["release_label"] = "synthetic"
    values["config"]["final_output"] = "FO001"
    for key, value in config.items():
        values["config"][key] = value
    return values, gui.post("save", {"values": values, "revision": state["revision"], "replace": False})


def run_study_metadata_no(base: Path, donor: Path) -> None:
    ws = new_workspace(base, "study-metadata-no")
    gui = Gui(ws)
    try:
        record(not (ws / "inputs").exists(), "Open the GUI", "no input folders before saving")
        values, saved = save_form(gui)
        expected = ["inputs/fastq", "inputs/reference"]
        made = sorted(str(p.relative_to(ws)) for p in (ws / "inputs").iterdir()) if (ws / "inputs").is_dir() else []
        record("error" not in saved and made == expected, "① Create Configuration and Input Directories",
               f"input folders {made} (expected {expected})")
        record(folders(ws) == ["inputs"], "① ... nothing else created yet", f"workspace folders {folders(ws)}")
        record((ws / "config/config.yaml").is_file(), "① ... configuration saved", "config/config.yaml written")

        sent = []
        for path in sorted((donor / "inputs/fastq").iterdir()):
            sent.append(gui.post(f"upload-raw?filename={path.name}&role=raw", data=path.read_bytes()))
        for name in ("Pf3D7.fasta", "core_genome.bed"):
            path = donor / "inputs/reference" / name
            sent.append(gui.post(f"upload-raw?filename={name}&role=reference", data=path.read_bytes()))
        copied = all((ws / "inputs/fastq" / p.name).read_bytes() == p.read_bytes()
                     for p in (donor / "inputs/fastq").iterdir())
        record(copied and not any("error" in s for s in sent), "② Upload selected files",
               f"{len(sent)} files uploaded, contents identical")
        refused = gui.post("upload-raw?filename=study.tsv&role=metadata", data=b"sample_id\tcountry\nPF001\tX\n")
        record("error" in refused, "② Upload to metadata (Study Metadata: No)", f"refused: {refused.get('error')}")

        result = gui.post("validate-inputs", {})
        record(result.get("passed") is True, "③ Validate inputs", stages(result))

        result = gui.post("prepare-inputs", {})
        report = ws / "metadata/inputs_preparation.json"
        merge = json.loads(report.read_text())["study_metadata"]["merge"] if report.is_file() else None
        index = (ws / "prepared/reference/combined.fasta.bwt.2bit.64").is_file()
        record(result.get("passed") is True and (ws / "prepared/samples.tsv").is_file() and index,
               "④ Prepare inputs (Inputs Preparation)", f"{stages(result)}; sample metadata and bwa-mem2 index built")
        record(merge == "NOT_USED", "④ ... Study Metadata: No", f"merge = {merge} (continued without study metadata)")

        result = gui.post("preflight", {})
        record(result.get("passed") is True, "⑤ Run readiness checks", stages(result))

        result = gui.post("workflow-check", {})
        record(result.get("passed") is True and not (ws / "runs").exists(), "⑥ Check workflow",
               f"{stages(result)}; no run created")

        gui.post("workflow-run", {})
        status: dict[str, Any] = {}
        for _ in range(900):
            status = gui.get("workflow-status")
            if not status["running"] and status["exit_code"] is not None:
                break
            time.sleep(2)
        record(status.get("latest_status") == "COMPLETED", "⑥ Run workflow",
               f"run {status.get('latest_run')}: {status.get('latest_status')}")
        record("runs" in folders(ws), "⑥ ... run folders created on demand", f"workspace folders {folders(ws)}")
    finally:
        gui.stop()


def run_study_metadata_yes(base: Path, donor: Path) -> None:
    ws = new_workspace(base, "study-metadata-yes")
    gui = Gui(ws)
    try:
        values, saved = save_form(gui, study_metadata=True)
        record("error" not in saved and (ws / "inputs/metadata").is_dir(), "① Create ... (Study Metadata: Yes)",
               f"input folders {sorted(p.name for p in (ws / 'inputs').iterdir())}")
        for path in sorted((donor / "inputs/fastq").iterdir()):
            gui.post(f"upload-raw?filename={path.name}&role=raw", data=path.read_bytes())
        for name in ("Pf3D7.fasta", "core_genome.bed"):
            gui.post(f"upload-raw?filename={name}&role=reference", data=(donor / "inputs/reference" / name).read_bytes())
        table = "sample_id\tcountry\tcollection_date\nPF001\tEthiopia\t2025-01-15\n"
        sent = gui.post("upload-raw?filename=study_metadata.tsv&role=metadata", data=table.encode())
        record("error" not in sent and (ws / "inputs/metadata/study_metadata.tsv").is_file(),
               "② Upload to metadata (Study Metadata: Yes)", sent.get("path", sent.get("error", "")))
        result = gui.post("prepare-inputs", {})
        merged = json.loads((ws / "metadata/inputs_preparation.json").read_text())["study_metadata"]
        rows = (ws / "prepared/samples.tsv").read_text().splitlines()
        first = dict(zip(rows[0].split("\t"), rows[1].split("\t"), strict=False))
        record(result.get("passed") is True and merged["merge"] == "MERGED" and first.get("country") == "Ethiopia",
               "④ Prepare inputs ... study metadata merged",
               f"merge = {merged['merge']}, PF001 country = {first.get('country')!r}")
        state = gui.get("state?purpose=analysis")
        values = state["values"]
        values["config"].update(assay_design="targeted_umi", assay={**values["config"]["assay"],
                                "library_catalog": "inputs/metadata/library_catalog.yaml"})
        values["config"]["reference"].update(host_removal=True, host_fasta="inputs/reference/human.fasta")
        refused = gui.post("save", {"values": values, "revision": state["revision"], "replace": True})
        record("eligible only for assay_design wgs" in refused.get("error", ""),
               "① Create ... Human DNA removal True + targeted", "refused: wgs design only")
    finally:
        gui.stop()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--keep", action="store_true", help="keep the temporary workspaces and print where")
    args = parser.parse_args()
    base = Path(tempfile.mkdtemp(prefix="pfvarx-gui-buttons-"))
    try:
        donor = synthetic_inputs(base)
        run_study_metadata_no(base, donor)
        run_study_metadata_yes(base, donor)
    finally:
        failed = [r for r in RESULTS if not r[0]]
        print(f"\n{len(RESULTS) - len(failed)} of {len(RESULTS)} checks passed.")
        if args.keep:
            print(f"Workspaces kept in {base}")
        else:
            shutil.rmtree(base, ignore_errors=True)
    return 1 if failed or not RESULTS else 0


if __name__ == "__main__":
    raise SystemExit(main())
