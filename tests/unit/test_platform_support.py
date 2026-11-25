"""Native targets must select their own locks and retain resource safeguards."""

import runpy
import subprocess
from pathlib import Path

import pytest
import yaml

from pfvarx import installation
from pfvarx.planning import plan

ROOT = Path(__file__).resolve().parents[2]
BOOTSTRAP = runpy.run_path(str(ROOT / "bootstrap.py"))


@pytest.mark.parametrize("system,machine,target", [
    ("Linux", "x86_64", "linux-64"),
    ("Linux", "aarch64", "linux-aarch64"),
    ("Linux", "arm64", "linux-aarch64"),
    ("Darwin", "x86_64", "osx-64"),
    ("Darwin", "arm64", "osx-arm64"),
])
def test_native_platform_identifiers(system: str, machine: str, target: str) -> None:
    assert BOOTSTRAP["host_platform"](system, machine) == target


@pytest.mark.parametrize("system,machine", [("Windows", "AMD64"), ("Linux", "armv7l")])
def test_unsupported_hosts_have_an_actionable_error(system: str, machine: str) -> None:
    with pytest.raises(ValueError, match="WSL2"):
        BOOTSTRAP["host_platform"](system, machine)


def test_tool_versions_are_selected_per_platform(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    lock = {"packages": [], "environments": {"core": {"packages": {}}}}
    for target, version in [("linux-64", "1.21"), ("osx-arm64", "1.22")]:
        url = f"https://example.invalid/{target}/samtools-{version}-build.conda"
        lock["packages"].append({"conda": url})
        lock["environments"]["core"]["packages"][target] = [{"conda": url}]
    (tmp_path / "pixi.lock").write_text(yaml.safe_dump(lock))
    monkeypatch.setattr(installation, "host_platform", lambda root: "linux-64")
    assert installation.locked_versions(tmp_path)["core"]["samtools"] == "1.21"
    monkeypatch.setattr(installation, "host_platform", lambda root: "osx-arm64")
    assert installation.locked_versions(tmp_path)["core"]["samtools"] == "1.22"


def test_macos_total_memory_is_detected(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(plan.sys, "platform", "darwin")

    def sysctl(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        assert args == ["sysctl", "-n", "hw.memsize"]
        return subprocess.CompletedProcess(args, 0, stdout=str(16 * 1024**3), stderr="")

    monkeypatch.setattr(plan.subprocess, "run", sysctl)
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    assert plan._host_mem_mb() == 16384


def test_unreadable_macos_memory_is_unknown(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(plan.sys, "platform", "darwin")

    def unavailable(*args: object, **kwargs: object) -> None:
        raise OSError("sysctl unavailable")

    monkeypatch.setattr(plan.subprocess, "run", unavailable)
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    assert plan._host_mem_mb() is None


def test_every_environment_is_locked_for_every_native_target() -> None:
    import tomllib

    targets = {"linux-64", "linux-aarch64", "osx-64", "osx-arm64"}
    manifest = tomllib.loads((ROOT / "pixi.toml").read_text())
    assert set(manifest["workspace"]["platforms"]) == targets
    lock = yaml.load((ROOT / "pixi.lock").read_text(), Loader=yaml.CSafeLoader)
    assert {p["name"] for p in lock["platforms"]} == targets
    for environment in manifest["environments"]:
        assert set(lock["environments"][environment]["packages"]) == targets, environment


def test_missing_platform_lock_is_an_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "pixi.lock").write_text("environments: {core: {packages: {linux-64: []}}}\n")
    monkeypatch.setattr(installation, "host_platform", lambda root: "osx-arm64")
    with pytest.raises(ValueError, match="no lock for osx-arm64"):
        installation.locked_versions(tmp_path)


def test_spawn_workers_match_serial_pairhmm() -> None:
    """macOS uses spawn: worker functions must be importable without fork-inherited state."""
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    from pfvarx.caller.pairhmm import support

    hap = "ACGTTAGCATTTTTGACCAGTACGATCGA"
    alt = hap[:14] + "T" + hap[14:]
    read = alt[3:25]
    args = (read, [35] * len(read), hap, alt, 1e-3, 0.1)
    expected = support(*args)
    with ProcessPoolExecutor(2, mp_context=multiprocessing.get_context("spawn")) as pool:
        assert pool.submit(support, *args).result(timeout=30) == expected == 1
