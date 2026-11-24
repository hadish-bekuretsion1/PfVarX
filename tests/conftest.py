"""Shared fixtures: a temporary workspace that mirrors the repository layout."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture()
def workspace(tmp_path: Path) -> Path:
    for item in ("config", "workflow", "pixi.toml", "pyproject.toml"):
        src = REPO / item
        if src.is_dir():
            shutil.copytree(src, tmp_path / item)
        else:
            shutil.copyfile(src, tmp_path / item)
    for private in ("config.yaml", "quality.yaml", "compute.yaml", "profiles/resources.yaml"):
        (tmp_path / "config" / private).unlink(missing_ok=True)
    return tmp_path


@pytest.fixture(autouse=True)
def _private_history(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    """Tests never read or write this computer's resource history, and never open a browser."""
    monkeypatch.setenv("PFVARX_HISTORY_DIR", str(tmp_path_factory.mktemp("history")))
    monkeypatch.setenv("PFVARX_NO_BROWSER", "1")
