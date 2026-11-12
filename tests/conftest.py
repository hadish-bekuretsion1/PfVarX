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
