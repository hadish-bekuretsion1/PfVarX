"""One authoritative pipeline version; every current-release reference must agree with it.

Only statements about the *current* PfVarX release are checked. Dependency versions,
schema/parameter/policy/catalogue versions and earlier changelog sections are separate
version series and are deliberately not compared.
"""

from __future__ import annotations

import re
import subprocess
import sys
import tomllib
from pathlib import Path

import pfvarx

REPO = Path(__file__).resolve().parents[2]
VERSION = pfvarx.__version__
SEMVER = r"\d+\.\d+\.\d+(?:-[0-9A-Za-z.-]+)?"

# Documents that describe the current release, with the phrases that state its version.
CURRENT_RELEASE_DOCS = (
    "README.md", "docs/README.html", "docs/status.md",
    "docs/configuration-guide.md", "docs/configuration-guide.html",
    "docs/sop/PfVarX_SOP_Manual.md", "docs/sop/PfVarX_SOP_Manual.html",
)
PHRASES = (
    rf"PfVarX ({SEMVER})",           # "PfVarX 0.0.0 is ...", page titles
    rf"\bVersion ({SEMVER})",        # configuration guide
    rf"finalized (?:local )?({SEMVER})",
    rf"pfvarx-({SEMVER})-source",    # archive names must be derived, never written out
)


def test_version_is_semantic() -> None:
    assert re.fullmatch(SEMVER, VERSION), VERSION


def test_package_metadata_reads_the_single_source() -> None:
    project = tomllib.loads((REPO / "pyproject.toml").read_text())
    assert "version" not in project["project"], "pyproject.toml must not declare its own version"
    assert "version" in project["project"]["dynamic"]
    assert project["tool"]["hatch"]["version"]["path"] == "src/pfvarx/__init__.py"
    workspace = tomllib.loads((REPO / "pixi.toml").read_text())["workspace"]
    assert "version" not in workspace, "pixi.toml must not declare a second pipeline version"


def test_cli_reports_the_authoritative_version() -> None:
    out = subprocess.run([sys.executable, "-m", "pfvarx.cli", "version"], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == VERSION


def test_source_archive_name_uses_the_authoritative_version() -> None:
    import runpy

    packager = runpy.run_path(str(REPO / "scripts/package.py"))
    assert packager["pipeline_version"](REPO) == VERSION
    ci = (REPO / ".github/workflows/ci.yaml").read_text()
    assert not re.search(rf"pfvarx-{SEMVER}-source", ci), "CI must derive the archive name"


def test_current_release_documents_agree() -> None:
    stale = []
    for name in CURRENT_RELEASE_DOCS:
        text = (REPO / name).read_text(encoding="utf-8")
        for phrase in PHRASES:
            stale += [f"{name}: {m.group(0)}" for m in re.finditer(phrase, text) if m.group(1) != VERSION]
    assert not stale, "current-release references disagree with __version__: " + "; ".join(stale)


def test_no_other_document_names_a_current_version() -> None:
    """Other docs must not restate the release (the earlier 0.4.3 contradiction)."""
    offenders = []
    for path in sorted((REPO / "docs").rglob("*.md")):
        rel = path.relative_to(REPO).as_posix()
        if rel in CURRENT_RELEASE_DOCS:
            continue
        text = path.read_text(encoding="utf-8")
        offenders += [f"{rel}: {m.group(0)}" for m in re.finditer(
            rf"(?:CLI|package|pipeline|release)[/A-Za-z]* version is ({SEMVER})|PfVarX ({SEMVER})|— ({SEMVER})$",
            text, re.MULTILINE)]
    assert not offenders, offenders


def test_changelog_has_one_section_for_the_current_version() -> None:
    headings = re.findall(rf"^## ({SEMVER})\b", (REPO / "CHANGELOG.md").read_text(), re.MULTILINE)
    assert headings and headings[0] == VERSION, f"first changelog section must be {VERSION}"
    assert headings.count(VERSION) == 1, f"{VERSION} appears in {headings.count(VERSION)} changelog headings"
