"""Software archives must exclude study data, user configuration and environments."""
import runpy
import shutil
import zipfile
from pathlib import Path
from uuid import uuid4

import pytest

ROOT = Path(__file__).resolve().parents[2]
PACKAGER = runpy.run_path(str(ROOT / "scripts/package.py"))


def _copy_source(destination: Path) -> Path:
    destination.mkdir()
    for path in PACKAGER["source_files"](ROOT):
        target = destination / path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    return destination


def test_archive_excludes_runtime_and_is_reproducible(tmp_path: Path) -> None:
    source = _copy_source(tmp_path / "source")
    sentinel = "PRIVATE_" + uuid4().hex
    for name in ("inputs/fastq/private.fastq.gz", "results/private.vcf.gz", "runs/private.json",
                 ".pixi/private.txt", "docs/generated-dashboard.json", "config/config.yaml", "config/compute.yaml",
                 "config/quality.yaml", "config/profiles/resources.yaml", "src/pfvarx/__pycache__/private.pyc"):
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(sentinel)
    first = PACKAGER["build_archive"](source, tmp_path / "first")
    second = PACKAGER["build_archive"](source, tmp_path / "second")
    assert first["sha256"] == second["sha256"]
    with zipfile.ZipFile(first["archive"]) as bundle:
        assert "pfvarx/src/pfvarx/cli.py" in bundle.namelist()
        for asset in ("dashboard.css", "core.js", "exports.js", "tables.js", "charts.js", "views.js", "dashboard.schema.json"):
            assert "pfvarx/src/pfvarx/reporting/" + asset in bundle.namelist()
        assert "pfvarx/docs/reporting/synthetic-demo.json" not in bundle.namelist()
        assert "pfvarx/docs/reporting/README.md" in bundle.namelist()
        assert "pfvarx/docs/analysis/README.md" in bundle.namelist()
        assert "pfvarx/docs/benchmarking/README.md" in bundle.namelist()
        assert all(sentinel.encode() not in bundle.read(name) for name in bundle.namelist())
        assert not any(name.startswith(("pfvarx/inputs/", "pfvarx/runs/", "pfvarx/results/", "pfvarx/.pixi/"))
                       for name in bundle.namelist())
    with pytest.raises(ValueError, match="not pristine"):
        PACKAGER["check_pristine"](source)


def test_pristine_archive_roundtrip_and_symlink_rejection(tmp_path: Path) -> None:
    source = _copy_source(tmp_path / "source")
    assert PACKAGER["check_pristine"](source)["runtime_files"] == 0
    result = PACKAGER["build_archive"](source, tmp_path / "archive")
    with zipfile.ZipFile(result["archive"]) as bundle:
        bundle.extractall(tmp_path / "unpacked")
    assert PACKAGER["check_pristine"](tmp_path / "unpacked/pfvarx")["runtime_files"] == 0
    external = tmp_path / "private.txt"
    external.write_text("PRIVATE_STUDY_DATA")
    (source / "src/pfvarx/unsafe.py").symlink_to(external)
    with pytest.raises(ValueError, match="symlinks"):
        PACKAGER["build_archive"](source, tmp_path / "rejected")


def test_development_content_is_not_released(tmp_path: Path) -> None:
    """The development guide and other development/ material stay out of the archive."""
    assert (ROOT / "development/PfVarX_Development_Guide.html").is_file()
    source = _copy_source(tmp_path / "source")
    shutil.copytree(ROOT / "development", source / "development")
    result = PACKAGER["build_archive"](source, tmp_path / "archive")
    with zipfile.ZipFile(result["archive"]) as bundle:
        assert not any(name.startswith("pfvarx/development/") for name in bundle.namelist())
        assert not any("Development_Guide" in name for name in bundle.namelist())
