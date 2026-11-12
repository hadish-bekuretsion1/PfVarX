"""pytest wrapper: runs the browser checks when PFVARX_TEST_REPORT points at a built report."""

import os
from pathlib import Path

import pytest

playwright = pytest.importorskip("playwright")


@pytest.mark.skipif(not os.environ.get("PFVARX_TEST_REPORT"), reason="set PFVARX_TEST_REPORT=<run>/results/report")
def test_dashboard(tmp_path: Path) -> None:
    from check_dashboard import main  # type: ignore[import-not-found]

    r = main(os.environ["PFVARX_TEST_REPORT"], str(tmp_path))
    assert r["passed"], r
