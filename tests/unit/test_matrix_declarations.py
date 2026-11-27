"""Each analysis declares whether it reads the cohort evidence matrix, and the code must agree.

Layer 5 skips loading the matrix for analyses that never read it, which is the difference
between gigabytes and nothing for a large cohort. A wrong declaration would be silent and
dangerous: the analysis would run against an empty matrix and report no evidence rather than
fail. So the declaration is checked against the source, following imports inside
``pfvarx.downstream``, since an analysis can reach the matrix through a helper.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from pfvarx.constants import DOWNSTREAM_IDS
from pfvarx.downstream.run import analysis_module, needs_matrix

DOWNSTREAM = Path(__file__).resolve().parents[2] / "src" / "pfvarx" / "downstream"
READS_MATRIX = ("ctx.matrix", ".matrix[")


def local_imports(module: str) -> set[str]:
    """Other ``pfvarx.downstream`` modules this one imports."""
    tree = ast.parse((DOWNSTREAM / f"{module}.py").read_text())
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("pfvarx.downstream."):
            found.add((node.module or "").rsplit(".", 1)[1])
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("pfvarx.downstream."):
                    found.add(alias.name.rsplit(".", 1)[1])
    return found - {"common", "run"}


def reaches_matrix(module: str, seen: set[str] | None = None) -> bool:
    """Whether the module, or anything it imports from this package, reads the matrix."""
    seen = seen if seen is not None else set()
    if module in seen:
        return False
    seen.add(module)
    path = DOWNSTREAM / f"{module}.py"
    if not path.is_file():
        return False
    if any(marker in path.read_text() for marker in READS_MATRIX):
        return True
    return any(reaches_matrix(imported, seen) for imported in local_imports(module))


@pytest.mark.parametrize("da", DOWNSTREAM_IDS)
def test_declaration_matches_the_code(da: str) -> None:
    module = analysis_module(da)
    declared = needs_matrix(da)
    actual = reaches_matrix(module)
    assert declared == actual, (
        f"{da} ({module}.py) declares USES_COHORT_MATRIX = {declared} but the code "
        f"{'reads' if actual else 'never reads'} the cohort matrix. A wrong declaration gives the "
        "analysis an empty matrix instead of failing.")


@pytest.mark.parametrize("da", DOWNSTREAM_IDS)
def test_every_analysis_declares_itself(da: str) -> None:
    source = (DOWNSTREAM / f"{analysis_module(da)}.py").read_text()
    assert "USES_COHORT_MATRIX" in source, f"{da} does not declare whether it reads the cohort matrix"


def test_the_analyses_that_work_from_alignments_skip_the_matrix() -> None:
    # The point of the declaration: these four never pay for the matrix.
    assert [da for da in DOWNSTREAM_IDS if not needs_matrix(da)] == ["DA002", "DA007", "DA008", "DA009"]
