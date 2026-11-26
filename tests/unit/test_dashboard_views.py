"""The dashboard's analysis tabs must build their charts in the function that renders them.

``views.js`` has two per-analysis functions: ``daKpis`` returns the headline numbers and
returns early for each analysis, while ``daViews`` builds the charts and tables. A visual
block written inside ``daKpis`` is unreachable, so the tab silently shows only its manifest
panel. That happened to DA008 and is checked here for every analysis, because no Python test
executes the dashboard JavaScript.
"""

from pathlib import Path

import pytest

from pfvarx.constants import DOWNSTREAM_IDS

VIEWS = Path(__file__).resolve().parents[2] / "src" / "pfvarx" / "reporting" / "views.js"
SOURCE = VIEWS.read_text(encoding="utf-8")
BUILDERS = ("visualGroup(", "recordsTable(", "chartPanel(", "out.push(")


def body(name: str) -> str:
    """Source of a top-level JavaScript function, by brace matching."""
    start = SOURCE.index(f"function {name}(")
    depth = 0
    for i in range(SOURCE.index("{", start), len(SOURCE)):
        if SOURCE[i] == "{":
            depth += 1
        elif SOURCE[i] == "}":
            depth -= 1
            if depth == 0:
                return SOURCE[start:i + 1]
    raise AssertionError(f"unbalanced braces in {name}")


@pytest.mark.parametrize("builder", BUILDERS)
def test_kpi_function_builds_no_views(builder: str) -> None:
    assert builder not in body("daKpis"), (
        f"daKpis calls {builder}: it returns early per analysis, so the block never runs. "
        "Move it into daViews.")


def test_every_analysis_visual_group_is_in_the_view_function() -> None:
    views = body("daViews")
    for da in DOWNSTREAM_IDS:
        groups = [part for part in SOURCE.split('id:"') if part.startswith(f"{da}-")]
        for group in groups:
            marker = 'id:"' + group.split('"')[0] + '"'
            assert marker in views, f"{marker} is built outside daViews, so the {da} tab cannot show it"


def test_each_analysis_with_results_has_a_tab_and_a_view_block() -> None:
    kpis, views = body("daKpis"), body("daViews")
    for da in DOWNSTREAM_IDS:
        assert f'{da}:{{tab:' in SOURCE.replace(" ", ""), f"{da} has no tab in DA_META"
        assert f'id==="{da}"' in kpis, f"{da} has no headline numbers in daKpis"
        assert f'id==="{da}"' in views, f"{da} has no charts or tables in daViews"
