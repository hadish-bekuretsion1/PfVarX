"""The SOP manual's appendices E, F and G are copies of guides; they must not drift apart."""

import runpy
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SYNC = runpy.run_path(str(REPO / "scripts/sync_sop.py"))


def test_sop_appendices_match_the_guides_they_mirror() -> None:
    stale = SYNC["stale"](REPO)
    assert not stale, f"{stale} differ from docs/caller.md, configuration.md or sequencing.md; run python scripts/sync_sop.py"


def test_every_mirrored_appendix_exists_in_both_formats() -> None:
    md = (REPO / "docs/sop/PfVarX_SOP_Manual.md").read_text()
    html = (REPO / "docs/sop/PfVarX_SOP_Manual.html").read_text()
    for letter, guide in SYNC["MIRRORS"].items():
        title = (REPO / guide).read_text().splitlines()[0][2:]
        assert f"# Appendix {letter} — {title}" in md
        assert f'<h1 id="{SYNC["slug"](f"Appendix {letter} — {title}")}">' in html
