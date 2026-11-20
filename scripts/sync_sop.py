"""Keep the SOP manual's mirrored appendices identical to the guides they copy.

Appendices E, F and G of ``docs/sop/PfVarX_SOP_Manual.md`` (and its ``.html``) are copies of
``docs/caller.md``, ``docs/configuration.md`` and ``docs/sequencing.md``. This rewrites each
appendix from its guide: the guide's title becomes the appendix title, links to another
mirrored guide point to its appendix, other relative links gain ``../``, and the HTML is
rendered with markdown-it (CommonMark with tables), headings carrying the same ids as the
rest of the manual.

    pixi run python scripts/sync_sop.py           # rewrite the appendices
    pixi run python scripts/sync_sop.py --check   # report only; exit 1 if out of date

``tests/unit/test_sop_mirrors.py`` runs the check, so a guide cannot change without the SOP.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOP_MD = ROOT / "docs/sop/PfVarX_SOP_Manual.md"
SOP_HTML = ROOT / "docs/sop/PfVarX_SOP_Manual.html"
MIRRORS = {"E": "docs/caller.md", "F": "docs/configuration.md", "G": "docs/sequencing.md"}


def slug(text: str) -> str:
    text = re.sub(r"<[^>]+>", "", text).lower()
    text = re.sub(r"[^a-z0-9\s-]", "", text)
    return re.sub(r"[\s-]+", "-", text).strip("-")


def appendix(letter: str, root: Path) -> tuple[str, str]:
    """(title, markdown body) of the appendix that mirrors one guide."""
    lines = (root / MIRRORS[letter]).read_text().splitlines()
    assert lines[0].startswith("# "), f"{MIRRORS[letter]} must start with a level-1 title"
    title = f"Appendix {letter} — {lines[0][2:].strip()}"
    body = "\n".join(lines[1:]).strip("\n") + "\n"
    anchors = {Path(path).name: slug(f"Appendix {other} — " + (root / path).read_text().splitlines()[0][2:])
               for other, path in MIRRORS.items()}

    def link(match: re.Match[str]) -> str:
        target = match.group(2)
        if re.match(r"[a-z]+:|#|/", target):
            return match.group(0)
        name, _, fragment = target.partition("#")
        if name in anchors:
            return f"]({'#' + (fragment or anchors[name])})"
        return f"](../{target})"

    return title, re.sub(r"\]\(()([^)\s]+)\)", link, body)


def render(markdown: str) -> str:
    from markdown_it import MarkdownIt

    html = MarkdownIt("commonmark", {"html": True}).enable("table").render(markdown)

    def heading(match: re.Match[str]) -> str:
        return f'<h{match.group(1)} id="{slug(match.group(2))}">{match.group(2)}</h{match.group(1)}>'

    return re.sub(r"<h([1-6])>(.*?)</h\1>", heading, html)


def next_heading(md: str, start: int) -> int:
    """Offset of the next level-1 heading after *start*, ignoring '# ' lines inside code fences."""
    offset, fenced = start, False
    for line in md[start:].splitlines(keepends=True):
        if line.startswith("```"):
            fenced = not fenced
        elif not fenced and line.startswith("# "):
            return offset
        offset += len(line)
    return len(md)


def updated(root: Path) -> dict[Path, str]:
    md, html = (root / SOP_MD.relative_to(ROOT)).read_text(), (root / SOP_HTML.relative_to(ROOT)).read_text()
    for letter in MIRRORS:
        title, body = appendix(letter, root)
        section = f"# {title}\n\n{body}\n"
        md_start = re.search(rf"(?m)^# Appendix {letter} — .*$", md)
        assert md_start, f"SOP markdown has no Appendix {letter}"
        end = next_heading(md, md_start.end())
        md = md[:md_start.start()] + section + md[end:]
        html_start = re.search(rf'<h1 id="appendix-{letter.lower()}-[^"]*">', html)
        assert html_start, f"SOP HTML has no Appendix {letter}"
        html_end = re.compile(r"<h1 id=").search(html, html_start.end())
        tail = html[html_end.start():] if html_end else html[html.index("</main>"):]
        html = html[:html_start.start()] + render(section) + tail
    return {root / SOP_MD.relative_to(ROOT): md, root / SOP_HTML.relative_to(ROOT): html}


def stale(root: Path = ROOT) -> list[str]:
    return [str(path.relative_to(root)) for path, text in updated(root).items() if path.read_text() != text]


def main() -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--check", action="store_true", help="report only; exit 1 if an appendix is out of date")
    args = parser.parse_args()
    if args.check:
        out = stale()
        print("SOP appendices out of date: " + ", ".join(out) if out else "SOP appendices are current")
        return 1 if out else 0
    changed = [path for path, text in updated(ROOT).items() if path.read_text() != text]
    for path in changed:
        path.write_text(updated(ROOT)[path])
    print(f"updated {len(changed)} file(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
