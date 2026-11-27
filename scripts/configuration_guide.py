"""Keep the configuration guide's setting tables identical to the GUI field help.

``src/pfvarx/configuration/gui/field-guide.json`` is the single source for what each setting
means; the GUI shows it next to every field. This script rewrites the example and help cells
of the setting tables in ``docs/configuration-guide.md`` and ``.html`` from it, so terminal
users read exactly what GUI users read.

    python scripts/configuration_guide.py          # rewrite both guide files
    python scripts/configuration_guide.py --check  # exit 1 if either file is out of date
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIELD_GUIDE = ROOT / "src/pfvarx/configuration/gui/field-guide.json"
GUIDE_MD = ROOT / "docs/configuration-guide.md"
GUIDE_HTML = ROOT / "docs/configuration-guide.html"
MD_ROW = re.compile(r"^\| `(?P<key>[a-z_]+(?:\.[a-z_0-9]+)*)` \| .*? \| .* \|$", re.MULTILINE)
HTML_ROW = re.compile(r'<tr><th scope="row"><code>(?P<key>[a-z_]+(?:\.[a-z_0-9]+)*)</code></th>.*?</tr>')


def guide_key(doc_key: str) -> str:
    """Docs drop the 'config.' prefix of project settings; quality and compute keep theirs."""
    return doc_key if doc_key.split(".")[0] in ("quality", "compute") else "config." + doc_key


def render(root: Path = ROOT) -> dict[Path, str]:
    """Return the guide files as they should read, keyed by path."""
    guide = json.loads((root / FIELD_GUIDE.relative_to(ROOT)).read_text())

    def entry(key: str) -> dict[str, str]:
        if guide_key(key) not in guide:
            raise KeyError(f"{key}: no entry in field-guide.json")
        return guide[guide_key(key)]

    def md(match: re.Match[str]) -> str:
        e = entry(match["key"])
        cell = lambda text: text.replace("|", "\\|")  # noqa: E731
        return f"| `{match['key']}` | `{cell(e['example'])}` | {cell(e['help'])} |"

    def page(match: re.Match[str]) -> str:
        e = entry(match["key"])
        return (f'<tr><th scope="row"><code>{match["key"]}</code></th><td><code>{html.escape(e["example"], False)}'
                f'</code></td><td>{html.escape(e["help"], False)}</td></tr>')

    md_path, html_path = root / GUIDE_MD.relative_to(ROOT), root / GUIDE_HTML.relative_to(ROOT)
    return {md_path: MD_ROW.sub(md, md_path.read_text()), html_path: HTML_ROW.sub(page, html_path.read_text())}


def documented_keys(root: Path = ROOT) -> dict[str, set[str]]:
    """Setting keys (in field-guide form) listed by each guide file."""
    return {
        "md": {guide_key(m["key"]) for m in MD_ROW.finditer((root / GUIDE_MD.relative_to(ROOT)).read_text())},
        "html": {guide_key(m["key"]) for m in HTML_ROW.finditer((root / GUIDE_HTML.relative_to(ROOT)).read_text())},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="report out-of-date files without writing")
    args = parser.parse_args()
    stale = [path for path, text in render().items() if path.read_text() != text]
    if args.check:
        for path in stale:
            print(f"out of date: {path.relative_to(ROOT)} (run python scripts/configuration_guide.py)")
        return 1 if stale else 0
    for path, text in render().items():
        path.write_text(text)
    print(f"updated {len(stale)} file(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
