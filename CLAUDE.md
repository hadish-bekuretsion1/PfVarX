# Working on PfVarX

## CLI and GUI changes go together

The CLI (`src/pfvarx/cli.py`, YAML templates in `config/templates/`) and the configuration
GUI (`src/pfvarx/configuration/web.py`, `gui/`) do the same configuration and run the same
stages. Any change to one must be made to the other in the same change, following the table
in `docs/architecture.md` ("CLI and GUI parity"):

- New or changed setting: model, `gui/field-guide.json`, the matching templates, then
  `python scripts/configuration_guide.py` to regenerate the configuration guide.
- Setting help text lives only in `gui/field-guide.json`.
- New or changed GUI tab or stage: matching CLI command and "Terminal equivalent" text,
  and the reverse.
- Shared behaviour belongs in Python used by both, not in `app.js` alone.

`pixi run pytest -q tests/unit/test_cli_gui_parity.py` must pass.

## Documentation goes with every change

Update all documentation in the same change: README.md, docs/*.md and docs/*.html (quick
start, configuration guide, SOP manual), docs/analysis, docs/benchmarking, docs/reporting and
CHANGELOG.md. Grep the md and html files for every name, path, command or setting the change
replaces. Appendices E, F and G of the SOP copy docs/caller.md, docs/configuration.md and
docs/sequencing.md: run `python scripts/sync_sop.py` after editing those guides
(`tests/unit/test_sop_mirrors.py` must pass).

## Development material stays out of releases

Material only developers need (the development guide, notes, experiment scripts, comparison
results) goes in `development/`. `scripts/package.py` never packages it and
`tests/unit/test_distribution.py` checks that, so the release users receive contains the
analysis software only. User documentation must not link to `development/`.

## Before committing

`pixi run lint`, `pixi run typecheck` and `pixi run pytest -q tests` (about 13 minutes).
The pipeline version is set only in `src/pfvarx/__init__.py`.
