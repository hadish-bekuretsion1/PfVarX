# Development material (not released)

This folder holds material for developing and maintaining PfVarX. It is kept in the
repository for the developers and is **left out of the release archive** that users receive
(`scripts/package.py`; checked by `tests/unit/test_distribution.py`).

| File | Contents |
| --- | --- |
| `local_ci.sh` | Local CI: the same checks as `.github/workflows/ci.yaml` (lint, types, all tests, the release archive), run before every merge and release; `--quick` skips the integration tests |
| `PfVarX_Development_Guide_0.0.0.html` | Step-by-step development reference for 0.0.0: branch and release model, repository and module map, workflow rules and their tools, tool versions, project rules, the testing ladder with expected output, how to add settings, parameters, rules and analyses, packaging, the release procedure and bringing fixes to `main` |
| `PfVarX_Development_Guide.html` | How to maintain PfVarX on Ubuntu or WSL2: every tool, package, step, algorithm, parameter and QC gate with its path and purpose; testing parameter values and QC thresholds; measuring, replacing and upgrading tools; the testing ladder; documentation rules; and all Git and GitHub commands for this project |

Open the guide with `xdg-open development/PfVarX_Development_Guide.html` (in WSL2:
`explorer.exe development/PfVarX_Development_Guide.html`).

Its inventory tables are a snapshot of `main` at the date in its header. After tools,
parameters or QC gates change, check `pixi.toml`, `config/parameters/` and
`config/registry/` directly, or regenerate the tables.

Anything else that only developers need (notes, experiment scripts, comparison results)
belongs here too, so that it never reaches a release.
