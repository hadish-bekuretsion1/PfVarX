# Development, verification and packaging

Install the included locked environment before checking source:

```bash
pixi install --locked
pixi run lint
pixi run typecheck
pixi run pytest -q tests/unit tests/scientific tests/integration
```

CI enforces lint, type checking, CLI smoke, unit/scientific tests and integration
preflight. Integration tests copy current source and do not require a Git commit.
Tests use temporary workspaces; no study dataset is bundled. A normal local install
may create caches/environments; those are never exported by the packager.

## What the automated suite covers

`pixi run pytest -q tests` takes about 13 minutes on four cores. Besides unit and
scientific tests, it runs real pipelines on the synthetic fixture:

| Integration test | What it proves |
| --- | --- |
| Dry runs (WGS, targeted non-UMI, targeted UMI) | Each design builds its workflow and creates no run or results |
| Annotated WGS cohort to DA002 | Alignment, calling, cohort merge, BGZF annotation and one analysis |
| Every downstream analysis | An 11-sample WGS cohort runs DA001–DA009 to completion, clonality calls never contradict the simulated truth, DA008 finds no copy-number change where none was simulated, and DA009 calls every judged sample NO_DELETION with the WGS evidence class. `tests/unit/test_cnv.py` checks DA008 on a constructed BAM with a known duplication and deletion; `tests/unit/test_hrp.py` checks DA009's joint states, both evidence models (including a dropout with failed controls), prevalence, co-occurrence and trend |
| Targeted UMI to DA002 (depth 60) | The UMI route reaches cohort, annotation and resistance analysis |
| No sample passing QC (UMI, depth 15) | The run completes with an explicit empty-cohort outcome instead of a crash |
| UMI with short reads | Fixed read structures and adapter dimers do not stop fgbio |

Consistency checks that fail when one part changes without another:

| Test | Fails when |
| --- | --- |
| `tests/unit/test_version.py` | A current-release reference disagrees with `__version__` in `src/pfvarx/__init__.py` |
| `tests/unit/test_cli_gui_parity.py` | A setting lacks GUI help, a template entry or a guide row; guide text differs from the GUI; or a GUI tab or stage has no CLI command |
| `tests/unit/test_validation_guidance.py` | A validation failure stops explaining its cause and fix, or a gene list passes as a core mask |
| `tests/unit/test_no_hardcoded_parameters.py` | Code repeats a maintainer parameter as an inline fallback, a parameter name is used in two layer files, a number in a parameter file is read as text, or `legacy_parameters.yaml` names an unknown parameter |
| `tests/unit/test_resource_plan.py` | A workflow rule writes its threads or memory as a number instead of reading `config/registry/job_resources.yaml`, or a rule with requests is missing from that file |
| `tests/unit/test_microhap.py` (settings) | A bare `off` in a hand-edited configuration stops meaning "off" (YAML reads it as a boolean), a template leaves the value unquoted, or the default changes from `off` |
| `tests/unit/test_microhap.py` (scale) | Locus derivation starts scanning every site per window again: a whole genome tiles into about 160 000 windows, so the test fails long before a real run would hang |
| `tests/unit/test_microhap.py` (derivation) | Loci derived from the run's own panel targets or core-genome windows change: the site minimum, the window length that a fragment must span, the cap on derived loci, and the exclusion of indels and sites no sample passed |
| `tests/unit/test_hrp.py` (calibration) | The self-calibrating `background_percentile` misplaces a deleted, single-copy or duplicated gene among the sample's own windows, or the threshold-sensitivity grid stops marking the calls that depend on a provisional cut-off |
| `tests/unit/test_cohort_matrix.py` | The cohort evidence store returns values different from the per-cell mapping it replaced, rounds `wsaf` (which is compared against reporting floors), treats a missing record as present, or regresses towards per-cell dicts in memory |
| `tests/unit/test_matrix_declarations.py` | An analysis's `USES_COHORT_MATRIX` disagrees with its code, following imports inside the package: a wrong declaration would hand the analysis an empty matrix instead of failing |
| `tests/unit/test_microhap.py` | Fragment-level haplotype counting (including overlapping and disagreeing mates, and a low-quality base), the haplotype floors, locus resolution and territory rules, the whole `build` path, the locus-basis choice between microhaplotypes and biallelic SNVs, and that the shipped panel holds no invented loci |
| `tests/unit/test_dashboard_views.py` | A dashboard analysis tab builds its charts in `daKpis` (which returns early, so they never render), an analysis lacks a tab, headline numbers or a view block, or a visual group is built outside `daViews`. No Python test executes the dashboard JavaScript, so this structural check stands in for it |

`pixi run python scripts/check_installs.py` compares each kind of installation with the tools
its final output runs, from `pixi.lock` (see [deployment](deployment.md#check-what-an-installation-contains)).

After editing `docs/caller.md`, `docs/configuration.md` or `docs/sequencing.md`, run
`python scripts/sync_sop.py`: SOP appendices E, F and G copy them, and
`tests/unit/test_sop_mirrors.py` fails while they differ.

After changing setting help in `src/pfvarx/configuration/gui/field-guide.json`, run
`python scripts/configuration_guide.py` to regenerate the configuration guide
(`--check` only reports whether it is out of date).

The fixture command takes `--variant-scale N` to multiply the simulated variant sites
(default 1). The all-analyses test uses 4, so relatedness, structure and phylogeny have
enough informative sites on the 48 kb synthetic genome.

## End-to-end software verification

```bash
pixi run verify --output /absolute/path/to/new-verification
```

This runs WGS FO004, targeted/UMI FO002, benchmarking FO004, WGS FO001, WGS FO003 and targeted automatic FO004 in
separate newly created workspaces. Use scenario names to restrict the run:

```bash
pixi run verify fo001 fo003 --output /absolute/path/to/new-modes
pixi run verify wgs benchmark --output /absolute/path/to/new-science-workflows
```

Run `pixi run verify auto --output /absolute/path/to/new-auto-check` to check that
missing annotation and optional metadata block only the affected targeted analyses.

The destination/scenario directory must not exist. The runner checks input layout,
validation, preparation idempotency, readiness, dry-run isolation, expected alignment
and call files, cohort contracts, requested terminal manifests, frozen integrity and
benchmark development/held-out evaluation/refusal of repetition. Scientific eligibility
skips are expected outcomes and do not prove that the skipped estimator executed.
Fixtures and outputs are generated outside the pristine source when an external output
path is used. They are software tests, not biological validation.

## Dashboard browser checks

Browser dependencies are separate from the analysis lock. Use an external virtual
environment if the source tree must remain pristine:

```bash
pixi run python -m venv /absolute/path/to/browser-env
/absolute/path/to/browser-env/bin/pip install -r tests/dashboard/requirements.txt
/absolute/path/to/browser-env/bin/python -m playwright install chromium
/absolute/path/to/browser-env/bin/python tests/dashboard/check_dashboard.py /path/to/run/results/report /path/to/browser-checks
```

Alternatively, prefix the check with `PFVARX_BROWSER_EXECUTABLE=/path/to/chrome` to use
an existing Chrome/Chromium executable. Restricted systems need permission to launch
headless browser processes. The suite checks navigation, hover/drill-down, table search
and sort, CSV/XLSX/PNG exports, print-to-PDF, keyboard tooltips, downstream or
empty state, offline operation and console errors. Where present it also checks
every available PCA/FST, relatedness, map and tree visualization selector. Use WGS and sufficiently
informative benchmarking reports to exercise skipped and executed science paths.

Each check writes its JSON, screenshots and export artifacts to the specified output
folder. Network requests are blocked and fail the check. `test_dashboard.py` is a pytest
adapter requiring Playwright on its Python path and `PFVARX_TEST_REPORT` pointing at a
completed report. The table search derives the sample prefix from the report, so real
accession IDs are supported.

## Fault injection

For tests only, `PFVARX_FAULT_INJECT=pfvarx_call:SAMPLE_ID` deliberately fails the named
call job with exit 70. This tests independent-sample continuation and exclusion logs.
Never leave it set for a scientific run. Inspect `runtime_exclusions.json`, workflow
passes, terminal manifests and denominators; a partial workflow is not a clean run.

## Build a software-only archive

```bash
pixi run package --output /absolute/path/to/distribution
```

This writes `pfvarx-<version>-source.zip` (the version printed by `./pfvarx version`) and its `.sha256` sidecar. Run `sha256sum -c` from
the output directory to verify the archive. Export uses explicit software roots and
file types, never raw workspace copying. Active config YAML, input/reference/truth
data, profiles measured on a study, runs, results, environments, caches and symlinks
are excluded or rejected. Development-only material in `development/` (the development
guide, notes, experiment scripts) is never packaged. Fixed ZIP timestamps, modes and ordering make repeated
packages from identical source byte-identical on the same packaging runtime.

To assert that an unpacked release contains only software and no runtime directories:

```bash
python scripts/package.py --check
```

This requires only the Python standard library; no environment installation or project
initialization is needed. Installing Pixi or running tests in that same directory will
make it a working checkout rather than pristine source. Package it using the allowlist,
or unpack the source ZIP again, to obtain a fresh pristine copy.

The distribution tests check runtime/private-data exclusion, reproducibility, unpacked
roundtrip integrity and source-symlink rejection. Keep validation outputs outside the
release, then verify the archive and its extracted tree before distributing it.

The dashboard suite also checks local JSON import, schema errors, missing optional data,
missing-versus-zero values, view persistence, pagination and export escaping. Pillow and
openpyxl inspect actual PNG and XLSX downloads. See [reporting](reporting/README.md).

## Orchestration regression checks

From full source/development installation:

```bash
pixi run pytest -q tests/unit/test_orchestration.py
python3 bootstrap.py --final-output FO002 --workspace /absolute/path/to/new-sample
cd /absolute/path/to/new-sample
./pfvarx fixture --n-mono 1 --n-poly 0 --depth 20
# Set an explicitly synthetic reference release label before this smoke run.
./pfvarx orchestrate
./pfvarx verify-prepared
```

Use only a disposable empty workspace for fixture generation. The unit suite tests
selective closure/export, reference/link registration, discovery ambiguities, failed
HTTPS checksums/no overwrite, actual FAI/dictionary/hash corruption, capability errors,
locks and conservative retention. The full six-case runner uses the current source
and tests both omitted and executed scientific paths. Container context generation
can be checked without Apptainer, but is not image/E2E qualification.

## Configuration editor browser test

Use a browser environment with Playwright and the pipeline's installed Python.

```bash
PFVARX_GUI_PYTHON=/absolute/path/to/pixi/environment/bin/python \
PFVARX_BROWSER_EXECUTABLE=/usr/bin/google-chrome \
  /absolute/path/to/browser-env/bin/python tests/dashboard/check_configure.py
```

The test creates a disposable workspace, changes assay and analysis choices,
checks field visibility and exclusivity, verifies the input browser and saves
the configuration twice. It checks that no sample sheet or preparation was created.

## Dry-run-only verification

From the release source root with the locked development environment available:

```bash
pixi run python scripts/verify.py --dry-run-only --output "/absolute/path/PfVarX preflight"
```

The verifier prepares synthetic input resources in disposable workspaces. It
requires exit zero from Snakemake for each selected scenario and checks that no
run or result directory was created. The per-scenario JSON links to its stage
logs. Jobs behind runtime checkpoints are checked during execution.
