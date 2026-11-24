# Interactive scientific reporting

PfVarX produces a standalone offline HTML report from existing run artifacts. Reporting
reads published scientific results; it does not call variants, refit models, change QC
thresholds, or change input/output paths. No server, network service, CDN or database is
required. Use a current Chrome, Chromium, Edge, Firefox or Safari browser; automated
acceptance checks currently exercise Chromium.

## Open or refresh a report

After a run finishes, open `results/report/pfvarx_dashboard.html` in your browser.

**Where the results are saved** (top of Overview) gives the run folder's full path with
**Copy path** and **Open results folder**, a button for each output folder (Alignments,
Quality control, Variant calls, Cohort, Annotation, Eligibility, Downstream analyses, Report,
logs, benchmarks, frozen settings, software versions), and every final output file grouped by
folder, with its size and what it is. Nothing is downloaded: **Open folder** shows the file's
folder in the browser, and **Copy path** copies its full path. File names are links only for
files the browser displays (HTML, JSON, TXT, logs, images, PDF); BAM, VCF, TSV, BED and other
files open with their own program (IGV, a spreadsheet, a text editor), so paste the copied
path there or into the file manager (Ctrl+L). Search finds a file by name, sample or type.
Links are relative to `results/report/`, so they work when the dashboard is opened from its
run folder on the computer that ran PfVarX; a copy sent elsewhere keeps the list but not the
files. Temporary work files are not listed, and tool byproducts that no step reads are not kept.
The same list is in `dashboard_data.json` under `outputs`.
The existing `results` link selects the current run. Each run retains the same paths
beneath `runs/<run-id>/results/`:

| Artifact | Contract |
| --- | --- |
| `report/pfvarx_dashboard.html` | Single HTML file containing styles, scripts, data and bundled map geometry |
| `report/dashboard_data.json` | Local data backing the report; schema `pfvarx.dashboard/2` |
| `qc/centralized_quality_matrix.json` | Centralized QC gate evidence and sample states |

Refresh reporting only from the workspace root with the existing command:

```bash
./pfvarx report --run-id YOUR_RECORDED_RUN_ID
```

This overwrites those three report artifacts, leaving scientific files and frozen run
configuration untouched. Archive a previous report first if you need its exact bytes.
New reporting assets must be installed in that workspace's software package before
refreshing. Selected workspaces use `./pfvarx`; full-source operation remains supported.
Refreshing a report does not silently upgrade a frozen scientific snapshot;
the report records the current renderer's version and SHA-256 separately under
`software.report_renderer`, alongside the original run provenance. The pipeline version
shown is the one recorded in the run's `snapshot.json` when the run was made, not the
version installed when the report is rebuilt. The software section lists the probed
command-line tools and the exact version and build of every scientific package in the
run's environment (fgbio, Picard, SnpEff, IQ-TREE, hmmIBD, MOIRE, Dcifer, mgcv and
others), marked "(package)".

For an in-progress run, missing final artifacts display an unfinished state. These are
artifact-derived states, not a live scheduler subscription: reload a freshly regenerated
report to observe progress. Final interpretation should use a finalized run.

Each run's dashboard carries its own data, so open that run's `pfvarx_dashboard.html` to
see another run; the embedded report remains usable without its adjacent JSON file. Invalid required fields produce a clear
error; malformed optional containers become unavailable with provenance warnings.
Dashboard JSON /2 uses the new analysis IDs. Keep the original renderer for /1 reports.
The schema is shipped at `src/pfvarx/reporting/dashboard.schema.json`.
Client validation checks the required identity and structural types; it is not a
scientific validity assessment or a complete JSON Schema engine.

## Clonality report

Every run with variant calls writes `results/report/clonality_report.txt`, one readable
line per sample, and the same rows in `clonality_report.tsv`; the dashboard overview
shows them as a table:

```text
Sample X is monoclonal
Sample Y is polyclonal with 3 clones: Clone 1: 50.0%, Clone 2: 30.0%, Clone 3: 20.0%
```

The number of clones and the percentages come from the PfVarX-Call mixture model: the
number of DNA components it selects and their weights, largest first. When DA001 ran,
its MOIRE monoclonal/polyclonal call is combined with it; MOIRE's own clone count is not
used, because one sample or a small cohort cannot support it. The report says so when the
two disagree (INDETERMINATE), when DA001 sees more clones than the mixture model resolves
("not resolved"), and when the largest number of clones tested was selected ("at least N").
Clones below 2% of the DNA cannot be resolved. Percentages are shares of parasite DNA,
which approximate but are not exactly the share of parasites, and stay uncalibrated until
the caller is benchmarked. `./pfvarx report --run-id ID` rebuilds the report for an
earlier run without rerunning any analysis.

## Sections and scientific views

The dashboard is titled **PfVarX Genomic Analysis Dashboard** (navy header with the teal
PfVarX mark, centred; on phones and tablets the **Menu** button sits on the same row and
the section list slides in from the left). The header has no buttons. Nothing else sits above the sections. The project, run, final output, design and purpose,
the data kind (SYNTHETIC DEMO DATA for fixtures) and the calibration status appear once,
at the top of **Overview**. The sidebar
contains eight sections, in the order a run is reviewed: **Overview** (run
summary), **Inputs** (validation and libraries), **Quality Control** (gates and metrics; for FO003 and FO004 also the cohort QC gate G5 and each sample's cohort QC verdict),
**Variants** (records and annotations), **Downstream Analyses** (DA001–DA007),
**Resources** (execution and utilization), **Software** (versions and provenance) and
**Logs** (diagnostics and exclusions). Each section opens with one heading, such as **PfVarX · Inputs**; Downstream Analyses
shows the name of the open analysis instead. Each is a card with an icon; the open section is
dark with a teal bar. On a phone, **Menu** opens the same cards.
**Downstream Analyses** has one tab per analysis across the top: Infection Complexity
(DA001), Resistance Markers (DA002), Genetic Diversity (DA003), Genetic Relatedness (DA004),
Population Structure (DA005), Phylogenetic Analysis (DA006) and Spatiotemporal Analysis
(DA007). The open tab is navy with a teal bar. Each tab
opens the analysis's own page: its full name as the title (for example *Antimalarial
Resistance Marker Analysis*, *DA002 · Marker calls, callability and prevalence*), summary
tiles (for DA002: samples, markers, callable denominator, UNKNOWN), and its charts and
evidence tables side by side under one set of filters. DA002 adds callability by marker;
DA001 adds P(COI > 1) by sample and the classification share; DA003 adds nucleotide
diversity; DA004 adds estimates with 95% confidence intervals; DA005 adds explained
variance. The run status and denominators of the analysis follow at the bottom.

In the other sections, every dataset with multiple views has one **Choose visualization**
selector. Exactly one view is visible per dataset there. Choices persist during navigation and the browser
session for that run; closing the tab restores the defaults. Different scientifically
distinct datasets have separate groups and explicitly identified denominators.

| Analysis | Default and available evidence |
| --- | --- |
| DA002 | Prevalence forest plot with supplied or explicitly calculated Wilson intervals; bar, state heatmap, table; map only with actual coordinates |
| DA003 | Supplied population diversity comparisons, supplied point-estimate distributions and tables; FWS is a separate dataset |
| DA001 | COI posterior medians with supplied bounds; distributions, class counts and table; technical mixture K is not COI |
| DA004 | Network for small cohorts, heatmap for larger cohorts; actual pairs, method/status/sample filters and optional display threshold |
| DA005 | PCA scatter with actual components and supplied explained variance; coordinates table; FST has its own matrix/table |
| DA006 | Supplied Newick tree, collapse/expand and support table; distance matrix only if IQ-TREE published `tree.mldist` |
| DA007 | Actual-coordinate grid map/table; separate country-month heatmap, time series and table for the corresponding source summary |

**Overview** and **Resources** show when the run started and finished, its total time and
the time expected before it ran. **Resources** then compares, step by step, the estimate
frozen with the plan (learned from earlier real-data runs on this computer, see
`docs/orchestration.md`) with what the run measured: time, peak RAM and the size of each
results folder. A run without an estimate says so.
While a run works, its live page `runs/<RUN_ID>/progress.html` shows progress; the report
dashboard is written when the run ends.

**Resources** continues with **Resources by step: tool, time and RAM**, one row per workflow
step (fastp read QC, bwa-mem2 alignment, Picard duplicate marking, the PfVarX caller,
bcftools, SnpEff, each analysis DA001–DA007, and the others): its layer, the tools and
versions that ran it, the number of jobs (samples, lanes or analyses), total and longest
job time, CPU time, peak and mean RAM, average CPU cores and its share of the summed time.
Charts show time, peak RAM or CPU time per step. **Resources by job** below gives the same
measurements for every sample, lane or analysis. Every step records its measurements in
`runs/<RUN_ID>/benchmarks/<step>/<job>.tsv`; time and RAM are Snakemake's measurements
(wall seconds and maximum resident memory). Steps run in parallel can overlap, so summed
time can exceed the run's elapsed time.

Filters select existing records and never refit an estimator or alter eligibility.
DA002 defaults to the ALL population group where present; DA004 starts with one recorded
method and requires a single method for a matrix/network; DA007 starts with one actual
marker. Reset filters to inspect all source rows. Spatial grid records lack country and
month metadata, so those controls apply only to the country-month summary. Heatmaps
require the relevant marker selection when multiple markers would give ambiguous cells.

Hover or keyboard-focus a chart mark for its actual identity, value, uncertainty and
available denominator. Axis titles identify units. Maps use the bundled Natural Earth
geometry with attribution. Missing coordinates never receive invented locations.
Missing branch lengths trigger a clearly labelled topological tree layout. Missing
measurements, skipped methods and unavailable uncertainty remain distinct from zero.
No violin density, posterior samples, historical control chart or execution timeline
is fabricated from summary-only inputs. Distribution quartiles describe the supplied
point estimates, not confidence intervals.

## Colours and what they mean

The whole dashboard uses one set of colours, each with an icon so colour is never the only cue:

| Look | Meaning | Examples |
| --- | --- | --- |
| green ✔ | accepted / done | PASS, COMPLETED, ELIGIBLE |
| amber ▲ | accepted with a warning | WARNING, LOW_CONFIDENCE, CONDITIONALLY_ELIGIBLE, COMPLETED_WITH_EXCLUSIONS |
| red ✖ | rejected / failed / blocked | FAIL, FAILED, INELIGIBLE, BLOCKED_…, FILTERED |
| orange ◆ | a finding (biology, not quality) | MUTANT, ALT_DOMINANT |
| blue ○ | a category (biology, not quality) | MONOCLONAL, MIXED, WILD_TYPE, REF_DOMINANT |
| gray – | not applicable, not evaluated or not callable | NOT_APPLICABLE, NOT_EVALUATED, INSUFFICIENT_DEPTH, UNCALLABLE |

In tables, a QC metric (for example `median_depth_core`, `duplicate_fraction`,
`cohort_missing_fraction`) is coloured against the run's own QC policy limits: green within,
amber ▲ past the warning limit, red ✖ past the failure limit. Other fractions between 0 and 1
(WSAF and its interval, prevalence, P(COI > 1)) show a small blue bar behind the number for
magnitude, with no judgement. Missing values stay a plain dash. Overview tiles take the same
colours: failed samples red when any, QC warnings amber when any (read/alignment or cohort
QC). Charts use the same status colours; clonality uses blue (monoclonal) and orange
(polyclonal), a pair checked for colour-vision deficiency.

## Filters, tables and exports

Tables provide text search, available categorical filters, keyboard-accessible sorting,
50-row pagination, reset, and explicit result counts. Charts have bounded display sizes
with visible truncation notes; table/CSV/XLSX exports retain all filtered embedded rows.
The backend's configured variant embedding cap is unchanged. Full variant totals and
category counts are harvested before that cap; detailed variant charts use the labelled
embedded subset. Consult the evidence TSV paths for omitted records.

- **CSV** exports all filtered rows, quotes commas/newlines correctly, and preserves
  missing cells. Provenance lines identify project, run, dataset and filter settings.
- **XLSX** writes a genuine Excel workbook, with data and provenance sheets, numeric
  cells where numeric values are supplied, and empty cells for unavailable data.
- **PNG** exports the current SVG view on white, at least 1200 pixels wide, with title
  and project/run caption. Zoom affects the exported visible view; reset zoom for full
  extents. Export feedback states whether a download was prepared or failed.
- **PDF**: print from the browser (Ctrl+P, or Share → Print on a phone) and choose Save
  as PDF. All eight sections
  print, with the open downstream analysis (all its charts and tables) and selected chart views. Tables print
  their current page. Use CSV/XLSX for complete tables and print each analysis separately
  when a DA001–DA007 supplement is required. Browser headers/footers are user settings.

**Pseudonymize sample/location fields in table exports** affects matching columns only.
It does not anonymize the HTML, plots, free text, input paths or every identifier.
Logs contain a bounded tail (up to 60 files, 64 KiB/file and 512 KiB total), with displayed
truncation indicators. Secret-named values, common token assignments, Bearer tokens,
URL credentials and private-key blocks are redacted before embedding. External log
symlinks are not followed. Review reports before sharing; paths and scientific
identifiers remain present. No process environment is harvested by the reporter.

## Synthetic example and validation

The source contains no generated demonstration JSON or historical reports. To create
an explicitly synthetic example outside the source tree, use a full/development
installation containing the verification script and all selected-method dependencies:

```bash
pixi run verify benchmark --output /absolute/path/to/new-verification
```

The resulting `benchmark/results/report/pfvarx_dashboard.html` shows **SYNTHETIC DEMO
DATA** at the top of Overview. The complete
synthetic generator and verification code are retained; generated data is not shipped.
This is a software fixture, not an experimental finding or accuracy claim.

Run core tests from a full/development installation with `pixi run test`. Browser dependencies are separate from the analysis
lock; see [testing](../testing.md). The acceptance runner exercises all available views,
navigation and persistence, true-data tooltips, filters, pagination, schema errors,
missing values, local JSON import, CSV/XLSX/PNG/PDF exports, and offline operation.
Use both an informative synthetic FO004 benchmark and FO001/FO002/skipped-analysis
reports. Scientific methods and all pre-existing input/output contracts are outside
this reporting update. Analysis and benchmarking instructions remain in their separate
[analysis](../analysis/README.md) and [benchmarking](../benchmarking/README.md) manuals.

## Implementation map

`build.py` harvests existing artifacts and inlines the frontend at build time.
`dashboard_template.html` supplies the shell; `dashboard.css` styles screen and print;
`core.js` validates input and manages state/filter/view groups; `tables.js` handles
tables; `charts.js` renders SVG; `views.js` connects datasets to the eight sections;
`exports.js` implements downloads and print. These assets ship with the source and
Python package; generated HTML has no external runtime assets.
