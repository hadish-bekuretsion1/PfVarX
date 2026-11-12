# Interactive scientific reporting

PfVarX produces a standalone offline HTML report from existing run artifacts. Reporting
reads published scientific results; it does not call variants, refit models, change QC
thresholds, or change input/output paths. No server, network service, CDN or database is
required. Use a current Chrome, Chromium, Edge, Firefox or Safari browser; automated
acceptance checks currently exercise Chromium.

## Open or refresh a report

After a run finishes, open `results/report/pfvarx_dashboard.html` in your browser.
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
the report records the current renderer's SHA-256 separately under
`software.report_renderer`, alongside the original run provenance.

For an in-progress run, missing final artifacts display an unfinished state. These are
artifact-derived states, not a live scheduler subscription: reload a freshly regenerated
report to observe progress. Final interpretation should use a finalized run.

Use **Open report JSON** to import another local `dashboard_data.json` without a server.
The embedded report remains usable without its adjacent JSON file. Import uses the
browser file picker, not filesystem fetch. Invalid required fields produce a clear
error; malformed optional containers become unavailable with provenance warnings.
Dashboard JSON /2 uses the new analysis IDs. Keep the original renderer for /1 reports.
The schema is shipped at `src/pfvarx/reporting/dashboard.schema.json`.
Client validation checks the required identity and structural types; it is not a
scientific validity assessment or a complete JSON Schema engine.

## Sections and scientific views

The sidebar contains eight sections: overview, downstream analyses, variant calling and
annotation, QC, resources, input validation, software/provenance, and troubleshooting.
Every dataset with multiple views has one **Choose visualization** selector. Exactly
one view is visible per dataset. Choices persist during navigation and the browser
session for that run; **Reset dashboard** restores defaults. Different scientifically
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
- **Export PDF** opens the browser print dialog. Select Save as PDF. All eight sections
  print, with the current downstream analysis and selected chart views. Tables print
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

The resulting `benchmark/results/report/dashboard_data.json` can be loaded with
**Open report JSON**. Its **SYNTHETIC DEMO DATA** banner remains visible. The complete
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
