# Orchestration and command reference

Run commands from the selected workspace root. The Git/source directory is the export
source; the workspace is the installation and execution location. Plain Git cannot
run setup automatically. The explicit `python3 bootstrap.py ...` operation is the
post-clone installation entry point. Review it before running it.

## Select the smallest valid workspace

```bash
# Run these from the source checkout, each into a NEW directory.
python3 bootstrap.py --final-output FO001 --workspace ../alignment-qc  # alignment and QC tools only
python3 bootstrap.py --final-output FO002 --workspace ../sample-wgs
python3 bootstrap.py --final-output FO003 --workspace ../cohort-wgs
python3 bootstrap.py --final-output FO004 --analyses DA005 --workspace ../population
python3 bootstrap.py --final-output FO002 --design targeted_umi --workspace ../targeted-umi
python3 bootstrap.py --profile benchmark --workspace ../benchmark
python3 bootstrap.py --profile benchmark --analyses DA002 DA001 DA007 --rtg --workspace ../benchmark-full
```

DA005 automatically includes DA001 and its shared helpers; DA007 includes DA002 as
its outcome source. Selected algorithms remain scientifically identical. Unselected
method modules/R scripts, cohort/annotation layers where not needed, and map geometry
where unused are omitted from the workspace. Shared policy/parameter metadata is
retained because validation, provenance and the control plane read it. Sample profiles
omit population analysis implementation and its tool dependencies. Locked transitive
dependencies required by selected tools are retained. The manifest/lock names multiple
possible environments but only the selected environment is installed.

`config/installation.json` records the capability boundary/environment, source Git
revision and tracked-dirty flag when exported from Git (null for an archive). Frozen
file hashes remain the exact content identity for the run. Do not delete
or edit it to add capabilities. To move from a sample-only to a cohort/downstream
objective, bootstrap a NEW appropriate workspace and register the same read-only input
paths. Higher output modes, unavailable analyses or changed assay profiles are rejected
with an actionable configuration error. A development/full-source workspace has all
capabilities. `--no-install` exports software and capability metadata only; install the
reported environment then run its `workspace-init` before use.

Bootstrap starts targeted UMI configuration with an **illustrative** `8M+T 8M+T`
structure. Replace it from real wet-lab design before validation. No generated library
identifier or assay setting is automatically established biological truth.

## CLI stages and contracts

Use `./pfvarx COMMAND`; use `./pfvarx COMMAND --help` for options. JSON goes to stdout;
Pixi/tool diagnostics go to stderr. `PFVARX_JSON_OUT` can save command JSON atomically.

| Stage/command | Required state | What is written / success criteria |
| --- | --- | --- |
| `workspace-init` | Installed profile or full source | Three missing editable YAML files and, unless --config-only is set, input folders and an initialization report; never overwrite user config |
| `config-init --purpose analysis|benchmarking` | Full source/compatible selection | Missing config/quality/compute templates only; legacy explicit setup |
| `config-migrate [--apply]` | Workspace with a schema-3 `config/config.yaml` | Reports, then with `--apply` performs, the move to schema 4: `assay_design`, `inputs/metadata/library_catalog.yaml` (targeted panel and UMI protocol) and the manifest columns; keeps `.before-config4` backups |
| `join-metadata --table PATH [--apply]` | Workspace with a sample sheet | Reports, then with `--apply` performs, the join of a tab-separated study table into `samples.tsv` by `sample_id`: only empty cells are filled, file and library columns are never taken, conflicts and unmatched IDs are listed, and coordinates and dates are checked. Keeps a `.before-join-<time>` backup. Exits 1 when conflicts remain |
| `init` | Valid config | Assay/purpose-specific input directories only |
| `discover-inputs --directory DIR` | Actual FASTQ folder | `metadata/discovered_samples.tsv`; REVIEW_REQUIRED, complete unambiguous pairs only |
| `register-inputs --manifest FILE --mode reference|link` | Reviewed TSV and valid config | Registered sample sheet, optional organized symlinks, metadata/registration.json; refuses nonempty sheet overwrite |
| `fetch-input --url URL --sha256 HASH --role ROLE --name NAME` | Explicit HTTPS source/trusted checksum | Only requested inputs/ROLE/NAME plus acquisition metadata; atomic checksum gate and no overwrite |
| `check-layout` | Config and supplied files | Read-only layout report; PASS required |
| `validate` | Registered inputs/resources | Detailed validation JSON, input quality facts TSV and effective configuration/parameter summary; PASS or reviewed PARTIAL. The printed `summary` explains each problem or warning and the steps that fix it |
| `prepare --dry-run` | Config/resources | Read-only list of missing/stale/damaged derived artifacts |
| `prepare` | Accepted input validation | Only necessary reference derivatives, artifact hashes/seal and preparation report; verification PASS |
| `verify-prepared` | Prepared resources and current registry | Input validation plus source/derived checksums, FASTA/FAI/dictionary compatibility; mandatory failures stop execution |
| `ready` | Valid inputs/current preparation | Read-only resource/eligibility pre-plan; ready must be true |
| `dry-run` | Same mandatory preflight | Snakemake DAG check in external temporary directory; no runs/results created. Also prints (on stderr) the expected time, peak memory and storage of each step for these inputs and this final output (see "Time, memory and storage estimates") |
| `pipeline --mode NEW` | Accepted preflight | Frozen run, actual workflow/QC/scientific products, report and policy cleanup. Writes and refreshes the live progress page `runs/<RUN_ID>/progress.html` (opened in the browser when run from a desktop terminal) |
| `resources-learn --run-id ID` | Existing COMPLETED real-data run | Adds the run's measured step times, peak memory and result sizes to this computer's resource history (done automatically after each completed run; use it for runs made before estimates existed). Synthetic runs are skipped |
| `pipeline --mode RESUME --run-id ID` | Unchanged interrupted run | Restart only with compatible source/config; preserved attempt logs |
| `pipeline --mode EXTEND --run-id ID` | Declared compatible extension | New run, conservatively reused artifacts where reuse key permits |
| `orchestrate --through STAGE` | User-reviewed configuration and registered actual inputs | Ordered initialize/validate/prepare/verify/preflight/execute; stop at selected boundary or mandatory failure |
| `report --run-id ID` | Existing stored results | Same HTML/JSON/QC report paths only; no calling/analysis rerun |
| `run-verify --run-id ID` | Existing frozen run | Frozen-file integrity and workspace-drift report |
| `cleanup --run-id ID [--apply]` | COMPLETED, inactive run | Preview or audited deletion of safe regular files under that run's results/work only |
| `benchmark-partition`, `benchmark-evaluate`, `benchmark-resources` | Appropriate benchmark configuration/truth/run | Separate benchmarking procedures in its manual |
| `status` | Any workspace | Effective configuration hash and the latest run's status |
| `version` | Any installation | The pipeline release version |
| `fixture` | Disposable empty workspace | Synthetic software-test reads, reference, GFF, sample sheet and test catalogue (`--design`, `--depth`, `--n-mono`, `--n-poly`, `--variant-scale`); never biological data |

### Commands that match the configuration GUI

Each GUI tab runs the same code as one of these commands, so you can switch between the
GUI and the terminal at any stage. `tests/unit/test_cli_gui_parity.py` keeps them matched.

| GUI tab or button | Command | What it does |
| --- | --- | --- |
| (open the GUI) | `configure --gui` | Local configuration form at 127.0.0.1 with a session token; the YAML files stay usable |
| Input directories | `input-directories` | Create the input folders for the saved configuration (after editing the YAML) |
| Input directories: Install required tools | `install-tools [--check]` | On-demand workspaces: install the locked environment the saved configuration needs (final output, design, input type, analyses and their prerequisites) and record it; nothing when it is installed. `--check` reports whether it is (exit 1 if not) and lists each tool and package with its locked version, purpose and installed state |
| Input files | `input-files`, `add-input --file PATH`, `fetch-input --url URL --name NAME` | List input folders; copy files from disk or download over HTTPS into them, with the same checks as the GUI |
| Input Files Validation | `validate-input-files` | Run `validate` on the saved configuration |
| Inputs Preparation | `inputs-preparation [--metadata-table PATH]` | Run `prepare-samples`, `prepare`, `verify-prepared` and `ready` in order |
| Inputs Preparation (steps 1–2 only) | `prepare-samples [--metadata-table PATH]` | Build the sample metadata `prepared/samples.tsv` from the input files if it is missing; with Study Metadata: Yes, merge the study table from `study_metadata_dir` (none: continue with a warning; several: stop and ask); with No, continue with the sample metadata only |
| (earlier name) | `prepare-reference` | Run `prepare`, `verify-prepared` and `ready` (references only, no sample sheet or merge) |
| Workflow Readiness Checks | `workflow-readiness-checks` | Run `prepare-samples`, `validate`, `prepare`, `verify-prepared` and `ready` in order |
| Workflow Check and Run: Check workflow | `workflow-check` | Run `dry-run`; no outputs are created. Shows the expected time, peak memory and storage per step |
| Workflow Check and Run: Run workflow | `workflow-run` | Check the workflow, then run it as a new analysis. The GUI opens the live progress page in a new tab; it opens the report dashboard when the run ends |

Initialization creates no placeholder runs, fake logs or empty output trees. Run
creation establishes run execution/log/result roots before processing; individual
rules create their required subdirectories. This is intentional demand-driven setup.
Sequencing quality facts before alignment are labelled input-validation observations;
post-alignment QC is computed only when the tools actually run.

## Run progress and time, memory and storage estimates

**Live progress.** While a run works, PfVarX rewrites `runs/<RUN_ID>/progress.html` (and
`progress.json`) every 15 seconds from the Snakemake log. The page shows when the run started,
how long it has run, the steps running now (with sample, start time, time so far and expected
time), a count of finished, running and failed jobs per step, the expected finish time and
the time remaining. The page refreshes itself; when the run ends it shows the final status
and total time and then opens the report dashboard. The progress page opens by itself when
you start a run from a desktop terminal (`./pfvarx workflow-run` or `pipeline`), or in a new
browser tab when you press **Run workflow** in the GUI (the GUI also has an **Open run
progress** link). Set `PFVARX_NO_BROWSER=1` to never open a browser; the page path is
always printed. The page only reads run files and never changes the run.

**Estimates.** After every COMPLETED real-data run, PfVarX records each step's wall time,
CPU time and peak memory with the input size of the sample it processed, and the size of
each results folder, in a history file for this computer
(`~/.local/share/pfvarx/history/<host>.jsonl`, or `$PFVARX_HISTORY_DIR`). A new plan scales
those measurements by its own input sizes (median seconds per input GB, preferring earlier
runs with the same design and input type) and packs samples `max_jobs` at a time. The result
is frozen with the plan (`frozen_plan.json` → `resources.estimate`), printed by
`workflow-check` and the GUI's **Check workflow**, used by the progress page, and compared
with what the run measured in the report's Resources section. Nothing is invented: a step
never run on this computer is UNKNOWN, so the total time is unknown until every step has
history; synthetic test runs are never learned from. Estimates are advice for planning only
and never change results. For runs made before this feature, add them once with
`./pfvarx resources-learn --run-id RUN_ID`.

**Memory requests.** When a job class used more memory in an earlier real-data run on this
computer than its fixed request, the plan raises the request to that peak plus 20 %
(`memory_basis: MEASURED_HISTORY`; at present the variant caller reads it), within
`per_job.max_mem_mb`. Snakemake then starts no more of those jobs at once than `max_mem_mb`
holds. If the peak needs more than `per_job.max_mem_mb`, the check warns
(`RESOURCE_MEMORY_TIGHT`). It also warns (`RESOURCE_MEMORY_IN_USE`) when the memory free at
check time cannot hold `max_jobs` of the largest job, and suggests a `max_jobs` value: other
programs (a browser, an editor) count, and when memory runs out the system can stop the run.

**Stopped runs.** A run records its process; pressing Ctrl+C or a polite kill marks it
INTERRUPTED. A run killed outright (for example by the out-of-memory killer) or ended by a
restart cannot write that itself, so `status`, the GUI and the progress page report it as
INTERRUPTED once its process is gone. Continue it with
`./pfvarx pipeline --mode RESUME --run-id RUN_ID`.

## Organize without duplicating originals

Discovery accepts `SAMPLE_R1.fastq[.gz]`, `SAMPLE_R2.fastq[.gz]`, optional `_L001` lane
and optional `_001` before the extension; `.fq` is also recognized. Duplicate candidates
or orphan mates are rejected. Nonmatching files are reported, not guessed. Review every
sample/lane/library and optional population/location/date field. For other naming
conventions provide an explicit TSV using the ordinary sample-sheet columns.

```bash
./pfvarx discover-inputs --directory /absolute/path/to/fastq
nano metadata/discovered_samples.tsv
./pfvarx register-inputs --manifest metadata/discovered_samples.tsv --mode link
```

Relative input paths in an imported manifest resolve relative to **that manifest's
parent directory**. Registered reference mode uses absolute original paths; link mode
creates workspace input symlinks and links existing BAM/CRAM sibling indexes. Originals
are not copied/rewritten. Missing or invalid fields fail before publishing the sheet.
A populated registered sheet is never overwritten: review/edit it explicitly and
revalidate, or use a new workspace.

References/targets/catalogues are declared in config and may be absolute read-only
paths. Explicit download example (replace all placeholders with verified values):

```bash
./pfvarx fetch-input --url https://provider.example/verified-reference.fasta \
  --sha256 TRUSTED_64_HEX_SHA256 --role reference --name reference.fasta
```

This is an example URL, not an actual supplied reference. HTTPS credentials in the URL
are rejected; HTTPS downgrades, missing checksums and existing destinations are refused.
Query strings are omitted from saved provenance. No guessed reference/truth release is
fetched. Use existing local resources when no suitable public download is provided.

## Preparation verification and recovery

`prepared/manifest.json` now seals required derived artifacts as schema
`pfvarx.prepared/2`. FASTQ entry needs combined-reference BWA-MEM2 index files. Aligned
BAM/CRAM entry needs no BWA index, but does require the relevant FASTA/FAI/dictionary.
Source fingerprint changes, absent/empty outputs and hash changes invalidate prep.
FASTA index and dictionary name/length compatibility is checked separately. Existing
unsealed manifests require `prepare` to seal/rebuild as necessary before execution.

```bash
./pfvarx prepare --dry-run
./pfvarx prepare
./pfvarx verify-prepared
```

Corruption is a failure, not an advisory warning. The repair plan rebuilds the affected
derivative; originals remain read-only. A forged manifest is not an authenticity trust
root; retain independent input checksums and reviewed software. Additional SnpEff
resources are built only when FO004 annotation executes and are reference-checked there.

## Retain what is needed, clean only proven intermediates

```yaml
retention:
  intermediates: successful   # default; or keep for debugging
```

Snakemake always runs with `--scheduler greedy`: it starts whatever jobs fit the limits in
`config/compute.yaml`. PfVarX's jobs are few per sample with fixed limits, so Snakemake's ILP
optimiser adds nothing and could otherwise pause 10 seconds per scheduling decision ("Failed to
solve scheduling problem with ILP solver in time") before falling back to the same greedy
choice. Results do not depend on the scheduler.

`keep` adds Snakemake `--notemp`. `successful` lets Snakemake delete declared temporary
outputs after their dependencies finish, then removes safe remaining results/work
regular files only after a COMPLETED run. A failed run never receives the final bulk
cleanup. Essential BAM/VCF/evidence/callability/QC/frozen/log/benchmark/report outputs
live outside that temporary subtree and are retained. Original configured input paths
and symlinks are excluded. `cleanup.json` records applied removals and bytes. Explicit
cleanup is preview-first and never removes reference preparation or final data.
Workspace locking rejects overlapping runs/cleanup; do not remove a live lock manually.

## Recovery rules

Start from the validation `summary`, which says what failed, why and how to fix it. Then
inspect `metadata/validation_report.json`, preparation reports,
`runs/ID/status.json`, `logs/snakemake.pass*.log`, events and per-tool logs. A locked
workspace means work is active, not permission to start a duplicate. Fix mandatory
input/reference/index/resource failures before execution. Use unchanged RESUME only
for a genuine interruption. Scientific input/mask/reference changes require NEW and
independent checksums; not all external file content is covered by the reuse key.
Never relax scientific gates merely to obtain a preferred output.

Actual reference-tool stdout/stderr is recorded under `metadata/preparation_logs/`
for each executed indexing/dictionary action. Inspect these logs when preparation fails;
failed/empty artifacts are not accepted as a verified preparation.

See the [configuration editor and input procedure](configuration.md) and the
[sequencing and downstream eligibility matrix](sequencing.md).

Inputs Preparation builds the missing sample metadata (`prepared/samples.tsv`) from the input files, or a blank template for BAM/CRAM input, and with Study Metadata: Yes merges the study table into it. Existing files are preserved.
