# Orchestration and command reference — 0.0.0

Run commands from the selected workspace root. The Git/source directory is the export
source; the workspace is the installation and execution location. Plain Git cannot
run setup automatically. The explicit `python3 bootstrap.py ...` operation is the
post-clone installation entry point. Review it before running it.

## Select the smallest valid workspace

```bash
# Run these from the source checkout, each into a NEW directory.
python3 bootstrap.py --profile sample --workspace ../sample-wgs
python3 bootstrap.py --profile cohort --workspace ../cohort-wgs
python3 bootstrap.py --profile downstream --analyses DA005 --workspace ../population
python3 bootstrap.py --profile sample --design TARGETED_UMI --workspace ../targeted-umi
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
| `workspace-init` | Installed profile or full source | Three missing editable YAML files, config-required input folders and metadata/init report; never overwrite user config |
| `config-init --purpose analysis|benchmarking` | Full source/compatible selection | Missing config/quality/compute templates only; legacy explicit setup |
| `init` | Valid config | Assay/purpose-specific input directories and blank sample-sheet header only |
| `discover-inputs --directory DIR` | Actual FASTQ folder | `metadata/discovered_samples.tsv`; REVIEW_REQUIRED, complete unambiguous pairs only |
| `register-inputs --manifest FILE --mode reference|link` | Reviewed TSV and valid config | Registered sample sheet, optional organized symlinks, metadata/registration.json; refuses nonempty sheet overwrite |
| `fetch-input --url URL --sha256 HASH --role ROLE --name NAME` | Explicit HTTPS source/trusted checksum | Only requested inputs/ROLE/NAME plus acquisition metadata; atomic checksum gate and no overwrite |
| `check-layout` | Config and supplied files | Read-only layout report; PASS required |
| `validate` | Registered inputs/resources | Detailed validation JSON, input quality facts TSV and effective configuration/parameter summary; PASS or reviewed PARTIAL |
| `prepare --dry-run` | Config/resources | Read-only list of missing/stale/damaged derived artifacts |
| `prepare` | Accepted input validation | Only necessary reference derivatives, artifact hashes/seal and preparation report; verification PASS |
| `verify-prepared` | Prepared resources and current registry | Input validation plus source/derived checksums, FASTA/FAI/dictionary compatibility; mandatory failures stop execution |
| `ready` | Valid inputs/current preparation | Read-only resource/eligibility pre-plan; ready must be true |
| `dry-run` | Same mandatory preflight | Snakemake DAG check in external temporary directory; no runs/results created |
| `pipeline --mode NEW` | Accepted preflight | Frozen run, actual workflow/QC/scientific products, report and policy cleanup |
| `pipeline --mode RESUME --run-id ID` | Unchanged interrupted run | Restart only with compatible source/config; preserved attempt logs |
| `pipeline --mode EXTEND --run-id ID` | Declared compatible extension | New run, conservatively reused artifacts where reuse key permits |
| `orchestrate --through STAGE` | User-reviewed configuration and registered actual inputs | Ordered initialize/validate/prepare/verify/preflight/execute; stop at selected boundary or mandatory failure |
| `report --run-id ID` | Existing stored results | Same HTML/JSON/QC report paths only; no calling/analysis rerun |
| `run-verify --run-id ID` | Existing frozen run | Frozen-file integrity and workspace-drift report |
| `cleanup --run-id ID [--apply]` | COMPLETED, inactive run | Preview or audited deletion of safe regular files under that run's results/work only |
| `benchmark-partition`, `benchmark-evaluate`, `benchmark-resources` | Appropriate benchmark configuration/truth/run | Separate benchmarking procedures in its manual |

Initialization creates no placeholder runs, fake logs or empty output trees. Run
creation establishes run execution/log/result roots before processing; individual
rules create their required subdirectories. This is intentional demand-driven setup.
Sequencing quality facts before alignment are labelled input-validation observations;
post-alignment QC is computed only when the tools actually run.

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

`keep` adds Snakemake `--notemp`. `successful` lets Snakemake delete declared temporary
outputs after their dependencies finish, then removes safe remaining results/work
regular files only after a COMPLETED run. A failed run never receives the final bulk
cleanup. Essential BAM/VCF/evidence/callability/QC/frozen/log/benchmark/report outputs
live outside that temporary subtree and are retained. Original configured input paths
and symlinks are excluded. `cleanup.json` records applied removals and bytes. Explicit
cleanup is preview-first and never removes reference preparation or final data.
Workspace locking rejects overlapping runs/cleanup; do not remove a live lock manually.

## Recovery rules

Inspect printed JSON problems, `metadata/validation_report.json`, preparation reports,
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
