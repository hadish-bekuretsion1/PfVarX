# Configure inputs, validate and prepare

Run commands from the installed workspace root. The GUI is optional and runs on
loopback using the same Pydantic models and installation capability checks as YAML.
It saves configuration and creates input folders. Use CLI commands for all data
inspection and preparation.

## Open the configuration editor

```bash
./pfvarx configure --gui
```

For a GUI-first workflow, install the locked environment before initializing the
workspace. In a full source copy run `pixi install --locked`. With a selective export,
use bootstrap's `--no-install` option, then run its printed Pixi install command
and launch the GUI in place of `workspace-init`.

Choose the analysis or benchmarking template in a fresh workspace. If YAML exists,
the form loads it. Complete project ID, input type, assay, reference release and
paths. The form shows targeted BED/primer fields and the UMI read structure only
where applicable. Choose an explicit analysis set or automatic selection for FO004.
The form validates the compute envelope, caller options and dependent fields using
the Python backend. The summary shows the YAML-equivalent settings.

Enter a workspace-relative or absolute file path. The form checks path configuration syntax.
The `validate` command checks the actual files. Select **Create Configuration and
Input Directories** when ready. Existing YAML with changes requires **Save with
backup**. Backups are stored under `config/backups/REVISION/`. Refresh if someone
edits the YAML while the form is open. Repeated identical saves retain the files.

On a remote Linux machine, start the editor with a fixed loopback port:

```bash
./pfvarx configure --gui --no-browser --port 8765
```

In a local terminal, forward that port and open the exact URL printed by the
server. Treat the URL as a private session link.

```bash
ssh -L 8765:127.0.0.1:8765 USER@HOST
```

For headless configuration, use `./pfvarx config-init --purpose analysis`, edit
`config/config.yaml`, `config/quality.yaml` and `config/compute.yaml`, then use
`./pfvarx init` to create input folders. GUI saves also create only folders.

Both GUI save and CLI initialization use the same workspace-relative input
folder layout. The GUI preview and the CLI `init` report show directory paths
with a trailing slash. File fields still require filenames.

## Place and register input data

| Folder | Place or reference these files |
| --- | --- |
| `inputs/reference/` | Parasite FASTA, and core BED. Add a matching GFF for explicit FO004 or to enable markers in automatic FO004. |
| `inputs/fastq/` | Paired FASTQ for FASTQ entry. Register read pairs with their verified sample and lane identities. |
| `inputs/alignments/` | BAM or CRAM with the appropriate index for aligned entry. |
| `inputs/metadata/` | Reviewed sample TSV and optional population and geographic fields. |
| `inputs/targets/` | Reference-matched target BED and primer TSV for targeted assays. |
| `inputs/truth/`, `inputs/partitions/` | Benchmark truth, confidence regions and partition resources for benchmark workspaces. |

For conventionally named FASTQs, use the following and review the proposed sheet
before registration. For other inputs, supply a reviewed manifest.

```bash
./pfvarx discover-inputs --directory inputs/fastq
nano metadata/discovered_samples.tsv
./pfvarx register-inputs --manifest metadata/discovered_samples.tsv --mode reference
```

## Transfer raw data and run preflight

After saving, open **Input files** to upload reads or alignments, or download
them using an HTTPS link and a destination filename. A trusted SHA-256 is
optional. When provided, PfVarX checks it before saving. Otherwise, PfVarX
records the computed checksum and marks the download unverified. Put reference
files in `inputs/reference/`. Register your reads or alignments before running
**Workflow Readiness Checks → Run readiness checks**. The tab runs `validate`, `prepare`,
`verify-prepared`, and `ready` in that order and stops on failure.

## Validate, prepare and execute

```bash
./pfvarx validate
./pfvarx prepare --dry-run
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline --mode NEW
```

`./pfvarx validate` (and the GUI Input Files Validation tab) prints a `summary`: one
entry per problem or warning, naming the input, why it failed, why the input matters and
numbered steps that fix it. For example, a missing core-genome mask reports the expected
path, notes when `config/config.yaml` gives only a folder, lists files of the same type
already in that folder, and explains how to build the mask from the standard core-genome
regions. A core mask that looks like a gene list (gene names, thousands of or overlapping
intervals, or coverage far from the usual 89% of 3D7) passes with a `CORE_MASK_SUSPECT`
warning, because the format is valid but the content is probably wrong.
The same fields (`file`, `used_for`, `why`, `actions`) are stored for every finding in
`metadata/validation_report.json`, with the sample-level statuses. Repair
mandatory failures before proceeding. Under the existing partial-admission policy,
only validated samples enter execution. Prepare rechecks validation, builds the
necessary reference derivatives under `prepared/reference/`, and records tool logs
in `metadata/preparation_logs/`. The `prepared/` tree is a sibling of `inputs/`.
Current sealed derivatives are reused. A changed or damaged derivative must pass
preparation and verification again. Originals stay at their registered paths.
Execution requires verified preparation and records a frozen run, QC decisions,
selected outputs and an offline dashboard.

On a SLURM cluster, the supported approach is to allocate a single node and run
the local Snakemake executor inside that allocation. Set `max_cores`, `max_mem_mb`
and `max_jobs` within the allocation, and use a shared filesystem for the workspace
and inputs. Example from the workspace root:

```bash
sbatch --cpus-per-task=8 --mem=24G --time=08:00:00 \
  --wrap='./pfvarx pipeline --mode NEW'
```

Adjust resources and partition to your site. Per-rule SLURM scheduling has no
qualified execution profile in this release. The single-node example requires a
site-specific smoke test.

Select the input type before the final output. FASTQ and FASTQ.GZ enter at
Layer 1 and can stop at FO001. BAM and CRAM enter at Layer 2 and start at
FO002. Bootstrap a matching input profile when changing entry types in a
selective installation.

The executor menu records Local, HPC, Server or Cloud. These use the same
single-host workflow engine. Submit HPC runs inside one allocated node using
the deployment guide. Server and Cloud refer to Linux hosts and cloud VMs.

For older workspaces, remove the retired host reference setting from the reference block. Run preparation again before starting a new run.

## Add study metadata

Keep study metadata under `inputs/metadata/`. Join a separate metadata TSV
to the configured sample registry by `sample_id` before validation. Preserve
the read or alignment paths, lane and library IDs. Repeat consistent
sample-level values across lanes.

| Column | Example | Use |
| --- | --- | --- |
| `country` | Ethiopia | DA007 |
| `latitude` | 13.50 | DA007, decimal degrees |
| `longitude` | 39.47 | DA007, decimal degrees |
| `collection_date` | 2025-07-18 | DA007, YYYY-MM-DD |
| `population` | Site_A | Grouped analyses, or use the configured grouping column |
| `age_group`, `study_arm` | Adult, Arm_A | Optional study records |

Review the coordinates included in the dashboard before sharing it.

The standard folder prefixes shown for the sample registry, reference FASTA,
core BED, GFF, target BED and primer TSV are accepted in GUI and YAML. They
resolve to the corresponding template filenames. Enter a complete path when
your filename differs. Validation still checks that each required file exists.

The **Input directories** tab saves configuration and creates folders. In
**Input files**, select a file category to upload or download into its folder.
Use the filenames in your saved configuration. Register sample files, then
use **Input Files Validation** followed by **Reference Preparation**. Preparation revalidates inputs,
generates FASTA indexes and the sequence dictionary under `prepared/reference/`,
and builds BWA-MEM2 indexes for FASTQ entry. Assay-specific read processing and
SnpEff annotation remain workflow stages. **Workflow Readiness Checks** runs the checks together.

## Terminal workflow

From your study workspace, edit the three configuration YAML files and run
the corresponding commands for the GUI stages.

| GUI stage | CLI equivalent |
| --- | --- |
| Input directories | `./pfvarx input-directories` |
| Input files: upload from your computer | `./pfvarx add-input --file /path/to/file` (repeat `--file` for several) |
| Input files: download from an HTTPS link | `./pfvarx fetch-input --url URL --name FILENAME` |
| Input Files Validation | `./pfvarx validate-input-files` |
| Reference Preparation | `./pfvarx prepare-reference` |
| Workflow Readiness Checks | `./pfvarx workflow-readiness-checks` |
| Workflow Check and Run: Check workflow | `./pfvarx workflow-check` |
| Workflow Check and Run: Run workflow | `./pfvarx workflow-run` |

Each GUI stage and its command run the same code and apply the same rules, so
you can switch between them at any stage. Both refuse to create input folders
while `reference.release_label` is still the template placeholder.

`add-input` and `fetch-input` place reads or alignments by default. Add
`--role reference`, `--role targets`, `--role metadata` or `--role truth` for the
other folders, as the GUI file categories do. `./pfvarx input-files` lists the
folders. Both routes need the folder to exist, accept only files whose content belongs
in that folder (below), and never overwrite a file. `add-input` copies the file and leaves the original
unchanged; to use large files in place without copying, register them by path
with `register-inputs`.

```bash
./pfvarx add-input --file /data/run1/sample_R1.fastq.gz --file /data/run1/sample_R2.fastq.gz
./pfvarx add-input --role reference --file /data/ref/Pf3D7.fasta --file /data/ref/core_genome.bed
./pfvarx fetch-input --name sample_R1.fastq.gz --url 'https://example.org/sample_R1.fastq.gz'
```

Replace the example paths, URL and filename with your actual inputs. Add
`--sha256` to `fetch-input` when you have a trusted checksum. For FASTQ input,
validation writes the sample sheet from the paired files in `inputs/fastq/`
if none exists; review it before analysis. `workflow-run` repeats the workflow
check and starts a new analysis only if it passes, exactly as the GUI **Run
workflow** button does. `./pfvarx dry-run` and `./pfvarx pipeline --mode NEW`
remain available as the separate steps.

### Supported input content

PfVarX recognises every input file by its content and structure, never by its name or
extension. A file is accepted when its bytes show that it is the right kind (FASTA
header lines, FASTQ records, BAM/CRAM headers, BED or GFF columns, a sample-sheet,
primer or catalogue header), and refused, with an explanation, when they do not.

| Folder | Accepted content |
| --- | --- |
| inputs/reference/ | FASTA sequence, BED intervals (core mask), GFF3 annotation |
| inputs/fastq/ | FASTQ reads, plain or gzip-compressed |
| inputs/alignments/ | BAM or CRAM alignment (as configured) and its index (.bai, .csi or .crai content) |
| inputs/targets/ | BED intervals (targets) and the primer table |
| inputs/metadata/ | Sample sheet (header with `sample_id`), marker catalogue (header with `marker_id`), other tab-separated tables |
| inputs/truth/ | VCF variants, BED intervals, JSON |
| inputs/partitions/ | partition.json, generated by the benchmark partition command |

Two naming rules remain because tools depend on them. fastp decides from the name
whether reads are compressed, so a gzip-compressed FASTQ must end in `.gz` and a plain
one must not; otherwise fastp reads no records without an error. Reference, BED, GFF,
primer, sample-sheet and catalogue files must be uncompressed, because PfVarX reads
them as plain text.

When a setting in `config/config.yaml` names only a folder (for example
`core_mask_bed: inputs/reference/`), PfVarX uses the one file in that folder whose
content suits the setting. With no such file, or several, validation reports the
problem and lists the candidates. Validation also checks the structure of each file:
FASTA sequences, BED intervals within the contigs, GFF3 columns, contigs and gene/CDS
features, the primer table, the catalogue columns and status, FASTQ record pairing and
quality encoding, and alignment headers, sorting and read groups. The right content
alone does not establish scientific suitability: check that the files match one
reference release.

The Input files categories follow the saved configuration and existing input
folders. After changing input type, save before returning to Input files.
Folders from earlier settings can remain on disk while their tabs are hidden.

Reference Preparation creates a blank sample sheet if it is missing. Fill it or register a reviewed manifest, then rerun preparation. Existing sheets are preserved. To use combined readiness checks, register samples first.
