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

The editor shows the six workflow steps as numbered tabs, always visible: 1 Input
directories, 2 Input files, 3 Input Files Validation, 4 the preparation step, 5 Workflow
Readiness Checks and 6 Workflow Check and Run (one row on a wide window, rows of three or two
on narrower ones). Each tab names its terminal command. On the first tab every setting is a
row of a form table, grouped in cards (config, quality, compute): the name on the left, with a
red mark when it is required, the input on the right and the help and an example underneath.
**Create Configuration and Input Directories** and **Save with backup** sit together in the
action area below the folder list. Printed or saved as PDF, the page keeps its colours.

Choose the analysis or benchmarking template in a fresh workspace. If YAML exists,
the form loads it. Complete project ID, input type, assay, reference release and
paths. The form shows targeted BED/primer fields and the UMI read structure only
where applicable. Choose an explicit analysis set or automatic selection for FO004.
The form validates the compute envelope, caller options and dependent fields using
the Python backend. The summary shows the YAML-equivalent settings.

In an on-demand workspace (`bootstrap.py` without `--final-output`), click **Install
required tools** after saving. It installs only the locked environment the saved
configuration needs and shows which one; the status line under it says whether the tools for
the current configuration are installed. Inputs Preparation, readiness, the dry run and the
pipeline stop with `TOOLS_NOT_INSTALLED` and this instruction until they are. After changing
the final output, design, input type or analyses, click it again: it installs only the newly
needed environment. Workspaces installed up front (`--final-output`) do not show the button.

The **Tools for this configuration** table under the button lists every tool and package the
form's current choices need: its name, the exact version locked in `pixi.lock`, what PfVarX
uses it for and whether it is installed yet. It follows each change of final output, design,
input type and analyses before you save; FO001 for example lists Python, pysam, Snakemake,
fastp, BWA-MEM2, SAMtools, Picard and mosdepth, and FO002 adds NumPy, SciPy and BCFtools.
`./pfvarx install-tools --check` reports the same list (`tools`) for the saved configuration.

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
| `inputs/reference/` | Parasite FASTA, and core BED. Add a matching GFF for explicit FO004 or to enable markers in automatic FO004. With Human DNA removal (`reference.host_removal: true`, wgs only), also the human genome FASTA. |
| `inputs/fastq/` | Paired FASTQ for FASTQ entry. Register read pairs with their verified sample and lane identities. |
| `inputs/alignments/` | BAM or CRAM with the appropriate index for aligned entry. |
| `inputs/metadata/` | Created only with Study Metadata: Yes, or for a targeted design's library catalog: the study metadata table (tab-separated, with `sample_id`: population, country, latitude, longitude, collection_date, ...), the library catalog YAML and an optional marker catalogue. |
| `inputs/targets/` | Reference-matched target BED and primer TSV for targeted assays. |
| `inputs/truth/`, `inputs/partitions/` | Benchmark truth, confidence regions and partition resources for benchmark workspaces. |

`inputs/` holds only the files you supply to start a run; PfVarX reads them and never
writes there. Everything it generates lives beside `inputs/`, in folders created on demand by
the step that first writes them:

| Folder | Written by | Contents | Safe to delete |
| --- | --- | --- | --- |
| `prepared/` | Inputs Preparation (`./pfvarx inputs-preparation`; `prepare` for the references) | `samples.tsv` (the sample metadata, built from the input files; review it and keep your edits), `reference/` (combined and parasite FASTA, `.fai`, `.dict`, bwa-mem2 index), `snpeff/` database, `manifest.json` sealing them | `reference/` and `snpeff/` yes, `prepare` rebuilds them; keep a reviewed `samples.tsv` |
| `metadata/` | `validate`, `prepare`, `ready` | validation and preparation reports, preparation logs | yes |
| `runs/RUN_ID/` | `./pfvarx pipeline` | frozen configuration and plan, logs, `results/` | only when archiving |
| `results` | `./pfvarx pipeline` | link to the latest run's `results/` | yes |

Workspaces prepared by earlier versions kept the reference derivatives in
`inputs/reference/prepared/`; the next `./pfvarx prepare` moves them to
`prepared/reference/` without rebuilding.

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

## Human DNA removal (optional, whole-genome sequencing only)

Reads from blood or dried blood spots carry human DNA. With `reference.host_removal: true`
(GUI: **Human DNA removal**) and `reference.host_fasta` set to a human primary-assembly
FASTA without ALT contigs, each FASTQ read pair is aligned once with bwa-mem2 to the
parasite and human genomes together. Pairs whose best (primary) alignment of either mate is
on the human genome are removed before duplicate marking; secondary and supplementary hits
on human contigs, and human entries in `SA`/`XA` tags, are removed from the pairs kept. The
BAM keeps the parasite contigs only, so calling, cohort and annotation are unchanged. Each
sample's QC records `host_templates_removed`, `templates_before_host_removal` and
`host_fraction` (gate `G2.host_fraction`, reported, not gated). BAM/CRAM input is used as it
is. With `host_removal: false` (the default) reads align to the parasite genome only and
`host_fasta` is ignored.

| `assay_design` | Human DNA removal |
| --- | --- |
| `wgs` | Eligible (optional) |
| `targeted_no_umi` | Ineligible |
| `targeted_umi` | Ineligible |

Targeted libraries are built from parasite-specific amplicons or probes, so human reads
are a small off-target share that the on-target QC already measures. `host_removal: true`
with a targeted design is refused when the configuration is loaded (CLI and GUI), with the
message to set it to false.

PfVarX does not build the bwa-mem2 index of a human-sized reference: building one is
commonly reported to need about 64 GB of memory, and aligning against it 10–16 GB.

```bash
./pfvarx prepare                       # writes prepared/reference/combined.fasta, then stops:
                                       # HOST_INDEX_REQUIRED, with these steps
# on a machine with enough memory, in a folder holding a copy of combined.fasta:
bwa-mem2 index combined.fasta          # bwa-mem2 2.2.1, as pinned here
# copy combined.fasta.0123 .amb .ann .bwt.2bit.64 .pac back into prepared/reference/
./pfvarx prepare                       # checks the index matches combined.fasta, then seals it
./pfvarx verify-prepared && ./pfvarx ready
```

`ready` adds `layer1.bwa_mem2_host_extra_mem_mb` (12 000 MB, provisional) to the bwa-mem2
request and reports `RESOURCE_MEMORY_EXCEEDED` when `compute.yaml` allows less. Files of
1 GiB or more under `prepared/reference/` are sealed by size and modification time instead
of SHA-256.

The earlier `host_fasta` decoy (retired in November 2026) appended the host genome without
removing reads; a workspace that still has `host_fasta` set keeps aligning to the parasite
genome only until `host_removal` is set to true.

## Add study metadata

Study metadata is optional. Set **Study Metadata: Yes** (`study_metadata: true`) in the
configuration, put the study table (tab-separated, with a `sample_id` column) in
`inputs/metadata/`, and **Inputs Preparation** merges it into the sample metadata
(`prepared/samples.tsv`) it builds from the input files. With **No** it continues with the
sample metadata alone. You can also type study values into `prepared/samples.tsv` as extra
columns, one value per sample, repeated on every lane row. To join a table by hand, or to
preview a join:

```bash
./pfvarx join-metadata --table inputs/metadata/study_metadata.tsv           # report only
./pfvarx join-metadata --table inputs/metadata/study_metadata.tsv --apply   # write
```

Only empty cells are filled, on every lane of the sample. A cell that already holds a
different value is reported and kept. File, lane, library and library-description columns
are never taken from the study table. Sample IDs found in only one of the files are
listed, and coordinates (decimal degrees) and dates (`YYYY-MM-DD`) are checked. The sheet is
copied to `samples.tsv.before-join-<UTC time>` first. With Study Metadata: Yes, **Inputs
Preparation** does this merge automatically, with the same rules.

| Column | Example | Use |
| --- | --- | --- |
| `country` | Ethiopia | DA007 |
| `latitude` | 13.50 | DA007, decimal degrees |
| `longitude` | 39.47 | DA007, decimal degrees |
| `collection_date` | 2025-07-18 | DA007, YYYY-MM-DD |
| `population` | Site_A | Grouped analyses, or use the configured grouping column |
| `age_group`, `study_arm` | Adult, Arm_A | Optional study records |

Review the coordinates included in the dashboard before sharing it.

### Naming standard

Files are recognised by content, so names never decide anything, but a common standard keeps
studies readable. `./pfvarx validate` warns (`NAMING_CONVENTION`, never a failure) when a row
departs from it:

| Item | Standard | Example |
| --- | --- | --- |
| Sample ID | letters, digits, `_`, `-` | `ETH001` |
| Lane | `L` and 3 digits | `L001` |
| Library ID | `<sample_id>_lib<n>` (recommended, not checked) | `ETH001_lib1` |
| FASTQ | `<sample_id>_<lane>_R1.fastq.gz` and `_R2.fastq.gz`; names that start with the sample ID and name the mate (Illumina `ETH001_S1_L001_R1_001.fastq.gz`) are accepted | `ETH001_L001_R1.fastq.gz` |
| Alignment | `<sample_id>.bam` (`.cram`) with its index | `ETH001.bam` |
| Reference | `<genome>_PlasmoDB-<release>.fasta`, `.gff`, `_core.bed` (recommended) | `Pf3D7_PlasmoDB-68.fasta` |
| Panel files | `<panel_id>_v<version>_targets.bed`, `_primers.tsv`, `_probes.bed` (recommended) | `MAD4HATTER_v1_targets.bed` |

### Library records

`assay_design` (`wgs`, `targeted_umi`, `targeted_no_umi`) is the only run setting that
selects the processing branch. What each library is goes in the manifest columns
`technique`, `amplification_route`, `panel_id`, `panel_version`,
`library_protocol_id` and, for imported BAM/CRAM files, `processing_state`; panels and
library protocols are defined in the library catalog (`assay.library_catalog`). The
values, the agreement rules between technique labels and designs, and examples are in
[assay design and library records](sequencing.md).

### Migrating a schema-3 configuration

Configurations written before this release (`schema_version: pfvarx.config/3`) are
refused with an instruction to migrate. `./pfvarx config-migrate` reports the change;
`./pfvarx config-migrate --apply` writes it: `assay.design` becomes `assay_design`,
`sequencing_type` and `amplification_route` become manifest columns, the targets and
primers become the catalog panel `LEGACY_PANEL` (version `config3`), and the UMI read
structure becomes the protocol `LEGACY_UMI_PROTOCOL`. Imported targeted BAM/CRAM rows
get the processing state schema 3 assumed (`primers_trimmed`, or
`consensus_primers_trimmed` for UMI). Only empty manifest cells are filled, comment lines
in the manifest are not kept, and `config.yaml` and the manifest are copied to
`*.before-config4` first. Saving the configuration in the GUI performs the same
migration. Afterwards rename `LEGACY_PANEL` to the real panel name and version if you
like, in the catalog and the manifest together.

The standard folder prefixes shown for the sample registry, reference FASTA,
core BED, GFF and library catalog are accepted in GUI and YAML. They
resolve to the corresponding template filenames. Enter a complete path when
your filename differs. Validation still checks that each required file exists.

The **Input directories** tab saves configuration and creates folders. In
**Input files**, select a file category to upload or download into its folder.
Use the filenames in your saved configuration. Register sample files, then
use **Input Files Validation** followed by **Inputs Preparation** (below). Assay-specific read
processing and SnpEff annotation remain workflow stages. **Workflow Readiness Checks** runs
the checks together.

## Inputs Preparation

GUI tab **Inputs Preparation**; terminal `./pfvarx inputs-preparation [--metadata-table PATH]`.

**Input directories** creates only the folders you put input files into (by upload,
download or browsing): the reads or alignments folder and `inputs/reference/`, plus
`inputs/targets/` for targeted designs, `inputs/truth/` and `inputs/partitions/` for
benchmarking, and the Study Metadata folder (`inputs/metadata/`) only when **Study
Metadata: Yes** or a targeted design needs its library catalog there. Everything else
(`prepared/`, `metadata/`, `runs/` and their subfolders) is created on demand by the step
that writes it, from Inputs Preparation to the end of the pipeline.

Two settings in the form control the metadata:

| Setting (GUI) | YAML | Values |
| --- | --- | --- |
| Study Metadata | `study_metadata` | No (`false`, default) or Yes (`true`) |
| Study Metadata folder (shown when Yes) | `study_metadata_dir` | `inputs/metadata/` |

Inputs Preparation runs four stages in order and stops at the first failure:

1. **Sample metadata** (`prepare-samples`). When `prepared/samples.tsv` (the fixed
   `samples_tsv` setting) does not exist, it is built from the input files, with its folder:
   paired FASTQ names give one row per sample and lane, and `library_id` defaults to
   `sample_id`, so review them. BAM/CRAM names do not give sample identities, so a blank
   template is written and preparation stops until it is filled in. An existing file is used
   as it is and never overwritten.
2. **Study metadata** (`prepare-samples`).

   | Study Metadata | Study Metadata folder holds | Inputs Preparation |
   | --- | --- | --- |
   | No | (not created; ignored if present) | continues with the sample metadata only (`NOT_USED`) |
   | Yes | one study table (tab-separated, with `sample_id`) | merges it into the sample metadata (`MERGED`): empty cells are filled, existing values are never overwritten, differences and unmatched IDs are listed |
   | Yes | no study table yet | continues with the sample metadata only and warns (`NO_STUDY_TABLE`); add the table and run it again |
   | Yes | several study tables | stops and lists them; keep one there or choose it with `--metadata-table PATH` |

   The library catalog, the marker catalogue, the sample metadata and earlier merge
   backups are never taken for study tables. Repeating the merge writes nothing new. The
   result is in `metadata/inputs_preparation.json`.
3. **References** (`prepare`). Revalidates the inputs and builds what the input files and
   final output need under `prepared/reference/`: the FASTA index and sequence dictionary,
   and the BWA-MEM2 index when reads are aligned (FASTQ input and UMI BAM input). Current
   files are reused.
4. **Checks** (`verify-prepared`, `ready`).

## Terminal workflow

From your study workspace, edit the three configuration YAML files and run
the corresponding commands for the GUI stages.

| GUI stage | CLI equivalent |
| --- | --- |
| Input directories | `./pfvarx input-directories` |
| Input directories: Install required tools | `./pfvarx install-tools` (`--check` reports only) |
| Input files: upload from your computer | `./pfvarx add-input --file /path/to/file` (repeat `--file` for several) |
| Input files: download from an HTTPS link | `./pfvarx fetch-input --url URL --name FILENAME` |
| Input Files Validation | `./pfvarx validate-input-files` |
| Inputs Preparation | `./pfvarx inputs-preparation` |
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
Inputs Preparation (or validation, if it runs first) writes the sample metadata
`prepared/samples.tsv` from the paired files in `inputs/fastq/` if none exists; review it
before analysis. `workflow-run` repeats the workflow
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
| inputs/targets/ | BED intervals (targets, probes) and the primer table, named by a panel in the library catalog |
| inputs/metadata/ | Study metadata table (header with `sample_id`), marker catalogue (header with `marker_id`), library catalog (YAML with `schema_version: pfvarx.library_catalog/1`), other tab-separated tables |
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

Inputs Preparation builds a missing sample sheet from the input files (a blank template for BAM/CRAM input). Fill or review it, or register a reviewed manifest, then rerun preparation. Existing sheets are preserved.
