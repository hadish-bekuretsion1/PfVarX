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
`./pfvarx init` if you want a blank sample-sheet template. GUI saves create no sheet.

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

Inspect `metadata/validation_report.json` and the sample-level statuses. Repair
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
| Input files | `./pfvarx input-files` lists destinations. Copy files there, use `fetch-input`, or register originals by path. |
| Input Files Validation | `./pfvarx validate-input-files` |
| Reference Preparation | `./pfvarx prepare-reference` |
| Workflow Readiness Checks | `./pfvarx workflow-readiness-checks` |

For example, download a file from a direct HTTPS link:

```bash
./pfvarx fetch-input --role fastq --name sample_R1.fastq.gz --url 'https://example.org/sample_R1.fastq.gz'
```

Replace the example URL and filename with your actual input. Add `--sha256`
when you have a trusted checksum. Register the sample inputs before validating.
Use either the separate validation and preparation commands or the combined
readiness command. Both stop on failure. Run `./pfvarx dry-run` to check the
workflow DAG before launching an analysis.

### Supported input filenames

| Folder | Accepted extensions |
| --- | --- |
| inputs/reference/ | .fa, .fasta, .fna, .bed, .gff, .gff3 |
| inputs/fastq/ | .fastq, .fq, .fastq.gz, .fq.gz |
| inputs/alignments/ | .bam and .bai/.csi for BAM. CLI CRAM supports .cram and .crai/.csi. |
| inputs/targets/ | .bed, .tsv |
| inputs/metadata/ | .tsv |
| inputs/truth/ | .truth.vcf, .bed, .json |
| inputs/partitions/ | partition.json, generated by the benchmark partition command |

The GUI file picker filters by extension. Upload and download commands reject
unsupported names. For files copied manually, validation checks registered
input extensions and runs the existing file checks. A valid extension alone
does not establish file integrity or reference compatibility.
