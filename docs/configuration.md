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
| `inputs/reference/` | Parasite FASTA, core BED and optional host FASTA. Add a matching GFF for explicit FO004 or to enable markers in automatic FO004. |
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
