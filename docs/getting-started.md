# Getting started: install, configure and run PfVarX

This guide takes you from a new computer to a finished run and its report, one step at a
time. Every step can be done in the configuration editor (GUI) or in the terminal; both use
the same files, so you can switch between them at any point.

## 1. What you need

| Item | Requirement |
| --- | --- |
| Computer | Linux x86-64; newly enabled experimental targets: Linux ARM64 and macOS 13+ (Intel/Apple silicon). See platform support in the deployment guide. Windows uses WSL2 with supported Linux. |
| Software | Git, Python 3.11 or later, a web browser for the GUI |
| Network | Internet access to install Pixi and the tools (once per tools environment) |
| Memory | About 8 GB for parasite-only analyses; about 16 GB free with Human DNA removal |
| Disk | 1–3 GB per tools environment (FO001 about 1 GB; analyses with R more), plus your reads, references and results |
| GitHub | Access to the private repository `hadish-bekuretsion1/PfVarX` |

## 2. Install Pixi (once per computer)

Pixi installs PfVarX's tools in locked, reproducible environments. Install it once:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
exec bash                 # reload the shell so that `pixi` is found
pixi --version            # prints the version, for example pixi 0.83.0
```

If `pixi` is not found afterwards, add `~/.pixi/bin` to your `PATH` in `~/.bashrc` and open a
new terminal.

## 3. Get PfVarX

```bash
mkdir -p ~/pfvarx && cd ~/pfvarx
git clone https://github.com/hadish-bekuretsion1/PfVarX.git PfVarX-source
```

For the private repository, log in first with the GitHub CLI (`gh auth login`, then
`gh auth setup-git`), or use a personal access token when Git asks for a password.

## 4. Create a study workspace

A workspace is the folder of one study: its configuration, input files, prepared references
and runs. Create a new one for each study.

```bash
cd ~/pfvarx
python3 PfVarX-source/bootstrap.py --workspace ~/pfvarx/my-study
cd ~/pfvarx/my-study
./pfvarx version
```

This installs only the small **control** environment (configuration, GUI and validation).
The pipeline tools are installed in step 6, for exactly what your configuration needs.

On a machine that will have no internet access later (for example an HPC compute node),
install the tools now instead, for the final output you need:

```bash
python3 PfVarX-source/bootstrap.py --workspace ~/pfvarx/my-study --final-output FO002
```

| Final output | What it produces |
| --- | --- |
| FO001 | Alignment and QC (FASTQ input only) |
| FO002 | Per-sample variant calls |
| FO003 | Cohort VCF across samples |
| FO004 | Annotation and downstream analyses (add `--analyses DA001 DA002 ...` or `--auto-analyses`) |

## 5. Configure the study

### With the GUI

```bash
./pfvarx configure --gui
```

Your browser opens the configuration editor (if not, open the link the terminal prints,
including its token). On tab **1 Input directories** set:

| Setting | What to choose |
| --- | --- |
| project id | Your study name: letters, digits, `_` or `-` |
| input type | FASTQ (paired reads) or BAM (aligned reads) |
| final output | FO001 to FO004, from the table above |
| Assay design | Whole-genome sequencing, Targeted without UMI or Targeted with UMI |
| fasta, Core-genome mask (BED) | Your reference files in `inputs/reference/` |
| Human DNA removal | False, or True (WGS only) with the Human genome reference |
| Reference release label | The exact release of your FASTA and GFF, for example PlasmoDB 68, or Other to type it |
| Study Metadata | Yes if you have a study table (country, dates, population, ...), else No |
| compute | Cores, memory and disk PfVarX may use on this computer |

Then click **Create Configuration and Input Directories**. If the workspace already has a
configuration (bootstrap writes one), tick **Save with backup if existing YAML differs** first.
Only the folders you put input files into are created, for example `inputs/fastq/` and
`inputs/reference/`.

### Or by editing the files

The same settings are in `config/config.yaml`, `config/quality.yaml` and
`config/compute.yaml`. Edit the lines marked `EDIT` or `REVIEW`, leave `FIXED` lines as they
are, then create the folders:

```bash
nano config/config.yaml           # or: code config/config.yaml
./pfvarx input-directories
```

Every setting, its values and an example are listed in the
[configuration guide](configuration-guide.md).

## 6. Install the required tools

Click **Install required tools** under the save area (terminal: `./pfvarx install-tools`).
It installs only what the saved final output, design, input type and analyses need; the
table under the button lists each tool, its version and what it is used for. Change the
configuration later and click it again to add what is newly needed.

```bash
./pfvarx install-tools --check    # are the tools for the saved configuration installed?
```

Workspaces installed with `--final-output` (step 4) already have their tools.

## 7. Add your input files

On tab **2 Input files**, choose a category and upload files from your computer or download
them by HTTPS link. In the terminal:

```bash
./pfvarx add-input --file /data/run1/S1_L001_R1.fastq.gz --file /data/run1/S1_L001_R2.fastq.gz
./pfvarx add-input --role reference --file /data/ref/Pf3D7.fasta --file /data/ref/core_genome.bed
```

| Folder | Put here |
| --- | --- |
| `inputs/fastq/` | Paired reads named like `S1_L001_R1.fastq.gz` and `S1_L001_R2.fastq.gz` |
| `inputs/bam/` | BAM or CRAM files with their index (BAM/CRAM input). The pipeline's own alignments are written to `results/alignments/` |
| `inputs/reference/` | Parasite FASTA, core-genome BED, GFF for FO004, human FASTA with Human DNA removal |
| `inputs/metadata/` | Your study table with a `sample_id` column (Study Metadata: Yes); library catalog for targeted designs |
| `inputs/targets/` | Target BED and primer table (targeted designs) |

## 8. Check, prepare and run

Work through the tabs in order; each names its terminal command.

| Tab | Button | Terminal | What it does |
| --- | --- | --- | --- |
| 3 Input Files Validation | Validate inputs | `./pfvarx validate` | Checks every file; each problem comes with numbered steps to fix it |
| 4 Inputs Preparation | Prepare inputs | `./pfvarx inputs-preparation` | Builds the sample metadata `prepared/samples.tsv` from your files, merges the study table, prepares the reference indexes and checks readiness |
| 5 Workflow Readiness Checks | Run readiness checks | `./pfvarx workflow-readiness-checks` | Repeats every check before a run |
| 6 Workflow Check and Run | Check workflow, then Run workflow | `./pfvarx workflow-check`, `./pfvarx workflow-run` | Plans the run without outputs and shows the expected time, memory and storage of each step, then runs it and opens the live progress page |

Review `prepared/samples.tsv` after step 4: sample and lane IDs come from the file names.
A run continues in the background if you close the browser; reopen tab 6 to see its status,
or **Open run progress** for the progress page: when it started, the steps running now,
the expected finish and the time remaining. When the run ends the page opens the report
dashboard. Expected times come from earlier real-data runs on this computer, so the first
run shows them as unknown.

## 9. Review the results

Each run has its own folder `runs/<RUN_ID>/`. Open its report dashboard in a browser:

```bash
xdg-open runs/<RUN_ID>/results/report/pfvarx_dashboard.html
```

The dashboard is a single file that works offline. At the top of **Overview**, **Where the
results are saved** shows the run folder (with **Copy path** and **Open results folder**),
buttons for each output folder, and every output file with its size and what it is, with
**Open folder** and **Copy path** (nothing is downloaded). Check the QC section and any exclusions before interpreting calls. See the [reporting guide](reporting/README.md).

## 10. Later: update PfVarX or start another study

```bash
cd ~/pfvarx/PfVarX-source && git pull      # newer PfVarX
python3 bootstrap.py --workspace ~/pfvarx/next-study
```

Start a new workspace for a new study or a new PfVarX version; existing workspaces keep the
version they were created with, so their runs stay reproducible.

## 11. When something goes wrong

| Message or problem | What to do |
| --- | --- |
| `pixi: command not found` | Repeat step 2, then open a new terminal |
| The GUI does not open | Copy the link the terminal printed (with its token) into the browser |
| `Save with backup` is asked | Tick **Save with backup if existing YAML differs**; the old file goes to `config/backups/` |
| Create Configuration stays greyed out | Choose the **Reference release label** and fill every field marked with a red `*` |
| `TOOLS_NOT_INSTALLED` | Click **Install required tools** (`./pfvarx install-tools`) |
| Installing tools fails | It needs internet access to the conda channels; see `metadata/install_tools.log` |
| `RESOURCE_MEMORY_EXCEEDED` in readiness | Lower the job memory in `config/compute.yaml`, or use a machine with more memory |
| `HOST_INDEX_REQUIRED` | Human DNA removal needs a bwa-mem2 index you build on a large machine; `./pfvarx prepare` lists the steps |
| A validation problem | Follow the numbered steps it prints; the full record is in `metadata/validation_report.json` |

For every command and option, see the [orchestration guide](orchestration.md); for the full
procedures, the [SOP manual](sop/PfVarX_SOP_Manual.md).

