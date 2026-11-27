# PfVarX

PfVarX 0.0.0 is a configuration-driven pipeline for paired-end *Plasmodium falciparum*
SNP/small-indel analysis, with a custom caller, QC, cohort integration, annotation,
eligibility-gated downstream methods and one offline scientific dashboard.

This source release contains software, locked environment definitions, templates,
policies and tests. It contains no study reads, reference genome, truth dataset, run
or historical report. Methods remain **scientifically unvalidated and uncalibrated**.
See [scientific status](docs/status.md) before interpreting results.

## 1. Install

New to PfVarX? [Getting started](docs/getting-started.md) walks through every step below,
from installing Pixi to opening the report, with a troubleshooting table.

Native installation targets are Linux x86-64, Linux ARM64, and macOS 13+ on Intel or
Apple silicon. You need Git, Python 3.11+ for bootstrap, and Pixi on PATH. On Windows, use WSL2 with a supported Linux distribution.
macOS and Linux ARM64 are newly enabled, experimental targets pending successful native
CI runs; a solved lock is not runtime or scientific validation. See
[platform support](docs/platforms.md). This example creates a workspace whose tools are installed on demand.
Run each step in order and resolve any errors before continuing.

```bash
mkdir -p ~/demo
cd ~/demo
git clone https://github.com/hadish-bekuretsion1/PfVarX.git PfVarX-source
python3 PfVarX-source/bootstrap.py
cd pfvarx-analysis
./pfvarx version
```

```text
~/demo/
├── PfVarX-source/     # GitHub clone
└── pfvarx-analysis/   # study workspace
```

Bootstrap installs the locked control environment (or, with `--final-output`, the tools too)
and writes configuration templates.
If the source is already cloned, reuse it. If the study workspace is already
initialized, enter it directly. Bootstrap requires a new destination.

The tools are installed for the final output you choose, in one of two ways:

- **On demand (default).** `bootstrap.py` without options installs only a small control
  environment (configuration, GUI, validation; no pipeline tools). After saving the
  configuration, click **Install required tools** on the Input directories tab (or run
  `./pfvarx install-tools`): it installs exactly what the saved final output, design, input
  type and analyses need. The **Tools for this configuration** table under it lists each tool
  and package your current choices need, with its locked version, what it is used for and
  whether it is installed, and follows every change in the form. Change the configuration
  later and click it again to add what is newly needed. Steps that run tools stop with that instruction until the tools are
  installed. This needs internet access (or a filled Pixi cache) at that step.
- **Up front.** `--final-output FO00N` installs the tools at bootstrap, for machines without
  internet access later, such as HPC compute nodes. A later final output then needs a new
  workspace; the configuration and the GUI refuse one above the installed final output and
  name the bootstrap option.

Either way each final output gets only the tools it uses:

| Final output | Bootstrap option (up front) | Installed beyond the shared base* |
| --- | --- | --- |
| FO001 alignment and QC | `--final-output FO001` | fastp, bwa-mem2 (and fgbio for targeted designs, as `fgbio-minimal` without R). No variant caller, bcftools, numpy, scipy, numba or R |
| FO002 per-sample calls | `--final-output FO002` | adds the variant caller (numpy, scipy, numba) and bcftools |
| FO003 cohort | `--final-output FO003` | same tools as FO002, plus the cohort step |
| FO004 annotation and analyses | `--final-output FO004 --analyses DA003` (or `--auto-analyses`) | adds SnpEff and only the requested analyses' tools: R with MOIRE (DA001), Dcifer and hmmIBD (DA004), scikit-allel (DA005), IQ-TREE (DA006), mgcv (DA007); DA002 and DA008 need no extra tool |
| Benchmarking | `--profile benchmark` | as FO002, plus RTG Tools with `--rtg` |
| Full development environment | `--profile development` | everything, plus the test tools |

\* Snakemake, samtools, Picard (`picard-slim`, without R) and mosdepth, with the Python,
pysam and htslib of the control environment: every final output uses them. To see what each installation gets, run
`pixi run python scripts/check_installs.py` ([details](docs/deployment.md#check-what-an-installation-contains)). The earlier `--profile sample`, `--profile cohort` and
`--profile downstream` options still work and match FO002, FO003 and FO004.

For aligned input, add `--input-type BAM` or `--input-type CRAM`. For targeted
libraries, add `--design targeted_no_umi` or `--design targeted_umi`. The design
(`assay_design` in `config.yaml`) selects the processing branch; each library's
technique, amplification route, panel, protocol and imported-BAM state are recorded in the
sample manifest and library catalog ([assay design and library records](docs/sequencing.md)).

## 2. Configure

Run the remaining commands from the study workspace. Set the actual assay,
reference FASTA, core-genome mask BED, reference release and compute budget. The core
mask is the accessible core of the genome (the Core regions of Miles et al. 2016, as used
by MalariaGEN), not a gene list.
Choose one configuration route:

- **CLI:** Edit `config/config.yaml`, `config/quality.yaml`, and
  `config/compute.yaml`. Save the files, then run `./pfvarx init`.
- **Optional GUI:** Run `./pfvarx configure --gui`. If needed, open the local
  URL printed in the terminal. Review your settings and click **Create
  Configuration and Input Directories**.

Either route creates the input folders after configuration. See the
[Configuration guide](docs/configuration-guide.md) for field examples.

## 3. Register inputs

For conventionally named paired FASTQs, replace the example directory and run:

```bash
./pfvarx discover-inputs --directory /absolute/path/to/fastq
nano metadata/discovered_samples.tsv
./pfvarx register-inputs --manifest metadata/discovered_samples.tsv --mode reference
```

Confirm sample and library IDs and R1/R2 pairing before registering. For BAM,
CRAM or other naming conventions, use a reviewed TSV manifest instead.
Reference mode records paths to the original files.

After saving, open **Input files** to upload reads (`inputs/fastq/`) or BAM/CRAM files (`inputs/bam/`), or download
them using an HTTPS link and a destination filename. A trusted SHA-256 is
optional. When provided, PfVarX checks it before saving. Otherwise, PfVarX
records the computed checksum and marks the download unverified. Put reference
files in `inputs/reference/`. Register your reads or BAM/CRAM files before running
**Workflow Readiness Checks → Run readiness checks**. The tab runs `prepare-samples`,
`validate`, `prepare`, `verify-prepared`, and `ready` in that order and stops on failure.

Study metadata is optional. Set **Study Metadata: Yes** (`study_metadata: true`) if you
have a study table with a `sample_id` column (country, latitude, longitude,
collection_date, population, ...); Input directories then creates `inputs/metadata/` for
it. **Inputs Preparation** (`./pfvarx inputs-preparation`) builds the sample metadata
`prepared/samples.tsv` from your input files and merges the study table into it. With
Study Metadata: No it continues with the sample metadata only. See
[Inputs Preparation](docs/configuration.md#inputs-preparation).

## 4. Check and run

The [terminal workflow](docs/configuration.md#terminal-workflow) lists CLI
equivalents for the GUI tabs.

```bash
./pfvarx validate &&
./pfvarx inputs-preparation &&
./pfvarx dry-run
```

If validation fails, its `summary` names each problem, why it failed and the numbered
steps that fix it; the full record is in `metadata/validation_report.json`. When
readiness is true and the dry run reports `OK`, start the analysis:

```bash
./pfvarx pipeline --mode NEW
```

`./pfvarx workflow-run` performs the dry run and the start in one step. In the GUI,
the **Workflow Check and Run** tab does the same. Every GUI tab has an equivalent
command; see the [terminal workflow](docs/configuration.md#terminal-workflow).

## 5. Review results

Record the run ID and open `results/report/pfvarx_dashboard.html`. Check sample
QC and exclusion reasons before interpreting calls. If no sample passes QC, the run
still completes and `results/cohort/cohort_status.json` explains which samples were
blocked and why. Each run retains its
outputs, logs and frozen configuration under `runs/<RUN_ID>/`.

| Mode | Final products |
| --- | --- |
| FO001 | FASTQ entry only. Processed BAM/BAI, QC and dashboard |
| FO002 | Per-sample VCF/TBI, evidence, callability, caller QC and dashboard |
| FO003 | Above plus variant-only cohort VCF/evidence/missingness |
| FO004 | Cohort outputs, SnpEff annotation from the configured GFF, and the requested downstream analyses that are eligible. DA001 (COI) runs when requested or when another analysis needs it. |

## Guides

- [Getting started: install Pixi, install, configure and run](docs/getting-started.md) ([HTML](docs/getting-started.html))
- [Step-by-step HTML quick start](docs/README.html)
- [HTML Operational SOP](docs/sop/PfVarX_SOP_Manual.html) and [GitHub-readable SOP](docs/sop/PfVarX_SOP_Manual.md)
- [Complete orchestration/command guide](docs/orchestration.md)
- [Sample analysis](docs/analysis/README.md) and [benchmarking](docs/benchmarking/README.md); step-by-step HTML guides
  with commands and expected output: [sample analysis](docs/PfVarX_Sample_Analysis_Guide.html) and
  [benchmarking](docs/PfVarX_Benchmarking_Guide.html)
- [Reporting](docs/reporting/README.md), [deployment](docs/deployment.md), [testing](docs/testing.md)
- [Configuration guide](docs/configuration-guide.md) and [thresholds and their sources](docs/thresholds.md);
  every threshold, setting and compute request is in a configuration file, none in the code
  ([where values are kept](docs/thresholds.md#where-values-are-kept))
- [Caller contract](docs/caller.md) and [sequencing and downstream eligibility](docs/sequencing.md)
- [Architecture](docs/architecture.md), including the CLI and GUI parity rule for contributors,
  and [scientific/operational limits](docs/status.md)


