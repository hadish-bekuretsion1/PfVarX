# PfVarX

PfVarX 0.0.0 is a configuration-driven pipeline for paired-end *Plasmodium falciparum*
SNP/small-indel analysis, with a custom caller, QC, cohort integration, annotation,
eligibility-gated downstream methods and one offline scientific dashboard.

This source release contains software, locked environment definitions, templates,
policies and tests. It contains no study reads, reference genome, truth dataset, run
or historical report. Methods remain **scientifically unvalidated and uncalibrated**.
See [scientific status](docs/status.md) before interpreting results.

## 1. Install

You need Linux x86-64, Git, Python 3.11+, Pixi on PATH, and access to the private
GitHub repository. On Windows, use WSL2 with a Linux distribution; macOS and Linux ARM64
are not supported. This example installs the paired-end WGS sample profile.
Run each step in order and resolve any errors before continuing.

```bash
mkdir -p ~/demo
cd ~/demo
git clone https://github.com/hadish-bekuretsion1/PfVarX.git PfVarX-source
python3 PfVarX-source/bootstrap.py --profile sample
cd pfvarx-analysis
./pfvarx version
```

```text
~/demo/
├── PfVarX-source/     # GitHub clone
└── pfvarx-analysis/   # study workspace
```

Bootstrap installs the locked environment and writes configuration templates.
If the source is already cloned, reuse it. If the study workspace is already
initialized, enter it directly. Bootstrap requires a new destination.

| Goal | Bootstrap option |
| --- | --- |
| Alignment/QC or per-sample calls | `--profile sample` |
| Cohort integration | `--profile cohort` |
| Downstream analyses | `--profile downstream --analyses DA003` for diversity, or `--profile downstream --auto-analyses` |
| Benchmarking | `--profile benchmark` |
| Full development environment | `--profile development` |

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

After saving, open **Input files** to upload reads or alignments, or download
them using an HTTPS link and a destination filename. A trusted SHA-256 is
optional. When provided, PfVarX checks it before saving. Otherwise, PfVarX
records the computed checksum and marks the download unverified. Put reference
files in `inputs/reference/`. Register your reads or alignments before running
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

- [Step-by-step HTML quick start](docs/README.html)
- [HTML Operational SOP](docs/sop/PfVarX_SOP_Manual.html) and [GitHub-readable SOP](docs/sop/PfVarX_SOP_Manual.md)
- [Complete orchestration/command guide](docs/orchestration.md)
- [Sample analysis](docs/analysis/README.md) and [benchmarking](docs/benchmarking/README.md)
- [Reporting](docs/reporting/README.md), [deployment](docs/deployment.md), [testing](docs/testing.md)
- [Configuration guide](docs/configuration-guide.md) and [thresholds and their sources](docs/thresholds.md)
- [Caller contract](docs/caller.md) and [sequencing and downstream eligibility](docs/sequencing.md)
- [Architecture](docs/architecture.md), including the CLI and GUI parity rule for contributors,
  and [scientific/operational limits](docs/status.md)

