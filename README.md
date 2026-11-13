# PfVarX

PfVarX 0.0.0 is a configuration-driven pipeline for paired-end *Plasmodium falciparum*
SNP/small-indel analysis, with a custom caller, QC, cohort integration, annotation,
eligibility-gated downstream methods and one offline scientific dashboard.

This source release contains software, locked environment definitions, templates,
policies and tests. It contains no study reads, reference genome, truth dataset, run
or historical report. Methods remain **scientifically unvalidated and uncalibrated**.
See [scientific status](docs/status.md) before interpreting results.

## Clone and initialize only the selected work

Linux x86-64, Python 3.11+ for bootstrap, Git and Pixi are required. An authorized
GitHub login is needed while this repository is private. A clone never executes code
by itself: run the explicit reviewed bootstrap command after cloning.

The source and study workspace sit side by side. For example:

```text
~/demo/
├── PfVarX-source/     # GitHub clone
└── pfvarx-analysis/   # study workspace
```

```bash
mkdir -p ~/demo
cd ~/demo
git clone --filter=blob:none --sparse https://github.com/hadish-bekuretsion1/PfVarX.git PfVarX-source
python3 PfVarX-source/bootstrap.py --profile sample
cd pfvarx-analysis
```

Bootstrap installs the locked environment and writes editable configuration templates.
Choose either configuration route from the study workspace:

- **CLI:** Edit `config/config.yaml`, `config/quality.yaml`, and
  `config/compute.yaml`, then run `./pfvarx init` to create the input folders.
- **Optional GUI:** Run `./pfvarx configure --gui`. Open the local URL printed in
  the terminal if needed. Review your settings and click **Create Configuration
  and Input Directories**.

The [HTML configuration guide](docs/configuration%20guide.html) explains the
settings for both routes. If the study workspace already exists, enter it directly.

| Profile | Scope | Installation examples |
| --- | --- | --- |
| `sample` | FO001/FO002 processing and per-sample calls | `--profile sample` |
| `cohort` | Through FO003 cohort integration | `--profile cohort` |
| `downstream` | FO004, explicit selected DA methods and prerequisites | `--profile downstream --analyses DA003 DA001` |
| `benchmark` | Separate benchmarking purpose; caller-only unless DAs selected | `--profile benchmark`; optional `--rtg` |
| `development` | Full source/tests and all development/method dependencies | `--profile development` |

Add `--input-type BAM|CRAM` for aligned entry. Add `--design TARGETED_NON_UMI` or `--design TARGETED_UMI` for the real targeted assay.
Only those profiles install fgbio. R/MOIRE, relatedness tools, scikit-allel, IQ-TREE,
spatial tools, SnpEff and RTG are selected from the analysis closure, not installed for
a sample-only WGS workspace. Shared validation/policy metadata remains available.

Use `./pfvarx` in the workspace root: it selects the installed locked environment.
Use `./pfvarx tool samtools ...` for a tool in that environment. Plain `pixi run` uses
the full default environment and can install unnecessary packages in a selected workspace.

## Configure → organize → validate → prepare → verify → execute

```bash
# Review actual assay, reference/mask/GFF, outputs, analyses and compute budget.
nano config/config.yaml config/quality.yaml config/compute.yaml

# Optional candidate discovery; inspect identities, library IDs and all metadata.
./pfvarx discover-inputs --directory /absolute/path/to/fastq
nano metadata/discovered_samples.tsv
./pfvarx register-inputs --manifest metadata/discovered_samples.tsv --mode reference

./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline
```

Alternatively, after configuring and registering actual inputs:

```bash
./pfvarx orchestrate --through execute
```

This initializes demand-driven folders, validates, prepares, verifies integrity,
checks the workflow and executes/report/cleans according to retention settings.
It stops at mandatory failures. Use `--through initialize|validate|prepare|verify` to
stop at a preparation boundary. Existing `pipeline --mode NEW|RESUME|EXTEND` remains
available for explicit run management. No scientific estimator is changed by bootstrap.

Inputs are referenced or symlinked, never copied by registration. Explicit HTTPS
acquisition requires a trusted SHA-256; no reference release is downloaded implicitly.
Discovery generates a **review-required** candidate sheet, not verified lab identities.
`verify-prepared` checks original validation, derived hashes, indexes and dictionary
compatibility; execution also enforces it. BAM/CRAM entry does not build unused BWA indexes.

## Outputs and retention

| Mode | Final products |
| --- | --- |
| FO001 | Processed BAM/BAI, QC and dashboard |
| FO002 | Per-sample VCF/TBI, evidence, callability, caller QC and dashboard |
| FO003 | Above plus variant-only cohort VCF/evidence/missingness |
| FO004 | Above plus reference-matched annotation and selected eligible DA outputs |

Existing output paths remain under `runs/<RUN_ID>/results/`; `results` selects the
latest run. Report filename remains `report/pfvarx_dashboard.html`. Missing values,
UNKNOWN markers and scientific skips retain their distinct meanings.

`retention.intermediates: successful` uses dependency-complete Snakemake temporary
cleanup and removes safe remaining `results/work` files after a completed run. Set
`keep` before running for debugging. Explicit cleanup defaults to preview:

```bash
./pfvarx cleanup --run-id YOUR_RUN_ID
./pfvarx cleanup --run-id YOUR_RUN_ID --apply
```

Failed/interrupted runs, originals, references, final results, frozen records, logs
and QC are preserved. Workspace locking prevents concurrent execution/cleanup.

## Manuals and source layout

- [Fifty operational SOPs](docs/sop/PfVarX_SOP_Manual.md)
- [Complete orchestration/command guide](docs/orchestration.md)
- [Sample analysis](docs/analysis/README.md) and [benchmarking](docs/benchmarking/README.md)
- [Reporting](docs/reporting/README.md), [deployment](docs/deployment.md), [testing](docs/testing.md)
- [Architecture](docs/architecture.md) and [scientific/operational limits](docs/status.md)

`bootstrap.py` exports/installs selected workspaces; `pfvarx` selects their environment.
`src/pfvarx/` owns control and science; `workflow/` owns scheduling; `config/templates/`,
`parameters/`, `registry/`, `geo/` own configuration/static resources. `scripts/` owns
packaging/verification/optional image building. `tests/` and `.github/` own verification.
Installed/runtime folders are created on demand and are excluded from source packaging.

Use the [caller procedure](docs/caller.md) to configure and qualify a Layer 2 replacement.

See the [configuration editor and input procedure](docs/configuration.md) and the
[sequencing and downstream eligibility matrix](docs/sequencing.md).

Open the [step-by-step HTML quick start](docs/README.html) after downloading the source.

See the [HTML configuration guide](docs/configuration                   0uide.html) for field-by-field choices and examples.
