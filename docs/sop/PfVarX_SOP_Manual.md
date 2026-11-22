# PfVarX 0.0.0 — Operational SOP manual

**Document:** PFVARX-SOP-001 · **Revision:** 3.2 · **Prepared:** 12 November 2026.
**Audience:** Linux bioinformatics researchers, developers and benchmark analysts.
**Authority:** finalized local 0.0.0 implementation, workflow rules, parameter/policy
registries and the included operating manuals. This is an operational companion, not
a change to the software, calibrated settings or scientific claims. A laboratory must
approve its study-specific validation protocol before routine surveillance decisions.

Read this manual in order within the relevant part. Part I is for maintainers; routine
users start with Part II after the common setup. Part III uses a separate workspace.
Every SOP has the same seven sections so an operator can record a formal completion.

## How to use commands safely and consistently

Run commands in Bash. Stop when a command fails; inspect its error before continuing.
Commands below are issued from a named **workspace root** containing `pixi.toml`,
`src/`, `config/` and `workflow/`. Do not execute jobs from `src/` or a data directory.
All paths containing `$...` are shell variables set below; replace any UPPERCASE
sample/run identifier or `/absolute/path/...` resource with its real verified value.
Examples are not biological data. A command block is executable only after its
prerequisites have been satisfied. `nano` is an editor example; use your usual editor.

Set these variables once in a terminal. If the SOP bundle is moved, update `SOP_HOME`.
Set the source path to the reviewed checkout or a fresh source-archive extraction.

```bash
export PFVARX_SOURCE=/absolute/path/to/PfVarX-source
export SOP_HOME="$PFVARX_SOURCE/docs/sop"
export LAB="$HOME/PfVarX-lab"
export DEV="$LAB/development"
export ANALYSIS="$LAB/analysis"
export BENCH="$LAB/benchmarking"
export RECORDS="$LAB/records"
BOOTSTRAP_ARGS=(--final-output FO002 --design WGS --input-type FASTQ)
mkdir -p "$LAB" "$RECORDS"
test -f "$PFVARX_SOURCE/workflow/Snakefile"
test -f "$SOP_HOME/PfVarX_SOP_Manual.md"
```

Keep pristine source read-only in practice. The commands create workspaces only on
request. Do not copy `.pixi/` between locations; install separately in each workspace.
Use a full-source distribution or selected bootstrap export, not a wheel alone.
`./pfvarx` selects the profile environment; `./pfvarx tool TOOL` runs a utility there.
Plain `pixi run` installs the full default and is reserved here for full/development
checks. Set BOOTSTRAP_ARGS before R01 for the actual profile/design/input type.
For FO004 use `--final-output FO004 --analyses ...` with the exact requested methods;
for aligned entry add `--input-type BAM|CRAM`. Export a new workspace for a changed
capability set. Browser QA dependencies are separate.
Unix utilities used include Bash, Git, `find`, `sha256sum`, `gzip`, `awk`, `sed`, `grep`,
`tee`, `df`, `free`, `nproc` and a terminal editor. Confirm utilities are present.
Use `set -o pipefail` before pipelines so an upstream failure is not hidden by `tee`.

## Procedures included in this manual

| Supported procedure | Required input/state |
| --- | --- |
| Paired-end WGS FASTQ processing | Direct or actual SWGA library, matching reference and core mask |
| Targeted non-UMI / UMI processing | Real target/primer definitions; genuine inline UMI structure where used |
| WGS BAM and CRAM input | Valid reference identity, coordinate sort, index and read groups |
| Existing targeted aligned input | Genuine completed primer clipping; UMI consensus provenance where applicable |
| FO001–FO004 outputs | Mode-appropriate validation, annotation and scientific eligibility |
| DA001–DA007 | Only registry-permitted designs, scales, metadata and prerequisites |
| Built-in benchmarking | Declared pure/mixture truth, confidence BED, frozen partition, completed run |
| Resource, provenance and offline report operations | Existing measured artifacts and immutable run definitions |
| Development/test/package/release operations | Existing implementation and actual tests; isolated generated work |

Only the implemented Pixi deployment and supported input routes are documented as
execution procedures. Existing VCF/joint-VCF entry, container deployment, external
caller automation and proposed methods are omitted. CRAM and targeted aligned entry
have implemented contracts but less extensive verification than WGS BAM; validate the
actual input carefully and record the evidence. No procedure here establishes field
accuracy or clinical validation by software completion alone.

**Biological handling rule:** use the implemented haploid dominant-allele representation;
retain WSAF and uncertainty for mixtures. GT=0 at a retained low-frequency candidate
can mean a reference dominant allele while an alternative allele is present. It does
not make the alternative record a reference call. Missing VCF positions are not known
reference genotypes. UNKNOWN marker states are not wild type. Technical DNA-mixture K,
MOIRE biological COI and FWS are distinct quantities. AT-rich/low-complexity regions,
homopolymers, amplification bias and uneven target coverage require stratified review.
Do not apply human diploid defaults or generic low-complexity filtering automatically.

## Operator record required for every SOP

In `$RECORDS`, keep an execution sheet with SOP ID, operator, date, workspace, source
commit/version, input checksums, exact commands, exit codes, run ID, deviations, output
paths, QC decision, exclusion reason and reviewer sign-off. Example header:

```bash
printf 'sop_id\toperator\tdate\tworkspace\trun_id\tstatus\tdeviation\treviewer\n' > "$RECORDS/execution-template.tsv"
```

This creates a blank form, not evidence that any procedure ran. Do not log tokens,
passwords or the complete environment. For JSON commands keep stdout separate from
stderr. `pipeline` returns nonzero for incomplete/failure outcomes; QC exclusion and
scientific ineligibility must be reviewed independently of process completion.

## Procedure index

- [D01 — Initialize an isolated development checkout](#procedure-d01)
- [D02 — Install the locked Pixi environment](#procedure-d02)
- [D03 — Maintain configuration schemas and parameters](#procedure-d03)
- [D04 — Maintain validators and select entry routes](#procedure-d04)
- [D05 — Develop FASTQ preprocessing and QC](#procedure-d05)
- [D06 — Develop reference and annotation preparation](#procedure-d06)
- [D07 — Maintain alignment and BAM/CRAM processing](#procedure-d07)
- [D08 — Develop the custom SNP and indel caller](#procedure-d08)
- [D09 — Maintain filtering, status classes and callability](#procedure-d09)
- [D10 — Develop diversity, local haplotypes and clonality](#procedure-d10)
- [D11 — Maintain resistance marker analysis](#procedure-d11)
- [D12 — Maintain population structure, relatedness, trees and spatial analysis](#procedure-d12)
- [D13 — Integrate and execute Snakemake changes](#procedure-d13)
- [D14 — Maintain diagnostics and reproducibility](#procedure-d14)
- [D15 — Maintain the offline reporting dashboard](#procedure-d15)
- [D16 — Test, package and prepare a release](#procedure-d16)
- [R01 — Create and install a routine-analysis workspace](#procedure-r01)
- [R02 — Register one sample or a cohort](#procedure-r02)
- [R03 — Register and prepare the exact reference resources](#procedure-r03)
- [R04 — Preflight resources, validation and the workflow](#procedure-r04)
- [R05 — Execute paired-end WGS FASTQ analysis](#procedure-r05)
- [R06 — Execute targeted non-UMI FASTQ analysis](#procedure-r06)
- [R07 — Execute targeted UMI FASTQ analysis](#procedure-r07)
- [R08 — Enter an existing WGS BAM](#procedure-r08)
- [R09 — Qualify existing targeted alignments](#procedure-r09)
- [R10 — Qualify existing CRAM input](#procedure-r10)
- [R11 — Select one or multiple downstream analyses before launch](#procedure-r11)
- [R12 — Validate every final-output mode and review QC](#procedure-r12)
- [R13 — Diagnose, resume or extend safely](#procedure-r13)
- [R14 — Operate and export the reporting dashboard](#procedure-r14)
- [R15 — Archive a study and maintain the software boundary](#procedure-r15)
- [R16 — Discover and register inputs without copying originals](#procedure-r16)
- [R17 — Verify prepared artifacts before execution](#procedure-r17)
- [R18 — Orchestrate the complete chain and apply retention safely](#procedure-r18)
- [B01 — Create a separate benchmark workspace and select inputs](#procedure-b01)
- [B02 — Run the supported synthetic software smoke test](#procedure-b02)
- [B03 — Prepare independent truth and confidence regions](#procedure-b03)
- [B04 — Freeze the protocol and lock the development/held-out partition](#procedure-b04)
- [B05 — Execute the benchmark pipeline and score development once](#procedure-b05)
- [B06 — Review SNP detection accuracy](#procedure-b06)
- [B07 — Review indel accuracy and complex-region limitations](#procedure-b07)
- [B08 — Review low-frequency detection and WSAF calibration](#procedure-b08)
- [B09 — Evaluate technical mixture K agreement](#procedure-b09)
- [B10 — Review implemented diversity and clonality outputs in benchmarks](#procedure-b10)
- [B11 — Audit exact matching, confidence and callable denominators](#procedure-b11)
- [B12 — Review built-in summary metrics and uncertainty](#procedure-b12)
- [B13 — Measure resources and generate the supported profile](#procedure-b13)
- [B14 — Measure sample-count scaling using ordinary supported runs](#procedure-b14)
- [B15 — Verify repeatability and frozen identity](#procedure-b15)
- [B16 — Perform one held-out evaluation and archive the experiment](#procedure-b16)

# Section 4 — Unified downstream eligibility and output matrix

**E** means eligible under the listed baseline requirements. **CE** means an additional analysis-specific condition must be met. **NE** means the current method excludes the category or the workflow lacks that route. Read the clonality and sequencing columns separately. Combine both with the input, metadata, QC and evidence gates before interpreting any result. Every analysis records a `status_manifest.json` in its output directory.

The table abbreviates `runs/RUN_ID/results/` as **R**. The cohort evidence inputs below are **M** = `R/cohort/evidence_matrix.tsv.gz`, **S** = `R/cohort/sites.tsv`, and **C** = `R/cohort/cohort_samples.txt`. **Q** = `R/qc/samples/SAMPLE/qc_decision.json`. **L** = `R/downstream/DA001/clonality.tsv`. All methods share the upstream admission and sample identity checks. Under automatic selection, an unset GFF omits annotation and blocks marker and spatial outcomes.

| DA identifier | Downstream analysis name | Scientific objective | Required input files | Required metadata | Input file formats | Monoclonal eligibility | Polyclonal eligibility | WGS eligibility | sWGS eligibility | Targeted non-UMI eligibility | Targeted UMI eligibility | Eligibility conditions and scientific justification | Expected output files | Output formats | Output directory paths | Output descriptions and interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DA001 | COI and clonal status | Estimate infection multiplicity and assign a supported clonal category | M, S, C, Q | Sample ID | TSV.gz, TSV, text, JSON | E | E | CE | CE | CE | CE | MOIRE uses observed alleles at thinned informative SNP loci. Depth or panel informativeness can leave too few qualifying loci. sWGS uses the WGS gate. UMI evidence must first pass the consensus route. Previous clonality is not required. | `clonality.tsv`, `moire_coi.tsv`, `moire_allele_freqs.tsv` | TSV | R/downstream/DA001/ | Posterior summaries and clonal categories. Convergence is not assessed by this wrapper. Technical caller K and WSAF cannot establish clone count or clone proportions. |
| DA002 | Drug-resistance markers | Summarize callable codon states and marker prevalence | `R/annotation/cohort.ann.vcf.gz`, reference FASTA/GFF, versioned catalogue, `R/alignments/SAMPLE.bam`, Q | Sample ID; population for grouped summaries | VCF.gz, FASTA, GFF3, TSV, BAM, JSON | E | E | CE | CE | CE | CE | Require a verified matching catalogue and read-backed callable codons. Shallow depth increases UNKNOWN states. Targeted panels must cover each codon, with qualified UMI consensus evidence where applicable. | `marker_calls.tsv`, `prevalence.tsv`, `catalog_provenance.json` | TSV, JSON | R/downstream/DA002/ | Callable denominator and marker states, including mixed states. UNKNOWN does not mean wild type. This is a molecular marker outcome. |
| DA003 | Genetic diversity | Describe FWS and diversity over the analysed callable scope | M, S, C, Q, `R/variants/samples/SAMPLE.callable.bed` | Sample ID; group labels for group analyses | TSV.gz, TSV, text, JSON, BED | E | E | E | CE | CE | CE | Baseline callable coverage and sufficient variation are required. sWGS must retain interpretable callable denominators. Targeted results describe the panel and require panel-specific interpretation. Group sizes constrain group analyses. | `within_host_fws.tsv`, `population_diversity.tsv` | TSV | R/downstream/DA003/ | Descriptive FWS, heterozygosity and diversity. FWS does not assign COI. |
| DA004 | Relatedness and IBD | Estimate genetic relatedness for eligible pairs | M, S, C, Q, L, `R/eligibility/cohort_decisions.json` | Sample ID | TSV.gz, TSV, text, JSON | CE | CE | CE | CE | CE | CE | Require assessed clonality, at least one eligible pair and adequate shared informative markers. hmmIBD is used for monoclonal pairs. Dcifer uses COI for mixed-infection pairs. Shallow or sparse panel data can fail the shared-marker gate: WGS pairs need 200 informative markers, targeted pairs 20, and pairs below 200 are reported as low precision without a relatedness label. | `pairwise_relatedness.tsv` | TSV | R/downstream/DA004/ | Pair estimates, method and support. Tool estimands differ. Genetic relatedness does not prove direct transmission. |
| DA005 | Population structure | Describe haploid PCA and Hudson FST across qualifying groups | M, S, C, Q, L, `R/eligibility/cohort_decisions.json` | Population labels and sample ID | TSV.gz, TSV, text, JSON | CE | CE | CE | CE | NE | NE | Registry restricts this method to the WGS branch. Require assessed clonality, enough complete variable SNPs and qualifying samples or groups. sWGS must meet the same completeness gate. Mixed samples use dominant alleles and can shift group summaries. | `pca.tsv`, `pca_model.json`, `fst.tsv` | TSV, JSON | R/downstream/DA005/ | PCA positions and eligible FST estimates. Dominant alleles represent one approximation for mixed infections. |
| DA006 | Single-lineage phylogeny | Infer a tree from qualifying single-lineage samples | M, S, C, Q, L, `R/eligibility/cohort_decisions.json`, callable BEDs | Sample ID | TSV.gz, TSV, text, JSON, BED | CE | NE | CE | CE | NE | NE | Requires DA001 monoclonal status, WGS branch, adequate sample count, low missingness and informative SNPs. The alignment represents one lineage per sample. Mixed infections and targeted panels are excluded by the implemented method. | `alignment.fasta`, `tree.treefile`, `runtime_gates.json` | FASTA, Newick, JSON | R/downstream/DA006/ | Ascertainment-corrected tree and support when the gates pass. Branches do not establish transmission history. |
| DA007 | Spatial and temporal surveillance | Describe marker outcomes over supplied places and dates | `R/downstream/DA002/status_manifest.json`, `R/downstream/DA002/marker_calls.tsv`, sample sheet, Q | Country, latitude, longitude, collection date | JSON, TSV | E | E | CE | CE | CE | CE | Require completed DA002 outcomes and valid geography/date. Every sequencing branch inherits DA002 callable denominators. GAM fitting additionally depends on support in the supplied data. | `country_month_summary.tsv`, `grid_summary.json`, conditional `gam_fit.tsv` | TSV, JSON | R/downstream/DA007/ | Descriptive summaries and optional model output. Review ascertainment and geographic sampling bias. |

TBS with UMI follows the targeted UMI column when its FASTQs, target BED, primer TSV and UMI layout match that route. TBS without UMI follows targeted non-UMI under the corresponding contract. The implementation has no barcode-table entry, independent sWGS or TBS preset, or validated haplotype-level clone-proportion output. Consult the actual gate and terminal manifest for each run. E and CE describe eligibility, while execution and empirical evaluation have separate status.

# PART I — Development and implementation

These procedures maintain the existing implementation. They do not instruct a researcher
to replace the working caller with placeholders or a different variant caller. D01
creates the entire real directory/module tree by copying source; subsequent SOPs identify
which existing files to edit, how to test a change, and when to reject it.


<a name="procedure-d01"></a>

## D01 — Initialize an isolated development checkout

### 1. Procedure identification

**Module:** repository and packaging. **Workflow:** development.

**Purpose and scope:** Create the actual source/module/workflow tree without bringing study outputs into development.

### 2. Prerequisites and requirements

Common setup completed; finalized source available; Git and Python 3 available. `$DEV` must not exist. Linux x86-64 is the supported host.

### 3. Detailed operational procedure

**Step 1. Create the development directory tree**

```bash
test ! -e "$DEV"
cp -a "$PFVARX_SOURCE" "$DEV"
cd "$DEV"
python3 scripts/package.py --check
find src/pfvarx workflow config tests -type f | sort
```

**Step 2. Initialize local development history if this is a source extraction**

```bash
git init -b main
git config user.name "YOUR_REVIEWED_COMMIT_NAME"
git config user.email "YOUR_VERIFIED_COMMIT_EMAIL"
git add .
git diff --cached --stat
git commit -m "Import verified PfVarX source"
git switch -c codex/development
```

### 4. Expected outputs

`$DEV/{src/pfvarx,workflow/rules,workflow/scripts,config/templates,config/parameters,config/registry,config/geo,tests,docs,scripts}` and local Git metadata. No `inputs`, `runs` or `results` initially.

### 5. Verification and quality control

Before Git initialization, the source check must report PASS and zero runtime files. Inspect staged files; raw reads, reference genomes, caches and reports must be absent.

### 6. Error handling and troubleshooting

Existing `$DEV`: stop and choose a new directory; do not copy inside it. Missing source components: re-extract the reviewed source archive. Do not copy a full analysis checkout as a substitute.

### 7. Completion criteria

Record source archive SHA-256 or reviewed commit, branch and clean initial tree. Continue to D02 only when the imported tree is complete.


<a name="procedure-d02"></a>

## D02 — Install the locked Pixi environment

### 1. Procedure identification

**Module:** Pixi and declared Bioconda/Conda dependencies. **Workflow:** development.

**Purpose and scope:** Install the release-supported dependency environment without altering the lockfile.

### 2. Prerequisites and requirements

D01 complete; Linux x86-64; Pixi available; network or cached packages; disk sufficient for both installation and separate test data. The full default environment is required for all implemented methods.

### 3. Detailed operational procedure

**Step 1. Inspect hardware and install**

```bash
cd "$DEV"
pixi --version
uname -m
nproc
free -h
df -h "$DEV" /tmp
pixi install --locked
pixi run pfvarx version
pixi run pfvarx --help
```

**Step 2. Check executables and installed R libraries**

```bash
pixi run samtools --version
pixi run bcftools --version
pixi run snakemake --version
pixi run Rscript -e 'for (p in c("moire","dcifer","mgcv","jsonlite")) { stopifnot(requireNamespace(p, quietly=TRUE)); cat(p, as.character(packageVersion(p)), "\n") }'
```

### 4. Expected outputs

Location-specific `.pixi/envs/default/` and installed tools specified by pixi.toml/pixi.lock. Save version outputs in the execution record; no environment is included in source distribution.

### 5. Verification and quality control

Locked installation and CLI checks must succeed. Confirm actual RAM and disk availability, not just planner estimates. Version declarations must remain consistent.

### 6. Error handling and troubleshooting

Outdated lock: maintainers review the dependency change and matching lock; routine users do not regenerate it to bypass an error. Disk quota: relocate a new workspace or remove only disposable test work. Never move a live installed prefix.

### 7. Completion criteria

Full default Pixi environment installed, essential executables/packages present, and source/lock identity recorded.


<a name="procedure-d03"></a>

## D03 — Maintain configuration schemas and parameters

### 1. Procedure identification

**Module:** configuration/models.py, loader.py, config/. **Workflow:** development.

**Purpose and scope:** Develop validated user configuration while preserving protected scientific policy.

### 2. Prerequisites and requirements

D02 complete; current schema, tests/unit/test_configuration.py and template files. Change user schema only with matching validation and migration/documentation changes.

### 3. Detailed operational procedure

**Step 1. Inspect authoritative definitions**

```bash
cd "$DEV"
sed -n '1,240p' src/pfvarx/configuration/models.py
cat config/templates/analysis.yaml config/templates/compute.yaml
cat config/parameters/layer2.yaml
```

**Step 2. Edit schema, templates and GUI help together**

The CLI (YAML templates) and the configuration GUI configure the same settings, so a
setting change touches both in the same commit (see "CLI and GUI parity" in
`docs/architecture.md`). Setting help lives only in `gui/field-guide.json`; regenerate the
configuration guide from it.

```bash
nano src/pfvarx/configuration/models.py config/templates/analysis.yaml config/templates/benchmarking.yaml
nano src/pfvarx/configuration/gui/field-guide.json
python scripts/configuration_guide.py
nano tests/unit/test_configuration.py
pixi run pytest -q tests/unit/test_configuration.py tests/unit/test_cli_gui_parity.py
pixi run lint
pixi run typecheck
```

**Step 3. Inspect the effective parameter provenance on a test workspace**

Use the isolated fixture workspace from D16; inspect its `runs/<RUN_ID>/frozen_config.json`. Confirm each numerical parameter retains a value and source class. Do not relabel `UNCALIBRATED_PROVISIONAL` as calibrated without independent experimental evidence.

### 4. Expected outputs

Updated schema/template/test source, and test records outside pristine source. No study configuration or empirical result is a release default.

### 5. Verification and quality control

Unknown keys and invalid assay combinations must be rejected by existing/new regression tests. Retain both analysis and benchmarking template validation.

### 6. Error handling and troubleshooting

Validation accepting an invalid combination: add a failing regression before correcting the owning model. Template/schema disagreement: reject the change. Do not permit user configuration to disable fixed gates.

### 7. Completion criteria

Configuration tests, lint and type checks pass; documentation and parameter source classes reflect the actual change.


<a name="procedure-d04"></a>

## D04 — Maintain validators and select entry routes

### 1. Procedure identification

**Module:** inputs/registry.py, validate.py, alignment_entry.py. **Workflow:** development.

**Purpose and scope:** Add or correct validation using malformed/valid fixtures, without pretending VCF entry is implemented.

### 2. Prerequisites and requirements

D02–D03; FASTQ pair, BAM/CRAM and reference validation code; tests/unit/test_inputs.py. Test input is generated in pytest temporary directories.

### 3. Detailed operational procedure

**Step 1. Inspect existing validation and routing**

```bash
cd "$DEV"
sed -n '1,150p' src/pfvarx/inputs/registry.py
sed -n '250,360p' src/pfvarx/inputs/validate.py
cat src/pfvarx/inputs/alignment_entry.py
```

**Step 2. Create a focused regression in the existing validator test module**

```bash
nano tests/unit/test_inputs.py
pixi run pytest -q tests/unit/test_inputs.py
nano src/pfvarx/inputs/validate.py
pixi run pytest -q tests/unit/test_inputs.py tests/unit/test_read_evidence.py
```

### 4. Expected outputs

Validated source change and a regression covering real input acceptance/rejection. FASTQ pair names/counts, reference contigs, sort/index/read-group checks remain enforced.

### 5. Verification and quality control

Expect malformed records, wrong REF, missing indexes and wrong SM to fail. VCF entry is explicitly rejected at configuration validation; this manual covers only implemented input routes.

### 6. Error handling and troubleshooting

Never add synthetic @PG labels to claim preprocessing occurred. Do not make unsupported VCF downstream routing pass by dropping its required evidence.

### 7. Completion criteria

All affected validator tests pass and the supported route matrix is updated if a real new entry implementation is added.


<a name="procedure-d05"></a>

## D05 — Develop FASTQ preprocessing and QC

### 1. Procedure identification

**Module:** workflow/rules/l1_fastq_qc.smk, the fastp options in `fastp_args` (workflow/rules/common.smk) and layer1.yaml. **Workflow:** development.

**Purpose and scope:** Maintain fastp behavior without erasing AT-rich sequence or inline UMIs.

### 2. Prerequisites and requirements

D02; original paired reads represented by synthetic branch fixtures; existing fastp rule and layer1 parameter registry.

### 3. Detailed operational procedure

**Step 1. Inspect branch-specific options**

```bash
cd "$DEV"
cat workflow/rules/l1_fastq_qc.smk workflow/rules/common.smk
cat config/parameters/layer1.yaml
```

**Step 2. Edit only the owning rule/parameter and exercise assay branches**

```bash
nano workflow/rules/l1_fastq_qc.smk
pixi run verify wgs targeted umi --output "$LAB/verification-D05"
```

**Step 3. Review generated fastp reports**

In each new scenario workspace inspect `results/qc/samples/<S>/<S>.<LANE>.fastp.json` and `.html`. UMI extraction must precede any trimming that would remove inline molecular tags. Preserve the AT-rich low-complexity policy; poly-G handling follows the declared chemistry.

### 4. Expected outputs

Updated rule and per-lane fastp JSON/HTML in disposable verification workspaces. Trimmed FASTQ lives under run `results/work/trimmed/` and may be deleted by Snakemake `temp()` semantics.

### 5. Verification and quality control

Compare paired retained-read counts, length distributions and adapter/quality summaries against the predeclared fixture expectations; do not require a particular retained percentage for every assay.

### 6. Error handling and troubleshooting

UMI structure failure after trimming: inspect extraction order. Unexpected AT-rich depletion: review any newly enabled complexity filter. Revert unqualified changes rather than relaxing QC.

### 7. Completion criteria

All three branches run correctly and preprocessing differences are documented with output comparisons.


<a name="procedure-d06"></a>

## D06 — Develop reference and annotation preparation

### 1. Procedure identification

**Module:** references/prepare.py, annotation/snpeff_db.py. **Workflow:** development.

**Purpose and scope:** Maintain content-aware derived indexes and reference-matched annotation.

### 2. Prerequisites and requirements

D02; FASTA/GFF/BED fixture resources; tests and separate verification output directory.

### 3. Detailed operational procedure

**Step 1. Inspect preparation and annotation code**

```bash
cd "$DEV"
sed -n '1,240p' src/pfvarx/references/prepare.py
sed -n '1,220p' src/pfvarx/annotation/snpeff_db.py
```

**Step 2. Implement a regression for a verified preparation defect**

```bash
nano tests/unit/test_inputs.py
pixi run pytest -q tests/unit/test_inputs.py
nano src/pfvarx/references/prepare.py
pixi run verify wgs --output "$LAB/verification-D06"
```

**Step 3. Check the preparation contract**

The runner checks that `prepare --dry-run` reports no actions after unchanged preparation. Verify outputs under `prepared/`; original FASTA/GFF/BED hashes must remain unchanged. Inspect annotation QC and cohort ANN fields for FO004.

### 4. Expected outputs

Derived parasite/combined FASTA, indexes and preparation manifest in the test workspace; SnpEff-derived annotation only after matching-reference checks.

### 5. Verification and quality control

Contig IDs/lengths and GFF coordinates must agree; wrong resources must fail. Changed resource checksums must invalidate stale derivatives.

### 6. Error handling and troubleshooting

Never repair a reference mismatch by renaming contigs or substituting another release without documented mapping/validation. Do not index in the original protected data folder.

### 7. Completion criteria

Idempotency, stale detection and annotation-reference checks pass; record source resource identities.


<a name="procedure-d07"></a>

## D07 — Maintain alignment and BAM/CRAM processing

### 1. Procedure identification

**Module:** l1_align.smk, l1_umi.smk, alignment_entry.py. **Workflow:** development.

**Purpose and scope:** Preserve read groups, WGS duplicate marking, targeted primer clipping and UMI consensus branches.

### 2. Prerequisites and requirements

D02/D05/D06; lane/library metadata, matching targets/primers where required; existing alignment validation contracts.

### 3. Detailed operational procedure

**Step 1. Inspect existing rules and entry contract**

```bash
cd "$DEV"
cat workflow/rules/l1_align.smk
cat workflow/rules/l1_umi.smk
cat src/pfvarx/inputs/alignment_entry.py
```

**Step 2. Edit the relevant branch and verify it end to end**

```bash
nano workflow/rules/l1_align.smk
pixi run verify fo001 targeted umi --output "$LAB/verification-D07"
```

**Step 3. Inspect an actual generated alignment**

```bash
cd "$LAB/verification-D07/fo001"
SAMPLE=PF001
pixi run samtools quickcheck -v "results/alignments/$SAMPLE.bam"
pixi run samtools view -H "results/alignments/$SAMPLE.bam"
pixi run samtools flagstat "results/alignments/$SAMPLE.bam"
```

### 4. Expected outputs

Coordinate-sorted BAM/BAI with correct sample/library/read-group provenance and branch QC. Existing CRAM is converted to a processed BAM; CRAM is not the standard final output.

### 5. Verification and quality control

Require correct @SQ, @RG SM and real @PG processing records. WGS duplicate flags are not equivalent to removing all duplicate reads. Preserve target primer semantics.

### 6. Error handling and troubleshooting

Missing or mixed read-group SM: correct upstream metadata/processing in a new derived file, not by faking a label. CRAM reference failures require its exact encoding reference.

### 7. Completion criteria

Fixture branches and BAM integrity checks pass; qualify real CRAM/targeted-alignment entry separately before routine use.


<a name="procedure-d08"></a>

## D08 — Develop the custom SNP and indel caller

### 1. Procedure identification

**Module:** caller/call.py, stats.py, pairhmm.py. **Workflow:** development.

**Purpose and scope:** Maintain the custom implementation using independent mathematical and read-evidence tests.

### 2. Prerequisites and requirements

D02/D06/D07; numerical tests, candidate/evidence schema and frozen-configuration contract. No external caller substitution is permitted.

### 3. Detailed operational procedure

**Step 1. Run the current scientific checks before editing**

```bash
cd "$DEV"
pixi run pytest -q tests/scientific/test_caller_math.py tests/unit/test_read_evidence.py
pixi run python -m pfvarx.caller.call --help
```

**Step 2. Implement a regression then a minimal fix**

```bash
nano tests/scientific/test_caller_math.py tests/unit/test_read_evidence.py
pixi run pytest -q tests/scientific/test_caller_math.py tests/unit/test_read_evidence.py
nano src/pfvarx/caller/call.py src/pfvarx/caller/pairhmm.py
pixi run pytest -q tests/scientific tests/unit/test_read_evidence.py
pixi run verify wgs fo003 --output "$LAB/verification-D08"
```

**Step 3. Review numerical and output contracts**

Use SciPy beta-binomial calculations as independent reference as in existing tests. Check pair-HMM preference on genuine REF/ALT read fixtures, strand/read orientation, overlapping-pair handling, indel homopolymers and finite/missing WSAF intervals. Production entry remains `pixi run pipeline`, not ad hoc caller invocation that bypasses gates.

### 4. Expected outputs

Custom caller source/tests; sample VCF, evidence TSV, callability/focal-loss BED and caller QC JSON from isolated fixtures.

### 5. Verification and quality control

Every altered mathematical branch must have an independent expected result. QUAL/SCR are screening-derived scores, not an automatically calibrated posterior probability. Preserve GT/AD/DP/PST/WSAF/interval/LHID semantics.

### 6. Error handling and troubleshooting

Numerical failure: isolate likelihood, boundary or normalization issue before changing tolerances. Do not change thresholds to hide inaccurate test results.

### 7. Completion criteria

Mathematical/read evidence tests and end-to-end products pass; unresolved scientific-model changes require a separate reviewed protocol.


<a name="procedure-d09"></a>

## D09 — Maintain filtering, status classes and callability

### 1. Procedure identification

**Module:** quality/, caller/, l2_caller.smk. **Workflow:** development.

**Purpose and scope:** Preserve scientific missingness and status-specific filtering.

### 2. Prerequisites and requirements

D08; QC policy and protected parameter definitions; candidate and coverage fixtures.

### 3. Detailed operational procedure

**Step 1. Inspect status and policy ownership**

```bash
cd "$DEV"
cat src/pfvarx/constants.py
cat config/registry/qc_policy_pfvarx_qc_v1.yaml
cat workflow/rules/l2_caller.smk
```

**Step 2. Edit only the confirmed owner and run gate regressions**

```bash
nano tests/unit/test_quality_eligibility.py
pixi run pytest -q tests/unit/test_quality_eligibility.py
pixi run verify benchmark --output "$LAB/verification-D09"
```

**Step 3. Review filter counts on the fixture**

```bash
cd "$LAB/verification-D09/benchmark"
SAMPLE=PF001
pixi run bcftools query -f '%FILTER\n' "results/variants/samples/$SAMPLE.vcf.gz" | sort | uniq -c
cat "results/qc/samples/$SAMPLE/qc_decision.json"
```

### 4. Expected outputs

Preserved PASS/LOW_CONFIDENCE and other declared statuses, callable intervals and centralized QC evidence. Non-PASS records remain available for review.

### 5. Verification and quality control

Uncallable positions must not be recoded as reference. NOT_EVALUATED and NOT_APPLICABLE are distinct from failures. Examine low-depth exclusions without forcing inclusion.

### 6. Error handling and troubleshooting

Inconsistent FILTER versus FORMAT/PST: trace the owning normalization/status code and regression. Do not manufacture all-PASS VCFs.

### 7. Completion criteria

Gate/filter tests pass and all status changes, if intentionally introduced, have a reviewed compatibility rationale.


<a name="procedure-d10"></a>

## D10 — Develop diversity, local haplotypes and clonality

### 1. Procedure identification

**Module:** downstream/diversity.py, coi.py, caller/mixture.py. **Workflow:** development.

**Purpose and scope:** Maintain distinct within-host diversity, local read-backed linkage, technical mixture K and MOIRE COI outputs.

### 2. Prerequisites and requirements

D08–D09; allele evidence, callable denominators, MOIRE installed; shared marker sufficiency. Full clone-resolved haplotypes are not implemented.

### 3. Detailed operational procedure

**Step 1. Inspect Python and R owners**

```bash
cd "$DEV"
cat workflow/scripts/da001_moire.R
sed -n '1,180p' src/pfvarx/downstream/coi.py
pixi run pytest -q tests/scientific/test_caller_math.py
```

**Step 2. Edit an identified component and exercise the benchmark cohort**

```bash
nano src/pfvarx/downstream/coi.py src/pfvarx/downstream/diversity.py
pixi run verify benchmark --output "$LAB/verification-D10"
```

**Step 3. Inspect actual COI/diversity evidence**

```bash
cd "$LAB/verification-D10/benchmark"
head results/downstream/DA001/clonality.tsv
cat results/downstream/DA001/status_manifest.json
head results/downstream/DA003/within_host_fws.tsv
head results/downstream/DA003/population_diversity.tsv
```

### 4. Expected outputs

Clonality TSV, MOIRE summaries, FWS and population diversity. VCF LHID/LHWIN identifies local read-backed linkage scope, not reconstructed complete clones.

### 5. Verification and quality control

Review informative loci, posterior bounds/probabilities, convergence evidence actually supplied and sample exclusions. No diagnostic absent from MOIRE output is assumed successful.

### 6. Error handling and troubleshooting

Few loci/unstable estimates: retain INDETERMINATE or skip and report uncertainty. Technical K agreement alone does not validate COI.

### 7. Completion criteria

Distinct estimands and uncertainty are preserved; proceed only when the implementation and docs do not overstate identifiability.


<a name="procedure-d11"></a>

## D11 — Maintain resistance marker analysis

### 1. Procedure identification

**Module:** downstream/resistance.py and catalogue registry. **Workflow:** development.

**Purpose and scope:** Maintain codon-based marker calls with callable denominators and independently verified catalogue provenance.

### 2. Prerequisites and requirements

D06/D09; reference-matched FASTA/GFF, verified marker catalogue or explicitly TEST_FIXTURE for software tests only.

### 3. Detailed operational procedure

**Step 1. Inspect the guard and genotype logic**

```bash
cd "$DEV"
head -12 config/registry/resistance_catalog_pfvarx_v1.tsv
sed -n '1,220p' src/pfvarx/downstream/resistance.py
pixi run pytest -q tests/unit/test_benchmark_and_resistance.py
```

**Step 2. Run the existing labelled marker fixture test**

```bash
pixi run verify benchmark --output "$LAB/verification-D11"
cd "$LAB/verification-D11/benchmark"
head results/downstream/DA002/marker_calls.tsv
head results/downstream/DA002/prevalence.tsv
cat results/downstream/DA002/catalog_provenance.json
```

**Step 3. Review any real catalogue update**

Verify each gene/codon/reference residue/alternative and evidence source against the chosen assembly and an authoritative marker publication. Store independent review evidence and catalogue checksum. Only after that review may a maintainer supply a genuinely verified catalogue; changing a status string alone is not verification.

### 4. Expected outputs

Marker calls, prevalence, Wilson uncertainty and catalogue provenance when eligible. Draft catalogue execution is blocked as designed.

### 5. Verification and quality control

UNKNOWN never contributes a wild-type observation; verify numerator and callable denominator. Marker presence is a molecular finding, not a guaranteed patient treatment outcome.

### 6. Error handling and troubleshooting

Draft/mismatched catalogue: stop DA002, obtain proper independent verification. Never relabel the draft or a TEST_FIXTURE catalogue for real surveillance.

### 7. Completion criteria

Guard, reference/codon checks and denominator tests pass; real catalogue promotion includes evidence and review sign-off.


<a name="procedure-d12"></a>

## D12 — Maintain population structure, relatedness, trees and spatial analysis

### 1. Procedure identification

**Module:** downstream/{relatedness,structure,phylogeny,spatial}.py. **Workflow:** development.

**Purpose and scope:** Integrate actual haploid PCA/FST, eligibility-specific relatedness, monoclonal phylogeny and supported spatial outcomes.

### 2. Prerequisites and requirements

D10–D11; population metadata for DA005, monoclonal evidence for DA006, actual country/date/coordinates and source outcome for DA007.

### 3. Detailed operational procedure

**Step 1. Inspect method ownership and prerequisites**

```bash
cd "$DEV"
cat config/registry/downstream_analyses.yaml
cat config/parameters/layer5.yaml
cat workflow/scripts/da004_dcifer.R workflow/scripts/da007_gam.R
```

**Step 2. Exercise all implemented methods on the informative fixture**

```bash
pixi run verify benchmark --output "$LAB/verification-D12"
cd "$LAB/verification-D12/benchmark"
find results/downstream -name status_manifest.json -print
cat results/downstream/DA006/runtime_gates.json
head results/downstream/DA005/pca.tsv
head results/downstream/DA005/fst.tsv
```

### 4. Expected outputs

Pairwise relatedness TSV, PCA/FST summaries, tree only when gates permit, and country-month/grid summaries; conditional GAM output where supported.

### 5. Verification and quality control

Check method-specific eligibility, haploid representation and sample counts. DA006 requires sufficient monoclonal samples/informative sites; pairwise methods cannot be mixed into one common estimand. Review spatial sampling bias.

### 6. Error handling and troubleshooting

Unsupported targeted PCA/tree request: reject it. Missing metadata: correct verified metadata before a NEW run. No informative tree: report skip, never a placeholder tree.

### 7. Completion criteria

All requested methods have truthful terminal manifests and reproducible outputs or justified skip/block states.


<a name="procedure-d13"></a>

## D13 — Integrate and execute Snakemake changes

### 1. Procedure identification

**Module:** workflow/Snakefile, workflow/rules/ (rules; helper functions in common.smk), execution/runner.py. **Workflow:** development.

**Purpose and scope:** Modify scheduling without bypassing Layer 0 plans, QC checkpoints or failure isolation.

### 2. Prerequisites and requirements

D02/D03 and affected component SOPs; Linux local executor; review per-job/global memory, disk and thread envelopes.

### 3. Detailed operational procedure

**Step 1. Inspect DAG and rule naming**

```bash
cd "$DEV"
cat workflow/Snakefile
pixi run pytest -q tests/integration
```

**Step 2. Edit the owning rule and verify the real checkpoint path**

```bash
nano workflow/rules/l2_caller.smk
pixi run verify fo001 fo003 --output "$LAB/verification-D13"
```

**Step 3. Verify dry-run isolation and resource declarations**

In a new fixture workspace run `pixi run ready` and `pixi run dry-run`. It must not create `runs/` or `results`. Runtime checkpoints can change later scheduling; a dry run is not evidence that all later scientific methods execute.

### 4. Expected outputs

Updated existing rule or runner with declared inputs/outputs/logs/benchmarks and a tested dependency chain. Runtime completion marker is owned by the real workflow.

### 5. Verification and quality control

Inputs and outputs must use existing run-relative contracts. Global/per-job budgets must fit. Keep local execution; SLURM is unsupported in this release.

### 6. Error handling and troubleshooting

Missing input/cyclic graph: fix dependencies, not by touching fake outputs. Scientific exclusion is not a transient tool error and must not trigger indiscriminate retries.

### 7. Completion criteria

Integration dry runs and affected end-to-end cases pass, with checkpoint exclusion/failure behavior preserved.


<a name="procedure-d14"></a>

## D14 — Maintain diagnostics and reproducibility

### 1. Procedure identification

**Module:** diagnostics/, provenance/, runs/, io.py. **Workflow:** development.

**Purpose and scope:** Preserve frozen runs, content hashes and useful redacted diagnostics.

### 2. Prerequisites and requirements

D13; run integrity contracts, error-handling tests and no secret-bearing environment capture.

### 3. Detailed operational procedure

**Step 1. Inspect persistence and diagnostic owners**

```bash
cd "$DEV"
cat src/pfvarx/provenance/snapshot.py
cat src/pfvarx/diagnostics/events.py
pixi run pytest -q tests/unit/test_reporting.py
```

**Step 2. Exercise frozen integrity on an isolated completed run**

```bash
cd "$LAB/verification-D13/fo003"
pixi run status
# Set RUN_ID to the exact completed identifier printed above.
export RUN_ID=YOUR_COMPLETED_RUN_ID
pixi run run-verify --run-id "$RUN_ID"
cat "runs/$RUN_ID/status.json"
```

### 4. Expected outputs

`snapshot.json`, frozen configuration/plan/code, status and events/tool logs; separate renderer provenance for refreshed reports.

### 5. Verification and quality control

Frozen-file integrity must PASS. Workspace drift is a separate signal. External reads/truth/masks are not all frozen or reuse-key-covered; maintain independent checksums.

### 6. Error handling and troubleshooting

Frozen mismatch: stop reuse, recover correct source or start NEW with documented inputs. Do not edit the snapshot to erase drift. Do not log `env` output or credentials.

### 7. Completion criteria

Failure causes are distinguished from QC/eligibility; reproducibility records are sufficient to identify the exact experiment.


<a name="procedure-d15"></a>

## D15 — Maintain the offline reporting dashboard

### 1. Procedure identification

**Module:** reporting/ assets and build.py. **Workflow:** development.

**Purpose and scope:** Modify the real modular offline dashboard without changing scientific estimators or paths.

### 2. Prerequisites and requirements

D02/D14; existing report schema/template/CSS/core/tables/charts/views/exports and browser QA runner. Browser dependencies require a separate environment.

### 3. Detailed operational procedure

**Step 1. Inspect and edit the existing frontend**

```bash
cd "$DEV"
ls src/pfvarx/reporting
nano src/pfvarx/reporting/views.js src/pfvarx/reporting/dashboard.css
pixi run pytest -q tests/unit/test_reporting.py
```

**Step 2. Prepare browser QA and generate a synthetic report externally**

```bash
pixi run python -m venv "$LAB/browser-env"
"$LAB/browser-env/bin/pip" install -r tests/dashboard/requirements.txt
"$LAB/browser-env/bin/python" -m playwright install chromium
pixi run verify benchmark --output "$LAB/verification-D15"
"$LAB/browser-env/bin/python" tests/dashboard/check_dashboard.py "$LAB/verification-D15/benchmark/results/report" "$LAB/browser-check-D15"
```

**Step 3. Inspect screen and exported artifacts**

Open the generated HTML locally; test all eight sections and every available visualization, exact hover values, missing states, filters, pagination and export feedback. Inspect actual CSV/XLSX contents, PNG labels/white background and printed selected views. No external network request is permitted.

### 4. Expected outputs

One self-contained HTML, dashboard JSON and centralized QC JSON at the existing paths in the test run; browser checks and temporary export artifacts outside source.

### 5. Verification and quality control

Synthetic banner remains explicit. Missing/null never becomes zero. Geographic positions/uncertainty/distances require real source data. CSV/XLSX must honor searches and active filters.

### 6. Error handling and troubleshooting

JavaScript/schema/export failure: inspect browser console and source field mapping. Do not substitute demonstration data for absent scientific results.

### 7. Completion criteria

Backend and browser QA pass; documentation reflects current table/PDF scope; reporting paths and scientific code remain unchanged.


<a name="procedure-d16"></a>

## D16 — Test, package and prepare a release

### 1. Procedure identification

**Module:** tests/, scripts/verify.py, package.py and CI. **Workflow:** development.

**Purpose and scope:** Produce a reviewed software-only release from tested source.

### 2. Prerequisites and requirements

D01–D15 as applicable; temporary external space; private-data exclusion rules; reviewed Git identity and authorized remote for any publication.

### 3. Detailed operational procedure

**Step 1. Run source and workflow checks**

```bash
cd "$DEV"
pixi run lint
pixi run typecheck
pixi run pytest -q tests/unit tests/scientific tests/integration
pixi run verify --output "$LAB/verification-D16"
```

**Step 2. Build and inspect the software-only archive**

```bash
pixi run package --output "$LAB/release"
cd "$LAB/release"
sha256sum -c pfvarx-*-source.zip.sha256
python3 -m zipfile -e pfvarx-*-source.zip extracted
python3 extracted/pfvarx/scripts/package.py --check
```

**Step 3. Prepare publication after review**

```bash
cd "$DEV"
git status --short
git diff
git diff --check
# Review and commit only the approved source change.
git add src workflow config tests docs scripts .github bootstrap.py pfvarx README.md CHANGELOG.md pixi.toml pixi.lock pyproject.toml .gitignore LICENSE
git diff --cached --stat
```

Publication/commit/tag actions require the repository owner’s release decision. The release version is set only in `src/pfvarx/__init__.py`; `tests/unit/test_version.py` checks every current-release reference against it; regenerate/review the lock only when appropriate. Do not retag a published release to conceal changes. Run hosted CI before a new tag.

### 4. Expected outputs

Source ZIP and SHA-256 sidecar containing no generated datasets/reports/environments, with maintained manuals and tests. Full test outputs remain outside source.

### 5. Verification and quality control

All executed checks must have actual exit-status evidence. Optional browser pytest skipping is not a browser pass; run D15 separately. Clean extraction must report zero runtime files.

### 6. Error handling and troubleshooting

Test failure: identify product/test/environment cause and preserve the failure record. Do not delete a failing test or disable a method to obtain green CI.

### 7. Completion criteria

Reviewed source and archive match, tests are recorded, limitations are disclosed and publication is separately authorized. Remove disposable generated test work after evidence capture.


# PART II — Routine sample analysis

Follow R01–R04, then select **one assay/entry branch** (R05–R10). Branches are alternatives,
not a sequence that reprocesses one study under conflicting assay labels. Use a new
workspace for each assay. Configure downstream selections in R11 **before** the NEW run
if they are needed; otherwise FO002 is the first caller pilot. Then perform R12–R15.
Do not launch twice: R11 only configures/preflights downstream selections; the chosen
entry branch launches the run. Afterward return to R11 Step 3 for manifest review.
For BAM/CRAM, register the alignment column and input_type in the aligned-entry SOP;
use the reference-registration portions of R03 and defer validation until those fields
are correct. Do not register fictional FASTQ paths for an aligned-input experiment.
Registered sample/library/lane identifiers are mandatory; country/population/date/location
metadata are optional only when the requested analyses do not require them.


<a name="procedure-r01"></a>

## R01 — Create and install a routine-analysis workspace

### 1. Procedure identification

**Module:** initialization and dependency management. **Workflow:** routine analysis.

**Purpose and scope:** Create a writable workspace independent of pristine source and benchmarking.

### 2. Prerequisites and requirements

Common variables set; Linux x86-64; Pixi installed; matching resources to be supplied. `$ANALYSIS` must not already exist. Choose an assay-specific directory if running a targeted/BAM branch.

### 3. Detailed operational procedure

**Step 1. Create an isolated analysis workspace**

```bash
test ! -e "$ANALYSIS"
python3 "$PFVARX_SOURCE/bootstrap.py" "${BOOTSTRAP_ARGS[@]}" --workspace "$ANALYSIS"
cd "$ANALYSIS"
./pfvarx version
./pfvarx --help
```

**Step 2. Configure the workspace**

Run `./pfvarx configure --gui` and review the assay, reference release, input paths,
caller and compute budget. Use the backup option to update YAML created by bootstrap.
For a headless session, edit the three YAML files directly. See Appendix F for the
GUI-first sequence and remote access.

**Step 3. Record hardware and installation identity**

```bash
mkdir -p "$RECORDS/analysis"
nproc > "$RECORDS/analysis/cpu-count.txt"
free -h > "$RECORDS/analysis/memory.txt"
df -h "$ANALYSIS" /tmp > "$RECORDS/analysis/disk.txt"
sha256sum pixi.lock > "$RECORDS/analysis/pixi-lock.sha256"
cat config/compute.yaml
```

### 4. Expected outputs

Bootstrap installs `.pixi` and writes three editable configuration files. Save the GUI configuration or run `./pfvarx init` after YAML edits to create input folders.

### 5. Verification and quality control

CLI version must identify the chosen release; `config-init` must not overwrite existing config. 4 cores/6500 MiB workflow RAM/20 GB temporary disk in the template are planner examples, not a whole-genome sizing guarantee. Set available budgets honestly.

### 6. Error handling and troubleshooting

Existing workspace: stop and inspect it, or choose another new directory. Insufficient quota/permission: relocate the new workspace or obtain resources. Do not reduce scientific thresholds to fit hardware.

### 7. Completion criteria

Record source/lock identity and hardware; installation complete; proceed to R02.


<a name="procedure-r02"></a>

## R02 — Register one sample or a cohort

### 1. Procedure identification

**Module:** prepared/samples.tsv and config/config.yaml. **Workflow:** routine analysis.

**Purpose and scope:** Create a tab-separated manifest with explicit read pairing and library provenance.

### 2. Prerequisites and requirements

R01; actual paired FASTQ files and laboratory sample/lane/library mapping. gzip FASTQ must be intact; population/location values must come from verified metadata.

### 3. Detailed operational procedure

**Step 1. Edit an explicit registration list in the workspace**

```bash
cd "$ANALYSIS"
nano prepared/samples.tsv
```

Use the generated columns, real tabs and one row per sample/lane. Do not infer identity from a filename alone. For a one-sample pilot write one data row; for a cohort add all verified rows. The following Python example creates a real tab-separated file; **replace the sample IDs, paths and metadata before executing**:

```bash
./pfvarx tool python - <<'PYCODE'
import csv
from pathlib import Path
columns=['sample_id','lane','fastq_1','fastq_2','bam','vcf','library_id','flowcell_type','population','country','latitude','longitude','collection_date']
rows=[dict(sample_id='S01',lane='L001',fastq_1='/absolute/path/S01_R1.fastq.gz',fastq_2='/absolute/path/S01_R2.fastq.gz',library_id='LIB01'),
      dict(sample_id='S02',lane='L001',fastq_1='/absolute/path/S02_R1.fastq.gz',fastq_2='/absolute/path/S02_R2.fastq.gz',library_id='LIB02')]
# Delete the second row for a one-sample pilot. Add verified metadata as needed.
with Path('prepared/samples.tsv').open('w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=columns,delimiter='	');w.writeheader();w.writerows(rows)
PYCODE
```

**Step 2. Review the exact registered files**

```bash
column -t -s $'\t' prepared/samples.tsv
./pfvarx tool python - <<'PYCODE'
import csv
from pathlib import Path
for r in csv.DictReader(open("prepared/samples.tsv"), delimiter="\t"):
    for k in ("fastq_1","fastq_2"):
        p=Path(r[k]); print(r["sample_id"],r["lane"],k,p,p.is_file())
PYCODE
```

If `column` is unavailable, use `cat` instead. All path existence checks must be true before validation. IDs use schema-permitted letters/digits/underscore/hyphen; sample/lane pairs must be unique. Multiple lanes share a sample ID and retain accurate library IDs. For WGS, all lanes of one sample must share one `flowcell_type` (duplicates are marked once per sample with one optical distance); validation rejects a sample that mixes patterned and unpatterned lanes.

### 4. Expected outputs

`prepared/samples.tsv`, containing required sample/lane/library fields and branch-specific file columns. Blank unused fields are preserved.

### 5. Verification and quality control

Confirm R1/R2 mapping with the wet-lab manifest. Do not use the analysis output as a source of truth for sample identity or population labels.

### 6. Error handling and troubleshooting

Missing/duplicate row or wrong pair: correct the manifest before running. Missing required population/date/coordinate fields: supply verified values or do not request the dependent analysis.

### 7. Completion criteria

Operator signs off sample identities, paths, library design, optional metadata and expected cohort size.


<a name="procedure-r03"></a>

## R03 — Register and prepare the exact reference resources

### 1. Procedure identification

**Module:** reference configuration and references/. **Workflow:** routine analysis.

**Purpose and scope:** Supply matching FASTA, core mask and annotation, derive indexes without modifying original resources.

### 2. Prerequisites and requirements

R01–R02; independently justified reference FASTA and core-genome mask BED (the accessible core: Core regions of Miles et al. 2016, as used by MalariaGEN; not a gene or feature BED); matching GFF for explicit FO004 or automatic FO004 marker outcomes. BED is zero-based half-open; GFF/VCF positions are one-based.

### 3. Detailed operational procedure

**Step 1. Set the actual reference and initial analysis mode**

```bash
cd "$ANALYSIS"
nano config/config.yaml
```

Edit the generated file without deleting schema/other keys. A valid edit fragment is:

```yaml
project_id: STUDY_A
final_output: FO002
input_type: FASTQ
reference:
  fasta: /absolute/path/to/verified/Pf3D7.fasta
  gff: null
  core_mask_bed: /absolute/path/to/verified/core_genome.bed
  release_label: YOUR_EXACT_REFERENCE_RELEASE
  snpeff_genome_id: null
```

Use actual filenames and release text. For annotation, set `gff` to that release's matching resource. Never substitute the earlier simulated pilot mask as a validated study core definition.

**Step 2. Validate and derive reference products**

```bash
./pfvarx check-layout
./pfvarx validate
./pfvarx prepare --dry-run
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx prepare --dry-run
find prepared -maxdepth 3 -type f | sort
```

**Step 3. Record original resource checksums**

```bash
sha256sum /absolute/path/to/verified/Pf3D7.fasta /absolute/path/to/verified/core_genome.bed > "$RECORDS/analysis/reference.sha256"
```

Add the actual GFF, optional host, targets/primers and verified catalogue to this manifest when used. For long filenames retain proper quoting. Record download/release provenance and licensing separately.

### 4. Expected outputs

`prepared/reference/parasite.fasta`, combined reference and their applicable indexes; `prepared/manifest.json`. Required source files stay read-only.

### 5. Verification and quality control

The second dry preparation plan should have no actions for unchanged resources. All contig names, lengths and intervals must match; inspect `metadata/validation_report.json` for source facts and failures. A `CORE_MASK_SUSPECT` warning means the core BED looks like a gene list (gene names, more than 1 000 or overlapping intervals, or coverage outside 70–97% of the reference); replace it with the real core mask before preparation.

### 6. Error handling and troubleshooting

Mismatched FASTA/GFF/BED: stop. Obtain consistent resources; do not silently rename contigs or regenerate a whole-genome mask to remove errors. Automatic FO004 can leave GFF unset. Layer 4 and its dependent marker/spatial outputs then remain unavailable.

### 7. Completion criteria

Input validation accepted, reference preparation current, checksums recorded and the biological core mask justified independently.


<a name="procedure-r04"></a>

## R04 — Preflight resources, validation and the workflow

### 1. Procedure identification

**Module:** CLI ready/dry-run and config/compute.yaml. **Workflow:** routine analysis.

**Purpose and scope:** Make an explicit launch decision before creating a run.

### 2. Prerequisites and requirements

R01–R03; chosen assay configuration and sample sheet; any R11 selections already applied. If branch/config/input changes afterward, repeat this SOP.

### 3. Detailed operational procedure

**Step 1. Set a realistic compute envelope**

```bash
cd "$ANALYSIS"
nano config/compute.yaml
cat config/compute.yaml
```

Use available total cores/RAM/jobs/disk and per-job ceilings. The current template is 4 cores, 6500 MiB total/per-job RAM, one concurrent job and 20000 MB temporary disk. These are illustrative planning limits. Leave RAM/disk headroom and use measured profiles from comparable studies; requests do not impose universal OS limits.

**Step 2. Run each preflight command and inspect its JSON**

```bash
set -o pipefail
./pfvarx check-layout
./pfvarx validate > "$RECORDS/analysis/validation-command.json"
cat "$RECORDS/analysis/validation-command.json"
./pfvarx prepare --dry-run
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready > "$RECORDS/analysis/ready.json"
cat "$RECORDS/analysis/ready.json"
./pfvarx dry-run > "$RECORDS/analysis/dry-run.json"
cat "$RECORDS/analysis/dry-run.json"
```

### 4. Expected outputs

Validation report, current preparation, `ready: true` and `dry-run` status `OK`. On a new workspace the dry run does not create `runs/` or `results`.

### 5. Verification and quality control

Review every failed sample, warning and deferred gate. The validation `summary` states for each problem or warning the input, why it failed, why the input matters and numbered steps to fix it. PARTIAL permits only individually valid samples under continue-independent policy; record exclusions explicitly. Reference/global validation failure must be resolved.

### 6. Error handling and troubleshooting

False readiness: inspect JSON resources and prep status. Parameter UNSET_BLOCKING/unsupported metadata: correct the owning requirement, not a generated plan. A dry-run pass does not guarantee later QC/clonality eligibility.

### 7. Completion criteria

Launch only with accepted validation, `ready: true`, dry-run OK and signed-off exclusions/resources.


<a name="procedure-r05"></a>

## R05 — Execute paired-end WGS FASTQ analysis

### 1. Procedure identification

**Module:** WGS FASTQ branch, layers 1–2 or selected later mode. **Workflow:** routine analysis.

**Purpose and scope:** Process direct or genuine SWGA gDNA reads using the implemented WGS branch.

### 2. Prerequisites and requirements

R01–R04 complete; paired FASTQ and matching resources. Choose FO002 for the first caller pilot or configure FO003/FO004 and R11 beforehand. Do not choose SWGA for unamplified gDNA.

### 3. Detailed operational procedure

**Step 1. Confirm the assay settings before repeating preflight if needed**

```bash
cd "$ANALYSIS"
nano config/config.yaml
```

```yaml
assay_design: wgs
assay:
  library_catalog: null
  two_colour_chemistry: false
  platform: ILLUMINA
```

These are example settings for direct non-UMI gDNA; change chemistry/amplification only from real library metadata. Do not replace the full configuration with this fragment.

**Step 2. Execute one NEW run and capture the result**

```bash
./pfvarx pipeline --mode NEW > "$RECORDS/analysis/pipeline.json"
cat "$RECORDS/analysis/pipeline.json"
./pfvarx status
# Copy the exact returned run identifier, not this placeholder.
export RUN_ID=YOUR_COMPLETED_RUN_ID
./pfvarx run-verify --run-id "$RUN_ID"
```

**Step 3. Inspect WGS output integrity**

```bash
SAMPLE=S01
./pfvarx tool samtools quickcheck -v "results/alignments/$SAMPLE.bam"
./pfvarx tool samtools flagstat "results/alignments/$SAMPLE.bam"
cat "results/qc/samples/$SAMPLE/qc_decision.json"
cat "results/qc/samples/$SAMPLE/$SAMPLE.caller_qc.json"
```

### 4. Expected outputs

Per-sample BAM/BAI, fastp/duplicate/coverage QC and FO002+ VCF/evidence/callability outputs; run logs/benchmarks and the offline dashboard. One sample produces one processed sample set; a multi-row cohort produces one set per admitted sample.

### 5. Verification and quality control

Run status and frozen integrity must pass; inspect mapping, duplicate burden, median depth, breadth/callability and warnings with their source thresholds. No universal accepted coverage value is invented by this manual.

### 6. Error handling and troubleshooting

Caller is CPU-active but silent: inspect process/resource use and log timestamps; do not launch a duplicate job. Failures/exclusions: follow R12–R13. A completed process can still contain QC-blocked samples.

### 7. Completion criteria

Record exact run ID and all per-sample acceptance/exclusion decisions before interpretation or downstream promotion.


<a name="procedure-r06"></a>

## R06 — Execute targeted non-UMI FASTQ analysis

### 1. Procedure identification

**Module:** targeted_no_umi branch (assay_design). **Workflow:** routine analysis.

**Purpose and scope:** Apply amplicon alignment/primer clipping without imposing WGS duplicate assumptions.

### 2. Prerequisites and requirements

Use a new `$LAB/analysis-targeted-non-umi` workspace: set `ANALYSIS` to it and complete R01–R03 there. Actual target BED and fgbio primer TSV required; no UMIs in the library.

### 3. Detailed operational procedure

**Step 1. Declare the real targeted assay**

```bash
export ANALYSIS="$LAB/analysis-targeted-non-umi"
BOOTSTRAP_ARGS=(--final-output FO002 --design targeted_no_umi --input-type FASTQ)
# Complete R01–R03 in this new workspace before the next command.
cd "$ANALYSIS"
nano config/config.yaml
```

```yaml
assay_design: targeted_no_umi
assay:
  library_catalog: inputs/metadata/library_catalog.yaml
  two_colour_chemistry: false
  platform: ILLUMINA
```

Run `./pfvarx init` after setting the assay design; it places the library catalog template at `inputs/metadata/library_catalog.yaml`. Define the panel there (`panel_id`, `panel_version`, `targets_bed`, `primers_tsv`), supply the actual panel files, and give every manifest row its `panel_id` and `panel_version` (optional `technique`: `AMPLICON` or `TBS_NON_UMI`). Targets are zero-based half-open. Primer header must be `chrom`, `left_start`, `left_end`, `right_start`, `right_end` with one-based inclusive coordinates and `left_start <= left_end < right_start <= right_end`. Example primer rows are not provided as substitute biology.

**Step 2. Validate the panel and execute after R04**

```bash
./pfvarx init
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline --mode NEW
./pfvarx status
```

**Step 3. Review the primer-clipping path**

```bash
SAMPLE=S01
./pfvarx tool samtools view -H "results/alignments/$SAMPLE.bam"
cat "results/qc/samples/$SAMPLE/qc_decision.json"
find runs -path "*/logs/trim_primers/*" -type f -print
```

### 4. Expected outputs

Primer-clipped BAM/BAI, preclip/coverage QC, sample variants and evidence for callable assayed regions, and report. Panel outputs are not whole-genome observations.

### 5. Verification and quality control

Require primer coordinate validity and real TrimPrimers execution. Review target-specific callability/uneven depth and amplicon bias. DA005 and DA006 run on a panel only when it has at least 10 separate target regions; their results then describe the panel loci, not the genome.

### 6. Error handling and troubleshooting

Panel/reference mismatch or absent primers: stop. Do not reuse a WGS workspace merely by changing its assay after a run. Do not deduplicate amplicon reads using WGS assumptions.

### 7. Completion criteria

Panel and preprocessing provenance verified; accepted targets/codons and inaccessible regions are explicitly documented.


<a name="procedure-r07"></a>

## R07 — Execute targeted UMI FASTQ analysis

### 1. Procedure identification

**Module:** targeted_umi branch (assay_design). **Workflow:** routine analysis.

**Purpose and scope:** Preserve inline tags and run genuine molecular consensus processing before calling.

### 2. Prerequisites and requirements

New `$LAB/analysis-targeted-umi` workspace with R01–R03; actual read structure, target BED and primer TSV; original FASTQ retains inline UMIs.

### 3. Detailed operational procedure

**Step 1. Declare UMI structure from the library design**

```bash
export ANALYSIS="$LAB/analysis-targeted-umi"
BOOTSTRAP_ARGS=(--final-output FO002 --design targeted_umi --input-type FASTQ)
# Complete R01–R03 before continuing.
cd "$ANALYSIS"
nano config/config.yaml
```

```yaml
assay_design: targeted_umi
assay:
  library_catalog: inputs/metadata/library_catalog.yaml
  two_colour_chemistry: false
  platform: ILLUMINA
```

Define the panel and a library protocol in `inputs/metadata/library_catalog.yaml`; the protocol's `umi.read_structure` (for example `8M+T`) and any consensus settings that differ from `layer1.yaml` are used for every manifest row that names it in `library_protocol_id`. `8M+T` is illustrative only for exactly eight inline molecular bases at each read start followed by template. Use actual fgbio structures; never infer them from read filenames or copy this setting for a different library. Allele counts and depths of this design count consensus molecules (`count_unit` in QC and the VCF header).

**Step 2. Validate and execute the UMI branch**

```bash
./pfvarx init
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline --mode NEW
./pfvarx status
```

**Step 3. Inspect consensus and clipping evidence**

```bash
SAMPLE=S01
./pfvarx tool samtools view -H "results/alignments/$SAMPLE.bam"
find runs -path "*/logs/*" -type f | sort
ls "results/qc/samples/$SAMPLE"
cat "results/qc/samples/$SAMPLE/qc_decision.json"
```

### 4. Expected outputs

Consensus-derived realigned primer-clipped BAM/BAI, UMI/QC records and caller products. Raw read depth is not molecular depth; review actual consensus/family metrics supplied by the tools.

### 5. Verification and quality control

Check ExtractUmis/GroupReads/consensus/filter/realignment/clipping logs and read-group/tag integrity. Low-depth small synthetic fixtures establish branch behavior only, not an adequate minimum molecular coverage.

### 6. Error handling and troubleshooting

UMI layout mismatch, removed tags or incorrect first read lengths: stop and obtain original data/design. Do not fabricate RX tags or rerun as non-UMI to conceal failed consensus.

### 7. Completion criteria

Correct library structure and end-to-end molecular processing verified; molecular evidence and callable scope accepted.


<a name="procedure-r08"></a>

## R08 — Enter an existing WGS BAM

### 1. Procedure identification

**Module:** alignment_entry.py and WGS caller branch. **Workflow:** routine analysis.

**Purpose and scope:** Use an already aligned sample while enforcing alignment identity and processing contracts.

### 2. Prerequisites and requirements

New `$LAB/analysis-bam` workspace with R01 and R03 resources; coordinate-sorted BAM/BAI and correct @SQ/@RG. Register one row per actual sample, not per BAM containing unrelated samples.

### 3. Detailed operational procedure

**Step 1. Set the workspace input route and BAM sample row**

```bash
export ANALYSIS="$LAB/analysis-bam"
BOOTSTRAP_ARGS=(--final-output FO002 --design WGS --input-type BAM)
# Complete R01 and reference setup in this new workspace.
cd "$ANALYSIS"
nano config/config.yaml
nano prepared/samples.tsv
```

Set `input_type: BAM`, `assay_design: wgs`, and FO002 or later. In the TSV fill `bam`, `sample_id`, `lane`, `library_id`, the actual `amplification_route` (`DIRECT` or `SWGA`) and `processing_state` (`aligned` to mark duplicates, `duplicates_marked` to use the file as is; empty marks duplicates unless the header records it); leave FASTQ/VCF fields empty. Every read-group SM must equal the registered sample ID. All read groups need ID/SM/LB/PL.

**Step 2. Inspect the supplied BAM read-only**

```bash
export INPUT_BAM=/absolute/path/to/S01.bam
./pfvarx tool samtools quickcheck -v "$INPUT_BAM"
./pfvarx tool samtools view -H "$INPUT_BAM"
./pfvarx tool samtools idxstats "$INPUT_BAM"
```

**Step 3. Resolve index/sort problems in a new derived file, then validate**

If the producer can regenerate a valid file, prefer that. If sorting/indexing is the only problem, preserve the original and create a registered derivative:

```bash
mkdir -p prepared/entry
./pfvarx tool samtools sort -o prepared/entry/S01.sorted.bam "$INPUT_BAM"
./pfvarx tool samtools index prepared/entry/S01.sorted.bam
# Register prepared/entry/S01.sorted.bam in the sample sheet.
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline --mode NEW
```

### 4. Expected outputs

Entry contract JSON, processed BAM/BAI and downstream caller products. WGS entry runs duplicate marking if no corresponding actual program record exists.

### 5. Verification and quality control

Accept only exact matching reference contigs/lengths/checksums where present, coordinate sort, readable index and correct sample read groups. Review entry actions, not just a label.

### 6. Error handling and troubleshooting

Wrong assembly or sample SM: stop; obtain correct input or a scientifically justified new derived file with documented provenance. Do not add fake @PG records to bypass checks.

### 7. Completion criteria

Alignment validation and entry contract accepted; duplicate/coverage/caller QC reviewed and run integrity recorded.


<a name="procedure-r09"></a>

## R09 — Qualify existing targeted alignments

### 1. Procedure identification

**Module:** targeted alignment entry. **Workflow:** routine analysis with qualification.

**Purpose and scope:** Use primer-clipped targeted BAM only when its true upstream state is established.

### 2. Prerequisites and requirements

R08 alignment checks plus correct target/primer definitions and genuine upstream fgbio TrimPrimers @PG. Targeted UMI entry also requires independently documented consensus/read-structure provenance; a header alone is insufficient proof.

### 3. Detailed operational procedure

**Step 1. Declare targeted BAM configuration and inspect provenance**

```bash
export ANALYSIS="$LAB/analysis-targeted-bam"
BOOTSTRAP_ARGS=(--final-output FO002 --design targeted_no_umi --input-type BAM)
# Complete R01 in this new directory and register real reference/panel paths.
cd "$ANALYSIS"
nano config/config.yaml
./pfvarx tool samtools view -H /absolute/path/to/targeted.S01.bam
```

Set BAM plus the actual targeted `assay_design`, the library catalog panel, and each row's `processing_state`: `aligned` or `primers_trimmed` (`targeted_no_umi`); `aligned` (raw reads, UMI in `RX`), `consensus` (consensus is not generated again; primers are trimmed) or `consensus_primers_trimmed` (`targeted_umi`). Validation rejects a state the reads contradict (fgbio `cD` consensus tags, `RX` UMI tags, a TrimPrimers `@PG` record). Use a dedicated workspace, not the WGS BAM workspace. Inspect upstream producer logs and retained molecular/primer evidence.

**Step 2. Apply the gate before executing**

```bash
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
```

The entry adapter rejects targeted alignments without a real TrimPrimers record. Qualification must also confirm the panel coordinates and prior consensus behavior. Only after approval run `./pfvarx pipeline --mode NEW` and inspect the sample entry contract.

### 4. Expected outputs

Accepted entry contract and sample products only after route-specific qualification. Unsupported processing state produces rejection, not an automatically clipped rescue alignment.

### 5. Verification and quality control

Compare a representative sample against the original FASTQ targeted route, including primer boundary/callability and molecular evidence, before approving routine reuse.

### 6. Error handling and troubleshooting

Absent clipping/consensus provenance: stop; use original FASTQ with R06/R07. Never edit @PG to claim a tool ran.

### 7. Completion criteria

Mark the route QUALIFIED only with comparison evidence; otherwise record NOT QUALIFIED and do not interpret its downstream results.


<a name="procedure-r10"></a>

## R10 — Qualify existing CRAM input

### 1. Procedure identification

**Module:** CRAM entry conversion and validation. **Workflow:** routine analysis with qualification.

**Purpose and scope:** Decode a reference-dependent alignment and verify the BAM derivative.

### 2. Prerequisites and requirements

Dedicated CRAM workspace; R08/R09 branch-specific prerequisites; exact CRAM decoding reference and an index. CRAM entry is implemented but less extensively verified.

### 3. Detailed operational procedure

**Step 1. Check decoding with the correct reference**

```bash
export ANALYSIS="$LAB/analysis-cram"
BOOTSTRAP_ARGS=(--final-output FO002 --design WGS --input-type CRAM)
# Complete R01 in this new directory and register real reference paths.
cd "$ANALYSIS"
export INPUT_CRAM=/absolute/path/to/S01.cram
export CRAM_REF=/absolute/path/to/exact-encoding-reference.fasta
./pfvarx tool samtools view -T "$CRAM_REF" -H "$INPUT_CRAM"
./pfvarx tool samtools view -T "$CRAM_REF" -c "$INPUT_CRAM"
```

**Step 2. Register and validate the CRAM route**

Set `input_type: CRAM` and put the CRAM path in the TSV `bam` column. Confirm coordinate sort, CRAI/index and read groups. Follow R03–R04, then run:

```bash
./pfvarx pipeline --mode NEW
SAMPLE=S01
cat "results/qc/samples/$SAMPLE/$SAMPLE.entry_contract.json"
./pfvarx tool samtools quickcheck -v "results/alignments/$SAMPLE.bam"
```

### 4. Expected outputs

BAM/BAI derivative and entry JSON; the original CRAM remains unchanged. The configured parasite/combined reference must support actual decoding identity.

### 5. Verification and quality control

Inspect decoded record counts, contig/reference identity and qualification agreement with a producer-validated BAM equivalent. Standard final output is BAM, not newly generated CRAM.

### 6. Error handling and troubleshooting

Reference MD5/sequence mismatch: stop and obtain the exact decoding reference. Do not substitute the nearest assembly or disable checking. Targeted CRAM also must meet R09 requirements.

### 7. Completion criteria

Reference-dependent conversion and branch qualification accepted and recorded before routine use.


<a name="procedure-r11"></a>

## R11 — Select automatic or explicit downstream analyses

### 1. Procedure identification

**Module:** FO004 analysis configuration and eligibility. **Workflow:** routine analysis.

**Purpose and scope:** Request actual supported scientific methods with required metadata and prerequisites.

### 2. Prerequisites and requirements

R01–R03 and chosen branch. Explicit FO004 needs a matching GFF. In automatic FO004, the GFF is required for annotation and marker outcomes. DA002 requires a genuinely verified catalogue; DA007 currently implements DA002-derived outcomes. DA005/DA006 require WGS.

### 3. Detailed operational procedure

**Step 1. Choose automatic or explicit execution**

For automatic execution, bootstrap a matching downstream workspace with
`--auto-analyses`. Use `selection: auto` and leave `requested` and `scales` empty.
Appendix G lists the candidates, gates, expected artifacts and skip reasons.
Use explicit selection for a subset or a particular analysis scale.

**Step 2. Configure explicit selections when needed**

Before R01 choose, for example, `BOOTSTRAP_ARGS=(--final-output FO004 --design wgs --analyses DA003 DA001 DA004 DA005 DA006)`. For DA003 only, select just DA003. If a sample-only workspace already exists, bootstrap a new matching downstream workspace and register the same original read paths; do not edit installation metadata to bypass its capability guard.

From the relevant workspace, preserve all generated keys and edit:

```bash
cd "$ANALYSIS"
nano config/config.yaml
```

For one diversity analysis:

```yaml
final_output: FO004
analyses:
  selection: explicit
  requested: [DA003]
  scales: {DA003: [POPULATION]}
  population_column: population
  resistance_catalog: null
  spatial_outcome_source: null
```

For combined WGS diversity/COI/relatedness/structure/tree work:

```yaml
analyses:
  selection: explicit
  requested: [DA003, DA001, DA004, DA005, DA006]
  scales:
    DA003: [POPULATION]
    DA004: [PAIRWISE]
    DA005: [BETWEEN_POPULATIONS]
    DA006: [COHORT]
  population_column: population
  resistance_catalog: null
  spatial_outcome_source: null
```

Also set the real matching `reference.gff`. For resistance/spatial work use `[DA002, DA007]`, a verified `resistance_catalog`, `spatial_outcome_source: DA002` and verified country/date/latitude/longitude metadata. On targeted data omit DA005/DA006. DA001 is optional for explicit selection; it is inserted automatically as an internal prerequisite when another selected method needs it.

**Step 3. Validate selections; launch only in the chosen entry SOP**

```bash
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
```

**Step 4. Inspect every selected method**

After the selected entry branch executes its single NEW run, list manifests with `find results/downstream -name status_manifest.json -print`. Read `results/eligibility/cohort_decisions.json`, `sample_analysis_eligibility.tsv` and each DA status manifest. Record requested, eligible, analysed, excluded and reasons. Do not equate requested count with analysed count.

### 4. Expected outputs

Section 4 lists outputs and their interpretation. Each requested method records a terminal manifest.

### 5. Verification and quality control

Check current registry/scales and parameter thresholds rather than inventing acceptance values. DA006 currently requires at least four eligible monoclonal samples and ten informative sites, subject to further missingness gates; these software floors do not establish a scientifically adequate study size. DA004 needs 200 informative markers per WGS pair (Taylor 2019) and 20 per targeted pair; pairs with fewer than 200 are `ESTIMATED_LOW_PRECISION` and carry no relatedness label, and DA004 reports SKIPPED_INELIGIBLE with `INSUFFICIENT_INFORMATIVE_SITES` when no pair can be estimated. If no sample passes QC into the cohort, every requested method reports SKIPPED_INELIGIBLE with `NO_ELIGIBLE_SAMPLES` and `NO_COHORT`.

### 6. Error handling and troubleshooting

BLOCKED_MISSING_PREREQUISITE: supply the real required resource or change the study request before a new run. SKIPPED_INELIGIBLE: report it; do not reduce gates merely to force a tree or comparison.

### 7. Completion criteria

Every requested method has an accepted truthful terminal state, with source outputs/uncertainty and all exclusions documented.



<a name="procedure-r12"></a>

## R12 — Validate every final-output mode and review QC

### 1. Procedure identification

**Module:** results/ and run status. **Workflow:** routine analysis.

**Purpose and scope:** Accept files and biological interpretation only within the selected final-output contract.

### 2. Prerequisites and requirements

Completed or terminal run from the chosen branch; exact RUN_ID; sample list and expected FO001–FO004. Use actual existing run paths, not a different latest symlink accidentally.

### 3. Detailed operational procedure

**Step 1. Resolve the selected run and inspect status**

```bash
cd "$ANALYSIS"
export RUN_ID=YOUR_RECORDED_RUN_ID
export RUN="$(pwd)/runs/$RUN_ID"
export RES="$RUN/results"
test -d "$RES"
cat "$RUN/status.json"
./pfvarx run-verify --run-id "$RUN_ID"
find "$RES" -maxdepth 3 -type f | sort
```

**Step 2. Check each admitted sample, including the caller status distribution**

```bash
SAMPLE=S01
./pfvarx tool samtools quickcheck -v "$RES/alignments/$SAMPLE.bam"
cat "$RES/qc/samples/$SAMPLE/qc_decision.json"
# The following checks apply only to FO002+ samples admitted to calling.
./pfvarx tool bcftools view -h "$RES/variants/samples/$SAMPLE.vcf.gz"
./pfvarx tool bcftools query -f '%FILTER\n' "$RES/variants/samples/$SAMPLE.vcf.gz" | sort | uniq -c
./pfvarx tool bcftools view -f PASS -H "$RES/variants/samples/$SAMPLE.vcf.gz" | wc -l
head "$RES/variants/samples/$SAMPLE.evidence.tsv"
```

**Step 3. Review cohort/annotation layers only when selected**

```bash
# FO003/FO004 only:
cat "$RES/cohort/cohort_status.json"   # COMPLETED, or SKIPPED_NO_ELIGIBLE_SAMPLES with reasons
./pfvarx tool bcftools query -l "$RES/cohort/cohort.vcf.gz"
gzip -t "$RES/cohort/evidence_matrix.tsv.gz"
head "$RES/cohort/missingness.tsv"
# FO004 only:
cat "$RES/annotation/annotation_qc.json"
```

The annotated `cohort.ann.vcf.gz` may be ordinary gzip rather than BGZF; do not assume a tabix index exists. Read it with gzip-aware text/appropriate tools. Inspect ANN and upstream QC instead of indexing original annotations ad hoc.

### 4. Expected outputs

The complete runtime layout is listed in Appendix A. Mode-specific final products:

| Mode | Required scope |
| --- | --- |
| FO001 | BAM/BAI, per-sample Layer 1 QC, run status/logs and report; no caller VCF is expected |
| FO002 | FO001 processing plus admitted sample VCF/TBI, evidence, callability/focal loss and caller QC |
| FO003 | FO002 plus variant-only cohort VCF, sample/site lists, forced evidence matrix and missingness |
| FO004 | FO003 plus matching-reference annotation, eligibility and requested/internal downstream terminal manifests; scientific outputs conditional on eligibility |

VCF variants include non-PASS statuses. Count accepted calls separately. Candidate totals are not accepted-SNP totals. The VCF GT describes the dominant allele; WSAF may indicate a minor alternative even when GT is 0.

### 5. Verification and quality control

Accept technical integrity and scientific QC separately. PASS, WARNING, FAIL, NOT_APPLICABLE, NOT_EVALUATED and UNKNOWN are distinct. Report failed/blocked samples and method skips; do not hide denominator changes.

### 6. Error handling and troubleshooting

Missing artifact: first check mode, QC eligibility and run ID; do not interpret missing VCF as zero variants. FAILED_EXECUTION requires log diagnosis; scientific skip can be legitimate. Reject unexplained integrity failure.

### 7. Completion criteria

Exact output contract and per-sample/method QC reviewed; all limitations, source hashes and decisions recorded before release of findings.


<a name="procedure-r13"></a>

## R13 — Diagnose, resume or extend safely

### 1. Procedure identification

**Module:** runs/ and pipeline modes. **Workflow:** routine analysis.

**Purpose and scope:** Recover an interrupted unchanged experiment or create a new justified extension without altering prior evidence.

### 2. Prerequisites and requirements

R12; original run ID and input/source checksum records; no concurrent duplicate job. Scientific changes require NEW rather than assuming cache/reuse validity.

### 3. Detailed operational procedure

**Step 1. Inspect the cause and preserve logs**

```bash
cd "$ANALYSIS"
export RUN_ID=YOUR_RECORDED_RUN_ID
cat "runs/$RUN_ID/status.json"
ls "runs/$RUN_ID/logs"
cat "runs/$RUN_ID/logs/events.jsonl"
find "runs/$RUN_ID/logs" -name '*.log' -type f -print
```

**Step 2. Select the correct recovery mode**

For the **same** source/configuration/inputs after a genuine interruption:

```bash
./pfvarx run-verify --run-id "$RUN_ID"
./pfvarx pipeline --mode RESUME --run-id "$RUN_ID"
```

For a declared cohort extension with unchanged scientific meaning, preserve old IDs/inputs, register only justified additional samples, validate and inspect the reuse basis before:

```bash
./pfvarx pipeline --mode EXTEND --run-id "$RUN_ID"
```

For changed reference, mask, reads, primers, policy or scientific preprocessing, use a new workspace/NEW experiment; do not rely on incomplete external-input reuse keys.

### 4. Expected outputs

Additional attempt logs/status for RESUME; a new run for EXTEND/NEW. Older evidence remains retained.

### 5. Verification and quality control

Restore a transient resource/system issue only after diagnosis. Verify frozen integrity and original external input checksums. Check all final products again using R12.

### 6. Error handling and troubleshooting

Unchanged caller log does not itself prove a hang. Check CPU/RAM/disk and workflow status before interruption. Do not delete lockfiles while a job is still running or edit frozen manifests.

### 7. Completion criteria

Recovery basis and exact commands documented; original record retained; recovered/new results pass R12.


<a name="procedure-r14"></a>

## R14 — Operate and export the reporting dashboard

### 1. Procedure identification

**Module:** reporting/ and report CLI. **Workflow:** routine analysis.

**Purpose and scope:** Review actual results offline and export selected data without changing scientific analyses.

### 2. Prerequisites and requirements

R12 accepted or explicitly incomplete run; modern local browser; current reporting assets in the workspace. No server/CDN required.

### 3. Detailed operational procedure

**Step 1. Refresh only reporting and open the HTML**

```bash
cd "$ANALYSIS"
./pfvarx report --run-id "$RUN_ID"
printf '%s\n' "$ANALYSIS/runs/$RUN_ID/results/report/pfvarx_dashboard.html"
```

Open the printed file path in the browser. This overwrites the report HTML/JSON and centralized QC JSON only; archive their previous bytes first if required for your audit.

**Step 2. Review the eight sections and export**

Use overview, downstream DA001–DA007, variant/annotation, QC, resources, inputs, software/provenance and logs. Choose one visualization per group. Filter actual fields; inspect active counts and click/focus a sample/mark for its source evidence. CSV/XLSX exports include all filtered embedded rows; tables display 50 rows/page. PNG exports the actual selected chart. Export PDF opens native print; all sections print, but tables use the visible page and downstream uses the current selected analysis. Print each desired analysis for a complete supplement.

**Step 3. Check provenance and sharing scope**

The variant embedding cap is currently 200000; full summaries precede that cap, detailed charts use the labelled embedded subset. Missing fields stay missing. Synthetic fixtures display SYNTHETIC DEMO DATA. Redacted log tails are bounded. Table pseudonymization is not full report anonymization; review source paths, free text and identifiers before sharing.

### 4. Expected outputs

`report/pfvarx_dashboard.html`, `report/dashboard_data.json`, `qc/centralized_quality_matrix.json`; user-selected CSV/XLSX/PNG/PDF and redacted-log downloads outside source.

### 5. Verification and quality control

Confirm run/reference/sample identities, QC source metrics, true denominators, view/export filters and retained renderer hash. Resource wall-time sums are not elapsed run duration; absent measurements are not zero.

### 6. Error handling and troubleshooting

Invalid local JSON: use the actual schema-compatible dashboard_data.json; do not insert fake values. A skipped method cannot be made present by switching views. Failed export: inspect visible feedback/browser restrictions.

### 7. Completion criteria

Report and exports match the selected run and intended subset; sharing review and report generation time recorded.


<a name="procedure-r15"></a>

## R15 — Archive a study and maintain the software boundary

### 1. Procedure identification

**Module:** study records, input checksums and frozen runs. **Workflow:** routine analysis.

**Purpose and scope:** Preserve reproducible scientific evidence separately from pristine distributable software.

### 2. Prerequisites and requirements

R12–R14; storage/backup policy, data-access authorization and actual input provenance. Do not remove study evidence merely to clean a software repository.

### 3. Detailed operational procedure

**Step 1. Record checksum evidence and integrity before archiving**

```bash
cd "$ANALYSIS"
./pfvarx run-verify --run-id "$RUN_ID" > "$RECORDS/analysis/frozen-integrity.json"
sha256sum pixi.lock config/config.yaml config/quality.yaml config/compute.yaml > "$RECORDS/analysis/final-config.sha256"
sha256sum "$RES/report/pfvarx_dashboard.html" "$RES/report/dashboard_data.json" > "$RECORDS/analysis/report.sha256"
```

**Step 2. Create a deliberate run archive outside source**

```bash
tar -czf "$RECORDS/analysis/run-$RUN_ID.tar.gz" "runs/$RUN_ID"
sha256sum "$RECORDS/analysis/run-$RUN_ID.tar.gz" > "$RECORDS/analysis/run-$RUN_ID.tar.gz.sha256"
sha256sum -c "$RECORDS/analysis/run-$RUN_ID.tar.gz.sha256"
```

Archive original external input/truth/mask checksum manifests and the wet-lab metadata separately; the frozen source does not guarantee their preservation. Protect identifiable location/sample metadata under the study policy.

### 4. Expected outputs

Run archive/checksum and study protocol, configuration, input identities, acceptance/exclusion records and report provenance, all outside the source release.

### 5. Verification and quality control

Test the archive inventory with `tar -tzf`; verify backup restoration according to the site policy. Confirm raw inputs still exist at approved locations or are independently archived.

### 6. Error handling and troubleshooting

Checksum/backup failure: stop disposal and repair backup. Never clean the only copy of a real sample or frozen run. Software packaging must use scripts/package.py, not tar of a live study workspace.

### 7. Completion criteria

Archive/restore and access controls accepted; only independently backed-up disposable intermediates may be removed under site policy.


<a name="procedure-r16"></a>

## R16 — Discover and register inputs without copying originals

### 1. Procedure identification

**Module:** inputs/organize.py. **Workflow:** routine analysis/benchmarking. Generate
reviewed manifests and reference/link registrations using the implemented CLI.

### 2. Prerequisites and requirements

Installed selected workspace; real unambiguous R1/R2 filenames or explicit TSV;
valid assay/input config. Original files must exist. Register only reviewed wet-lab IDs.

### 3. Detailed operational procedure

**Step 1. Discover candidates or supply an explicit sheet.**

```bash
cd "$ANALYSIS"
./pfvarx discover-inputs --directory /absolute/path/to/fastq
nano metadata/discovered_samples.tsv
```

Discovery supports SAMPLE_R1/R2, optional _L001 and _001, FASTQ/FQ with optional gzip.
It rejects orphan/ambiguous pairs and reports ignored names. Sample/library IDs inferred
from filenames are candidates only; replace them from actual lab records and add required metadata.

**Step 2. Publish reviewed registration without duplicating reads.**

```bash
./pfvarx register-inputs --manifest metadata/discovered_samples.tsv --mode reference
./pfvarx validate
```

Use `--mode link` for organized workspace symlinks. Imported relative paths resolve
from the registration manifest's parent. A populated sheet is not overwritten.

### 4. Expected outputs

Candidate TSV, final configured sample sheet, optional symlinks and
metadata/registration.json. Originals are neither copied nor modified.

### 5. Verification and quality control

Review sample/lane/library mapping and all file paths. Validation must accept current
pairs/reference/metadata. Candidate generation alone is not validation or identity proof.

### 6. Error handling and troubleshooting

Unmatched/duplicate mate: resolve delivery/naming or supply an explicit TSV. Existing
registered sheet: review/edit it deliberately or create a new workspace. Do not force overwrite.

### 7. Completion criteria

Reviewed registration and accepted validation recorded, with explicit original paths
and any generated-identifier corrections. Continue to preparation only after review.

<a name="procedure-r17"></a>

## R17 — Verify prepared artifacts before execution

### 1. Procedure identification

**Module:** references/prepare.py and execution/preflight.py. **Workflow:** all execution.
Detect missing/stale/damaged required derivatives before creating a run.

### 2. Prerequisites and requirements

Configured/registered real inputs; current preparation; native locked tools. FASTQ needs
BWA indexes; aligned input does not. Original resource fingerprints must be available.

### 3. Detailed operational procedure

**Step 1. Prepare only required derivatives.**

```bash
cd "$ANALYSIS"
./pfvarx prepare --dry-run
./pfvarx prepare
./pfvarx verify-prepared
```

**Step 2. Inspect the actual verification record and launch decision.**

```bash
cat metadata/prepared_verification.json
cat metadata/preparation_report.json
./pfvarx ready
./pfvarx dry-run
```

A non-PASS preparation verification stops execution. `prepare` repairs missing/damaged
artifacts; no operator should edit hashes to claim a damaged file is valid.

### 4. Expected outputs

Sealed prepared manifest schema pfvarx.prepared/2, source/derived hashes, compatibility
result, current input-validation/quality/configuration reports. No scientific run yet.

### 5. Verification and quality control

Confirm nonempty derivatives, checksum agreement and FASTA/FAI/dictionary names/lengths.
A checksum seal is not an independent authenticity certificate; preserve reviewed source inputs.

### 6. Error handling and troubleshooting

Corruption/stale source: inspect reported paths and rerun preparation after resolving
underlying storage/source issues. Never bypass mandatory checks or delete original references.

### 7. Completion criteria

verify-prepared PASS, ready true and dry-run OK, with any PARTIAL sample exclusions
explicitly accepted. Execution independently enforces these preparation requirements.

<a name="procedure-r18"></a>

## R18 — Orchestrate the complete chain and apply retention safely

### 1. Procedure identification

**Module:** execution/preflight.py, runner.py and retention.py. **Workflow:** all native
analyses. Run ordered controls and retain/delete only the configured intermediate scope.

### 2. Prerequisites and requirements

Selected installed capabilities, reviewed config/registration and real resources.
No active overlapping execution. Choose keep versus successful retention before launch.

### 3. Detailed operational procedure

**Step 1. Select policy and optional preparation-only boundary.**

```yaml
retention:
  intermediates: successful
```

Use `keep` for debugging. `successful` is the default. Then:

```bash
cd "$ANALYSIS"
./pfvarx orchestrate --through verify
```

**Step 2. Execute the full ordered chain once.**

```bash
./pfvarx orchestrate --through execute
./pfvarx status
```

This performs initialize → validate → prepare → verify → preflight → execute/report.
It stops at mandatory failures; do not also launch a redundant pipeline NEW for the same intended run.

**Step 3. Preview explicit completed-run cleanup and apply only after review.**

```bash
./pfvarx cleanup --run-id YOUR_COMPLETED_RUN_ID
./pfvarx cleanup --run-id YOUR_COMPLETED_RUN_ID --apply
```

### 4. Expected outputs

The ordinary frozen run/scientific/report paths, stage metadata and an applied
cleanup.json with removed regular work files/bytes. keep retains Snakemake temporary files.

### 5. Verification and quality control

Workspace locking rejects overlapping execution/cleanup. Bulk cleanup requires
COMPLETED and stays under that run's results/work; final BAM/VCF/evidence/QC/frozen/log/report
and original input/reference paths are preserved. Symlinks are excluded.

### 6. Error handling and troubleshooting

Failed/interrupted run: preserve work/logs and diagnose/resume under R13; no final bulk
cleanup. Active workspace: wait, do not delete a live lock or launch duplicates.

### 7. Completion criteria

Run/verification/report accepted; retention outcome reviewed against policy. Original
and essential result identities remain unchanged and all deletion evidence is retained.

# PART III — Built-in benchmarking and performance evaluation

Use a completely separate workspace. These SOPs describe the implemented evaluator,
its exact-normalized pure-sample matching, mixture calibration, protected partition,
resource summaries and reproducibility checks. They do not add external callers,
new truth adapters, unsupported statistical tests or a biological-MOI validator.
Use the same scientific caller implementation as routine analysis.


<a name="procedure-b01"></a>

## B01 — Create a separate benchmark workspace and select inputs

### 1. Procedure identification

**Module:** benchmark initialization. **Workflow:** benchmarking.

**Purpose and scope:** Declare an independently justified dataset and install the same implemented software in a separate location.

### 2. Prerequisites and requirements

Common variables; pristine source; original sample/reference provenance; independent truth to be supplied. `$BENCH` must be new. Required reads, reference/GFF/core resources follow R02–R03.

### 3. Detailed operational procedure

**Step 1. Create and initialize the benchmark copy**

```bash
test ! -e "$BENCH"
BENCH_BOOTSTRAP_ARGS=(--profile benchmark --design WGS --input-type FASTQ)
# Add --analyses DA003 DA001 ... BEFORE bootstrap when FO004 is intended.
python3 "$PFVARX_SOURCE/bootstrap.py" "${BENCH_BOOTSTRAP_ARGS[@]}" --workspace "$BENCH"
cd "$BENCH"
nano config/config.yaml
./pfvarx workspace-init
```

**Step 2. Choose the experiment scope before acquisition**

Set a unique project ID, real assay/amplification route, FO002 for caller-only experiments or FO004 with a matching GFF and R11 selections for downstream evaluation. Choose independently characterized pure samples and mixtures spanning the intended depth, minor-fraction, homopolymer/low-complexity and WGS/targeted contexts. Record actual supplier/accession/release/license/checksum facts; no benchmark accession is invented by this manual.

**Step 3. Register acquired resources**

Use R02–R03 with `$BENCH` as the workspace (`export ANALYSIS="$BENCH"` for those commands). Original study inputs may be registered read-only by absolute path. Do not copy or run a fixture over real data. Save a protocol and source checksum manifest outside the source distribution.

### 4. Expected outputs

Benchmarking-purpose config, sample registry and runtime `inputs/truth/` and `inputs/partitions/` folders; no truth/calls yet.

### 5. Verification and quality control

Check assay, actual sample IDs and resource identity. A previously run field pilot without independent truth cannot be converted into an accuracy reference by declaring its calls true.

### 6. Error handling and troubleshooting

Missing independent truth or uncertain mixture composition: stop the corresponding performance claim. Synthetic smoke data is permitted only under B02 and must remain labelled.

### 7. Completion criteria

Dataset/source identities and intended performance domain recorded; workspace is separate from routine analysis.


<a name="procedure-b02"></a>

## B02 — Run the supported synthetic software smoke test

### 1. Procedure identification

**Module:** scripts/verify.py / fixtures/synthetic.py. **Workflow:** benchmarking software verification.

**Purpose and scope:** Verify pipeline integration using on-demand synthetic data without calling it biological accuracy.

### 2. Prerequisites and requirements

Installed full source/Pixi tools; new external output directory and adequate free disk; no real input file will be overwritten.

### 3. Detailed operational procedure

**Step 1. Run the isolated built-in verification scenario**

```bash
# Verification of all scientific methods uses a separate full/development tool copy.
python3 "$PFVARX_SOURCE/bootstrap.py" --profile development --workspace "$LAB/verification-tools"
cd "$LAB/verification-tools"
./pfvarx tool python scripts/verify.py benchmark --output "$LAB/benchmark-smoke"
```

**Step 2. Inspect its acceptance and expected exclusion**

```bash
cat "$LAB/benchmark-smoke/benchmark/verification.json"
cat "$LAB/benchmark-smoke/benchmark/benchmarks/development/benchmark_report.json"
find "$LAB/benchmark-smoke/benchmark/results/qc/samples" -name qc_decision.json -print
```

The runner creates an informative WGS fixture, a low-depth exclusion, all requested DA manifests, development/held-out evaluation and the expected refusal of repeated held-out evaluation. Its test manifest must report PASSED.

### 4. Expected outputs

Only the new external smoke workspace contains synthetic FASTQ/reference/truth/calls/reports and evaluation outputs. Production source and real benchmark inputs remain unchanged.

### 5. Verification and quality control

Check the SYNTHETIC DEMO DATA report banner, expected low-depth sample blocking and all method terminal states. A skipped method is not evidence it executed; inspect the informative fixture manifests.

### 6. Error handling and troubleshooting

Output already exists: choose another new directory. Tool or data failure: inspect the numbered command log and run tool logs; do not overwrite the failed evidence with a successful-looking fixture.

### 7. Completion criteria

Software integration verified; record expressly that scores are synthetic plumbing checks, not field/clinical accuracy. Clean disposable smoke data only after needed evidence is saved.


<a name="procedure-b03"></a>

## B03 — Prepare independent truth and confidence regions

### 1. Procedure identification

**Module:** inputs/truth/ and reference validation. **Workflow:** benchmarking.

**Purpose and scope:** Supply actual truth in the evaluator format and independently inspect its identity/completeness.

### 2. Prerequisites and requirements

B01; actual pure-sample VCFs and, for mixtures, independently established clone variants/DNA fractions; exact assembly FASTA; confidence BED; registered sample IDs.

### 3. Detailed operational procedure

**Step 1. Populate required truth paths and inspect format**

```bash
cd "$BENCH"
ls inputs/truth
python3 -m json.tool inputs/truth/truth.json
head inputs/truth/confident_regions.bed
```

Required layout: `truth.json`, `confident_regions.bed`, and `<SAMPLE>.truth.vcf` for every pure sample. Pure truth VCFs are uncompressed `.vcf` in this adapter. `truth.json` is a JSON object with `samples` and `clones`. Each sample has `clonality`, `true_k`, and `components` containing `clone`/`proportion`. Each clone variant has `chrom`, `pos`, `ref`, `alt`. For example, the following is **format illustration only**, not usable truth:

```json
{"samples":{"PURE_A":{"clonality":"MONOCLONAL","true_k":1,"components":[{"clone":"cloneA","proportion":1.0}]}},"clones":{"cloneA":[]}}
```

An empty clone list does not prove a real clone is reference everywhere. Mixture scoring requires complete constituent-clone variant lists and actual proportions summing to one. Use VCF/JSON one-based positions and BED zero-based half-open coordinates.

**Step 2. Validate source identity without changing truth**

```bash
./pfvarx validate
./pfvarx prepare
TRUTH_SAMPLE=PURE_A
./pfvarx tool bcftools view -h "inputs/truth/$TRUTH_SAMPLE.truth.vcf"
./pfvarx tool bcftools norm --check-ref e -f prepared/reference/parasite.fasta "inputs/truth/$TRUTH_SAMPLE.truth.vcf" -Ou -o /dev/null
sha256sum inputs/truth/truth.json inputs/truth/confident_regions.bed "inputs/truth/$TRUTH_SAMPLE.truth.vcf" > "$RECORDS/truth.sha256"
```

**Step 3. Perform independent truth review**

Review every sample ID, clone identifier, DNA fraction, true K definition, VCF REF/ALT and complete confidence interval set. Check BED contig membership, nonnegative/in-range bounds and sorted nonoverlapping regions. The built-in pipeline validation is not a complete truth validator. Require a second reviewer or documented independent truth protocol before scoring.

### 4. Expected outputs

Actual supplied truth resources at the conventional paths with provenance/checksums and review record; no inferred truth from PfVarX output.

### 5. Verification and quality control

Every pure truth sample has a real matching VCF; mixture components are fully declared and sample membership is valid. The evaluator uses one shared confidence BED; experiments requiring incompatible per-sample masks must be separated.

### 6. Error handling and troubleshooting

REF mismatch, incomplete truth or uncertain proportions: stop. Do not use `--check-ref s` or swap reference alleles to hide disagreement. Unknown true K must not be guessed from the caller.

### 7. Completion criteria

Independent truth/reference/coordinate/completeness review accepted; inputs are frozen for partitioning.


<a name="procedure-b04"></a>

## B04 — Freeze the protocol and lock the development/held-out partition

### 1. Procedure identification

**Module:** benchmarking/partitions.py. **Workflow:** benchmarking.

**Purpose and scope:** Prevent tuning on held-out outcomes and record the experiment before calling/evaluation.

### 2. Prerequisites and requirements

B01/B03; enough independent samples for the proposed inference; documented lineage/replicate relationships and predeclared success criteria.

### 3. Detailed operational procedure

**Step 1. Write the protocol before seeing results**

In `$RECORDS`, record assay, dataset, source/lock, reference/core/confidence masks, truth, preprocessing, caller parameter source classes, acceptable error rates/intervals for the intended use, missing-sample policy and reporting strata. Acceptance limits must be justified by the study; this manual supplies no universal clinical sensitivity threshold.

**Step 2. Lock partition once**

```bash
cd "$BENCH"
./pfvarx benchmark-partition --holdout-fraction 0.25 --seed 20261106
cat inputs/partitions/partition.json
sha256sum inputs/partitions/partition.json > "$RECORDS/partition-file.sha256"
```

**Step 3. Review membership before tuning**

The 0.25 fraction and seed are explicit illustrative choices and must be predeclared, not searched for favourable results. The built-in split is sample-random, not family/lineage-aware. Prevent biological leakage through sample selection and independent study design; different IDs do not make technical replicates independent. Holdout access control is a protocol responsibility.

### 4. Expected outputs

Locked `inputs/partitions/partition.json` with development (`benchmarking`) and held-out IDs, seed/content hash.

### 5. Verification and quality control

Check both groups have the planned composition and adequate independent truth support. Verify there was no prior outcome-driven split selection.

### 6. Error handling and troubleshooting

Existing lock or hash mismatch: stop and audit experiment identity; do not delete/edit it to unlock another split. If design is invalid, document a genuinely new experiment before exposing outcomes.

### 7. Completion criteria

Protocol and partition approved/frozen; held-out outcomes remain unseen while development proceeds.


<a name="procedure-b05"></a>

## B05 — Execute the benchmark pipeline and score development once

### 1. Procedure identification

**Module:** pipeline and benchmark-evaluate. **Workflow:** benchmarking.

**Purpose and scope:** Call without truth-informed tuning and produce the implemented development reports.

### 2. Prerequisites and requirements

B01/B03/B04; selected mode/analyses configured before launch; prepared resources and realistic compute budget. Use explicit run IDs when multiple runs exist.

### 3. Detailed operational procedure

**Step 1. Preflight and run**

```bash
cd "$BENCH"
./pfvarx check-layout
./pfvarx validate
./pfvarx prepare --dry-run
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline --mode NEW
./pfvarx status
export BENCH_RUN_ID=YOUR_COMPLETED_BENCHMARK_RUN_ID
./pfvarx run-verify --run-id "$BENCH_RUN_ID"
```

**Step 2. Review sample QC before scoring**

Apply R12 in this workspace. Record input/Layer 1 exclusions, missing VCFs and caller status distributions. Accuracy scoring is for pure MONOCLONAL truth samples; mixture calibration is separate.

**Step 3. Generate development evaluation outputs**

```bash
./pfvarx benchmark-evaluate --run-id "$BENCH_RUN_ID" --out benchmarks/development
cat benchmarks/development/benchmark_report.json
head benchmarks/development/per_sample_accuracy.tsv
head benchmarks/development/mixture_sites.tsv
head benchmarks/development/k_selection.tsv
```

### 4. Expected outputs

`benchmarks/development/{benchmark_report.json,per_sample_accuracy.tsv,mixture_sites.tsv,k_selection.tsv}` and normalized call/truth VCF intermediates. Downstream/report outputs remain under the selected run.

### 5. Verification and quality control

Report partition must be `benchmarking` and run_id must match. Review `NO_CALLS` rows separately; omitted/missing samples are not perfect negative predictions. Subsequent SOPs inspect these same files rather than repeatedly overwriting them.

### 6. Error handling and troubleshooting

Incomplete run or missing expected pure VCF/truth: diagnose before scoring. Preserve evaluation evidence for each configuration in a distinct predeclared directory; do not repurpose held-out folders for development.

### 7. Completion criteria

Completed run verified, development outputs created and denominators/exclusions documented.


<a name="procedure-b06"></a>

## B06 — Review SNP detection accuracy

### 1. Procedure identification

**Module:** built-in exact-normalized SNV evaluator. **Workflow:** benchmarking.

**Purpose and scope:** Inspect pure-sample SNV agreement using actual TP/FP/FN and ratio intervals.

### 2. Prerequisites and requirements

B05 development evaluation; truth confidence scope accepted. This SOP reads existing outputs and does not re-expose held-out data.

### 3. Detailed operational procedure

**Step 1. Inspect pure SNV rows and summary**

```bash
cd "$BENCH"
awk -F '\t' 'NR==1 || $2=="SNV"' benchmarks/development/per_sample_accuracy.tsv
./pfvarx tool python - <<'PYCODE'
import json
d=json.load(open("benchmarks/development/benchmark_report.json"))
print(json.dumps(d["accuracy"]["SNV"],indent=2))
PYCODE
```

**Step 2. Review mismatch representation**

Inspect normalized `.calls.norm.vcf.gz` and `.truth.norm.vcf.gz` for a selected scored sample using `./pfvarx tool bcftools view -H`. Check sample identity, alleles, PASS filtering and confidence-region anchors before interpreting discordances. Stratify interpretation by the declared assay/depth/context; the built-in summary does not automatically implement every proposed stratum.

### 4. Expected outputs

Existing per-sample SNV TP/FP/FN, sensitivity, precision, FP/callable-Mb and summary uncertainty; interpretation record linked to source rows.

### 5. Verification and quality control

Sensitivity is TP/(TP+FN); precision TP/(TP+FP). Read actual uncertainty metadata. Zero denominators are unavailable, not zero. Exact normalized identity is the implemented match rule.

### 6. Error handling and troubleshooting

Unexpected perfect/poor scores: verify truth/reference, sample correspondence, masks and exclusions first. Do not remove difficult loci after looking at held-out performance.

### 7. Completion criteria

SNV metrics accepted only for the recorded truth/assay/mask domain, with explicit sample and site denominators.


<a name="procedure-b07"></a>

## B07 — Review indel accuracy and complex-region limitations

### 1. Procedure identification

**Module:** built-in exact-normalized INDEL evaluator. **Workflow:** benchmarking.

**Purpose and scope:** Inspect small-indel agreement without treating exact identity as universal haplotype equivalence.

### 2. Prerequisites and requirements

B05; INDEL truth and confidence domains with independently reviewed repeat/homopolymer content.

### 3. Detailed operational procedure

**Step 1. Inspect indel rows and summary**

```bash
cd "$BENCH"
awk -F '\t' 'NR==1 || $2=="INDEL"' benchmarks/development/per_sample_accuracy.tsv
./pfvarx tool python - <<'PYCODE'
import json
d=json.load(open("benchmarks/development/benchmark_report.json"))
print(json.dumps(d["accuracy"]["INDEL"],indent=2))
PYCODE
```

**Step 2. Check normalization source files**

```bash
TRUTH_SAMPLE=PURE_A
./pfvarx tool bcftools view -H "benchmarks/development/$TRUTH_SAMPLE.calls.norm.vcf.gz"
./pfvarx tool bcftools view -H "benchmarks/development/$TRUTH_SAMPLE.truth.norm.vcf.gz"
```

Review left alignment/decomposition and representative homopolymer/repeat discordances against the exact reference. The evaluator groups non-SNV alleles as INDEL by allele lengths; inspect the declared small-variant scope so complex/multi-base events are not silently misinterpreted.

### 4. Expected outputs

Existing per-sample and pooled indel metrics plus mismatch review notes. No new haplotype-aware performance claim is generated by these commands.

### 5. Verification and quality control

Report exact normalized allele matching and its limitations for complex equivalent representations. Confidence anchoring and callability denominator definitions must accompany indel performance.

### 6. Error handling and troubleshooting

Reference/normalization warning: stop and investigate. Do not manually change indel positions to improve apparent concordance.

### 7. Completion criteria

Indel scope, normalization and difficult-region limitations documented separately from SNV findings.


<a name="procedure-b08"></a>

## B08 — Review low-frequency detection and WSAF calibration

### 1. Procedure identification

**Module:** mixture_sites.tsv and benchmark report. **Workflow:** benchmarking.

**Purpose and scope:** Assess declared mixture SNV fractions using implemented detection bins and observed intervals.

### 2. Prerequisites and requirements

B03/B05 with independently characterized mixture clone variants/DNA proportions. Detection scope is mixture SNVs, not a validated arbitrary low-frequency indel assay.

### 3. Detailed operational procedure

**Step 1. Inspect site-level mixture evidence and actual summary bins**

```bash
cd "$BENCH"
head -30 benchmarks/development/mixture_sites.tsv
./pfvarx tool python - <<'PYCODE'
import json
d=json.load(open("benchmarks/development/benchmark_report.json"))
for key in ["wsaf_calibration","detection_by_expected_fraction"]:
    print(key,json.dumps(d[key],indent=2))
PYCODE
```

**Step 2. Verify denominator and uncertainty interpretation**

Current expected-fraction bins are [0,0.1), [0.1,0.2), [0.2,0.35), [0.35,0.5), [0.5,0.999). A site counts detected-pass only when its actual status is PASS. Absent candidates are NOT_DETECTED, not WSAF=0. The calibration slope/intercept/RMSE/interval coverage uses observed sites only and requires at least three; report detection misses separately. Bias correction is not implemented; WSAF is a DNA fraction, not a cell-count fraction.

### 4. Expected outputs

Mixture-site expected/observed WSAF, interval and status table; bin counts/detection rates; conditional calibration summaries.

### 5. Verification and quality control

Compare expected fraction/depth/context with site evidence. Interval coverage conditional on observed sites is not unconditional coverage of every truth site. Do not infer a limit of detection from a small fixture or one successful minor allele.

### 6. Error handling and troubleshooting

Empty calibration: inspect eligible mixture sites and declared clone truth; do not fill missing coefficients. Unequal amplification/representation: record the limitation, do not reinterpret it as known correction.

### 7. Completion criteria

Detection and calibration results include missing/detected denominators, tested fraction/depth range and assay-specific limitations.


<a name="procedure-b09"></a>

## B09 — Evaluate technical mixture K agreement

### 1. Procedure identification

**Module:** k_selection.tsv and caller mixture evidence. **Workflow:** benchmarking.

**Purpose and scope:** Compare technical caller K against independently declared synthetic/experimental mixture composition.

### 2. Prerequisites and requirements

B03/B05 with actual true_k/composition. This is not a biological MOI/COI validation procedure.

### 3. Detailed operational procedure

**Step 1. Inspect declared and selected K plus fit status**

```bash
cd "$BENCH"
cat benchmarks/development/k_selection.tsv
./pfvarx tool python - <<'PYCODE'
import json
d=json.load(open("benchmarks/development/benchmark_report.json"))
print(json.dumps(d["k_selection"],indent=2))
PYCODE
```

**Step 2. Review the underlying fitted evidence**

```bash
SAMPLE=MIX_AB
cat "results/qc/samples/$SAMPLE/$SAMPLE.caller_qc.json"
```

Replace MIX_AB with an actual registered sample. Review selected K, fit stability, component weights, available status and missing/unidentified fits. The summary correct/total count is not a substitute for the per-sample status table.

### 4. Expected outputs

Existing true_k/technical_k/status/weights table and aggregate agreement; no genome-wide clone reconstruction or biological-MOI truth inference.

### 5. Verification and quality control

Check constituent-clone identifiability, sufficient informative loci and low minor fraction. Distinct laboratory clones can be difficult to distinguish at the assayed loci; do not assume agreement must be perfect.

### 6. Error handling and troubleshooting

Unavailable/unstable fit: report it rather than assigning K=1. Do not change true_k to match the estimate or reinterpret K as MOIRE COI.

### 7. Completion criteria

Agreement denominators, non-estimable cases and tested mixture design recorded explicitly.


<a name="procedure-b10"></a>

## B10 — Review implemented diversity and clonality outputs in benchmarks

### 1. Procedure identification

**Module:** FO004 DA003/DA001. **Workflow:** benchmarking.

**Purpose and scope:** Check the actual downstream numerical/eligibility outputs when the experiment requested them.

### 2. Prerequisites and requirements

B05 configured as FO004 with GFF and DA003/DA001 selected before NEW execution. For caller-only FO002 these files are not expected.

### 3. Detailed operational procedure

**Step 1. Inspect analysis outcomes and source tables**

```bash
cd "$BENCH"
cat results/downstream/DA003/status_manifest.json
cat results/downstream/DA001/status_manifest.json
head results/downstream/DA003/within_host_fws.tsv
head results/downstream/DA003/population_diversity.tsv
head results/downstream/DA001/clonality.tsv
```

**Step 2. Review uncertainty and informative-data scope**

Record actual eligible/analysed samples, FWS/diversity availability, COI posterior medians/bounds, classification and informative loci. Inspect conditional MOIRE log/output when executed. This is supported output QA; the built-in benchmark report does not provide a complete independent biological-MOI accuracy analysis. Do not claim that validation from technical K agreement.

### 4. Expected outputs

DA003/DA001 manifests and conditional tables already generated by the same workflow; benchmark QA notes linked to the actual run.

### 5. Verification and quality control

Few loci, low group size or poor identifiability can legitimately yield missing/indeterminate results. Check biological interpretability independently of software execution.

### 6. Error handling and troubleshooting

FO002 missing downstream files: expected, not failure. If a separate downstream benchmark is required, configure FO004 before a NEW experiment and retain its own run identity.

### 7. Completion criteria

Output QA completed without unsupported MOI accuracy claims; all missing/excluded cases retained.


<a name="procedure-b11"></a>

## B11 — Audit exact matching, confidence and callable denominators

### 1. Procedure identification

**Module:** benchmark matching intermediates and BEDs. **Workflow:** benchmarking.

**Purpose and scope:** Verify the implemented concordance rule and denominator choices before comparing performance claims.

### 2. Prerequisites and requirements

B05; sample truth/calls/confidence BED/callability BED; exact normalized VCF intermediates.

### 3. Detailed operational procedure

**Step 1. Inspect confidence anchors and callable span**

```bash
cd "$BENCH"
SAMPLE=PURE_A
head inputs/truth/confident_regions.bed
head "results/variants/samples/$SAMPLE.callable.bed"
awk 'BEGIN{s=0} !/^#/ && NF>=3 {s+=$3-$2} END{print "callable_bp",s; print "callable_Mb",s/1000000}' "results/variants/samples/$SAMPLE.callable.bed"
```

**Step 2. Record the exact evaluator semantics**

Calls are restricted to PASS and variants are compared as `(chrom,pos,ref,alt)` after bcftools normalization/decomposition. Confidence selection uses the variant anchor. Truth outside caller-callable intervals can still count as FN. FP-per-callable-Mb divides by the whole caller callable span, not automatically its intersection with confidence regions. Record these definitions beside every summary.

### 4. Expected outputs

Denominator audit notes and recorded callable bp/Mb, linked to exact VCFs/BEDs and the existing metric rows.

### 5. Verification and quality control

Recalculate a sample FP/callable-Mb from actual table FP and measured span. Require exact agreement within displayed rounding; do not silently switch denominator definitions between experiments.

### 6. Error handling and troubleshooting

Unexpected discrepancy: verify selected sample/run, BED overlap/length validity and metric source. Do not change the mask after seeing final results.

### 7. Completion criteria

Matching and all denominator scopes explicitly accepted and reproducible.


<a name="procedure-b12"></a>

## B12 — Review built-in summary metrics and uncertainty

### 1. Procedure identification

**Module:** benchmark_report.json and linearized-ratio wild cluster intervals. **Workflow:** benchmarking.

**Purpose and scope:** Interpret actual sensitivity/precision/F1 summaries and their supported uncertainty.

### 2. Prerequisites and requirements

B05–B11; enough independent pure truth samples; predeclared study acceptance criteria. No new competing-configuration statistical model is introduced.

### 3. Detailed operational procedure

**Step 1. Read summary uncertainty and independently check pooled count formulas**

```bash
cd "$BENCH"
./pfvarx tool python - <<'PYCODE'
import csv,json
rows=list(csv.DictReader(open("benchmarks/development/per_sample_accuracy.tsv"),delimiter="\t"))
for kind in ("SNV","INDEL"):
    selected=[r for r in rows if r["status"]=="SCORED" and r["variant_type"]==kind]
    if not selected:
        print(kind,"NO SCORED SAMPLES");continue
    tp,fp,fn=[sum(int(r[k]) for r in selected) for k in ("tp","fp","fn")]
    print(kind,dict(tp=tp,fp=fp,fn=fn,sensitivity=tp/(tp+fn) if tp+fn else None,
        precision=tp/(tp+fp) if tp+fp else None,f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None))
print(json.dumps(json.load(open("benchmarks/development/benchmark_report.json"))["accuracy"],indent=2))
PYCODE
```

**Step 2. Interpret the supplied intervals correctly**

The built-in method is a wild cluster bootstrap of the linearized ratio estimator over pure truth samples: Rademacher weights with at least ten clusters, Webb below that, and enumeration where small enough. It is not a paired superiority test between callers. Small samples limit inference. Missing metrics are not zero; if displayed F1 is absent, independently inspect TP/FP/FN before interpreting it. A positive denominator with TP=0 implies arithmetic F1=0, while an empty denominator remains unavailable.

### 4. Expected outputs

Existing ratio estimates/interval metadata and a transparent console check of actual pooled counts. No specificity value is computed by the built-in evaluator.

### 5. Verification and quality control

Inspect counts, number of scored samples, interval method and exclusions. Predeclared acceptance should consider uncertainty and intended assay domain, not only a point estimate.

### 6. Error handling and troubleshooting

Nonfinite/undefined summary: inspect missing samples, zero denominators and confidence scope. Do not invent true-negative counts or statistical significance.

### 7. Completion criteria

Summary, count arithmetic and uncertainty limitations recorded; study reviewer determines acceptance against the frozen protocol.


<a name="procedure-b13"></a>

## B13 — Measure resources and generate the supported profile

### 1. Procedure identification

**Module:** Snakemake benchmark directives / benchmark-resources. **Workflow:** benchmarking.

**Purpose and scope:** Summarize measured job resources and review the planner profile before reuse.

### 2. Prerequisites and requirements

B05 completed with available instrumented benchmark TSVs; known hardware, sample count, assay/depth and preparation/cache state.

### 3. Detailed operational procedure

**Step 1. Inspect raw measurements and write the profile**

```bash
cd "$BENCH"
find "runs/$BENCH_RUN_ID/benchmarks" -name '*.tsv' -type f -print
./pfvarx benchmark-resources --run-id "$BENCH_RUN_ID"
cat config/profiles/resources.yaml
```

**Step 2. Review measurements in the dashboard**

Open the corresponding run dashboard and select measured wall time, RSS, CPU time/average active cores or the table. Preserve source benchmark TSVs. Record hardware and workload basis with the generated profile; it is not automatically portable to another machine or assay.

### 4. Expected outputs

`config/profiles/resources.yaml` generated at runtime plus the original per-job benchmark TSVs. Available values include wall seconds, CPU seconds, peak RSS and IO traffic.

### 5. Verification and quality control

Sum of job wall times can exceed elapsed run time because jobs overlap. Peak job RSS is not a measured concurrent whole-workflow RAM peak. IO traffic is not peak disk usage. Missing measurements remain missing.

### 6. Error handling and troubleshooting

No benchmark records: report the unmeasured scope; do not estimate a measured peak. OOM/quota: inspect actual hardware/budgets and only adjust legitimate computational settings before a new run.

### 7. Completion criteria

Resource profile accepted only for a comparable stated dataset/hardware basis; no invented execution timeline or disk peak.


<a name="procedure-b14"></a>

## B14 — Measure sample-count scaling using ordinary supported runs

### 1. Procedure identification

**Module:** separate benchmarking workspaces and sample registry. **Workflow:** benchmarking performance.

**Purpose and scope:** Run predeclared cohort sizes without adding a new scalability engine.

### 2. Prerequisites and requirements

B01–B05; sufficient real samples, fixed software/scientific settings and compute budget; retain all lanes per selected sample. Use new workspace/experiment IDs and measured resource profiles for every size.

### 3. Detailed operational procedure

**Step 1. Predeclare cohort sizes and copy source per size**

Select actual available sample counts before measuring, e.g. 1, 4, 8 only if those cohorts exist. Keep the same assay, sample-depth distribution, reference/mask and compute budget. For each N set a new BENCH path and repeat B01/R02/R03 with exactly N unique sample IDs and every lane for those IDs. Use absolute resource/read paths so new workspace copies resolve the original files correctly.

```bash
export BENCH="$LAB/benchmark-scale-N4"
# Repeat B01/R02/R03 for this new directory with exactly four actual samples.
cd "$BENCH"
./pfvarx tool python - <<'PYCODE'
import csv
rows=list(csv.DictReader(open("prepared/samples.tsv"),delimiter="\t"))
ids=sorted({r["sample_id"] for r in rows})
assert len(ids)==4, ids
print("samples",ids,"lanes",len(rows))
PYCODE
```

**Step 2. Execute and collect actual resources for each N**

```bash
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
date -u +%FT%TZ
./pfvarx pipeline --mode NEW
date -u +%FT%TZ
./pfvarx status
# Set BENCH_RUN_ID from this run, then:
./pfvarx benchmark-resources --run-id "$BENCH_RUN_ID"
```

Record start/end timestamps and measured profile separately for every size. Declare whether reference preparation and Pixi startup are inside the timed interval, and whether caches are warm. Use actual QC-admitted as well as registered sample counts.

### 4. Expected outputs

One distinct frozen run and actual resource records per cohort size; an operator-maintained comparison table of N, admitted N, hardware, elapsed scope, measured job CPU/RSS and completion state.

### 5. Verification and quality control

Compare only like-for-like budgets/assays. Do not extrapolate a two-sample synthetic run into a whole-genome production guarantee. Report QC/admission changes as possible causes of scaling differences.

### 6. Error handling and troubleshooting

Unavailable N, sample overlap error or changed settings: reject that comparison and correct experiment registration before rerunning. Do not overwrite earlier runs/profiles as if all belonged to one size.

### 7. Completion criteria

Every size has its own validated run and workload/resource basis; scalability conclusions are limited to measured conditions.


<a name="procedure-b15"></a>

## B15 — Verify repeatability and frozen identity

### 1. Procedure identification

**Module:** NEW/RESUME modes and normalized scientific fields. **Workflow:** benchmarking reproducibility.

**Purpose and scope:** Check scientific result stability with fixed source, inputs and seed, distinguishing it from byte-identical timestamps.

### 2. Prerequisites and requirements

Completed run A with original input/reference/mask/truth checksums; fixed random_seed and source/lock. Independent repeat B uses NEW in an identical approved setup; RESUME instead tests recovery.

### 3. Detailed operational procedure

**Step 1. Record identity and execute the predeclared repeat**

```bash
cd "$BENCH"
./pfvarx run-verify --run-id "$BENCH_RUN_ID"
sha256sum pixi.lock config/config.yaml config/quality.yaml config/compute.yaml
./pfvarx pipeline --mode NEW
./pfvarx status
export RUN_A=FIRST_COMPLETED_RUN_ID
export RUN_B=REPEATED_COMPLETED_RUN_ID
```

**Step 2. Compare defined scientific fields, not generated timestamps**

```bash
SAMPLE=PURE_A
mkdir -p benchmarks/repeatability
for r in "$RUN_A" "$RUN_B"; do
  ./pfvarx tool bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%PST\t%WSAF\t%WSAFLO\t%WSAFHI]\n' "runs/$r/results/variants/samples/$SAMPLE.vcf.gz" | LC_ALL=C sort > "benchmarks/repeatability/$r.$SAMPLE.tsv"
done
diff -u "benchmarks/repeatability/$RUN_A.$SAMPLE.tsv" "benchmarks/repeatability/$RUN_B.$SAMPLE.tsv"
```

Repeat for every relevant admitted sample and review downstream fields under a predeclared tolerance if numerical output is not byte-identical. This does not authorize post hoc tolerance selection.

### 4. Expected outputs

Frozen identity checks and scientific-field comparison tables/diffs; original run archives remain intact.

### 5. Verification and quality control

Identical scientific fields are expected only under the specified deterministic scope. gzip headers, report timestamps and run IDs can legitimately differ. Check random seed, scheduling/tool versions and numerical conditions if science fields differ.

### 6. Error handling and troubleshooting

Unexplained differences: stop promotion and diagnose source/input/seed/tool identity. Do not overwrite the first run or manually reconcile differences.

### 7. Completion criteria

Repeatability outcome and exact comparison scope/tolerances documented; unresolved scientific differences are not hidden.


<a name="procedure-b16"></a>

## B16 — Perform one held-out evaluation and archive the experiment

### 1. Procedure identification

**Module:** held-out guard, evaluator and reporting. **Workflow:** benchmarking final evaluation.

**Purpose and scope:** Evaluate the frozen experiment once and preserve all decision/attempt evidence.

### 2. Prerequisites and requirements

All development analysis completed; protocol/settings/source/truth/masks frozen; held-out outcomes unexposed; exact completed final run ID. The output-directory-based guard is not an access-control system.

### 3. Detailed operational procedure

**Step 1. Confirm the final evaluation basis**

```bash
cd "$BENCH"
./pfvarx run-verify --run-id "$BENCH_RUN_ID"
cat inputs/partitions/partition.json
sha256sum inputs/truth/truth.json inputs/truth/confident_regions.bed pixi.lock
test ! -e benchmarks/heldout/holdout_attempts.jsonl
```

**Step 2. Execute the held-out command once**

```bash
./pfvarx benchmark-evaluate --run-id "$BENCH_RUN_ID" --holdout --out benchmarks/heldout
cat benchmarks/heldout/benchmark_report.json
cat benchmarks/heldout/holdout_attempts.jsonl
```

**Step 3. Review and archive the supported final outputs**

Apply B06–B12 to the held-out files without tuning/repeating them. Report development and held-out metrics separately, actual scored/excluded samples, masks/matching rules, uncertainty, low-frequency scope and limitations. Archive truth/checksum/partition/attempt records, frozen run, evaluation JSON/TSVs, normalization intermediates and resource profiles outside pristine source. Use R14 only for actual workflow/dashboard datasets; the dashboard does not automatically create a comparative benchmark-accuracy plot.

### 4. Expected outputs

`benchmarks/heldout/benchmark_report.json`, accuracy/mixture/K TSVs, normalized intermediates and the one-attempt JSONL. A failure can still consume the held-out attempt; preserve that evidence.

### 5. Verification and quality control

Run/partition identity must match. Evaluate frozen protocol acceptance with uncertainty and independent review. Never describe agreement with declared truth as universal clinical accuracy.

### 6. Error handling and troubleshooting

Second attempt refused: expected protection; do not delete the attempt log or change output directory to bypass it. An exposed/failed held-out experiment requires a documented new experiment and appropriate unseen material, not relabelling old samples.

### 7. Completion criteria

One final evaluation and all exclusions/attempts archived; reviewer sign-off recorded; disposable verification work cleaned only after preserving required evidence.


# Appendix A — Runtime output map and acceptance checklist

The source distribution ships none of the following study data. Paths are relative to
the chosen analysis/benchmark workspace. Optional outputs exist only for the selected
mode, assay, executed tool and scientifically eligible method. Tool versions can add
auxiliary logs/index files; use `find runs/<RUN_ID> -type f | sort` for the literal
inventory of that run. Snakemake temporary intermediates may disappear after success.

```text
WORKSPACE/
├── config/
│   ├── config.yaml                         editable study configuration
│   ├── quality.yaml                        selected policy
│   ├── compute.yaml                        computational envelope
│   └── profiles/resources.yaml             only after benchmark-resources
├── inputs/
│   ├── reference/                          optional local resources; absolute originals may be elsewhere
│   ├── fastq/ OR alignments/                applicable entry route
│   ├── metadata/                           optional study metadata and catalogues
│   ├── targets/                            targeted assays only
│   ├── truth/                              benchmarking only
│   │   ├── truth.json
│   │   ├── confident_regions.bed
│   │   └── <PURE_SAMPLE>.truth.vcf
│   └── partitions/partition.json           locked benchmark membership
├── prepared/
│   ├── manifest.json
│   ├── reference/                          derived parasite/combined FASTA and indexes
│   └── snpeff/                              conditional reference-matched SnpEff derivatives
├── metadata/validation_report.json
├── runs/<RUN_ID>/
│   ├── frozen/                             copied configuration/software
│   ├── frozen_config.json
│   ├── frozen_plan.json
│   ├── snapshot.json
│   ├── status.json
│   ├── runtime_exclusions.json             when exclusions are recorded
│   ├── workflow.complete                   real successful workflow completion marker
│   ├── logs/
│   │   ├── events.jsonl
│   │   ├── snakemake.pass*.log
│   │   └── <RULE>/<SAMPLE_OR_JOB>.log
│   ├── benchmarks/<RULE>/<SAMPLE_OR_JOB>.tsv measured instrumented jobs
│   └── results/
│       ├── alignments/<SAMPLE>.bam[.bai]
│       ├── qc/
│       │   ├── centralized_quality_matrix.json
│       │   └── samples/<SAMPLE>/
│       │       ├── qc_decision.json
│       │       ├── <SAMPLE>.caller_qc.json  FO002+
│       │       ├── <SAMPLE>.<LANE>.fastp.json/.html  FASTQ only
│       │       ├── <SAMPLE>.flagstat.json
│       │       ├── <SAMPLE>.markdup.txt     applicable WGS processing
│       │       ├── core.mosdepth.*         WGS coverage
│       │       ├── targets.*               targeted coverage
│       │       ├── <SAMPLE>.umi_family_sizes.txt   UMI branch
│       │       └── <SAMPLE>.entry_contract.json   aligned entry
│       ├── variants/samples/               FO002+
│       │   ├── <SAMPLE>.vcf.gz[.tbi]
│       │   ├── <SAMPLE>.evidence.tsv
│       │   ├── <SAMPLE>.callable.bed
│       │   └── <SAMPLE>.focal_loss.bed
│       ├── cohort/                         FO003+
│       │   ├── cohort_status.json          cohort formed, or SKIPPED_NO_ELIGIBLE_SAMPLES
│       │   ├── cohort.vcf.gz
│       │   ├── cohort_samples.txt
│       │   ├── sites.tsv
│       │   ├── evidence_matrix.tsv.gz
│       │   └── missingness.tsv
│       ├── annotation/                     FO004
│       │   ├── cohort.ann.vcf.gz
│       │   └── annotation_qc.json
│       ├── eligibility/                    requested downstream work
│       │   ├── cohort_decisions.json
│       │   ├── sample_analysis_eligibility.tsv
│       │   └── sample_analysis_eligibility.json
│       ├── downstream/DA001–DA007/         selected/internal analyses only
│       │   ├── status_manifest.json
│       │   └── method-specific files listed in R11
│       ├── report/
│       │   ├── pfvarx_dashboard.html
│       │   ├── dashboard_data.json
│       │   ├── clonality_report.txt        one line per sample: monoclonal, or polyclonal with clones and %
│       │   └── clonality_report.tsv
│       └── work/                           transient tool/caller intermediates
├── results -> runs/<LATEST_RUN_ID>/results
└── benchmarks/                            evaluation output, outside individual run
    ├── development/
    └── heldout/
        ├── benchmark_report.json
        ├── per_sample_accuracy.tsv
        ├── mixture_sites.tsv
        ├── k_selection.tsv
        ├── <SAMPLE>.calls.norm.vcf.gz
        ├── <SAMPLE>.truth.norm.vcf.gz
        └── holdout_attempts.jsonl          held-out only
```

The `development/` evaluation contains the same applicable metrics/intermediates,
without a held-out attempt record. Names such as `<SAMPLE>` are actual IDs, not literal
filenames. The literal derived SnpEff database tree is `prepared/snpeff/data/<GENOME_ID>/`,
with build configuration/check records under `prepared/snpeff/`. It is created by
the annotation workflow when FO004 executes, not by ordinary reference preparation.

# Appendix B — Current parameters: inspect before use

User configuration selects mode/design/resources and analyses; it does not calibrate
scientific thresholds. Current examples/defaults below are recorded implementation
values, not universal recommendations or evidence of scientific validity.

| Owner | Current value / operational interpretation |
| --- | --- |
| `layer1.yaml` | fastp quality threshold 15, maximum unqualified percentage 40, N-base limit 5, minimum length 50; cut-right false (window 4, mean quality 20 when enabled); fgbio consensus filter: minimum reads 3, base quality 45, maximum base error 0.2 |
| `layer1.yaml` | low-complexity filter false because indiscriminate filtering can bias AT-rich sequence; UMI branch preserves tag extraction order |
| `layer2.yaml` | min MAPQ 20, base quality 13, usable depth 5 (MalariaGEN Pf7/Pf8); candidates need 2 reads; pair-HMM caps 200 reads with 60 bp flanks; these remain source-class-labelled |
| `layer2.yaml` | screening -log10 threshold 6, PASS minimum fragments 2 and minimum score 10; not a calibrated detection guarantee; the mixture model skips SNVs within 10 bp of a called indel |
| `layer5.yaml` | DA001 minimum informative loci 10, MOIRE only with at least 10 informative samples (fewer: caller mixture test, no COI count or Fws), burn-in 10000 and 1000 samples, loci thinned to 500 bp; DA004 minimum markers 200 for WGS and 20 for targeted panels (pairs below 200 are low precision, unlabelled); Dcifer level 0.05 |
| `layer5.yaml` | DA006 minimum monoclonal samples 4 and informative sites 10; sample/site missingness ceilings 0.2; UFBoot/SH-aLRT 1000 each |
| `layer5.yaml` | DA003/DA005 population minimum 3; DA003 bootstrap 200; DA002 minimum callable/group 3, reporting floor 0.05 and MIXED below 0.95; DA007 GAM after 10 samples; provisional settings require scientific review |
| `layer6.yaml` | detailed variant embedding cap 200000; full totals remain distinguished from embedded details |
| `compute.yaml` | example total cores 4, RAM 6500 MiB, jobs 1, temp disk 20000 MB, per-job runtime estimate/ceiling field 720 minutes; `tmp_dir` null (system temporary directory) |

Every value, its source and its reference is listed in `docs/thresholds.md`. Confirm the actual release/configuration values with:

```bash
cd "$ANALYSIS"
cat config/parameters/layer1.yaml config/parameters/layer2.yaml
cat config/parameters/layer5.yaml config/parameters/layer6.yaml
cat config/compute.yaml
```

Routine users do not alter scientific defaults to force acceptance. Maintainers change
parameters only through an independently justified, versioned and reviewed development
experiment. Do not relabel provisional values as calibrated merely because tests pass.

# Appendix C — Concise troubleshooting decision table

| Observed condition | Inspect | Required action / stop condition |
| --- | --- | --- |
| Must run from workspace root | `pwd`, pixi.toml, workflow/Snakefile | Change to the intended installed workspace; do not run from a submodule |
| CONFIG_ERROR / unknown field | printed problems, generated schema/templates | Fix the actual typo/assay combination; do not remove validation |
| Global reference/GFF/BED failure | validation `summary`, validation_report.json, source identities | Stop; follow the numbered steps in the summary; obtain mutually matching resources |
| `CORE_MASK_SUSPECT` warning | the core BED's column 4, interval count and coverage | Replace a gene/feature BED with the real core-genome mask; record the reason if a custom mask is intended |
| Mixed flow-cell types for one sample | sample sheet `flowcell_type` per lane | Register lanes from different flow-cell types as separate samples |
| FASTQ mismatch/truncated gzip | pair validation and original delivery checksums | Obtain correct original pairs; do not trim records to make counts match |
| PARTIAL validation | failed sample IDs and continue-independent policy | Explicitly review/record exclusions; fix global failures before launch |
| Missing BAM index / wrong sort / wrong SM | @HD/@SQ/@RG, index, registry | Use a documented derived repair or correct producer file; never fake provenance |
| Targeted entry state rejected | genuine TrimPrimers/consensus evidence | Use original FASTQ route when upstream state cannot be established |
| Resource infeasible/OOM/quota | ready JSON, free/df, job log | Correct actual resources/envelope; preserve failed attempt; do not alter scientific gates |
| Silent long caller stage | job/process CPU, memory/disk and log timestamps | Allow active work; no reliable percentage estimate; avoid duplicate launches |
| Method skipped | eligibility/status_manifest/runtime_gates | Report the actual reason and denominator; do not create a placeholder result |
| `cohort_status.json` is SKIPPED_NO_ELIGIBLE_SAMPLES | per-sample `qc_decision.json`, blocked-sample reasons | No sample passed QC into the cohort; analyses report NO_ELIGIBLE_SAMPLES. Review QC; never relax gates to force a cohort |
| Resistance blocked | catalogue provenance and matching reference | Obtain independent verified catalogue; no status-only bypass |
| Frozen integrity failed | run-verify/snapshot/source identity | Stop reuse; restore correct software or run a documented NEW experiment |
| NO_CALLS in scoring | sample QC, run/mode, missing VCF | Report exclusion/missing result; never call it perfect reference prediction |
| Held-out attempt refused | partition and holdout_attempts.jsonl | Preserve attempt; do not delete/change directory to peek again |
| Export/schema problem | report warnings/console, source JSON/filter scope | Correct rendering/data contract; never substitute fabricated science |

# Appendix D — Documentation basis and verification of this manual

Commands and output contracts were checked against the finalized 0.0.0 local CLI,
configuration models/templates, Snakemake rules, caller/benchmark/downstream source and
analysis/benchmarking/reporting/deployment manuals. The provided user CI transcript records the earlier 0.2.0 release. The 0.0.0
orchestration update requires its own reviewed commit/hosted CI run; no new hosted-CI
result is inferred from that earlier transcript.
This document itself does not claim a real study or new benchmark was executed.

Authoritative companion sources in the release are `docs/analysis/README.md`,
`docs/benchmarking/README.md`, `docs/reporting/README.md`, `docs/deployment.md`,
`docs/testing.md` and `docs/status.md`. For tool command semantics consult the primary
[Pixi lock documentation](https://pixi.prefix.dev/latest/workspace/lock_file/) and
[BCFtools manual](https://samtools.github.io/bcftools/bcftools). The software's included
lock controls installed versions; a current web manual is not permission to silently
upgrade tools. Unsupported technologies/proposed methods from the original request
have deliberately not been developed into SOPs.


# Appendix E — Select, qualify and replace the Layer 2 caller

## Choose the caller for a new run

The project configuration defaults to `caller.backend: native`, including older YAML files
that omit the caller block. FO001 stops after alignment and QC. FO002–FO004 require
a caller that implements `pfvarx.layer2/1`. For an existing study, keep its prior frozen
run and use NEW when upgrading the pipeline or changing the backend.

1. Bootstrap the profile required by the desired final output. Review the assay,
   reference, core mask, targets when applicable, and sample sheet.
2. For the native caller, keep the default. For a replacement, install the caller and
   its runtime in a separately managed environment. Write a reviewed executable
   adapter that translates between the tool and the contract below. Pin the tool,
   adapter and dependency versions and record their checksums in the study protocol.
3. Run `sha256sum /absolute/path/to/adapter`. Set the following in
   `config/config.yaml`, using the actual digest.

```yaml
caller:
  backend: external
  executable: /absolute/path/to/adapter
  sha256: "<64 lowercase hexadecimal characters>"
```

4. In the workspace root, run the following after registering inputs.

```bash
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline --mode NEW
```

5. Check the Layer 2 logs in `runs/RUN_ID/logs/pfvarx_call/`, the per-sample caller QC,
   and the dashboard. Confirm the expected sample set, decisions, callability and
   downstream eligibility. Run `./pfvarx run-verify --run-id RUN_ID` for frozen files.
   External executable hashes are checked at planning and invocation. Interpreter,
   tool binaries and resources referenced by the executable require independent
   version and checksum records.

## Implement the adapter contract

The executable receives separate command-line arguments. It runs once per
QC-eligible sample, with the workspace root as its working directory. Use the supplied
paths verbatim and write all five outputs before returning zero. A nonzero exit or
invalid output fails the caller job. The runner validates external output in a staging
directory before publishing it to the existing result paths.

| Input option | Contract |
| --- | --- |
| `--bam` | Layer 1 coordinate-sorted alignment and adjacent index. Apply the assay-specific evidence policy to the parasite contigs. |
| `--reference` | Indexed parasite FASTA matching the alignment contigs. |
| `--core-mask` | Study mask in zero-based, half-open BED coordinates. |
| `--context-unreliable` | Pipeline context exclusion BED. It may contain no intervals. |
| `--frozen-config` | JSON containing effective assay, reference, policy and parameter configuration. Native Layer 2 parameters apply to the native implementation. |
| `--sample` | Exact sample ID for VCF and QC. |

| Output option | Contract |
| --- | --- |
| `--out-vcf` | Single-sample, variant-only VCF with reference-matched explicit biallelic SNP/small-indel alleles. Declare GT, AD, DP, PST, WSAF, WSAFLO and WSAFHI. Encode the dominant allele as haploid GT, observed ref/alt support as AD, usable depth as DP and the decision as PST. Use PASS only for accepted records. An absent record does not establish a reference genotype. |
| `--out-evidence` | Tab-separated header in the order shown below. POS is one-based. Record measured counts, decision and allele-fraction estimates with uncertainty. |
| `--out-callable` | Three-column BED of callable regions, ordered and nonoverlapping within contigs. Callability needs independent usable-depth evidence. |
| `--out-focal` | Three-column BED of measured focal coverage loss, ordered and nonoverlapping within contigs. |
| `--out-qc` | JSON with schema `pfvarx.caller_qc/1`, exact sample, caller version, calibration state, core and callable bases, record and status counts, and a `mixture` object with a status. Use an explicit unavailable status when technical mixture estimation is unavailable. |

```text
chrom pos ref alt kind y n alt_f alt_r ref_f ref_r mean_alt_bq tract screen status wsaf wsaf_lo wsaf_hi lhid pairhmm_reads
```

The separators above represent tabs. Use the status values in
`src/pfvarx/constants.py`. Use finite fractions between zero and one and ordered
intervals. Existing benchmarking consumes numeric WSAF estimates and bounds in
the evidence table. If the replacement lacks the needed evidence, qualify a real
calibration and evidence component before using it for that comparison. A VCF-only
wrapper cannot supply valid callability or uncertainty by filling in constants.
Use the native writer and QC schema as implementation references.

The external route loads shared interval and read-filter utilities independently of
the native caller package. After validating the external route in a disposable
workspace, the native `src/pfvarx/caller/` package can be omitted from that workspace.
Keep `caller_adapter.py` and `genomic_evidence.py`. The distributed package retains
the native caller as the default. External tool dependencies are managed separately.

## Run a fair comparison

1. Before inspecting held-out results, record the candidate versions, target assay,
   thresholds, resource budget, acceptance metrics and tie-breaking rule. Freeze
   independent truth, confidence intervals and development/held-out assignments.
2. Use separate benchmark workspaces with the same reads, reference, masks, truth,
   partition, QC and compute budget. Use development data to select parameters.
   Verify that each tool actually processed the same evaluable samples and regions.
3. Run the supported evaluation procedure for each backend. Compare SNP and indel
   precision/recall, callable coverage, minor-allele recovery and WSAF calibration.
   Review wall time and peak memory. Stratify by depth, assay, clonality, AT-rich and
   low-complexity context, including homopolymers. Report unsupported endpoints as
   unavailable and keep denominators visible.
4. Apply the predefined criteria to the held-out evaluation. Retain or replace the
   native caller only after reviewing the evidence for the intended assay. Follow
   the separate benchmarking guide for the available scoring and partition commands.

Native calling uses context-stratified error fitting, probabilistic screening, local
pair-HMM indel rescoring, read-backed linkage and adaptive technical mixture estimation.
Full local haplotype assembly and validated assay-specific calibration remain open.
Technical mixture count describes DNA evidence. Parasite genomes remain haploid.

The implemented assay branches are WGS, targeted non-UMI and targeted UMI. Targeted
barcode data can enter a targeted branch only when its inputs and protocol satisfy
that branch's contract. Shallow WGS uses the WGS branch with applicable QC and
callability gates. Neither TBS nor sWGS has an independently validated preset here.
No established third-party caller is bundled or prequalified by the adapter interface.

# Appendix F — Configure inputs, validate and prepare

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

**Reference release label** is a drop-down: PlasmoDB 68, Pf3D7 - release unknown, Custom
reference - study version 1, or Reference source and release unknown. Choose PlasmoDB 68 only
if you verified that release. **Other release (type it)** opens a text box for any other exact
release, such as another PlasmoDB release or a study's own assembly. In `config/config.yaml`
the setting stays `reference.release_label`, any text; the template placeholder
`SET_REFERENCE_RELEASE` must be replaced before input directories are created.

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
[assay design and library records](#appendix-g-assay-design-library-records-and-conditional-analysis).

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

# Appendix G — Assay design, library records and conditional analysis

Paired-end short-read genomic DNA is fixed. One run setting selects the processing
branch: `assay_design` (`wgs`, `targeted_umi` or `targeted_no_umi`). Each design accepts
FASTQ and BAM/CRAM entry (`input_type`). Read depth and filenames alone cannot identify a
laboratory protocol, so what each library *is* is recorded per sample in the manifest and
in the library catalog, and the pipeline checks it against the design.

## Designs and their branches

| `assay_design` | Libraries | FASTQ entry | BAM/CRAM entry (by `processing_state`) | Counts |
| --- | --- | --- | --- | --- |
| `wgs` | Direct WGS (`amplification_route: DIRECT`) and selective WGS (`SWGA`) share this branch | fastp, bwa-mem2, optional Human DNA removal (`reference.host_removal`, wgs only), merge lanes, duplicate marking (Picard) | `aligned`: duplicates marked; `duplicates_marked`: used as is; empty: duplicates marked unless the header records it | reads |
| `targeted_no_umi` | Amplicon or panel libraries without UMIs (`AMPLICON`) | fastp, bwa-mem2, merge lanes, primer trimming (fgbio TrimPrimers) | `aligned`: primers trimmed; `primers_trimmed`: used as is | reads |
| `targeted_umi` | Amplicon libraries with UMIs (`AMPLICON`) | fastp, UMI extraction (FastqToBam), alignment, grouping, consensus, filtering, realignment, primer trimming | `aligned` (raw reads, UMI in RX): grouping, consensus, realignment, trimming; `consensus`: trimming only, consensus is **not** generated again; `consensus_primers_trimmed`: used as is | UMI consensus molecules |

Primer trimming needs primers: a panel with only probe definitions (capture) is not
trimmed. A documented state is compared with the file before anything runs: consensus
reads carry the fgbio `cD` tag and raw UMI reads the `RX` tag, so a state the reads
contradict is rejected. A file whose header records fgbio TrimPrimers cannot be
documented as untrimmed (trimming twice would remove template bases). A claimed trim or
duplicate marking that the header does not record gives a warning.

Allele counts and depths of `targeted_umi` runs count consensus molecules, not sequenced
reads. The unit is written as `##PfVarX_count_unit=` in every per-sample VCF, as
`count_unit` in the caller QC and QC decision, and in DA002's method record.

## The sample manifest

Required columns stay `sample_id`, `lane`, `library_id` and the input files (`fastq_1` and
`fastq_2`, or `bam`). Optional columns describe the library:

| Column | Values | Rule |
| --- | --- | --- |
| `technique` | `WGS`, `SELECTIVE_WGS`, `AMPLICON`, `TBS_UMI`, `TBS_NON_UMI`; empty, `null`, `NA` or `.` = unspecified | `WGS` and `SELECTIVE_WGS` only with `wgs` (`SELECTIVE_WGS` records the SWGA route); `AMPLICON` with either targeted design; `TBS_UMI` only with `targeted_umi`; `TBS_NON_UMI` only with `targeted_no_umi` |
| `amplification_route` | `DIRECT`, `SWGA`, `AMPLICON` | `DIRECT` or `SWGA` for `wgs`, `AMPLICON` for targeted designs. Recorded as provenance: `DIRECT` only when documented, `SWGA` when documented or implied by `SELECTIVE_WGS`, `AMPLICON` when documented or implied by the panel's primers |
| `panel_id`, `panel_version` | IDs in the library catalog | Required for targeted designs; one panel per run |
| `library_protocol_id` | ID in the library catalog | Required for `targeted_umi` FASTQ (UMI extraction); optional otherwise |
| `processing_state` | see the table above | BAM/CRAM entry only; required for targeted designs |
| `adapter_kit` | `TRUSEQ`, `NEXTERA`, `AUTO`, `CUSTOM`; empty = unspecified | The 3' adapters fastp trims (FASTQ entry). Empty uses the library protocol's `adapters`, else `layer1.fastp_adapter_kit_default` (`TRUSEQ`). `CUSTOM` needs a protocol with the sequences. `AUTO` uses fastp's detection, which on AT-rich reads can pick its poly-T `FlowCell1` sequence and trim real T runs; use it only when the kit is unknown |

The technique label, route, protocol and panel must agree with each other and with the
design. A conflict is reported at input validation as
`sample S1 (library LIB1): technique TBS_UMI requires assay_design targeted_umi, but the
run's assay_design is targeted_no_umi`; only the affected sample is excluded. The resolved
technique, route (with where it came from), panel, protocol, adapters (kit, sequences and
source), state, steps run on entry and count unit are frozen in the run plan and shown in the dashboard's input section.

Example (`prepared/samples.tsv`, tab-separated; BAM columns left empty):

```text
sample_id	lane	fastq_1	fastq_2	bam	library_id	technique	amplification_route	panel_id	panel_version	library_protocol_id	processing_state
PF001	L001	inputs/fastq/PF001_R1.fastq.gz	inputs/fastq/PF001_R2.fastq.gz		PF001_lib1	TBS_UMI	AMPLICON	MAD4HATTER	1	UMI8_DUAL
PF002	L001	inputs/fastq/PF002_R1.fastq.gz	inputs/fastq/PF002_R2.fastq.gz		PF002_lib1	null		MAD4HATTER	1	UMI8_DUAL
```

An imported consensus BAM for the same panel:

```text
sample_id	lane	bam	library_id	technique	panel_id	panel_version	processing_state
PF003	L001	inputs/alignments/PF003.consensus.bam	PF003_lib1	TBS_UMI	MAD4HATTER	1	consensus
```

## The library catalog

`assay.library_catalog` (required for targeted designs) names a YAML file. The template
is `config/templates/library_catalog.yaml`; `./pfvarx init` copies it into place for
targeted workspaces.

```yaml
schema_version: pfvarx.library_catalog/1
panels:
  - panel_id: MAD4HATTER
    panel_version: "1"
    targets_bed: inputs/targets/targets.bed     # targeted intervals: target QC, territory
    primers_tsv: inputs/targets/primers.tsv     # amplicon primers: trimmed
    probes_bed: null                            # capture probes: recorded only
library_protocols:
  - library_protocol_id: UMI8_DUAL
    technique: TBS_UMI                          # optional default label
    amplification_route: AMPLICON               # optional default route
    adapters:                                   # optional: 3' adapters fastp trims
      kit: CUSTOM                               # TRUSEQ | NEXTERA | AUTO | CUSTOM
      read1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  # CUSTOM only
      read2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  # CUSTOM only
    umi:
      read_structure: "8M+T 8M+T"               # UMI extraction
      filter_min_reads: 3                       # optional; omitted settings use layer1.yaml
```

UMI consensus settings a protocol does not give use the maintainer defaults in
`config/parameters/layer1.yaml`; each resolved value is recorded with its source
(`LIBRARY_PROTOCOL` or `PARAMETER_DEFAULT`).

Adapters are resolved the same way for every design: the manifest `adapter_kit`, else the
protocol's `adapters`, else `layer1.fastp_adapter_kit_default`. The two sources must agree,
and all libraries of a sample must use one kit. UMI libraries whose read structure ends in a
fixed segment keep adapter trimming off, as before.

## Downstream eligibility from the panel and callable evidence

Every analysis applies to every design; whether it runs is decided from what was
sequenced (the **territory**: the core genome for `wgs`, the panel's targets for
targeted designs) and from what is callable. Decisions record the territory facts.

| ID | Analysis | Panel and callable-evidence rule | Further run-time gates |
| --- | --- | --- | --- |
| DA001 | Clone inference (COI, clonality) | — | ≥ 10 informative loci per sample; MOIRE needs ≥ 10 informative samples, otherwise the caller's mixture test |
| DA002 | Resistance markers | The territory must contain a whole catalogue codon (`PANEL_LACKS_MARKERS`); each sample needs at least one whole catalogue codon in its callable regions (`NO_CALLABLE_MARKERS`) | Per-marker callability (`UNKNOWN`), WHO K13 flag needs ≥ 3 callable infections |
| DA003 | Diversity | — | Group sizes; Fws needs ≥ 10 samples |
| DA004 | Relatedness and IBD | — | 200 informative markers per WGS pair, 20 per targeted pair (low precision below 200) |
| DA005 | Population structure | A panel needs ≥ 10 separate target regions (`da005_min_independent_loci`, `TOO_FEW_INDEPENDENT_LOCI`) | Population groups; ≥ 3 informative sites |
| DA006 | Phylogenetics | A panel needs ≥ 10 separate target regions (`da006_min_independent_loci`) | ≥ 4 monoclonal samples, ≥ 10 parsimony-informative sites, missingness limits |
| DA007 | Spatial and temporal | Inherits DA002 | Coordinates, dates, country |

A panel's loci are its separate target regions: sites inside one amplicon are linked, so
DA005 and DA006 results from a panel carry a caveat that they describe the panel loci,
not the genome. The 10-region minimums are provisional (no published standard).

Outputs live under `runs/RUN_ID/results/downstream/DAxxx/`. Review each
`status_manifest.json`, the eligibility matrix and cohort decisions before using
scientific output.

## Select automatic execution

For a fresh targeted non-UMI workspace with the compatible dependencies:

```bash
python3 bootstrap.py --workspace ../targeted-study --final-output FO004 \
  --design targeted_no_umi --auto-analyses
cd ../targeted-study
./pfvarx configure --gui
```

Choose FO004 and automatic selection in the form, or edit the equivalent YAML:

```yaml
schema_version: pfvarx.config/4
final_output: FO004
assay_design: targeted_no_umi
assay:
  library_catalog: inputs/metadata/library_catalog.yaml
analyses:
  selection: auto
  requested: []
  scales: {}
  population_column: population
```

Automatic selection requests every analysis; eligibility then records which ones the
panel and evidence support. Supply a matching GFF to enable annotation and marker
outcomes; without one, Layer 4 is omitted and DA002 and DA007 are blocked with recorded
reasons. Keep `CONTINUE_INDEPENDENT` to let qualified independent analyses continue.

## Upgrade from configuration schema 3

Schema 3 kept the panel, UMI and technique settings in `config.yaml`. Run:

```bash
./pfvarx config-migrate           # report what would change
./pfvarx config-migrate --apply   # write it, keeping .before-config4 backups
```

| Schema 3 (`config.yaml`) | Schema 4 |
| --- | --- |
| `assay.design: WGS / TARGETED_NON_UMI / TARGETED_UMI` | `assay_design: wgs / targeted_no_umi / targeted_umi` |
| `assay.sequencing_type` (`SWGS` becomes `SELECTIVE_WGS`) | manifest `technique` |
| `assay.amplification_route` | manifest `amplification_route` |
| `assay.targets_bed`, `assay.primers_tsv` | catalog panel `LEGACY_PANEL` version `config3`; manifest `panel_id`, `panel_version` |
| `assay.read_structure` | catalog protocol `LEGACY_UMI_PROTOCOL`; manifest `library_protocol_id` |
| targeted BAM/CRAM import (had to be trimmed, and for UMI consensus) | manifest `processing_state: primers_trimmed` or `consensus_primers_trimmed` |

Only empty manifest cells are filled. Saving a schema-3 workspace in the configuration
GUI runs the same migration first. Until it is migrated, a schema-3 `config.yaml` is
refused with this instruction. Runs frozen under schema 3 still resume, report and verify:
the workflow reads both spellings. Older identifier changes (config /2: DA001–DA003
renumbered) are still translated at load time.

# Appendix H — Reproduce the dry-run checks

Run from the installed development source root with its locked environment.

```bash
pixi run python scripts/verify.py --dry-run-only --output "/absolute/path/PfVarX preflight"
```

Choose a fresh destination. The verifier creates synthetic inputs, validates
and prepares them, checks readiness, and runs Snakemake preflight for each
scenario. Inspect each `verification.json` and its dry-run log. Acceptance
requires exit zero and no run or result directory. Checkpoints determine the
later QC and analytical jobs during execution. For an individual prepared study,
run `./pfvarx dry-run` from that workspace.

On RESUME, top-level Snakemake logs advance to the next
`logs/snakemake.passN.log`. Method-level logs can still be reused on retry.
Archive them before resuming if their original contents are needed.

## Transfer raw data and run preflight

After saving, open **Input files** to upload reads or alignments, or download
them using an HTTPS link and a destination filename. A trusted SHA-256 is
optional. When provided, PfVarX checks it before saving. Otherwise, PfVarX
records the computed checksum and marks the download unverified. Put reference
files in `inputs/reference/`. Register your reads or alignments before running
**Workflow Readiness Checks → Run readiness checks**. The tab runs `validate`, `prepare`,
`verify-prepared`, and `ready` in that order and stops on failure.
