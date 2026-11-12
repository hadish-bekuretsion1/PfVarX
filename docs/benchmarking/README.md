# Benchmarking

Benchmark the exact same caller and workflow used for analysis, in a **separate
workspace** initialized with `purpose: benchmarking`. This guide covers software
verification, real-truth input contracts, evaluation, interpretation and reproducibility.
Nothing in this distribution is a truth dataset, an analysis result, or a calibrated model.

## 1. Separate software verification from scientific benchmarking

| Activity | Inputs | What it establishes |
|---|---|---|
| Unit/scientific tests | Generated small numerical/test examples | Tested software behavior and mathematical checks |
| Synthetic end-to-end verification | `pfvarx fixture` | Tool integration, QC/eligibility, outputs and evaluation plumbing |
| Real-truth benchmarking | Independently established truth and confidence regions | Agreement with that truth in the tested assay/context |
| Calibration | Predeclared development experiments | Empirical settings with provenance; not automatically performed by the evaluator |
| Held-out evaluation | Locked unseen truth partition | Final performance estimate for the frozen experiment |

Do not describe synthetic scores as clinical/field accuracy. A successful pilot without
independent truth is not accuracy validation. Define the assay, variant scope, coverage
and mixture ranges, masks, truth provenance, exclusion rules and acceptance criteria
before inspecting held-out results. Baseline caller comparisons require separate,
reference-consistent runs; PfVarX does not automatically run alternative callers.

## 2. Bootstrap a separate benchmarking workspace

From the reviewed source checkout:

```bash
python3 bootstrap.py --profile benchmark --workspace ../pfvarx-benchmarking
cd ../pfvarx-benchmarking
```

This installs only caller/core dependencies by default. Add `--analyses DA002 DA001`
for the intended FO004 scientific methods and prerequisites, `--design` for the real
assay and `--input-type BAM|CRAM` for aligned entry. `--rtg` explicitly adds the optional
RTG cross-check tool; it is not required by the built-in evaluator. A sample/analysis
workspace is never reused as a benchmark copy. Bootstrap creates editable config,
benchmark truth/partition folders and a blank sample registry, not truth datasets.

Use `./pfvarx` and `./pfvarx tool TOOL` from this workspace to select the recorded
locked environment. Configure experiment ID, reference/mask/GFF, assays and modes before
registration. Use FO002 for caller-only scoring; selected FO004 analyses require a
matching GFF. Profile capabilities cannot be expanded in-place by editing installation
metadata; bootstrap a new matching workspace instead. See [orchestration](../orchestration.md).

## 3. Run a synthetic smoke experiment first

Use a newly initialized, empty benchmarking workspace. Do not run `fixture` over real
reads or reference resources; it writes synthetic data at the configured conventional
fixture paths. The CLI refuses an existing FASTQ fixture, but this is not a general
backup or overwrite-protection tool for all resource files.

For a WGS caller smoke experiment, set `reference.release_label` to
`SYNTHETIC_SOFTWARE_FIXTURE`, keep the template reference filenames and change
`reference.gff` to `inputs/reference/Pf3D7.gff` only if testing FO004. These data are
software fixtures, not a genome or mask for a real study.

```bash
./pfvarx fixture --design WGS --purpose benchmarking --n-mono 10 --n-poly 3 --low-depth-sample
./pfvarx benchmark-partition --holdout-fraction 0.25 --seed 20261106
./pfvarx check-layout
./pfvarx validate
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline
./pfvarx status
```

The deliberately low-depth sample should be blocked at its own QC gate; other samples
can continue. Do not change thresholds to force its inclusion. Record the completed
run ID. A small synthetic run tests behavior; its runtime is not a full-genome estimate.

For all assay/output scenarios, use a **development/full installation** with all tools
and scripts available, then run the isolated verification runner instead:

```bash
pixi run verify --output /absolute/path/to/new-verification
```

Default scenarios are `wgs targeted umi benchmark fo001 fo003`. A subset is supported in that full/development installation:

```bash
pixi run verify benchmark --output /absolute/path/to/new-benchmark-verification
```

The runner copies software into new scenario workspaces, creates fixtures there,
checks preparation idempotency and dry-run isolation, runs the pipeline, reads output
contracts, checks frozen integrity and requested terminal manifests, and evaluates
both development and held-out partitions. It also checks refusal of repeat held-out
evaluation. Output locations must be new; no data is generated in the pristine release.

## 4. Supply real truth resources

Use a fresh benchmarking workspace for real experiments. Populate the standard sample
sheet and read/reference resources described in the analysis guide, then supply:

```text
inputs/truth/
  truth.json
  confident_regions.bed
  <MONOCLONAL_SAMPLE_ID>.truth.vcf
```

These files are required at execution time and are deliberately absent from the source
release. Truth VCFs must use the identical assembly, alleles and contig identifiers as
the analysis reference. The evaluator expects uncompressed `.truth.vcf` filenames for
monoclonal samples. Confident regions are a zero-based, half-open BED on that assembly.
The built-in evaluator uses one shared confidence BED; assays needing different truth
regions per sample require carefully separated experiments or external evaluation.

`truth.json` supplies sample composition and clone variants. This is the structural
format, not usable biological truth; the example positions/alleles are illustrative:

```json
{
  "samples": {
    "PURE_A": {
      "clonality": "MONOCLONAL",
      "true_k": 1,
      "components": [{"clone": "cloneA", "proportion": 1.0}]
    },
    "PURE_B": {
      "clonality": "MONOCLONAL",
      "true_k": 1,
      "components": [{"clone": "cloneB", "proportion": 1.0}]
    },
    "MIX_AB": {
      "clonality": "POLYCLONAL",
      "true_k": 2,
      "components": [
        {"clone": "cloneA", "proportion": 0.7},
        {"clone": "cloneB", "proportion": 0.3}
      ]
    }
  },
  "clones": {
    "cloneA": [{"chrom": "Pf3D7_01_v3", "pos": 12345, "ref": "A", "alt": "G"}],
    "cloneB": []
  }
}
```

Every truth sample ID must match the registered study sample ID. For pure samples,
provide the corresponding VCF. For mixtures, use independently established clone
variants, actual DNA proportions summing to one, and all constituent clones. Variant
positions in JSON/VCF are one-based; BED coordinates are zero-based. An empty `clones`
object is structurally sufficient for an all-monoclonal experiment but provides no
mixture truth. Do not infer or fill unknown true K/composition from PfVarX's own output.
Technical mixture K and biological COI are different estimands.

The current truth/evaluation adapter is not a complete biological-truth validator.
Independently verify reference matching, VCF validity, confident-region bounds, truth
completeness, clone identifiers, sample membership and proportions before evaluation.
Store checksums and provenance for FASTA, BED, reads, VCFs, JSON and catalogues.

## 5. Lock development and held-out membership before tuning

```bash
./pfvarx benchmark-partition --holdout-fraction 0.25 --seed 20261106
```

This writes `inputs/partitions/partition.json` with development IDs (`benchmarking`),
held-out IDs, seed and content hash. A second lock is refused; a changed partition hash
is rejected during evaluation. Use a fraction strictly between zero and one and enough
samples for meaningful development and held-out estimates. The current split is random
by sample ID, not stratified or family/lineage-aware. Prevent biological leakage through
study design before using it; technical replicates or related isolates are not
independent merely because they have different IDs.

Keep held-out labels/results inaccessible during development as an experiment protocol;
the software does not provide an access-control boundary. Freeze preprocessing,
parameters, reference/masks, caller version and acceptance criteria before the held-out
step. If scientific settings are calibrated, maintain their parameter source/provenance
entries and code review; the evaluation command does not tune or relabel them for you.

## 6. Execute calls without truth and verify the run

```bash
./pfvarx validate
./pfvarx prepare --dry-run
./pfvarx prepare
./pfvarx verify-prepared
./pfvarx ready
./pfvarx dry-run
./pfvarx pipeline
./pfvarx run-verify --run-id YOUR_COMPLETED_RUN_ID
```

Review per-sample QC and failures before scoring. A QC-excluded sample is not silently
made callable. Record exclusion reasons and report eligible/analysed denominators.
Scientific methods can skip despite successful processing. Do not evaluate an incomplete
run or reinterpret missing VCFs as perfect reference calls.

## 7. Score development and held-out sets into separate directories

Always specify the completed run ID when multiple runs exist:

```bash
./pfvarx benchmark-evaluate --run-id YOUR_COMPLETED_RUN_ID --out benchmarks/development
```

Inspect and finish the development experiment. After freezing the protocol, perform
one final held-out evaluation:

```bash
./pfvarx benchmark-evaluate --run-id YOUR_COMPLETED_RUN_ID --holdout --out benchmarks/heldout
```

Separate output directories preserve both reports; the evaluator otherwise overwrites
filenames in its output directory. Keep `benchmarks/heldout/holdout_attempts.jsonl` with
the experiment. A second held-out attempt in that same output directory is refused,
including after a failed attempt. The guard is output-directory-based, not a global
experiment lock; changing `--out` is not permission to repeatedly inspect held-out data.
Do not delete or bypass the attempt record. Log a new experiment if a protocol genuinely
changes, and do not reuse exposed held-out samples as if they were unseen.

## 8. Understand the metrics and files

Each evaluation output directory contains:

| File | Contents |
|---|---|
| `benchmark_report.json` | Summary accuracy, uncertainty, mixture calibration and partition identity |
| `per_sample_accuracy.tsv` | SNV/indel TP, FP, FN, sensitivity, precision and FP per callable Mb |
| `mixture_sites.tsv` | Expected/observed WSAF, interval bounds, detection and coverage |
| `k_selection.tsv` | Declared true K versus technical caller K |
| `*.calls.norm.vcf.gz`, `*.truth.norm.vcf.gz` | Normalized comparison intermediates |
| `holdout_attempts.jsonl` | One-attempt record, held-out evaluation only |

The implemented matcher compares exact `(chrom, pos, ref, alt)` after `bcftools norm`
left alignment and multiallelic decomposition. Only PASS calls are scored for variant
accuracy, with variant anchors restricted to the declared confidence BED. Truth variants
outside the callability BED can still be counted as false negatives; document this
interpretation. FP-per-callable-Mb uses the caller's callable span, not automatically
the intersection with confidence regions. These denominator choices must be reported
when comparing external tools. Reference mismatch warnings from normalization must be
investigated; the evaluator is not a substitute for strict truth/reference validation.

Sensitivity is TP/(TP+FN), precision TP/(TP+FP), and F1 their harmonic mean. Uncertainty
uses a wild cluster bootstrap of the linearized ratio estimator over pure truth samples:
Rademacher weights with at least ten clusters, Webb weights below, enumeration when
small enough. Very small truth sets cannot establish general performance.

Mixture summaries include WSAF slope/intercept/RMSE/interval coverage, detection by
expected-frequency bin and technical K agreement. Check source tables and sample
exclusions; a summary alone does not establish mixture identifiability or clinical COI.

RTG Tools is installed only with `--rtg` or the full/development environment. The executable is detected and reported as available, but RTG
`vcfeval` is **not automatically executed**. For a pure haploid truth sample, a manual cross-check can use the following commands
(adapt the sample ID, paths, thread budget and experiment directory):

```bash
TRUTH_SAMPLE=PURE_A
mkdir -p benchmarks/rtg
./pfvarx tool rtg format -o benchmarks/rtg/reference.sdf prepared/reference/parasite.fasta
./pfvarx tool bcftools norm -f prepared/reference/parasite.fasta -m -any --check-ref e -Oz -o benchmarks/rtg/truth.vcf.gz "inputs/truth/${TRUTH_SAMPLE}.truth.vcf"
./pfvarx tool bcftools view -f PASS -Oz -o benchmarks/rtg/calls.vcf.gz "results/variants/samples/${TRUTH_SAMPLE}.vcf.gz"
./pfvarx tool bcftools index -t benchmarks/rtg/truth.vcf.gz
./pfvarx tool bcftools index -t benchmarks/rtg/calls.vcf.gz
./pfvarx tool rtg vcfeval --baseline benchmarks/rtg/truth.vcf.gz --calls benchmarks/rtg/calls.vcf.gz --template benchmarks/rtg/reference.sdf --evaluation-regions inputs/truth/confident_regions.bed --sample "$TRUTH_SAMPLE" --sample-ploidy 1 --vcf-score-field QUAL --at-score 0 --threads 2 --output benchmarks/rtg/evaluation
```

Use a new RTG output directory for each comparison. `QUAL` is explicit because PfVarX
does not supply the default RTG `GQ` field; ploidy one is appropriate for this pure
haploid comparison. The score threshold is predeclared to include the passed call set,
rather than selecting a best threshold on held-out truth. Archive commands, versions
and all settings. Review `rtg ... --help` for the installed release and design an
appropriate protocol separately for mixtures. hap.py and automatic baseline-caller
comparisons are not included. Exact normalized matching can disagree with haplotype-aware
scoring on complex indels; report that limitation explicitly.

## 9. Measure resources and preserve reproducibility

```bash
./pfvarx benchmark-resources --run-id YOUR_COMPLETED_RUN_ID
```

This summarizes available Snakemake benchmark records into
`config/profiles/resources.yaml`, generated only at runtime. It records dataset/hardware
basis; coverage is limited to instrumented rule classes. Review measurements before
using the profile on another assay, machine, or cohort. The planner may use it for
future estimates; verify basis compatibility yourself rather than assuming automatic
hardware compatibility enforcement.

Archive the frozen run, input/truth/mask checksums, partition and attempt records,
software/lockfile identity, command/configuration changes, all evaluation directories,
QC exclusions, resource measurements and the experiment protocol outside the pristine
software distribution. `run-verify` checks frozen software/configuration integrity;
external truth/read files need their own checksum manifest. Keep the routine analysis
copy separate and promote calibrated settings only after documented scientific review.

## 10. Operational additions in 0.3.0

Discovery/registration and explicit checksum-verified HTTPS fetching work identically
to analysis mode. No truth/reference release is downloaded implicitly.
`verify-prepared` must pass before execution, and pipeline enforces it independently.
`orchestrate` can execute the preparation/control chain after real inputs are configured.
Retention removes only completed dependency intermediates/safe run work files; benchmark
evaluation outputs, held-out attempts, truth, partitions, frozen records and scientific
products are retained. Keep them under the independent experiment archival protocol.

For native versus external comparisons, use the [caller qualification procedure](../caller.md)
and preserve the same truth, confidence regions, partition and resource budget across workspaces.
