# Changelog

## 0.0.0 (unreleased)

The release version was reset to 0.0.0 after 0.4.2, as the baseline for the first reviewed
release. The sections below 0.0.0 are the earlier development numbering, kept as history.

### Clone counts not inflated by indel artefacts

- SNVs within 10 bp of a called indel are left out of the mixture model
  (`mixture_indel_exclusion_bp`). Misaligned reads beside true indels produced false
  low-fraction SNVs that made single-clone samples look polyclonal. On the synthetic 40×
  cohort, single-clone samples estimated as one clone rose from 3 of 6 to 6 of 6 (PF002
  78.7/21.3, PF004 85/15 and PF006 62.7/29.5/7.8 are now monoclonal); mixed samples and all
  calls are unchanged. Caller QC records `mixture.snvs_excluded_near_indels`.

### Error model robust to mixed infections

- In mixed infections, minor-clone variants among the "reference" positions inflated the
  per-sample error model (dispersion phi fell from thousands to the bound of 10), so real
  variants were screened out. The fit now leaves out positions with a non-reference allele
  on at least 2 reads and at least 5% of reads, with a truncated likelihood
  (`error_fit_exclude_fraction`). On synthetic 40× cohorts with known truth, minor-clone
  recall rose from 25 of 839 to 715 of 839, majority-clone recall in mixed samples recovered
  (PF008: 55 to 383 of 395), precision rose from 95.8% to 96.8%, and the mixture weights moved
  to the truth (70/30: 79.5/20.5 to 71/29; 50/30/20: 91.8/8.2 to 54.7/25.4/19.9). Results
  change for every sample; rerun earlier studies with this version.
- Mixture-constrained calling is implemented and off by default
  (`mixture_constrained_calling`, `mixture_constrained_min_log10_lr`); when on, records
  carry MCF, MCLR and MCP. The 500 bp mixture-site spacing is now a parameter.

### Faster caller with identical results

- Profiling a real 29× sample showed the mixture bootstrap (58%), read evidence (22%), indel
  rescoring (10%) and error-model fit dominated a 66-minute caller run. PfVarX-Call now runs
  read evidence in 250 kb pieces, indel rescoring per contig, and bootstrap replicates and EM
  restarts in worker processes (`--threads`, set from the rule's threads), and evaluates the
  error-model likelihood once per distinct (errors, depth) pair. The same sample takes 9.7
  minutes with 24 workers, and its evidence table, BED files, caller QC and VCF are
  byte-identical to the serial run. Bootstrap randomness is drawn in the original order, and
  the stream is rewound where the serial stopping rule would stop.
- The caller log reports the time of each stage.
- Fixed while testing: parallel indel rescoring matched candidates back by allele, so two
  deletions sharing an anchor base (TTA>T and TTATA>T) collided; candidates are now matched by
  index, with a regression test.

### Shipped resistance catalogue repaired

- `config/registry/resistance_catalog_pfvarx_v1.tsv` had 12 fields per marker row under an
  11-column header, so every value after `gene_name` sat one or more columns off (the codon
  under `codon_start`, the wild-type amino acid under `alt_aa`). The 26 rows are rewritten into
  their declared columns with unchanged values; the status stays DRAFT_REQUIRES_VERIFICATION.
- Catalogue validation now fails a row whose field count differs from the header, a codon that
  is not a whole number, or an amino acid that is not a three-letter code. A test keeps the
  shipped catalogue well formed.

### Per-sample clonality report

- Every run with variant calls writes `results/report/clonality_report.txt` and `.tsv`, and
  the dashboard overview shows a clonality table: each sample is monoclonal, or polyclonal
  with the number of clones and each clone's percentage of the parasite DNA, from the
  PfVarX-Call mixture model, combined with DA001's call when DA001 ran. Disagreements,
  unresolved minor clones and the largest tested clone count are stated explicitly.

### SnpEff database ID safe for any release label

- The SnpEff database ID derived from `reference.release_label` now keeps only letters,
  digits, dots and underscores (`Pf3D7 - release unknown` becomes
  `pfvarx_Pf3D7_release_unknown`). A label with spaces made SnpEff fail to read its own
  configuration, so annotation and every FO004 analysis failed. A configured
  `reference.snpeff_genome_id` must use letters, digits, dots, underscores and hyphens;
  the GUI help says so. The annotated integration test now uses a label with spaces.

### Inputs recognised by content, not name

- Every input file is identified from its bytes and structure (FASTA, FASTQ, BAM, CRAM,
  alignment index, BED, GFF3, VCF, sample sheet, primer table, marker catalogue, JSON) by
  the new `src/pfvarx/inputs/content.py`. Validation, GUI upload and download, `add-input`,
  `fetch-input` and registration by link all use it; file extensions no longer decide
  anything, and `src/pfvarx/inputs/formats.py` is removed. A file with the wrong content is
  refused with `CONTENT_MISMATCH`, naming what it contains and what was needed.
- Two name rules stay because tools need them: a gzip-compressed FASTQ must end in `.gz`
  (fastp otherwise reads zero records without an error), and reference, BED, GFF, primer,
  sample-sheet and catalogue files must be uncompressed.
- A configuration setting that names only a folder now uses the one file in it whose content
  suits the setting; with none or several, validation lists the candidates.
- New structure checks: GFF3 (nine columns, contigs in the FASTA, coordinates inside them,
  gene-type and CDS features; PlasmoDB gene types are accepted) and the marker catalogue
  (required columns, marker rows, status line). Text inputs are decompressed by content.
- The GUI Input files tab shows the accepted content for each folder instead of filtering
  the file picker by extension.

### Documentation brought up to date

- README, HTML quick start, analysis, orchestration, testing, reporting, status, sequencing
  and SOP manual (Markdown and HTML) describe the current behaviour: validation summaries and
  the `CORE_MASK_SUSPECT` warning, the empty-cohort outcome and `cohort_status.json`, the
  design-specific DA004 marker rule, the flow-cell rule, `compute.tmp_dir`, report and package
  provenance, current parameter values (usable depth 5, DA004 200/20, MOIRE settings), and
  platform support (Linux x86-64; WSL2 on Windows).
- The orchestration guide documents every CLI command, with a table pairing each GUI tab with
  its command; `tests/unit/test_cli_gui_parity.py` now fails if a command is undocumented.
- The testing guide lists what the automated suite covers; `docs/thresholds.md` is linked
  from the README.
- Corrected stale statements: the default compute budget is 6,500 MiB (not 8,000 MB), and
  DA001 runs for FO004 only when requested or needed by another analysis.

### CLI and GUI parity

- The CLI and the configuration GUI now change together: `tests/unit/test_cli_gui_parity.py`
  fails when a setting lacks GUI help, a template entry or a configuration-guide row, when the
  guide text differs from the GUI help, or when a GUI tab or stage has no CLI command. The rule
  is written in `docs/architecture.md` and `CLAUDE.md`.
- `gui/field-guide.json` is the single source of setting help; `scripts/configuration_guide.py`
  regenerates the setting tables of `docs/configuration guide.md` and `.html` from it.
- Fixed help that disagreed between the two: DA001 is not forced for FO004 (it runs internally
  when another analysis needs it), and the sample-sheet path is not read-only in the GUI.
- The project templates now list `assay.sequencing_type`, `caller.executable` and
  `caller.sha256`, which only the GUI showed before.

### Validation explains failures

- Every failed or warned input check now says which input it is, what the pipeline uses it
  for, why it failed and numbered steps that fix it. `./pfvarx validate` and the GUI print
  this as a `summary`, and `metadata/validation_report.json` stores it per finding. A
  missing file names its expected path, notes a folder-only setting in `config/config.yaml`
  and lists files of the same type already present.
- New `CORE_MASK_SUSPECT` warning: a core-genome mask that looks like a gene list (gene
  names, more than 1 000 or overlapping intervals, or coverage outside 70-97% of the
  reference) is flagged, with instructions to build the standard core mask.

### Provenance, parameters and portability audit fixes

- Reports keep the pipeline version recorded when the run was made; the renderer version is
  listed separately under `report_renderer`. Run snapshots record the exact version and
  build of every scientific package (fgbio, Picard, SnpEff, IQ-TREE, hmmIBD, MOIRE, Dcifer,
  mgcv and others) from the environment's package records.
- Constants that were written into code are now parameters with a source class, at the same
  values: fastp sliding-window trim, caller candidate admission, pair-HMM read cap and flank,
  focal callability length, DA003 bootstrap replicates, Dcifer level and grid, and the DA007
  GAM sample minimum, basis size and interval level.
- New optional `compute.tmp_dir` sends Snakemake, Picard and fgbio temporary files to a chosen
  directory such as node-local scratch. The ZipperBams heap follows the job's memory instead
  of a fixed 2 GB.
- Validation rejects a WGS sample whose lanes mix patterned and unpatterned flow cells, since
  duplicates are marked once per sample with one optical distance.
- Workflow: cohort evidence declares its BAMs and callability files as inputs; target
  coverage has threads, memory and a benchmark; shell commands no longer read globals.
- `bootstrap.py` runs workspace initialisation with `--locked`; the container builder takes
  `--base` to pin the base image by digest and records whether it was pinned.
- MOIRE COI summaries are written to two decimals. Platform support (Linux x86-64; WSL2 on
  Windows) is documented, and confirmed MalariaGEN Pf8 citations are added.

### Version, cohort and analysis audit fixes

- One pipeline version: `__version__` in `src/pfvarx/__init__.py` is the only declaration.
  `pyproject.toml` reads it (`[tool.hatch.version]`), `pixi.toml` no longer declares one, and
  the archive name in CI and the SOP is derived. `docs/deployment.md` no longer states a
  conflicting 0.4.3 version. `tests/unit/test_version.py` checks every current-release
  reference.
- No sample passing QC no longer crashes the run: the cohort, annotation and analyses are
  skipped with `results/cohort/cohort_status.json` and `SKIPPED_INELIGIBLE` manifests
  (`NO_ELIGIBLE_SAMPLES`, `NO_COHORT`).
- DA004 marker rule depends on the assay: WGS needs 200 informative markers (Taylor 2019);
  targeted panels need `da004_min_informative_markers_targeted` (20, provisional), and pairs
  below 200 are `ESTIMATED_LOW_PRECISION` without a relatedness label. Dcifer pairs follow the
  same rule. DA004 reports `SKIPPED_INELIGIBLE` when no pair can be estimated.
- QC and analysis thresholds cite published standards where they exist (`docs/thresholds.md`).
- Fixture `--variant-scale` adds variant sites. New integration tests run all seven
  analyses on a WGS cohort, the targeted UMI route through DA002, and a run where every
  sample fails QC.

### GUI and CLI parity

- One rule set for placing input files: GUI upload/download and the CLI now share it. New
  `add-input` copies files from disk like a GUI upload; `fetch-input` gained the GUI's checks
  (folder must belong to the saved configuration and exist; extensions follow the input type).
- The GUI gained a **Workflow Check and Run** tab; `workflow-check` and `workflow-run` are the
  matching commands. Each GUI tab names its command.
- `init`/`input-directories` refuse the placeholder reference release label, as the GUI does.
- Sample sheets discovered inside the workspace record workspace-relative paths.

### End-to-end execution audit

Found by running real pipelines on the synthetic fixture for every design and entry route.

- UMI route: accept a fully fixed read structure such as `8M143T 8M143T` (validation rejected
  every one), pass it to fgbio as `8M+T` because adapter trimming shortens reads, and remove
  pairs too short to hold the UMI before fgbio FastqToBam, which stopped on the first such read.
- Annotation: write `cohort.ann.vcf.gz` as BGZF with a `.tbi` index instead of plain gzip,
  verify the GFF checksum as well as the FASTA, and use the protozoan mitochondrial codon
  table for the mitochondrial contig.
- Portability: the SnpEff configuration and the prepared-reference manifest no longer store
  absolute paths, so a workspace can be moved without re-indexing; DA007 reads country
  boundaries from the workspace rather than from the source tree.
- RESUME: retry samples that failed at run time; a finished run was reported as up to date and
  the failed sample stayed excluded.
- Input validation: an empty FASTQ read no longer marks the whole file malformed.
- Tests: two real end-to-end runs (annotated WGS cohort, UMI with short reads) join the dry runs.
- Lint: restore import order in `execution/preflight.py`, which had been failing CI.

### CPU deployment budgets

- Default to one job and a 6500 MiB scheduling budget.
- Check Linux CPU affinity, cgroup v2 and SLURM per-node limits during planning.
- Document PC, server, VM and single-node SLURM operation.

## 0.4.2 — Dry-run and path audit

- Quote paths in Snakemake shell commands and preserve list-valued arguments.
- Retain distinct top-level Snakemake pass logs across resumed invocations.
- Add an isolated dry-run-only verifier and refresh stale method messages.

## 0.4.1 — Audit corrections

- Isolate failed downstream methods under CONTINUE_INDEPENDENT and retry recorded failures on RESUME.
- Gate DA001 estimates on informative loci and suppress failed result payloads in the dashboard.
- Check external VCF depth and allele-fraction evidence and preserve defined zero F1 scores.
- Unify the SOP downstream matrix and document the audit evidence.

## 0.4.0 — Sequencing-aware selection and input configuration

- Added a loopback-only configuration editor using the CLI schemas.
- Added declared sequencing protocols, automatic method selection, and explicit exclusion reporting.
- Promoted COI to DA001 and translated legacy project identifiers on load.
- Updated the SOP, dashboard and caller integration guidance.

## 0.3.1 — Replaceable Layer 2 caller

- Added native and checksum-pinned external backend selection, a versioned output contract,
  and external artifact validation before publication.
- Moved interval and alignment conventions into a shared module and documented caller qualification.

## 0.3.0 — Configuration-driven orchestration

- Added explicit post-clone bootstrap with selective source export and locked
  assay/input/analysis-specific environments; selected methods include prerequisites.
- Added profile-aware workspace initialization, reviewed FASTQ discovery, zero-copy
  TSV registration and checksum-required/no-overwrite HTTPS resource acquisition.
- Added preparation artifact seals and mandatory derived-hash/FAI/dictionary checks;
  BAM/CRAM entry omits unused BWA indexing.
- Added ordered orchestrate boundaries and explicit verify-prepared command.
- Added keep/successful intermediate retention, preview-first safe cleanup and workspace
  locking; preserved scientific algorithms, original inputs and final output paths.
- Added optional Apptainer context/image builder; image execution remains unqualified.
- Updated CLI, manuals/SOPs, CI packaging and source archive assets.

## Software-only distribution — 2026-11-06

This historical 0.2.0 distribution retained its existing runtime version and tool lock.

- Removed study inputs, simulated masks, pilot launchers, generated analyses, results,
  historical evidence files, embedded source archives, installed environments and caches
  from the distribution. Existing working projects are preserved outside it.
- Consolidated duplicate purpose-specific compute/quality templates into one copy each.
  Purpose configurations are `config/templates/analysis.yaml` and `benchmarking.yaml`.
  `config-init` still creates the same three active user YAML filenames.
- Added separate analysis and benchmarking operating readmes. Consolidated implementation
  status/limitations and replaced documentation that referenced removed artifacts.
- Retained all implemented control, WGS/targeted/UMI/alignment-entry, calling, cohort,
  annotation, downstream, reporting and benchmarking components.
- Retained type-check fixes and explicit handling of absent optional BAM fields.
- Added deterministic allowlist-based packaging, archive-exclusion/roundtrip tests,
  isolated six-scenario verification and separate development/held-out evaluation outputs.
- Removed unverified reference-release defaults; GFF defaults to null for FO002 and is
  explicitly supplied for annotation workflows; FO004 now rejects missing GFF during
  configuration validation. Updated documentation and verification.

## 0.2.0

- Seven-layer typed control plane and Snakemake workflow with per-sample QC and downstream
  eligibility checkpoints, frozen run provenance, NEW/RESUME/EXTEND and offline reporting.
- Assay-specific WGS, targeted non-UMI and UMI read processing; existing-alignment entry.
- PfVarX-Call beta-binomial screening, local pair-HMM indel rescoring, linkage blocks,
  technical clone-mixture estimation, WSAF intervals and callability/evidence output.
- Variant-only cohort integration, reference-matched SnpEff and implemented DA001–DA007.
- Benchmark partitions, exact-normalized truth comparison, mixture evaluation, clustered
  uncertainty and measured resource profiles. Scientific calibration remains outstanding.

## Reporting dashboard refresh

- Modular offline frontend inlined into the same HTML output, with exclusive persistent
  visualization selectors, actual-value tooltips, filters and paginated tables.
- Added local JSON validation/import, synthetic labelling, measured resource views,
  input facts, full variant summary counts, actual optional VCF/PCA/tree fields and redacted logs.
- Verified CSV/XLSX/PNG/PDF export controls; documented available-data limits and renderer provenance.
- Reporting paths, scientific algorithms, workflow, input configuration and upstream output contracts unchanged.

## Production source finalization

- Removed generated synthetic report JSON and historical verification documentation;
  retained fixture generation and all dashboard source.
- Restricted packaged documentation to Markdown, added deployment/dependency guidance,
  and added CI source-archive validation. Scientific algorithms and output contracts unchanged.

