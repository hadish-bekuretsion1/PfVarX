# Changelog

## 0.0.0 — End-to-end execution audit

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

## 0.0.0 — CPU deployment budgets

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

