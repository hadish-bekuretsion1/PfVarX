# Changelog

## 0.0.0 (unreleased)

The release version was reset to 0.0.0 after 0.4.2, as the baseline for the first reviewed
release. The sections below 0.0.0 are the earlier development numbering, kept as history.

### Greedy job scheduling

- Snakemake runs with `--scheduler greedy`. The default ILP optimiser added nothing for
  PfVarX's jobs and could pause 10 seconds per scheduling decision before falling back
  ("Failed to solve scheduling problem with ILP solver in time"). Results are unchanged.

### Output files open in place

- In "Where the results are saved", nothing is downloaded any more: each file has **Open
  folder** (its folder in the browser) and **Copy path** (full path for IGV, a spreadsheet or
  the file manager). File names are links only for files the browser displays (HTML, JSON,
  TXT, logs, images, PDF); BAM, VCF, TSV and BED are no longer links.

### Output files in the dashboard

- Overview starts with **Where the results are saved**: the run folder's full path (Copy path,
  Open results folder), a button for each output folder, logs, benchmarks, frozen settings and
  software versions, and every final output file grouped by folder with its size, what it is
  and a link that opens it, with search. Temporary work files are left out. The list is also
  in `dashboard_data.json` (`outputs`).

### No unused files

- Tool byproducts that no step reads and that add nothing are no longer kept: mosdepth's
  whole-genome distribution and per-region depth tables for the core genome
  (`core.mosdepth.global.dist.txt`, `core.regions.bed.gz` and its index), the unused mosdepth
  files of target depth, and IQ-TREE's working files (`tree.ckp.gz`, `tree.bionj`,
  `tree.model.gz`, `tree.splits.nex`). The files the QC gates, analyses and report read, and
  the summaries, logs and inputs kept for checking and reproducing a run, are unchanged.

### CI uses half the GitHub Actions minutes

- CI runs on pull requests and on pushes to `main` and `develop`, not on every branch push,
  so each change is tested once instead of twice. A newer push to the same pull request
  cancels the older run. It can also be started by hand (Actions → Run workflow).

### One page per downstream analysis

- Downstream Analyses no longer uses an analysis drop-down: the seven analyses are tabs
  across the top (Infection Complexity, Resistance Markers, Genetic Diversity, Genetic
  Relatedness, Population Structure, Phylogenetic Analysis, Spatiotemporal Analysis). Each opens a page
  titled with the analysis name, with its own summary tiles and its charts and tables shown
  together in a grid under shared filters, instead of one chosen view at a time.
- New views: callability by marker (DA002), P(COI > 1) and classification share (DA001),
  nucleotide diversity (DA003), relatedness estimates with 95% intervals (DA004),
  explained variance (DA005).

### Time and RAM for every step

- Every workflow step now records its time and memory (Snakemake `benchmark` on all 29
  rules; before, only 5 did), in `runs/<RUN_ID>/benchmarks/<step>/`.
- The report's Resources section starts with **Resources by step: tool, time and RAM**: per
  step, its layer, the tools and versions used, jobs, total and longest time, CPU time,
  peak and mean RAM, average CPU cores and share of time, as a table or charts. Each
  measurement in `dashboard_data.json` carries `step`, `layer` and `tools`.

### Report dashboard redesign

- The report dashboard has a new look: a navy header with the teal PfVarX mark and
  the title "PfVarX Genomic Analysis Dashboard", centred, with the Menu button on the same
  row on phones and tablets; the Open report JSON, Reset dashboard and Export PDF buttons
  are removed (print from the browser for a PDF); run details, the synthetic-data label and
  calibration status appear once at the top of Overview; switching sections and analysis
  tabs is smoother because pages are prepared in the background, each section opens with one large
  heading such as "PfVarX · Inputs" (Downstream Analyses keeps its analysis titles), navigation cards with icons, the open section shown dark with a teal bar,
  KPI tiles with an accent bar, and cleaner tables, panels and status badges. Phones get
  the same cards through **Menu**.
- The sections follow the order a run is reviewed: Overview, Inputs, Quality Control,
  Variants, Downstream Analyses, Resources, Software, Logs. Contents, filters and exports
  are unchanged.

### Development guide, kept out of releases

- New `development/` folder with `PfVarX_Development_Guide.html`: every tool, package,
  pipeline step, algorithm, parameter and QC gate with its path and purpose; how to test
  parameter values and QC thresholds and choose the best; how to measure, replace or upgrade
  a tool; the testing ladder; documentation rules; and the Git and GitHub commands for this
  project on Ubuntu or WSL2.
- `scripts/package.py` refuses to package anything under `development/`; the release
  archive contains the analysis software only (`tests/unit/test_distribution.py`).

### Tools for this configuration

- Under **Install required tools**, the configuration editor lists every tool and package the
  form's current choices need: name, version locked in `pixi.lock`, what PfVarX uses it for
  and whether it is installed. It follows each change before saving (endpoint `tools-plan`).
  `./pfvarx install-tools --check` includes the same list (`tools`) for the saved
  configuration.

### Reference release label as a drop-down

- In the configuration editor the release setting is now **Reference release label**, a
  drop-down of PlasmoDB 68, Pf3D7 - release unknown, Custom reference - study version 1 and
  Reference source and release unknown, plus **Other release (type it)** for any other exact
  release. The YAML setting `reference.release_label` is unchanged and still takes any text.
- Form messages name the field by its label instead of its internal ID.

### Getting started guide

- New `docs/getting-started.md` (and `.html`): installing Pixi, getting PfVarX, creating a
  workspace (tools on demand or up front), configuring in the GUI or the YAML files,
  installing the required tools, adding input files, the steps from validation to the run,
  reviewing the report, updating, and a troubleshooting table. Linked from the README and
  the HTML quick start.

### Optional Human DNA removal for whole-genome sequencing (`reference.host_removal: true | false`)

- New true/false setting `reference.host_removal` (GUI: **Human DNA removal**, default
  false) with `reference.host_fasta`. When true, FASTQ read pairs are aligned once to the
  parasite plus human genome (`prepared/reference/combined.fasta`, human contigs prefixed
  `host_`), and pairs whose primary alignment is on the human genome are removed before
  duplicate marking. The BAM header keeps the parasite contigs only; later steps use
  `parasite.fasta` (with a new `parasite.dict`).
- Eligible only for `assay_design: wgs`. `targeted_no_umi` and `targeted_umi` are
  ineligible: `host_removal: true` with a targeted design is refused when the
  configuration is loaded, in the CLI and the GUI.
- GUI: for a targeted design the Human DNA removal checkbox is switched off and greyed out
  with the reason shown; the Human genome FASTA field appears only when removal is on.
- QC records `host_templates_removed`, `templates_before_host_removal` and `host_fraction`
  (descriptive gate `G2.host_fraction`; QC policy 1.2.0). `G2.mapped` is measured after removal.
- The human-sized bwa-mem2 index is supplied, not built: `prepare` stops with
  `HOST_INDEX_REQUIRED` and the steps, then checks a supplied index against
  `combined.fasta` (`.ann` against `.fai`) and seals it. `ready` adds
  `layer1.bwa_mem2_host_extra_mem_mb` (12 000 MB, provisional) to the bwa-mem2 request.
- Files of 1 GiB or more under `prepared/reference/` are sealed by size and modification time.

### SOP manual kept in step with the guides

- SOP appendices E, F and G copy `docs/caller.md`, `docs/configuration.md` and
  `docs/sequencing.md`. They had fallen behind (the adapter kits, Human DNA removal and most of
  the current configuration chapter were missing) and are now regenerated, in `.md` and
  `.html`, by `scripts/sync_sop.py`; `tests/unit/test_sop_mirrors.py` fails when they differ.
- Install examples across the README, guides, quick start and SOP use `--final-output`.

### Install required tools on demand

- `bootstrap.py` without `--final-output` or `--profile` now creates an on-demand workspace:
  all pipeline code, but only the small `control` environment (Python, pydantic, typer,
  pyyaml, pysam, htslib; no Snakemake, samtools, Picard or Java). Snakemake, samtools, Picard
  and mosdepth moved to the pixi feature `pipeline`, part of every tools environment.
- **Install required tools** on the Input directories tab, and `./pfvarx install-tools
  [--check]`, install the locked environment the saved configuration needs (final output,
  design, input type, analyses with prerequisites; the same rules as `--final-output`) and
  record it; later steps run in it. Inputs Preparation, readiness, dry run and pipeline stop
  with `TOOLS_NOT_INSTALLED` until it is installed; after a configuration change it installs
  only the newly needed environment. Up-front workspaces (`--final-output`, `--profile`) work
  as before.

### Installation follows the final output

- `bootstrap.py --final-output FO001|FO002|FO003|FO004` installs only the tools that final
  output uses and records it as the installation's maximum final output. FO001 (alignment
  and QC) gets the new `alignment-qc` environments: no variant caller code, bcftools, numpy
  or scipy (now the pixi feature `calling`, added to every FO002+ environment). The earlier
  `--profile` options still work; the two must agree when both are given.
- Picard comes from `picard-slim` (the same 3.5.0, without its R dependency; PfVarX runs only
  MarkDuplicates and CreateSequenceDictionary) and fgbio from `fgbio-minimal` (the same 4.1.1,
  without R, ggplot2 and scales; PfVarX runs no plotting fgbio tool), so R is installed only for
  the analyses that use it. An FO001 environment is about 1.0 GB instead of 1.9 GB.
- `scripts/check_installs.py` compares each kind of installation with the tools its final
  output runs, from `pixi.lock`, and reports anything missing or extra.
- A configuration above the installed final output, design, input type or analyses is refused
  with the reason and the bootstrap option to use; the GUI greys out final outputs that need
  a new installation.

### Adapters are set per library, not detected

- fastp's adapter detection (`--detect_adapter_for_pe`) was on for every library. On AT-rich
  reads it chose fastp's built-in poly-T `FlowCell1` sequence: on the synthetic WGS fixture,
  which has no adapters, 19% of reads lost about 10 bases of real sequence at T runs.
- Each library now names its kit: the manifest column `adapter_kit` (`TRUSEQ`, `NEXTERA`,
  `AUTO`, `CUSTOM`) or the library protocol's `adapters` (with `read1`/`read2` for `CUSTOM`),
  else `layer1.fastp_adapter_kit_default` (`TRUSEQ`, replacing
  `layer1.fastp_detect_adapter_for_pe`). fastp gets the sequences explicitly; `AUTO` keeps
  the old detection for libraries whose kit is unknown. Conflicting sources fail input
  validation for that sample (`ADAPTER_UNRESOLVED`). The resolved kit, sequences and source
  are frozen in the plan, written to each sample's QC decision and shown in the dashboard.

### Inputs Preparation: sample metadata from the input files, optional Study Metadata

- The GUI tab **Reference Preparation** is now **Inputs Preparation**, with the terminal
  command `./pfvarx inputs-preparation [--metadata-table PATH]`. It runs `prepare-samples`,
  `prepare`, `verify-prepared` and `ready` in order. `prepare` and `prepare-reference` work as
  before.
- `prepare-samples` builds the sample metadata from the input files when it is missing (a
  blank template for BAM/CRAM input) and creates its folder. New settings **Study
  Metadata** (`study_metadata`, No/Yes, default No) and **Study Metadata folder**
  (`study_metadata_dir`, default `inputs/metadata`): with Yes the study table there is merged
  (empty cells filled, nothing overwritten; several tables stop and ask; none continues with
  a warning), with No preparation continues with the sample metadata only. The result is
  written to `metadata/inputs_preparation.json`.
- The sample metadata now defaults to `prepared/samples.tsv`, a fixed path that is no longer
  shown in the configuration form. Existing workspaces keep the `samples_tsv` path in their
  configuration.
- **Input directories** creates only folders for files the user supplies; `inputs/metadata/`
  only with Study Metadata: Yes or a targeted design's library catalog. Later folders are
  created on demand.
- GUI: **Human DNA removal** is a False/True drop-down with **Human genome reference** shown
  only for True, and **Study Metadata** a No/Yes drop-down with its folder shown only for Yes.
- Workflow Readiness Checks now starts with `prepare-samples`.

### Configuration editor: professional layout

- The six workflow steps are numbered tabs (1-6), always visible as an even grid of horizontal
  buttons (one row of six on wide screens, 3 x 2 on laptops, 2 x 3 on narrow windows); labels
  no longer break inside a tab and the bar no longer scrolls. The tabs are the only step
  navigation (no duplicate drop-down).
- Every settings group is a card with a header; inside it each setting is a form-table row:
  name on the left (with the required mark), input on the right, help and the example value
  underneath. Narrow windows stack the name above the input.
- A header band and one action area (backup option, Create Configuration button, result
  message). One navy colour scheme defined once in `style.css`: primary and outline buttons,
  a clear disabled state, focus outlines; colours are kept when the page is printed.

### Study metadata join and naming standard

- `./pfvarx join-metadata --table PATH [--apply]` (and **Join study metadata** in the GUI's
  Input files tab) fills empty study columns of the sample sheet from a study table, matched
  by `sample_id`, on every lane. Existing values are never overwritten (conflicts are
  listed), file and library columns are never taken, unmatched IDs are listed, coordinates
  and dates are checked, and a backup is kept.
- `./pfvarx validate` warns (`NAMING_CONVENTION`) when lanes, FASTQ or alignment names depart
  from the documented naming standard; files are still recognised by content.

### `inputs/` holds only supplied files

- Reference derivatives (combined and parasite FASTA, `.fai`, `.dict`, bwa-mem2 index) are
  written to `prepared/reference/` beside their manifest and the SnpEff database, as the
  documentation already stated; they used to go to `inputs/reference/prepared/`, so the
  documented `prepared/reference/parasite.fasta` commands failed. `./pfvarx prepare` moves
  an existing `inputs/reference/prepared/` once, without rebuilding the indexes; read-only
  checks only report that the move is needed.

### Assay design selects the branch; libraries are described per sample (config schema 4)

- `assay_design` (`wgs`, `targeted_umi`, `targeted_no_umi`) is the one run setting that
  selects the processing branch, for FASTQ and for BAM/CRAM entry. `assay.sequencing_type`,
  `assay.amplification_route`, `assay.read_structure`, `assay.targets_bed` and
  `assay.primers_tsv` are no longer run settings: paired-end short-read genomic DNA is fixed.
- New optional manifest columns: `technique` (`WGS`, `SELECTIVE_WGS`, `AMPLICON`,
  `TBS_UMI`, `TBS_NON_UMI`; empty or `null` = unspecified), `amplification_route`
  (`DIRECT`, `SWGA`, `AMPLICON`, provenance), `panel_id` + `panel_version`,
  `library_protocol_id` and, for imported alignments, `processing_state`. Labels, routes,
  protocols and panels must agree with each other and the design; conflicts fail input
  validation for the affected sample, naming the sample and library. `SELECTIVE_WGS`
  records `SWGA`; direct and sWGA WGS share the WGS branch.
- New library catalog (`assay.library_catalog`, `pfvarx.library_catalog/1`): panels (targets
  BED, primers or probes) and library protocols (technique, route, UMI read structure and
  consensus settings, falling back to `layer1.yaml`). Required for targeted designs; one
  panel per run. Template: `config/templates/library_catalog.yaml`, placed by `init`.
- Imported BAM/CRAM files run only the steps their documented state still needs: WGS
  `aligned` marks duplicates; `targeted_no_umi` `aligned` trims primers; `targeted_umi`
  `aligned` (raw reads with RX UMIs) groups, builds consensus, realigns and trims, while
  `consensus` only trims (consensus is not generated again). Before schema 4, targeted
  imports without a TrimPrimers record were rejected and UMI imports were taken as finished.
  The documented state is checked against the reads (fgbio `cD` consensus tags, `RX` UMI
  tags) and `@PG`; a contradiction is rejected. UMI BAM entry now gets the aligner index.
- Read counts and UMI consensus-molecule counts are distinguished: `count_unit` in caller QC,
  QC decisions and DA002, and `##PfVarX_count_unit=` in each per-sample VCF.
- Eligibility comes from the panel and callable evidence (rule version 1.1.0, registry
  2.1.0): DA002 needs a catalogue codon in the territory (`PANEL_LACKS_MARKERS`) and one
  callable per sample (`NO_CALLABLE_MARKERS`); DA005 and DA006 are open to panels with at
  least 10 separate target regions (`da005_min_independent_loci`,
  `da006_min_independent_loci`, provisional), with a caveat that panel loci are linked.
- Migration: `./pfvarx config-migrate [--apply]` (and saving in the GUI) moves a schema-3
  workspace to schema 4, writing the catalog and filling empty manifest cells, with
  `.before-config4` backups. Schema-3 configs are refused with that instruction; runs frozen
  under schema 3 still resume and report. 24 missing Pixi environments (targeted with
  structure or phylogeny) were added; the 91 existing environments are unchanged.

### Second end-to-end audit

Scientific corrections:

- DA002: the WHO PfKelch13 flag (`who_k13_alert`) needs at least 3 callable infections in a
  group (`da002_min_callable_per_group`, previously defined but unused); fewer give
  `INSUFFICIENT_CALLABLE`. One mutant in one sample was a 100% prevalence and raised the flag
  for a whole population.
- DA001: MOIRE is run only with at least 10 samples that have enough informative loci
  (`population_frequency_min_samples`). It estimates allele frequencies and COI jointly, and in
  small cohorts neither is identified: on the 9-sample synthetic truth cohort all 6 single-clone
  samples were called polyclonal (P(COI > 1) 0.99-1.00), and across five seeds their P ranged
  from 0.00 to 1.00. One real sample got a COI of 9. Below the minimum, the
  monoclonal/polyclonal call comes from the caller's per-sample mixture test (new `basis`
  column, `CALLER_MIXTURE_LRT`; 6/6 single-clone and 3/3 mixed samples correct on that cohort),
  and the COI count and Fws are withheld. The clonality report does not count such a call as
  separate DA001 evidence.
- DA001 and DA003: Fws is withheld below the same minimum. It uses allele frequencies estimated
  from the analysed samples, in which each sample is 1/n of its own reference; one sample gave
  Fws 0 by construction.
- DA004: Dcifer no longer assumes COI 1 for a sample without a MOIRE count; it uses Dcifer's
  naive estimate (`getCOI`, `lrank = 2`). Monoclonal samples stay at 1.
- Targeted designs: target depth for the QC gate counts only reads at the caller's mapping
  quality (`-Q`), as the WGS core depth already did. The on-target read fraction is unchanged.
- The caller admits no indel longer than the declared scope (`indel_length_bound_bp`, 50 bp).

Consistency and streamlining:

- Parameters that changed nothing were removed: `quality_encoding_offset` (Phred+33 is an
  input requirement), `bootstrap_seed` (`random_seed` is used), `bwa_mem2_scoring`,
  `k_bootstrap_replicates`, `haplotypes_per_window_max`, `merge_missing_to_ref`,
  `vep_cross_check_fraction` and `control_chart_min_runs`. `fastp_low_complexity_filter` and
  `fastp_dedup` now set the fastp options (both stay off).
- The preflight resource plan disagreed with the rules (caller 1 thread and 2 000 MB against
  24 threads; UMI jobs 2 threads against 4). It now matches, and a test compares the two.
- The caller requests 8 threads and 5 000 MB (was 24 and 4 000). Measured on a 29x WGS sample:
  579 s and 3.6 GB with 8 workers, 474 s and 3.8 GB with 24; the old request was below the
  measured peak, and 8 threads let three samples run at once on a 24-core host.
- Workflow helper functions moved to `workflow/rules/common.smk`; `cohort_merge` is a shell
  rule. Unused code was removed. `docs/configuration guide.*` is renamed
  `docs/configuration-guide.*`.
- A stray test output (`0.95`) committed at the repository root was removed.
- CI: timeout raised to 60 minutes (runs took up to 28 of 30), and the integration step is
  named for what it runs.

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

