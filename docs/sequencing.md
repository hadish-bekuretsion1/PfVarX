# Assay design, library records and conditional analysis

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
PF003	L001	inputs/bam/PF003.consensus.bam	PF003_lib1	TBS_UMI	MAD4HATTER	1	consensus
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
| DA001 | Clone inference (COI, clonality) | Multiallelic microhaplotype loci when the cohort has ≥ `microhap_min_loci_for_basis` of them (`analyses.microhap_panel`), else biallelic SNVs; the basis is recorded | ≥ 10 informative loci per sample; MOIRE needs ≥ 10 informative samples, otherwise the caller's mixture test |
| DA002 | Resistance markers | The territory must contain a whole catalogue codon (`PANEL_LACKS_MARKERS`); each sample needs at least one whole catalogue codon in its callable regions (`NO_CALLABLE_MARKERS`) | Per-marker callability (`UNKNOWN`), WHO K13 flag needs ≥ 3 callable infections |
| DA003 | Diversity | — | Group sizes; Fws needs ≥ 10 samples |
| DA004 | Relatedness and IBD | Dcifer uses the same microhaplotype loci as DA001 when available; hmmIBD always uses biallelic SNVs | 200 informative markers per WGS pair, 20 per targeted pair (low precision below 200). The marker thresholds were derived for biallelic markers (Taylor 2019) and are applied to the locus count, so a low-precision label on a microhaplotype panel is conservative |
| DA005 | Population structure | A panel needs ≥ 10 separate target regions (`da005_min_independent_loci`, `TOO_FEW_INDEPENDENT_LOCI`) | Population groups; ≥ 3 informative sites |
| DA006 | Phylogenetics | A panel needs ≥ 10 separate target regions (`da006_min_independent_loci`) | ≥ 4 monoclonal samples, ≥ 10 parsimony-informative sites, missingness limits |
| DA007 | Spatial and temporal | Inherits DA002 | Coordinates, dates, country |
| DA008 | Copy-number markers | Whole-genome sequencing only (amplicon depth does not follow copy number) | Median core depth ≥ 20; verified copy-number catalogue; GFF |
| DA009 | pfhrp2/pfhrp3 joint deletion | Every design, with a different evidence model each: WGS by depth ratio (as DA008), targeted panels by amplicon presence against control amplicons. A panel must carry a pfhrp2 or pfhrp3 target (`PANEL_LACKS_HRP_TARGETS`) and at least `da009_min_control_amplicons` other targets (`CONTROL_AMPLICONS_INSUFFICIENT`) | WGS: median core depth ≥ 20. Targeted: control-amplicon depth ≥ `da009_min_control_depth`, else every gene is UNKNOWN. Verified copy-number catalogue; GFF |

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
outcomes; without one, Layer 4 is omitted and DA002, DA007, DA008 and DA009 are blocked with recorded
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
