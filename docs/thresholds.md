# Thresholds and their scientific basis

Every numeric decision PfVarX makes is listed here with its three outcomes — accepted,
accepted with warning, rejected or filtered out — and where the value comes from.
The values live in `config/registry/qc_policy_pfvarx_qc_v1.yaml` (policy version 1.3.0) and
`config/parameters/layer*.yaml`; no threshold or setting is written in the program code (see
[Where values are kept](#where-values-are-kept)). Each entry there carries a source class:

| Source class | Meaning |
| --- | --- |
| `LITERATURE_SUPPORTED` | A published standard or the tool authors' stated recommendation. The reference is given below. |
| `SOFTWARE_DEFAULT` | The pinned tool's own default. |
| `FRAMEWORK_POLICY` | A design decision of this pipeline, not an empirical claim. |
| `UNCALIBRATED_PROVISIONAL` | No published standard exists for the metric. The value is a working value; outputs that depend on it are labelled uncalibrated. |
| `UNSET_BLOCKING` | No defensible value exists. The metric is reported and the gate is `NOT_EVALUATED`, never `PASS`. |

A literature value is a standard used by others on their data and their caller. It has not
been calibrated against truth data for PfVarX-Call. The [benchmarking
protocol](benchmarking/README.md) is still the route to calibrated values.

## How the three outcomes work

- **Accepted (PASS):** the sample continues to every stage.
- **Accepted with warning (WARNING):** the sample continues only into the stages the gate
  permits, and every affected result carries the warning.
- **Rejected (FAIL):** the sample stops at that gate and stays in the denominators as an
  exclusion. A gate marked "warning only" never rejects.

## Stage 0 — inputs that must be supplied

| Stage | Required information |
| --- | --- |
| All runs | Reference FASTA; core-genome BED for the same assembly; exact reference release label; sample sheet with `sample_id`, `lane`, `library_id`. |
| FASTQ entry | Paired `fastq_1`/`fastq_2`, Phred+33. For WGS, `flowcell_type` (`patterned` or `unpatterned`) sets the optical-duplicate distance. |
| BAM/CRAM entry | Coordinate-sorted, indexed, `@RG` with ID, SM, LB, PL; contigs matching the reference. Each targeted row documents its `processing_state`; consensus states need fgbio `cD` tags in the reads, raw UMI reads need `RX` tags. |
| Targeted designs | A library catalog with the panel (target BED; primer table, 1-based inclusive, or probe BED) named by `panel_id` + `panel_version` in every row; UMI FASTQ also needs a library protocol with the read structure. |
| Annotation and DA002 | GFF matching the FASTA; a marker catalogue with status `MAINTAINER_VERIFIED`. |
| DA003, DA005 by population | Population column in the sample sheet. |
| DA007 | `latitude`, `longitude`, `collection_date`, `country`. |
| DA008 | GFF matching the FASTA; a copy-number catalogue with status `MAINTAINER_VERIFIED`; whole-genome sequencing. |

The core genome should follow the accessible-genome definition of Miles et al. 2016,
which excludes subtelomeric and internal hypervariable regions.

## Stage 1 — read QC, alignment and coverage

| Gate | Metric | Accepted | Accepted with warning | Rejected | Source |
| --- | --- | --- | --- | --- | --- |
| G1.read_retention | Fraction of reads kept by fastp | ≥ 0.70 | 0.50–0.70 | < 0.50 | Provisional. No published standard. |
| G1.mean_length | Mean read length after trimming | ≥ 70 bp | < 70 bp | never (warning only) | Li 2013: BWA-MEM is designed for reads of 70 bp and longer. |
| G1.q30 | Fraction of bases ≥ Q30 | reported | — | — | Not gated: instrument specifications differ by platform and read length. |
| G2.mapped | Fraction of reads mapped | ≥ 0.80 | 0.50–0.80 | < 0.50 | Provisional. Host DNA content varies with sample type (Oyola 2016). With host-read removal, measured after host pairs are removed. |
| G2.host_fraction | Fraction of read pairs removed as host (human) | reported | — | — | Only with `reference.host_fasta`. Not gated: the depth and callability gates decide whether enough parasite reads remain. |
| G2.properly_paired | Fraction properly paired | reported | — | — | Unset. No published standard. |
| G2.duplicates (WGS) | Duplicate fraction | ≤ 0.30 | 0.30–0.60 | > 0.60 | Provisional. Selective whole-genome amplification raises duplication. |
| G2.depth (WGS) | Median depth over the core genome | ≥ 20 | 5–20 | < 5 | Rejection level: MalariaGEN Pf7 and Pf8. Warning level: framework policy. |
| G2.callable (WGS) | Fraction of the core genome with ≥ 5 reads | ≥ 0.80 | 0.50–0.80 | < 0.50 | Rejection level: MalariaGEN Pf7 and Pf8. Warning level: framework policy. |
| G3.on_target (targeted) | Fraction of mapped reads on target | ≥ 0.70 | 0.50–0.70 | < 0.50 | Provisional. Panel specific. |
| G3.dropout (targeted) | Fraction of amplicons below usable depth | ≤ 0.20 | 0.20–0.50 | > 0.50 | Provisional. Panel specific. |
| G3.amplicon_depth (targeted, no UMI) | Median depth per amplicon | ≥ 100 | < 100 | never (warning only) | Early 2019; MAD4HatTeR guidance. |
| G3.umi_singletons (UMI) | Fraction of single-read UMI families | reported | — | — | Unset. No published standard. |

## Stage 3 — cohort QC (FO003 and FO004)

Applied after the cohort evidence step to every sample in the cohort
(`results/cohort/cohort_qc.json`). A rejected sample keeps its own VCF and stays in the
cohort VCF, but is excluded from the downstream analyses (reason `COHORT_QC_FAIL`); a sample
accepted with warning continues into them (`COHORT_QC_WARNING` is recorded only when a
stage is not permitted).

| Gate | Metric | Accepted | Accepted with warning | Rejected | Source |
| --- | --- | --- | --- | --- | --- |
| G5.cohort_missingness | Fraction of cohort sites where the sample is MISSING_LOW_DEPTH or UNCALLABLE | ≤ 0.20 | 0.20–0.50 | > 0.50 | Provisional. No published standard; set from the phylogeny missingness limit (≤ 20 %) and the Pf7/Pf8 callability rejection level (50 %). |
| G5.cohort_mixed_fraction | Fraction of cohort sites with MIXED evidence | reported | — | — | Descriptive: mixed infection or artefact (for example paralog mismapping). |

The same file reports how many cohort sites are called (ALT, MIXED or REF) in at least
80 % of the samples (`layer3.cohort_site_call_rate_min`); that site completeness is reported,
not gated.

Basis for the literature-supported rows:

- **Callability and depth.** MalariaGEN Pf7 computed callability "with a minimum depth
  of 5" and removed samples "with callability <50%". A median core depth below 5 means
  fewer than half the positions reach 5 reads, so the depth and callability rejection
  levels are the same criterion seen two ways. Oyola 2016 used the same yardstick
  (">5× coverage for 50% of the genome") for dried blood spots.
- **Amplicon depth.** Early 2019 found reduced sensitivity and precision in samples with
  "low read count (< 100 reads per amplicon)" and 46% accuracy at "< 25 reads per
  amplicon". The MAD4HatTeR documentation recommends 100 reads or more per amplicon.
  The rule is stated for raw reads, so it is not applied to UMI consensus reads. MalariaGEN Pf8 (2025) applies the same rule: samples "with <50% of their genome callable at 5X" fail QC. The read-fraction rule (fewer than 2 reads or under 5%) is cited from Pf7, whose methods state it; the Pf8 main text defers those details to its supplementary methods.

Fixed read-processing settings:

| Setting | Value | Source |
| --- | --- | --- |
| fastp base quality counted as qualified | Q15 | fastp default (Chen 2018). |
| fastp unqualified bases allowed per read | 40% | fastp default. |
| fastp N bases allowed per read | 5 | fastp default. |
| fastp minimum read length | 50 bp | Provisional. |
| fastp 3′ sliding-window trim (`--cut_right`) | window 4, mean quality 20 | fastp defaults. |
| Optical duplicate distance | 100 unpatterned, 2500 patterned | Picard documentation. |
| UMI grouping | directed adjacency, 1 edit | Smith 2017; fgbio best-practice pipeline. |
| UMI consensus: minimum input base quality | 20 | fgbio best-practice pipeline. |
| UMI consensus filter: minimum reads per molecule | 3 | fgbio best-practice pipeline. |
| UMI consensus filter: minimum consensus base quality | 45 | fgbio best-practice pipeline. |
| UMI consensus filter: maximum base error rate | 0.2 | fgbio best-practice pipeline. |

The fgbio authors state that their filter values "aim for a balance between generating
high quality consensus reads and yield" and need tuning for the library preparation.
Molecules seen fewer than three times are filtered out; a library with small UMI
families will lose depth at this step.

## Stage 2 — variant calling

| Decision | Accepted | Filtered out | Source |
| --- | --- | --- | --- |
| Usable depth at a position | ≥ 5 reads | < 5 reads: `INSUFFICIENT_DEPTH`, position not callable | MalariaGEN Pf7 and Pf8. |
| Reads supporting an allele | ≥ 2 | < 2: allele treated as absent | MalariaGEN Pf7. |
| Mapping quality | ≥ 20 | < 20: read not used | Provisional. |
| Base quality | ≥ 13 | < 13: base not used | Provisional. |

MalariaGEN Pf7 sets the allele fraction "to 0 if there are < 2 or <5% alternative
allele reads" and does not compute it "when there were less than 5 reads in total".

The remaining caller settings — screening and PASS operating points, strand-bias limit,
error-model defaults, pair-HMM gap penalties, mixture-model limits — are specific to
PfVarX-Call. No external publication can supply them. They stay
`UNCALIBRATED_PROVISIONAL` until the benchmarking experiments named in
`config/parameters/layer2.yaml` are run. The same file also holds the candidate admission count
(2 reads), the pair-HMM read cap (200) and flank (60 bp), and the shortest focal callability
loss reported (50 bp).

The error model leaves out of its fit any position where a non-reference allele has at least
2 reads and at least 5% of reads (`error_fit_exclude_fraction`), and corrects for this with a
truncated likelihood. Without it, minor-clone variants in mixed infections inflated the error
model and were screened out (on synthetic 40× mixtures, 3% of minor-clone variants were
called before and 85% after). Mixture-constrained calling (`mixture_constrained_calling`) is
implemented but off: after the error-model fix it added 34–75 true and 4–17 false calls on
the same mixtures, depending on its threshold.

The mixture model leaves out SNVs within 10 bp of a called indel, including a deletion's span
(`mixture_indel_exclusion_bp`). Reads across an indel are sometimes aligned with mismatches
instead of a gap, which creates false low-fraction SNVs beside true indels; a few of them were
enough to make a single-clone sample fit two or three clones. On the synthetic 40× cohort, 3 of
6 single-clone samples were estimated as 2–3 clones before and all 6 as one clone after; the
mixed samples' estimates did not change. Calls are not affected. The number of sites left out
is recorded in the caller QC (`mixture.snvs_excluded_near_indels`).

## Stage 3 — cohort evidence

| State at a cohort site | Rule | Source |
| --- | --- | --- |
| `ALT_DOMINANT` | alternative allele ≥ 95% of reads | MalariaGEN Pf7 (the other allele is below 5%). |
| `MIXED` | between 5% and 95% | MalariaGEN Pf7. |
| `REF_DOMINANT` | alternative allele ≤ 5% of reads | MalariaGEN Pf7. |
| `MISSING_LOW_DEPTH` | fewer than 5 usable reads | MalariaGEN Pf7. |
| `UNCALLABLE` | outside the sample's callable intervals | Framework policy. |

A missing record is never counted as a reference call.

## Stage 5 — downstream analyses

| Analysis | Decision | Value | Source |
| --- | --- | --- | --- |
| DA001 COI | Allele counted as present | ≥ 5% of reads | MalariaGEN Pf7. |
| DA001 COI | Polyclonal / monoclonal call | P(COI > 1) ≥ 0.95 / ≤ 0.05; otherwise indeterminate | Provisional. |
| DA001 COI | Minimum informative loci | 10; below: indeterminate | Provisional. |
| DA001 COI | MOIRE burn-in / samples per chain | 10 000 / 1 000 | MOIRE defaults (Murphy and Greenhouse 2024). |
| DA001 COI | Fws comparison line | Fws ≥ 0.95 read as clonal | Manske 2012; Auburn 2012; MalariaGEN Pf7 and Pf8. Descriptive only. |
| DA001 / DA003 | Samples needed for MOIRE and Fws | ≥ 10 samples with enough informative loci (`population_frequency_min_samples`). Fewer: MOIRE is not run; the monoclonal/polyclonal call comes from the caller's per-sample mixture test (`basis` = `CALLER_MIXTURE_LRT`); COI count and Fws are left blank | Provisional. MOIRE estimates allele frequencies and COI jointly: on a 9-sample synthetic cohort, P(COI > 1) of truly single-clone samples ranged from 0.00 to 1.00 across five random seeds (mixtures stayed at 1.00). Fws uses frequencies in which each sample is 1/n of its own reference (one sample: Fws 0 by construction). |
| DA004 relatedness | COI given to Dcifer | DA001's MOIRE count; monoclonal samples 1; otherwise Dcifer's naive estimate (`getCOI`, `lrank = 2`, `da004_dcifer_coi_lrank`) | Gerlovina 2022. |
| DA002 resistance | Marker callable | ≥ 5 reads spanning the codon; fewer: `UNKNOWN` | MalariaGEN Pf7. |
| DA002 resistance | `MUTANT` | mutant fraction ≥ 95% | MalariaGEN Pf7. |
| DA002 resistance | `MIXED` | 5% to 95% | MalariaGEN Pf7. |
| DA002 resistance | `BELOW_REPORTING_FLOOR` | mutant reads present but under 5% | MalariaGEN Pf7. |
| DA002 resistance | `who_k13_alert` in `prevalence.tsv` | validated PfKelch13 marker in more than 5% of callable infections | WHO (see below). |
| DA002 resistance | Callable infections needed for the K13 flag | ≥ 3 per group; fewer: `INSUFFICIENT_CALLABLE` (`da002_min_callable_per_group`) | Provisional. One mutant in one or two samples would otherwise raise the WHO flag. |
| DA003 diversity | Minimum group size | 3 | Provisional. |
| DA004 relatedness (WGS) | Markers needed for an estimate | ≥ 200 biallelic; fewer: `EXCLUDED_FEW_MARKERS` | Taylor 2019 (about 200 biallelic markers keep the error below 0.1). |
| DA004 relatedness (targeted) | Markers needed for an estimate | ≥ 20; 20–199: `ESTIMATED_LOW_PRECISION`, reported without a label; fewer: `EXCLUDED_FEW_MARKERS` | Provisional. No published minimum exists for SNPs that are few and linked within amplicons; the 200-marker precision level is Taylor 2019. |
| DA004 relatedness | `CLOSELY_RELATED` / `NEAR_IDENTICAL` labels | IBD ≥ 0.5 / ≥ 0.9 | Reporting convention (e.g. Pereus 2026). |
| DA003 diversity | Bootstrap replicates for intervals | 200 | Framework policy. |
| DA004 relatedness | Dcifer significance level | 0.05 | Dcifer usage (Gerlovina 2022). |
| DA004 relatedness | Dcifer estimate grid (`nr`) | 1 000 | Dcifer default. |
| DA005 structure | LD pruning window / step / r² | 100 / 20 / 0.1 | scikit-allel defaults. |
| DA005 structure | Minimum group size | 3 | Provisional. |
| DA005 / DA006 on a panel | Separate target regions (independent loci) | ≥ 10 (`da005_min_independent_loci`, `da006_min_independent_loci`); fewer: `TOO_FEW_INDEPENDENT_LOCI` | Provisional: the same floor as DA006 informative sites and DA001 informative loci. The core genome is not limited. |
| DA002 resistance | Territory and callable evidence | The core genome or panel contains ≥ 1 whole catalogue codon (else `PANEL_LACKS_MARKERS`); each sample has ≥ 1 whole catalogue codon in its callable regions (else `NO_CALLABLE_MARKERS`) | Framework policy: a marker outside the sequenced territory has no denominator. |
| DA006 phylogeny | Minimum monoclonal samples | 4 | Mathematical floor for an informative unrooted tree. |
| DA006 phylogeny | Sample / site missingness allowed | ≤ 20% | Provisional. |
| DA006 phylogeny | Bootstrap replicates | UFBoot 1 000, SH-aLRT 1 000 | IQ-TREE. |
| DA006 phylogeny | Branch considered supported | UFBoot ≥ 95 and SH-aLRT ≥ 80 | IQ-TREE guidance (Hoang 2018; Guindon 2010). |
| DA007 spatial | Minimum samples per unit | 3 | Provisional. |
| DA007 spatial | Samples with place and date before a GAM is fitted | 10 | Provisional. |
| DA007 spatial | GAM basis dimension / interval level | k ≤ 4 / 95% | Framework policy. |
| DA008 copy number | Amplified: depth ratio / lower bootstrap bound | ≥ 1.5 / ≥ 1.25 | Provisional (halfway between one and two copies). |
| DA008 copy number | Deleted: depth ratio, with flank ratio | ≤ 0.1, flanks ≥ 0.3 | Provisional. |
| DA008 copy number | Present (deletion loci) / single-copy floor | ratio ≥ 0.3; between 0.1 and 0.3 is UNKNOWN | Provisional. |
| DA008 copy number | Median core-genome depth | ≥ 20 | Provisional; below it every locus is UNKNOWN. |
| DA008 copy number | Windows / flanks / GC model / bootstrap | 200 bp / 2 000 bp / 3 000 core windows / 1 000 replicates | Framework policy. |

Notes on interpretation:

- **Low-frequency alleles.** The 5% floor is the whole-genome convention. Amplicon panels
  report lower limits under validated conditions: MAD4HatTeR detected minor alleles at
  within-sample frequencies "as low as 1%" in high-density samples (Aranda-Díaz 2025), and
  Katairo 2025 reported 100% sensitivity at 2% for samples of 1 000 parasites/µL or more.
  Those limits belong to those assays; PfVarX keeps 5% for every design until its own
  detection limit is measured.
- **Relatedness labels.** Taylor 2019 recommends "∼100 polyallelic or 200 biallelic
  markers" to keep the error of a relatedness estimate below 0.1. The 0.5 and 0.9 labels
  are a reporting convention, not tests; relatedness is continuous, and Mehra 2025 shows
  that estimates are biased downwards on sparse data while high thresholds stay robust.
- **Branch support.** IQ-TREE: "you should only start to rely on a branch if its support
  is >= 95%" for UFBoot, and "one would be more confident if a clade has its SH-aLRT >=
  80% and UFboot >= 95%". The developers state this for single-gene trees.

## Public-health reference levels

These are not pipeline parameters, except the first, which sets a flag. PfVarX reports
marker prevalence; it never declares resistance.

| Reference level | Value | Source |
| --- | --- | --- |
| Validated PfKelch13 marker prevalence that prompts follow-up | more than 5% of patients | WHO Malaria Policy Advisory Group background document, March 2024. |
| Evidence needed to confirm artemisinin partial resistance | the marker level above, plus delayed parasite clearance in a therapeutic efficacy study | WHO. |
| Criteria for a PfKelch13 marker to be "validated" | association with clearance half-life > 5 hours, and > 1% survival in the ring-stage survival assay in at least five isolates | WHO. |
| Validated PfKelch13 markers | F446I, N458Y, C469Y, M476I, Y493H, R539T, I543T, P553L, R561H, P574L, C580Y, R622I, A675V | WHO. These are the 13 kelch13 entries in the shipped catalogue. |

The shipped catalogue remains `DRAFT_REQUIRES_VERIFICATION`. Check it against the current
WHO list before use; the list changes.

## Where values are kept

Every threshold, analysis setting and compute request is read from a configuration file; the
program code holds none of them, and `tests/unit/test_no_hardcoded_parameters.py` fails if a
parameter is repeated in the code as a fallback.

| File | Holds | Who edits it |
| --- | --- | --- |
| `config/analysis.yaml`, `config/quality.yaml`, `config/compute.yaml` | the project's choices: inputs, design, outputs, analyses, QC policy, computer limits | the user (CLI or GUI) |
| `config/registry/qc_policy_*.yaml` | QC gate thresholds | maintainer |
| `config/parameters/layer0.yaml` … `layer6.yaml` | every numeric setting of input checks, read processing, the caller, cohort, annotation, analyses and report, with its source class | maintainer |
| `config/registry/job_resources.yaml` | threads and memory of each workflow rule, and the planning factors (memory headroom, temporary-disk factor, progress refresh) | maintainer |
| `config/registry/*.tsv`, `downstream_analyses.yaml` | marker catalogues and analysis contracts | maintainer |

A run freezes the parameters it used (`frozen_config.json`). A run frozen before a parameter
was added takes the value it was frozen with from `src/pfvarx/configuration/legacy_parameters.yaml`;
that file records history and is not a place to change behaviour.

These settings of the numerical methods are design decisions (`FRAMEWORK_POLICY`, or
`SOFTWARE_DEFAULT` for a tool's own default) and are listed here for completeness:

| Parameter | Value | What it does |
| --- | --- | --- |
| `layer0.core_mask_fraction_min` / `_max`, `core_mask_max_intervals` | 0.70 / 0.97, 1000 | a core mask outside these limits is flagged `CORE_MASK_SUSPECT` |
| `layer0.alignment_state_reads_checked` | 2000 | reads read from an imported BAM/CRAM to check its processing state |
| `layer2.homopolymer_tract_classes` | 1, 2–4, 5–9, ≥ 10 bp | homopolymer classes with their own indel error model |
| `layer2.indel_error_init_divisor`, `indel_error_floor` | 5, 10⁻⁶ | starting value (error_mean_default / 5) and floor of the indel error rate |
| `layer2.error_rate_cap` | 0.25 | largest error rate used for WSAF intervals and constrained scoring |
| `layer2.error_fit_min_sites`, `error_fit_min_nonzero_for_phi` | 50, 50 | sites needed to fit the error model and its overdispersion |
| `layer2.overdispersion_phi_bounds` | 10 – 10⁷ | range of the fitted overdispersion |
| `layer2.wsaf_profile_grid_points` | 2001 | grid of the WSAF profile likelihood (step 0.0005) |
| `layer2.indel_anchor_bp` | 5 | reference bases a read must cover on each side of an indel to be re-scored |
| `layer2.pileup_max_depth` | 1 000 000 | pileup depth cap, far above any real depth |
| `layer2.evidence_chunk_bp` | 250 000 | work unit for parallel evidence collection; results do not depend on it |
| `layer2.k_bootstrap_max_iterations`, `k_bootstrap_restart_divisor`, `k_bootstrap_min_restarts` | 200, 10, 2 | EM limits of the bootstrap replicate fits in the clone-number test |
| `layer5.da001_moire_pt_chains` | 1 | MOIRE parallel-tempering chains |
| `layer5.da003_interval_level`, `da003_min_bootstrap_values` | 0.95, 10 | heterozygosity bootstrap interval |
| `layer5.da004_dcifer_afreq_tol` | 10⁻⁵ | Dcifer allele-frequency tolerance |
| `layer5.da005_min_samples`, `da005_min_sites`, `da005_pca_components` | 3, 3, 4 | PCA minimums and components reported |
| `layer5.da007_gam_grid_points`, `da007_time_origin` | 20, 2000-01-01 | GAM curve points and time axis origin |
| `layer5.da008_interval_level`, `da008_gc_min_windows_per_bin` | 0.95, 20 | depth-ratio interval and GC-matched windows per expected depth |

What stays in code is not a setting: file-format facts (VCF, GFF and BED columns, BAM flags,
the Phred+33 offset), mathematical constants (three alternative bases, the normal quantile
computed from a level), unit conversions, numerical-safety bounds of the optimisers, and the
status and reason vocabularies in `src/pfvarx/constants.py`.

## What still has no published basis

These values stay provisional because no standard exists, not because they were
overlooked: read retention, mapped fraction, duplicate fraction, on-target fraction,
amplicon dropout, the warning levels for depth and callability, minimum mapping and base
quality, minimum group sizes, phylogeny missingness limits, the clonality credibility
level, and every PfVarX-Call model setting. Properly-paired fraction and UMI singleton
fraction have no value at all and are reported only.

## References

- Aranda-Díaz A, et al. Sensitive and modular amplicon sequencing of *Plasmodium falciparum* diversity and resistance for research and public health. *Scientific Reports* 2025.
- Auburn S, et al. Characterization of within-host *Plasmodium falciparum* diversity using next-generation sequence data. *PLoS ONE* 2012.
- Chen S, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 2018.
- Early AM, et al. Detection of low-density *Plasmodium falciparum* infections using amplicon deep sequencing. *Malaria Journal* 2019.
- Fulcrum Genomics. fgbio best practice consensus pipeline. <https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md>
- Gerlovina I, et al. Dcifer: an IBD-based method to calculate genetic distance between polyclonal infections. *Genetics* 2022.
- Guindon S, et al. New algorithms and methods to estimate maximum-likelihood phylogenies. *Systematic Biology* 2010.
- Hoang DT, et al. UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 2018.
- IQ-TREE. Frequently asked questions. <https://iqtree.github.io/doc/Frequently-Asked-Questions>
- Katairo T, et al. Performance of molecular inversion probe DR23K and Paragon MAD4HatTeR amplicon sequencing panels. *Malaria Journal* 2025.
- Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997, 2013.
- MalariaGEN, et al. Pf7: an open dataset of *Plasmodium falciparum* genome variation in 20,000 worldwide samples. *Wellcome Open Research* 2023. Methods: <https://www.malariagen.net/wp-content/uploads/2023/11/Pf7-Details-of-bioinformatics-methods.pdf>
- MalariaGEN, et al. Pf8: an open dataset of *Plasmodium falciparum* genome variation in 33,325 worldwide samples. *Wellcome Open Research* 10:325, 2025. doi:10.12688/wellcomeopenres.24031.1
- MAD4HatTeR pipeline documentation. <https://eppicenter.github.io/mad4hatter/>
- Manske M, et al. Analysis of *Plasmodium falciparum* diversity in natural infections by deep sequencing. *Nature* 2012.
- Mehra S, et al. Systematic bias in malaria parasite relatedness estimation. *G3* 2025.
- Miles A, et al. Indels, structural variation, and recombination drive genomic diversity in *Plasmodium falciparum*. *Genome Research* 2016.
- Murphy M, Greenhouse B. MOIRE: a software package for the estimation of allele frequencies and effective multiplicity of infection from polyallelic data. *Bioinformatics* 2024.
- Oyola SO, et al. Whole genome sequencing of *Plasmodium falciparum* from dried blood spots using selective whole genome amplification. *Malaria Journal* 2016.
- Pereus D, et al. Genetic metrics decode *Plasmodium falciparum* diversity: complexity of infection, parasite connectivity, and transmission intensity in Mainland Tanzania's diverse regions. *Frontiers in Genetics* 2026.
- Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors in unique molecular identifiers to improve quantification accuracy. *Genome Research* 2017.
- Taylor AR, et al. Estimating relatedness between malaria parasites. *Genetics* 2019.
- World Health Organization. Malaria: artemisinin partial resistance (questions and answers). <https://www.who.int/news-room/questions-and-answers/item/artemisinin-resistance>
- World Health Organization Malaria Policy Advisory Group. Antimalarial drug resistance in Africa, background document, March 2024. <https://cdn.who.int/media/docs/default-source/malaria/mpac-documentation/mpag-march2024-session5-antimalarial-drug-resistance-africa.pdf>
