# Implementation and scientific status

This distribution contains the complete **implemented** PfVarX 1.0.0 software, including
subsequent type-safety and source-distribution cleanup. It is ready for controlled
benchmark experiments. It is not a claim that every proposed research method exists,
that empirical thresholds are calibrated, or that biological accuracy is established.
No datasets, pilot results or historical analysis evidence are bundled.

Apptainer build/context generation is implemented but container execution has not been
qualified here because Apptainer is unavailable. Only native Linux/Pixi execution is
verified; no automatic cluster/container deployment claim is made.

## Implemented components

- Explicit post-clone selective bootstrap, locked analysis/input-specific environments,
  minimal source exports, reviewed discovery/zero-copy registration, checksum-gated
  HTTPS acquisition, integrity-enforced reference verification and safe retention.
- Typed YAML and optional local HTML configuration, kept in step by a parity test;
  automatic assay candidate selection, sample/reference validation that explains each
  failure and its fix, conditional preparation, QC policy with cited thresholds
  ([thresholds](thresholds.md)), resource planning, scientific eligibility, frozen runs,
  diagnostics and orchestration.
- WGS, targeted non-UMI and targeted UMI branches selected by `assay_design`, each with FASTQ and
  BAM/CRAM entry; imported alignments run only the steps their documented `processing_state`
  still needs (verified end to end on the synthetic fixture for every design and state).
- PfVarX-Call read evidence, reduced per-sample beta-binomial error models, candidate
  screening, local pair-HMM indel rescoring, read-backed linkage, technical mixture/K
  selection, WSAF intervals, callability and six site-status classes.
- Variant-only cohort merge and forced per-sample cohort-site evidence, with explicit
  missing/uncallable states and an explicit outcome when no sample enters the cohort;
  reference-matched SnpEff annotation.
- DA001 MOIRE; DA002 codon marker genotyping; DA003 diversity; DA004 hmmIBD/Dcifer;
  DA005 haploid PCA/Hudson FST; DA006 IQ-TREE; DA007 descriptive/GAM spatial-temporal work;
  DA008 GC-corrected read-depth copy number at catalogue loci (WGS); DA009 pfhrp2/pfhrp3
  joint deletion on every design, by depth ratio for WGS and by amplicon presence against
  control amplicons for targeted panels.
- Offline dashboard, exports, benchmark partitions, exact-normalized accuracy/mixture
  evaluation, wild-cluster uncertainty and resource-profile generation.
- Typed-source checks, unit/scientific/integration tests (including a synthetic cohort
  that runs all nine analyses), version and CLI/GUI consistency tests, browser checks and
  clean-distribution tests, run during development; that test suite is not part of this
  distribution. Isolated synthetic end-to-end verification ships as `scripts/verify.py`.

## Limits that matter before scientific benchmarking

| Area | Current implementation / limitation |
|---|---|
| Empirical settings | `UNCALIBRATED_PROVISIONAL`; associated outputs remain UNCALIBRATED. UNSET_BLOCKING optional metrics are NOT_EVALUATED, not PASS. |
| Resistance catalogue | Shipped catalogue is DRAFT_REQUIRES_VERIFICATION; real DA002 is blocked until a maintainer independently verifies a matching catalogue. TEST_FIXTURE catalogues are only for software tests. |
| Caller error model | Reduced substitution/homopolymer classes and fitted dispersion; full covariate/random-effect model is not implemented. |
| Haplotypes | Local reference/alternative pair-HMM and linkage blocks; bounded combinatorial enumeration and pass-two mixture factor remain proposed extensions. |
| Mixture K | Technical DNA-mixture count, not biological COI. SNVs within 10 bp of a called indel are left out of the mixture sites, because alignment artefacts beside indels made single-clone samples fit 2–3 components (3 of 6 synthetic single-clone samples before, 0 of 6 after). An 85/15 synthetic mixture is still estimated as 3 components. Analysis fallback bootstrap count is provisional; precomputed critical-value surface is not implemented. |
| WSAF | Profile-likelihood intervals; representation-bias correction remains unimplemented. |
| Clonality | MOIRE implemented, descriptive FWS retained; THE REAL McCOIL and DEploid remain unimplemented alternatives. MCMC settings/convergence require scientific review. With fewer than 10 informative samples MOIRE is not run (its allele frequencies and COI are not identified; on a 9-sample truth cohort single-clone samples got P(COI>1) from 0.00 to 1.00 depending on the seed). DA001 then takes the monoclonal/polyclonal call from the caller's per-sample mixture test and withholds the COI count and FWS. Whether 10 samples is enough for MOIRE has not been tested. |
| Relatedness | hmmIBD/Dcifer implemented; isoRelate is not included. WGS pairs need 200 informative markers (Taylor 2019); the targeted minimum of 20 is provisional, and pairs below 200 are reported as low precision without a label. |
| Annotation | SnpEff implemented; VEP cross-check and BCFtools/csq consequences are not included. |
| Exact benchmark matcher | Normalized allele identity, not haplotype-aware equivalence. RTG availability is reported but vcfeval is manual; hap.py/baseline pipelines are not installed/run automatically. |
| Confidence/callability | Shared confidence BED; compare documented metric denominators carefully. Truth validation is partly an external experiment responsibility. |
| Held-out protection | Partition hash and refusal of repeated attempts in one evaluation directory; no global access-control/experiment lock or family-aware splitting. |
| Resampling intervals | Unrestricted wild cluster bootstrap of the linearized ratio estimator; not restricted wild cluster bootstrap-t. |
| Execution | Local executor with locks for Linux x86-64, Linux ARM64 and macOS 13+ (Intel/Apple silicon); macOS and Linux ARM64 are experimental pending native qualification (see [platform support](platforms.md)). Linux can run inside a single SLURM node allocation; the `executor` setting is recorded but every value runs on one host. Per-rule cluster or cloud scheduling and object storage are not implemented. `tmp_dir` can point job temporary files at node-local scratch. CRAM and targeted existing-alignment entry have less validation than WGS BAM. VCF entry is explicitly rejected; use FASTQ/BAM/CRAM. |
| Cohort scale (Layer 5) | The analyses that read cohort evidence hold one cell per sample per cohort site: measured at 25 bytes per cell (473 before the columnar store, an 18.8x reduction), so 200 samples over 5 000 sites is about 25 MB and 500 over 50 000 about 630 MB. DA002, DA007, DA008 and DA009 work from the alignments and never load it. The eligibility checkpoint warns (`RESOURCE_MEMORY_TIGHT`) before the analyses run when the estimate passes half of `per_job.max_mem_mb`. Evidence is still read whole into memory, not streamed, so a cohort of several thousand samples over a dense site set remains beyond one host. |
| Accelerators | Graphics processors are detected and recorded in provenance; nothing computes on them. No CUDA, Metal or other GPU backend, and no result depends on an accelerator. |
| Resources/performance | Resource requests/heuristic estimates are not OS limits. The caller runs read evidence (250 kb pieces, `layer2.evidence_chunk_bp`), the error-model fits, indel rescoring (one candidate per task), linkage windows, bootstrap replicates and EM restarts in parallel worker processes, with output identical to a serial run. Measured on a real 29× WGS sample (20.8 Mb core), before the Numba pair-HMM and the faster evidence collection: 66 min serial, 9.7 min with 8 workers and 7.9 min with 24; memory 3.6–3.8 GB whatever the worker count. The rule requests 8 threads and 5 000 MB, raised to the peak measured on the same computer in earlier real-data runs plus 20 % (within `per_job.max_mem_mb`); `workflow-check` warns when the memory free at check time cannot hold `max_jobs` jobs. The log reports each stage. Time, memory and storage estimates per step are learned from earlier real-data runs on the same computer and are UNKNOWN until such a run exists; on three public WGS samples (5.5 GB of BAM) the estimate was 2 h 20 min and the runs took 2 h 48–50 min, mostly in variant calling. The remaining time is mostly the mixture bootstrap. Measured profiles need manual basis review. |
| Reuse/provenance | Frozen software/configuration integrity is checked; not all external reads, masks or truth files are frozen or covered by reuse keys. Maintain checksum manifests; use NEW for changed scientific inputs. |
| Reporting | A live progress page (`runs/<RUN_ID>/progress.html`) follows a running workflow; the report dashboard is written when the run ends, so use finalized runs for interpretation. Missing upstream distributions/timestamps stay unavailable. Historical control charts are not implemented. PDF uses browser printing and visible table pages. |
| Spatial interpretation | Descriptive aggregation/GAM only; preferential sampling is not modelled. |
| Copy number (DA008) | Gene-level amplification and deletion at six catalogue loci from GC-corrected read depth, WGS only; the shipped copy-number catalogue is DRAFT_REQUIRES_VERIFICATION and blocks DA008 until verified. Thresholds are UNCALIBRATED_PROVISIONAL and not yet benchmarked against samples of known copy number; mapping bias in subtelomeric genes (pfhrp2/3) is not corrected; amplicon breakpoints are not resolved (junction read pairs are counted as support only). |
| Microhaplotypes | Read-backed haplotypes at catalogue loci, or at loci derived from the run's own panel targets or core-genome windows (`analyses.microhap_loci`). Derived loci are specific to that definition, so frequencies compare only with cohorts defined the same way; a published panel's own catalogue remains the better input where one exists. Read-backed haplotypes at supplied catalogue loci, used by DA001 and DA004 as multiallelic markers. Haplotypes are limited to this cohort's PASS SNVs inside each locus, so the panel's own documented variable positions are not independently verified; a fragment must resolve every position or it is incomplete. The floors (fragment depth, haplotype reads and share) are UNCALIBRATED_PROVISIONAL. The shipped panel holds no loci: coordinates belong to a published panel and a reference release, so a site supplies and verifies its own. No phasing across loci and no haplotype-level clone reconstruction. |
| pfhrp2/3 thresholds | The targeted dropout and whole-genome deletion cut-offs are UNCALIBRATED_PROVISIONAL and not yet benchmarked against samples of known pfhrp2/3 status. Two things stand without them: `background_percentile`, which places the gene's depth ratio among the sample's own core-genome windows and carries no threshold - read across samples at the same locus, since the background is core sequence while these genes are subtelomeric - and `threshold_sensitivity.tsv`, which recomputes every state across a grid of cut-offs and counts the calls that move. Neither replaces calibration against truth samples. |
| pfhrp2/3 deletion (DA009) | Joint single and dual deletion, their co-occurrence against independence and a country-month trend. The WGS evidence is the DA008 depth measurement; a targeted panel gives amplicon presence against control amplicons, which is not a copy number, and a mutated primer site or a failed PCR produces the same dropout as a deletion. The two evidence classes are recorded per row and prevalence is published pooled and stratified; results from runs with different designs are comparable only with that composition in view. Thresholds are UNCALIBRATED_PROVISIONAL and the catalogue must be verified. Not a patient-level diagnosis and not evidence of selection. |
| Unsupported scope | Genome-wide structural-variant and copy-number discovery, clone-resolved trees and a validated ML classifier are not delivered by this implementation. |

The source release retains these limits explicitly instead of shipping simulated
resources or successful-looking historical outputs as evidence of readiness. Read the
[benchmarking protocol](benchmarking/README.md) to establish real truth, freeze decisions,
and assess performance before promoting calibrated settings into routine analysis.

