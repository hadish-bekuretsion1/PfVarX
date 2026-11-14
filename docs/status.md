# Implementation and scientific status

This distribution contains the complete **implemented** PfVarX 0.0.0 software, including
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
- WGS, targeted non-UMI and targeted UMI branches, plus existing BAM/CRAM entry contracts.
- PfVarX-Call read evidence, reduced per-sample beta-binomial error models, candidate
  screening, local pair-HMM indel rescoring, read-backed linkage, technical mixture/K
  selection, WSAF intervals, callability and six site-status classes.
- Variant-only cohort merge and forced per-sample cohort-site evidence, with explicit
  missing/uncallable states and an explicit outcome when no sample enters the cohort;
  reference-matched SnpEff annotation.
- DA001 MOIRE; DA002 codon marker genotyping; DA003 diversity; DA004 hmmIBD/Dcifer;
  DA005 haploid PCA/Hudson FST; DA006 IQ-TREE; DA007 descriptive/GAM spatial-temporal work.
- Offline dashboard, exports, benchmark partitions, exact-normalized accuracy/mixture
  evaluation, wild-cluster uncertainty and resource-profile generation.
- Typed-source checks, unit/scientific/integration tests (including a synthetic cohort
  that runs all seven analyses), version and CLI/GUI consistency tests, isolated synthetic
  end-to-end verification, browser checks, and clean-distribution tests.

## Limits that matter before scientific benchmarking

| Area | Current implementation / limitation |
|---|---|
| Empirical settings | `UNCALIBRATED_PROVISIONAL`; associated outputs remain UNCALIBRATED. UNSET_BLOCKING optional metrics are NOT_EVALUATED, not PASS. |
| Resistance catalogue | Shipped catalogue is DRAFT_REQUIRES_VERIFICATION; real DA002 is blocked until a maintainer independently verifies a matching catalogue. TEST_FIXTURE catalogues are only for software tests. |
| Caller error model | Reduced substitution/homopolymer classes and fitted dispersion; full covariate/random-effect model is not implemented. |
| Haplotypes | Local reference/alternative pair-HMM and linkage blocks; bounded combinatorial enumeration and pass-two mixture factor remain proposed extensions. |
| Mixture K | Technical DNA-mixture count, not biological COI. Analysis fallback bootstrap count is provisional; precomputed critical-value surface is not implemented. |
| WSAF | Profile-likelihood intervals; representation-bias correction remains unimplemented. |
| Clonality | MOIRE implemented, descriptive FWS retained; THE REAL McCOIL and DEploid remain unimplemented alternatives. MCMC settings/convergence require scientific review. |
| Relatedness | hmmIBD/Dcifer implemented; isoRelate is not included. WGS pairs need 200 informative markers (Taylor 2019); the targeted minimum of 20 is provisional, and pairs below 200 are reported as low precision without a label. |
| Annotation | SnpEff implemented; VEP cross-check and BCFtools/csq consequences are not included. |
| Exact benchmark matcher | Normalized allele identity, not haplotype-aware equivalence. RTG availability is reported but vcfeval is manual; hap.py/baseline pipelines are not installed/run automatically. |
| Confidence/callability | Shared confidence BED; compare documented metric denominators carefully. Truth validation is partly an external experiment responsibility. |
| Held-out protection | Partition hash and refusal of repeated attempts in one evaluation directory; no global access-control/experiment lock or family-aware splitting. |
| Resampling intervals | Unrestricted wild cluster bootstrap of the linearized ratio estimator; not restricted wild cluster bootstrap-t. |
| Execution | Linux x86-64 local executor, including inside a single SLURM node allocation; the `executor` setting is recorded but every value runs on one host. Per-rule cluster or cloud scheduling and object storage are not implemented. `tmp_dir` can point job temporary files at node-local scratch. CRAM and targeted existing-alignment entry have less validation than WGS BAM. VCF entry is explicitly rejected; use FASTQ/BAM/CRAM. |
| Resources/performance | Resource requests/heuristic estimates are not OS limits. Caller is single-worker per sample and can take substantial full-genome time without intermediate logs: about 30 s for the 48 kb synthetic core at 40×, which extrapolates to hours per sample for the ~21 Mb 3D7 core (not yet measured on real data). Measured profiles need manual basis review. |
| Reuse/provenance | Frozen software/configuration integrity is checked; not all external reads, masks or truth files are frozen or covered by reuse keys. Maintain checksum manifests; use NEW for changed scientific inputs. |
| Reporting | Artifact-derived progress requires report regeneration; use finalized runs for interpretation. Missing upstream distributions/timestamps stay unavailable. Historical control charts are not implemented. PDF uses browser printing and visible table pages. |
| Spatial interpretation | Descriptive aggregation/GAM only; preferential sampling is not modelled. |
| Unsupported scope | Structural variants, copy-number/resistance amplification inference, clone-resolved trees and a validated ML classifier are not delivered by this implementation. |

The source release retains these limits explicitly instead of shipping simulated
resources or successful-looking historical outputs as evidence of readiness. Read the
[benchmarking protocol](benchmarking/README.md) to establish real truth, freeze decisions,
and assess performance before promoting calibrated settings into routine analysis.
