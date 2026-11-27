# Select, qualify and replace the Layer 2 caller

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
Its decision values (screening and PASS thresholds, the pair-HMM support ratio, the linkage
fragment minimum, the genotype cut-off and the EM settings) are read from
`config/parameters/layer2.yaml`. The pair-HMM runs as a Numba-compiled kernel when Numba is
installed and as the numpy reference otherwise, with the same read decisions; worker
processes share read evidence, error-model fits, indel rescoring, linkage windows and
bootstrap replicates, and the output does not depend on their number.
Full local haplotype assembly and validated assay-specific calibration remain open.
Technical mixture count describes DNA evidence. Parasite genomes remain haploid.

The implemented assay branches are WGS, targeted non-UMI and targeted UMI. Targeted
barcode data can enter a targeted branch only when its inputs and protocol satisfy
that branch's contract. Shallow WGS uses the WGS branch with applicable QC and
callability gates. Neither TBS nor sWGS has an independently validated preset here.
No established third-party caller is bundled or prequalified by the adapter interface.
