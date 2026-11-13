# Sequencing design and conditional analysis

Declare the sequencing protocol in `assay.sequencing_type` and choose its processing
branch in `assay.design`. The pipeline checks compatibility. Read depth and filenames
alone cannot identify a laboratory protocol. Older configurations can omit the protocol
label and continue to use their declared branch.

| Sequencing protocol | Processing branch and required evidence | DA001 COI and clonality | Applicable downstream analyses |
| --- | --- | --- | --- |
| WGS (`WGS`) | `WGS`. Declare direct or SWGA amplification, provide the matching reference and core mask, and review usable depth and breadth. | MOIRE can estimate COI when the cohort supplies enough informative loci and the fit qualifies. | DA002–DA007 are candidates. Each retains its marker, QC, callability, metadata and clonality gates. |
| Shallow WGS (`SWGS`) | `WGS`. Apply the WGS gates and review the observed low-depth exclusions. | Same model and evidence gates as WGS. Low information can leave clonality indeterminate or prevent estimation. | Same candidates as WGS. There is no separate validated shallow-depth preset. |
| Targeted amplicon (`AMPLICON`) | `TARGETED_NON_UMI` or `TARGETED_UMI` as actually sequenced. Provide target BED and primer TSV. For UMI reads, declare the actual read structure. | Conditional on informative callable loci and an adequate fit. Review panel ascertainment and linkage. | DA002, DA003, DA004 and DA007 are candidates. Report diversity and marker denominators within the measured panel. |
| Targeted barcode with UMI (`TBS_UMI`) | `TARGETED_UMI` when paired reads, panel coordinates, primer schema and UMI structure fit this workflow. Require usable consensus evidence. | Same targeted eligibility, with panel-specific qualification required. | DA002, DA003, DA004 and DA007 are candidates under the targeted contract. |
| Targeted barcode without UMI (`TBS_NON_UMI`) | `TARGETED_NON_UMI` when paired reads, panel coordinates and primer schema fit this workflow. Review amplification and usable-read evidence. | Same targeted eligibility. Review minor-clone sensitivity and panel information. | DA002, DA003, DA004 and DA007 are candidates under the targeted contract. |

TBS labels identify compatible targeted inputs. A barcode-only allele table has no
supported entry route. Sequencing type and clonality are distinct. COI is the
infection-level model result. The caller's technical mixture estimate is a separate
quantity. A monoclonal classification requires qualifying evidence. Polyclonal
samples can enter supported mixture-aware paths. Indeterminate or missing
clonality cannot satisfy a gate that requires an assessed category.

| ID | Analysis and prerequisite | Clonality and execution conditions | Expected result when eligible |
| --- | --- | --- | --- |
| DA001 | COI and clonal status. Needs cohort allele evidence, admitted samples and informative loci. | Estimates clonality. Review model convergence and information checks. | MOIRE COI summaries and `clonality.tsv`, or a recorded exclusion. |
| DA002 | Resistance markers. Needs matching annotation, callability and a verified marker catalogue. | Supports assessed and indeterminate infections. Coverage and marker-specific QC still apply. | Callable marker states and supported prevalence summaries. |
| DA003 | Genetic diversity. Needs cohort variants and callability. | Review mixed infections and panel scope. Group comparisons require qualifying population groups. | FWS and diversity summaries at the requested supported scale. |
| DA004 | Relatedness/IBD. Needs DA001 and an eligible sample pair with enough shared information. | hmmIBD for qualifying monoclonal pairs. Dcifer for qualifying pairs involving polyclonal samples. | Pair estimates and explicit exclusions. |
| DA005 | Population structure. Needs WGS, DA001, cohort variants and population labels. | Assessed mono/poly samples subject to the existing missingness, sample and group gates. | Haploid PCA and eligible Hudson FST summaries. |
| DA006 | Single-lineage phylogeny. Needs WGS, DA001, callability and sufficient informative sites. | Monoclonal samples only, with the configured minimum qualifying sample count. | IQ-TREE tree and support results. |
| DA007 | Spatial/temporal analysis. Needs DA002 outcome, country, coordinates and collection date. | Inherits source eligibility and applies geographic, date and model support checks. | Qualifying descriptive summaries and GAM results where supported. |

Outputs live under `runs/RUN_ID/results/downstream/DAxxx/`. Review each
`status_manifest.json`, the eligibility matrix and cohort decisions before using
scientific output. The dashboard reports excluded methods and terminal states.
Insufficient evidence, missing optional metadata and unsupported designs can
legitimately result in no estimate. Inspect the reason and exclusions.

## Select automatic execution

For a fresh targeted non-UMI workspace with the compatible dependencies:

```bash
python3 bootstrap.py --workspace ../targeted-study --profile downstream \
  --design TARGETED_NON_UMI --auto-analyses
cd ../targeted-study
./pfvarx configure --gui
```

Choose FO004 and automatic selection in the form, or edit the equivalent YAML.

```yaml
schema_version: pfvarx.config/3
final_output: FO004
analyses:
  selection: auto
  requested: []
  scales: {}
  population_column: population
```

This is a configuration excerpt. Complete the assay, input, reference and caller
fields using the generated template. Supply a matching GFF to enable annotation and marker outcomes. In automatic mode,
an unset GFF omits Layer 4 and blocks DA002 and DA007 with recorded reasons. An
explicit FO004 configuration still requires a matching GFF. The plan
selects all design-compatible methods, uses each registry default scale and resolves
DA007's source to DA002. Use explicit selection to request other supported scales.
An existing selective workspace must contain the dependencies for the full
candidate set. Bootstrap a matching workspace when its capability set is smaller.

With automatic selection, absent optional population and geographic fields do not
invalidate the whole sample registry. Eligibility and method checks record their
effect. Keep `CONTINUE_INDEPENDENT` to allow qualified independent analyses to
continue. `STRICT` retains its stop behavior. FO001–FO003 retain their respective
output boundaries. FO001 requires no caller. FO002–FO004 can use the native or a
qualified external caller through the existing Layer 2 artifact contract.

## Upgrade identifiers

| Method | Config /2 and archived runs | Config /3 and new runs |
| --- | --- | --- |
| COI and clonality | DA003 | DA001 |
| Resistance markers | DA001 | DA002 |
| Genetic diversity | DA002 | DA003 |
| Relatedness, population structure, phylogeny, spatial/temporal | DA004–DA007 | DA004–DA007 |

Config /2 YAML is translated at load time, including requested IDs, scale keys
and the spatial source. Leave the old schema marker in place unless you also
translate those IDs. Start a NEW run for upgraded results. Keep archived result
directories and dashboards with the original release. Dashboard JSON /2 uses the
new identifier order. Use the original renderer for archived dashboard JSON /1.

Explicit FO004 runs the selected methods and their declared prerequisites. The workflow evaluates COI first, and methods that require assessed clonality use that evidence. Independent methods can continue under CONTINUE_INDEPENDENT when COI cannot be estimated.
