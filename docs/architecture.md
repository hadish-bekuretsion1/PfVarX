# Architecture and ownership

Layer 2 exposes a versioned sample artifact contract. The workflow can load the native
caller or invoke a separately installed adapter. The same downstream consumers read
the resulting artifacts. The configuration editor, input validation, preparation
and execution are separate entry points.

| Boundary | Owner | Contract |
| --- | --- | --- |
| Post-clone export/install | `bootstrap.py` | Standard-library/Python 3.11+ entry point; new destination only; full/sparse source support; selected source/locked environment |
| Environment dispatch | root `pfvarx` launcher | Run from workspace root; use installation metadata to choose Pixi environment; `tool` dispatch for installed utilities |
| Configuration/capabilities | `configuration/`, templates, installation metadata | Strict typed config; reject unavailable modes/analyses/assays/input types and unsupported VCF entry |
| Optional configuration editor | `configuration/web.py`, `configuration/gui/` | Loopback form uses the same schema and capability checks, saves YAML with backup protection, and creates input folders |
| Layout | `initialization/layout.py`, `workspace.py` | Config-driven input folders, three editable YAMLs, sample header and init metadata; no fake runs/results |
| Registration/acquisition | `inputs/organize.py` | Review-required discovery, reviewed TSV imports, references/symlinks, HTTPS checksum/no-clobber publication |
| Validation/preparation | `inputs/`, `references/`, `execution/preflight.py` | Actual input facts, detailed reports, effective config/parameters, necessary derived reference artifacts and verification |
| Scheduling/execution | `execution/runner.py`, Snakemake | Mandatory prep integrity/resource gates, frozen run, checkpoint scheduling, bounded local concurrency, failure isolation |
| Retention | `execution/retention.py` | Workspace lock; keep/notemp or dependency cleanup; completed-run safe work-only cleanup; preview-first CLI |
| Caller boundary | `caller_adapter.py`, `caller/`, `genomic_evidence.py` | Native or external backend emits VCF, evidence, callable/focal BED and QC under the same contract |
| Downstream eligibility | `eligibility/`, method registry and Layer 5 checkpoint | Assay candidate selection, QC admission, DA001 clonality, metadata, artifact and method gates with reason codes |
| Provenance | `runs/`, `provenance/`, diagnostics | Frozen config/plan/software plus launcher/installation metadata, attempts/events/logs; external file checksums remain a study responsibility |
| Reporting | `reporting/` | Offline HTML, dashboard JSON /2 and QC artifacts; renderer hash separate from frozen scientific provenance |

## Selection and dependency closure

FO001–FO003 installations need no downstream implementation (FO001 also leaves out the caller and its tools). Selected FO004 analyses
add annotation and requested modules plus genuine prerequisites/shared helpers. DA004,
DA005 and DA006 depend on DA001; DA007 currently uses DA002 outcomes; DA008 and DA009 have no
prerequisite, and DA009 reuses DA008's depth model for whole-genome evidence. Registry/policy
metadata remains shared because control/provenance code reads it. Tool feature groups
are reads, targeted, annotation, coi, related, structure, phylogeny, spatial, bench and
dev. Valid feature combinations are locked in the single Pixi manifest/lock. Core is
FASTQ sample/cohort; aligned skips FASTQ alignment tools. Optional features are installed
only in the selected environment; a full/development environment retains every tool.

The source export omits unselected modules/rules/R adapters/unused geometry. Snakemake
includes cohort/annotation/downstream rules only for layers in the validated plan, and
UMI rules only for UMI FASTQ. The downstream dispatcher imports only the invoked method.
Changing objectives beyond a selected export requires a NEW appropriate workspace;
installation metadata is a capability declaration, not a security/access-control system.

## Preparation and execution lifecycle

Clone → explicit bootstrap → review configuration → discover/register/fetch actual
inputs → validate → prepare → verify-prepared → ready/dry-run → execute → report →
retention/explicit cleanup. `orchestrate` implements the ordered chain after inputs are
configured and can stop at initialize/validate/prepare/verify. It holds a workspace lock
through mutation/execution. Explicit pipeline entry also locks and re-enforces preflight.

Reference manifest v2 stores source fingerprints and nonempty derived artifact hashes.
Separate semantic checks compare derived FASTA index/dictionary names and lengths.
Corruption/missing/stale artifacts prevent run creation; `prepare` repairs affected
resources. BWA indexes are required only for FASTQ entry. Annotation-specific resources
are deferred to the actual FO004 layer and checked against its reference/GFF there.

No setup relies on clone hooks or implicit execution. No reference/truth resource is
chosen/downloaded without explicit user input. No genotype, assay chemistry, geographic
coordinate or wet-lab library identity is inferred as scientific truth by discovery.

## Data boundaries

Runtime paths are `inputs/`, `prepared/`, `metadata/`, `runs/`, `results`, benchmarks,
installed `.pixi` and the workspace lock; they are excluded from source archives.
The completed run retains BAM/VCF/evidence/callability/QC/frozen/log/measurement/report
artifacts. Only safe regular files under `runs/ID/results/work` are bulk cleanup
candidates; originals, symlinks and failed-run evidence are excluded. An explicit
successful retention policy controls Snakemake temporary-output deletion as well.

Only native Linux local execution has been qualified here. The optional Apptainer
context/image builder is a deployment implementation with an explicit qualification
gap, not evidence of a tested container image or cluster deployment.

Layer 2 selects a backend behind the [caller contract](caller.md). Cohort, reporting,
and benchmark consumers use the existing sample artifact paths. Shared read-filter
and interval conventions live outside the native caller package.

See the [configuration editor and input procedure](configuration.md) and the
[sequencing and downstream eligibility matrix](sequencing.md).

## Cohort evidence in memory

Layer 5 analyses that need a genotype at every cohort site read `cohort/evidence_matrix.tsv.gz`
into memory. The cohort's sample order is known before the file is read, so each site holds
columns addressed by the sample's position: one pass, constant-time lookup and no intermediate
copy. That costs a measured 25 bytes per cell against 473 for a dict per cell, and `wsaf` is
kept as a float64 exactly as the file holds it, so no comparison against a reporting floor can
move. Each analysis declares whether it reads the matrix (`USES_COHORT_MATRIX`); the four that
work from the alignments never pay for it, and an analysis that does not declare itself is given
the matrix, so a new one cannot lose its evidence by omission
(`tests/unit/test_matrix_declarations.py`).

## Microhaplotype loci

Layer 3 writes one multiallelic marker per catalogue locus (`cohort/microhaplotypes.tsv.gz`)
from the reads of each sample, using this cohort's PASS SNVs inside the locus: no variant is
discovered there, so a haplotype cannot contain a site the caller did not pass. DA001 and DA004
prefer those loci, because MOIRE and Dcifer are built for multiallelic markers, and fall back to
biallelic SNVs when the cohort has too few loci. The choice, the locus count and the reason are
recorded in each analysis manifest (`locus_basis`), never left implicit.

## Where values are kept

No threshold, analysis setting or compute request is written in the code. User choices are in
the three user YAML files; maintainer values are in `config/parameters/layer*.yaml` (each with
its source class), `config/registry/qc_policy_*.yaml` and `config/registry/job_resources.yaml`
(threads and memory per workflow rule, read by both the rules and the preflight plan). Code
reads a run's values from its frozen configuration; a run frozen before a parameter existed
uses `configuration/legacy_parameters.yaml` through `loader.frozen_param`. To add a setting,
add it to the matching file with a source class and a note, read it with `frozen_param`, and
list it in [thresholds](thresholds.md#where-values-are-kept).
`tests/unit/test_no_hardcoded_parameters.py` and `tests/unit/test_resource_plan.py` enforce it.

## CLI and GUI parity

The terminal commands and the local configuration GUI are two front ends for the same
workspace: the same YAML files, the same validation and the same pipeline stages. A change
to one must be made to the other in the same commit.

| When you change | Also change |
| --- | --- |
| A setting in `configuration/models.py` | Its entry in `gui/field-guide.json`, both project templates (`config/templates/analysis.yaml`, `benchmarking.yaml`) or the quality/compute template, then run `python scripts/configuration_guide.py` to refresh `docs/configuration-guide.md` and `.html` |
| The help text of a setting | `gui/field-guide.json` only (it is the single source), then run `python scripts/configuration_guide.py` |
| A GUI tab or button | Its CLI command in `cli.py` and the "Terminal equivalent" line in `gui/index.html` |
| A CLI configuration or preflight command | The GUI tab or stage that runs it (`configuration/web.py`, `gui/app.js`) |
| Shared rules (paths, transfers, validation) | Put them in one Python module both use, never in the GUI JavaScript alone |

`tests/unit/test_cli_gui_parity.py` fails when the two sides drift: a setting missing from
the GUI help, the templates or the guide; guide text that differs from the GUI; a GUI tab
naming a command that does not exist; or a GUI stage the CLI lacks.
