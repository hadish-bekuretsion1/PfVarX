# Installation, deployment and maintenance

## Supported native installation

Linux x86-64 / local executor is supported. Git/source access, Python 3.11+ for bootstrap
and Pixi are prerequisites. The lockfile covers only `linux-64`: on Windows use WSL2 with a
Linux distribution; macOS (Intel or Apple silicon) and Linux ARM64 are not supported, because
several Bioconda tools have no builds there and none are locked. Private GitHub repositories require authorized login. Clone
alone cannot run setup; use an explicit post-clone command:

```bash
git clone --filter=blob:none --sparse https://github.com/hadish-bekuretsion1/PfVarX.git PfVarX-source
python3 PfVarX-source/bootstrap.py --profile sample
cd pfvarx-analysis
./pfvarx --help
```

Full extraction/clone also works. Bootstrap fetches/exports required blobs and installs
the matching locked environment. A selected workspace contains required source,
templates/static metadata, its launcher and installation capability declaration. It
initializes actual input/metadata folders, not artificial empty run/output trees.
Choose profile/design/input type/analyses before installation. See
[orchestration](orchestration.md) for every option, stage and input/output contract.

Use `--profile development` for all methods/tests/scripts. A standalone Python wheel
is not a complete pipeline: workflow, policy, resource and tool dependencies are
required. Do not copy `.pixi` between filesystem prefixes. Installing creates a local
environment; source distribution itself contains none. Use `./pfvarx` in selected
workspaces rather than accidentally installing the full default with plain `pixi run`.

## Dependencies and scientific resources

`pixi.toml` declares direct feature requirements; `pixi.lock` pins Linux package URLs,
versions and hashes. The one source manifest declares valid combinations; only the
selected environment is installed. Base dependencies include the Python control plane,
NumPy/SciPy/pysam, Snakemake, SAMtools/BCFtools/HTSlib, Picard and mosdepth. FASTQ adds
BWA-MEM2/fastp; targeted FASTQ adds fgbio. FO004 adds SnpEff, then exactly selected
method feature closures: MOIRE/R/jsonlite, Dcifer/hmmIBD, scikit-allel, IQ-TREE and/or
mgcv/R. Optional `--rtg` adds RTG; dev adds testing/lint/type tools. Transitive packages
needed by selected tools remain. Browser QA uses separate Playwright/openpyxl/Pillow.

The software supplies no FASTQ/BAM, genome, GFF, core mask, primers or truth. Obtain the
correct original release, record identifiers/URLs/licenses/checksums and independently
verify scientific compatibility. Local absolute paths avoid duplication. Discovery
creates candidates requiring human metadata review; registration references or links
files. `fetch-input` downloads only an explicitly requested HTTPS file with a trusted
checksum. No automatic "latest genome" or guessed mask is used.

## Optional Apptainer builder — not yet image-qualified

The source contains a real selected-profile build/context generator. Apptainer is not
available in the verification host, so no successful image build/container E2E is
claimed. Use native Pixi for verified deployment. A site with Apptainer may review and
qualify the recipe:

```bash
# From the full source directory, with the pinned Pixi binary on PATH:
python3 scripts/build_container.py --profile sample --emit-context /absolute/path/new-context
# Review new-context/Apptainer.def and build_manifest.json.
python3 scripts/build_container.py --profile sample --image /absolute/path/pfvarx.sif
# Run using a matching source workspace mounted read-write at /workspace:
apptainer run --cleanenv --bind /absolute/path/study:/workspace /absolute/path/pfvarx.sif --help
```

The image provides tools under `/opt/pfvarx`; the mounted workspace supplies current
source/config/runtime, with `PYTHONPATH=/workspace/src`. Mount external original inputs
read-only at the exact declared paths as needed. Use matching reviewed source and image
lock identities. Do not assume a build/help test proves native/container equivalence:
run all six end-to-end and browser scenarios and compare scientific outputs. Record
image SHA-256, Pixi hash, source/lock hashes, host/version and resolved base OCI digest.
The emitted base tag `ubuntu:24.04` must be pinned/recorded by digest for strict image
reproducibility before qualification. See the single-node SLURM procedure below.

## GitHub release and upgrades

Publish the contents of this pristine source folder as the repository root, not an
outer study workspace. `.gitignore` excludes runtime paths, mutable study config,
installation metadata, workspace locks, caches/environments and build artifacts.
### Versioning

The pipeline release version has one source: `__version__` in `src/pfvarx/__init__.py`.
The package metadata reads it from there (`[tool.hatch.version]` in `pyproject.toml`);
`./pfvarx version`, VCF headers, run snapshots, reports and the source archive name all
take it from the same place. `pixi.toml` deliberately has no version. To release, change
`__version__`, add a changelog section with that version, run the tests
(`tests/unit/test_version.py` fails on any disagreeing current-release reference) and tag
the reviewed commit `v<version>`. Schema, parameter-file, QC-policy, registry and
catalogue versions are independent and are not changed by a release.

Before publication run [testing](testing.md), isolated workflow checks, source package
checks and checksum/extraction validation. The CI pins an explicit Ubuntu runner and
validates the software-only archive for the current version. Attach the final source
ZIP/checksum to the reviewed release. Never attach a study workspace or a guessed
scientific-validation claim.

Keep the previous source/lock for rollback. Use NEW workspaces for capability/assay/input
changes and independent input checksums for changed external resources. Old frozen
runs remain immutable. An old reference manifest is resealed/repaired via `prepare`
before new execution; RESUME still requires unchanged effective source/configuration.
Back up scientific evidence before any manual cleanup; the automated policy never
removes originals/final products or bulk-cleans failed-run work.

## CPU deployment profiles

For an 8 GB Linux PC, start with the supplied template (one concurrent job,
6500 MiB workflow budget, up to four threads). On existing workspaces, edit
both memory limits and set max_jobs to 1 in config/compute.yaml. Set max_cores
and per_job.max_threads to the available CPU count if it is lower. Run
`./pfvarx ready` and `./pfvarx dry-run` before launching a study. The budget
reserves space for the OS. Snakemake resource requests govern scheduling,
while the OS or scheduler enforces actual memory limits. Deep libraries and
large cohorts can require a larger machine.

For a Linux server or cloud VM, install the same locked environment and set
the compute limits for the VM or container allocation. Use a persistent
filesystem for the workspace and registered inputs. The planner checks
physical memory, the standard Linux cgroup v2 hierarchy and SLURM per-node
memory when supplied. CPU reporting uses process affinity on Linux.

For HPC, submit the local executor inside one SLURM node allocation. From the
prepared workspace root, with Pixi available to batch jobs:

```bash
sbatch --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=12:00:00 \
  --wrap='./pfvarx ready && ./pfvarx dry-run && ./pfvarx pipeline --mode NEW'
```

Match the batch request and compute YAML to the study. Set `tmp_dir` in
`config/compute.yaml` to node-local scratch (an absolute path such as `/local/scratch/alice`,
created before the job; YAML does not expand `$USER`) so Snakemake, Picard and fgbio temporary files do not fill a small `/tmp`.
This uses one node per workflow run. Per-rule distributed scheduling and managed cloud batch
services require separate executor integration and qualification.
