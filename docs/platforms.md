# Native platforms and qualification

PfVarX uses the same terminal commands, configuration GUI and CPU caller on each native
installation target. A terminal is an interface, not an operating-system compatibility layer.

| Target | Pixi platform | Installation status |
| --- | --- | --- |
| Linux x86-64 | `linux-64` | Existing baseline; previous package selections retained |
| Linux ARM64 (AArch64) | `linux-aarch64` | Dependencies locked; experimental pending native CI |
| macOS Intel | `osx-64` | Dependencies locked; experimental pending native CI |
| macOS Apple silicon | `osx-arm64` | Dependencies locked; experimental pending native CI |

The locks declare macOS 13+, or Linux kernel 4.18+ with glibc 2.28+. Use a native
terminal, Python and Pixi for your architecture. On Apple silicon, avoid running the
installation under a Rosetta Intel terminal if you intend to use the ARM64 build.
32-bit ARM, native Windows and other architectures are outside these native targets.
Windows users can use WSL2 with one of the supported Linux architectures.

## Install from a terminal

Use the [getting-started guide](getting-started.md) unchanged: install the reviewed Pixi
version, clone/extract the source, then run `python3 PfVarX-source/bootstrap.py`.
Bootstrap needs Python 3.11+; Pixi installs the locked Python 3.12 environment.
Pixi chooses the host platform automatically; no aligner or scientific model is substituted.
`./pfvarx configure --gui`, `./pfvarx install-tools` and the pipeline commands are identical.

On every new computer, export a fresh workspace and install from its lock. Do not copy
`.pixi` between computers, architectures or installation prefixes. Source and input files
can be transferred separately. Check reference indexes and prepared files again on the
destination rather than assuming architecture-specific indexes are interchangeable.

The GUI's tool table and CLI installation checks select versions from the host's own
lock. Keep the full lock with the run: dependency builds and floating-point libraries can
differ by platform. The existing provenance snapshot records platform and package versions.

## Qualification and limitations

All 118 environments resolve for all four targets. That establishes package availability,
not a successful pipeline run or equivalent scientific results on new hardware.
macOS/ARM targets remain experimental until native integration checks pass.

The manually dispatched **PfVarX CI** workflow runs the four native platforms independently:
locked installation, platform identity, lint, typing, unit/scientific tests, full synthetic
integration runs, and source-archive checks. It also tests the spawn process model used by
macOS. Review the results on the exact source revision before claiming runtime support.
Synthetic checks do not replace independent scientific calibration or held-out validation.
CPU architectures can differ numerically, especially near decision thresholds; no bitwise
cross-platform equivalence or unchanged accuracy is claimed from dependency solving.

macOS total memory is detected for resource limits. The current available-memory warning
uses Linux `MemAvailable` and is unavailable on macOS; leave headroom for the OS and other
applications when setting `max_mem_mb`. Apple silicon CPU and GPU share physical memory.

These additions provide CPU execution targets. They add no CUDA, Metal or other GPU
backend. A run does record the graphics processors it found on its machine (`nvidia-smi` on
Linux and Windows hosts, `system_profiler` on macOS) in its provenance, with
`used_for_compute: false`, and the report says so in the Software section: PfVarX runs
correctly on a machine with a GPU, uses its CPUs, and no result depends on an accelerator.
Every pinned tool (fastp, bwa-mem2, samtools, fgbio, mosdepth, SnpEff, IQ-TREE, hmmIBD and the
R packages) is a CPU program, and the native caller uses NumPy and SciPy. Cloud/server execution still uses one host; HPC examples still run inside one
allocated Linux node. The optional Apptainer image builder is a Linux deployment route,
not a native macOS installer or a qualified multi-architecture container distribution.
