"""Layer 3 variant-only cohort merge that works at any cohort size.

``bcftools merge`` keeps every input VCF open at once, so a cohort larger than the process's
open-file limit fails ("Too many open files"): about 1 000 samples under the common Linux
limit of 1 024, about 250 under the macOS default of 256. This raises the limit for the merge
to what the cohort needs, up to the system's hard limit, and then runs the same merge as
before, so every cohort that fitted gives the same file. Only a cohort beyond the hard limit
is merged in batches: each batch to a temporary BCF, then the batches. The batched cohort VCF
holds the same records, genotypes and sample order; records that share a position may be
listed in a different order. The strategy used is written to the log.

The file list is passed to bcftools in a file (``-l``), never on the command line, so the
number of samples is not limited by the operating system's command-line length either.
"""

from __future__ import annotations

import argparse
import math
import resource
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

MERGE = ["bcftools", "merge", "-m", "none", "-i", "-"]


def open_file_limit(needed: int) -> int:
    """Raise this process's soft open-file limit towards *needed*; return the limit in force.

    Children (bcftools) inherit it. The hard limit is never exceeded.
    """
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    ceiling = needed if hard == resource.RLIM_INFINITY else min(needed, hard)
    if ceiling > soft:
        try:
            resource.setrlimit(resource.RLIMIT_NOFILE, (ceiling, hard))
        except (ValueError, OSError):  # some systems refuse values below the hard limit too
            pass
    return resource.getrlimit(resource.RLIMIT_NOFILE)[0]


def _merge(files: list[Path], out: Path, output_type: str, workdir: Path) -> None:
    """One bcftools merge of *files* (copied when there is a single file)."""
    if len(files) == 1:
        subprocess.run(["bcftools", "view", "-O", output_type, "-o", str(out), str(files[0])], check=True)
    else:
        listing = workdir / f"{out.name}.list"
        listing.write_text("".join(f"{f}\n" for f in files))
        subprocess.run([*MERGE, "-O", output_type, "-o", str(out), "-l", str(listing)], check=True)
    subprocess.run(["bcftools", "index", "-f", *(["-t"] if output_type == "z" else []), str(out)], check=True)


def merge(vcfs: list[Path], out: Path, headroom: int) -> str:
    """Merge per-sample VCFs (in the given order) into *out* (bgzipped, tabix-indexed).

    *headroom* is ``planning.cohort_merge_file_headroom``: open files the merge needs besides
    its inputs (output, standard streams, libraries). Returns ``copy``, ``direct`` or ``batched:<size>``.
    """
    if len(vcfs) == 1:  # one sample: the cohort is its VCF, copied unchanged
        shutil.copyfile(vcfs[0], out)
        subprocess.run(["bcftools", "index", "-t", "-f", str(out)], check=True)
        return "copy"
    limit = open_file_limit(len(vcfs) + headroom)
    with tempfile.TemporaryDirectory(prefix=".merge-", dir=out.parent) as tmp:
        workdir = Path(tmp)
        if len(vcfs) + headroom <= limit:
            _merge(vcfs, out, "z", workdir)
            return "direct"
        size = limit - headroom
        if size < 2:
            raise SystemExit(f"open-file limit {limit} leaves no room to merge (headroom {headroom}); "
                             "raise the limit with 'ulimit -n'")
        batches = [vcfs[i:i + size] for i in range(0, len(vcfs), size)]
        if len(batches) > size:
            raise SystemExit(f"{len(vcfs)} samples need more than two merge levels under an open-file limit "
                             f"of {limit}; raise the limit with 'ulimit -n'")
        width = len(str(len(batches)))
        parts = []
        for number, batch in enumerate(batches):
            part = workdir / f"batch_{number:0{width}d}.bcf"
            _merge(batch, part, "b", workdir)
            parts.append(part)
        _merge(parts, out, "z", workdir)
        return f"batched:{size}"


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--list", required=True, help="file with one per-sample VCF per line, in cohort order")
    ap.add_argument("--out", required=True)
    ap.add_argument("--headroom", type=int, required=True)
    a = ap.parse_args(argv)
    vcfs = [Path(line) for line in Path(a.list).read_text().splitlines() if line.strip()]
    if not vcfs:
        raise SystemExit("no per-sample VCFs to merge")
    if shutil.which("bcftools") is None:
        raise SystemExit("bcftools is not installed")
    strategy = merge(vcfs, Path(a.out), a.headroom)
    soft = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    print(f"cohort merge: {len(vcfs)} sample(s), strategy {strategy}, open-file limit {soft}",
          file=sys.stderr, flush=True)
    if strategy.startswith("batched"):
        print(f"cohort merge: {math.ceil(len(vcfs) / int(strategy.split(':')[1]))} batch(es); records sharing a "
              "position may be listed in a different order than a single merge would give", file=sys.stderr)


if __name__ == "__main__":
    main()
