"""Build a SnpEff database from the pinned FASTA+GFF, check it, and annotate (Layer 4).

The genome identifier, reference checksum and GFF checksum are written to a build-check
record; annotation refuses to run when the record does not match the files.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from pfvarx.io import sha256


def _snpeff() -> str:
    exe = shutil.which("snpEff")
    if not exe:
        raise SystemExit("TOOL_NOT_FOUND: snpEff (install the 'full' pixi environment)")
    return exe


def build(fasta: Path, gff: Path, genome: str, data_dir: Path, check: Path) -> None:
    gdir = data_dir / "data" / genome
    gdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(fasta, gdir / "sequences.fa")
    shutil.copyfile(gff, gdir / "genes.gff")
    cfg = data_dir / "snpEff.config"
    cfg.write_text(f"data.dir = {data_dir / 'data'}\n{genome}.genome : {genome}\n"
                   f"{genome}.codonTable : Standard\n")
    r = subprocess.run([_snpeff(), "build", "-gff3", "-noCheckCds", "-noCheckProtein", "-c", str(cfg),
                        "-v", genome], capture_output=True, text=True)
    print(r.stdout[-4000:], r.stderr[-4000:])
    if r.returncode != 0:
        raise SystemExit(f"snpEff build failed with exit code {r.returncode}")
    check.write_text(json.dumps({"genome": genome, "fasta_sha256": sha256(fasta), "gff_sha256": sha256(gff),
                                 "snpeff_version": subprocess.run([_snpeff(), "-version"], capture_output=True,
                                                                  text=True).stdout.strip(),
                                 "status": "BUILT"}, indent=2))


def annotate(vcf: Path, genome: str, data_dir: Path, check: Path, ud: int, out_vcf: Path, out_qc: Path) -> None:
    info = json.loads(check.read_text())
    gdir = data_dir / "data" / genome
    if info["genome"] != genome or sha256(gdir / "sequences.fa") != info["fasta_sha256"]:
        raise SystemExit("REFERENCE_MISMATCH: SnpEff database does not match the pinned reference (FAIL)")
    cmd = [_snpeff(), "ann", "-c", str(data_dir / "snpEff.config"), "-noStats", "-ud", str(ud),
           "-noLog", genome, str(vcf)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr)
        raise SystemExit(f"snpEff ann failed with exit code {r.returncode}")
    with gzip.open(out_vcf, "wt") as h:
        h.write(r.stdout)
    effects: Counter[str] = Counter()
    records = 0
    for line in r.stdout.splitlines():
        if line.startswith("#"):
            continue
        records += 1
        for field in line.split("\t")[7].split(";"):
            if field.startswith("ANN="):
                for ann in field[4:].split(","):
                    effects[ann.split("|")[1]] += 1
    out_qc.write_text(json.dumps({"records": records, "consequence_counts": dict(effects),
                                  "genome": genome, "database_check": info,
                                  "cross_annotator": "VEP cross-check not run (SPECIFIED)",
                                  "haplotype_aware_csq": "BCFtools/csq not run (SPECIFIED)"},
                                 indent=2))


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--annotate", action="store_true")
    ap.add_argument("--fasta")
    ap.add_argument("--gff")
    ap.add_argument("--vcf")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--check", required=True)
    ap.add_argument("--ud", type=int, default=0)
    ap.add_argument("--out-vcf")
    ap.add_argument("--out-qc")
    a = ap.parse_args(argv)
    if a.annotate:
        annotate(Path(a.vcf), a.genome, Path(a.data_dir), Path(a.check), a.ud, Path(a.out_vcf), Path(a.out_qc))
    else:
        build(Path(a.fasta), Path(a.gff), a.genome, Path(a.data_dir), Path(a.check))


if __name__ == "__main__":
    main()
