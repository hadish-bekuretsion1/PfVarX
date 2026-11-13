"""Build a SnpEff database from the pinned FASTA+GFF, check it, and annotate (Layer 4).

The genome identifier, reference checksum and GFF checksum are written to a build-check
record; annotation refuses to run when the record does not match the files.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import pysam

from pfvarx.io import sha256

# NCBI translation table 4 as shipped in SnpEff's own snpEff.config ('+' = start, '*' = stop).
# Plasmodium mitochondrial genes (cox1, cox3, cytb) read TGA as Trp; the standard table would
# report a premature stop there.
MITOCHONDRIAL_TABLE = "Protozoan_Mitochondrial"
MITOCHONDRIAL_CODONS = (
    "TTT/F, TTC/F, TTA/L+, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, "
    "TGT/C, TGC/C, TGA/W, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, "
    "CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I+, ATC/I+, ATA/I+, ATG/M+, "
    "ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, "
    "GTT/V, GTC/V, GTA/V, GTG/V+, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, "
    "GGT/G, GGC/G, GGA/G, GGG/G"
)
# PlasmoDB names the mitochondrial contig Pf3D7_MIT_v3 (older releases: M76611).
MITOCHONDRIAL_CONTIG = re.compile(r"(^|_)(MIT|MITO|MT)(_|$)|^M76611$|^chrM$", re.IGNORECASE)


def _snpeff() -> str:
    exe = shutil.which("snpEff")
    if not exe:
        raise SystemExit("TOOL_NOT_FOUND: snpEff (install an environment with the 'annotation' feature)")
    return exe


def mitochondrial_contigs(fasta: Path) -> list[str]:
    names = []
    with fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if MITOCHONDRIAL_CONTIG.search(name):
                    names.append(name)
    return names


def write_config(data_dir: Path, genome: str, mitochondrial: list[str]) -> Path:
    """Write snpEff.config with a data directory relative to the file, so a workspace can be moved."""
    lines = ["data.dir = ./data/", f"{genome}.genome : {genome}"]
    if mitochondrial:
        lines.append(f"codon.{MITOCHONDRIAL_TABLE} : {MITOCHONDRIAL_CODONS}")
        lines += [f"{genome}.{name}.codonTable : {MITOCHONDRIAL_TABLE}" for name in mitochondrial]
    cfg = data_dir / "snpEff.config"
    cfg.write_text("\n".join(lines) + "\n")
    return cfg


def build(fasta: Path, gff: Path, genome: str, data_dir: Path, check: Path) -> None:
    gdir = data_dir / "data" / genome
    gdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(fasta, gdir / "sequences.fa")
    shutil.copyfile(gff, gdir / "genes.gff")
    mitochondrial = mitochondrial_contigs(fasta)
    cfg = write_config(data_dir, genome, mitochondrial)
    r = subprocess.run([_snpeff(), "build", "-gff3", "-noCheckCds", "-noCheckProtein", "-c", str(cfg),
                        "-v", genome], capture_output=True, text=True)
    print(r.stdout[-4000:], r.stderr[-4000:])
    if r.returncode != 0:
        raise SystemExit(f"snpEff build failed with exit code {r.returncode}")
    check.write_text(json.dumps({"genome": genome, "fasta_sha256": sha256(fasta), "gff_sha256": sha256(gff),
                                 "snpeff_version": subprocess.run([_snpeff(), "-version"], capture_output=True,
                                                                  text=True).stdout.strip(),
                                 "codon_tables": {name: MITOCHONDRIAL_TABLE for name in mitochondrial},
                                 "status": "BUILT"}, indent=2))


def annotate(vcf: Path, genome: str, data_dir: Path, check: Path, ud: int, out_vcf: Path, out_qc: Path) -> None:
    info = json.loads(check.read_text())
    gdir = data_dir / "data" / genome
    if (info["genome"] != genome or sha256(gdir / "sequences.fa") != info["fasta_sha256"]
            or sha256(gdir / "genes.gff") != info["gff_sha256"]):
        raise SystemExit("REFERENCE_MISMATCH: SnpEff database does not match the pinned reference/GFF (FAIL)")
    # Rewrite the configuration so a database built before the workspace moved stays usable.
    cfg = write_config(data_dir, genome, sorted(info.get("codon_tables", {})))
    cmd = [_snpeff(), "ann", "-c", str(cfg), "-noStats", "-ud", str(ud), "-noLog", genome, str(vcf)]
    plain = out_vcf.with_name(out_vcf.name + ".tmp.vcf")
    try:
        with plain.open("w") as handle:
            r = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            print(r.stderr)
            raise SystemExit(f"snpEff ann failed with exit code {r.returncode}")
        effects: Counter[str] = Counter()
        records = 0
        # BGZF (not plain gzip) so the annotated VCF can be indexed and queried by region.
        with plain.open() as source, pysam.BGZFile(str(out_vcf), "wb", index=None) as target:
            for line in source:
                target.write(line.encode())
                if line.startswith("#"):
                    continue
                records += 1
                for field in line.split("\t")[7].split(";"):
                    if field.startswith("ANN="):
                        for ann in field[4:].split(","):
                            effects[ann.split("|")[1]] += 1
    finally:
        plain.unlink(missing_ok=True)
    pysam.tabix_index(str(out_vcf), preset="vcf", force=True)
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
