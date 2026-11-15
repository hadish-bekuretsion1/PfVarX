"""Dispatch Layer 2 through the versioned sample output contract."""

from __future__ import annotations

import argparse
import csv
import importlib
import importlib.util
import os
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import pysam

from pfvarx.configuration.models import CallerConfig
from pfvarx.constants import CallStatus
from pfvarx.io import read_json, sha256, write_json

EVIDENCE_COLUMNS = (
    "chrom pos ref alt kind y n alt_f alt_r ref_f ref_r mean_alt_bq tract "
    "screen status wsaf wsaf_lo wsaf_hi lhid pairhmm_reads"
).split()
INPUTS = ("bam", "reference", "core-mask", "frozen-config", "sample", "context-unreliable")
OUTPUTS = ("vcf", "evidence", "callable", "focal", "qc")


def check_backend(config: CallerConfig) -> None:
    if config.backend == "native":
        if importlib.util.find_spec("pfvarx.caller") is None or importlib.util.find_spec("pfvarx.caller.call") is None:
            raise ValueError("Native caller is unavailable. Configure an external caller.")
        return
    executable = config.executable
    if executable is None or not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError("External caller executable is missing or is not executable")
    if sha256(executable) != config.sha256:
        raise ValueError("External caller executable checksum mismatch")


def validate_outputs(paths: dict[str, Path], sample: str, reference: Path) -> None:
    with pysam.FastaFile(str(reference)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths, strict=True))
        with pysam.VariantFile(str(paths["vcf"])) as vcf:
            if list(vcf.header.samples) != [sample]:
                raise ValueError("Caller VCF must contain the requested sample only")
            required_format = {"GT", "AD", "DP", "PST", "WSAF", "WSAFLO", "WSAFHI"}
            if not required_format.issubset(vcf.header.formats):
                raise ValueError("Caller VCF is missing required FORMAT declarations")
            for record in vcf:
                if record.contig not in lengths or record.start < 0 or record.stop > lengths[record.contig]:
                    raise ValueError("Caller VCF contains an out-of-reference record")
                ref = record.ref
                if ref is None or fasta.fetch(record.contig, record.start, record.start + len(ref)).upper() != ref.upper():
                    raise ValueError("Caller VCF REF does not match the reference")
                fields = record.samples[sample]
                gt = fields.get("GT")
                if not gt or len(gt) != 1 or gt[0] not in (None, 0, 1):
                    raise ValueError("Caller VCF requires a haploid GT per biallelic record")
                if len(record.alts or ()) != 1 or not all(
                    allele and set(allele.upper()) <= set("ACGTN")
                    for allele in (record.alleles or ())
                ):
                    raise ValueError("Caller VCF requires explicit biallelic SNP or indel alleles")
                depth, counts = fields.get("DP"), fields.get("AD")
                if (
                    not isinstance(depth, int) or depth < 0 or not counts or len(counts) != 2
                    or any(not isinstance(count, int) or count < 0 for count in counts)
                    or sum(counts) > depth
                ):
                    raise ValueError("Caller VCF requires nonnegative AD and DP with AD sum no greater than DP")
                lo, value, hi = (fields.get(key) for key in ("WSAFLO", "WSAF", "WSAFHI"))
                if lo is None or value is None or hi is None or not 0 <= lo <= value <= hi <= 1:
                    raise ValueError("Caller VCF contains missing or invalid allele-fraction bounds")
                if fields.get("PST") not in set(CallStatus):
                    raise ValueError("Caller VCF contains an invalid decision status")
        for kind in ("callable", "focal"):
            with paths[kind].open() as handle:
                previous: dict[str, int] = {}
                for line in handle:
                    chrom, start, end = line.rstrip("\n").split("\t")
                    lo, hi = int(start), int(end)
                    if not previous.get(chrom, 0) <= lo < hi <= lengths.get(chrom, 0):
                        raise ValueError(f"Invalid or overlapping {kind} interval")
                    previous[chrom] = hi
    with paths["evidence"].open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if rows.fieldnames != EVIDENCE_COLUMNS:
            raise ValueError("Caller evidence TSV has an incompatible header")
        for row in rows:
            if (
                None in row or None in row.values()
                or not 1 <= int(row["pos"]) <= lengths.get(row["chrom"], 0)
                or row["status"] not in set(CallStatus)
                or not 0 <= int(row["y"]) <= int(row["n"])
                or not 0 <= float(row["wsaf_lo"]) <= float(row["wsaf"]) <= float(row["wsaf_hi"]) <= 1
            ):
                raise ValueError("Caller evidence contains an invalid row")
    qc = read_json(paths["qc"])
    if (
        qc.get("schema") != "pfvarx.caller_qc/1" or qc.get("sample") != sample
        or not isinstance(qc.get("mixture"), dict) or not qc["mixture"].get("status")
        or not isinstance(qc.get("status_counts"), dict)
    ):
        raise ValueError("Caller QC JSON has an incompatible schema or sample")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in INPUTS:
        parser.add_argument("--" + name, required=True)
    for name in OUTPUTS:
        parser.add_argument("--out-" + name, required=True, type=Path)
    # Worker processes for the native caller only; the external contract is unchanged.
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    config = CallerConfig.model_validate(read_json(Path(args.frozen_config))["config"].get("caller", {}))
    check_backend(config)
    if config.backend == "native":
        importlib.import_module("pfvarx.caller.call").run(args)
        return
    for name in OUTPUTS:
        getattr(args, "out_" + name).parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".caller-", dir=args.out_vcf.parent) as directory:
        paths = {name: Path(directory) / f"sample.{name}" for name in OUTPUTS}
        command = [str(config.executable)]
        for name in INPUTS:
            command.extend(["--" + name, getattr(args, name.replace("-", "_"))])
        for name, path in paths.items():
            command.extend(["--out-" + name, str(path)])
        subprocess.run(command, check=True)
        validate_outputs(paths, args.sample, Path(args.reference))
        qc = read_json(paths["qc"])
        qc["adapter"] = {"contract": "pfvarx.layer2/1", "backend": "external", "sha256": config.sha256}
        write_json(paths["qc"], qc)
        for name, path in paths.items():
            path.replace(getattr(args, "out_" + name))


if __name__ == "__main__":
    main()
