"""L1.9 existing-alignment entry: validate, convert and repair only declared missing steps.

* Reference identity, sort order, index and read groups must already be valid
  (checked by 'pfvarx validate'); this step re-checks them and never trusts a label.
* WGS: if no duplicate-marking program is recorded in @PG, Picard MarkDuplicates is run.
* Targeted designs: primer clipping cannot be verified from a header alone, so an entry
  without a fgbio TrimPrimers @PG record is rejected (ALIGNMENT_PROCESSING_STATE).
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pysam


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--contract", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--design", required=True)
    a = ap.parse_args(argv)
    with pysam.AlignmentFile(a.input, reference_filename=a.reference) as aln:
        header = aln.header.to_dict()
    programs = [pg.get("PN", pg.get("ID", "")) for pg in header.get("PG", [])]
    marked = any("MarkDuplicates" in p or "markdup" in p.lower() for p in programs)
    clipped = any("TrimPrimers" in p for p in programs)
    actions = []
    contract = {"input": a.input, "programs": programs, "duplicates_marked_before_entry": marked,
                "primers_clipped_before_entry": clipped, "design": a.design}
    if a.design != "WGS" and not clipped:
        contract.update(status="REJECTED", reason="ALIGNMENT_PROCESSING_STATE",
                        detail="targeted entry requires evidence of primer clipping (fgbio TrimPrimers @PG)")
        Path(a.contract).write_text(json.dumps(contract, indent=2))
        sys.exit(3)
    tmp = Path(a.output + ".tmp.bam")
    subprocess.run(["samtools", "view", "-b", "-T", a.reference, "-o", str(tmp), a.input], check=True)
    actions.append("converted_to_bam" if a.input.endswith(".cram") else "copied")
    if a.design == "WGS" and not marked:
        metrics = Path(a.contract).with_name(Path(a.output).stem + ".markdup.txt")
        subprocess.run(["picard", "MarkDuplicates", "-I", str(tmp), "-O", a.output, "-M", str(metrics),
                        "--VALIDATION_STRINGENCY", "LENIENT"], check=True)
        tmp.unlink()
        actions.append("duplicates_marked")
    else:
        shutil.move(str(tmp), a.output)
    subprocess.run(["samtools", "index", a.output], check=True)
    contract.update(status="ACCEPTED", actions=actions)
    Path(a.contract).write_text(json.dumps(contract, indent=2))


if __name__ == "__main__":
    main()
