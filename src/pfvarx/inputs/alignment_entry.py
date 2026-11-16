"""L1.9 imported-alignment entry: re-check the file and do the steps its documented state needs.

The sample manifest documents each imported BAM/CRAM's ``processing_state``; validation has
already compared it with the reads (fgbio consensus and UMI tags) and the @PG records.
This step re-checks reference identity through samtools, converts CRAM to BAM, and then:

* ``wgs``: marks duplicates (Picard) when the state is ``aligned``, or when no state is
  documented and no duplicate-marking program is recorded in @PG (the schema-3 behaviour);
* ``targeted_umi`` with raw reads (``aligned``): adds mate tags (``samtools fixmate -m``) so
  the reads can be template-coordinate sorted and grouped by UMI; consensus follows in the
  workflow;
* otherwise the alignment is passed on unchanged. Primer trimming, when the state needs it,
  is a separate workflow step.

The contract JSON records the state, the steps planned and the actions taken.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pysam


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--contract", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--state", default="")
    ap.add_argument("--steps", default="")
    a = ap.parse_args(argv)
    steps = [s for s in a.steps.split(",") if s]
    with pysam.AlignmentFile(a.input, reference_filename=a.reference) as aln:
        header = aln.header.to_dict()
    programs = [pg.get("PN", pg.get("ID", "")) for pg in header.get("PG", [])]
    marked = any("MarkDuplicates" in p or "markdup" in p.lower() for p in programs)
    contract = {"input": a.input, "design": a.design, "processing_state": a.state or None,
                "planned_steps": steps, "programs": programs, "duplicates_marked_before_entry": marked,
                "primers_clipped_before_entry": any("TrimPrimers" in p for p in programs)}
    actions: list[str] = []
    out = Path(a.output)
    tmp = Path(a.output + ".tmp.bam")
    subprocess.run(["samtools", "view", "-b", "-T", a.reference, "-o", str(tmp), a.input], check=True)
    actions.append("converted_to_bam" if a.input.lower().endswith(".cram") else "copied")
    mark = "mark_duplicates" in steps or ("mark_duplicates_if_unrecorded" in steps and not marked)
    if a.design == "wgs" and mark:
        metrics = Path(a.contract).with_name(out.stem + ".markdup.txt")
        subprocess.run(["picard", "MarkDuplicates", "-I", str(tmp), "-O", str(out), "-M", str(metrics),
                        "--VALIDATION_STRINGENCY", "LENIENT"], check=True)
        tmp.unlink()
        actions.append("duplicates_marked")
    elif a.design == "targeted_umi" and "umi_consensus" in steps:
        # Template-coordinate sorting for GroupReadsByUmi needs mate CIGAR (MC) tags.
        named = Path(a.output + ".qn.bam")
        subprocess.run(["samtools", "sort", "-n", "-o", str(named), str(tmp)], check=True)
        fixmate = subprocess.Popen(["samtools", "fixmate", "-m", "-u", str(named), "-"], stdout=subprocess.PIPE)
        subprocess.run(["samtools", "sort", "-o", str(out), "-"], stdin=fixmate.stdout, check=True)
        if fixmate.stdout is not None:
            fixmate.stdout.close()
        if fixmate.wait() != 0:
            raise subprocess.CalledProcessError(fixmate.returncode, "samtools fixmate")
        tmp.unlink()
        named.unlink()
        actions.append("mate_tags_added_for_umi_grouping")
    else:
        shutil.move(str(tmp), out)
    subprocess.run(["samtools", "index", str(out)], check=True)
    contract.update(status="ACCEPTED", actions=actions)
    Path(a.contract).write_text(json.dumps(contract, indent=2))


if __name__ == "__main__":
    main()
