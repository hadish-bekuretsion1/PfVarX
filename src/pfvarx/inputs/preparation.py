"""Inputs Preparation, first part: the sample metadata and the optional study metadata merge.

``./pfvarx inputs-preparation`` (GUI: **Inputs Preparation**) runs, in order:

1. **Sample metadata.** When the configured sheet (``samples_tsv``, by default
   ``prepared/samples.tsv``) does not exist, it is built from the input files: paired FASTQ
   names give one row per sample and lane. Its folder is created here, when it is needed. For
   BAM/CRAM input a blank template is written and preparation stops, because sample
   identities cannot be read from alignment file names. An existing sheet is used as it is
   and never overwritten.
2. **Study metadata (Study Metadata: Yes).** With ``study_metadata: false`` nothing is merged
   and preparation continues with the sample metadata alone. With ``true``, the study table in
   ``study_metadata_dir`` (tab-separated, with a ``sample_id`` column) is merged into the
   sample metadata: empty cells are filled, existing values are never overwritten and
   differences are listed. If the folder holds no study table yet, preparation continues with
   the sample metadata alone and says so; with several, it stops and asks which one to use
   (``--metadata-table``).
3. Reference preparation, verification and readiness follow as separate stages
   (``prepare``, ``verify-prepared``, ``ready``). Every folder they write is created when
   it is first needed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pfvarx.configuration.loader import load_bundle
from pfvarx.initialization.layout import SAMPLE_COLUMNS
from pfvarx.inputs.content import detect
from pfvarx.inputs.metadata_join import join_metadata
from pfvarx.inputs.organize import discover
from pfvarx.io import write_json, write_tsv

REPORT = "metadata/inputs_preparation.json"


def _relative(root: Path, path: Path) -> str:
    return str(path.relative_to(root)) if path.is_relative_to(root) else str(path)


def study_tables(folder: Path, exclude: set[Path]) -> list[Path]:
    """Tab-separated tables with a sample_id column in *folder*, except *exclude*."""
    if not folder.is_dir():
        return []
    found = []
    for path in sorted(folder.iterdir()):
        if path.name.startswith(".") or not path.is_file() or path.resolve() in exclude:
            continue
        if ".before-join-" in path.name:  # backups written by an earlier merge
            continue
        if detect(path) == "SAMPLE_SHEET":  # TSV whose header has sample_id
            found.append(path)
    return found


def build_sample_sheet(root: Path) -> dict[str, Any]:
    """Step 1: use the existing sample metadata, or build it from the input files."""
    cfg = load_bundle(root).config
    sheet = root / cfg.samples_tsv
    if sheet.exists() or sheet.is_symlink():
        return {"status": "EXISTING", "sample_sheet": _relative(root, sheet),
                "detail": "The existing sample metadata is used as it is (never overwritten)."}
    if cfg.input_type.value == "FASTQ":
        found = discover(root, root / "inputs/fastq", sheet)
        return {"status": "CREATED", "sample_sheet": _relative(root, sheet), "samples": found["samples"],
                "lanes": found["lanes"], "ignored_files": found["ignored_files"],
                "detail": "Built from the paired FASTQ names in inputs/fastq/. Sample and lane IDs come from the "
                          "file names and library_id defaults to sample_id: review them."}
    write_tsv(sheet, [], SAMPLE_COLUMNS)
    raise ValueError(f"Created a blank sample metadata sheet at {_relative(root, sheet)}: {cfg.input_type.value} "
                     "file names do not give sample identities. Enter one row per alignment (sample_id, lane, bam, "
                     "library_id), then run Inputs Preparation again.")


def merge_study_metadata(root: Path, table: str | None = None) -> dict[str, Any]:
    """Step 2: merge the study metadata into the sample metadata when Study Metadata is Yes."""
    cfg = load_bundle(root).config
    sheet = root / cfg.samples_tsv
    if table is not None:
        path = Path(table).expanduser()
        path = path if path.is_absolute() else root / path
        if not path.is_file():
            raise ValueError(f"Study metadata table not found: {table}")
        return {**join_metadata(root, sheet, path, apply=True), "merge": "MERGED_CHOSEN_TABLE"}
    if not cfg.study_metadata:
        return {"merge": "NOT_USED", "detail": "Study Metadata: No. Continued with the sample metadata built from "
                "the input files only."}
    folder = root / cfg.study_metadata_dir
    own = {sheet.resolve()} | {(root / p).resolve() for p in (cfg.assay.library_catalog,
                                                              cfg.analyses.resistance_catalog) if p is not None}
    tables = study_tables(folder, own)
    if not tables:
        return {"merge": "NO_STUDY_TABLE", "warning": f"Study Metadata: Yes, but {_relative(root, folder)}/ holds no "
                "study table with a sample_id column yet. Continued with the sample metadata only; add the table "
                "and run Inputs Preparation again to merge it."}
    if len(tables) > 1:
        names = ", ".join(_relative(root, t) for t in tables)
        raise ValueError(f"Several study metadata tables in {_relative(root, folder)}/: {names}. Keep only the one "
                         "to merge in that folder, or choose it with ./pfvarx inputs-preparation --metadata-table PATH.")
    return {**join_metadata(root, sheet, tables[0], apply=True), "merge": "MERGED"}


def prepare_samples(root: Path, table: str | None = None) -> dict[str, Any]:
    """Steps 1 and 2 of Inputs Preparation; the result is also written to metadata/inputs_preparation.json."""
    result = {"sample_sheet": build_sample_sheet(root), "study_metadata": merge_study_metadata(root, table)}
    write_json(root / REPORT, result)
    return result
