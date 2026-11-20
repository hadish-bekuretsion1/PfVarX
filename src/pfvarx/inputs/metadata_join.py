"""Join a study-metadata table into the sample sheet by ``sample_id`` (``./pfvarx join-metadata``).

Study teams usually keep population, location and dates in their own table. Typing them
into ``samples.tsv`` by hand is the most error-prone step before a run, so this module
merges them, shared by the CLI and the configuration GUI:

* rows are matched on ``sample_id``; every lane row of a sample gets the same values;
* only **empty** cells are filled. A cell that already holds a different value is a
  conflict: it is reported and left unchanged, never overwritten;
* input-file and identity columns (paths, lane, library) are never taken from the study
  table;
* sample IDs found in only one of the two files are reported;
* ``latitude``, ``longitude`` and ``collection_date`` are checked (decimal degrees in range,
  ``YYYY-MM-DD``) and problems are reported as warnings.

Without ``apply`` nothing is written. With it, the sample sheet is first copied to
``<name>.before-join-<UTC timestamp>``.
"""

from __future__ import annotations

import datetime as dt
import shutil
from pathlib import Path
from typing import Any

from pfvarx.inputs.content import content_problem
from pfvarx.io import read_tsv, write_tsv

KEY = "sample_id"
# Never filled from a study table: they describe the files and libraries, not the study.
PROTECTED = {"sample_id", "lane", "library_id", "fastq_1", "fastq_2", "bam", "vcf", "flowcell_type",
             "technique", "amplification_route", "panel_id", "panel_version", "library_protocol_id",
             "processing_state", "adapter_kit"}
NULLS = {"", "na", "n/a", "null", "none", "."}


def _blank(value: str | None) -> bool:
    return value is None or value.strip().lower() in NULLS


def check_value(column: str, value: str) -> str | None:
    """Why a study value is not usable by the analyses, or None."""
    if _blank(value):
        return None
    if column in ("latitude", "longitude"):
        try:
            number = float(value)
        except ValueError:
            return f"{column} {value!r} is not a decimal number (decimal degrees)"
        limit = 90 if column == "latitude" else 180
        if not -limit <= number <= limit:
            return f"{column} {value} is outside -{limit} to {limit}"
    if column == "collection_date":
        try:
            dt.date.fromisoformat(value)
        except ValueError:
            return f"collection_date {value!r} is not YYYY-MM-DD"
    return None


def join_metadata(root: Path, sheet: Path, table: Path, apply: bool = False) -> dict[str, Any]:
    """Report (and with *apply*, write) the join of *table* into *sheet*."""
    problem = content_problem(table, ("SAMPLE_SHEET", "TSV"))
    if problem:
        raise ValueError(f"{problem}. The study table must be tab-separated with a {KEY} column")
    study = read_tsv(table)
    if not study or KEY not in study[0]:
        raise ValueError(f"{table.name}: the study table needs a header with a {KEY} column")
    rows = read_tsv(sheet) if sheet.is_file() else []
    if not rows:
        raise ValueError(f"{sheet} has no samples yet; register the input files first")
    by_id: dict[str, dict[str, str]] = {}
    duplicates: list[str] = []
    for record in study:
        sid = record.get(KEY, "").strip()
        if not sid:
            continue
        if sid in by_id and by_id[sid] != record:
            duplicates.append(sid)
        by_id[sid] = record
    if duplicates:
        raise ValueError(f"{table.name}: {KEY} values with differing rows: {sorted(set(duplicates))}; "
                         "keep one row per sample")
    columns = [c for c in study[0] if c != KEY]
    skipped = sorted(c for c in columns if c in PROTECTED)
    joined = [c for c in columns if c not in PROTECTED]
    sheet_ids = {r[KEY] for r in rows}
    existing = list(rows[0])  # sheet columns before the join
    filled = 0
    conflicts: list[dict[str, str]] = []
    warnings: list[str] = []
    for row in rows:
        match = by_id.get(row[KEY])
        if match is None:
            continue
        for column in joined:
            value = (match.get(column) or "").strip()
            if _blank(value):
                continue
            current = (row.get(column) or "").strip()
            if _blank(current):
                row[column] = value
                filled += 1
            elif current != value:
                conflicts.append({"sample_id": row[KEY], "lane": row.get("lane", ""), "column": column,
                                  "sample_sheet": current, "study_table": value})
    for sid, record in sorted(by_id.items()):
        if sid in sheet_ids:
            for column in joined:
                why = check_value(column, record.get(column) or "")
                if why:
                    warnings.append(f"sample {sid}: {why}")
    report: dict[str, Any] = {
        "status": "APPLIED" if apply else "PLANNED",
        "sample_sheet": str(sheet.relative_to(root)) if sheet.is_relative_to(root) else str(sheet),
        "study_table": str(table),
        "columns_joined": joined, "columns_skipped_protected": skipped,
        "new_columns": [c for c in joined if c not in existing],
        "cells_filled": filled,
        "conflicts_not_overwritten": conflicts,
        "samples_without_study_row": sorted(sheet_ids - set(by_id)),
        "study_rows_without_sample": sorted(set(by_id) - sheet_ids),
        "warnings": warnings,
    }
    if apply and filled:
        stamp = dt.datetime.now(dt.UTC).strftime("%Y%m%dT%H%M%SZ")
        backup = sheet.with_name(f"{sheet.name}.before-join-{stamp}")
        shutil.copy2(sheet, backup)
        order = existing + [c for c in joined if c not in existing]
        write_tsv(sheet, rows, order)
        report["backup"] = str(backup.relative_to(root)) if backup.is_relative_to(root) else str(backup)
    elif apply:
        report["status"] = "NOTHING_TO_WRITE"
    return report


def join_for_workspace(root: Path, table: str, apply: bool = False) -> dict[str, Any]:
    """Join into the configured sample sheet; *table* may be absolute or workspace-relative."""
    from pfvarx.configuration.loader import load_bundle

    cfg = load_bundle(root).config
    path = Path(table).expanduser()
    path = path if path.is_absolute() else root / path
    if not path.is_file():
        raise ValueError(f"Study table not found: {table}")
    return join_metadata(root, root / cfg.samples_tsv, path, apply)
