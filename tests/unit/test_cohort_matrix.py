"""The cohort evidence matrix: identical values to the mapping it replaces, at constant cost.

Layer 5 holds every sample's evidence at every cohort site. One dict per cell costs about 474
bytes (measured), so a few hundred samples over tens of thousands of sites runs into gigabytes
before any analysis starts. The columns in ``common.SiteRow`` cost about 30 bytes per cell.

What must not change is the numbers: ``wsaf`` is compared against reporting floors (0.05, 0.95),
so it is stored as a float64 exactly as the file holds it. These tests pin that, and the
reading interface every analysis and test fixture uses.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from pfvarx.downstream.common import Cell, SiteRow, load_context, load_matrix

STATES = ("ALT_DOMINANT", "MIXED", "REF_DOMINANT", "MISSING_LOW_DEPTH", "UNCALLABLE")
# Values that sit exactly on the reporting floors, where a rounding change would flip a call.
EDGE_WSAF = ("0.05", "0.95", "0.0500", "0.949999", "0.0499999", "1.0", "0.0", "0.3333333333")


def _dict_cell(row: dict[str, str]) -> dict[str, object]:
    """Exactly the cell the previous implementation built, as the reference."""
    return {"state": row["state"], "ref": int(row["ref_count"]), "alt": int(row["alt_count"]),
            "depth": int(row["depth"]), "wsaf": float(row["wsaf"]) if row["wsaf"] else None,
            "caller_status": row["caller_status"]}


def _write(run: Path, samples: list[str], rows: list[dict[str, str]]) -> None:
    cohort = run / "results" / "cohort"
    cohort.mkdir(parents=True, exist_ok=True)
    (cohort / "cohort_samples.txt").write_text("\n".join(samples) + "\n")
    columns = ["site_id", "sample_id", "state", "ref_count", "alt_count", "depth", "wsaf", "caller_status"]
    with gzip.open(cohort / "evidence_matrix.tsv.gz", "wt") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row[c] for c in columns) + "\n")


def _context(tmp_path: Path, samples: list[str], rows: list[dict[str, str]]):
    run = tmp_path / "run"
    _write(run, samples, rows)
    (run / "frozen_config.json").write_text(json.dumps(
        {"config": {"samples_tsv": "samples.tsv"}, "parameters": {}}))
    (run / "frozen_plan.json").write_text(json.dumps({"samples": [{"sample_id": s} for s in samples]}))
    (tmp_path / "samples.tsv").write_text("sample_id\n" + "\n".join(samples) + "\n")
    ctx = load_context(tmp_path, run, with_matrix=False)
    skipped = load_matrix(ctx, run / "results" / "cohort" / "evidence_matrix.tsv.gz")
    return ctx, skipped


def _rows(samples: list[str], sites: int) -> list[dict[str, str]]:
    out = []
    for i in range(sites):
        for j, sample in enumerate(samples):
            out.append({"site_id": f"ctg_{i}", "sample_id": sample, "state": STATES[(i + j) % len(STATES)],
                        "ref_count": str(10 + j), "alt_count": str(j), "depth": str(10 + 2 * j),
                        "wsaf": EDGE_WSAF[(i + j) % len(EDGE_WSAF)], "caller_status": "PASS"})
    return out


def test_every_cell_matches_the_mapping_it_replaces(tmp_path: Path) -> None:
    samples = [f"S{i:03d}" for i in range(7)]
    rows = _rows(samples, 11)
    ctx, skipped = _context(tmp_path, samples, rows)
    assert skipped == {"unknown_sample": 0, "unknown_state": 0}
    reference = {(r["site_id"], r["sample_id"]): _dict_cell(r) for r in rows}
    assert len(reference) == 77
    for (site, sample), expected in reference.items():
        cell = ctx.matrix[site].get(sample)
        assert cell is not None, (site, sample)
        for key, value in expected.items():
            assert cell[key] == value, (site, sample, key, cell[key], value)
            # float64 round trip: the stored wsaf is the same number the file holds
            if key == "wsaf" and value is not None:
                assert repr(cell[key]) == repr(value)


def test_the_reading_interface_is_unchanged(tmp_path: Path) -> None:
    samples = ["A", "B"]
    rows = _rows(samples, 2)
    ctx, _ = _context(tmp_path, samples, rows)
    row = ctx.matrix["ctg_0"]
    assert isinstance(row, SiteRow)
    assert row.get("A")["state"] == rows[0]["state"]          # mapping-style read
    assert row["A"]["ref"] == int(rows[0]["ref_count"])        # subscript
    assert row.get("ABSENT") is None and "ABSENT" not in row   # missing sample, as a dict gave
    with pytest.raises(KeyError):
        row["ABSENT"]
    assert sorted(row) == ["A", "B"] and len(row) == 2          # iteration and length
    cell = row.get("A")
    assert cell.get("state") == cell["state"] and cell.get("nope", "fallback") == "fallback"
    assert "wsaf=" in repr(cell)


def test_a_missing_record_is_absent_not_a_zero(tmp_path: Path) -> None:
    samples = ["A", "B", "C"]
    rows = [r for r in _rows(samples, 1) if r["sample_id"] != "B"]
    ctx, _ = _context(tmp_path, samples, rows)
    row = ctx.matrix["ctg_0"]
    assert row.get("B") is None, "a sample with no record must not read as depth 0"
    assert sorted(row) == ["A", "C"]
    assert row.get("A")["depth"] == int(rows[0]["depth"])


def test_an_empty_wsaf_stays_none_and_is_never_compared_as_a_number(tmp_path: Path) -> None:
    samples = ["A"]
    rows = _rows(samples, 1)
    rows[0]["wsaf"] = ""
    ctx, _ = _context(tmp_path, samples, rows)
    assert ctx.matrix["ctg_0"].get("A")["wsaf"] is None


def test_rows_outside_the_cohort_or_vocabulary_are_skipped_and_counted(tmp_path: Path) -> None:
    samples = ["A"]
    rows = _rows(samples, 1)
    rows.append({**rows[0], "sample_id": "NOT_IN_COHORT"})
    rows.append({**rows[0], "site_id": "ctg_1", "state": "SOMETHING_NEW"})
    ctx, skipped = _context(tmp_path, samples, rows)
    assert skipped == {"unknown_sample": 1, "unknown_state": 1}
    assert "ctg_1" not in ctx.matrix


def test_each_cell_costs_far_less_than_a_dict(tmp_path: Path) -> None:
    """A regression guard on the memory, measured on a small matrix (~30 B/cell against ~474)."""
    import tracemalloc

    samples = [f"S{i:03d}" for i in range(40)]
    rows = _rows(samples, 300)
    run = tmp_path / "run"
    _write(run, samples, rows)
    (run / "frozen_config.json").write_text(json.dumps(
        {"config": {"samples_tsv": "samples.tsv"}, "parameters": {}}))
    (run / "frozen_plan.json").write_text(json.dumps({"samples": []}))
    (tmp_path / "samples.tsv").write_text("sample_id\n" + "\n".join(samples) + "\n")
    ctx = load_context(tmp_path, run, with_matrix=False)
    tracemalloc.start()
    load_matrix(ctx, run / "results" / "cohort" / "evidence_matrix.tsv.gz")
    held, _ = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    per_cell = held / (len(samples) * 300)
    assert per_cell < 120, f"{per_cell:.0f} bytes per cell: the columns regressed towards per-cell dicts"


def test_cell_is_constructed_directly_for_fixtures() -> None:
    cell = Cell("MIXED", 7, 3, 10, 0.3, "PASS")
    assert cell["state"] == "MIXED" and cell["ref"] == 7 and cell["wsaf"] == 0.3
