"""Study metadata joins into the sample sheet by sample_id; names get advisory warnings."""

from __future__ import annotations

from pathlib import Path

import pytest

from pfvarx.constants import InputType
from pfvarx.inputs.metadata_join import join_metadata
from pfvarx.inputs.validate import naming_problems
from pfvarx.io import read_tsv

SHEET = ("sample_id\tlane\tfastq_1\tfastq_2\tlibrary_id\tpopulation\n"
         "ETH001\tL001\tinputs/fastq/ETH001_L001_R1.fastq.gz\tinputs/fastq/ETH001_L001_R2.fastq.gz\tETH001_lib1\t\n"
         "ETH001\tL002\tinputs/fastq/ETH001_L002_R1.fastq.gz\tinputs/fastq/ETH001_L002_R2.fastq.gz\tETH001_lib1\t\n"
         "ETH002\tL001\tinputs/fastq/ETH002_L001_R1.fastq.gz\tinputs/fastq/ETH002_L001_R2.fastq.gz\tETH002_lib1\tBahirDar\n"
         "ETH003\tL001\tinputs/fastq/ETH003_L001_R1.fastq.gz\tinputs/fastq/ETH003_L001_R2.fastq.gz\tETH003_lib1\t\n")
STUDY = ("sample_id\tpopulation\tcountry\tlatitude\tlongitude\tcollection_date\tfastq_1\tage_group\n"
         "ETH001\tGondar\tEthiopia\t12.600\t37.467\t2025-07-18\t/elsewhere/x.fq.gz\tAdult\n"
         "ETH002\tGondar\tEthiopia\t11.59\t37.39\t18/07/2025\t\tChild\n"
         "ETH999\tGondar\tEthiopia\t12.6\t37.4\t2025-07-19\t\tAdult\n")


def files(tmp_path: Path) -> tuple[Path, Path]:
    sheet, study = tmp_path / "samples.tsv", tmp_path / "study.tsv"
    sheet.write_text(SHEET)
    study.write_text(STUDY)
    return sheet, study


def test_preview_changes_nothing_and_reports_everything(tmp_path: Path) -> None:
    sheet, study = files(tmp_path)
    report = join_metadata(tmp_path, sheet, study)
    assert report["status"] == "PLANNED" and sheet.read_text() == SHEET
    assert report["columns_joined"] == ["population", "country", "latitude", "longitude", "collection_date", "age_group"]
    assert report["new_columns"] == ["country", "latitude", "longitude", "collection_date", "age_group"]
    assert report["columns_skipped_protected"] == ["fastq_1"]  # file paths never come from a study table
    assert report["samples_without_study_row"] == ["ETH003"]
    assert report["study_rows_without_sample"] == ["ETH999"]
    assert report["conflicts_not_overwritten"] == [{"sample_id": "ETH002", "lane": "L001", "column": "population",
                                                    "sample_sheet": "BahirDar", "study_table": "Gondar"}]
    assert any("ETH002" in w and "YYYY-MM-DD" in w for w in report["warnings"])


def test_apply_fills_only_empty_cells_on_every_lane_and_keeps_a_backup(tmp_path: Path) -> None:
    sheet, study = files(tmp_path)
    report = join_metadata(tmp_path, sheet, study, apply=True)
    assert report["status"] == "APPLIED" and (tmp_path / report["backup"]).read_text() == SHEET
    rows = {(r["sample_id"], r["lane"]): r for r in read_tsv(sheet)}
    assert rows[("ETH001", "L001")]["latitude"] == rows[("ETH001", "L002")]["latitude"] == "12.600"
    assert rows[("ETH001", "L001")]["fastq_1"] == "inputs/fastq/ETH001_L001_R1.fastq.gz"
    assert rows[("ETH002", "L001")]["population"] == "BahirDar"  # conflict kept, not overwritten
    assert rows[("ETH003", "L001")]["country"] == ""
    assert join_metadata(tmp_path, sheet, study, apply=True)["status"] == "NOTHING_TO_WRITE"


def test_study_table_needs_one_row_per_sample(tmp_path: Path) -> None:
    sheet, study = files(tmp_path)
    study.write_text(STUDY + "ETH001\tAddis\tEthiopia\t9.0\t38.7\t2025-07-18\t\tAdult\n")
    with pytest.raises(ValueError, match="differing rows"):
        join_metadata(tmp_path, sheet, study)
    study.write_text("country\tlatitude\nEthiopia\t12\n")
    with pytest.raises(ValueError, match="sample_id"):
        join_metadata(tmp_path, sheet, study)


def test_names_outside_the_standard_are_warnings() -> None:
    good = {"sample_id": "ETH001", "lane": "L001", "fastq_1": "inputs/fastq/ETH001_L001_R1.fastq.gz",
            "fastq_2": "inputs/fastq/ETH001_L001_R2.fastq.gz", "bam": "inputs/bam/ETH001.bam"}
    assert naming_problems([good], InputType.FASTQ) == []
    assert naming_problems([good], InputType.BAM) == []
    illumina = {**good, "fastq_1": "x/ETH001_S1_L001_R1_001.fastq.gz", "fastq_2": "x/ETH001_S1_L001_R2_001.fastq.gz"}
    assert naming_problems([illumina], InputType.FASTQ) == []  # starts with the ID and names the mate
    odd = {**good, "lane": "1", "fastq_1": "reads_a.fq.gz", "bam": "aln.bam"}
    found = naming_problems([odd], InputType.FASTQ) + naming_problems([odd], InputType.BAM)
    assert any("lane '1'" in f for f in found) and any("fastq_1 'reads_a.fq.gz'" in f for f in found)
    assert any("ETH001.bam" in f for f in found)
