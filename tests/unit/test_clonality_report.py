"""Per-sample clonality report: monoclonal, or polyclonal with clone count and percentages."""

import json
from pathlib import Path

from pfvarx.reporting.clonality import classify, sentence, write_report


def _mix(k: int, weights: list[float], at_limit: bool = False) -> dict:
    tests = [{"k": i, "p_value": 0.025} for i in range(1, k)]
    # As PfVarX-Call writes it: the note sits on the largest K that was fitted.
    tests.append({"k": k, "p_value": None, "note": "K_MAX reached; larger K not tested"} if at_limit
                 else {"k": k, "p_value": 0.6})
    return {"status": "ESTIMATED", "k_selected": k, "fits": {str(k): {"weights": weights}}, "tests": tests,
            "min_component_weight": 0.02}


def test_monoclonal_and_polyclonal_sentences() -> None:
    mono = classify("X", _mix(1, [1.0]), None)
    assert mono["clonality"] == "MONOCLONAL" and sentence(mono) == "Sample X is monoclonal"
    poly = classify("Y", _mix(3, [0.2, 0.5, 0.3]), "POLYCLONAL")
    assert poly["number_of_clones"] == "3"
    assert poly["clone_proportions"] == "Clone 1: 50.0%; Clone 2: 30.0%; Clone 3: 20.0%"  # largest first
    assert sentence(poly).startswith("Sample Y is polyclonal with 3 clones: Clone 1: 50.0%, Clone 2: 30.0%, Clone 3: 20.0%")


def test_limits_and_disagreements_are_stated() -> None:
    limit = classify("Z", _mix(3, [0.6, 0.3, 0.1], at_limit=True), None)
    assert limit["number_of_clones"] == "at least 3" and "more may be present" in limit["note"]
    hidden = classify("H", _mix(1, [1.0]), "POLYCLONAL")
    assert hidden["clonality"] == "POLYCLONAL" and hidden["number_of_clones"] == "not resolved"
    clash = classify("C", _mix(2, [0.9, 0.1]), "MONOCLONAL")
    assert clash["clonality"] == "INDETERMINATE" and "disagree" not in clash["note"] and "monoclonal" in clash["note"]
    missing = classify("M", None, None)
    assert missing["clonality"] == "NOT_DETERMINED" and "no caller mixture" in missing["note"]
    da001_only = classify("D", {"status": "INSUFFICIENT_INFORMATIVE_SITES"}, "MONOCLONAL")
    assert da001_only["clonality"] == "MONOCLONAL" and da001_only["basis"] == "DA001 (MOIRE)"


def test_report_files(tmp_path: Path) -> None:
    results = tmp_path / "results"
    for sample, mix in (("S1", _mix(1, [1.0])), ("S2", _mix(2, [0.1014, 0.8986]))):
        qc = results / "qc/samples" / sample
        qc.mkdir(parents=True)
        (qc / f"{sample}.caller_qc.json").write_text(json.dumps({"mixture": mix}))
    rows = write_report(results, ["S1", "S2"], "RUN1")
    assert [r["clonality"] for r in rows] == ["MONOCLONAL", "POLYCLONAL"]
    text = (results / "report/clonality_report.txt").read_text()
    assert "Sample S1 is monoclonal" in text
    assert "Sample S2 is polyclonal with 2 clones: Clone 1: 89.9%, Clone 2: 10.1%" in text
    assert (results / "report/clonality_report.tsv").read_text().startswith("sample_id\tclonality\tnumber_of_clones")


def test_da001_calls_from_the_caller_mixture_are_not_counted_twice(tmp_path) -> None:
    from pfvarx.reporting.clonality import _da001

    d = tmp_path / "downstream" / "DA001"
    d.mkdir(parents=True)
    (d / "clonality.tsv").write_text("sample_id\tclonality\tbasis\nA\tPOLYCLONAL\tMOIRE\n"
                                     "B\tPOLYCLONAL\tCALLER_MIXTURE_LRT\n")
    assert _da001(tmp_path) == {"A": "POLYCLONAL"}
    (d / "clonality.tsv").write_text("sample_id\tclonality\nC\tMONOCLONAL\n")  # runs before the basis column
    assert _da001(tmp_path) == {"C": "MONOCLONAL"}
