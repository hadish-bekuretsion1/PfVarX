"""assay_design selects the branch; the manifest and library catalog describe each library."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml

from pfvarx.configuration.library import CATALOG_SCHEMA, LibraryCatalog, umi_settings
from pfvarx.configuration.loader import ConfigError, load_parameters, validate_bundle
from pfvarx.configuration.migrate import migrate_workspace, plan_v3
from pfvarx.constants import Design, InputType
from pfvarx.inputs import assay
from pfvarx.inputs.content import detect

REPO = Path(__file__).resolve().parents[2]
PARAMS = load_parameters(REPO)
CATALOG = LibraryCatalog.model_validate({
    "schema_version": CATALOG_SCHEMA,
    "panels": [{"panel_id": "P", "panel_version": "1", "targets_bed": "t.bed", "primers_tsv": "p.tsv"},
               {"panel_id": "P", "panel_version": "2", "targets_bed": "t.bed", "primers_tsv": "p.tsv"},
               {"panel_id": "CAP", "panel_version": "1", "targets_bed": "t.bed", "probes_bed": "probes.bed"}],
    "library_protocols": [
        {"library_protocol_id": "UMI8", "technique": "TBS_UMI",
         "umi": {"read_structure": "8M+T 8M+T", "filter_min_reads": 2}},
        {"library_protocol_id": "SWGA_KIT", "technique": "SELECTIVE_WGS"},
    ],
})


def row(sample: str = "S1", **cells: str) -> dict[str, str]:
    return {"sample_id": sample, "lane": "L1", "library_id": f"{sample}_lib", **cells}


def resolve(design: Design, rows: list[dict[str, str]], input_type: InputType = InputType.FASTQ
            ) -> tuple[list[str], dict[str, Any], Any]:
    rep, facts, panel = assay.resolve(design, input_type, rows, CATALOG, PARAMS)
    return [f.detail for f in rep.findings if f.level == "FAIL"], facts, panel


@pytest.mark.parametrize(("design", "label", "ok"), [
    (Design.WGS, "WGS", True), (Design.WGS, "SELECTIVE_WGS", True), (Design.WGS, "AMPLICON", False),
    (Design.TARGETED_NO_UMI, "AMPLICON", True), (Design.TARGETED_UMI, "AMPLICON", True),
    (Design.TARGETED_UMI, "TBS_UMI", True), (Design.TARGETED_NO_UMI, "TBS_UMI", False),
    (Design.TARGETED_NO_UMI, "TBS_NON_UMI", True), (Design.TARGETED_UMI, "TBS_NON_UMI", False),
    (Design.TARGETED_NO_UMI, "WGS", False), (Design.WGS, "TBS_UMI", False),
])
def test_technique_label_must_match_the_design(design: Design, label: str, ok: bool) -> None:
    # Targeted UMI cases import a documented consensus BAM, which needs no UMI protocol.
    cells = {} if design == Design.WGS else {"panel_id": "P", "panel_version": "1"}
    input_type = InputType.FASTQ
    if design == Design.TARGETED_UMI:
        cells["processing_state"] = "consensus"
        input_type = InputType.BAM
    failures, _, _ = resolve(design, [row(technique=label, **cells)], input_type)
    conflicts = [f for f in failures if "technique" in f and label in f]
    assert (not conflicts) == ok, failures
    if not ok:
        assert conflicts[0].startswith("sample S1 (library S1_lib)")


@pytest.mark.parametrize("empty", ["", "null", "NULL", "NA", ".", "None"])
def test_empty_or_null_technique_is_unspecified(empty: str) -> None:
    failures, facts, _ = resolve(Design.WGS, [row(technique=empty)])
    assert failures == [] and facts["S1"].technique is None and facts["S1"].amplification_route is None


def test_routes_are_recorded_from_the_label_protocol_or_panel() -> None:
    _, facts, _ = resolve(Design.WGS, [row(technique="SELECTIVE_WGS")])
    assert (facts["S1"].amplification_route, facts["S1"].route_source) == ("SWGA", "technique SELECTIVE_WGS")
    _, facts, _ = resolve(Design.WGS, [row(amplification_route="DIRECT")])
    assert facts["S1"].amplification_route == "DIRECT"  # confirmed direct WGS
    _, facts, _ = resolve(Design.WGS, [row(technique="WGS")])
    assert facts["S1"].amplification_route is None  # a WGS label alone does not confirm DIRECT
    _, facts, _ = resolve(Design.WGS, [row(library_protocol_id="SWGA_KIT")])
    assert facts["S1"].amplification_route == "SWGA" and facts["S1"].technique == "SELECTIVE_WGS"
    failures, facts, panel = resolve(Design.TARGETED_NO_UMI, [row(panel_id="P", panel_version="1")])
    assert failures == [] and facts["S1"].amplification_route == "AMPLICON" and panel.key == ("P", "1")


def test_route_conflicts_name_the_sample_and_library() -> None:
    failures, _, _ = resolve(Design.WGS, [row(technique="SELECTIVE_WGS", amplification_route="DIRECT")])
    assert any(f.startswith("sample S1 (library S1_lib)") and "disagree" in f for f in failures)
    failures, _, _ = resolve(Design.WGS, [row(amplification_route="AMPLICON")])
    assert any("does not fit assay_design wgs" in f for f in failures)
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row(panel_id="P", panel_version="1", amplification_route="SWGA")])
    assert any("SWGA" in f for f in failures)
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row(panel_id="CAP", panel_version="1", amplification_route="AMPLICON")])
    assert any("defines no primers" in f for f in failures)


def test_panel_and_protocol_are_resolved_from_the_catalog() -> None:
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row()])
    assert any("need panel_id and panel_version" in f for f in failures)
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row(panel_id="P", panel_version="9")])
    assert any("panel P version 9 is not in the library catalog" in f for f in failures)
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row("S1", panel_id="P", panel_version="1"),
                                                      row("S2", panel_id="P", panel_version="2")])
    assert any(f.startswith("one panel per run") for f in failures)
    failures, facts, _ = resolve(Design.TARGETED_UMI, [row(panel_id="P", panel_version="1", library_protocol_id="UMI8")])
    assert failures == []
    umi = facts["S1"].umi
    assert umi["read_structure"]["value"] == "8M+T 8M+T"
    assert umi["filter_min_reads"] == {"value": 2, "source": "LIBRARY_PROTOCOL"}
    assert umi["consensus_min_reads"]["source"] == "PARAMETER_DEFAULT"
    assert facts["S1"].count_unit == "consensus_molecules"
    failures, _, _ = resolve(Design.TARGETED_UMI, [row(panel_id="P", panel_version="1")])
    assert any("umi.read_structure" in f for f in failures)  # FASTQ needs UMI extraction settings
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row(panel_id="P", panel_version="1", library_protocol_id="UMI8")])
    assert any("defines UMIs" in f for f in failures)


def test_imported_alignments_use_the_documented_state() -> None:
    steps = {(d, s): assay.entry_steps(d, s) for d in Design for s in assay.DESIGN_STATES[d]}
    assert steps[(Design.WGS, "aligned")] == ["mark_duplicates"]
    assert steps[(Design.WGS, "duplicates_marked")] == []
    assert steps[(Design.TARGETED_NO_UMI, "aligned")] == ["trim_primers"]
    assert steps[(Design.TARGETED_NO_UMI, "primers_trimmed")] == []
    assert steps[(Design.TARGETED_UMI, "aligned")] == ["umi_consensus", "trim_primers"]
    assert steps[(Design.TARGETED_UMI, "consensus")] == ["trim_primers"]  # consensus is not generated again
    assert steps[(Design.TARGETED_UMI, "consensus_primers_trimmed")] == []
    assert assay.entry_steps(Design.WGS, None) == ["mark_duplicates_if_unrecorded"]  # schema-3 behaviour
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row(panel_id="P", panel_version="1")], InputType.BAM)
    assert any("need processing_state" in f for f in failures)
    failures, _, _ = resolve(Design.TARGETED_NO_UMI, [row(panel_id="P", panel_version="1",
                                                          processing_state="consensus")], InputType.BAM)
    assert any("not valid for assay_design targeted_no_umi" in f for f in failures)


def evidence(**kw: Any) -> dict[str, Any]:
    return {"programs": [], "reads_inspected": 100, "consensus_tagged": 0, "umi_tagged": 0,
            "trim_primers_recorded": False, "consensus_recorded": False, "duplicates_marked_recorded": False, **kw}


@pytest.mark.parametrize(("design", "state", "seen", "fails"), [
    (Design.TARGETED_UMI, "consensus", evidence(), True),  # documented consensus, no cD tags
    (Design.TARGETED_UMI, "consensus", evidence(consensus_tagged=90), False),
    (Design.TARGETED_UMI, "aligned", evidence(consensus_tagged=90, umi_tagged=90), True),  # already consensus
    (Design.TARGETED_UMI, "aligned", evidence(), True),  # raw reads without RX
    (Design.TARGETED_UMI, "aligned", evidence(umi_tagged=100), False),
    (Design.TARGETED_NO_UMI, "aligned", evidence(trim_primers_recorded=True), True),  # would trim twice
    (Design.TARGETED_NO_UMI, "primers_trimmed", evidence(), False),  # warning only
    (Design.WGS, "aligned", evidence(consensus_tagged=5), True),
])
def test_documented_state_is_checked_against_the_reads(design: Design, state: str, seen: dict[str, Any],
                                                       fails: bool) -> None:
    from pfvarx.inputs.validate import FileReport

    fact = assay.LibraryFacts(sample_id="S1", library_ids=["L"], design=design.value, processing_state=state,
                              entry_steps=assay.entry_steps(design, assay.ProcessingState(state)))
    rep = FileReport(path="x.bam", kind="ALIGNMENT")
    assay.check_alignment_state(fact, seen, rep)
    assert (rep.status == "FAIL") == fails, rep.findings


def schema3(design: str, input_type: str = "FASTQ", **assay_cfg: Any) -> dict[str, Any]:
    cfg = yaml.safe_load((REPO / "config/templates/analysis.yaml").read_text())
    cfg["schema_version"] = "pfvarx.config/3"
    cfg["samples_tsv"] = "inputs/metadata/samples.tsv"  # where schema-3 workspaces kept the sheet
    cfg["input_type"] = input_type
    cfg.pop("assay_design")
    cfg["assay"] = {"design": design, "sequencing_type": None, "amplification_route": "DIRECT",
                    "read_structure": None, "targets_bed": None, "primers_tsv": None,
                    "two_colour_chemistry": False, "platform": "ILLUMINA", **assay_cfg}
    return cfg


def test_schema3_is_planned_into_assay_design_catalog_and_manifest_columns() -> None:
    new, catalog, defaults = plan_v3(schema3("WGS", sequencing_type="SELECTIVE_WGS", amplification_route="SWGA"))
    assert new["schema_version"] == "pfvarx.config/4" and new["assay_design"] == "wgs"
    assert set(new["assay"]) == {"two_colour_chemistry", "platform"} and catalog is None
    assert defaults == {"technique": "SELECTIVE_WGS", "amplification_route": "SWGA"}
    new, catalog, defaults = plan_v3(schema3(
        "TARGETED_UMI", "BAM", amplification_route="AMPLICON", read_structure="8M+T 8M+T",
        targets_bed="inputs/targets/t.bed", primers_tsv="inputs/targets/p.tsv"))
    assert new["assay_design"] == "targeted_umi" and new["assay"]["library_catalog"]
    LibraryCatalog.model_validate(catalog)
    assert catalog["library_protocols"][0]["umi"]["read_structure"] == "8M+T 8M+T"
    assert defaults["processing_state"] == "consensus_primers_trimmed"  # schema 3 took imports as finished
    assert plan_v3(schema3("TARGETED_NON_UMI", targets_bed="t", primers_tsv="p"))[0]["assay_design"] == "targeted_no_umi"


def test_schema3_config_is_refused_with_the_migration_command() -> None:
    parsed = {"config": schema3("WGS"), "quality": {}, "compute": {}}
    with pytest.raises(ConfigError, match="config-migrate --apply"):
        validate_bundle(REPO, parsed)


def test_migration_fills_only_empty_cells_and_keeps_backups(tmp_path: Path) -> None:
    (tmp_path / "config").mkdir()
    cfg = schema3("TARGETED_NON_UMI", amplification_route="AMPLICON", targets_bed="inputs/targets/t.bed",
                  primers_tsv="inputs/targets/p.tsv")
    (tmp_path / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    sheet = tmp_path / "inputs/metadata/samples.tsv"
    sheet.parent.mkdir(parents=True)
    sheet.write_text("sample_id\tlane\tlibrary_id\ttechnique\nS1\tL1\tA\tTBS_NON_UMI\nS2\tL1\tB\t\n")
    assert migrate_workspace(tmp_path)["status"] == "PLANNED"
    assert yaml.safe_load((tmp_path / "config/config.yaml").read_text())["schema_version"] == "pfvarx.config/3"
    report = migrate_workspace(tmp_path, apply=True)
    assert report["status"] == "APPLIED" and (tmp_path / "config/config.yaml.before-config4").is_file()
    assert yaml.safe_load((tmp_path / "config/config.yaml").read_text())["assay_design"] == "targeted_no_umi"
    rows = [line.split("\t") for line in sheet.read_text().splitlines()]
    header = rows[0]
    assert rows[1][header.index("technique")] == "TBS_NON_UMI"  # existing value kept
    assert rows[1][header.index("panel_id")] == rows[2][header.index("panel_id")] == "LEGACY_PANEL"
    assert detect(tmp_path / "inputs/metadata/library_catalog.yaml") == "LIBRARY_CATALOG"
    assert migrate_workspace(tmp_path)["status"] == "NOT_NEEDED"


def test_umi_settings_fall_back_to_the_parameter_defaults() -> None:
    settings = umi_settings(None, PARAMS)
    assert settings["read_structure"]["value"] is None
    assert settings["group_strategy"] == {"value": PARAMS["layer1.fgbio_group_strategy"]["value"],
                                          "source": "PARAMETER_DEFAULT"}


def test_migration_without_a_manifest_reports_the_values_to_record(tmp_path: Path) -> None:
    (tmp_path / "config").mkdir()
    cfg = schema3("WGS", amplification_route="SWGA")
    (tmp_path / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    report = migrate_workspace(tmp_path, apply=True)
    assert report["status"] == "APPLIED" and report["manifest_cells_filled"] == 0
    assert "'amplification_route': 'SWGA'" in report["manifest_not_found"]


ADAPTER_CATALOG = LibraryCatalog.model_validate({
    "schema_version": CATALOG_SCHEMA,
    "library_protocols": [
        {"library_protocol_id": "NXT", "adapters": {"kit": "NEXTERA"}},
        {"library_protocol_id": "OWN", "adapters": {"kit": "CUSTOM", "read1": "ACGTACGTAC", "read2": "TGCATGCATG"}},
        {"library_protocol_id": "BARE"},
    ],
})


def resolve_adapters(rows: list[dict[str, str]]) -> tuple[list[str], dict[str, Any]]:
    rep, facts, _ = assay.resolve(Design.WGS, InputType.FASTQ, rows, ADAPTER_CATALOG, PARAMS)
    return [f.detail for f in rep.findings if f.level == "FAIL"], facts


def test_default_adapter_kit_is_explicit_not_detected() -> None:
    # fastp detection picks its poly-T FlowCell1 adapter on AT-rich reads; never the default.
    assert PARAMS["layer1.fastp_adapter_kit_default"]["value"] == "TRUSEQ"
    failures, facts = resolve_adapters([row()])
    assert failures == []
    assert facts["S1"].adapters == {"kit": "TRUSEQ", "read1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                                    "read2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "source": "PARAMETER_DEFAULT"}


def test_adapters_come_from_the_manifest_or_the_protocol() -> None:
    _, facts = resolve_adapters([row(adapter_kit="nextera")])
    assert facts["S1"].adapters["kit"] == "NEXTERA" and facts["S1"].adapters["source"] == "manifest"
    _, facts = resolve_adapters([row(library_protocol_id="OWN")])
    assert facts["S1"].adapters == {"kit": "CUSTOM", "read1": "ACGTACGTAC", "read2": "TGCATGCATG",
                                    "source": "library protocol"}
    _, facts = resolve_adapters([row(adapter_kit="AUTO")])
    assert facts["S1"].adapters == {"kit": "AUTO", "read1": None, "read2": None, "source": "manifest"}
    _, facts = resolve_adapters([row(library_protocol_id="BARE", adapter_kit=".")])
    assert facts["S1"].adapters["source"] == "PARAMETER_DEFAULT"


@pytest.mark.parametrize("cells, message", [
    ({"adapter_kit": "SMARTER"}, "unknown adapter_kit 'SMARTER'"),
    ({"adapter_kit": "TRUSEQ", "library_protocol_id": "NXT"}, "adapter kits disagree"),
    ({"adapter_kit": "CUSTOM"}, "adapter_kit CUSTOM needs a library_protocol_id"),
])
def test_adapter_problems_fail_the_sample(cells: dict[str, str], message: str) -> None:
    failures, facts = resolve_adapters([row(**cells)])
    assert any(message in f and f.startswith("sample S1 (library S1_lib)") for f in failures), failures
    assert facts["S1"].adapters is None


def test_libraries_of_one_sample_must_share_a_kit() -> None:
    failures, _ = resolve_adapters([row(adapter_kit="TRUSEQ"), {**row(adapter_kit="NEXTERA"), "lane": "L2"}])
    assert any("libraries disagree on adapter_kit" in f for f in failures)


@pytest.mark.parametrize("adapters, message", [
    ({"kit": "CUSTOM", "read1": "ACGTACGTAC"}, "needs read1 and read2"),
    ({"kit": "TRUSEQ", "read1": "ACGTACGTAC"}, "only for kit CUSTOM"),
    ({"kit": "CUSTOM", "read1": "ACGTNNNNAC", "read2": "ACGTACGTAC"}, "pattern"),
])
def test_catalog_adapter_definitions_are_checked(adapters: dict[str, str], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        LibraryCatalog.model_validate({"schema_version": CATALOG_SCHEMA,
                                       "library_protocols": [{"library_protocol_id": "X", "adapters": adapters}]})
