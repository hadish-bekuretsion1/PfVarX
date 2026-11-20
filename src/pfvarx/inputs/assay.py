"""Resolve what each library is, from the sample manifest and the library catalog.

``assay_design`` (the run setting) selects the processing branch. Each manifest row then
documents its library with optional columns:

``technique``            WGS | SELECTIVE_WGS | AMPLICON | TBS_UMI | TBS_NON_UMI (empty, null,
                         NA or . = unspecified)
``amplification_route``  DIRECT | SWGA | AMPLICON (provenance)
``panel_id``, ``panel_version``  the targeted panel (targeted designs; one panel per run)
``library_protocol_id``  a protocol in the library catalog (UMI extraction and consensus)
``processing_state``     imported BAM/CRAM only: what has already been done
``adapter_kit``          TRUSEQ | NEXTERA | AUTO | CUSTOM: the 3' adapters fastp trims (empty = the
                         library protocol's ``adapters``, else ``layer1.fastp_adapter_kit_default``)

Every conflict between the design, a technique label, the route, the protocol and the panel
is reported with the sample and library ID, as an input-validation failure for that sample.
The resolved facts go into the frozen plan, where the workflow reads them.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from pfvarx.configuration.library import LibraryCatalog, Panel, adapter_settings, umi_settings
from pfvarx.constants import (
    DESIGN_STATES,
    TECHNIQUE_DESIGNS,
    TECHNIQUE_ROUTE,
    AdapterKit,
    CountUnit,
    Design,
    InputType,
    ProcessingState,
    Reason,
    Route,
    Technique,
)
from pfvarx.inputs.validate import FileReport

NULL_VALUES = {"", "null", "none", "na", "n/a", "."}
DESIGN_ROUTES = {Design.WGS: (Route.DIRECT, Route.SWGA),
                 Design.TARGETED_UMI: (Route.AMPLICON,), Design.TARGETED_NO_UMI: (Route.AMPLICON,)}


def blank(value: str | None) -> bool:
    return value is None or value.strip().lower() in NULL_VALUES


@dataclass
class LibraryFacts:
    """Resolved assay facts for one sample (all of its libraries)."""

    sample_id: str
    library_ids: list[str]
    design: str
    technique: str | None = None
    technique_source: str | None = None
    amplification_route: str | None = None
    route_source: str | None = None
    panel_id: str | None = None
    panel_version: str | None = None
    library_protocol_id: str | None = None
    processing_state: str | None = None
    entry_steps: list[str] = field(default_factory=list)
    count_unit: str = CountUnit.READS.value
    umi: dict[str, Any] | None = None
    adapters: dict[str, Any] | None = None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def entry_steps(design: Design, state: ProcessingState | None) -> list[str]:
    """Steps still needed for an imported alignment in its documented state."""
    if design == Design.WGS:
        return ["mark_duplicates"] if state == ProcessingState.ALIGNED else \
            [] if state == ProcessingState.DUPLICATES_MARKED else ["mark_duplicates_if_unrecorded"]
    if design == Design.TARGETED_NO_UMI:
        return ["trim_primers"] if state == ProcessingState.ALIGNED else []
    return {ProcessingState.ALIGNED: ["umi_consensus", "trim_primers"],
            ProcessingState.CONSENSUS: ["trim_primers"]}.get(state, []) if state else []


def resolve(design: Design, input_type: InputType, rows: list[dict[str, str]],
            catalog: LibraryCatalog | None, params: dict[str, dict[str, Any]]
            ) -> tuple[FileReport, dict[str, LibraryFacts], Panel | None]:
    """Resolve every sample; findings name the sample and library and the conflicting values."""
    rep = FileReport(path="samples_tsv + assay.library_catalog", kind="LIBRARY_ASSAY")
    by_sample: dict[str, list[tuple[int, dict[str, str]]]] = {}
    for line, row in enumerate(rows, start=2):
        if row.get("sample_id"):
            by_sample.setdefault(row["sample_id"], []).append((line, row))
    facts: dict[str, LibraryFacts] = {}
    for sid, members in by_sample.items():
        libraries = sorted({r.get("library_id", "") for _, r in members if r.get("library_id")})
        fact = LibraryFacts(sample_id=sid, library_ids=libraries, design=design.value)

        def fail(code: Reason, text: str, line: int, library: str | None = None) -> None:
            where = f"sample {sid}" + (f" (library {library})" if library else "")
            rep.fail(code, f"{where}: {text}", line)

        for column in ("technique", "amplification_route", "panel_id", "panel_version",
                       "library_protocol_id", "processing_state", "adapter_kit"):
            values = {(r.get(column) or "").strip() for _, r in members if not blank(r.get(column))}
            if len(values) > 1:
                fail(Reason.TECHNIQUE_CONFLICT, f"libraries disagree on {column}: {sorted(values)}; "
                     "register libraries that differ as separate samples", members[0][0])
        line, row = members[0]
        library = row.get("library_id") or None
        # Library protocol
        protocol = None
        protocol_id = next((r["library_protocol_id"].strip() for _, r in members
                            if not blank(r.get("library_protocol_id"))), None)
        if protocol_id:
            protocol = catalog.protocol(protocol_id) if catalog else None
            if protocol is None:
                fail(Reason.PROTOCOL_UNRESOLVED, f"library_protocol_id {protocol_id!r} is not in the library "
                     "catalog", line, library)
            fact.library_protocol_id = protocol_id
        # Technique label: manifest, else protocol; both must suit the design and agree.
        labels: dict[str, str] = {}
        for _, r in members:
            value = r.get("technique")
            if not blank(value):
                labels["manifest"] = str(value).strip().upper()
        if protocol and protocol.technique:
            labels["library protocol"] = protocol.technique.value
        for source, label in labels.items():
            try:
                technique = Technique(label)
            except ValueError:
                fail(Reason.TECHNIQUE_CONFLICT, f"unknown technique {label!r} ({source}); use "
                     "WGS, SELECTIVE_WGS, AMPLICON, TBS_UMI, TBS_NON_UMI or leave it empty", line, library)
                continue
            if design not in TECHNIQUE_DESIGNS[technique]:
                allowed = " or ".join(d.value for d in TECHNIQUE_DESIGNS[technique])
                fail(Reason.TECHNIQUE_CONFLICT, f"technique {label} ({source}) requires assay_design {allowed}, "
                     f"but the run's assay_design is {design.value}", line, library)
        if len(set(labels.values())) > 1:
            fail(Reason.TECHNIQUE_CONFLICT, f"technique labels disagree: {labels}", line, library)
        if labels:
            source, label = next(iter(labels.items()))
            fact.technique, fact.technique_source = label, source
        # Panel (targeted designs; one panel for the whole run, checked below).
        panel_id = next((r["panel_id"].strip() for _, r in members if not blank(r.get("panel_id"))), None)
        version = next((r["panel_version"].strip() for _, r in members if not blank(r.get("panel_version"))), None)
        panel = None
        if design.targeted:
            if not panel_id or not version:
                fail(Reason.PANEL_UNRESOLVED, "targeted designs need panel_id and panel_version in the manifest",
                     line, library)
            else:
                panel = catalog.panel(panel_id, version) if catalog else None
                if panel is None:
                    fail(Reason.PANEL_UNRESOLVED, f"panel {panel_id} version {version} is not in the library "
                         "catalog", line, library)
            fact.panel_id, fact.panel_version = panel_id, version
        elif panel_id or version:
            fail(Reason.TECHNIQUE_CONFLICT, "panel_id/panel_version given but assay_design is wgs", line, library)
        # Amplification route: every documented source must agree and suit the design.
        routes: dict[str, str] = {}
        given = next((r["amplification_route"].strip().upper() for _, r in members
                      if not blank(r.get("amplification_route"))), None)
        if given:
            routes["manifest"] = given
        if protocol and protocol.amplification_route:
            routes["library protocol"] = protocol.amplification_route.value
        if fact.technique and Technique(fact.technique) in TECHNIQUE_ROUTE:
            routes[f"technique {fact.technique}"] = TECHNIQUE_ROUTE[Technique(fact.technique)].value
        if panel is not None and panel.primers_tsv is not None:
            routes["panel primers"] = Route.AMPLICON.value
        for source, value in routes.items():
            if value not in {r.value for r in Route}:
                fail(Reason.TECHNIQUE_CONFLICT, f"unknown amplification_route {value!r} ({source}); use DIRECT, "
                     "SWGA or AMPLICON", line, library)
            elif Route(value) not in DESIGN_ROUTES[design]:
                allowed = " or ".join(r.value for r in DESIGN_ROUTES[design])
                fail(Reason.TECHNIQUE_CONFLICT, f"amplification_route {value} ({source}) does not fit "
                     f"assay_design {design.value} (use {allowed})", line, library)
        if len(set(routes.values())) > 1:
            fail(Reason.TECHNIQUE_CONFLICT, f"amplification routes disagree: {routes}", line, library)
        if routes:
            source, value = next(iter(routes.items()))
            fact.amplification_route, fact.route_source = value, source
        if panel is not None and panel.primers_tsv is None and fact.amplification_route == Route.AMPLICON:
            fail(Reason.TECHNIQUE_CONFLICT, f"amplicon library but panel {panel_id} {version} defines no primers "
                 "(primers_tsv) to trim", line, library)
        # Adapters: manifest, else protocol, else the parameter default; sources must agree.
        kits: dict[str, str] = {}
        given_kit = next((r["adapter_kit"].strip().upper() for _, r in members
                          if not blank(r.get("adapter_kit"))), None)
        if given_kit:
            kits["manifest"] = given_kit
        if protocol and protocol.adapters:
            kits["library protocol"] = protocol.adapters.kit.value
        known = {k.value for k in AdapterKit}
        for source, value in kits.items():
            if value not in known:
                fail(Reason.ADAPTER_UNRESOLVED, f"unknown adapter_kit {value!r} ({source}); use TRUSEQ, NEXTERA, "
                     "AUTO, or CUSTOM with sequences in the library protocol", line, library)
        if len(set(kits.values())) > 1:
            fail(Reason.ADAPTER_UNRESOLVED, f"adapter kits disagree: {kits}", line, library)
        elif kits and next(iter(kits.values())) in known:
            source, value = next(iter(kits.items()))
            kit = AdapterKit(value)
            if kit == AdapterKit.CUSTOM and not (protocol and protocol.adapters):
                fail(Reason.ADAPTER_UNRESOLVED, "adapter_kit CUSTOM needs a library_protocol_id whose protocol "
                     "gives adapters.read1 and adapters.read2", line, library)
            else:
                fact.adapters = adapter_settings(kit, source, protocol)
        elif not kits:
            default = AdapterKit(str(params["layer1.fastp_adapter_kit_default"]["value"]).upper())
            fact.adapters = adapter_settings(default, "PARAMETER_DEFAULT")
        # UMI settings and the count unit.
        if design == Design.TARGETED_UMI:
            fact.count_unit = CountUnit.CONSENSUS_MOLECULES.value
            fact.umi = umi_settings(protocol, params)
        elif protocol and protocol.umi:
            fail(Reason.TECHNIQUE_CONFLICT, f"library protocol {protocol_id} defines UMIs but assay_design is "
                 f"{design.value}", line, library)
        # Imported alignments: the documented processing state.
        state_text = next((r["processing_state"].strip().lower() for _, r in members
                           if not blank(r.get("processing_state"))), None)
        state = None
        if input_type in (InputType.BAM, InputType.CRAM):
            if state_text:
                try:
                    state = ProcessingState(state_text)
                except ValueError:
                    state = None
                if state is None or state not in DESIGN_STATES[design]:
                    allowed = ", ".join(s.value for s in DESIGN_STATES[design])
                    fail(Reason.ALIGNMENT_PROCESSING_STATE, f"processing_state {state_text!r} is not valid for "
                         f"assay_design {design.value} (use {allowed})", line, library)
                    state = None
            elif design.targeted:
                allowed = ", ".join(s.value for s in DESIGN_STATES[design])
                fail(Reason.ALIGNMENT_PROCESSING_STATE, "imported targeted alignments need processing_state "
                     f"({allowed})", line, library)
            fact.processing_state = state.value if state else None
            fact.entry_steps = entry_steps(design, state)
        elif state_text:
            rep.warn(Reason.ALIGNMENT_PROCESSING_STATE,
                     f"sample {sid}: processing_state applies only to imported BAM/CRAM; ignored for FASTQ")
        needs_structure = design == Design.TARGETED_UMI and input_type == InputType.FASTQ
        if needs_structure and not (protocol and protocol.umi):
            fail(Reason.PROTOCOL_UNRESOLVED, "UMI FASTQ input needs a library_protocol_id whose protocol has a "
                 "umi.read_structure (how to extract the UMIs)", line, library)
        facts[sid] = fact
    run_panel = None
    if design.targeted:
        panels = {(f.panel_id, f.panel_version) for f in facts.values() if f.panel_id and f.panel_version}
        if len(panels) > 1:
            rep.fail(Reason.PANEL_UNRESOLVED, "one panel per run: samples use "
                     + "; ".join(f"{p[0]} {p[1]} ({', '.join(s for s, f in sorted(facts.items()) if (f.panel_id, f.panel_version) == p)})"
                                 for p in sorted(panels)))
        elif panels and catalog:
            run_panel = catalog.panel(*next(iter(panels)))
    rep.facts = {"design": design.value, "samples": len(facts),
                 "panel": None if run_panel is None else {"panel_id": run_panel.panel_id,
                                                          "panel_version": run_panel.panel_version}}
    return rep, facts, run_panel


def sample_failures(rep: FileReport, rows: list[dict[str, str]]) -> set[str]:
    """Samples with a FAIL finding (findings carry the manifest line of the sample's first row)."""
    out: set[str] = set()
    for finding in rep.findings:
        if finding.level == "FAIL" and finding.record is not None and 2 <= finding.record <= len(rows) + 1:
            sid = rows[finding.record - 2].get("sample_id")
            if sid:
                out.add(sid)
    return out


def check_alignment_state(fact: LibraryFacts, evidence: dict[str, Any], rep: FileReport) -> None:
    """Compare the documented processing state with what the alignment itself shows.

    Content that contradicts the documented state is a failure (consensus reads documented as
    raw reads would be collapsed again; raw reads documented as consensus would skip
    consensus). A trimming or duplicate-marking claim that the header does not record is a
    warning: it cannot be checked from the reads, so the documented state is used."""
    sid = fact.sample_id
    state = fact.processing_state
    consensus = evidence["consensus_tagged"] > 0
    rep.facts["processing_evidence"] = evidence
    if fact.design == Design.TARGETED_UMI.value:
        if state == ProcessingState.ALIGNED and consensus:
            rep.fail(Reason.ALIGNMENT_PROCESSING_STATE, f"sample {sid}: processing_state is aligned (raw UMI reads) "
                     "but the reads carry fgbio consensus tags (cD); use consensus or consensus_primers_trimmed")
        if state == ProcessingState.ALIGNED and evidence["umi_tagged"] == 0:
            rep.fail(Reason.ALIGNMENT_PROCESSING_STATE, f"sample {sid}: raw UMI reads need the UMI in the RX tag; "
                     "none of the inspected reads has one")
        if state in (ProcessingState.CONSENSUS, ProcessingState.CONSENSUS_PRIMERS_TRIMMED) and not consensus:
            rep.fail(Reason.ALIGNMENT_PROCESSING_STATE, f"sample {sid}: processing_state is {state} but no "
                     "inspected read carries an fgbio consensus tag (cD)")
    elif consensus:
        rep.fail(Reason.ALIGNMENT_PROCESSING_STATE, f"sample {sid}: reads carry UMI consensus tags (cD) but "
                 f"assay_design is {fact.design}; use targeted_umi")
    trimmed_claim = state in (ProcessingState.PRIMERS_TRIMMED, ProcessingState.CONSENSUS_PRIMERS_TRIMMED)
    if "trim_primers" in fact.entry_steps and evidence["trim_primers_recorded"]:
        rep.fail(Reason.ALIGNMENT_PROCESSING_STATE, f"sample {sid}: the header records fgbio TrimPrimers but "
                 f"processing_state {state} says primers are not trimmed; trimming twice would remove template "
                 "bases. Correct processing_state")
    if trimmed_claim and not evidence["trim_primers_recorded"]:
        rep.warn(Reason.ALIGNMENT_PROCESSING_STATE, f"sample {sid}: processing_state {state} documents primer "
                 "trimming that the header does not record (no TrimPrimers @PG); the documented state is used")
    if state == ProcessingState.DUPLICATES_MARKED and not evidence["duplicates_marked_recorded"]:
        rep.warn(Reason.ALIGNMENT_PROCESSING_STATE, f"sample {sid}: duplicates_marked is documented but no "
                 "duplicate-marking program is recorded in @PG; the documented state is used")
