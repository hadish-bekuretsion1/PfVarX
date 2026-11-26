"""Centralised, versioned identifiers, status vocabularies and reason codes.

Every status that crosses a layer boundary is defined here so that Layer 0,
the Snakemake rules, the downstream analyses and the dashboard share one
vocabulary. Changing a value is a schema change and requires a version bump.
"""

from __future__ import annotations

from enum import StrEnum

STATUS_SCHEMA_VERSION = "pfvarx.status/1"

# Default random seed of a project (config random_seed) and of the benchmark partition.
DEFAULT_RANDOM_SEED = 20261106


class Purpose(StrEnum):
    ANALYSIS = "analysis"
    BENCHMARKING = "benchmarking"


class FinalOutput(StrEnum):
    FO001 = "FO001"  # caller-ready alignment, index and QC
    FO002 = "FO002"  # filtered per-sample variant files, evidence and QC
    FO003 = "FO003"  # compatible variant-only cohort VCF/BCF and missingness
    FO004 = "FO004"  # user-selected eligible DA001-DA008 outputs


class InputType(StrEnum):
    FASTQ = "FASTQ"
    BAM = "BAM"
    CRAM = "CRAM"
    VCF = "VCF"


class Design(StrEnum):
    """Assay design (config ``assay_design``): the setting that selects the processing branch."""

    WGS = "wgs"
    TARGETED_UMI = "targeted_umi"
    TARGETED_NO_UMI = "targeted_no_umi"

    @classmethod
    def parse(cls, value: object) -> Design:
        """Accept the current values and the spellings of config schema 3 and earlier runs."""
        text = str(value.value if isinstance(value, StrEnum) else value).strip().lower()
        return cls(LEGACY_DESIGNS.get(text, text))

    @property
    def targeted(self) -> bool:
        return self is not Design.WGS


LEGACY_DESIGNS = {"targeted_non_umi": "targeted_no_umi"}


class Route(StrEnum):
    """Amplification route, recorded per library as provenance (not a run setting)."""

    DIRECT = "DIRECT"  # confirmed WGS without selective amplification
    SWGA = "SWGA"  # documented selective whole-genome amplification
    AMPLICON = "AMPLICON"  # targeted library prepared by amplicon PCR


class Technique(StrEnum):
    """Optional sequencing-technique label in the sample manifest or library protocol."""

    WGS = "WGS"
    SELECTIVE_WGS = "SELECTIVE_WGS"
    AMPLICON = "AMPLICON"
    TBS_UMI = "TBS_UMI"
    TBS_NON_UMI = "TBS_NON_UMI"


class AdapterKit(StrEnum):
    """Adapter family of a library, given to fastp explicitly (per library, not per run).

    fastp's own adapter detection is not used by default: on AT-rich *P. falciparum* reads it
    can pick its built-in poly-T ``FlowCell1`` sequence and trim genuine T runs.
    """

    TRUSEQ = "TRUSEQ"  # Illumina TruSeq / most ligation kits (AGATCGGAAGAGC...)
    NEXTERA = "NEXTERA"  # Illumina Nextera / tagmentation kits (CTGTCTCTTATACACATCT)
    CUSTOM = "CUSTOM"  # sequences given in the library protocol (adapters.read1 / read2)
    AUTO = "AUTO"  # fastp --detect_adapter_for_pe (only when the kit is unknown)


# Read-1 and read-2 3' adapter sequences fastp trims for each kit.
ADAPTER_SEQUENCES: dict[AdapterKit, tuple[str, str]] = {
    AdapterKit.TRUSEQ: ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
    AdapterKit.NEXTERA: ("CTGTCTCTTATACACATCT", "CTGTCTCTTATACACATCT"),
}


# Designs each technique label may be used with, and the route a label documents.
TECHNIQUE_DESIGNS = {
    Technique.WGS: (Design.WGS,),
    Technique.SELECTIVE_WGS: (Design.WGS,),
    Technique.AMPLICON: (Design.TARGETED_UMI, Design.TARGETED_NO_UMI),
    Technique.TBS_UMI: (Design.TARGETED_UMI,),
    Technique.TBS_NON_UMI: (Design.TARGETED_NO_UMI,),
}
TECHNIQUE_ROUTE = {Technique.SELECTIVE_WGS: Route.SWGA, Technique.AMPLICON: Route.AMPLICON}


class ProcessingState(StrEnum):
    """Documented processing state of an imported BAM/CRAM (sample manifest ``processing_state``)."""

    ALIGNED = "aligned"  # aligned reads only (UMI designs: raw reads with UMIs in the RX tag)
    DUPLICATES_MARKED = "duplicates_marked"  # wgs
    PRIMERS_TRIMMED = "primers_trimmed"  # targeted_no_umi
    CONSENSUS = "consensus"  # targeted_umi: consensus reads, primers not yet trimmed
    CONSENSUS_PRIMERS_TRIMMED = "consensus_primers_trimmed"  # targeted_umi: consensus and trimmed


DESIGN_STATES = {
    Design.WGS: (ProcessingState.ALIGNED, ProcessingState.DUPLICATES_MARKED),
    Design.TARGETED_NO_UMI: (ProcessingState.ALIGNED, ProcessingState.PRIMERS_TRIMMED),
    Design.TARGETED_UMI: (ProcessingState.ALIGNED, ProcessingState.CONSENSUS,
                          ProcessingState.CONSENSUS_PRIMERS_TRIMMED),
}


class CountUnit(StrEnum):
    """What allele counts and depths count: sequenced reads, or UMI consensus molecules."""

    READS = "reads"
    CONSENSUS_MOLECULES = "consensus_molecules"


class Scale(StrEnum):
    SAMPLE = "SAMPLE"
    PAIRWISE = "PAIRWISE"
    COHORT = "COHORT"
    POPULATION = "POPULATION"
    BETWEEN_POPULATIONS = "BETWEEN_POPULATIONS"


class QC(StrEnum):
    PASS = "PASS"
    WARNING = "WARNING"
    FAIL = "FAIL"
    NOT_EVALUATED = "NOT_EVALUATED"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class Eligibility(StrEnum):
    ELIGIBLE = "ELIGIBLE"
    CONDITIONALLY_ELIGIBLE = "CONDITIONALLY_ELIGIBLE"
    INELIGIBLE = "INELIGIBLE"
    INDETERMINATE = "INDETERMINATE"
    UNSUPPORTED = "UNSUPPORTED"


class Terminal(StrEnum):
    COMPLETED = "COMPLETED"
    COMPLETED_WITH_EXCLUSIONS = "COMPLETED_WITH_EXCLUSIONS"
    SKIPPED_INELIGIBLE = "SKIPPED_INELIGIBLE"
    BLOCKED_MISSING_PREREQUISITE = "BLOCKED_MISSING_PREREQUISITE"
    FAILED_EXECUTION = "FAILED_EXECUTION"


class Clonality(StrEnum):
    MONOCLONAL = "MONOCLONAL"
    POLYCLONAL = "POLYCLONAL"
    INDETERMINATE = "INDETERMINATE"


class CallStatus(StrEnum):
    """Per-site output statuses of the PfVarX-Call decision core (Table 3.7)."""

    PASS = "PASS"
    FILTERED = "FILTERED"
    LOW_CONFIDENCE = "LOW_CONFIDENCE"
    INSUFFICIENT_DEPTH = "INSUFFICIENT_DEPTH"
    FOCAL_CALLABILITY_LOSS = "FOCAL_CALLABILITY_LOSS"
    CONTEXT_UNRELIABLE = "CONTEXT_UNRELIABLE"


class MarkerState(StrEnum):
    MUTANT = "MUTANT"
    MIXED = "MIXED"
    WILD_TYPE = "WILD_TYPE"
    BELOW_REPORTING_FLOOR = "BELOW_REPORTING_FLOOR"
    UNKNOWN = "UNKNOWN"


class ParameterSource(StrEnum):
    """Recorded source class for every numerical setting (construction prompt 1.3.7)."""

    SOFTWARE_DEFAULT = "SOFTWARE_DEFAULT"
    LITERATURE_SUPPORTED = "LITERATURE_SUPPORTED"
    ASSAY_VALIDATED = "ASSAY_VALIDATED"
    EMPIRICALLY_CALIBRATED = "EMPIRICALLY_CALIBRATED"
    UNSET_BLOCKING = "UNSET_BLOCKING"
    # Resolution of discrepancy D-03 (see docs/status.md): a value that
    # permits technical execution before calibration. Every output that depends on it
    # is labelled UNCALIBRATED; it is never reported as validated.
    UNCALIBRATED_PROVISIONAL = "UNCALIBRATED_PROVISIONAL"
    DECLARED_SCOPE = "DECLARED_SCOPE"
    FRAMEWORK_POLICY = "FRAMEWORK_POLICY"


class FailureClass(StrEnum):
    SCIENTIFIC_EXCLUSION = "SCIENTIFIC_EXCLUSION"
    QC_FAILURE = "QC_FAILURE"
    RESOURCE_FAILURE = "RESOURCE_FAILURE"
    SOFTWARE_FAILURE = "SOFTWARE_FAILURE"
    INPUT_FAILURE = "INPUT_FAILURE"


class Severity(StrEnum):
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class Reason(StrEnum):
    """Stable machine-readable reason codes."""

    # inputs
    FASTQ_TRUNCATED = "FASTQ_TRUNCATED"
    FASTQ_MATE_COUNT_MISMATCH = "FASTQ_MATE_COUNT_MISMATCH"
    FASTQ_MATE_NAME_MISMATCH = "FASTQ_MATE_NAME_MISMATCH"
    FASTQ_QUALITY_ENCODING = "FASTQ_QUALITY_ENCODING"
    FASTQ_MALFORMED_RECORD = "FASTQ_MALFORMED_RECORD"
    FILE_MISSING = "FILE_MISSING"
    REFERENCE_MISMATCH = "REFERENCE_MISMATCH"
    REFERENCE_MALFORMED = "REFERENCE_MALFORMED"
    BED_MALFORMED = "BED_MALFORMED"
    CORE_MASK_SUSPECT = "CORE_MASK_SUSPECT"
    CONTENT_MISMATCH = "CONTENT_MISMATCH"
    GFF_MALFORMED = "GFF_MALFORMED"
    CATALOG_MALFORMED = "CATALOG_MALFORMED"
    PRIMER_TABLE_MALFORMED = "PRIMER_TABLE_MALFORMED"
    METADATA_MISSING_FIELD = "METADATA_MISSING_FIELD"
    METADATA_DUPLICATE_SAMPLE = "METADATA_DUPLICATE_SAMPLE"
    ALIGNMENT_UNSORTED = "ALIGNMENT_UNSORTED"
    ALIGNMENT_NO_INDEX = "ALIGNMENT_NO_INDEX"
    ALIGNMENT_READ_GROUP = "ALIGNMENT_READ_GROUP"
    ALIGNMENT_PROCESSING_STATE = "ALIGNMENT_PROCESSING_STATE"
    UMI_READ_STRUCTURE = "UMI_READ_STRUCTURE"
    TECHNIQUE_CONFLICT = "TECHNIQUE_CONFLICT"
    NAMING_CONVENTION = "NAMING_CONVENTION"
    LIBRARY_CATALOG_INVALID = "LIBRARY_CATALOG_INVALID"
    PANEL_UNRESOLVED = "PANEL_UNRESOLVED"
    PROTOCOL_UNRESOLVED = "PROTOCOL_UNRESOLVED"
    ADAPTER_UNRESOLVED = "ADAPTER_UNRESOLVED"
    VCF_INCOMPATIBLE = "VCF_INCOMPATIBLE"
    # QC
    QC_GATE_FAIL = "QC_GATE_FAIL"
    QC_GATE_WARNING = "QC_GATE_WARNING"
    QC_GATE_UNCALIBRATED = "QC_GATE_UNCALIBRATED"
    QC_UPSTREAM_FAILED = "QC_UPSTREAM_FAILED"
    COHORT_QC_FAIL = "COHORT_QC_FAIL"
    COHORT_QC_WARNING = "COHORT_QC_WARNING"
    # eligibility
    DESIGN_NOT_SUPPORTED = "DESIGN_NOT_SUPPORTED"
    SCALE_NOT_SUPPORTED = "SCALE_NOT_SUPPORTED"
    CLONALITY_POLYCLONAL = "CLONALITY_POLYCLONAL"
    CLONALITY_INDETERMINATE = "CLONALITY_INDETERMINATE"
    CLONALITY_EVIDENCE_MISSING = "CLONALITY_EVIDENCE_MISSING"
    MONOCLONAL_COUNT_BELOW_MINIMUM = "MONOCLONAL_COUNT_BELOW_MINIMUM"
    INSUFFICIENT_INFORMATIVE_SITES = "INSUFFICIENT_INFORMATIVE_SITES"
    MISSINGNESS_ABOVE_MAXIMUM = "MISSINGNESS_ABOVE_MAXIMUM"
    GROUP_ASSIGNMENT_MISSING = "GROUP_ASSIGNMENT_MISSING"
    GROUP_TOO_SMALL = "GROUP_TOO_SMALL"
    METADATA_GEOGRAPHY_INVALID = "METADATA_GEOGRAPHY_INVALID"
    METADATA_DATE_INVALID = "METADATA_DATE_INVALID"
    SOURCE_OUTCOME_INCOMPLETE = "SOURCE_OUTCOME_INCOMPLETE"
    CATALOG_NOT_VERIFIED = "CATALOG_NOT_VERIFIED"
    METHOD_NOT_INSTALLED = "METHOD_NOT_INSTALLED"
    METHOD_NOT_VALIDATED = "METHOD_NOT_VALIDATED"
    NO_ELIGIBLE_PAIR = "NO_ELIGIBLE_PAIR"
    NO_ELIGIBLE_SAMPLES = "NO_ELIGIBLE_SAMPLES"
    NO_COHORT = "NO_COHORT"
    PANEL_LACKS_MARKERS = "PANEL_LACKS_MARKERS"
    NO_CALLABLE_MARKERS = "NO_CALLABLE_MARKERS"
    TOO_FEW_INDEPENDENT_LOCI = "TOO_FEW_INDEPENDENT_LOCI"
    PARAMETER_UNSET = "PARAMETER_UNSET"
    # resources
    RESOURCE_MEMORY_EXCEEDED = "RESOURCE_MEMORY_EXCEEDED"
    RESOURCE_TIME_LIMIT_EXCEEDED = "RESOURCE_TIME_LIMIT_EXCEEDED"
    RESOURCE_DISK_EXHAUSTED = "RESOURCE_DISK_EXHAUSTED"
    RESOURCE_CPU_LIMIT = "RESOURCE_CPU_LIMIT"
    RESOURCE_SCHEDULER_REJECTED = "RESOURCE_SCHEDULER_REJECTED"
    RESOURCE_ESTIMATE_UNAVAILABLE = "RESOURCE_ESTIMATE_UNAVAILABLE"
    # software
    TOOL_NOT_FOUND = "TOOL_NOT_FOUND"
    TOOL_NONZERO_EXIT = "TOOL_NONZERO_EXIT"
    RUN_REUSE_INCOMPATIBLE = "RUN_REUSE_INCOMPATIBLE"
    STRICT_RUN_STOP = "STRICT_RUN_STOP"


DOWNSTREAM_IDS = ("DA002", "DA003", "DA001", "DA004", "DA005", "DA006", "DA007", "DA008")
