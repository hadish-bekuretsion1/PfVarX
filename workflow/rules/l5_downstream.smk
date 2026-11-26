# Layer 5: downstream analyses with two-phase eligibility (checkpoint), then dynamic expansion.

# REQUESTED, INTERNAL, DOWN and ELIG and the input functions are in common.smk.


rule da001_clonality:
    """Clonal-status evidence (DA001), run when requested or as an internal prerequisite."""
    input:
        unpack(common_inputs),
        gates=lambda wc: [str(R / "qc" / "samples" / s / "qc_decision.json") for s in sorted(SAMPLES) if s not in EXCLUDED],
    output:
        manifest=str(DOWN / "DA001" / "status_manifest.json"),
        clonality=str(DOWN / "DA001" / "clonality.tsv"),
    params:
        run=str(RUN),
        root=str(ROOT),
    resources:
        mem_mb=mem("da001_clonality"),
        runtime=RUNTIME,
    log:
        str(LOGS / "downstream" / "DA001.log"),
    benchmark:
        str(BENCH / "downstream" / "DA001.tsv")
    shell:
        "python -m pfvarx.downstream.run --da DA001 --run-dir {params.run:q} --root {params.root:q} > {log:q} 2>&1"


checkpoint da_eligibility:
    input:
        unpack(eligibility_inputs),
    output:
        decisions=str(ELIG / "cohort_decisions.json"),
        matrix=str(ELIG / "sample_analysis_eligibility.tsv"),
    params:
        run=str(RUN),
        root=str(ROOT),
    benchmark:
        str(BENCH / "da_eligibility" / "da_eligibility.tsv")
    log:
        str(LOGS / "downstream" / "eligibility.log"),
    shell:
        "python -m pfvarx.downstream.run --eligibility --run-dir {params.run:q} --root {params.root:q} > {log:q} 2>&1"


rule da_execute:
    input:
        unpack(common_inputs),
        decisions=str(ELIG / "cohort_decisions.json"),
        source=source_manifest,
    output:
        manifest=str(DOWN / "{da}" / "status_manifest.json"),
    params:
        run=str(RUN),
        root=str(ROOT),
    wildcard_constraints:
        da="DA00[2345678]",
    threads: threads("da_execute")
    resources:
        mem_mb=mem("da_execute"),
        runtime=RUNTIME,
    log:
        str(LOGS / "downstream" / "{da}.log"),
    benchmark:
        str(BENCH / "downstream" / "{da}.tsv")
    shell:
        "python -m pfvarx.downstream.run --da {wildcards.da} --run-dir {params.run:q} --root {params.root:q} "
        "--threads {threads} > {log:q} 2>&1"


