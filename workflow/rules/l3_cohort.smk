# Layer 3: variant-only merge (NOT joint genotyping) and the cohort evidence matrix.


rule cohort_merge:
    input:
        vcfs=lambda wc: [str(R / "variants" / "samples" / f"{s}.vcf.gz") for s in passing("L3")],
    output:
        vcf=str(R / "cohort" / "cohort.vcf.gz"),
        tbi=str(R / "cohort" / "cohort.vcf.gz.tbi"),
        samples=str(R / "cohort" / "cohort_samples.txt"),
    params:
        # Sample IDs in merge order (IDs never contain dots, so the file stem is the ID).
        ids=lambda wc, input: " ".join(Path(v).name.split(".")[0] for v in input.vcfs),
        n=lambda wc, input: len(input.vcfs),
    log:
        str(LOGS / "cohort" / "merge.log"),
    # final_targets requests the cohort only when a sample is permitted into Layer 3. One
    # sample is copied (bcftools merge needs two files); several are merged variant-only, with
    # missing genotypes left missing (no -0: absence of a call is not a reference call).
    shell:
        "(printf '%s\\n' {params.ids} > {output.samples:q} "
        "&& if [ {params.n} -eq 1 ]; "
        "then cp {input.vcfs:q} {output.vcf:q}; "
        "else bcftools merge -m none -i - -Oz -o {output.vcf:q} {input.vcfs:q}; fi "
        "&& bcftools index -t -f {output.vcf:q}) > {log:q} 2>&1"


rule cohort_status:
    """Record whether a cohort was formed; with none, mark requested analyses as skipped."""
    input:
        unpack(cohort_status_inputs),
        gates=lambda wc: [str(R / "qc" / "samples" / s / "qc_decision.json") for s in sorted(SAMPLES) if s not in EXCLUDED],
    output:
        str(R / "cohort" / "cohort_status.json"),
    params:
        empty=lambda wc: "--empty" if not passing("L3") else "",
        run=str(RUN),
        root=str(ROOT),
    log:
        str(LOGS / "cohort" / "status.log"),
    shell:
        "python -m pfvarx.downstream.run --cohort-status {params.empty} --run-dir {params.run:q} "
        "--root {params.root:q} > {log:q} 2>&1"


rule cohort_evidence:
    input:
        vcf=str(R / "cohort" / "cohort.vcf.gz"),
        samples=str(R / "cohort" / "cohort_samples.txt"),
        bams=lambda wc: [str(R / "alignments" / f"{s}.bam") for s in passing("L3")],
        callable=lambda wc: [str(R / "variants" / "samples" / f"{s}.callable.bed") for s in passing("L3")],
    output:
        matrix=str(R / "cohort" / "evidence_matrix.tsv.gz"),
        missing=str(R / "cohort" / "missingness.tsv"),
        sites=str(R / "cohort" / "sites.tsv"),
    params:
        frozen=FROZEN,
        ref=PREF,
    resources:
        mem_mb=mem(2000),
        runtime=RUNTIME,
    log:
        str(LOGS / "cohort" / "evidence.log"),
    shell:
        "python -m pfvarx.cohort.evidence --vcf {input.vcf:q} --samples {input.samples:q} --bams {input.bams:q} "
        "--callable {input.callable:q} --reference {params.ref:q} --frozen-config {params.frozen:q} "
        "--out-matrix {output.matrix:q} --out-missing {output.missing:q} --out-sites {output.sites:q} > {log:q} 2>&1"
