# Layer 3: variant-only merge (NOT joint genotyping) and the cohort evidence matrix.


rule cohort_merge:
    input:
        vcfs=lambda wc: [str(R / "variants" / "samples" / f"{s}.vcf.gz") for s in passing("L3")],
    output:
        vcf=str(R / "cohort" / "cohort.vcf.gz"),
        tbi=str(R / "cohort" / "cohort.vcf.gz.tbi"),
        samples=str(R / "cohort" / "cohort_samples.txt"),
    log:
        str(LOGS / "cohort" / "merge.log"),
    run:
        if not input.vcfs:
            # final_targets never requests the cohort without an admitted sample; fail clearly
            # rather than with a bcftools usage message if that contract is ever broken.
            raise ValueError("cohort_merge needs at least one sample permitted into Layer 3")
        Path(output.samples).write_text("\n".join(Path(v).name.split(".")[0] for v in input.vcfs) + "\n")
        if len(input.vcfs) == 1:
            shell("(cp {input.vcfs[0]:q} {output.vcf:q} && bcftools index -t -f {output.vcf:q}) > {log:q} 2>&1")
        else:
            shell("(bcftools merge -m none -i - -Oz -o {output.vcf:q} {input.vcfs:q} "
                  "&& bcftools index -t {output.vcf:q}) > {log:q} 2>&1")


def cohort_status_inputs(wc):
    if passing("L3"):
        return {"samples": str(R / "cohort" / "cohort_samples.txt")}
    return {}


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
