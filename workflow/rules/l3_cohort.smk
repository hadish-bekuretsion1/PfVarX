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
        Path(output.samples).write_text("\n".join(Path(v).name.split(".")[0] for v in input.vcfs) + "\n")
        if len(input.vcfs) == 1:
            shell("(cp {input.vcfs[0]:q} {output.vcf:q} && bcftools index -t -f {output.vcf:q}) > {log:q} 2>&1")
        else:
            shell("(bcftools merge -m none -i - -Oz -o {output.vcf:q} {input.vcfs:q} "
                  "&& bcftools index -t {output.vcf:q}) > {log:q} 2>&1")


rule cohort_evidence:
    input:
        vcf=str(R / "cohort" / "cohort.vcf.gz"),
        samples=str(R / "cohort" / "cohort_samples.txt"),
    output:
        matrix=str(R / "cohort" / "evidence_matrix.tsv.gz"),
        missing=str(R / "cohort" / "missingness.tsv"),
        sites=str(R / "cohort" / "sites.tsv"),
    params:
        bams=lambda wc: [str(R / "alignments" / f"{s}.bam") for s in passing("L3")],
        callable=lambda wc: [str(R / "variants" / "samples" / f"{s}.callable.bed") for s in passing("L3")],
        frozen=FROZEN,
        ref=PREF,
    resources:
        mem_mb=mem(2000),
        runtime=RUNTIME,
    log:
        str(LOGS / "cohort" / "evidence.log"),
    shell:
        "python -m pfvarx.cohort.evidence --vcf {input.vcf:q} --samples {input.samples:q} --bams {params.bams:q} "
        "--callable {params.callable:q} --reference {params.ref:q} --frozen-config {params.frozen:q} "
        "--out-matrix {output.matrix:q} --out-missing {output.missing:q} --out-sites {output.sites:q} > {log:q} 2>&1"
