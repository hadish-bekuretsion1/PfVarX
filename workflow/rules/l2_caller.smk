# Layer 2: PfVarX-Call, then reference-consistent normalisation (left-align, decompose).


rule pfvarx_call:
    input:
        bam=str(R / "alignments" / "{sample}.bam"),
    output:
        vcf=temp(str(R / "work" / "calls" / "{sample}.raw.vcf")),
        evidence=str(R / "variants" / "samples" / "{sample}.evidence.tsv"),
        callable=str(R / "variants" / "samples" / "{sample}.callable.bed"),
        focal=str(R / "variants" / "samples" / "{sample}.focal_loss.bed"),
        qc=str(R / "qc" / "samples" / "{sample}" / "{sample}.caller_qc.json"),
    params:
        ref=PREF,
        core=CORE,
        frozen=FROZEN,
        unreliable=str(ROOT / "config" / "registry" / "context_unreliable.bed"),
        # What AD/DP count: sequenced reads, or UMI consensus molecules (targeted_umi).
        unit=lambda wc: count_unit(wc.sample),
    # Chromosomes, indel rescoring and bootstrap replicates run in parallel; results do not
    # depend on the thread count. Measured on a 29x WGS sample (20.8 Mb core): 579 s and 3.6 GB
    # with 8 workers, 474 s and 3.8 GB with 24. Memory is set by the genome, not the workers,
    # and beyond 8 workers the gain is small, so 8 lets several samples run at once.
    threads: threads(8)
    resources:
        mem_mb=mem(5000),
        runtime=RUNTIME,
    retries: COMPUTE["max_retries"]
    log:
        str(LOGS / "pfvarx_call" / "{sample}.log"),
    benchmark:
        str(BENCH / "pfvarx_call" / "{sample}.tsv")
    shell:
        "python -m pfvarx.caller_adapter --bam {input.bam:q} --reference {params.ref:q} --core-mask {params.core:q} "
        "--frozen-config {params.frozen:q} --sample {wildcards.sample} --context-unreliable {params.unreliable:q} "
        "--out-vcf {output.vcf:q} --out-evidence {output.evidence:q} --out-callable {output.callable:q} "
        "--out-focal {output.focal:q} --out-qc {output.qc:q} --threads {threads} "
        "--count-unit {params.unit} > {log:q} 2>&1"


rule normalise_sample:
    input:
        str(R / "work" / "calls" / "{sample}.raw.vcf"),
    output:
        vcf=str(R / "variants" / "samples" / "{sample}.vcf.gz"),
        tbi=str(R / "variants" / "samples" / "{sample}.vcf.gz.tbi"),
    params:
        ref=PREF,
        m=P["layer3.norm_multiallelic"],
        c=P["layer3.norm_check_ref"],
        o=P["layer3.norm_multi_overlaps"],
        unit=lambda wc: count_unit(wc.sample),
    benchmark:
        str(BENCH / "normalise_sample" / "{sample}.tsv")
    log:
        str(LOGS / "norm" / "{sample}.log"),
    shell:
        "(bcftools norm -f {params.ref:q} -m {params.m} --check-ref {params.c} --multi-overlaps {params.o} -Ou {input:q} "
        "| bcftools annotate --header-line '##PfVarX_count_unit={params.unit}' -Oz -o {output.vcf:q} "
        "&& bcftools index -t {output.vcf:q}) > {log:q} 2>&1"
