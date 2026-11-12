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
    threads: 1
    resources:
        mem_mb=mem(2000),
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
        "--out-focal {output.focal:q} --out-qc {output.qc:q} > {log:q} 2>&1"


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
    log:
        str(LOGS / "norm" / "{sample}.log"),
    shell:
        "(bcftools norm -f {params.ref:q} -m {params.m} --check-ref {params.c} --multi-overlaps {params.o} "
        "-Oz -o {output.vcf:q} {input:q} && bcftools index -t {output.vcf:q}) > {log:q} 2>&1"
