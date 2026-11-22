# L1 QC metrics and the per-sample Layer 1 gate (checkpoint).

# TARGETS (the run panel's intervals) is defined in the Snakefile.


rule flagstat:
    input:
        str(R / "alignments" / "{sample}.bam"),
    output:
        str(R / "qc" / "samples" / "{sample}" / "{sample}.flagstat.json"),
    benchmark:
        str(BENCH / "flagstat" / "{sample}.tsv")
    log:
        str(LOGS / "flagstat" / "{sample}.log"),
    shell:
        "samtools flagstat -O json {input:q} > {output:q} 2> {log:q}"


rule mosdepth_core:
    input:
        str(R / "alignments" / "{sample}.bam"),
    output:
        dist=str(R / "qc" / "samples" / "{sample}" / "core.mosdepth.region.dist.txt"),
        summary=str(R / "qc" / "samples" / "{sample}" / "core.mosdepth.summary.txt"),
    params:
        prefix=lambda wc: str(R / "qc" / "samples" / wc.sample / "core"),
        bed=CORE,
        flag=P["layer1.mosdepth_exclude_flag"],
        mapq=P["layer2.min_mapping_quality"],
    threads: threads(2)
    resources:
        mem_mb=mem(1000),
        runtime=RUNTIME,
    benchmark:
        str(BENCH / "mosdepth_core" / "{sample}.tsv")
    log:
        str(LOGS / "mosdepth" / "{sample}.core.log"),
    shell:
        "mosdepth -n -t {threads} -F {params.flag} -Q {params.mapq} -b {params.bed:q} {params.prefix:q} {input:q} > {log:q} 2>&1"


rule mosdepth_targets:
    input:
        str(R / "alignments" / "{sample}.bam"),
    output:
        regions=str(R / "qc" / "samples" / "{sample}" / "targets.regions.bed.gz"),
        ontarget=str(R / "qc" / "samples" / "{sample}" / "{sample}.ontarget.txt"),
    params:
        prefix=lambda wc: str(R / "qc" / "samples" / wc.sample / "targets"),
        bed=TARGETS,
        flag=P["layer1.mosdepth_exclude_flag"],
        # Same mapping-quality floor as the core-genome depth and the caller: target depth
        # counts only reads the caller can use. The on-target read count stays unfiltered.
        mapq=P["layer2.min_mapping_quality"],
    threads: threads(2)
    resources:
        mem_mb=mem(1000),
        runtime=RUNTIME,
    log:
        str(LOGS / "mosdepth" / "{sample}.targets.log"),
    benchmark:
        str(BENCH / "mosdepth" / "{sample}.targets.tsv")
    shell:
        "(mosdepth -n -t {threads} -F {params.flag} -Q {params.mapq} -b {params.bed:q} {params.prefix:q} {input:q} "
        "&& samtools view -c -F {params.flag} -L {params.bed:q} {input:q} > {output.ontarget:q}) > {log:q} 2>&1"


checkpoint l1_gate:
    input:
        unpack(gate_inputs),
    output:
        str(R / "qc" / "samples" / "{sample}" / "qc_decision.json"),
    params:
        frozen=FROZEN,
        root=str(ROOT),
    benchmark:
        str(BENCH / "l1_gate" / "{sample}.tsv")
    log:
        str(LOGS / "gate" / "{sample}.log"),
    script:
        "../scripts/l1_gate.py"
