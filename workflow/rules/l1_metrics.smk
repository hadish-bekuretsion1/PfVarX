# L1 QC metrics and the per-sample Layer 1 gate (checkpoint).

TARGETS = str(ROOT / CFG["assay"]["targets_bed"]) if CFG["assay"]["targets_bed"] else ""


rule flagstat:
    input:
        str(R / "alignments" / "{sample}.bam"),
    output:
        str(R / "qc" / "samples" / "{sample}" / "{sample}.flagstat.json"),
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
    log:
        str(LOGS / "mosdepth" / "{sample}.targets.log"),
    shell:
        "(mosdepth -n -F {params.flag} -b {params.bed:q} {params.prefix:q} {input:q} "
        "&& samtools view -c -F {params.flag} -L {params.bed:q} {input:q} > {output.ontarget:q}) > {log:q} 2>&1"


def gate_inputs(wc):
    d = R / "qc" / "samples" / wc.sample
    files = {"flagstat": str(d / f"{wc.sample}.flagstat.json")}
    if PLAN["input_type"] == "FASTQ":
        files["fastp"] = [str(d / f"{wc.sample}.{x['lane']}.fastp.json") for x in SAMPLES[wc.sample]["lanes"]]
    if DESIGN == "WGS":
        files["mosdepth"] = str(d / "core.mosdepth.region.dist.txt")
        if PLAN["input_type"] == "FASTQ":
            files["markdup"] = str(d / f"{wc.sample}.markdup.txt")
    else:
        files["targets"] = str(d / "targets.regions.bed.gz")
        files["ontarget"] = str(d / f"{wc.sample}.ontarget.txt")
    if DESIGN == "TARGETED_UMI" and PLAN["input_type"] == "FASTQ":
        files["umi"] = str(d / f"{wc.sample}.umi_family_sizes.txt")
    return files


checkpoint l1_gate:
    input:
        unpack(gate_inputs),
    output:
        str(R / "qc" / "samples" / "{sample}" / "qc_decision.json"),
    params:
        frozen=FROZEN,
        root=str(ROOT),
    log:
        str(LOGS / "gate" / "{sample}.log"),
    script:
        "../scripts/l1_gate.py"
