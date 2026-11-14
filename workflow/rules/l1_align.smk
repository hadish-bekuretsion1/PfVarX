# L1.2-L1.5, L1.9: alignment, duplicate marking (WGS), primer clipping (targeted, no UMI)
# and existing-alignment entry.


def rg(wc):
    s = SAMPLES[wc.sample]
    ln = next(x for x in s["lanes"] if x["lane"] == wc.lane)
    return (f"@RG\\tID:{wc.sample}.{wc.lane}\\tSM:{wc.sample}\\tLB:{ln['library_id']}"
            f"\\tPL:{CFG['assay']['platform']}\\tPU:{wc.sample}.{wc.lane}")


rule bwa_mem2_lane:
    input:
        r1=str(R / "work" / "trimmed" / "{sample}.{lane}.R1.fastq.gz"),
        r2=str(R / "work" / "trimmed" / "{sample}.{lane}.R2.fastq.gz"),
    output:
        temp(str(R / "work" / "lanes" / "{sample}.{lane}.bam")),
    params:
        rg=rg,
        ref=REF,
        K=P["layer1.bwa_mem2_fixed_batch_K"],
    threads: threads(4)
    resources:
        mem_mb=mem(6000),
        runtime=RUNTIME,
    log:
        str(LOGS / "bwa_mem2" / "{sample}.{lane}.log"),
    benchmark:
        str(BENCH / "bwa_mem2" / "{sample}.{lane}.tsv")
    shell:
        "(bwa-mem2 mem -t {threads} -K {params.K} -R {params.rg:q} {params.ref:q} {input.r1:q} {input.r2:q} "
        "| samtools sort -@ {threads} -o {output:q} -) 2> {log:q}"


rule merge_lanes:
    input:
        lambda wc: [str(R / "work" / "lanes" / f"{wc.sample}.{x['lane']}.bam") for x in SAMPLES[wc.sample]["lanes"]],
    output:
        temp(str(R / "work" / "merged" / "{sample}.bam")),
    threads: threads(2)
    resources:
        mem_mb=mem(1000),
        runtime=RUNTIME,
    log:
        str(LOGS / "merge" / "{sample}.log"),
    shell:
        "samtools merge -f -@ {threads} -o {output:q} {input:q} 2> {log:q}"


if DESIGN == "WGS" and PLAN["input_type"] == "FASTQ":

    rule mark_duplicates:
        input:
            str(R / "work" / "merged" / "{sample}.bam"),
        output:
            bam=str(R / "alignments" / "{sample}.bam"),
            metrics=str(R / "qc" / "samples" / "{sample}" / "{sample}.markdup.txt"),
        params:
            pixel=lambda wc: (P["layer1.picard_optical_pixel_distance_patterned"]
                              if SAMPLES[wc.sample]["lanes"][0].get("flowcell_type") == "patterned"
                              else P["layer1.picard_optical_pixel_distance_unpatterned"]),
        threads: 1
        resources:
            mem_mb=mem(4000),
            runtime=RUNTIME,
        retries: COMPUTE["max_retries"]
        log:
            str(LOGS / "markdup" / "{sample}.log"),
        benchmark:
            str(BENCH / "markdup" / "{sample}.tsv")
        shell:
            "picard -Djava.io.tmpdir={resources.tmpdir:q} -Xmx{resources.mem_mb}m MarkDuplicates --TMP_DIR {resources.tmpdir:q} -I {input:q} -O {output.bam:q} -M {output.metrics:q} "
            "--OPTICAL_DUPLICATE_PIXEL_DISTANCE {params.pixel} --VALIDATION_STRINGENCY LENIENT > {log:q} 2>&1 "
            "&& samtools index {output.bam:q} 2>> {log:q}"


if DESIGN == "TARGETED_NON_UMI" and PLAN["input_type"] == "FASTQ":

    rule insert_metrics_preclip:
        input:
            str(R / "work" / "merged" / "{sample}.bam"),
        output:
            str(R / "qc" / "samples" / "{sample}" / "{sample}.preclip.stats.txt"),
        log:
            str(LOGS / "preclip" / "{sample}.log"),
        shell:
            "samtools stats {input:q} > {output:q} 2> {log:q}"

    rule trim_primers:
        input:
            bam=str(R / "work" / "merged" / "{sample}.bam"),
            stats=str(R / "qc" / "samples" / "{sample}" / "{sample}.preclip.stats.txt"),
        output:
            str(R / "alignments" / "{sample}.bam"),
        params:
            primers=str(ROOT / CFG["assay"]["primers_tsv"]) if CFG["assay"]["primers_tsv"] else "",
            slop=P["layer1.fgbio_trimprimers_slop"],
            ref=REF,
        resources:
            mem_mb=mem(3000),
            runtime=RUNTIME,
        log:
            str(LOGS / "trim_primers" / "{sample}.log"),
        benchmark:
            str(BENCH / "trim_primers" / "{sample}.tsv")
        shell:
            "fgbio -Djava.io.tmpdir={resources.tmpdir:q} -Xmx{resources.mem_mb}m TrimPrimers -i {input.bam:q} -o {output:q} -p {params.primers:q} "
            "-S {params.slop} -r {params.ref:q} --first-of-pair=false > {log:q} 2>&1 "
            "&& samtools index {output:q} 2>> {log:q}"


if PLAN["input_type"] in ("BAM", "CRAM"):

    rule alignment_entry:
        """L1.9: re-use a validated caller-ready alignment; never silently skip missing steps."""
        input:
            lambda wc: str(ROOT / SAMPLES[wc.sample]["bam"]),
        output:
            bam=str(R / "alignments" / "{sample}.bam"),
            contract=str(R / "qc" / "samples" / "{sample}" / "{sample}.entry_contract.json"),
        params:
            ref=REF,
            design=DESIGN,
        resources:
            mem_mb=mem(4000),
            runtime=RUNTIME,
        log:
            str(LOGS / "entry" / "{sample}.log"),
        shell:
            "python -m pfvarx.inputs.alignment_entry --input {input:q} --output {output.bam:q} "
            "--contract {output.contract:q} --reference {params.ref:q} --design {params.design} > {log:q} 2>&1"
