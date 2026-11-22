# L1.2-L1.5, L1.9: alignment, duplicate marking (WGS), primer clipping (targeted, no UMI)
# and existing-alignment entry.


# Human DNA removal (reference.host_removal: true, wgs only): reads are aligned to the
# parasite+host reference and the host filter removes pairs aligned to the host before
# sorting (counts in <sample>.<lane>.host.json).
HOST_COUNTS = str(R / "qc" / "samples" / "{sample}" / "{sample}.{lane}.host.json")


rule bwa_mem2_lane:
    input:
        r1=str(R / "work" / "trimmed" / "{sample}.{lane}.R1.fastq.gz"),
        r2=str(R / "work" / "trimmed" / "{sample}.{lane}.R2.fastq.gz"),
    output:
        bam=temp(str(R / "work" / "lanes" / "{sample}.{lane}.bam")),
        **({"host": HOST_COUNTS} if HOST else {}),
    params:
        rg=rg,
        ref=ALIGN_REF,
        K=P["layer1.bwa_mem2_fixed_batch_K"],
        host=host_filter(HOST_COUNTS),
    threads: threads(4)
    resources:
        mem_mb=mem(6000, host=True),
        runtime=RUNTIME,
    log:
        str(LOGS / "bwa_mem2" / "{sample}.{lane}.log"),
    benchmark:
        str(BENCH / "bwa_mem2" / "{sample}.{lane}.tsv")
    shell:
        "(bwa-mem2 mem -t {threads} -K {params.K} -R {params.rg:q} {params.ref:q} {input.r1:q} {input.r2:q} "
        "{params.host}| samtools sort -@ {threads} -o {output.bam:q} -) 2> {log:q}"


rule merge_lanes:
    input:
        lambda wc: [str(R / "work" / "lanes" / f"{wc.sample}.{x['lane']}.bam") for x in SAMPLES[wc.sample]["lanes"]],
    output:
        temp(str(R / "work" / "merged" / "{sample}.bam")),
    threads: threads(2)
    resources:
        mem_mb=mem(1000),
        runtime=RUNTIME,
    benchmark:
        str(BENCH / "merge_lanes" / "{sample}.tsv")
    log:
        str(LOGS / "merge" / "{sample}.log"),
    shell:
        "samtools merge -f -@ {threads} -o {output:q} {input:q} 2> {log:q}"


if DESIGN == "wgs" and PLAN["input_type"] == "FASTQ":

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


# Primer trimming (fgbio TrimPrimers) when the panel defines primers; a probe-only (capture)
# panel has nothing to trim, so the reads pass through unchanged.
TRIM_PRIMERS = (
    "if [ -n {params.primers:q} ]; then "
    "fgbio -Djava.io.tmpdir={resources.tmpdir:q} -Xmx{resources.mem_mb}m TrimPrimers -i {input.bam:q} "
    "-o {output.bam:q} -p {params.primers:q} -S {params.slop} -r {params.ref:q} --first-of-pair=false; "
    "else cp {input.bam:q} {output.bam:q}; fi"
)


if DESIGN == "targeted_no_umi" and PLAN["input_type"] == "FASTQ":

    rule insert_metrics_preclip:
        input:
            str(R / "work" / "merged" / "{sample}.bam"),
        output:
            str(R / "qc" / "samples" / "{sample}" / "{sample}.preclip.stats.txt"),
        benchmark:
            str(BENCH / "insert_metrics_preclip" / "{sample}.tsv")
        log:
            str(LOGS / "preclip" / "{sample}.log"),
        shell:
            "samtools stats {input:q} > {output:q} 2> {log:q}"

    rule trim_primers:
        input:
            bam=str(R / "work" / "merged" / "{sample}.bam"),
            stats=str(R / "qc" / "samples" / "{sample}" / "{sample}.preclip.stats.txt"),
        output:
            bam=str(R / "alignments" / "{sample}.bam"),
        params:
            primers=PRIMERS,
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
            "(" + TRIM_PRIMERS + " && samtools index {output.bam:q}) > {log:q} 2>&1"


if PLAN["input_type"] in ("BAM", "CRAM"):

    rule alignment_entry:
        """L1.9: re-check an imported alignment and do what its documented processing state
        still needs before trimming or consensus (WGS: duplicate marking; raw UMI reads: mate
        tags for UMI grouping). The state was checked against the reads at validation."""
        input:
            lambda wc: str(ROOT / SAMPLES[wc.sample]["bam"]),
        output:
            bam=temp(str(R / "work" / "entry" / "{sample}.bam")),
            contract=str(R / "qc" / "samples" / "{sample}" / "{sample}.entry_contract.json"),
        params:
            ref=REF,
            design=DESIGN,
            state=lambda wc: library(wc.sample).get("processing_state") or "",
            steps=lambda wc: ",".join(entry_steps(wc.sample)),
        resources:
            mem_mb=mem(4000),
            runtime=RUNTIME,
        benchmark:
            str(BENCH / "alignment_entry" / "{sample}.tsv")
        log:
            str(LOGS / "entry" / "{sample}.log"),
        shell:
            "python -m pfvarx.inputs.alignment_entry --input {input:q} --output {output.bam:q} "
            "--contract {output.contract:q} --reference {params.ref:q} --design {params.design} "
            "--state={params.state} --steps={params.steps} > {log:q} 2>&1"

    rule entry_trim_primers:
        """Primer trimming for imported alignments documented as not yet trimmed."""
        input:
            bam=entry_trim_source,
        output:
            bam=temp(str(R / "work" / "entry" / "{sample}.trimmed.bam")),
            stats=str(R / "qc" / "samples" / "{sample}" / "{sample}.preclip.stats.txt"),
        params:
            primers=PRIMERS,
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
            "(samtools stats {input.bam:q} > {output.stats:q} && " + TRIM_PRIMERS + ") > {log:q} 2>&1"

    rule entry_finalize:
        """The caller-ready alignment: the last step the documented state needed, coordinate-sorted
        and indexed (trimming moves alignment starts)."""
        input:
            entry_final_source,
        output:
            str(R / "alignments" / "{sample}.bam"),
        benchmark:
            str(BENCH / "entry_finalize" / "{sample}.tsv")
        log:
            str(LOGS / "entry" / "{sample}.finalize.log"),
        shell:
            "(samtools sort -@ 1 -o {output:q}.sorted.bam {input:q} && mv {output:q}.sorted.bam {output:q} "
            "&& samtools index {output:q}) > {log:q} 2>&1"
