# L1.1 FASTQ QC and trimming with fastp (per sample and lane). Raw FASTQ is read-only.


rule fastp:
    input:
        r1=lambda wc: str(ROOT / lane(wc)["fastq_1"]),
        r2=lambda wc: str(ROOT / lane(wc)["fastq_2"]),
    output:
        r1=temp(str(R / "work" / "trimmed" / "{sample}.{lane}.R1.fastq.gz")),
        r2=temp(str(R / "work" / "trimmed" / "{sample}.{lane}.R2.fastq.gz")),
        json=str(R / "qc" / "samples" / "{sample}" / "{sample}.{lane}.fastp.json"),
        html=str(R / "qc" / "samples" / "{sample}" / "{sample}.{lane}.fastp.html"),
    params:
        args=lambda wc: fastp_args(wc.sample),
    threads: threads("fastp")
    resources:
        mem_mb=mem("fastp"),
        runtime=RUNTIME,
    log:
        str(LOGS / "fastp" / "{sample}.{lane}.log"),
    benchmark:
        str(BENCH / "fastp" / "{sample}.{lane}.tsv")
    shell:
        "fastp -i {input.r1:q} -I {input.r2:q} -o {output.r1:q} -O {output.r2:q} "
        "-j {output.json:q} -h {output.html:q} -w {threads} {params.args} > {log:q} 2>&1"
