# L1.1 FASTQ QC and trimming with fastp (per sample and lane). Raw FASTQ is read-only.


def lane(wc):
    return next(x for x in SAMPLES[wc.sample]["lanes"] if x["lane"] == wc.lane)


def fastp_args():
    a = []
    if P["layer1.fastp_detect_adapter_for_pe"]:
        a.append("--detect_adapter_for_pe")
    a += [f"--qualified_quality_phred {P['layer1.fastp_qualified_quality_phred']}",
          f"--unqualified_percent_limit {P['layer1.fastp_unqualified_percent_limit']}",
          f"--n_base_limit {P['layer1.fastp_n_base_limit']}",
          f"--length_required {P['layer1.fastp_length_required']}"]
    if P["layer1.fastp_cut_right"]:
        a.append("--cut_right --cut_window_size 4 --cut_mean_quality 20")
    if CFG["assay"]["two_colour_chemistry"]:
        a.append("--trim_poly_g")
    else:
        a.append("--disable_trim_poly_g")
    if DESIGN == "TARGETED_UMI":
        # UMIs are inline at the 5' end: no front or quality trimming before extraction.
        a = [x for x in a if not x.startswith("--cut_right")] + ["--disable_quality_filtering",
                                                                   "--disable_length_filtering"]
    return " ".join(a)


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
        args=fastp_args(),
    threads: threads(4)
    resources:
        mem_mb=mem(2000),
        runtime=RUNTIME,
    log:
        str(LOGS / "fastp" / "{sample}.{lane}.log"),
    benchmark:
        str(BENCH / "fastp" / "{sample}.{lane}.tsv")
    shell:
        "fastp -i {input.r1:q} -I {input.r2:q} -o {output.r1:q} -O {output.r2:q} "
        "-j {output.json:q} -h {output.html:q} -w {threads} {params.args} > {log:q} 2>&1"
