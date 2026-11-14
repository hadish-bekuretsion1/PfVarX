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
        a.append(f"--cut_right --cut_window_size {P.get('layer1.fastp_cut_window_size', 4)} "
                 f"--cut_mean_quality {P.get('layer1.fastp_cut_mean_quality', 20)}")
    if CFG["assay"]["two_colour_chemistry"]:
        a.append("--trim_poly_g")
    else:
        a.append("--disable_trim_poly_g")
    if DESIGN == "TARGETED_UMI":
        # UMIs are inline at the 5' end: no quality trimming or filtering before extraction.
        # Pairs too short to hold the fixed read-structure segments (adapter dimers) are removed
        # here, because fgbio FastqToBam stops on the first such read.
        from pfvarx.inputs.validate import runtime_read_structures

        _, min_len, trim_safe = runtime_read_structures(CFG["assay"]["read_structure"])
        a = [x for x in a if not x.startswith(("--cut_right", "--length_required"))]
        a += ["--disable_quality_filtering", f"--length_required {min_len}"]
        if not trim_safe:
            # A fixed non-template segment ends the read: 3' trimming would corrupt it.
            a = [x for x in a if x not in ("--detect_adapter_for_pe", "--trim_poly_g", "--disable_trim_poly_g")]
            a += ["--disable_adapter_trimming", "--disable_trim_poly_g"]
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
