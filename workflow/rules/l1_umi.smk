# L1.6-L1.8: targeted sequencing with UMIs (fgbio). Order is fixed:
# FastqToBam -> align -> ZipperBams -> template-coordinate sort -> GroupReadsByUmi ->
# CallMolecularConsensusReads -> FilterConsensusReads -> realign -> ZipperBams ->
# coordinate sort -> pre-clipping metrics -> TrimPrimers -> index.

if DESIGN == "TARGETED_UMI" and PLAN["input_type"] == "FASTQ":

    rule umi_fastq_to_bam:
        input:
            r1=str(R / "work" / "trimmed" / "{sample}.{lane}.R1.fastq.gz"),
            r2=str(R / "work" / "trimmed" / "{sample}.{lane}.R2.fastq.gz"),
        output:
            temp(str(R / "work" / "umi" / "lanes" / "{sample}.{lane}.unmapped.bam")),
        params:
            rs=CFG["assay"]["read_structure"].split(),
            lib=lambda wc: next(x for x in SAMPLES[wc.sample]["lanes"] if x["lane"] == wc.lane)["library_id"],
        resources:
            mem_mb=mem(2000),
            runtime=RUNTIME,
        log:
            str(LOGS / "umi" / "{sample}.{lane}.fastqtobam.log"),
        shell:
            "fgbio -Xmx{resources.mem_mb}m FastqToBam -i {input.r1:q} {input.r2:q} -r {params.rs:q} "
            "-o {output:q} --sample {wildcards.sample} --library {params.lib:q} "
            "--read-group-id {wildcards.sample}.{wildcards.lane} --sort=true > {log:q} 2>&1"

    rule umi_align_zipper:
        input:
            str(R / "work" / "umi" / "lanes" / "{sample}.{lane}.unmapped.bam"),
        output:
            temp(str(R / "work" / "umi" / "lanes" / "{sample}.{lane}.mapped.bam")),
        params:
            ref=REF,
            K=P["layer1.bwa_mem2_fixed_batch_K"],
        threads: threads(4)
        resources:
            mem_mb=mem(6000),
            runtime=RUNTIME,
        log:
            str(LOGS / "umi" / "{sample}.{lane}.align.log"),
        shell:
            "(samtools fastq {input:q} | bwa-mem2 mem -p -t {threads} -K {params.K} -Y {params.ref:q} - "
            "| fgbio -Xmx2g --compression 0 ZipperBams -i /dev/stdin --unmapped {input:q} --ref {params.ref:q} "
            "--output {output:q}) 2> {log:q}"

    rule umi_group:
        input:
            lambda wc: [str(R / "work" / "umi" / "lanes" / f"{wc.sample}.{x['lane']}.mapped.bam") for x in SAMPLES[wc.sample]["lanes"]],
        output:
            bam=temp(str(R / "work" / "umi" / "{sample}.grouped.bam")),
            hist=str(R / "qc" / "samples" / "{sample}" / "{sample}.umi_family_sizes.txt"),
            metrics=str(R / "qc" / "samples" / "{sample}" / "{sample}.umi_grouping_metrics.txt"),
        params:
            strategy=P["layer1.fgbio_group_strategy"],
            edits=P["layer1.fgbio_group_edits"],
            mapq=P["layer1.fgbio_group_min_map_q"],
        threads: threads(2)
        resources:
            mem_mb=mem(4000),
            runtime=RUNTIME,
        log:
            str(LOGS / "umi" / "{sample}.group.log"),
        shell:
            "(samtools merge -f -u -o - {input:q} | samtools sort --template-coordinate -@ {threads} -o {output.bam:q}.tc.bam - "
            "&& fgbio -Xmx{resources.mem_mb}m GroupReadsByUmi -i {output.bam:q}.tc.bam -o {output.bam:q} "
            "-s {params.strategy} -e {params.edits} -m {params.mapq} -f {output.hist:q} -g {output.metrics:q} "
            "&& rm {output.bam:q}.tc.bam) > {log:q} 2>&1"

    rule umi_consensus:
        input:
            str(R / "work" / "umi" / "{sample}.grouped.bam"),
        output:
            temp(str(R / "work" / "umi" / "{sample}.consensus.unmapped.bam")),
        params:
            minr=P["layer1.fgbio_consensus_min_reads"],
            pre=P["layer1.fgbio_consensus_error_pre_umi"],
            post=P["layer1.fgbio_consensus_error_post_umi"],
            minq=P["layer1.fgbio_consensus_min_input_base_quality"],
            ref=REF,
            fminr=P["layer1.fgbio_filter_min_reads"],
            fminq=P["layer1.fgbio_filter_min_base_quality"],
            E=P["layer1.fgbio_filter_max_read_error_rate"],
            e=P["layer1.fgbio_filter_max_base_error_rate"],
            nc=P["layer1.fgbio_filter_max_no_call_fraction"],
        threads: threads(2)
        resources:
            mem_mb=mem(4000),
            runtime=RUNTIME,
        log:
            str(LOGS / "umi" / "{sample}.consensus.log"),
        shell:
            "(fgbio -Xmx{resources.mem_mb}m CallMolecularConsensusReads -i {input:q} -o {output:q}.raw.bam "
            "-M {params.minr} --error-rate-pre-umi={params.pre} --error-rate-post-umi={params.post} "
            "-m {params.minq} --threads={threads} "
            "&& fgbio -Xmx{resources.mem_mb}m FilterConsensusReads -i {output:q}.raw.bam -o {output:q} -r {params.ref:q} "
            "-M {params.fminr} -N {params.fminq} -E {params.E} -e {params.e} -n {params.nc} "
            "&& rm {output:q}.raw.bam) > {log:q} 2>&1"

    rule umi_realign:
        input:
            str(R / "work" / "umi" / "{sample}.consensus.unmapped.bam"),
        output:
            temp(str(R / "work" / "umi" / "{sample}.consensus.mapped.bam")),
        params:
            ref=REF,
            K=P["layer1.bwa_mem2_fixed_batch_K"],
        threads: threads(4)
        resources:
            mem_mb=mem(6000),
            runtime=RUNTIME,
        log:
            str(LOGS / "umi" / "{sample}.realign.log"),
        shell:
            "(samtools sort -n -u -o {output:q}.qn.bam {input:q} && samtools fastq {output:q}.qn.bam "
            "| bwa-mem2 mem -p -t {threads} -K {params.K} -Y {params.ref:q} - "
            "| fgbio -Xmx2g --compression 0 ZipperBams -i /dev/stdin --unmapped {output:q}.qn.bam --ref {params.ref:q} "
            "--tags-to-reverse Consensus --tags-to-revcomp Consensus "
            "| samtools sort -@ {threads} -o {output:q} - && rm {output:q}.qn.bam) 2> {log:q}"

    rule umi_preclip_and_trim:
        input:
            str(R / "work" / "umi" / "{sample}.consensus.mapped.bam"),
        output:
            bam=str(R / "alignments" / "{sample}.bam"),
            stats=str(R / "qc" / "samples" / "{sample}" / "{sample}.preclip.stats.txt"),
        params:
            primers=str(ROOT / CFG["assay"]["primers_tsv"]) if CFG["assay"]["primers_tsv"] else "",
            slop=P["layer1.fgbio_trimprimers_slop"],
            ref=REF,
        resources:
            mem_mb=mem(3000),
            runtime=RUNTIME,
        log:
            str(LOGS / "umi" / "{sample}.trim.log"),
        shell:
            "(samtools stats {input:q} > {output.stats:q} && fgbio -Xmx{resources.mem_mb}m TrimPrimers -i {input:q} "
            "-o {output.bam:q} -p {params.primers:q} -S {params.slop} -r {params.ref:q} --first-of-pair=false "
            "&& samtools index {output.bam:q}) > {log:q} 2>&1"
