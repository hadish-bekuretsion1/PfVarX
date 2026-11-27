# Layer 3: variant-only merge (NOT joint genotyping) and the cohort evidence matrix.


rule cohort_merge:
    input:
        vcfs=lambda wc: [str(R / "variants" / "samples" / f"{s}.vcf.gz") for s in passing("L3")],
    output:
        vcf=str(R / "cohort" / "cohort.vcf.gz"),
        tbi=str(R / "cohort" / "cohort.vcf.gz.tbi"),
        samples=str(R / "cohort" / "cohort_samples.txt"),
    params:
        # Sample IDs in merge order (IDs never contain dots, so the file stem is the ID).
        ids=lambda wc, input: " ".join(Path(v).name.split(".")[0] for v in input.vcfs),
        listing=str(R / "work" / "cohort_merge.vcfs.txt"),
        headroom=COHORT_MERGE_HEADROOM,
    benchmark:
        str(BENCH / "cohort_merge" / "cohort_merge.tsv")
    log:
        str(LOGS / "cohort" / "merge.log"),
    # final_targets requests the cohort only when a sample is permitted into Layer 3. One
    # sample is copied (bcftools merge needs two files); several are merged variant-only, with
    # missing genotypes left missing (no -0: absence of a call is not a reference call).
    # printf is a shell builtin, so neither list is limited by the command-line length; the
    # merge raises the open-file limit itself and merges in batches beyond it (cohort/merge.py).
    shell:
        "(printf '%s\\n' {params.ids} > {output.samples:q} "
        "&& mkdir -p $(dirname {params.listing:q}) && printf '%s\\n' {input.vcfs:q} > {params.listing:q} "
        "&& python -m pfvarx.cohort.merge --list {params.listing:q} --out {output.vcf:q} "
        "--headroom {params.headroom}) > {log:q} 2>&1"


rule cohort_status:
    """Record whether a cohort was formed; with none, mark requested analyses as skipped."""
    input:
        unpack(cohort_status_inputs),
        gates=lambda wc: [str(R / "qc" / "samples" / s / "qc_decision.json") for s in sorted(SAMPLES) if s not in EXCLUDED],
    output:
        str(R / "cohort" / "cohort_status.json"),
    params:
        empty=lambda wc: "--empty" if not passing("L3") else "",
        run=str(RUN),
        root=str(ROOT),
    benchmark:
        str(BENCH / "cohort_status" / "cohort_status.tsv")
    log:
        str(LOGS / "cohort" / "status.log"),
    shell:
        "python -m pfvarx.downstream.run --cohort-status {params.empty} --run-dir {params.run:q} "
        "--root {params.root:q} > {log:q} 2>&1"


rule cohort_evidence:
    input:
        vcf=str(R / "cohort" / "cohort.vcf.gz"),
        samples=str(R / "cohort" / "cohort_samples.txt"),
        bams=lambda wc: [str(R / "alignments" / f"{s}.bam") for s in passing("L3")],
        callable=lambda wc: [str(R / "variants" / "samples" / f"{s}.callable.bed") for s in passing("L3")],
    output:
        matrix=str(R / "cohort" / "evidence_matrix.tsv.gz"),
        missing=str(R / "cohort" / "missingness.tsv"),
        sites=str(R / "cohort" / "sites.tsv"),
    params:
        frozen=FROZEN,
        ref=PREF,
        bam_list=str(R / "work" / "cohort_evidence.bams.txt"),
        callable_list=str(R / "work" / "cohort_evidence.callable.txt"),
    # Samples are independent: they run in worker processes and are written in cohort order,
    # so the matrix does not depend on the thread count. The file lists are written with the
    # printf builtin and read from files, so cohort size is not limited by the command line.
    threads: threads("cohort_evidence")
    resources:
        mem_mb=mem("cohort_evidence"),
        runtime=RUNTIME,
    benchmark:
        str(BENCH / "cohort_evidence" / "cohort_evidence.tsv")
    log:
        str(LOGS / "cohort" / "evidence.log"),
    shell:
        "(mkdir -p $(dirname {params.bam_list:q}) && printf '%s\\n' {input.bams:q} > {params.bam_list:q} "
        "&& printf '%s\\n' {input.callable:q} > {params.callable_list:q} "
        "&& python -m pfvarx.cohort.evidence --vcf {input.vcf:q} --samples {input.samples:q} "
        "--bam-list {params.bam_list:q} --callable-list {params.callable_list:q} --threads {threads} "
        "--reference {params.ref:q} --frozen-config {params.frozen:q} "
        "--out-matrix {output.matrix:q} --out-missing {output.missing:q} --out-sites {output.sites:q}) > {log:q} 2>&1"


rule cohort_microhaplotypes:
    """Multiplexed microhaplotypes at catalogue loci: one multiallelic marker per locus for
    DA001 and DA004. Writes an empty artifact (with the reason) when no verified panel is
    configured, so the cohort step always has a recorded outcome."""
    input:
        sites=str(R / "cohort" / "sites.tsv"),
        samples=str(R / "cohort" / "cohort_samples.txt"),
        bams=lambda wc: [str(R / "alignments" / f"{s}.bam") for s in passing("L3")],
        callable=lambda wc: [str(R / "variants" / "samples" / f"{s}.callable.bed") for s in passing("L3")],
    output:
        haplotypes=str(R / "cohort" / "microhaplotypes.tsv.gz"),
        loci=str(R / "cohort" / "microhaplotype_loci.tsv"),
        summary=str(R / "cohort" / "microhaplotypes.json"),
    params:
        frozen=FROZEN,
        ref=PREF,
        root=str(ROOT),
        bam_list=str(R / "work" / "cohort_microhaplotypes.bams.txt"),
        callable_list=str(R / "work" / "cohort_microhaplotypes.callable.txt"),
    threads: threads("cohort_microhaplotypes")
    resources:
        mem_mb=mem("cohort_microhaplotypes"),
        runtime=RUNTIME,
    benchmark:
        str(BENCH / "cohort_microhaplotypes" / "cohort_microhaplotypes.tsv")
    log:
        str(LOGS / "cohort" / "microhaplotypes.log"),
    shell:
        "(mkdir -p $(dirname {params.bam_list:q}) && printf '%s\\n' {input.bams:q} > {params.bam_list:q} "
        "&& printf '%s\\n' {input.callable:q} > {params.callable_list:q} "
        "&& python -m pfvarx.cohort.microhap --sites {input.sites:q} --samples {input.samples:q} "
        "--bam-list {params.bam_list:q} --callable-list {params.callable_list:q} --threads {threads} "
        "--reference {params.ref:q} --root {params.root:q} --frozen-config {params.frozen:q} "
        "--out-haplotypes {output.haplotypes:q} --out-loci {output.loci:q} "
        "--out-summary {output.summary:q}) > {log:q} 2>&1"


rule cohort_qc:
    """Cohort QC (gate G5): each sample's completeness in the cohort; a FAIL excludes the
    sample from the downstream analyses."""
    input:
        missing=str(R / "cohort" / "missingness.tsv"),
        matrix=str(R / "cohort" / "evidence_matrix.tsv.gz"),
    output:
        str(R / "cohort" / "cohort_qc.json"),
    params:
        frozen=FROZEN,
        root=str(ROOT),
    benchmark:
        str(BENCH / "cohort_qc" / "cohort_qc.tsv")
    log:
        str(LOGS / "cohort" / "qc.log"),
    shell:
        "python -m pfvarx.quality.cohort --root {params.root:q} --frozen-config {params.frozen:q} "
        "--missingness {input.missing:q} --matrix {input.matrix:q} --out {output:q} > {log:q} 2>&1"
