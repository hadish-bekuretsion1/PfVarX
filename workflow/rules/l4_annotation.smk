# Layer 4: SnpEff with a database built from the pinned reference and GFF; mismatch is FAIL.

from pfvarx.annotation.snpeff_db import genome_id

GENOME_ID = genome_id(CFG["reference"].get("snpeff_genome_id"), CFG["reference"]["release_label"])
SNPEFF_DIR = ROOT / "prepared" / "snpeff"


rule snpeff_build:
    input:
        fasta=PREF,
        gff=str(ROOT / CFG["reference"]["gff"]) if CFG["reference"].get("gff") else [],
    output:
        db=str(SNPEFF_DIR / "data" / GENOME_ID / "snpEffectPredictor.bin"),
        check=str(SNPEFF_DIR / f"{GENOME_ID}.build_check.json"),
    params:
        genome=GENOME_ID,
        data=str(SNPEFF_DIR),
    resources:
        mem_mb=mem(4000),
        runtime=RUNTIME,
    benchmark:
        str(BENCH / "snpeff_build" / "snpeff_build.tsv")
    log:
        str(LOGS / "snpeff" / "build.log"),
    shell:
        "python -m pfvarx.annotation.snpeff_db --fasta {input.fasta:q} --gff {input.gff:q} --genome {params.genome:q} "
        "--data-dir {params.data:q} --check {output.check:q} > {log:q} 2>&1"


rule snpeff_annotate:
    input:
        vcf=str(R / "cohort" / "cohort.vcf.gz"),
        db=str(SNPEFF_DIR / "data" / GENOME_ID / "snpEffectPredictor.bin"),
        check=str(SNPEFF_DIR / f"{GENOME_ID}.build_check.json"),
    output:
        vcf=str(R / "annotation" / "cohort.ann.vcf.gz"),
        tbi=str(R / "annotation" / "cohort.ann.vcf.gz.tbi"),
        qc=str(R / "annotation" / "annotation_qc.json"),
    params:
        ud=P["layer4.snpeff_upstream_downstream_len"],
        genome=GENOME_ID,
        data=str(SNPEFF_DIR),
    resources:
        mem_mb=mem(4000),
        runtime=RUNTIME,
    benchmark:
        str(BENCH / "snpeff_annotate" / "snpeff_annotate.tsv")
    log:
        str(LOGS / "snpeff" / "annotate.log"),
    shell:
        "python -m pfvarx.annotation.snpeff_db --annotate --vcf {input.vcf:q} --genome {params.genome:q} "
        "--data-dir {params.data:q} --check {input.check:q} --ud {params.ud} --out-vcf {output.vcf:q} "
        "--out-qc {output.qc:q} > {log:q} 2>&1"
