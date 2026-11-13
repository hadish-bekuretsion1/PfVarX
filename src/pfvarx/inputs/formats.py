"""Filename formats supported by the input transfer interfaces."""

FORMATS = {
    "fasta": (".fasta", ".fa", ".fna"),
    "bed": (".bed",),
    "annotation": (".gff", ".gff3"),
    "metadata": (".tsv",),
    "fastq": (".fastq", ".fq", ".fastq.gz", ".fq.gz"),
    "bam": (".bam",),
    "cram": (".cram",),
}
FOLDER_FORMATS = {
    "reference": FORMATS["fasta"] + FORMATS["bed"] + FORMATS["annotation"],
    "targets": FORMATS["bed"] + FORMATS["metadata"],
    "metadata": FORMATS["metadata"],
    "fastq": FORMATS["fastq"],
    "alignments": FORMATS["bam"] + FORMATS["cram"] + (".bai", ".csi", ".crai"),
    "truth": (".truth.vcf", ".bed", ".json"),
}


def check_extension(name: str, extensions: tuple[str, ...], role: str) -> None:
    if not name.lower().endswith(extensions):
        raise ValueError(f"{name}: {role} accepts {', '.join(extensions)}")
