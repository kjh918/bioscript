# [METADATA]
# TOOL_NAME = fastqc
# VERSION = 0.12.1
# THREADS = 8

rule fastqc:
    input:
        SeqID = ""
        RawFastqDir = ""
        qcResDir = ""
    output:
        r1_html = "[qcResDir]/[SeqID]_R1_fastqc.html"
        r2_html = "[qcResDir]/[SeqID]_R2_fastqc.html"
        r1_zip = "[qcResDir]/[SeqID]_R1_fastqc.zip"
        r2_zip = "[qcResDir]/[SeqID]_R2_fastqc.zip"
        r1_dir = "[qcResDir]/[SeqID]_R1_fastqc"
        r2_dir = "[qcResDir]/[SeqID]_R2_fastqc"
    params:
        singularity_bin = "singularity"
        fastqc_bin = "fastqc"
        sif = "/storage/images/fastqc-0.12.1.sif"
        bind = "/storage,/data"
        threads = "8"
        fastqc_args = "--extract"
    threads: 8
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.fastqc_bin} {params.fastqc_args} --threads {threads} --outdir {input.qcResDir} {input.RawFastqDir}/{input.SeqID}_R1.fastq.gz {input.RawFastqDir}/{input.SeqID}_R2.fastq.gz
        """