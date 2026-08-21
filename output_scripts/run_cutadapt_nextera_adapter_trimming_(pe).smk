# [METADATA]
# TOOL_NAME = Cutadapt
# VERSION = 4.4
# THREADS = 4

rule cutadapt:
    input:
        SeqID = ""
        RawFastqDir = ""
        OutDir = ""
    output:
        r1_clean = "[OutDir]/[SeqID]_R1.trimmed.fastq.gz"
        r2_clean = "[OutDir]/[SeqID]_R2.trimmed.fastq.gz"
        trim_log = "[OutDir]/[SeqID]_cutadapt.log"
    params:
        InputSuffix = ".fastq.gz"
        AdapterR1 = "CTGTCTCTTATACACATCT"
        AdapterR2 = "CTGTCTCTTATACACATCT"
        QualityCutoff = "20"
        MinLength = "30"
        Threads = "4"
        cutadapt_bin = "/storage/home/jhkim/Apps/Python-3.11.13/python -m cutadapt"
    threads: 4
    shell:
        """
        {params.cutadapt_bin} -j {threads} -a {params.AdapterR1} -A {params.AdapterR2} -q {params.QualityCutoff} -m {params.MinLength} -o {output.r1_clean} -p {output.r2_clean} {input.RawFastqDir}/{input.SeqID}_R1{params.InputSuffix} {input.RawFastqDir}/{input.SeqID}_R2{params.InputSuffix} > {output.trim_log}
        """