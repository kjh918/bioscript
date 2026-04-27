# [METADATA]
# TOOL_NAME = fastp
# VERSION = 0.23.4
# THREADS = 1

rule fastp:
    input:
        SeqID = ""
        RawFastqDir = ""
        TrimFastqDir = ""
        qcResDir = ""
    output:
        out_read1 = "[TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz"
        out_read2 = "[TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz"
        json = "[qcResDir]/[SeqID].fastp.json"
        html = "[qcResDir]/[SeqID].fastp.html"
    params:
        singularity_bin = "singularity"
        fastp_bin = "fastp"
        sif = "/storage/images/fastp-0.23.4.sif"
        bind = "/storage,/data"
        Threads = "8"
        length_required = "100"
        average_qual = "10"
        qualified_quality_phred = "15"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.fastp_bin} --thread {threads} --in1 {input.RawFastqDir}/{input.SeqID}_R1.fastq.gz --in2 {input.RawFastqDir}/{input.SeqID}_R2.fastq.gz --out1 {output.out_read1} --out2 {output.out_read2} --json {output.json} --html {output.html} --trim_poly_g --detect_adapter_for_pe --length_required {params.length_required} --average_qual {params.average_qual} --qualified_quality_phred {params.qualified_quality_phred}
        """