# [METADATA]
# TOOL_NAME = samtools
# VERSION = 1.10
# THREADS = 1

rule samtools:
    input:
        SeqID = ""
        BamDir = ""
    output:
        sorted_bam = "[BamDir]/[SeqID].[InputSuffix].bam"
        sorted_bai = "[BamDir]/[SeqID].[OutputSuffix].bam.bai"
    params:
        InputSuffix = "primary"
        OutputSuffix = "sorted"
        samtools_bin = "samtools"
        Threads = "8"
    threads: 1
    shell:
        """
        {params.samtools_bin} sort -@ {threads} -o {output.sorted_bam} {input.BamDir}/{input.SeqID}.bam && {params.samtools_bin} index -b -@ {threads} {output.sorted_bam}
        """