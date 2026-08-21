# [METADATA]
# TOOL_NAME = samtools_flagstats
# VERSION = 1.10
# THREADS = 1

rule samtools_flagstats:
    input:
        SeqID = ""
        BamDir = ""
    output:
        flag_stats_txt = "[QcDir]/[SeqID][InputSuffix].[OutputSuffix].txt"
    params:
        InputSuffix = ".sorted"
        OutputSuffix = "filtered"
        samtools_bin = "samtools"
    threads: 1
    shell:
        """
        {params.samtools_bin} flagstats {input.BamDir}/{params.SeqId}{params.InputSuffix}.bam > {output.flag_stats_txt}
        """