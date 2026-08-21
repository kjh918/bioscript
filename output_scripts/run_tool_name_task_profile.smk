# [METADATA]
# TOOL_NAME = TOOL_NAME
# VERSION = VERSION
# THREADS = 8

rule tool_name:
    input:
        SeqID = ""
        InputDir = ""
        OutDir = ""
    output:
        output_file = "[OutDir]/[SeqID].out"
    params:
        tool_bin = "TOOL_COMMAND"
        InputSuffix = ".UseAnalysis.bam"
        OutputSuffix = ".txt"
        Threads = "8"
    threads: 8
    shell:
        """
        {params.tool_bin} INPUT_OPTIONS OUTPUT_OPTIONS
        """