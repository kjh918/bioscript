# [METADATA]
# TOOL_NAME = salmon
# VERSION = 1.10.0
# THREADS = 1

rule salmon:
    input:
        SeqID = ""
        BamDir = ""
        OutputDir = ""
        SalmonIndex = ""
    output:
        quant_sf = "[OutputDir]/[SeqID]_quant/quant.sf"
        quant_log = "[OutputDir]/[SeqID]_quant/logs/salmon_quant.log"
        lib_format = "[OutputDir]/[SeqID]_quant/lib_format_counts.json"
    params:
        Threads = "12"
        libType = "A"
        InputSuffix = ".Aligned.toTranscriptome.out.bam"
        extra_args = ""
        samlon_bin = "/storage/apps/salmon-1.10.0/bin/salmon"
    threads: 1
    shell:
        """
        {params.samlon_bin} quant -t {input.SalmonIndex} --threads {threads} --libType {params.libType} -a {input.BamDir}/{input.SeqID}{params.InputSuffix} -o {input.OutputDir} {params.extra_args}
        """