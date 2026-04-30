# [METADATA]
# TOOL_NAME = salmon_quant
# VERSION = 1.10.1
# THREADS = 1

rule salmon_quant:
    input:
        SeqID = ""
        FastqDir = ""
        OutputDir = ""
        SalmonIndex = ""
    output:
        quant_sf = "[OutputDir]/[SeqID]_quant/quant.sf"
        quant_log = "[OutputDir]/[SeqID]_quant/logs/salmon_quant.log"
        lib_format = "[OutputDir]/[SeqID]_quant/lib_format_counts.json"
    params:
        Threads = "12"
        libType = "A"
        InputSuffix = "fastq.gz"
        validateMappings = "--validateMappings"
        extra_args = ""
    threads: 1
    shell:
        """
        salmon quant --threads {threads} --index {input.SalmonIndex} --libType {params.libType} -1 {input.FastqDir}/{input.SeqID}_R1.{params.InputSuffix} -2 {input.FastqDir}/{input.SeqID}_R2.{params.InputSuffix} {params.validateMappings} -o {input.OutputDir}/{input.SeqID}_quant {params.extra_args}
        """