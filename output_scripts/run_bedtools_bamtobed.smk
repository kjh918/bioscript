# [METADATA]
# TOOL_NAME = bedtools
# VERSION = 2.27.1
# THREADS = 1

rule bedtools:
    input:
        SeqID = ""
        BamDir = ""
        BamToBedDir = ""
    output:
        bed_gz = "[BamDir]/[SeqID].[InputSuffix].bed.gz"
        final_bed = "[BamToBedDir]/[SeqID].bed.gz"
    params:
        InputSuffix = "analysisReady"
        bedtools_bin = "bedtools"
        bgzip_bin = "bgzip"
        ln_bin = "ln"
        cp_bin = "cp"
        Threads = "8"
    threads: 1
    shell:
        """
        {params.bedtools_bin} bamtobed -i {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam > {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bed && {params.bgzip_bin} -f -@ {threads} {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bed && {params.ln_bin} -Tsf {output.bed_gz} {input.BamDir}/{input.SeqID}.bed.gz && {params.cp_bin} {output.bed_gz} {output.final_bed}
        """