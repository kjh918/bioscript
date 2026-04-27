# [METADATA]
# TOOL_NAME = samtools_filter_index
# VERSION = 1.10
# THREADS = 1

rule samtools_filter_index:
    input:
        SeqID = ""
        BamDir = ""
    output:
        filtered_bam = "[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam"
        filtered_bai = "[BamDir]/[SeqID].[InputSuffix].[OutputSuffix].bam.bai"
        analysis_ready_bam = "[BamDir]/[SeqID].analysisReady.bam"
        analysis_ready_bai = "[BamDir]/[SeqID].analysisReady.bam.bai"
    params:
        InputSuffix = "recal"
        OutputSuffix = "filtered"
        samtools_bin = "samtools"
        ln_bin = "ln"
        Threads = "8"
        min_mapq = "20"
        include_flag = "0x2"
        exclude_flag = "0x104"
        expr_nm = "[NM] < 12"
    threads: 1
    shell:
        """
        {params.samtools_bin} view -b -h -q {params.min_mapq} -f {params.include_flag} -F {params.exclude_flag} -e '{params.expr_nm}' --threads {threads} {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam > {output.filtered_bam} && {params.samtools_bin} index --threads {threads} {output.filtered_bam} && {params.ln_bin} -Tsf {output.filtered_bam} {output.analysis_ready_bam} && {params.ln_bin} -Tsf {output.filtered_bai} {output.analysis_ready_bai}
        """