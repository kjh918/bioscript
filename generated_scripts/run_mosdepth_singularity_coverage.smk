# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 1

rule mosdepth:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
    output:
        prefix = "[qcResDir]/[SeqID].[InputSuffix]"
        regions_bed_gz = "[prefix].regions.bed.gz"
        regions_bed_gz_csi = "[prefix].regions.bed.gz.csi"
        summary_txt = "[prefix].mosdepth.summary.txt"
        global_dist_txt = "[prefix].global.dist.txt"
    params:
        InputSuffix = "analysisReady"
        singularity_bin = "singularity"
        mosdepth_bin = "/opt/mosdepth"
        sif = "/storage/images/mosdepth-0.3.6.sif"
        bind = "/storage,/data"
        Threads = "8"
        bin = "100000"
        mapq = "20"
        mosdepth_args = "--no-per-base --fast-mode"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.mosdepth_bin} --threads {threads} {params.mosdepth_args} --by {params.bin} --mapq {params.mapq} {output.prefix} {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam
        """