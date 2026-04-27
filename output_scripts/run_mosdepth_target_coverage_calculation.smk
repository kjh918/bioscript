# [METADATA]
# TOOL_NAME = mosdepth
# VERSION = 0.3.6
# THREADS = 1

rule mosdepth:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
        TargetBed = ""
    output:
        global_dist = "[qcResDir]/[SeqID].mosdepth.global.dist.txt"
        region_dist = "[qcResDir]/[SeqID].mosdepth.region.dist.txt"
        summary_txt = "[qcResDir]/[SeqID].mosdepth.summary.txt"
        regions_bed = "[qcResDir]/[SeqID].regions.bed.gz"
    params:
        InputSuffix = "analysisReady"
        singularity_bin = "singularity"
        mosdepth_sif = "/storage/images/mosdepth-0.3.6.sif"
        mosdepth_bin = "/opt/mosdepth"
        bind = "/storage,/data"
        Threads = "4"
        MapQ = "20"
        extra_args = ""
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.mosdepth_sif} {params.mosdepth_bin} --threads {threads} --by {input.TargetBed} --mapq {params.MapQ} {params.extra_args} {input.qcResDir}/{input.SeqID} {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam
        """