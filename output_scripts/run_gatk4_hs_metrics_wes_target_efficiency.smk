# [METADATA]
# TOOL_NAME = gatk4_hs_metrics
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4_hs_metrics:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
        ReferenceFasta = ""
        BaitIntervals = ""
        TargetIntervals = ""
        InputSuffix = "analysisReady"
    output:
        hs_metrics = "[qcResDir]/[SeqID].[InputSuffix].HS.metrics.txt"
    params:
        singularity_bin = "singularity"
        gatk4_sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        xmx_mb = "16384"
        Threads = "14"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.gatk4_sif} gatk CollectHsMetrics --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' --INPUT {input.BamDir}/{input.SeqID}.{input.InputSuffix}.bam --OUTPUT {output.hs_metrics} --REFERENCE_SEQUENCE {input.ReferenceFasta} --BAIT_INTERVALS {input.BaitIntervals} --TARGET_INTERVALS {input.TargetIntervals}
        """