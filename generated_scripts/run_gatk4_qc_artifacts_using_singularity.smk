# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1

rule gatk4:
    input:
        SeqID = ""
        BamDir = ""
        qcResDir = ""
        ReferenceFasta = ""
    output:
        artifacts_txt = "[qcResDir]/[SeqID].artifacts.txt"
    params:
        InputSuffix = "analysisReady"
        singularity_bin = "singularity"
        gatk_bin = "gatk"
        sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        Threads = "14"
        xmx_mb = "16384"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.gatk_bin} CollectSequencingArtifactMetrics --java-options '-XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m' --INPUT {input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam --OUTPUT {output.artifacts_txt} --FILE_EXTENSION .txt --REFERENCE_SEQUENCE {input.ReferenceFasta}
        """