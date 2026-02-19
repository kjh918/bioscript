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
        insert_size_metrics_txt = "[qcResDir]/[SeqID].insert_size.metrics.txt"
        insert_size_hist_pdf = "[qcResDir]/[SeqID].insert_size.histogram.pdf"
    params:
        InputSuffix = "analysisReady"
        singularity_bin = "singularity"
        java_bin = "java"
        gatk_jar = "/gatk/gatk-package-4.4.0.0-local.jar"
        sif = "/storage/images/gatk-4.4.0.0.sif"
        bind = "/storage,/data"
        Threads = "14"
        xmx_mb = "16384"
    threads: 1
    shell:
        """
        {params.singularity_bin} exec -B {params.bind} {params.sif} {params.java_bin} -XX:ParallelGCThreads={threads} -Xmx{params.xmx_mb}m -jar {params.gatk_jar} CollectInsertSizeMetrics --INPUT {input.BamDir}/{input.SeqID}.bam --OUTPUT {output.insert_size_metrics_txt} --Histogram_FILE {output.insert_size_hist_pdf} --REFERENCE_SEQUENCE {input.ReferenceFasta}
        """