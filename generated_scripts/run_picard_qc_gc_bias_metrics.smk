# [METADATA]
# TOOL_NAME = picard
# VERSION = 3.1.0
# THREADS = 1

rule picard:
    input:
        SeqID = ""
        BamDir = ""
        ReferenceFasta = ""
    output:
        gc_bias_metrics_txt = "[BamDir]/[SeqID].[InputSuffix].gc_bias_metrics.txt"
        gc_bias_summary_txt = "[BamDir]/[SeqID].[InputSuffix].gc_bias_summary.txt"
        gc_bias_chart_pdf = "[BamDir]/[SeqID].[InputSuffix].gc_bias_metrics.pdf"
    params:
        InputSuffix = "analysisReady"
        java_bin = "java"
        picard_jar = "/storage/apps/bin/picard-3.1.0.jar"
        gc_threads = "14"
        xmx_mb = "16384"
        min_gc = "0"
        max_gc = "100"
        window_size = "100"
    threads: 1
    shell:
        """
        {params.java_bin} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx_mb}m -jar {params.picard_jar} CollectGcBiasMetrics I={input.BamDir}/{input.SeqID}.{params.InputSuffix}.bam O={output.gc_bias_metrics_txt} S={output.gc_bias_summary_txt} CHART={output.gc_bias_chart_pdf} R={input.ReferenceFasta} MINIMUM_GC={params.min_gc} MAXIMUM_GC={params.max_gc} WINDOW_SIZE={params.window_size}
        """