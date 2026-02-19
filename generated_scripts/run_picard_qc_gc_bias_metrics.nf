// [METADATA]
// TOOL_NAME = picard
// THREADS = 1

process PICARD {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(ReferenceFasta)

    output:
    path "*", emit: gc_bias_metrics_txt
    path "*", emit: gc_bias_summary_txt
    path "*", emit: gc_bias_chart_pdf

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def java_bin = params.java_bin ?: "java"
    def picard_jar = params.picard_jar ?: "/storage/apps/bin/picard-3.1.0.jar"
    def gc_threads = params.gc_threads ?: "14"
    def xmx_mb = params.xmx_mb ?: "16384"
    def min_gc = params.min_gc ?: "0"
    def max_gc = params.max_gc ?: "100"
    def window_size = params.window_size ?: "100"
    """
    ${java_bin} -XX:ParallelGCThreads=${gc_threads} -Xmx${xmx_mb}m -jar ${picard_jar} CollectGcBiasMetrics I=${BamDir}/${SeqID}.${InputSuffix}.bam O=${gc_bias_metrics_txt} S=${gc_bias_summary_txt} CHART=${gc_bias_chart_pdf} R=${ReferenceFasta} MINIMUM_GC=${min_gc} MAXIMUM_GC=${max_gc} WINDOW_SIZE=${window_size}
    """
}