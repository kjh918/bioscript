// [METADATA]
// TOOL_NAME = gatk4_hs_metrics
// THREADS = 1

process GATK4_HS_METRICS {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir), val(ReferenceFasta), val(BaitIntervals), val(TargetIntervals), val(InputSuffix)

    output:
    path "*", emit: hs_metrics

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def BaitIntervals = params.BaitIntervals ?: ""
    def TargetIntervals = params.TargetIntervals ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk4_sif = params.gatk4_sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def xmx_mb = params.xmx_mb ?: "16384"
    def Threads = params.Threads ?: "14"
    """
    ${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk CollectHsMetrics --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --INPUT ${BamDir}/${SeqID}.${InputSuffix}.bam --OUTPUT ${hs_metrics} --REFERENCE_SEQUENCE ${ReferenceFasta} --BAIT_INTERVALS ${BaitIntervals} --TARGET_INTERVALS ${TargetIntervals}
    """
}