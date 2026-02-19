// [METADATA]
// TOOL_NAME = gatk4
// THREADS = 1

process GATK4 {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir), val(ReferenceFasta)

    output:
    path "*", emit: artifacts_txt

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk_bin = params.gatk_bin ?: "gatk"
    def sif = params.sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "14"
    def xmx_mb = params.xmx_mb ?: "16384"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} CollectSequencingArtifactMetrics --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --INPUT ${BamDir}/${SeqID}.${InputSuffix}.bam --OUTPUT ${artifacts_txt} --FILE_EXTENSION .txt --REFERENCE_SEQUENCE ${ReferenceFasta}
    """
}