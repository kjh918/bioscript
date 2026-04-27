// [METADATA]
// TOOL_NAME = gatk4_learn_read_orientation_model
// THREADS = 1

process GATK4_LEARN_READ_ORIENTATION_MODEL {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(qcResDir)

    output:
    path "*", emit: read_orientation_model

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def qcResDir = params.qcResDir ?: ""
    def InputSuffix = params.InputSuffix ?: ""
    def OutputSuffix = params.OutputSuffix ?: ""
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk4_sif = params.gatk4_sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def xmx_mb = params.xmx_mb ?: "8192"
    def Threads = params.Threads ?: "14"
    """
    ${singularity_bin} exec -B ${bind} ${gatk4_sif} gatk LearnReadOrientationModel --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' -I ${qcResDir}/${SeqID}.${InputSuffix}f1r2.tar.gz -O ${read_orientation_model}
    """
}