// [METADATA]
// TOOL_NAME = gatk4
// THREADS = 1

process GATK4 {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir), val(ReferenceFasta), val(KnownSnp), val(KnownIndel1), val(KnownIndel2)

    output:
    path "*", emit: recal_table
    path "*", emit: recal_bam
    path "*", emit: recal_bai

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def KnownSnp = params.KnownSnp ?: ""
    def KnownIndel1 = params.KnownIndel1 ?: ""
    def KnownIndel2 = params.KnownIndel2 ?: ""
    def InputSuffix = params.InputSuffix ?: "dedup"
    def OutputSuffix = params.OutputSuffix ?: "recal"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def gatk_bin = params.gatk_bin ?: "gatk"
    def sif = params.sif ?: "/storage/images/gatk-4.4.0.0.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "14"
    def xmx_mb = params.xmx_mb ?: "16384"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} BaseRecalibrator --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --input ${BamDir}/${SeqID}.${InputSuffix}.bam --reference ${ReferenceFasta} --output ${recal_table} --known-sites ${KnownSnp} --known-sites ${KnownIndel1} --known-sites ${KnownIndel2} && ${singularity_bin} exec -B ${bind} ${sif} ${gatk_bin} ApplyBQSR --java-options '-XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m' --input ${BamDir}/${SeqID}.${InputSuffix}.bam --bqsr-recal-file ${recal_table} --output ${recal_bam} --create-output-bam-index true
    """
}