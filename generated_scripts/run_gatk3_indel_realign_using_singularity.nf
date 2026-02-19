// [METADATA]
// TOOL_NAME = gatk3
// THREADS = 1

process GATK3 {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir), val(ReferenceFasta), val(KnownIndel1), val(KnownIndel2)

    output:
    path "*", emit: target_intervals
    path "*", emit: realigned_bam

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def KnownIndel1 = params.KnownIndel1 ?: ""
    def KnownIndel2 = params.KnownIndel2 ?: ""
    def InputSuffix = params.InputSuffix ?: "dedup"
    def OutputSuffix = params.OutputSuffix ?: "realign"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def java_bin = params.java_bin ?: "java"
    def gatk_jar = params.gatk_jar ?: "/usr/GenomeAnalysisTK.jar"
    def sif = params.sif ?: "/storage/images/gatk-3.8-1.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "8"
    def xmx_mb = params.xmx_mb ?: "16384"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${java_bin} -Xmx${xmx_mb}m -jar ${gatk_jar} -T RealignerTargetCreator -R ${ReferenceFasta} -I ${BamDir}/${SeqID}.${InputSuffix}.bam -o ${target_intervals} -known ${KnownIndel1} -known ${KnownIndel2} -nt ${Threads} && ${singularity_bin} exec -B ${bind} ${sif} ${java_bin} -Xmx${xmx_mb}m -jar ${gatk_jar} -T IndelRealigner -R ${ReferenceFasta} -targetIntervals ${target_intervals} -known ${KnownIndel1} -known ${KnownIndel2} -I ${BamDir}/${SeqID}.${InputSuffix}.bam -o ${realigned_bam}
    """
}