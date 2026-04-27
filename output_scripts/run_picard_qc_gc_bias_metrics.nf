// [METADATA]
// TOOL_NAME = picard
// THREADS = 1

process PICARD {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir)

    output:
    path "*", emit: out_bam
    path "*", emit: metrics

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def java_bin = params.java_bin ?: "java"
    def picard_jar = params.picard_jar ?: "/storage/apps/bin/picard.jar"
    def Threads = params.Threads ?: "14"
    def Memory = params.Memory ?: "16384m"
    def TmpDir = params.TmpDir ?: "/tmp"
    def create_index = params.create_index ?: "true"
    def remove_duplicates = params.remove_duplicates ?: "true"
    """
    ${java_bin} -XX:ParallelGCThreads=${Threads} -Xmx${Memory} -jar ${picard_jar} MarkDuplicates  INPUT=${BamDir}/${SeqID}.sorted.bam  OUTPUT=${out_bam}  METRICS_FILE=${metrics}  CREATE_INDEX=${create_index}  REMOVE_DUPLICATES=${remove_duplicates}  TMP_DIR=${TmpDir}
    """
}