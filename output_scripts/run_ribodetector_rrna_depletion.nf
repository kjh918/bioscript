// [METADATA]
// TOOL_NAME = ribodetector
// THREADS = 1

process RIBODETECTOR {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(FastqDir), val(OutputDir)

    output:
    path "*", emit: non_rrna_r1
    path "*", emit: non_rrna_r2

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def FastqDir = params.FastqDir ?: ""
    def OutputDir = params.OutputDir ?: ""
    def Threads = params.Threads ?: "10"
    def ReadLen = params.ReadLen ?: "100"
    def ChunkSize = params.ChunkSize ?: "256"
    def ExcludeMode = params.ExcludeMode ?: "rrna"
    def InputSuffix = params.InputSuffix ?: "fastq.gz"
    """
    ribodetector_cpu -t ${Threads} -l ${ReadLen} -i ${FastqDir}/${SeqID}_R1.${InputSuffix} ${FastqDir}/${SeqID}_R2.${InputSuffix} -e ${ExcludeMode} --chunk_size ${ChunkSize} -o ${non_rrna_r1} ${non_rrna_r2}
    """
}