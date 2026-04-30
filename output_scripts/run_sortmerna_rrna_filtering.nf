// [METADATA]
// TOOL_NAME = sortmerna
// THREADS = 1

process SORTMERNA {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(FastqDir), val(qcResDir), val(RefArgs), val(IndexDir)

    output:
    path "*", emit: non_rrna_r1
    path "*", emit: non_rrna_r2
    path "*", emit: sortmerna_log

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def FastqDir = params.FastqDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def RefArgs = params.RefArgs ?: ""
    def IndexDir = params.IndexDir ?: ""
    def Threads = params.Threads ?: "8"
    def InputSuffix = params.InputSuffix ?: "fastq.gz"
    def paired_cmd = params.paired_cmd ?: "--paired_in --out2"
    def extra_args = params.extra_args ?: ""
    """
    sortmerna ${RefArgs} --reads ${FastqDir}/${SeqID}_1.${InputSuffix} --reads ${FastqDir}/${SeqID}_2.${InputSuffix} --threads ${Threads} --workdir . --aligned ${qcResDir}/rRNA_reads --fastx --other ${qcResDir}/non_rRNA_reads ${paired_cmd} ${extra_args} && mv ${qcResDir}/non_rRNA_reads_fwd.f*q.gz ${non_rrna_r1} && mv ${qcResDir}/non_rRNA_reads_rev.f*q.gz ${non_rrna_r2} && mv ${qcResDir}/rRNA_reads.log ${sortmerna_log}
    """
}