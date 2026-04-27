// [METADATA]
// TOOL_NAME = mosdepth
// THREADS = 1

process MOSDEPTH {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir)

    output:
    path "*", emit: prefix
    path "*", emit: regions_bed_gz
    path "*", emit: regions_bed_gz_csi
    path "*", emit: summary_txt
    path "*", emit: global_dist_txt

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def mosdepth_bin = params.mosdepth_bin ?: "/opt/mosdepth"
    def sif = params.sif ?: "/storage/images/mosdepth-0.3.6.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "8"
    def bin = params.bin ?: "100000"
    def mapq = params.mapq ?: "20"
    def mosdepth_args = params.mosdepth_args ?: "--no-per-base --fast-mode"
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${mosdepth_bin} --threads ${Threads} ${mosdepth_args} --by ${bin} --mapq ${mapq} ${prefix} ${BamDir}/${SeqID}.${InputSuffix}.bam
    """
}