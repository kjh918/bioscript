// [METADATA]
// TOOL_NAME = mosdepth
// THREADS = 1

process MOSDEPTH {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir), val(TargetBed)

    output:
    path "*", emit: global_dist
    path "*", emit: region_dist
    path "*", emit: summary_txt
    path "*", emit: regions_bed

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def TargetBed = params.TargetBed ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def mosdepth_sif = params.mosdepth_sif ?: "/storage/images/mosdepth-0.3.6.sif"
    def mosdepth_bin = params.mosdepth_bin ?: "/opt/mosdepth"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "4"
    def MapQ = params.MapQ ?: "20"
    def extra_args = params.extra_args ?: ""
    """
    ${singularity_bin} exec -B ${bind} ${mosdepth_sif} ${mosdepth_bin} --threads ${Threads} --by ${TargetBed} --mapq ${MapQ} ${extra_args} ${qcResDir}/${SeqID} ${BamDir}/${SeqID}.${InputSuffix}.bam
    """
}