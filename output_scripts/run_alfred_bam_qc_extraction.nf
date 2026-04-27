// [METADATA]
// TOOL_NAME = alfred
// THREADS = 1

process ALFRED {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamDir), val(qcResDir), val(ReferenceFasta), val(TargetBed), val(InputSuffix)

    output:
    path "*", emit: alfred_raw_tsv
    path "*", emit: chr_map_stats
    path "*", emit: target_coverage

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamDir = params.BamDir ?: ""
    def qcResDir = params.qcResDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def TargetBed = params.TargetBed ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def alfred_sif = params.alfred_sif ?: "/storage/images/alfred-0.2.6.sif"
    def alfred_bin = params.alfred_bin ?: "/opt/alfred/bin/alfred"
    def bind = params.bind ?: "/storage,/data"
    """
    ${singularity_bin} exec -B ${bind} ${alfred_sif} ${alfred_bin} qc --reference ${ReferenceFasta} --bed ${TargetBed} --outfile ${alfred_raw_tsv} ${BamDir}/${SeqID}.${InputSuffix}.bam && zgrep '^CM' ${alfred_raw_tsv} | cut -f 2- > ${chr_map_stats} && zgrep '^TC' ${alfred_raw_tsv} | cut -f 2- > ${target_coverage}
    """
}