// [METADATA]
// TOOL_NAME = bbmap_bbsplit
// THREADS = 1

process BBMAP_BBSPLIT {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(Read1), val(Read2), val(PrimaryRef), val(OtherRefsArgs), val(OutputDir)

    output:

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def Read1 = params.Read1 ?: ""
    def Read2 = params.Read2 ?: ""
    def PrimaryRef = params.PrimaryRef ?: ""
    def OtherRefsArgs = params.OtherRefsArgs ?: ""
    def OutputDir = params.OutputDir ?: ""
    def MemoryMB = params.MemoryMB ?: ""
    def Threads = params.Threads ?: ""
    def BBMapPath = params.BBMapPath ?: ""
    def BaseNamePattern = params.BaseNamePattern ?: ""
    def ExtraArgs = params.ExtraArgs ?: "ambiguous2=all"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def bbmap_sif = params.bbmap_sif ?: "/storage/home/jhkim/Apps/bbmap_39.81.sif"
    def bind = params.bind ?: "/storage,/data"
    """
    ${singularity_bin} exec -B ${bind} ${bbmap_sif} bbsplit.sh -Xmx${MemoryMB}M ref_primary=${PrimaryRef} ${OtherRefsArgs} in1=${Read1} in2=${Read2} basename=${BaseNamePattern} refstats=${OutputDir}/${SeqID}.refstats.txt threads=${Threads} ${ExtraArgs}
    """
}