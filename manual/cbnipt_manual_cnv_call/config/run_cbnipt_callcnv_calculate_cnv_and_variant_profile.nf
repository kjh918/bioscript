// [METADATA]
// TOOL_NAME = cbNIPT_CallCNV
// THREADS = 4

process CBNIPT_CALLCNV {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamPath), val(AnnotatedBins)

    output:
    path "*", emit: OutDir

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamPath = params.BamPath ?: ""
    def AnnotatedBins = params.AnnotatedBins ?: ""
    def MinMapQ = params.MinMapQ ?: "20"
    def LowessFrac = params.LowessFrac ?: "0.2"
    def BaselinePloidy = params.BaselinePloidy ?: "2"
    def SmoothWindow = params.SmoothWindow ?: "5"
    def MinDepth = params.MinDepth ?: "2.0"
    def MinCoverage = params.MinCoverage ?: "0.5"
    def MaskLowerQ = params.MaskLowerQ ?: "0.01"
    def MaskUpperQ = params.MaskUpperQ ?: "0.99"
    def SegPenalty = params.SegPenalty ?: "10.0"
    def Threads = params.Threads ?: "4"
    def python_bin = params.python_bin ?: "python"
    def script_path = params.script_path ?: ""
    """
    ${python_bin} ${script_path} call-cnv --SeqID ${SeqID} --BamPath ${BamPath} --AnnotatedBins ${AnnotatedBins} --OutDir ${OutDir} --MinMapQ ${MinMapQ} --LowessFrac ${LowessFrac} --BaselinePloidy ${BaselinePloidy} --SmoothWindow ${SmoothWindow} --MinDepth ${MinDepth} --MinCoverage ${MinCoverage} --MaskLowerQ ${MaskLowerQ} --MaskUpperQ ${MaskUpperQ} --SegPenalty ${SegPenalty} --Threads ${Threads}
    """
}