// [METADATA]
// TOOL_NAME = cbnipt_cnv_call_cli
// THREADS = 4

process CBNIPT_CNV_CALL_CLI {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(BamPath), val(ReferenceFasta), val(OutDir)

    output:
    path "*", emit: OutDir

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def BamPath = params.BamPath ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def PythonBin = params.PythonBin ?: "/storage/home/jhkim/Apps/Python-3.11.13/python"
    def CliScript = params.CliScript ?: "/storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/scripts/cli.py"
    def Threads = params.Threads ?: "4"
    """
    ${PythonBin} ${CliScript} --SeqID ${SeqID} --BamPath ${BamPath} --ReferenceFasta ${ReferenceFasta} --OutDir ${OutDir} --Threads ${Threads}
    """
}