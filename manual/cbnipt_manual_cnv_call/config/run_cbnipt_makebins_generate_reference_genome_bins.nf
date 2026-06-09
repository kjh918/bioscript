// [METADATA]
// TOOL_NAME = cbNIPT_MakeBins
// THREADS = 1

process CBNIPT_MAKEBINS {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(ReferenceFasta), val(MappabilityBW)

    output:
    path "*", emit: OutBinFile

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def MappabilityBW = params.MappabilityBW ?: ""
    def BinSize = params.BinSize ?: "100000"
    def MinMappability = params.MinMappability ?: "0.9"
    def MinGC = params.MinGC ?: "0.3"
    def MaxGC = params.MaxGC ?: "0.7"
    def IncludeSexChrom = params.IncludeSexChrom ?: "--IncludeSexChrom"
    def python_bin = params.python_bin ?: "python"
    def script_path = params.script_path ?: ""
    """
    ${python_bin} ${script_path} make-bins --ReferenceFasta ${ReferenceFasta} --MappabilityBW ${MappabilityBW} --OutBinFile ${OutBinFile} --BinSize ${BinSize} --MinMappability ${MinMappability} --MinGC ${MinGC} --MaxGC ${MaxGC} ${IncludeSexChrom}
    """
}