// [METADATA]
// TOOL_NAME = CopyKit
// THREADS = 1

process COPYKIT {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(NGS_DataBaseDir), val(ResultBaseDir)

    output:
    path "*", emit: AnalysisRunDir

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def NGS_DataBaseDir = params.NGS_DataBaseDir ?: ""
    def ResultBaseDir = params.ResultBaseDir ?: ""
    def InputSuffix = params.InputSuffix ?: "analysisReady"
    def mkdir_bin = params.mkdir_bin ?: "mkdir"
    def rscript_bin = params.rscript_bin ?: "Rscript"
    def Rscript_path = params.Rscript_path ?: "/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R"
    def Threads = params.Threads ?: "1"
    def BinSize = params.BinSize ?: "220kb"
    def GenomeVersion = params.GenomeVersion ?: "hg38"
    def SamplePloidy = params.SamplePloidy ?: "2"
    """
    ${mkdir_bin} -p ${AnalysisRunDir} && ln -Tsf ${NGS_DataBaseDir}/${SeqID}.${InputSuffix}.bam ${AnalysisRunDir}/${SeqID}.bam && ln -Tsf ${NGS_DataBaseDir}/${SeqID}.${InputSuffix}.bam.bai ${AnalysisRunDir}/${SeqID}.bam.bai && ${rscript_bin} ${Rscript_path}  --SeqID ${SeqID}  --AnalysisRunDir ${AnalysisRunDir}  --BinSize ${BinSize}  --Ploidy ${SamplePloidy}  --GenomeVersion ${GenomeVersion}
    """
}