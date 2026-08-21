// [METADATA]
// TOOL_NAME = Cutadapt
// THREADS = 4

process CUTADAPT {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(RawFastqDir), val(OutDir)

    output:
    path "*", emit: r1_clean
    path "*", emit: r2_clean
    path "*", emit: trim_log

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def RawFastqDir = params.RawFastqDir ?: ""
    def OutDir = params.OutDir ?: ""
    def InputSuffix = params.InputSuffix ?: ".fastq.gz"
    def AdapterR1 = params.AdapterR1 ?: "CTGTCTCTTATACACATCT"
    def AdapterR2 = params.AdapterR2 ?: "CTGTCTCTTATACACATCT"
    def QualityCutoff = params.QualityCutoff ?: "20"
    def MinLength = params.MinLength ?: "30"
    def Threads = params.Threads ?: "4"
    def cutadapt_bin = params.cutadapt_bin ?: "/storage/home/jhkim/Apps/Python-3.11.13/python -m cutadapt"
    """
    ${cutadapt_bin} -j ${Threads} -a ${AdapterR1} -A ${AdapterR2} -q ${QualityCutoff} -m ${MinLength} -o ${r1_clean} -p ${r2_clean} ${RawFastqDir}/${SeqID}_R1${InputSuffix} ${RawFastqDir}/${SeqID}_R2${InputSuffix} > ${trim_log}
    """
}