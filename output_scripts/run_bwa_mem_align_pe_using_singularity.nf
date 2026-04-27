// [METADATA]
// TOOL_NAME = bwa_mem
// THREADS = 1

process BWA_MEM {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(TrimFastqDir), path(BamDir), val(ReferenceFasta), val(ReadGroupID), val(ReadGroupPlatform), val(ReadGroupLibrary), val(ReadGroupCenter)

    output:
    path "*", emit: primary_bam

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def TrimFastqDir = params.TrimFastqDir ?: ""
    def BamDir = params.BamDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def ReadGroupID = params.ReadGroupID ?: ""
    def ReadGroupPlatform = params.ReadGroupPlatform ?: ""
    def ReadGroupLibrary = params.ReadGroupLibrary ?: ""
    def ReadGroupCenter = params.ReadGroupCenter ?: ""
    def InputSuffix = params.InputSuffix ?: "trimmed"
    def OutputSuffix = params.OutputSuffix ?: "primary"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def bwa_bin = params.bwa_bin ?: "bwa"
    def samtools_bin = params.samtools_bin ?: "samtools"
    def sif = params.sif ?: "/storage/images/bwa-0.7.17.sif"
    def bind = params.bind ?: "/storage,/data"
    def Threads = params.Threads ?: "8"
    def mark_short_split = params.mark_short_split ?: "-M"
    def soft_clipping = params.soft_clipping ?: "-Y"
    def clipping_penalty = params.clipping_penalty ?: "-L 50,50"
    def other_args = params.other_args ?: ""
    """
    ${singularity_bin} exec -B ${bind} ${sif} ${bwa_bin} mem  ${mark_short_split} ${soft_clipping} ${clipping_penalty} ${other_args}  -t ${Threads} -R '@RG\tID:${ReadGroupID}\tPL:${ReadGroupPlatform}\tLB:${ReadGroupLibrary}\tSM:${SeqID}\tCN:${ReadGroupCenter}' ${ReferenceFasta} ${TrimFastqDir}/${SeqID}.${InputSuffix}_R1.fastq.gz ${TrimFastqDir}/${SeqID}.${InputSuffix}_R2.fastq.gz | ${samtools_bin} view -@ ${Threads} -bS -o ${primary_bam} -
    """
}