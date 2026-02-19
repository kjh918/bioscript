// [METADATA]
// TOOL_NAME = bwa_picard
// THREADS = 1

process BWA_PICARD {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(TrimFastqDir), path(BamDir), val(TmpDir), val(ReferenceFasta), val(ReadGroupID), val(ReadGroupPlatform), val(ReadGroupLibrary), val(ReadGroupCenter)

    output:
    path "*", emit: unmapped_bam
    path "*", emit: aligned_sam
    path "*", emit: primary_bam
    path "*", emit: primary_bai

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def TrimFastqDir = params.TrimFastqDir ?: ""
    def BamDir = params.BamDir ?: ""
    def TmpDir = params.TmpDir ?: ""
    def ReferenceFasta = params.ReferenceFasta ?: ""
    def ReadGroupID = params.ReadGroupID ?: ""
    def ReadGroupPlatform = params.ReadGroupPlatform ?: ""
    def ReadGroupLibrary = params.ReadGroupLibrary ?: ""
    def ReadGroupCenter = params.ReadGroupCenter ?: ""
    def InputSuffix = params.InputSuffix ?: "trimmed"
    def OutputSuffix = params.OutputSuffix ?: "primary"
    def java_bin = params.java_bin ?: "java"
    def picard_jar = params.picard_jar ?: "/storage/apps/bin/picard.jar"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def bwa_sif = params.bwa_sif ?: "/storage/images/bwa-0.7.17.sif"
    def bind = params.bind ?: "/storage,/data"
    def xmx_mb = params.xmx_mb ?: "16384"
    def Threads = params.Threads ?: "8"
    def bwa_bin = params.bwa_bin ?: "bwa"
    def mark_short_split = params.mark_short_split ?: "-M"
    def soft_clipping = params.soft_clipping ?: "-Y"
    def clipping_penalty = params.clipping_penalty ?: "-L 50,50"
    def other_args = params.other_args ?: ""
    def mba_strategy = params.mba_strategy ?: "MostDistant"
    def mba_attributes = params.mba_attributes ?: "XS"
    def mba_orientations = params.mba_orientations ?: "FR --EXPECTED_ORIENTATIONS RF"
    def mba_other_flags = params.mba_other_flags ?: "--CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false"
    """
    ${java_bin} -XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m -jar ${picard_jar} FastqToSam --FASTQ ${TrimFastqDir}/${SeqID}.${InputSuffix}_R1.fastq.gz --FASTQ2 ${TrimFastqDir}/${SeqID}.${InputSuffix}_R2.fastq.gz --SAMPLE_NAME ${SeqID} --OUTPUT ${unmapped_bam} --READ_GROUP_NAME ${ReadGroupID} --PLATFORM ${ReadGroupPlatform} --LIBRARY_NAME ${ReadGroupLibrary} --SEQUENCING_CENTER ${ReadGroupCenter} --TMP_DIR ${TmpDir} && ${singularity_bin} exec -B ${bind} ${bwa_sif} ${bwa_bin} mem  ${mark_short_split} ${soft_clipping} ${clipping_penalty} ${other_args}  -t ${Threads} -R '@RG\tID:${ReadGroupID}\tPL:${ReadGroupPlatform}\tLB:${ReadGroupLibrary}\tSM:${SeqID}\tCN:${ReadGroupCenter}' ${ReferenceFasta} ${TrimFastqDir}/${SeqID}.${InputSuffix}_R1.fastq.gz ${TrimFastqDir}/${SeqID}.${InputSuffix}_R2.fastq.gz > ${aligned_sam} && ${java_bin} -XX:ParallelGCThreads=${Threads} -Xmx${xmx_mb}m -jar ${picard_jar} MergeBamAlignment --UNMAPPED_BAM ${unmapped_bam} --ALIGNED_BAM ${aligned_sam} --REFERENCE_SEQUENCE ${ReferenceFasta} --OUTPUT ${primary_bam} --PRIMARY_ALIGNMENT_STRATEGY ${mba_strategy} --ATTRIBUTES_TO_RETAIN ${mba_attributes} --EXPECTED_ORIENTATIONS ${mba_orientations} ${mba_other_flags}
    """
}