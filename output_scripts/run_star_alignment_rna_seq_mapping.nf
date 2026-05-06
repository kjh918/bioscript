// [METADATA]
// TOOL_NAME = star_alignment
// THREADS = 1

process STAR_ALIGNMENT {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(FastqDir), path(BamDir), val(StarIndex)

    output:
    path "*", emit: aligned_bam
    path "*", emit: aligned_transcriptome_bam
    path "*", emit: gene_counts
    path "*", emit: mapping_log

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def FastqDir = params.FastqDir ?: ""
    def BamDir = params.BamDir ?: ""
    def StarIndex = params.StarIndex ?: ""
    def Threads = params.Threads ?: "12"
    def InputSuffix = params.InputSuffix ?: ".trimmed"
    def outSAMtype = params.outSAMtype ?: "BAM SortedByCoordinate"
    def quantMode = params.quantMode ?: "TranscriptomeSAM"
    def outSAMattributes = params.outSAMattributes ?: "NH HI AS nM XS NM"
    def chimSegmentMin = params.chimSegmentMin ?: "10"
    def twopassMode = params.twopassMode ?: "Basic"
    def outFilterMismatchNmax = params.outFilterMismatchNmax ?: "10"
    def outSAMunmapped = params.outSAMunmapped ?: "Within"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def star_sif = params.star_sif ?: "/storage/images/star-2.7.11.sif"
    def bind = params.bind ?: "/storage,/data"
    """
    ${singularity_bin} exec -B ${bind} ${star_sif} STAR --runThreadN ${Threads} --genomeDir ${StarIndex} --readFilesIn ${FastqDir}/${SeqID}${InputSuffix}_R1.fastq.gz ${FastqDir}/${SeqID}${InputSuffix}_R2.fastq.gz --readFilesCommand zcat --quantMode ${quantMode} --twopassMode ${twopassMode} --chimSegmentMin ${chimSegmentMin} --outFilterMismatchNmax ${outFilterMismatchNmax} --outFileNamePrefix ${BamDir}/${SeqID}. --outSAMtype ${outSAMtype} --outSAMattributes ${outSAMattributes} --outSAMunmapped ${outSAMunmapped}
    """
}