// [METADATA]
// TOOL_NAME = bedtools_genomecov_to_bigwig
// THREADS = 1

process BEDTOOLS_GENOMECOV_TO_BIGWIG {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), val(InputFile), val(GenomeSizes), val(Extension), val(OutputDir)

    output:

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def InputFile = params.InputFile ?: ""
    def GenomeSizes = params.GenomeSizes ?: ""
    def Extension = params.Extension ?: ""
    def OutputDir = params.OutputDir ?: ""
    def Threads = params.Threads ?: "4"
    def InputFlag = params.InputFlag ?: "-ibam"
    def GenomeSizesFlag = params.GenomeSizesFlag ?: ""
    def ScaleArgs = params.ScaleArgs ?: "-bg"
    def SortCmd = params.SortCmd ?: "| LC_ALL=C sort --buffer-size=8G -k1,1 -k2,2n"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def bedtools_sif = params.bedtools_sif ?: "/storage/images/bedtools-2.27.1.sif"
    def bigwig_sif = params.bigwig_sif ?: "/storage/images/ucsc-bedgraphtobigwig-445.sif"
    def bind = params.bind ?: "/storage,/data"
    """
    ${singularity_bin} exec -B ${bind} ${bedtools_sif} bedtools genomecov ${InputFlag} ${InputFile} ${GenomeSizesFlag} ${ScaleArgs} ${SortCmd} --parallel=${Threads} > ${OutputDir}/${SeqID}.${Extension} && ${singularity_bin} exec -B ${bind} ${bigwig_sif} bedGraphToBigWig ${OutputDir}/${SeqID}.${Extension} ${GenomeSizes} ${OutputDir}/${SeqID}.bw
    """
}