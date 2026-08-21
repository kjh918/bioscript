// [METADATA]
// TOOL_NAME = SRAToolkit
// THREADS = 8

process SRATOOLKIT {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SraAccession), val(OutDir), val(TmpDir)

    output:
    path "*", emit: sra_dir
    path "*", emit: r1_fastq
    path "*", emit: r2_fastq
    path "*", emit: r1_fastq_gz
    path "*", emit: r2_fastq_gz

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SraAccession = params.SraAccession ?: ""
    def OutDir = params.OutDir ?: ""
    def TmpDir = params.TmpDir ?: "[OutDir]/tmp"
    def prefetch_bin = params.prefetch_bin ?: "/storage/apps/sratoolkit.3.0.10/bin/prefetch"
    def fasterq_dump_bin = params.fasterq_dump_bin ?: "/storage/apps/sratoolkit.3.0.10/bin/fasterq-dump"
    def Threads = params.Threads ?: "8"
    def MaxSize = params.MaxSize ?: "100G"
    def SplitMode = params.SplitMode ?: "--split-files"
    """
    mkdir -p ${OutDir} ${TmpDir}/sra ${TmpDir}/fasterq && ${prefetch_bin} --max-size ${MaxSize} --output-directory ${TmpDir}/sra ${SraAccession} && ${fasterq_dump_bin} ${SplitMode} --threads ${Threads} --temp ${TmpDir}/fasterq --outdir ${OutDir} ${sra_dir} && mv ${OutDir}/${SraAccession}_1.fastq ${r1_fastq} && mv ${OutDir}/${SraAccession}_2.fastq ${r2_fastq} && gzip -f ${r1_fastq} && gzip -f ${r2_fastq}
    """
}