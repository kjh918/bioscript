// [METADATA]
// TOOL_NAME = Kraken2
// THREADS = 8

process KRAKEN2 {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(TrimFastqDir), val(KrakenDB), val(OutDir)

    output:
    path "*", emit: sub_r1
    path "*", emit: sub_r2
    path "*", emit: kraken_out
    path "*", emit: report_out
    path "*", emit: plot_out

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def TrimFastqDir = params.TrimFastqDir ?: ""
    def KrakenDB = params.KrakenDB ?: "/storage/home/jhkim/Apps/kraken2/bin/kraken2_db"
    def OutDir = params.OutDir ?: ""
    def seqkit_bin = params.seqkit_bin ?: "/storage/apps/seqkit-2.8.2/seqkit"
    def SampleSeed = params.SampleSeed ?: "100"
    def SampleNum = params.SampleNum ?: "100000"
    def R1_suffix = params.R1_suffix ?: ".trimmed_R1.fastq.gz"
    def R2_suffix = params.R2_suffix ?: ".trimmed_R2.fastq.gz"
    def kraken2_bin = params.kraken2_bin ?: "/storage/home/jhkim/Apps/kraken2/bin/kraken2"
    def krona_bin = params.krona_bin ?: "/storage/home/jhkim/Apps/Krona/KronaTools/bin/ktImportTaxonomy"
    def Threads = params.Threads ?: "8"
    """
    ${seqkit_bin} sample -s ${SampleSeed} -n ${SampleNum} ${TrimFastqDir}/${SeqID}${R1_suffix} | gzip > ${sub_r1} && ${seqkit_bin} sample -s ${SampleSeed} -n ${SampleNum} ${TrimFastqDir}/${SeqID}${R2_suffix} | gzip > ${sub_r2} && ${kraken2_bin} --db ${KrakenDB} --threads ${Threads} --report ${report_out} --paired ${sub_r1} ${sub_r2} > ${kraken_out} && ${krona_bin} -m 4 -t 3 ${kraken_out} -o ${plot_out}
    """
}