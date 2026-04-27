// [METADATA]
// TOOL_NAME = vep_ensembl
// THREADS = 1

process VEP_ENSEMBL {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(vcfDir), path(VcfTag), val(vepCacheDir), val(genomeFasta), path(vepBam), val(clinvarData), val(cosmicData), val(alphaMissense), val(caddSnp), val(caddIndel)

    output:
    path "*", emit: vep_vcf
    path "*", emit: vep_stats

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def vcfDir = params.vcfDir ?: ""
    def VcfTag = params.VcfTag ?: ""
    def vepCacheDir = params.vepCacheDir ?: ""
    def genomeFasta = params.genomeFasta ?: ""
    def vepBam = params.vepBam ?: ""
    def clinvarData = params.clinvarData ?: ""
    def cosmicData = params.cosmicData ?: ""
    def alphaMissense = params.alphaMissense ?: ""
    def caddSnp = params.caddSnp ?: ""
    def caddIndel = params.caddIndel ?: ""
    def assembly = params.assembly ?: "GRCh38"
    def species = params.species ?: "homo_sapiens"
    def Threads = params.Threads ?: "8"
    def buffer_size = params.buffer_size ?: "50000"
    def cosmic_args = params.cosmic_args ?: ""
    def clinvar_args = params.clinvar_args ?: "--custom file=[clinvarData],short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN"
    def singularity_bin = params.singularity_bin ?: "singularity"
    def vep_sif = params.vep_sif ?: "/storage/images/vep-110.sif"
    def bind = params.bind ?: "/storage,/data"
    """
    ${singularity_bin} exec -B ${bind} ${vep_sif} vep --force_overwrite --offline --cache --dir_cache ${vepCacheDir} --fasta ${genomeFasta} --bam ${vepBam} --species ${species} --assembly ${assembly} --input_file ${vcfDir}/${SeqID}.${VcfTag}.vcf --output_file ${vep_vcf} --vcf --stats_text --hgvs --hgvsg --canonical --exclude_predicted --ccds --uniprot --domains --fork ${Threads} --buffer_size ${buffer_size} ${clinvar_args} ${cosmic_args} --plugin AlphaMissense,file=${alphaMissense} --plugin CADD,snv=${caddSnp},indels=${caddIndel}
    """
}