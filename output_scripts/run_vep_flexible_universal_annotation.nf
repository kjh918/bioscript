// [METADATA]
// TOOL_NAME = vep_flexible
// THREADS = 1

process VEP_FLEXIBLE {
    tag "$sample_id"

    // YAML의 [qcResDir] 등을 Nextflow 변수 체계로 매핑
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(SeqID), path(vcfDir), path(VcfTag), val(vepCacheDir), val(vepPluginDir), val(genomeFasta), val(clinvarData), val(cosmicData), val(alphaMissense), val(caddSnp), val(caddIndel)

    output:
    path "*", emit: vep_vcf
    path "*", emit: vep_stats

    script:
    // 로컬 변수 정의 (YAML params 기반)
    def SeqID = params.SeqID ?: ""
    def vcfDir = params.vcfDir ?: ""
    def VcfTag = params.VcfTag ?: ""
    def vepCacheDir = params.vepCacheDir ?: ""
    def vepPluginDir = params.vepPluginDir ?: ""
    def genomeFasta = params.genomeFasta ?: ""
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
    ${singularity_bin} exec -B ${bind} ${vep_sif} vep --force_overwrite --offline --cache --dir_cache ${vepCacheDir} --dir_plugins ${vepPluginDir} --fasta ${genomeFasta} --species ${species} --assembly ${assembly} --input_file ${vcfDir}/${SeqID}.${VcfTag}.vcf --output_file ${vep_vcf} --vcf --stats_text --refseq --show_ref_allele --uploaded_allele --use_transcript_ref --variant_class --sift b --polyphen b --gene_phenotype --numbers --hgvs --hgvsg --symbol --canonical --biotype --regulatory --mirna --check_existing --max_af --af_1kg --af_gnomade --exclude_predicted --pick --flag_pick --fork ${Threads} --buffer_size ${buffer_size} ${clinvar_args} ${cosmic_args} --plugin AlphaMissense,file=${alphaMissense} --plugin CADD,snv=${caddSnp},indels=${caddIndel}
    """
}