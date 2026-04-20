
#!/bin/bash

# script version v0.5

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: VariantAnnotation.sh [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID (Used as sample's processing product folder name)
        [ -h | --help ]           : Print this usage 
        [ --baseDir ]             : Base Work Dir
        [ --nType ]               : type of normal sample. if organoid, 'ORG'. default = 'TS' 
        [ --threads ]             : N-CPU threads (fork  size)
        [ --assembly ]            : Genome assembly
        [ --variantType ]         : variant tpye 'somatic' or 'germline' (default = somatic)
        [ --variantCallMode ]     : variant call run mode 'tonly' or 'nt' (ignored when variantType = germline)
        [ --normalID ]            : normal sample seqID ( use only when variantCallMode = nt )
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    THREADS=15
    BASE_DIR=/data/wes
    GENOME_ASSEMBLY="hg19"
    NORMAL_TYPE="TS"
    VARIANT_TYPE="somatic"
    VARIANT_CALL_MODE="tonly"
    NORMAL_ID=null
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:l:h: --long seqFolder:,seqID:,baseDir:,threads:,nType:,assembly:,variantType:,variantCallMode:,normalID:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        -s | --seqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; shift ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --nType )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = TS" ; shift ;;
                *)  NORMAL_TYPE=$2 ; shift 2 ;;
            esac ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
            esac ;;
        --variantType )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = somatic" ; shift ;;
                *)  VARIANT_TYPE=$2 ; shift 2 ;;
            esac ;;
        --variantCallMode )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = tonly" ; shift ;;
                *)  VARIANT_CALL_MODE=$2 ; shift 2 ;;
            esac ;;
        --normalID )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = null" ; shift ;;
                *)  NORMAL_ID=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    VCF_DIR=${DATA_DIR}/${SEQ_ID}/vcf
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    if [ ! -d ${LOG_DIR} ]; then touch ${LOG_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_VEP="${RUN_SINGULARITY}/vep-110.1.sif /opt/vep/src/ensembl-vep/vep"
    RUN_VCF2MAF="perl /storage/apps/vcf2maf-1.6.21/vcf2maf.pl"
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.vep.variant.annotation.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
        VEP_CACHE_DIR=/storage/references_and_index/hg38/vep
        VEP_PLUGIN_DIR=/storage/references_and_index/hg19/vep/Plugins/110
        VEP_BAM=${VEP_CACHE_DIR}/bam/GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bam
        TOPLEVEL_FASTA=${VEP_CACHE_DIR}/homo_sapiens/110_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
        CLINVAR_DATA=/storage/references_and_index/hg38/annotation/clinvar/clinvar.vcf.gz
        ALPHA_MISSENSE=/storage/references_and_index/hg38/annotation/AlphaMissense/AlphaMissense_hg38.tsv.gz
        CADD_SNP=/storage/references_and_index/hg38/annotation/CADD/gnomad.genomes.r3.0.snv.tsv.gz
        CADD_INDEL=/storage/references_and_index/hg38/annotation/CADD/gnomad.genomes.r3.0.indel.tsv.gz
    else
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
        VEP_CACHE_DIR=/storage/references_and_index/hg19/vep
        VEP_PLUGIN_DIR=${VEP_CACHE_DIR}/Plugins/110
        VEP_BAM=${VEP_CACHE_DIR}/bam/GCF_000001405.25_GRCh37.p13_knownrefseq_alns.bam
        TOPLEVEL_FASTA=${VEP_CACHE_DIR}/homo_sapiens/110_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
        CLINVAR_DATA=/storage/references_and_index/hg19/annotation/clinvar/clinvar.vcf.gz
        COSMIC_DATA=/storage/references_and_index/hg19/annotation/cosmic/Cosmic_GenomeScreensMutant_v99_GRCh37.vcf.gz
        ALPHA_MISSENSE=/storage/references_and_index/hg19/annotation/AlphaMissense/AlphaMissense_b37.tsv.gz
        CADD_SNP=/storage/references_and_index/hg19/annotation/CADD/gnomad.genomes.r2.1.1.snv.tsv.gz
        CADD_INDEL=/storage/references_and_index/hg19/annotation/CADD/gnomad.genomes.r2.1.1.indel.tsv.gz
    fi
#--------------------------------------------------------------------------------------------------#

#---| VCF TAG |------------------------------------------------------------------------------------#
    if [ "${VARIANT_TYPE}" = "germline" ]; then
        VCF_TAG="deepvariant.germline.variant.filtered"
    else
        if [ "${VARIANT_CALL_MODE}" = "tonly" ]; then
            VCF_TAG="mutect2.bias.filtered"
        else         
            if [ "${NORMAL_TYPE}" = "ORG" ]; then
                VCF_TAG="mutect2.ORG.NT.bias.filtered"
            else
                VCF_TAG="mutect2.NT.bias.filtered"
            fi
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES VariantAnnotation log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " GENOME ASSEMBLY      : hg19" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT ANNOTATION START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VEP RefSeq  ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| VEP RefSeq |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" == "hg38" ]; then
        ${RUN_VEP} --force_overwrite \
            --input_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vcf \
            --vcf \
            --output_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.refseq.vcf \
            --stats_text \
            --show_ref_allele \
            --uploaded_allele \
            --species homo_sapiens \
            --assembly GRCh38 \
            --cache \
            --dir_cache ${VEP_CACHE_DIR} \
            --dir_plugins ${VEP_PLUGIN_DIR} \
            --refseq \
            --use_transcript_ref \
            --variant_class \
            --sift b \
            --polyphen b \
            --gene_phenotype \
            --numbers \
            --hgvs \
            --hgvsg \
            --symbol \
            --canonical \
            --biotype \
            --regulatory \
            --mirna \
            --check_existing \
            --max_af \
            --af_1kg \
            --af_gnomade \
            --exclude_predicted \
            --pick \
            --flag_pick \
            --custom file=${CLINVAR_DATA},short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN \
            --plugin AlphaMissense,file=${ALPHA_MISSENSE} \
            --plugin CADD,snv=${CADD_SNP},indels=${CADD_INDEL} \
            --fork ${THREADS} \
            --buffer_size 50000 \
            --bam ${VEP_BAM} \
            --fasta ${GENOME_FASTA} \
            --offline
    else
        ${RUN_VEP} --force_overwrite \
            --input_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vcf \
            --vcf \
            --output_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.refseq.vcf \
            --stats_text \
            --show_ref_allele \
            --uploaded_allele \
            --species homo_sapiens \
            --assembly GRCh37 \
            --cache \
            --dir_cache ${VEP_CACHE_DIR} \
            --dir_plugins ${VEP_PLUGIN_DIR} \
            --refseq \
            --use_transcript_ref \
            --variant_class \
            --sift b \
            --polyphen b \
            --gene_phenotype \
            --numbers \
            --hgvs \
            --hgvsg \
            --symbol \
            --canonical \
            --biotype \
            --regulatory \
            --mirna \
            --check_existing \
            --max_af \
            --af_1kg \
            --af_gnomade \
            --exclude_predicted \
            --pick \
            --flag_pick \
            --custom file=${CLINVAR_DATA},short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN \
            --custom file=${COSMIC_DATA},short_name=COSMIC,format=vcf,fields=HGVSC%LEGACY_ID%TIER \
            --plugin AlphaMissense,file=${ALPHA_MISSENSE} \
            --plugin CADD,snv=${CADD_SNP},indels=${CADD_INDEL} \
            --fork ${THREADS} \
            --buffer_size 50000 \
            --bam ${VEP_BAM} \
            --fasta ${GENOME_FASTA} \
            --offline
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VEP RefSeq  ] FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VEP Ensembl ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| VEP Ensembl |--------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        ${RUN_VEP} --force_overwrite \
            --input_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vcf \
            --vcf \
            --output_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.ensembl.vcf \
            --stats_text \
            --species homo_sapiens \
            --assembly GRCh38 \
            --cache \
            --dir_cache ${VEP_CACHE_DIR} \
            --ccds \
            --uniprot \
            --canonical \
            --hgvs \
            --hgvsg \
            --domains \
            --exclude_predicted \
            --fork ${THREADS} \
            --buffer_size 50000 \
            --bam ${VEP_BAM} \
            --fasta ${GENOME_FASTA} \
            --offline
    else
        ${RUN_VEP} --force_overwrite \
            --input_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vcf \
            --vcf \
            --output_file ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.ensembl.vcf \
            --stats_text \
            --species homo_sapiens \
            --assembly GRCh37 \
            --cache \
            --dir_cache ${VEP_CACHE_DIR} \
            --ccds \
            --uniprot \
            --canonical \
            --hgvs \
            --hgvsg \
            --domains \
            --exclude_predicted \
            --fork ${THREADS} \
            --buffer_size 50000 \
            --bam ${VEP_BAM} \
            --fasta ${GENOME_FASTA} \
            --offline
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VEP Ensembl ] FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VCF to MAF START." >> ${LOGFILE}
#==================================================================================================# 

#---| VCF to MAF |---------------------------------------------------------------------------------#
    if [ "${VARIANT_CALL_MODE}" = "nt" ]; then
        ## RefSeq
        ${RUN_VCF2MAF} \
            --inhibit-vep \
            --input-vcf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.refseq.vcf \
            --output-maf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.refseq.maf \
            --tumor-id ${SEQ_ID} \
            --normal-id ${NORMAL_ID} \
            --ref-fasta ${GENOME_FASTA} \
            --retain-info DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus \
            --retain-ann HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER
        ## Ensembl
        ${RUN_VCF2MAF} \
            --inhibit-vep \
            --input-vcf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.ensembl.vcf \
            --output-maf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.ensembl.maf \
            --tumor-id ${SEQ_ID} \
            --normal-id ${NORMAL_ID} \
            --ref-fasta ${GENOME_FASTA} \
            --retain-info FILTER \
            --retain-ann HGVSg
    else
        ## RefSeq
        ${RUN_VCF2MAF} \
            --inhibit-vep \
            --input-vcf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.refseq.vcf \
            --output-maf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.refseq.maf \
            --tumor-id ${SEQ_ID} \
            --ref-fasta ${GENOME_FASTA} \
            --retain-info DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus \
            --retain-ann HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER
        ## Ensembl
        ${RUN_VCF2MAF} \
            --inhibit-vep \
            --input-vcf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.ensembl.vcf \
            --output-maf ${VCF_DIR}/${SEQ_ID}.${VCF_TAG}.vep.ensembl.maf \
            --tumor-id ${SEQ_ID} \
            --ref-fasta ${GENOME_FASTA} \
            --retain-info FILTER \
            --retain-ann HGVSg
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VCF to MAF FINISHED." >> ${LOGFILE}
    echo " " >> ${LOGFILE}
    echo "WES VariantAnnotation DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

