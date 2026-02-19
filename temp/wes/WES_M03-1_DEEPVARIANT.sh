
#! /bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A09
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.Deepvariant_GermlineVariantCall.sh [ ARGUMENTS ]...
        [ -d | --seqFolder ]   : data processing base directory
        [ -s | --seqID ]       : sample ID (Used as sample's processing product folder name)
        [ -h | --help ]        : Print this usage 
        [ --baseDir ]          : Base Work Dir
        [ --bed ]              : target region bed file
        [ --modelType ]        : model type (WGS or WES)
        [ --threads ]          : threads 
        [ --assembly ]            : Genome assembly version 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    THREADS=4
    BASE_DIR=/data/wes
    TARGET_BED_FILE=/storage/references_and_index/hg19/bed/Twist_WES_2.0/hg19_twist.exome.2.0.target.bed
    GENOME_ASSEMBLY=hg19
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:l:h: --long seqFolder:,seqID:,baseDir:,bed:,modelType:,threads:,assembly:,help -- "$@" )
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
        --bed )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TARGET_BED_FILE=$2 ; shift 2 ;;
            esac ;;
        --modelType )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  MODEL_TYPE=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
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
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam
    GVCF_DIR=${DATA_DIR}/${SEQ_ID}/vcf
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    ##
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
    if [ ! -d ${GVCF_DIR} ]; then mkdir -p ${GVCF_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"   
    RUN_DEEPVARIANT="${RUN_SINGULARITY}/deepvariant-1.6.0.sif /opt/deepvariant/bin/run_deepvariant"
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.deepvariant.germline.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        # GENOME FASTA #
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    else
        # GENOME FASTA #
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta        
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES GermlineVariantCalling DeepVariant log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " GENOME ASSEMBLY      : ${GENOME_ASSEMBLY}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING START." >> ${LOGFILE}
#==================================================================================================# 

#---| RUN DEEPVARIANT |----------------------------------------------------------------------------#
    ${RUN_DEEPVARIANT} \
        --model_type ${MODEL_TYPE} \
        --sample_name ${SEQ_ID} \
        --num_shards ${THREADS} \
        --postprocess_cpus ${THREADS} \
        --ref ${GENOME_FASTA} \
        --regions ${TARGET_BED_FILE} \
        --reads ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
        --output_vcf ${GVCF_DIR}/${SEQ_ID}.deepvariant.germline.variant.vcf.gz
#--------------------------------------------------------------------------------------------------#

#---| VCF DECOMPRESS & FILTERING |-----------------------------------------------------------------#
    # DECOMPRESS
    bgzip -c -d ${GVCF_DIR}/${SEQ_ID}.deepvariant.germline.variant.vcf.gz > ${GVCF_DIR}/${SEQ_ID}.deepvariant.germline.variant.vcf
    #--------------------------------------------------------------------------#
    # FILTERING
    grep ^# ${GVCF_DIR}/${SEQ_ID}.deepvariant.germline.variant.vcf > ${GVCF_DIR}/${SEQ_ID}.deepvariant.germline.variant.filtered.vcf
    #--------------------------------------------------------------------------#
    bcftools query -i "FILTER='PASS' && QUAL>20 && FMT/DP>20" -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT" ${GVCF_DIR}/${SEQ_ID}.deepvariant.germline.variant.vcf >> ${GVCF_DIR}/${SEQ_ID}.deepvariant.germline.variant.filtered.vcf
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | VARIANT-CALLING FINISHED." >> ${LOGFILE}
    echo " " >> ${LOGFILE}
    echo "WES GermlineVariantCalling DeepVariant DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

