
#! /bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.MSI_Calculation.sh [ ARGUMENTS ]...
        [ --seqFolder ] : data processing base directory
        [ --tumorSeqID ]     : sample ID (Used as sample's processing product folder name)
        [ --normalSeqID ]    : normal sample SEQ_ID
        [ --targetBed ]      : target region bed file
        [ --baseDir ]        : work base folder
        [ --threads ]        : CPU core threads
        [ --runMode ]        : Run Mode. "tonly" or "nt"
        [ --assembly ]       : Genome assembly version
        [ -h | --help ]      : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR=/data/wes
    NORMAL_SEQ_ID=""
    RUN_MODE="tonly"
    THREADS=1
    TARGET_BED=/storage/references_and_index/hg19/bed/Twist_WES_2.0/hg19_twist.exome.2.0.target.bed
    GENOME_ASSEMBLY=hg19
    NORMAL_TYPE=TS
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long seqFolder:,tumorSeqID:,normalSeqID:,targetBed:,baseDir:,threads:,assembly:,runMode:,help -- "$@" )
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
        --tumorSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  TUMOR_SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        --normalSeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = NULL" ; shift ;;
                *)  NORMAL_SEQ_ID=$2 ; shift 2 ;;
            esac ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes" ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --targetBed )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = twist exome 2.0" ; shift ;;
                *)  TARGET_BED=$2 ; shift 2 ;;
            esac ;;
        --runMode )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = tumor only" ; shift ;;
                *)  RUN_MODE=$2 ; shift 2 ;;
            esac ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 1" ; shift ;;
                *)  THREADS=$2 ; shift 2 ;;
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
    TUMOR_BAM_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/bam
    NORMAL_BAM_DIR=${DATA_DIR}/${NORMAL_SEQ_ID}/bam
    LOG_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/log
    MSI_RES_DIR=${DATA_DIR}/${TUMOR_SEQ_ID}/msi
    #--------------------------------------------------------------------------#
    if [ ! -d ${MSI_RES_DIR} ]; then mkdir -p ${MSI_RES_DIR}; fi
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
#--------------------------------------------------------------------------------------------------#
    
#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_CONDA_MSI="${RUN_SINGULARITY}/conda_msi.sif"
    RUN_MSISENSOR2="/storage/apps/msisensor2-0.1/msisensor2 msi"
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${TUMOR_SEQ_ID}.MSI.calculation.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
        #----------------------------------------------------------------------#
        MANTIS_SITE_LIST=/storage/references_and_index/hg38/msi/hg38_microsatelittes_MANTIS.list
        MSISENSOR_MODEL=/storage/references_and_index/hg38/msi/hg38_msisensor2_models/
        MSISENSOR_SITE_LIST=/storage/references_and_index/hg38/msi/hg38_microsatelittes_MSISensor2.list
    else
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
        #----------------------------------------------------------------------#
        MANTIS_SITE_LIST=/storage/references_and_index/hg19/msi/hg19_microsatelittes_MANTIS.list
        MSISENSOR_MODEL=/storage/references_and_index/hg19/msi/hg19_msisensor2_models/
        MSISENSOR_SITE_LIST=/storage/references_and_index/hg19/msi/hg19_microsatelittes_MSISensor2.list
    fi        
#--------------------------------------------------------------------------------------------------#

#---| TAGS |---------------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" == "tonly" ]; then
        NORMAL_SAMPLE=""
        SCRIPT_RUN_MODE="Tumor Only"
    else
        NORMAL_SAMPLE=${NORMAL_SEQ_ID}
        SCRIPT_RUN_MODE="Matched Normal"
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES MSI Calculation log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " TUMOR_SEQ_ID         : ${TUMOR_SEQ_ID}" >> ${LOGFILE}
    echo " NORMAL_SEQ_ID        : ${NORMAL_SAMPLE}" >> ${LOGFILE}
    echo " GENOME ASSEMBLY      : ${GENOME_ASSEMBLY}" >> ${LOGFILE}
    echo " RUN MODE             : ${SCRIPT_RUN_MODE}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MSI CALCULATION START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSISensor2 ] START." >> ${LOGFILE}
#==================================================================================================# 

#---| MSISensor2 RUN |-----------------------------------------------------------------------------#
    if [ "${RUN_MODE}" == "tonly" ]; then
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSISensor2 ] Tumor-Only Mode START." >> ${LOGFILE}
        #==========================================================================================#
        ${RUN_MSISENSOR2} \
            -M ${MSISENSOR_MODEL} \
            -t ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            -e ${TARGET_BED} \
            -b ${THREADS} \
            -o ${MSI_RES_DIR}/${TUMOR_SEQ_ID}.msisensor2.msi
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSISensor2 ] Tumor-Only Mode FINISHED." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSISensor2 ] Matched-Normal Mode START." >> ${LOGFILE}
        #==========================================================================================#
        NORMAL_BAM_DIR=${DATA_DIR}/${NORMAL_SEQ_ID}/bam
        #----------------------------------------------------------------------#
        ${RUN_MSISENSOR2} \
            -d ${MSISENSOR_SITE_LIST} \
            -t ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            -n ${NORMAL_BAM_DIR}/${NORMAL_SEQ_ID}.analysisReady.bam \
            -e ${TARGET_BED} \
            -b ${THREADS} \
            -o ${MSI_RES_DIR}/${TUMOR_SEQ_ID}.msisensor2.msi  
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSISensor2 ] Matched-Normal Mode FINISHED." >> ${LOGFILE}
        #==========================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#

#---| MANTIS RUN |---------------------------------------------------------------------------------#
    if [ "${RUN_MODE}" == "tonly" ]; then
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MANTIS ] Run SKIPPED." >> ${LOGFILE}
        #==========================================================================================#
    else
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MANTIS ] Run START." >> ${LOGFILE}
        #==========================================================================================#
        ${RUN_CONDA_MSI} python /storage/apps/mantis-1.0.5/mantis.py \
            --threads ${THREADS} \
            --bedfile ${MANTIS_SITE_LIST} \
            --genome ${GENOME_FASTA} \
            -n ${NORMAL_BAM_DIR}/${NORMAL_SEQ_ID}.analysisReady.bam \
            -t ${TUMOR_BAM_DIR}/${TUMOR_SEQ_ID}.analysisReady.bam \
            -o ${MSI_RES_DIR}/${TUMOR_SEQ_ID}.mantis.msi
        #==========================================================================================#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MANTIS ] Run FINISHED." >> ${LOGFILE}
        #==========================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MSI CALCULATION FINISHED." >> ${LOGFILE}
    echo " " >> ${LOGFILE}
    echo "WES MSI Calculation DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 
