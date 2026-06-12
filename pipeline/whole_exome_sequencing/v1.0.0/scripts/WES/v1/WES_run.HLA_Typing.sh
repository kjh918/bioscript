#!/bin/bash
# [MODIFIED] Added strict error handling.
set -euo pipefail

# SCRIPT VERSION : v0.4
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com
# MODIFIED       : Applied strict resource controls (-tc for razers3) to prevent OOM server crashes.

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.HLA_Typing.sh  [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID
        [ -t | --inputType ]      : input NGS type (DNA or RNA)
        [ --baseDir ]             : Base Work Dir
        [ --threads ]             : N-CPU threads (WARNING: HLA-LA/razers3 use heavy RAM per thread. Recommend 4~8)
        [ -h | --help ]           : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    INPUT_TYPE=""
    BASE_DIR=""
    THREADS=""
    SEQ_FOLDER=""
    SEQ_ID=""

    ARGS=$(getopt -a -o d:s:t:h --long seqFolder:,seqID:,inputType:,baseDir:,threads:,help -- "$@" )
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
                -* | --* | "") echo "[ERROR] no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        -s | --seqID )
            case "$2" in
                -* | --* | "") echo "[ERROR] no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        -t | --inputType )
            case "$2" in
                -* | --* | "") echo "[ERROR] no value at $1" >&2 ; exit 2 ;;
                *)  INPUT_TYPE=$2 ; shift 2 ;;
            esac ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "[ERROR] no value at $1" >&2 ; exit 2 ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "[ERROR] no value at $1" >&2 ; exit 2 ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done

    if [ -z "$SEQ_FOLDER" ] || [ -z "$SEQ_ID" ] || [ -z "$INPUT_TYPE" ] || [ -z "$BASE_DIR" ] || [ -z "$THREADS" ]; then
        echo "[ERROR] Missing required arguments. All parameters must be explicitly provided." >&2
        usage >&2
        exit 1
    fi
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    HLA_RES_DIR=${DATA_DIR}/${SEQ_ID}/hla
    OPTITYPE_RES_DIR=${HLA_RES_DIR}/optitype
    HLA_LA_RES_DIR=${HLA_RES_DIR}/hla-la
    FASTQ_DIR=${DATA_DIR}/${SEQ_ID}/fastq
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam  
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    
    mkdir -p ${HLA_RES_DIR} ${OPTITYPE_RES_DIR} ${HLA_LA_RES_DIR} ${LOG_DIR}
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.HLA.typing.$(date '+%Y%m%d').log
    touch ${LOGFILE}
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_CONDA_HLA="${RUN_SINGULARITY}/conda_HLA.sif"
    RUN_HLA_LA="/storage/apps/HLA-LA-1.0.3/src/HLA-LA.pl"
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${INPUT_TYPE}" = "DNA" ]; then
        OPTITYPE_REF_FASTA=/storage/references_and_index/optitype_references/hla_reference_dna.fasta
    elif [ "${INPUT_TYPE}" = "RNA" ]; then
        OPTITYPE_REF_FASTA=/storage/references_and_index/optitype_references/hla_reference_rna.fasta
    else
        echo "[ERROR] Invalid inputType: ${INPUT_TYPE}. Must be DNA or RNA." >&2
        exit 1
    fi

    if [ ! -f "${OPTITYPE_REF_FASTA}" ]; then
        echo "[ERROR] Reference file not found: ${OPTITYPE_REF_FASTA}" >&2
        exit 1
    fi
#--------------------------------------------------------------------------------------------------#

    echo " " >> ${LOGFILE}
    echo "WES HLA-Typing log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " MAX THREADS          : ${THREADS}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}   

    echo "$(date '+%Y-%m-%d %H:%M:%S') | HLA-TYPING START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run OptiType START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - Fished-BAM generation Start (Threads: ${THREADS})." >> ${LOGFILE}


    OPTITYPE_REGION="6:28500000-33500000"

    ${RUN_CONDA_HLA} samtools view -@ 8 -b \
        ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
        ${OPTITYPE_REGION} | \
    ${RUN_CONDA_HLA} samtools fastq -@ 8 \
        -1 ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R1.fastq \
        -2 ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R2.fastq \
        -0 /dev/null \
        -s /dev/null \
        -n \
        -
#--------------------------------------------------------------------------------------------------#

    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - HLA-region FASTQ extraction Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - OptiType Start." >> ${LOGFILE}

#---| OPTITYPE : HLA-TYPING |----------------------------------------------------------------------#
    if [ "${INPUT_TYPE}" = "DNA" ]; then
        ${RUN_CONDA_HLA} OptiTypePipeline.py \
            --input ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R1.fastq ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R2.fastq \
            --dna \
            --outdir ${OPTITYPE_RES_DIR} \
            --prefix ${SEQ_ID}.optitype
    else
        ${RUN_CONDA_HLA} OptiTypePipeline.py \
            --input ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R1.fastq ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R2.fastq \
            --rna \
            --outdir ${OPTITYPE_RES_DIR} \
            --prefix ${SEQ_ID}.optitype
    fi
#--------------------------------------------------------------------------------------------------#    

    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run OptiType FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run HLA-LA START." >> ${LOGFILE}

#---| HLA-LA : HLA-TYPING |------------------------------------------------------------------------#
    ${RUN_HLA_LA} \
        --BAM ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID ${SEQ_ID} \
        --maxThreads ${THREADS} \
        --workingDir ${HLA_LA_RES_DIR} 
#--------------------------------------------------------------------------------------------------#    

    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run HLA-LA FINISHED." >> ${LOGFILE}
    echo " " >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | HLA-TYPING DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}