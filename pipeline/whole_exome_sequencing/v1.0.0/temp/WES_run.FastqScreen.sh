
#!/bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.FastqScreen.sh [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID (Used as sample's processing product folder name)
        [ -h | --help ]           : Print this usage 
        [ --baseDir ]             : Base Work Dir
        [ --threads ]             : N-CPU threads
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    THREADS=15
    BASE_DIR=/data/wes
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:l:h: --long seqFolder:,seqID:,baseDir:,threads:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1, default = 18" ; shift ;;
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
    RAW_FQ_DIR=${BASE_DIR}/fastq
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    QC_DIR=${DATA_DIR}/${SEQ_ID}/qcfiles
    ##
    if [ ! -d ${QC_DIR} ]; then 
        mkdir -p ${QC_DIR}
        mkdir -p ${QC_DIR}/qc_res
    fi
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    SCRIPT_DIR=/storage/home/kangsm/runScripts
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_FASTQ_SCREEN="${RUN_SINGULARITY}/fastqScreen-0.15.3.sif fastq_screen"
    RUN_CONFIG=${SCRIPT_DIR}/NGS_config.FastqScreen.conf
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.fastq.screen.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES FastqScreen log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " REFERENCE SPECIES    : HUMAN, MOUSE, RAT" >> ${LOGFILE}
    echo "                        E.Coli, PhiX, Pseudomonas, Mycoplasm" >> ${LOGFILE}
    echo "                        Adpaters, Vectors" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqScreen START." >> ${LOGFILE}
#==================================================================================================#

#---| RUN FASTQ-SCREEN |---------------------------------------------------------------------------#
    ${RUN_FASTQ_SCREEN} \
        --aligner bowtie2 \
        --conf ${RUN_CONFIG} \
        --outdir ${QC_DIR} \
        --threads ${THREADS} \
        --force \
        ${RAW_FQ_DIR}/${SEQ_ID}_R1.fastq.gz
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqScreen FINISHED." >> ${LOGFILE}
    echo " " >> ${LOGFILE}
    echo "WES FastqScreen PROCESS DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 


