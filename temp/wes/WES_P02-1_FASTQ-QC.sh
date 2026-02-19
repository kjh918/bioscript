
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A02
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.FastqQC.sh [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID (Used as sample's processing product folder name)
        [ --baseDir ]             : Base Work Dir
        [ --readLength ]          : NGS Read Length 
        [ --threads ]             : N-CPU threads
        [ --runMultiqc ]          : run MultiQC in this step or not ( default = false )
        [ -h | --help ]           : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    THREADS=15
    BASE_DIR=/data/wes
    RUN_MULTIQC="false"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:h: --long seqFolder:,seqID:,baseDir:,readLength:,threads:,runMultiqc:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1, default = /data/wes " ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --readLength )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  READ_LENGTH=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 14" ; shift ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        --runMultiqc )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  RUN_MULTIQC=$2 ; shift 2 ;;
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
    ## 
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    ##  
    FASTQ_DIR=${DATA_DIR}/${SEQ_ID}/fastq
    QC_DIR=${DATA_DIR}/${SEQ_ID}/qcfiles
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    ##
    if [ ! -d ${FASTQ_DIR} ]; then mkdir -p ${FASTQ_DIR}; fi
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
    if [ ! -d ${QC_DIR} ]; then 
        mkdir -p ${QC_DIR}
        mkdir -p ${QC_DIR}/qc_res
    fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.fastq.qc.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES FastqQC log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo " FASTQ Read Length    : ${READ_LENGTH}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqQC START." >> ${LOGFILE}
#==================================================================================================#

#---| RUN FASTQC |---------------------------------------------------------------------------------#    
    ${RUN_SINGULARITY}/fastqc-0.12.1.sif fastqc --extract \
        --threads ${THREADS} \
        --outdir ${QC_DIR} \
        ${RAW_FQ_DIR}/${SEQ_ID}_*.fastq.gz
#--------------------------------------------------------------------------------------------------#    
        
#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastQC FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | fastp TRIMMING START." >> ${LOGFILE}
#==================================================================================================#

#---| fastp CUT-OFF |------------------------------------------------------------------------------#
    LENGTH_CUTOFF=$(echo $READ_LENGTH | awk '$1 {print $1*.7}')

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ fastp ] min. length  cut-off = ${LENGTH_CUTOFF}" >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ fastp ] base quality cut-off = 20" >> ${LOGFILE}
#==================================================================================================# 

    ${RUN_SINGULARITY}/fastp-0.23.4.sif fastp \
        --in1 ${RAW_FQ_DIR}/${SEQ_ID}_R1.fastq.gz --in2 ${RAW_FQ_DIR}/${SEQ_ID}_R2.fastq.gz \
        --json ${QC_DIR}/${SEQ_ID}.fastp.json \
        --html ${FASTQ_DIR}/${SEQ_ID}.fastp.html \
        --length_required ${LENGTH_CUTOFF} \
        --qualified_quality_phred 20 \
        --thread ${THREADS} \
        --trim_poly_g \
        --detect_adapter_for_pe \
        --out1 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz --out2 ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | fastp TRIMMING FINISHED." >> ${LOGFILE}
#==================================================================================================# 

#--------------------------------------------------------------------------------------------------#    

#---| RUN MULTIQC OPTIONS |------------------------------------------------------------------------#    
    if [ "${RUN_MULTIQC}" = "true" ]; then
    #==============================================================================================#
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Run Option = ${RUN_MULTIQC}" >> ${LOGFILE}    
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Summary START. " >> ${LOGFILE}    
    #==============================================================================================#
    MULTIQC_CONFIG=/storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml
    ##
    ${RUN_SINGULARITY}/multiqc-1.18.sif multiqc --force \
        --filename ${SEQ_ID}_qc \
        --outdir ${QC_DIR}/qc_res \
        --config ${MULTIQC_CONFIG} \
        --data-dir ${QC_DIR}
    #==============================================================================================#
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Summary FINISHED. " >> ${LOGFILE}    
    #==============================================================================================#
    else
    #==============================================================================================#
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Run Option = ${RUN_MULTIQC}" >> ${LOGFILE}    
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Summary SKIPPED. " >> ${LOGFILE}    
    #==============================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES FastqQC PROCESS DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

