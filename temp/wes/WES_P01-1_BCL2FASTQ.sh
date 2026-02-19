
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A01
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.BCL2FASTQ.sh [ ARGUMENTS ]...
        --bcl  BCL data output folder
        --outputFastq  FASTQ output folder
        --baseDir  base work folder default=/data/wes
        --seqFolder  SeqFolder ID 
        --threadsR  N-CPU threads for reading. default=10
        --threadsP  N-CPU threads for processing. default=20
        --threadsW  N-CPU threads for writing. default=10
        -h | --help  print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    READ_THREADS=10
    PROCESS_THREADS=20
    WRITE_THREADS=10
    BASE_DIR=/data/wes
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long bcl:,outputFastq:,baseDir:,seqFolder:,threadsR:,threadsP:,threadsW:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --bcl )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BCL_DIR=$2 ; shift 2 ;;
            esac ;;
        --outputFastq )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  OUTPUT_FASTQ_DIR=$2 ; shift 2 ;;
            esac ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes " ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes " ; shift ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        --threadsR )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  READ_THREADS=$2 ; shift 2 ;;
            esac ;;
        --threadsP )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 20" ; shift ;;
                *)  PROCESS_THREADS=$2 ; shift 2 ;;
            esac ;;
        --threadsW )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 10" ; shift ;;
                *)  WRITE_THREADS=$2 ; shift 2 ;;
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
    RAW_FQ_TMP=${RAW_FQ_DIR}/tmp
    if [ ! -d ${RAW_FQ_TMP} ]; then mkdir -p ${RAW_FQ_TMP}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
#--------------------------------------------------------------------------------------------------#

#---| RUN BCL2FASTQ |---------------------------------------------------------------------------------#    
    ${RUN_SINGULARITY}/bcl2fastq-v2.20.sif bcl2fastq -R ${BCL_DIR} -o ${OUTPUT_FASTQ_DIR} -r ${READ_THREADS} -p ${PROCESS_THREADS} -w ${WRITE_THREADS} --no-lane-splitting
#--------------------------------------------------------------------------------------------------#    

#---| FASTQ RENAME |-------------------------------------------------------------------------------#
    echo "copy raw-fastq files to temp folder..."
    cp ${OUTPUT_FASTQ_DIR}/${SEQ_FOLDER}_*.fastq.gz ${RAW_FQ_TMP}/
    echo "rename and move fastq-files..."
    rename 's/_S[0-9][0-9]//' ${RAW_FQ_TMP}/${SEQ_FOLDER}_*.fastq.gz
    rename 's/_S[0-9]//' ${RAW_FQ_TMP}/${SEQ_FOLDER}_*.fastq.gz
    rename 's/_001.fastq.gz/.fastq.gz/' ${RAW_FQ_TMP}/${SEQ_FOLDER}_*.fastq.gz
#--------------------------------------------------------------------------------------------------#    

#---| MOVE FASTQ |---------------------------------------------------------------------------------#
    mv ${RAW_FQ_TMP}/${SEQ_FOLDER}_*.fastq.gz ${RAW_FQ_DIR}/
    echo "BCL To FASTQ Process DONE."
#--------------------------------------------------------------------------------------------------#    

