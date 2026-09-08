
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A03
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#


#---| Input Options |------------------------------------------------------------------------------#
    usage() {echo "
        Usage: WGS_M01-1_FASTQ-SCREEN.sh [ ARGUMENTS ]...
        
        --BASE_DIR      [ REQUIRED ] analysis base folder by NGS data type
        
        --BATCH_ID      [ REQUIRED ] dataset batch-id

        --SEQ_ID        [ REQURIED ] sample-id 

        --THREADS       parallel processing cpu cores 
                        default = 15 

        -h | --help     Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long BATCH_ID:,SEQ_ID:,BASE_DIR:,THREADS:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --BATCH_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        -s | --SEQ_ID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --BASE_DIR )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --THREADS )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| CHECK PARAMS and SET DEFAULT |----------------------------------------------------------------
    if [ -z "${Threads}" ]; then Threads=15; fi
#--------------------------------------------------------------------------------------------------#

#---| Manual Params |-------------------------------------------------------------------------------
    # BaseDir=
    # BatchID=
    # SeqID=
    # Threads=
#--------------------------------------------------------------------------------------------------- 

#---| FOLDERS |------------------------------------------------------------------------------------#
    RawFastqDir=${BaseDir}/${BatchID}/fastq
    DataDir=${BaseDir}/${BatchID}
    LogDir=${DataDir}/${SeqID}/log
    QcDir=${DataDir}/${SeqID}/qcfiles
    ##
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    Run_FastqScreen="singularity exec -B /storage,/data /storage/images/fastqScreen-0.15.3.sif fastq_screen"
    RunFastqScreenConfig=/storage/home/kangsm/runScripts/WGS/WGS_M01-2_FASTQ-SCREEN-CONFIG.conf
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.fastq.screen.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WGS FastqScreen log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${DataDir}" >> ${LogFile}
    echo " SEQ ID               : ${SeqID}" >> ${LogFile}
    echo " REFERENCE SPECIES    : HUMAN, MOUSE, RAT" >> ${LogFile}
    echo "                        E.Coli, PhiX, Pseudomonas, Mycoplasm" >> ${LogFile}
    echo "                        Adpaters, Vectors" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqScreen START." >> ${LogFile}
#==================================================================================================#

#---| RUN FASTQ-SCREEN |---------------------------------------------------------------------------#
    ${Run_FastqScreen} \
        --aligner bowtie2 \
        --conf ${RunFastqScreenConfig} \
        --outdir ${QcDir} \
        --THREADS ${Threads} \
        --force \
        ${RawFastqDir}/${SeqID}_R1.fastq.gz
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqScreen FINISHED." >> ${LogFile}
    echo " " >> ${LogFile}
    echo "WGS FastqScreen PROCESS DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 


