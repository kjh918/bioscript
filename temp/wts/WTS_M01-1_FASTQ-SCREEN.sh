
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A03
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#


#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WTS_M01-1_FASTQ-SCREEN.sh [ ARGUMENTS ]...
        -d | --batchID : data processing base directory
        -s | --SeqID   : sample ID (Used as sample's processing product folder name)
        --baseDir      : Base Work Dir
        --Threads       : N-CPU Threads
        -h | --help     : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:l:h: --long batchID:,SeqID:,baseDir:,Threads:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --batchID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        -s | --SeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --Threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 18" ; shift ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    RawFastqDir=${DataDir}/fastq
    LogDir=${DataDir}/${SeqID}/log
    QcDir=${DataDir}/${SeqID}/qcfiles
    ##
    if [ ! -d ${QcDir} ]; then 
        mkdir -p ${QcDir}/qc_res
    fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RunScriptDir=/storage/home/kangsm/runScripts/WTS
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    Run_FastqScreen="${RunSingularity}/fastqScreen-0.15.3.sif fastq_screen"
    RunFastqScreenConfig=${RunScriptDir}/WTS_M01-2_FASTQ-SCREEN-CONFIG.conf
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.fastq.screen.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WTS FastqScreen log" >> ${LogFile}
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
        --Threads ${Threads} \
        --force \
        ${RawFastqDir}/${SeqID}_R1.fastq.gz
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqScreen FINISHED." >> ${LogFile}
    echo " " >> ${LogFile}
    echo "WTS FastqScreen PROCESS DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 


