
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A03
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#


#---| Input Options |------------------------------------------------------------------------------#
    usage() {echo " Usage: WES_M01-1_FASTQ-SCREEN.sh [ ARGUMENTS ]...

        --baseDir       Base folder for processing and analysis. 
                        Generally specified by NGS data type.
                        default for WTS = '/data/wes'
      
        --batchID       [ REQUIRED ] NGS data batch-id. 
                        Batch-id also used as folder name of data batch.

        --seqID         [ REQUIRED ] Sample's seq-id. 
                        Seq-id also used as sample-data folder name.
        
        --threads       The number of cores to use
                        default = 10
                        
        -h | --help     Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,baseDir:,threads:,help -- "$@" )
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
        -s | --seqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --threads )
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
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RunScript=/storage/home/kangsm/runScripts
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunFastqScreen="${RunSingularity}/fastqScreen-0.15.3.sif fastq_screen"
    FastqScreenConfig=${RunScript}/NGS_config.FastqScreen.conf
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.fastq.screen.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WES FastqScreen log" >> ${LogFile}
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
    ${RunFastqScreen} \
        --aligner bowtie2 \
        --conf ${FastqScreenConfig} \
        --outdir ${QcDir} \
        --threads ${Threads} \
        --force \
        ${RawFastqDir}/${SeqID}_R1.fastq.gz
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqScreen FINISHED." >> ${LogFile}
    echo " " >> ${LogFile}
    echo "WES FastqScreen PROCESS DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 


