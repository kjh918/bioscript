
#!/bin/bash

#==============================================================================#
# ScriptID : GMC-NGS-B02-A02
# Author : kangsm
# version : v1.0 
# lastest update date : 2025.12.26
#==============================================================================#

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "
        Usage: WGS_P02-1_FASTQ-QC.sh [ ARGUMENTS ]...
        
        --BASE_DIR          [ REQUIRED ] analysis base folder by NGS data type

        --BATCH_ID          [ REQUIRED ] dataset batch-id
        
        --SEQ_ID            [ REQURIED ] sample-id 
        
        --THREADS           parallel processing cpu cores 
                            default = 15
        
        --RUN_MULTIQC       run MultiQC in this step or not
                            default = false
        
        -h | --help         Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#


#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:h: --long BATCH_ID:,SEQ_ID:,BASE_DIR:,THREADS:,RUN_MULTIQC:,help -- "$@" )
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
                    -* | --* | "") echo "no value at $1, default = /data/wes " ; shift ;;
                    *)  BaseDir=$2 ; shift 2 ;;
                esac ;;
            --THREADS )
                case "$2" in
                    -* | --* | "") echo "no value at $1, default = 14" ; shift ;;
                    *)  Threads=$2 ; shift 2 ;;
                esac ;;
            --RUN_MULTIQC )
                case "$2" in
                    -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                    *)  Run_Multiqc=$2 ; shift 2 ;;
                esac ;;
            -h | --help )
            usage >&2 ; exit 2 ;;   
            --) shift ; break ;;
            *)  usage >&2 ; exit 2 ;;
        esac
    done
#--------------------------------------------------------------------------------------------------#

#---| CHECK PARAMS and SET DEFAULT |----------------------------------------------------------------
    if [[ -z "${Threads}" ]]; then Threads=15; fi
    if [[ -z "${Run_Multiqc}" ]]; then Run_Multiqc="false"; fi
#--------------------------------------------------------------------------------------------------#

#---| Manual Params |------------------------------------------------------------------------------#
    # BaseDir="/data/wgs"
    # BatchID="WGS-26-01"
    # SeqID="WGS-26-01_02"
    # Threads=15
    # Run_Multiqc="false"
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    RawFastqDir=${DataDir}/fastq
    FastqDir=${DataDir}/${SeqID}/fastq
    QcDir=${DataDir}/${SeqID}/qcfiles
    LogDir=${DataDir}/${SeqID}/log
    ##
    mkdir -p ${FastqDir}
    mkdir -p ${LogDir}
    mkdir -p ${QcDir}/qc_res
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunFastqc="${RunSingularity}/fastqc-0.12.1.sif fastqc"
    RunFastp="${RunSingularity}/fastp-0.23.4.sif fastp"
    RunMultiqc="${RunSingularity}/multiqc-1.18.sif multiqc"
#--------------------------------------------------------------------------------------------------#

#---| RunLog |--------------------------------------------------------------------------------------
    RunLog=${LogDir}/${SeqID}.wgs.fastq-qc.$(date '+%Y%m%d').log
    echo " " > ${RunLog}
    echo " WGS FASTQ-QC RUN LOG" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " RUN DATE     $(date '+%Y-%m-%d %H:%M:%S')" >> ${RunLog}
    echo " BATCH ID     ${BatchID}" >> ${RunLog}
    echo " SEQ ID       ${SeqID}" >> ${RunLog}
    echo "----------------------------------------------------------------------" >> ${RunLog}
    echo " " >> ${RunLog}   
#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqQC START." >> ${LogFile}
#==================================================================================================#

#---| RUN FASTQC |---------------------------------------------------------------------------------#    

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run Fastqc ----------------------------------------------###" >> ${RunLog}

    export LC_ALL=C
    export LANG=C

    FASTQC_ARGS=()
    FASTQC_ARGS+=("--extract")
    FASTQC_ARGS+=("-t ${Threads}")
    FASTQC_ARGS+=("-o ${QcDir}")
    FASTQC_ARGS+=("${RawFastqDir}/${SeqID}_R1.fastq.gz")
    FASTQC_ARGS+=("${RawFastqDir}/${SeqID}_R2.fastq.gz")

    ${RunFastqc} ${FASTQC_ARGS[@]} 

    echo " " >> ${RunLog}   
#--------------------------------------------------------------------------------------------------#    
        
#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastQC FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | fastp TRIMMING START." >> ${LogFile}
#==================================================================================================#

#---| fastp CUT-OFF |------------------------------------------------------------------------------#

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run fastp -----------------------------------------------###" >> ${RunLog}

    FASTP_ARGS=()
    FASTP_ARGS+=("--in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz")
    FASTP_ARGS+=("--in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz")
    FASTP_ARGS+=("--json ${QcDir}/${SeqID}.fastp.json")
    FASTP_ARGS+=("--html ${QcDir}/${SeqID}.fastp.html")
    FASTP_ARGS+=("--out1 ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz")
    FASTP_ARGS+=("--out2 ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz")
    FASTP_ARGS+=("--thread ${Threads}")
    FASTP_ARGS+=("--length_required 50")
    FASTP_ARGS+=("--cut_tail")
    FASTP_ARGS+=("--cut_tail_window_size 4 ")
    FASTP_ARGS+=("--cut_tail_mean_quality 20")
    FASTP_ARGS+=("--trim_poly_g")
    FASTP_ARGS+=("--detect_adapter_for_pe")
    FASTP_ARGS+=("--n_base_limit 5 ")

    ${RunFastp} ${FASTP_ARGS[@]} 
    
    echo " " >> ${RunLog}   
#---------------------------------------------------------------------------------------------------

#---| RUN MULTIQC OPTIONS |------------------------------------------------------------------------#    

    if [ "${Run_Multiqc}" = "true" ]; then

        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Run Multi-QC --------------------------------------------###" >> ${RunLog}

        MULTIQC_ARGS=()
        MULTIQC_ARGS+=("--force")
        MULTIQC_ARGS+=("--filename ${SeqID}_qc")
        MULTIQC_ARGS+=("--outdir ${QcDir}/qc_res")
        MULTIQC_ARGS+=("--config /storage/home/kangsm/runScripts/WGS/WGS_P05-2_QC-POST_PROCESSING-MULTIQC-CONFIG.yaml")
        MULTIQC_ARGS+=("--data-dir ${QcDir}")

        ${RunMultiqc} ${MULTIQC_ARGS[@]} 

        echo " " >> ${RunLog}   
    
    else
        echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] SKIP Multi-QC Running" >> ${RunLog}
    fi
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] WGS FASTQ-QC PROCESS FINISHED." >> ${RunLog}
    echo " " >> ${RunLog}
#==================================================================================================# 

