
#!/bin/bash
 
# version : release-1.0

#---| Input Options |------------------------------------------------------------------------------#
    usage() { echo "    Usage: WTS_P01-1_FASTQ-QC.sh [ ARGUMENTS ]...

        --baseDir                   Base folder for processing and analysis. 
                                    Generally specified by NGS data type.
                                    default for WTS = '/data/wts'

        --batchID                   [ REQUIRED ] NGS data batch-id. 
                                    Batch-id also used as folder name of data batch.

        --seqID                     [ REQUIRED ] Sample's seq-id. 
                                    Seq-id also used as sample-data folder name.

        --readLength                The length of NGS reads in fastq file.
                                    default = 150

        --readLengthFilterCutoff    Minimum required read length after trimming.
                                    ( recommend : set to 70% of initial NGS data read length )

        --threads                   The number of cores to use
                                    default = 10

        --runMultiqc                Run Multi-QC to summarize FASTQ-QC results.
                                    It's also enable to run summarization of QC reuslts
                                    after all processing finished. (recommended)
                                    default = false

        --seqType                   NGS sequencing type. 
                                    'SE' for single-end and 'PE' for paired-end.
                                    default = PE

        -h | --help                 Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    # Threads=15
    # BaseDir=/data/wes
    # RunMultiqc="false"
    # SeqType="PE"
    # ReadLengthFilterCutoff
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,baseDir:,readLength:,readLengthFilterCutoff:,threads:,runMultiqc:,seqType:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes " ; shift ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
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
        --readLength )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  ReadLength=$2 ; shift 2 ;;
            esac ;;
        --readLengthFilterCutoff )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  ReadLengthFilterCutoff=$2 ; shift 2 ;;
            esac ;;

        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 14" ; shift ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --runMultiqc )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; shift ;;
                *)  RunMultiqc=$2 ; shift 2 ;;
            esac ;;
        --seqType )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = PE" ; shift ;;
                *)  SeqType=$2 ; shift 2 ;;
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
    FastqDir=${DataDir}/${SeqID}/fastq
    QcDir=${DataDir}/${SeqID}/qcfiles
    LogDir=${DataDir}/${SeqID}/log
    if [ ! -d ${FastqDir} ]; then mkdir -p ${FastqDir}; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RunSingularity="singularity exec -B /storage,/data /storage/images"
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.fastq.qc.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WTS FastqQC log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${DataDir}" >> ${LogFile}
    echo " SEQ ID               : ${SeqID}" >> ${LogFile}
    echo " FASTQ Read Length    : ${RreadLength}" >> ${LogFile}
    echo " SEQUENCING TYPE      : ${SeqType}" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqQC START." >> ${LogFile}
#==================================================================================================#

#---| RUN FASTQC |---------------------------------------------------------------------------------#    
    ${RunSingularity}/fastqc-0.12.1.sif fastqc --extract \
        --threads ${Threads} \
        --outdir ${QcDir} \
        ${RawFastqDir}/${SeqID}*.fastq.gz
#--------------------------------------------------------------------------------------------------#    
        
#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastQC FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | fastp TRIMMING START." >> ${LogFile}
#==================================================================================================#

#---| fastp CUT-OFF |------------------------------------------------------------------------------#
    #LENGTH_CUTOFF=$(echo $RreadLength | awk '$1 {print $1*.7}')
    # In RNAseq set this value as 50.
#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ fastp ] min. length  cut-off = ${LENGTH_CUTOFF}" >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | [ fastp ] base quality cut-off = 20" >> ${LogFile}
#==================================================================================================# 

    if [ "${SeqType}" = "SE" ]; then
        ${RunSingularity}/fastp-0.23.4.sif fastp \
            --in1 ${RawFastqDir}/${SeqID}.fastq.gz \
            --json ${QcDir}/${SeqID}.fastp.json \
            --html ${FastqDir}/${SeqID}.fastp.html \
            --length_required ${ReadLengthFilterCutoff} \
            --qualified_quality_phred 20 \
            --thread ${Threads} \
            --trim_poly_g \
            --out1 ${FastqDir}/${SeqID}.trimmed.fastq.gz
    else
        ${RunSingularity}/fastp-0.23.4.sif fastp \
            --in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz --in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz \
            --json ${QcDir}/${SeqID}.fastp.json \
            --html ${FastqDir}/${SeqID}.fastp.html \
            --length_required ${ReadLengthFilterCutoff} \
            --qualified_quality_phred 20 \
            --thread ${Threads} \
            --trim_poly_g \
            --detect_adapter_for_pe \
            --out1 ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz --out2 ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz
    fi
#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | fastp TRIMMING FINISHED." >> ${LogFile}
#==================================================================================================# 

#--------------------------------------------------------------------------------------------------#    

#---| RUN MULTIQC OPTIONS |------------------------------------------------------------------------#    
    if [ "${RunMultiqc}" = "true" ]; then
    #==============================================================================================#
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Run Option = ${RunMultiqc}" >> ${LogFile}    
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Summary START. " >> ${LogFile}    
    #==============================================================================================#
    MultiqcConfig=/storage/home/kangsm/runScripts/WTS/WTS_M01-2_QC-POST_PROCESSING-MULTIQC-CONFIG.yaml
    ##
    ${RunSingularity}/multiqc-1.18.sif multiqc --force \
        --filename ${SeqID}_qc \
        --outdir ${QcDir}/qc_res \
        --config ${MultiqcConfig} \
        --data-dir ${QcDir}
    #==============================================================================================#
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Summary FINISHED. " >> ${LogFile}    
    #==============================================================================================#
    else
    #==============================================================================================#
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Run Option = ${RunMultiqc}" >> ${LogFile}    
    echo "$(date '+%Y-%m-%d %H:%M:%S') | MultiQC Summary SKIPPED. " >> ${LogFile}    
    #==============================================================================================#
    fi
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WES FastqQC PROCESS DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 

