
#! /bin/bash



#---| SCRIPT USAGE |-------------------------------------------------------------------------------#
    usage() {
        echo "  Usage: WTS_P03-1_QUANTIFICATION-RSEM.sh [ARGUMENTS]...

        --baseDir           Base folder for processing and analysis. 
                            Generally specified by NGS data type.
                            default for WTS = '/data/wts'

        --batchID           [ REQUIRED ] NGS data batch-id. 
                            Batch-id also used as folder name of data batch.

        --seqID             [ REQUIRED ] Sample's seq-id. 
                            Seq-id also used as sample-data folder name.

        --threads           The number of cores to use
                            default = 10

        --readLength        The length of NGS reads in fastq file. 
                            '150' or '100' enable.
                            default = 150

        --stranded          NGS data produced via stranded-library preparation 
                            or unstranded-library preparation.
                            default = false

        --genomeAssembly    Reference genome assembly version.
                            'hg19' and 'hg38' enable.
                            default = hg19

        --paired            Sequencing type is paired-end sequencing or not.
                            default = true (paired-end)

        --xenograftPipeline Run as Xenograft pipeline. default = 'false'

        -h | --help       = Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARSE |-------------------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,threads:,baseDir:,readLength:,stranded:,genomeAssembly:,paired:,xenograftPipeline:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --BaseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wts" ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --BatchID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        --SeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; exit 2 ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --readLength )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 150" ; exit 2 ;;
                *)  ReadLength=$2 ; shift 2 ;;
            esac ;;        
        --standed )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = false" ; exit 2 ;;
                *)  Stranded=$2 ; shift 2 ;;
            esac ;; 
        --genomeAssembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; exit 2 ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;; 
        --paired )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  Paired=$2 ; shift 2 ;;
            esac ;; 
        --xenograftPipeline )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = true" ; exit 2 ;;
                *)  XenograftPipeline=$2 ; shift 2 ;;
            esac ;; 
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| MANUAL PARAMS |-----------------------------------------------------------------------------#
    # BaseDir="/data/wts"
    # BatchID=
    # SeqID=
    # Threads=15
    # ReadLength=100
    # Stranded="false"
    # GenomeAssembly="hg19"
    # Paired="true"
    # XenograftPipeline="false"
#--------------------------------------------------------------------------------------------------#

#---| ARGS CHECK |---------------------------------------------------------------------------------#
    if [ -z ${BatchID} ]; then 
        echo "NO BATCH-ID. Required."    
        exit 0 
    fi
    if [ -z ${SeqID} ]; then 
        echo "No SEQ-ID. Required."
        exit 0
    fi
    if [ -z ${BaseDir} ]; then BaseDir=/data/wts ; fi
    if [ -z ${Threads} ]; then Threads=10 ; fi
    if [ -z ${ReadLength} ]; then ReadLength=150 ; fi
    if [ -z ${Stranded} ]; then Stranded=false ; fi
    if [ -z ${GenomeAssembly} ]; then GenomeAssembly=hg19; fi
    if [ -z ${Paired} ]; then Paired=true; fi
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    BamDir=${DataDir}/${SeqID}/bam
    QuantDir=${DataDir}/${SeqID}/quant/rsem
    LogDir=${DataDir}/${SeqID}/log
    QcDir=${DataDir}/${SeqID}/qcfiles
    ##
    if [ ! -d ${QuantDir} ]; then mkdir -p ${QuantDir}; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES AND INDEX |-----------------------------------------------------------------------#
    RsemIndex=/storage/references_and_index/${GenomeAssembly}/star-rsem-index/PE-${ReadLength}/${GenomeAssembly}
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunRsem="${RunSingularity}/rsem-1.3.3.sif rsem-calculate-expression"
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.fastq.to.bam.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WTS RNAseq Quantification: RSEM log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${DataDir}" >> ${LogFile}
    echo " SEQ ID               : ${SeqID}" >> ${LogFile}
    echo " READ LENGTH          : ${ReadLength}" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | RSEM Quantification START." >> ${LogFile}
#==================================================================================================# 

#---| RUN STAR ALIGNMENT |-------------------------------------------------------------------------#
    if [ "${XenograftPipeline}" = "true" ]; then




    else
    
        if [ "${Stranded}" = "true" ]; then RsemStrandOption="reverse"; else RsemStrandOption="none"; fi
    
        if [ "${Paired}" = "true" ]; then
            ${RunRsem} --paired-end --no-bam-output --estimate-rspd \
                --strandedness ${RsemStrandOption} \
                --num-threads ${Threads} \
                --alignments ${BamDir}/${SeqID}.Aligned.toTranscriptome.out.bam \
                ${RsemIndex} \
                ${QuantDir}/${SeqID}     
        else
            ${RunRsem} --no-bam-output --estimate-rspd \
                --strandedness ${RsemStrandOption} \
                --num-threads ${Threads} \
                --alignments ${BamDir}/${SeqID}.Aligned.toTranscriptome.out.bam \
                ${RsemIndex} \
                ${QuantDir}/${SeqID}   
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | RSEM Quantification FINISHED." >> ${LogFile}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | RNAseq Quantification: RSEM DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 





