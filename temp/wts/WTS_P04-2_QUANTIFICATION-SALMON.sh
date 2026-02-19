
#! /bin/bash

#---| SCRIPT USAGE |-------------------------------------------------------------------------------#
    usage() {
        echo "  Usage: WTS_P04-2_QUANTIFICATION-SALMON.sh [ARGUMENTS]...

        --baseDir           Base folder for processing and analysis. 
                            Generally specified by NGS data type.
                            default for WTS = '/data/wts'

        --batchID           [ REQUIRED ] NGS data batch-id. 
                            Batch-id also used as folder name of data batch.

        --seqID             [ REQUIRED ] Sample's seq-id. 
                            Seq-id also used as sample-data folder name.

        --threads           The number of cores to use
                            default = 10

        --fastqPrefix       Input fastq file prefix
                            
        --stranded          NGS data produced via stranded-library preparation 
                            or unstranded-library preparation.
                            default = false

        --genomeAssembly    Reference genome assembly version.
                            'hg19' and 'hg38' enable.
                            default = hg19

        --paired            Sequencing type is paired-end sequencing or not.
                            default = true (paired-end)

        --outputPrefix      Output folder name prefix

        -h | --help         Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARSE |-------------------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,threads:,baseDir:,fastqPrefix:,stranded:,genomeAssembly:,paired:,outputPrefix:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1, default = /data/wts" ; exit 2 ;;
                *)  BaseDir=$2 ; shift 2 ;;
            esac ;;
        --batchID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BatchID=$2 ; shift 2 ;;
            esac ;;
        --seqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac  ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; exit 2 ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --fastqPrefix )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 150" ; exit 2 ;;
                *)  FastqPrefix=$2 ; shift 2 ;;
            esac ;;        
        --stranded )
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
        --outputPrefix )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  OutputPrefix=$2 ; shift 2 ;;
            esac ;; 
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| MANUAL PARAMS |------------------------------------------------------------------------------#
    # BaseDir="/data/wts"
    # BatchID="WTS_25_02"
    # SeqID="WTS_25_02_03"
    # Threads=15
    # Stranded="false"
    # GenomeAssembly="mm10"
    # Paired="true"
    # FastqPrefix=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed
    # OutputPrefix=${SeqID}.human
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
    if [ -z ${FastqPrefix} ]; then FastqPrefix=${BaseDir}/${BatchID}/${SeqID}/fastq/${SeqID}.trimmed ; fi
    if [ -z ${Stranded} ]; then Stranded=false ; fi
    if [ -z ${GenomeAssembly} ]; then GenomeAssembly=hg38; fi
    if [ -z ${Paired} ]; then Paired=true; fi
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    BamDir=${DataDir}/${SeqID}/bam
    QuantDir=${DataDir}/${SeqID}/quant/salmon
    LogDir=${DataDir}/${SeqID}/log
    QcDir=${DataDir}/${SeqID}/qcfiles
    ##
    if [ ! -d ${QuantDir} ]; then mkdir -p ${QuantDir}; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES AND INDEX |-----------------------------------------------------------------------#
    SalmonIndex=/storage/references_and_index/${GenomeAssembly}/salmon-index/decoy-aware
    GtfFile=/storage/references_and_index/${GenomeAssembly}/gtf/${GenomeAssembly}.gtf
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RunSalmon="/storage/apps/salmon-1.10.0/bin/salmon"
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.salmon.quant.${GenomeAssembly}.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WTS RNAseq Quantification: SALMON log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${DataDir}" >> ${LogFile}
    echo " SEQ ID               : ${SeqID}" >> ${LogFile}
    echo " READ LENGTH          : ${ReadLength}" >> ${LogFile}
    echo " REFERENCE GENOME     : ${GenomeAssembly}" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | SALMON Quantification START." >> ${LogFile}
#==================================================================================================# 

#---| RUN SALMON QUANTIFICATION |------------------------------------------------------------------#
    if [ "${Stranded}" = "true" ]; then 
        LibraryType="ISR"
    else
        LibraryType="IU"
    fi

    ${RunSalmon} quant \
        --index ${SalmonIndex} --geneMap ${GtfFile} --libType ${LibraryType} \
        --mates1 ${FastqPrefix}_R1.fastq.gz --mates2 ${FastqPrefix}_R2.fastq.gz \
        --threads ${Threads} --gcBias --seqBias --posBias \
        --output ${QuantDir}/${OutputPrefix}

#--------------------------------------------------------------------------------------------------#


#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | SALMON Quantification FINISHED." >> ${LogFile}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | RNAseq Quantification: SALMON DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 





