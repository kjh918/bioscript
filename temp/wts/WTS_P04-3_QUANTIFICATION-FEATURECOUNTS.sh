
#! /bin/bash

#---| SCRIPT USAGE |-------------------------------------------------------------------------------#
    usage() {
        echo "  Usage: WTS_P03-2_QUANTIFICATION-SALMON.sh [ARGUMENTS]...

        --baseDir           Base folder for processing and analysis. 
                            Generally specified by NGS data type.
                            default for WTS = '/data/wts'

        --batchID           [ REQUIRED ] NGS data batch-id. 
                            Batch-id also used as folder name of data batch.

        --seqID             [ REQUIRED ] Sample's seq-id. 
                            Seq-id also used as sample-data folder name.

        --bamFile           Input BAM file. 
                            If not specified, use Aligned.out.bam as default BAM file.

        --threads           The number of cores to use
                            default = 10

        --stranded          NGS data produced via stranded-library preparation 
                            or unstranded-library preparation.
                            default = false

        --genomeAssembly    Reference genome assembly version.
                            'hg19' and 'hg38' enable.
                            default = hg19

        --paired            Sequencing type is paired-end sequencing or not.
                            default = true (paired-end)

        --outputPrefix      Output file name prefix

        -h | --help        Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARSE |-------------------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,threads:,baseDir:,bamFile:,stranded:,genomeAssembly:,paired:,outputPrefix:,help -- "$@" )
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
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
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
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  Threads=$2 ; shift 2 ;;
            esac ;;
        --bamFile )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  BamFile=$2 ; shift 2 ;;
            esac ;;        
        --stranded )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  Stranded=$2 ; shift 2 ;;
            esac ;; 
        --genomeAssembly )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  GenomeAssembly=$2 ; shift 2 ;;
            esac ;; 
        --paired )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
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
    # BamFile=${BaseDir}/${BatchID}/${SeqID}/bam/${SeqID}.human.xenofilter.bam
    # Stranded="false"
    # GenomeAssembly="hg38"
    # Paired="true"
    # OutputPrefix=${SeqID}.human
#--------------------------------------------------------------------------------------------------#

#---| ARGS CHECK |---------------------------------------------------------------------------------#
    if [ -z ${BaseDir} ]; then BaseDir=/data/wts ; fi
    if [ -z ${BatchID} ]; then 
        echo "NO BATCH-ID. Required."    
        exit 0
    fi
    if [ -z ${SeqID} ]; then 
        echo "No SEQ-ID. Required."
        exit 0
    fi
    if [ -z ${OutputPrefix} ]; then OutputPrefix=${SeqID}; fi
    if [ -z ${Threads} ]; then Threads=10 ; fi
    if [ -z ${Stranded} ]; then Stranded=true ; fi
    if [ -z ${GenomeAssembly} ]; then GenomeAssembly=hg38; fi
    if [ -z ${Paired} ]; then Paired=true; fi

#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    BamDir=${DataDir}/${SeqID}/bam
    QuantDir=${DataDir}/${SeqID}/quant/featurecounts
    LogDir=${DataDir}/${SeqID}/log
    QcDir=${DataDir}/${SeqID}/qcfiles
    ##
    if [ ! -d ${QuantDir} ]; then mkdir -p ${QuantDir}; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
    if [ ! -d ${QcDir} ]; then mkdir -p ${QcDir}/qc_res; fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES AND INDEX |-----------------------------------------------------------------------#
    GtfFile=/storage/references_and_index/${GenomeAssembly}/gtf/${GenomeAssembly}.gtf
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RunFeatureCounts="/storage/apps/subread-2.1.1/bin/featureCounts"
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.featurecount.quant.${GenomeAssembly}.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WTS RNAseq Quantification: featureCount log" >> ${LogFile}
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
    echo "$(date '+%Y-%m-%d %H:%M:%S') | featureCount Quantification START." >> ${LogFile}
#==================================================================================================# 

#---| RUN FeatureCounts |--------------------------------------------------------------------------#
    ${RunFeatureCounts} \
        -F GTF -t exon -g gene_id -s 2 -p --countReadPairs -T ${Threads} \
        -a ${GtfFile} \
        -o ${QuantDir}/${OutputPrefix}.featurecounts \
        ${BamFile}
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | featureCount Quantification FINISHED." >> ${LogFile}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | RNAseq Quantification: featureCount DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 





