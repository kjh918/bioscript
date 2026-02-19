
#! /bin/bash

# version : release-1.0


#---| SCRIPT USAGE |-------------------------------------------------------------------------------#
    usage() { echo "  Usage: WTS_P05-1_XENOGRAFT-FASTQ-TO-BAM.sh [ARGUMENTS]...
        
        * IMPORTANT! 
        Xenograft RNAseq mapping supports ONLY 'hg38 and GRCm38 (mm10)' references and 'PAIRED-END' sequencing.
        And also NO option for strandness.

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
        
        -h | --help         Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARSE |-------------------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,threads:,baseDir:,readLength:,help -- "$@" )
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
        --readLength )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 150" ; exit 2 ;;
                *)  ReadLength=$2 ; shift 2 ;;
            esac ;;        
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| DEFAULT VALUES |-----------------------------------------------------------------------------#
    # BaseDir="/data/wts"
    # BatchID="WTS_25_02"
    # SeqID="WTS_25_02_03"
    # Threads=15
    # ReadLength=150
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
    if [ -z ${Threads} ]; then Threads=15 ; fi
    if [ -z ${ReadLength} ]; then ReadLength=150 ; fi
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DataDir=${BaseDir}/${BatchID}
    #---------------------------------------------------------------------------
    FastqDir=${DataDir}/${SeqID}/fastq
    BamDir=${DataDir}/${SeqID}/bam
    QcDir=${DataDir}/${SeqID}/qcfiles
    TmpDir=${DataDir}/${SeqID}/tmp
    LogDir=${DataDir}/${SeqID}/log
    #---------------------------------------------------------------------------
    if [ ! -d ${BamDir} ]; then mkdir -p ${BamDir}; fi
    if [ ! -d ${TmpDir} ]; then mkdir -p ${TmpDir}; fi
    if [ ! -d ${LogDir} ]; then mkdir -p ${LogDir}; fi
    if [ ! -d ${QcDir} ]; then 
        mkdir -p ${QcDir}
        mkdir -p ${QcDir}/qc_res
    fi
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES AND INDEX |-----------------------------------------------------------------------#
    HumanStarIndex=/storage/references_and_index/hg38/star-rsem-index/PE-${ReadLength}
    HumanGtfFile=/storage/references_and_index/hg38/gtf/hg38.gtf
    MouseStarIndex=/storage/references_and_index/mm10/star-rsem-index/PE-${ReadLength}
    MouseGtfFile=/storage/references_and_index/mm10/gtf/mm10.gtf
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunStar="${RunSingularity}/star-2.7.11.sif STAR"
    RunSamtools="${RunSingularity}/samtools-1.18.sif samtools"
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.xenograft.rnaseq.fastq.to.bam.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "Xenograft RNAseq Fastq-To-BAM log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${DataDir}" >> ${LogFile}
    echo " SEQ ID               : ${SeqID}" >> ${LogFile}
    echo " READ LENGTH          : ${ReadLength}" >> ${LogFile}
    echo " STAR MODE            : Two-Pass Mode" >> ${LogFile}
    echo " HUMAN GENOME         : hg38" >> ${LogFile}
    echo " MOUSE GENOME         : mm10" >> ${LogFile}
    echo " HUMAN/MOUSE FILTER   : xenofilteR" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Xenograft-FastqToBAM: STAR Alignment to 'HUMAN' START." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | STAR Alignment TwoPass Mode with HUMAN-REFERENCE" >> ${LogFile}
#==================================================================================================# 

#---| RUN STAR ALIGNMENT : HUMAN |-----------------------------------------------------------------#
    # mapping
    ${RunStar} \
        --genomeDir ${HumanStarIndex} \
        --readFilesIn ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
        --outFileNamePrefix ${BamDir}/${SeqID}.human. \
        --sjdbGTFfile ${HumanGtfFile} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --genomeLoad NoSharedMemory \
        --runThreadN ${Threads} \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NM NH HI AS nM XS \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outSAMunmapped Within \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimScoreJunctionNonGTAG 0 \
        --chimMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --outBAMcompression 0 \
        --chimOutType WithinBAM HardClip \
        --chimScoreDropMax 30 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3       

    # sort and index
    ${RunSamtools} sort -@ ${Threads} -o ${BamDir}/${SeqID}.human.Aligned.out.Sorted.bam ${BamDir}/${SeqID}.human.Aligned.out.bam
    ${RunSamtools} index  -@ ${Threads} -b ${BamDir}/${SeqID}.human.Aligned.out.Sorted.bam 
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | STAR Alignment TwoPass Mode with HUMAN-REFERENCE FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Xenograft-FastqToBAM: STAR Alignment to 'MOUSE' START." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | STAR Alignment TwoPass Mode with MOUSE-REFERENCE" >> ${LogFile}
#==================================================================================================# 

#---| RUN STAR ALIGNMENT : MOUSE |-----------------------------------------------------------------#
    # mapping
    ${RunStar} \
        --genomeDir ${MouseStarIndex} \
        --readFilesIn ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
        --outFileNamePrefix ${BamDir}/${SeqID}.mouse. \
        --sjdbGTFfile ${MouseGtfFile} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --genomeLoad NoSharedMemory \
        --runThreadN ${Threads} \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NM NH HI AS nM XS \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outSAMunmapped Within \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimScoreJunctionNonGTAG 0 \
        --chimMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --outBAMcompression 0 \
        --chimOutType WithinBAM HardClip \
        --chimScoreDropMax 30 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3       
    # sort and index
    ${RunSamtools} sort -@ ${Threads} -o ${BamDir}/${SeqID}.mouse.Aligned.out.Sorted.bam ${BamDir}/${SeqID}.mouse.Aligned.out.bam
    ${RunSamtools} index  -@ ${Threads} -b ${BamDir}/${SeqID}.mouse.Aligned.out.Sorted.bam 
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | STAR Alignment TwoPass Mode with MOUSE-REFERENCE FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run XenofilteR START." >> ${LogFile}
#==================================================================================================# 

#---| Run XenofilteR |-----------------------------------------------------------------------------#
    Rscript /storage/home/kangsm/runScripts/WTS/WTS_P05-2_RUN-XENOFILTER.R \
        --baseDir ${BaseDir} --batchID ${BatchID} --seqID ${SeqID} --threads ${Threads} 
#--------------------------------------------------------------------------------------------------#

#---| Post-Alignment Processing |------------------------------------------------------------------#
    cp ${BamDir}/${SeqID}.human.Log.final.out ${QcDir}/
    cp ${BamDir}/${SeqID}.mouse.Log.final.out ${QcDir}/
    cp ${BamDir}/${SeqID}.human.SJ.out.tab ${QcDir}/
    cp ${BamDir}/${SeqID}.mouse.SJ.out.tab ${QcDir}/
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Post-Processing FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqToBAM FINISHED." >> ${LogFile}
    echo " " >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqToBAM DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 





