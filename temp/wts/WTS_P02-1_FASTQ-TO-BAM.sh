
#! /bin/bash

# version : release-1.0

#---| SCRIPT USAGE |-------------------------------------------------------------------------------#
    usage() { echo "  Usage: WTS_run.FastqToBAM.sh [ARGUMENTS]...

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

        --genomeAssembly    Reference genome assembly version.
                            'hg19' and 'hg38' enable.
                            default = hg19

        --paired            Sequencing type is paired-end sequencing or not.
                            default = true (paired-end)
        
        -h | --help         Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| DEFAULT VALUES |-----------------------------------------------------------------------------#
    # BaseDir="/data/wts"
    # Threads=15
    # ReadLength=150
    # GenomeAssembly="hg19"
    # Paired="true"
#--------------------------------------------------------------------------------------------------#

#---| OPTION PARSE |-------------------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long batchID:,seqID:,threads:,baseDir:,readLength:,genomeAssembly:,paired:,help -- "$@" )
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
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
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
    if [ -z ${GenomeAssembly} ]; then GenomeAssembly=hg19; fi
    if [ -z ${Paired} ]; then Paired=true; fi
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
    StarIndex=/storage/references_and_index/${GenomeAssembly}/star-rsem-index/PE-${ReadLength}
    GtfFile=/storage/references_and_index/${GenomeAssembly}/gtf/${GenomeAssembly}.gtf
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#
    RunSingularity="singularity exec -B /storage,/data /storage/images"
    RunStar="${RunSingularity}/star-2.7.11.sif STAR"
    RunSamtools="${RunSingularity}/samtools-1.18.sif samtools"
#--------------------------------------------------------------------------------------------------#

#---| LogFile |------------------------------------------------------------------------------------#
    LogFile=${LogDir}/${SeqID}.fastq.to.bam.$(date '+%Y%m%d').log
    if [ ! -f ${LogFile} ]; then touch ${LogFile}; fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LogFile}
    echo "WTS FastqToBAM: STAR Alignment log" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LogFile}
    echo " BATCH ID             : ${DataDir}" >> ${LogFile}
    echo " SEQ ID               : ${SeqID}" >> ${LogFile}
    echo " READ LENGTH          : ${ReadLength}" >> ${LogFile}
    echo " STAR MODE            : Two-Pass Mode" >> ${LogFile}
    echo " OPTION COMPARTIBLE   : Arriba" >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqToBAM: STAR Alignment START." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | STAR Alignment TwoPass Mode with Arriba Compartoble Options START." >> ${LogFile}
#==================================================================================================# 

#---| RUN STAR ALIGNMENT |-------------------------------------------------------------------------#
    if [ "${ReadLength}" = "150" ]; then 
        MaxMismatch=15
    else
        MaxMismatch=10
    fi
    #--------------------------------------------------------------------------#
    if [ "${Paired}" = "true" ]; then
        ${RunStar} \
            --genomeDir ${StarIndex} \
            --readFilesIn ${FastqDir}/${SeqID}.trimmed_R1.fastq.gz ${FastqDir}/${SeqID}.trimmed_R2.fastq.gz \
            --outFileNamePrefix ${BamDir}/${SeqID}. \
            --sjdbGTFfile ${GtfFile} \
            --readFilesCommand zcat \
            --quantMode TranscriptomeSAM \
            --twopassMode Basic \
            --genomeLoad NoSharedMemory \
            --runThreadN ${Threads} \
            --outFilterMatchNminOverLread 0.33 \
            --outFilterScoreMinOverLread 0.33 \
            --readQualityScoreBase 20 \
            --outFilterMismatchNmax ${MaxMismatch} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outSAMattributes NH HI AS nM XS NM \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimScoreJunctionNonGTAG 0 \
            --chimMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --outBAMcompression 0 \
            --outFilterMultimapNmax 50 \
            --chimOutType WithinBAM HardClip \
            --chimScoreDropMax 30 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3       
    else
        ${RunStar} \
            --genomeDir ${StarIndex} \
            --readFilesIn ${FastqDir}/${SeqID}.trimmed.fastq.gz \
            --outFileNamePrefix ${BamDir}/${SeqID}. \
            --sjdbGTFfile ${GtfFile} \
            --readFilesCommand zcat \
            --quantMode TranscriptomeSAM \
            --twopassMode Basic \
            --genomeLoad NoSharedMemory \
            --runThreadN ${Threads} \
            --outFilterMatchNminOverLread 0.33 \
            --outFilterScoreMinOverLread 0.33 \
            --readQualityScoreBase 20 \
            --outFilterMismatchNmax ${MaxMismatch} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outSAMattributes NH HI AS nM XS NM \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimScoreJunctionNonGTAG 0 \
            --chimMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --outBAMcompression 0 \
            --outFilterMultimapNmax 50 \
            --chimOutType WithinBAM HardClip \
            --chimScoreDropMax 30 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3       
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | STAR Alignment TwoPass Mode with Arriba Compartoble Options FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Post-Processing START." >> ${LogFile}
#==================================================================================================# 

#---| Post-Alignment Processing |------------------------------------------------------------------#
    # SORTING ------------------------------------------------------------------
    ${RunSamtools} sort -@ ${Threads} -o ${BamDir}/${SeqID}.Aligned.out.Sorted.bam ${BamDir}/${SeqID}.Aligned.out.bam
    # CREATE INDEX -------------------------------------------------------------
    ${RunSamtools} index  -@ ${Threads} -b ${BamDir}/${SeqID}.Aligned.out.Sorted.bam 
    # QC FILES -----------------------------------------------------------------
    cp ${BamDir}/${SeqID}.Log.final.out ${QcDir}/
    cp ${BamDir}/${SeqID}.SJ.out.tab ${QcDir}/
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Post-Processing FINISHED." >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqToBAM FINISHED." >> ${LogFile}
    echo " " >> ${LogFile}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | FastqToBAM DONE." >> ${LogFile}
    echo "----------------------------------------------------------------------" >> ${LogFile}
    echo " " >> ${LogFile}
#==================================================================================================# 





