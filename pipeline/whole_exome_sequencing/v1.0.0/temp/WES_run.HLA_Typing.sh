
#!/bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.HLA_Typing.sh  [ ARGUMENTS ]...
        [ -d | --seqFolder ]      : data processing base directory
        [ -s | --seqID ]          : sample ID (Used as sample's processing product folder name)
        [ -t | --inputType ]      : input NGS type. DNA or RNA. (default = DNA)
        [ --baseDir ]             : Base Work Dir
        [ --threads ]             : N-CPU threads
        [ -h | --help ]           : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    INPUT_TYPE="DNA"
    IMPORT_DB="false"
    BASE_DIR=/data/wes
    THREADS=15
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o d:s:l:t:h: --long seqFolder:,seqID:,inputType:,baseDir:,threads:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        -d | --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        -s | --seqID )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_ID=$2 ; shift 2 ;;
            esac  ;;
        -t | --inputType )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = DNA" ; shift ;;
                *)  INPUT_TYPE=$2 ; shift 2 ;;
            esac ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes" ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac ;;
        --threads )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = 15" ; shift ;;
                *)  THREADS=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR=${BASE_DIR}/${SEQ_FOLDER}
    HLA_RES_DIR=${DATA_DIR}/${SEQ_ID}/hla
    if [ ! -d ${HLA_RES_DIR} ]; then mkdir -p ${HLA_RES_DIR}; fi
    ##
    OPTITYPE_RES_DIR=${HLA_RES_DIR}/optitype
    if [ ! -d ${OPTITYPE_RES_DIR} ]; then mkdir -p ${OPTITYPE_RES_DIR}; fi
    ##
    HLA_LA_RES_DIR=${HLA_RES_DIR}/hla-la
    if [ ! -d ${HLA_LA_RES_DIR} ]; then mkdir -p ${HLA_LA_RES_DIR}; fi
    ##
    FASTQ_DIR=${DATA_DIR}/${SEQ_ID}/fastq
    BAM_DIR=${DATA_DIR}/${SEQ_ID}/bam  
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/log
    if [ ! -d ${LOG_DIR} ]; then mkdir -p ${LOG_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/${SEQ_ID}.HLA.typing.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_CONDA_HLA="${RUN_SINGULARITY}/conda_HLA.sif"
    RUN_HLA_LA="/storage/apps/HLA-LA-1.0.3/src/HLA-LA.pl"
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${INPUT_TYPE}" = "DNA" ]; then
        OPTITYPE_REF_FASTA=/storage/references_and_index/optitype_references/hla_reference_dna.fasta
    else
        OPTITYPE_REF_FASTA=/storage/references_and_index/optitype_references/hla_reference_rna.fasta
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "WES HLA-Typing log" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " RUN DATE             : $(date '+%Y-%m-%d %H:%M:%S')" >> ${LOGFILE}
    echo " SEQ FOLDER           : ${DATA_DIR}" >> ${LOGFILE}
    echo " SEQ ID               : ${SEQ_ID}" >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}   
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | HLA-TYPING START." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run OptiType START." >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - Fished-BAM generation Start." >> ${LOGFILE}
#==================================================================================================# 

#---| OPTITYPE : FISHED-BAM |----------------------------------------------------------------------#
    READS="R1 R2"
    ##
    parallel -j 2 -k "${RUN_CONDA_HLA} razers3 \
        --percent-identity 95 \
        --max-hits 1 \
        --distance-range 0 \
        --output ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_fished_{}.bam \
        ${OPTITYPE_REF_FASTA} \
        ${FASTQ_DIR}/${SEQ_ID}.trimmed_{}.fastq.gz" ::: ${READS}
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - Fished-BAM generation Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - Revert-FASTQ generation Start." >> ${LOGFILE}
#==================================================================================================# 

#---| OPTITYPE : REVERT-FASTQ |----------------------------------------------------------------------#
    parallel -j 2 -k "${RUN_CONDA_HLA} samtools bam2fq \
        ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_fished_{}.bam > ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_{}.fastq" ::: ${READS}
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - Revert-FASTQ generation Finished." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') |   - OptiType Start." >> ${LOGFILE}
#==================================================================================================# 

#---| OPTITYPE : HLA-TYPING |----------------------------------------------------------------------#
    if [ "${INPUT_TYPE}" = "DNA" ]; then
        ${RUN_CONDA_HLA} OptiTypePipeline.py \
            --input ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R1.fastq ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R2.fastq \
            --dna \
            --outdir ${OPTITYPE_RES_DIR} \
            --prefix ${SEQ_ID}.optitype
    else
        ${RUN_CONDA_HLA} OptiTypePipeline.py \
            --input ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R1.fastq ${OPTITYPE_RES_DIR}/${SEQ_ID}.optitype_R2.fastq \
            --rna \
            --outdir ${OPTITYPE_RES_DIR} \
            --prefix ${SEQ_ID}.optitype
    fi 
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run OptiType FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run HLA-LA START." >> ${LOGFILE}
#==================================================================================================# 

#---| HLA-LA : HLA-TYPING |------------------------------------------------------------------------#
    ${RUN_HLA_LA} \
        --BAM ${BAM_DIR}/${SEQ_ID}.analysisReady.bam \
        --graph PRG_MHC_GRCh38_withIMGT \
        --sampleID ${SEQ_ID} \
        --maxThreads ${THREADS} \
        --workingDir ${HLA_LA_RES_DIR} 
#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run HLA-LA FINISHED." >> ${LOGFILE}
#==================================================================================================# 

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | HLA-TYPING DONE." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " " >> ${LOGFILE}
#==================================================================================================# 

