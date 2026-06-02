
#! /bin/bash

# SCRIPT VERSION : v0.1
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| Input Options |------------------------------------------------------------------------------#
    usage() {
        echo "Usage: WES_run.NGSCheckMate_Rev.sh [ ARGUMENTS ]...
        --baseDir         : work base folder
        --seqFolder       : data processing base directory
        --sampleListFile  : sample List (BAM file list)
        --assembly        : genome assembly version. default = 'hg19'
        --platform        : WES or TSO. default = 'WES' 
        -h | --help       : Print this usage 
        "
    }
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR=/data/wes
    GENOME_ASSEMBLY=hg19
    PLATFORM="hg19"
#--------------------------------------------------------------------------------------------------#

#---| Parse Option Arguments |---------------------------------------------------------------------#
    ARGS=$(getopt -a -o h: --long seqFolder:,baseDir:,assembly:,sampleListFile:,platform:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --seqFolder )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SEQ_FOLDER=$2 ; shift 2 ;;
            esac ;;
        --platform )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = WES" ; exit 2 ;;
                *)  PLATFORM=$2 ; shift 2 ;;
            esac  ;;
        --sampleListFile )
            case "$2" in
                -* | --* | "") echo "no value at $1" >&2 ; exit 2 ;;
                *)  SAMPLE_LIST_FILE=$2 ; shift 2 ;;
            esac  ;;
        --baseDir )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = /data/wes" ; shift ;;
                *)  BASE_DIR=$2 ; shift 2 ;;
            esac  ;;
        --assembly )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = hg19" ; shift ;;
                *)  GENOME_ASSEMBLY=$2 ; shift 2 ;;
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
    IND_MATCH_RES_DIR=${DATA_DIR}/meta/individual_matching
    IND_MATCH_VCF_DATA_DIR=${IND_MATCH_RES_DIR}/vcf
    LOG_DIR=${DATA_DIR}/${SEQ_ID}/meta/analysis_log
    ##
    if [ ! -d ${IND_MATCH_RES_DIR} ]; then mkdir -p ${IND_MATCH_RES_DIR}; fi
    if [ ! -d ${IND_MATCH_VCF_DATA_DIR} ]; then mkdir -p ${IND_MATCH_VCF_DATA_DIR}; fi
#--------------------------------------------------------------------------------------------------#

#---| LOGFILE |------------------------------------------------------------------------------------#
    LOGFILE=${LOG_DIR}/run.individual.matching.$(date '+%Y%m%d').log
    if [ ! -f ${LOGFILE} ]; then touch ${LOGFILE}; fi
#--------------------------------------------------------------------------------------------------#

#---| SWTOOLS RUN |--------------------------------------------------------------------------------#   
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_SAMTOOLS="${RUN_SINGULARITY}/ngscheckmate-1.0.1.sif /usr/local/bin/samtools"
    RUN_NGSCHECKMATE="${RUN_SINGULARITY}/ngscheckmate-1.0.1.sif python2 /opt/NGSCheckMate/ncm.py"
#--------------------------------------------------------------------------------------------------#

#---| REFERENCES |---------------------------------------------------------------------------------#
    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        export GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
        export BED=/opt/NGSCheckMate/SNP/SNP_GRCh38_hg38_wChr.bed
    else
        if [ "${PLATFORM}" = "TSO" ]; then
            export GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
            export BED=/opt/NGSCheckMate/SNP/SNP_GRCh37_hg19_wChr.bed
        else
            export GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
            export BED=//opt/NGSCheckMate/SNP/SNP_GRCh37_hg19_woChr.bed
        fi
    fi
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo " " >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo " INDIVIDUAL-MATCHING "  >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | SNP VCF Generation START." >> ${LOGFILE}
#==================================================================================================# 

#---| GENERATE INDIV-MATCH-SNP-VCF |---------------------------------------------------------------#
    cat $SAMPLE_LIST_FILE | awk '{print $1}' > ${IND_MATCH_RES_DIR}/indiv_match.list

    # for SEQ_ID in $(cat ${IND_MATCH_RES_DIR}/indiv_match.list)
    # do
    #     echo -e " >> Generate Individual-Matching-SNP VCF: ${SEQ_ID} ...."
    #     INPUT_BAM_FILE=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam/${SEQ_ID}.analysisReady.bam
    #     ${RUN_SAMTOOLS} mpileup -I -uf ${GENOME_FASTA} -l ${BED} ${INPUT_BAM_FILE} | bcftools call -c - > ${IND_MATCH_VCF_DATA_DIR}/${SEQ_ID}.indiv.match.vcf
    #     echo -e " >> Generate Individual-Matching-SNP VCF: ${SEQ_ID} ....DONE!"
    # done

    parallel -j 15 -k "${RUN_SAMTOOLS} mpileup -I -uf ${GENOME_FASTA} -l ${BED} ${BASE_DIR}/${SEQ_FOLDER}/{}/bam/{}.analysisReady.bam | bcftools call -c - > ${IND_MATCH_VCF_DATA_DIR}/{}.indiv.match.vcf" ::: $(cat ${IND_MATCH_RES_DIR}/indiv_match.list)

#--------------------------------------------------------------------------------------------------#    

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | SNP VCF Generation FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | INDIVIDUAL-MATCHING ANALYSIS START." >> ${LOGFILE}
#==================================================================================================# 

#---| RUN NGSCHECKMATE |---------------------------------------------------------------------------#
    touch ${IND_MATCH_RES_DIR}/indiv_match_analysis.list
    for SEQ_ID in $(cat ${IND_MATCH_RES_DIR}/indiv_match.list)
    do
        echo "${IND_MATCH_RES_DIR}/vcf/${SEQ_ID}.indiv.match.vcf" >> ${IND_MATCH_RES_DIR}/indiv_match_analysis.list
    done
    ${RUN_NGSCHECKMATE} -V -l ${IND_MATCH_RES_DIR}/indiv_match_analysis.list -bed ${BED} -O ${IND_MATCH_RES_DIR} -N ${SEQ_FOLDER}
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    echo "$(date '+%Y-%m-%d %H:%M:%S') | INDIVIDUAL-MATCHING ANALYSIS FINISHED." >> ${LOGFILE}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | Run NGS-CheckMate FINISHED." >> ${LOGFILE}
    echo "----------------------------------------------------------------------" >> ${LOGFILE}
#==================================================================================================# 



