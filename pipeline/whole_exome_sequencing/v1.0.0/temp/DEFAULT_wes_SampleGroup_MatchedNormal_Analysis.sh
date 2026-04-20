
#! /bin/bash

#$ -N JOB_NAME_HERE
#$ -S /bin/bash
#$ -pe mpi THREADS_N_HERE
#$ -t 1-TOTAL_TASK_HERE
#$ -tc JOB_N_HERE
#$ -e ERROR_LOG_HERE
#$ -o OUT_LOG_HERE

# >>>>> USER IMPUT START
#--------------------------------------------------------------------------------------------------#
BASE_DIR="BASE_DIR_HERE"
SEQ_FOLDER="BATCH_ID_HERE"
INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/matched_analysis.list
# NORMAL_SEQ_ID TUMOR_SEQ_ID SEQ_DEPTH NORMAL_TYPE #
#--------------------------------------------------------------------------------------------------#
THREADS=THREADS_N_HERE
GENOME_ASSEMBLY="hg19"
WES_LIB_KIT=NGS_LIB_HERE
MAF_DP=10
MAF_ALT=2
MAF_AF=0.01
# VCF_TAG="mutect2.NT"
#--------------------------------------------------------------------------------------------------#
RUN_MUTECT2_NT_GENERAL="ON"
RUN_MUTECT2_NT_KEEP_GERMLINE="ON"
RUN_MUTECT2_NT_ANNOTATION="ON"
RUN_MUTECT2_NT_MAF="ON"
RUN_TMB_NT="ON"
RUN_MSI_NT="ON"
RUN_CNV_ASCAT="ON"
RUN_CNV_PURECN_NT="ON"
RUN_CNV_CNVKIT_NT="OFF"
RUN_CNV_SUPERFREQ_NT="ON"
#--------------------------------------------------------------------------------------------------#
# <<<<< USER INPUT END
# NORMAL_SEQ_ID TUMOR_SEQ_ID DEPTH NORMAL_TYPE
#--------------------------------------------------------------------------------------------------#
SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
TUMOR_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $2}')
NORMAL_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
SEQ_DEPTH=$(echo ${SEQ_INFO} | awk '{print $3}')
NORMAL_TYPE=$(echo ${SEQ_INFO} | awk '{print $4}')
#--------------------------------------------------------------------------------------------------#
RUN_SCRIPT=/storage/home/kangsm/runScripts
BED_PREFIX=/storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
#--------------------------------------------------------------------------------------------------#

#---| PIPELINE RUN LOG |---------------------------------------------------------------------------#
    PIPELINE_LOG_DIR=${BASE_DIR}/${SEQ_FOLDER}/meta/analysis_log
    if [ ! -d ${PIPELINE_LOG_DIR} ]; then mkdir -p ${PIPELINE_LOG_DIR}; fi
    PIPELINE_LOG=${PIPELINE_LOG_DIR}/${TUMOR_SEQ_ID}.matched.analysis.$(date '+%Y-%m-%d').log
    if [ ! -f ${PIPELINE_LOG} ]; then touch ${PIPELINE_LOG}; fi
#--------------------------------------------------------------------------------------------------#

#---| BAM FILE CHECK |-----------------------------------------------------------------------------#
    TUMOR_BAM=${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/bam/${TUMOR_SEQ_ID}.analysisReady.bam
    NORMAL_BAM=${BASE_DIR}/${SEQ_FOLDER}/${NORMAL_SEQ_ID}/bam/${NORMAL_SEQ_ID}.analysisReady.bam
    if [ ! -f ${TUMOR_BAM} ]; then 
        echo "No TUMOR BAM FILE. Stop."
        echo "$(date '+%Y-%m-%d %H:%M:%S') | NO TUMOR BAM FILE. PIPELINE RUNNING STOPPED." >> ${PIPELINE_LOG}
        exit 0 
    fi
    if [ ! -f ${NORMAL_BAM} ]; then 
        echo "No matched NORMAL BAM FILE. Stop."
        echo "$(date '+%Y-%m-%d %H:%M:%S') | NO MATCHED NORMAL BAM FILE. PIPELINE RUNNING STOPPED." >> ${PIPELINE_LOG}
        exit 0 
    fi
#--------------------------------------------------------------------------------------------------#

#---| 1. MUTECT2 NT-MODE VARIANT CALLING |---------------------------------------------------------#
    if [ "${RUN_MUTECT2_NT_GENERAL}" = "ON" ]; then
        if [ "${RUN_MUTECT2_NT_KEEP_GERMLINE}" = "ON" ]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : KEEP GERMLINE ] START." >> ${PIPELINE_LOG}
            #------------------------------------------------------------------#
            ${RUN_SCRIPT}/WES_run.Mutect2_MatchedNormalMode.sh \
                --seqFolder ${SEQ_FOLDER} \
                --baseDir ${BASE_DIR} \
                --threads ${THREADS} \
                --tumorSeqID ${TUMOR_SEQ_ID} \
                --normalSeqID ${NORMAL_SEQ_ID} \
                --interval ${BED_PREFIX}.target.interval_list \
                --assembly ${GENOME_ASSEMBLY} \
                --seqDepth ${SEQ_DEPTH} \
                --keepGermline "true" \
                --nType ${NORMAL_TYPE}
            #------------------------------------------------------------------#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : KEEP GERMLINE ] FINISHED." >> ${PIPELINE_LOG}
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : GENERAL ] START." >> ${PIPELINE_LOG}
            #------------------------------------------------------------------#
            ${RUN_SCRIPT}/WES_run.Mutect2_MatchedNormalMode.sh \
                --seqFolder ${SEQ_FOLDER} \
                --baseDir ${BASE_DIR} \
                --threads ${THREADS} \
                --tumorSeqID ${TUMOR_SEQ_ID} \
                --normalSeqID ${NORMAL_SEQ_ID} \
                --interval ${BED_PREFIX}.target.interval_list \
                --assembly ${GENOME_ASSEMBLY} \
                --seqDepth ${SEQ_DEPTH} \
                --keepGermline "false" \
                --nType ${NORMAL_TYPE}
            #------------------------------------------------------------------#
            echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : GENERAL ] FINISHED." >> ${PIPELINE_LOG}
        fi        
    fi
#--------------------------------------------------------------------------------------------------#
#---| 2. MUTECT2 NT-MODE VARIANT ANNOTATION : GENERAL ONLY |---------------------------------------#
    if [ "${RUN_MUTECT2_NT_ANNOTATION}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT ANNOTATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.VariantAnnotation.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --seqID ${TUMOR_SEQ_ID} \
            --nType ${NORMAL_TYPE} \
            --assembly ${GENOME_ASSEMBLY} \
            --variantType "somatic" \
            --variantCallMode "nt"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT ANNOTATION ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 3. MUTECT2 NT-MODE VARIANT MAF-PROCESSING |--------------------------------------------------#
    if [ "${RUN_MUTECT2_NT_MAF}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.MAF_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${TUMOR_SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --nType ${NORMAL_TYPE} \
            --writeSomaticVarinatsOnlyMAF TRUE \
            --DP ${MAF_DP} \
            --ALT ${MAF_ALT} \
            --POPAF ${MAF_AF} \
            --writeDB TRUE \
            --variantCallMode "nt" \
            --GERMLINE FALSE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 4. TMB CALCULATION : NT-MODE |---------------------------------------------------------------#
    if [ "${RUN_TMB_NT}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.TMB_Calculation.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${TUMOR_SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --nType ${NORMAL_TYPE} \
            --useReferenceMethod TRUE \
            --useReferenceSize TRUE \
            --variantCallMode "nt" \
            --importDB TRUE
        #----------------------------------------------------------------------#
        TMB_RES_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/tmb

        if [ "${NORMAL_TYPE}" = "ORG" ]; then 
            RES_TAG="ORG.NT"
        else
            RES_TAG="NT"
        fi

        mv ${TMB_RES_DIR}/${TUMOR_SEQ_ID}.TMB.calculation.result.txt ${TMB_RES_DIR}/${TUMOR_SEQ_ID}.${RES_TAG}.TMB.calculation.result.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] FINISHED." >> ${PIPELINE_LOG}    
    fi
#--------------------------------------------------------------------------------------------------#
#---| 5. MSI CALCULATION : NT-MODE |---------------------------------------------------------------#
    if [ "${RUN_MSI_NT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSI CALCULATION NT-MODE ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.MSI_Calculation.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --runMode "nt" \
            --targetBed ${BED_PREFIX}.target.bed  \
            --assembly ${GENOME_ASSEMBLY}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.MSI_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${TUMOR_SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --runMode "nt" \
            --nType ${NORMAL_TYPE} \
            --importDB TRUE
        #----------------------------------------------------------------------#
        MSI_DIR=${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/msi

        if [ "${NORMAL_TYPE}" = "ORG" ]; then 
            RES_TAG="ORG.NT"
        else
            RES_TAG="NT"
        fi

        mv ${MSI_DIR}/${TUMOR_SEQ_ID}.MSI.status.calculation.results.txt ${MSI_DIR}/${TUMOR_SEQ_ID}.${RES_TAG}.MSI.status.calculation.results.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSI CALCULATION NT-MODE ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 6. ASCAT : CNV(1), PLOIDY, PURITY, CLONALITY(1) |--------------------------------------------#
    if [ "${RUN_CNV_ASCAT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ASCAT : CNV(1), PLOIDY, PURITY, CLONALITY(1) ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript /storage/home/kangsm/runScripts/WES_run.CNV_ASCAT.R \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID}  \
            --wesLibKit ${WES_LIB_KIT} \
            --assembly ${GENOME_ASSEMBLY} \
            --threads ${THREADS} \
            --nType ${NORMAL_TYPE}
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ASCAT : CNV(1), PLOIDY, PURITY, CLONALITY(1) ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 7. PureCN : PLOIDY & PURITY |----------------------------------------------------------------#
    if [ "${RUN_CNV_PURECN_NT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PureCN : PLOIDY & PURITY ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_PureCN.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --runMode "nt" \
            --assembly ${GENOME_ASSEMBLY} \
            --wesLibKit ${WES_LIB_KIT} \
            --callMinDepth 20 \
            --nType ${NORMAL_TYPE} \
            --forceCoverage "false"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PureCN : PLOIDY & PURITY ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 8. SUPERFREQ : CNV(2) & CLONALITY(2) |-------------------------------------------------------#
    if [ "${RUN_CNV_SUPERFREQ_NT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SUPERFREQ : CNV(2) & CLONALITY(2) ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.CNV_superFreq.R \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --assembly ${GENOME_ASSEMBLY} \
            --wesLibKit ${WES_LIB_KIT} \
            --mode "exome" \
            --threads ${THREADS} \
            --nType ${NORMAL_TYPE} \
            --setID ${TUMOR_SEQ_ID}.NormalTumor
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SUPERFREQ : CNV(2) & CLONALITY(2) ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 9. CNVKIT : CNV(3) |-------------------------------------------------------------------------#
    if [ "${RUN_CNV_CNVKIT_NT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNVKIT : CNV(3) ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_CNVkit.sh \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --runMode "nt" \
            --wesLibKit ${WES_LIB_KIT} \
            --assembly ${GENOME_ASSEMBLY} \
            --threads ${THREADS} \
            --nType ${NORMAL_TYPE} \
            --superfreqSetID NormalTumor
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNVKIT : CNV(3) ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#


