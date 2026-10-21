
#! /bin/bash

#$ -N hg19_cnv_paired
#$ -S /bin/bash
#$ -pe mpi 15
#$ -t 1-6
#$ -tc 3
#$ -e /storage/home/kangsm/sge_log/cnv_seqc2/paired.analysis.hg19.240325.02.err
#$ -o /storage/home/kangsm/sge_log/cnv_seqc2/paired.analysis.hg19.240325.02.out

#--------------------------------------------------------------------------------------------------#
BASE_DIR=/storage/home/kangsm/analysis/cnv_seqc2
SEQ_FOLDER=hg19
INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/NT_match.list
#--------------------------------------------------------------------------------------------------#
SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
TUMOR_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $2}')
NORMAL_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
SEQ_DEPTH=$(echo ${SEQ_INFO} | awk '{print $3}')
#--------------------------------------------------------------------------------------------------#
THREADS=15
GENOME_ASSEMBLY="hg19"
WES_LIB_KIT="sureselect.v6"
#----------------------------------------------------les----------------------------------------------#
RUN_SCRIPT=/storage/home/kangsm/runScripts
BED_PREFIX=/storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
MAF_DP=25
MAF_ALT=3
MAF_AF=0.01
#--------------------------------------------------------------------------------------------------#
RUN_MUTECT2_NT_GENERAL="OFF"
RUN_MUTECT2_NT_ANNOTATION="OFF"
RUN_MUTECT2_NT_MAF="OFF"
RUN_MUTECT2_NT_PURECN="OFF"
RUN_TMB_NT="OFF"
RUN_MSI_NT="OFF"
RUN_CNV_GATK="OFF"
RUN_CNV_PURECN="OFF"
RUN_CNV_CNVKIT="ON"
#--------------------------------------------------------------------------------------------------#
#---| PIPELINE RUN LOG |---------------------------------------------------------------------------#
    PIPELINE_LOG_DIR=${BASE_DIR}/${SEQ_FOLDER}/meta//analysis_run_log
    if [ ! -d ${PIPELINE_LOG_DIR} ]; then mkdir -p ${PIPELINE_LOG_DIR}; fi
    PIPELINE_LOG=${PIPELINE_LOG_DIR}/${TUMOR_SEQ_ID}.paired.analysis.hg19.$(date '+%Y-%m-%d').log
    if [ ! -f ${PIPELINE_LOG} ]; then touch ${PIPELINE_LOG}; fi
#--------------------------------------------------------------------------------------------------#

#---| 1. MUTECT2 NT-MODE VARIANT CALLING : GENERAL |-----------------------------------------------#
    if [ "${RUN_MUTECT2_NT_GENERAL}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : GENERAL ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.Mutect2_MatchedNormalMode.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --interval ${BED_PREFIX}.target.interval_list \
            --assembly ${GENOME_ASSEMBLY} \
            --seqDepth ${SEQ_DEPTH} \
            --purecn "false"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : GENERAL ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 2. MUTECT2 NT-MODE VARIANT ANNOTATION : GENERAL |--------------------------------------------#
    if [ "${RUN_MUTECT2_NT_ANNOTATION}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT ANNOTATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.VariantAnnotation.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --seqID ${TUMOR_SEQ_ID} \
            --varCallMode "nt" \
            --assembly ${GENOME_ASSEMBLY}
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
            --varCallMode "nt" \
            --writeSomaticVarinatsOnlyMAF TRUE \
            --DP ${MAF_DP} \
            --ALT ${MAF_ALT} \
            --POPAF ${MAF_AF} \
            --writeDB FALSE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 4. MUTECT2 NT-MODE VARIANT CALLING : PURECN |------------------------------------------------#
    if [ "${RUN_MUTECT2_NT_PURECN}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : PURECN ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.Mutect2_MatchedNormalMode.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --interval ${BED_PREFIX}.target.interval_list \
            --assembly ${GENOME_ASSEMBLY} \
            --seqDepth ${SEQ_DEPTH} \
            --purecn "true"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MUTECT2 NT-MODE VARIANT CALLING : PURECN ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 5. TMB CALCULATION : NT-MODE |---------------------------------------------------------------#
    if [ "${RUN_TMB_NT}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.TMB_Calculation.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${TUMOR_SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --maf ${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/vcf/${TUMOR_SEQ_ID}.mutect2.NT.filtered.somatic.variants.only.maf \
            --useReference TRUE \
            --useReferenceSize TRUE \
            --variantCallMode "nt"
        #----------------------------------------------------------------------#
        mv ${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/tmb/${TUMOR_SEQ_ID}.TMB.calculation.result.txt ${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/tmb/${TUMOR_SEQ_ID}.NT.TMB.calculation.result.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] FINISHED." >> ${PIPELINE_LOG}    
    fi
#--------------------------------------------------------------------------------------------------#
#---| 6. MSI CALCULATION : NT-MODE |---------------------------------------------------------------#
    if [ "${RUN_MSI_NT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSI CALCULATION NT-MODE ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.MSI_Calculation.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --runMode 'nt' \
            --targetBed ${BED_PREFIX}.target.bed  \
            --assembly ${GENOME_ASSEMBLY}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.MSI_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${TUMOR_SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --runMode 'nt' \
            --importDB FALSE
        #----------------------------------------------------------------------#
        MSI_DIR=${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/msi
        #----------------------------------------------------------------------#
        mv ${MSI_DIR}/${TUMOR_SEQ_ID}.MSI.status.calculation.results.txt ${MSI_DIR}/${TUMOR_SEQ_ID}.NT.MSI.status.calculation.results.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSI CALCULATION NT-MODE ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 7. CNV : GATK |------------------------------------------------------------------------------#
    if [ "${RUN_CNV_GATK}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : GATK ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_GATK.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --runNTmode "true" \
            --genomeAssembly ${GENOME_ASSEMBLY} \
            --wesLibKit ${WES_LIB_KIT}
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : GATK ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 8. CNV : PURECN |----------------------------------------------------------------------------#
    if [ "${RUN_CNV_PURECN}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : PURECN ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_PureCN.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --runMode 'nt' \
            --assembly ${GENOME_ASSEMBLY} \
            --wesLibKit ${WES_LIB_KIT} \
            --callMinDepth 20
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : PURECN ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#

#---| 8. CNV : CNVKIT |----------------------------------------------------------------------------#
    if [ "${RUN_CNV_CNVKIT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : CNVKIT ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_CNVkit.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --tumorSeqID ${TUMOR_SEQ_ID} \
            --normalSeqID ${NORMAL_SEQ_ID} \
            --runMode 'nt' \
            --assembly ${GENOME_ASSEMBLY} \
            --wesLibKit ${WES_LIB_KIT} \
            --threads ${THREADS}
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : CNVKIT ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#



