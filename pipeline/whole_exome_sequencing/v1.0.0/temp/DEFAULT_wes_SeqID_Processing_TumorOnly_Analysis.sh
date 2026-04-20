
#! /bin/bash

#$ -N  JOB_NAME_HERE
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
INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/sample.list
# SEQ_ID SEQ_DEPTH FLOWCELL_ID READ_LENGTH #
#------------------------------------------------------------------------------#
THREADS=THREADS_N_HERE
GENOME_ASSEMBLY="hg19"
WES_LIB_KIT="NGS_LIB_HERE"
MAF_DP=10
MAF_ALT=2
MAF_AF=0.01
#--------------------------------------------------------------------------------------------------#
RUN_FASTQ_SCREEN="ON"
RUN_FASTQ_QC="ON"
RUN_FASTQ_TO_BAM="ON"
RUN_BAM_QC="ON"
RUN_MUTECT2_TONLY="ON"
RUN_MUTECT2_TONLY_ANNOTATION="ON"
RUN_MUTECT2_TONLY_MAF="ON"
RUN_QC_SUMMARY="ON"
#------------------------------------------------------------------------------#
RUN_TMB_TONLY="ON"
RUN_MSI_TONLY="ON"
RUN_HLA="ON"
RUN_DEEPVARIANT="GERMLINE_VARS_RUN_OPTION"
RUN_DEEPVARIANT_ANNOTATION="GERMLINE_VARS_RUN_OPTION"
RUN_DEEPVARIANT_MAF="GERMLINE_VARS_RUN_OPTION"
RUN_CNV_ASCAT="ASCAT_CNV_RUN_OPTION"
RUN_CNV_PURECN="PURECN_CNV_RUN_OPTION"
RUN_CNV_SUPERFREQ="SUPERFREQ_CNV_RUN_OPTION"
RUN_CNV_CNVKIT_TONLY="OFF"
#--------------------------------------------------------------------------------------------------#
# <<<<< USER INPUT END

#--------------------------------------------------------------------------------------------------#
SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
SEQ_DEPTH=$(echo ${SEQ_INFO} | awk '{print $2}')
FLOWCELL_ID=$(echo ${SEQ_INFO} | awk '{print $3}')
READ_LENGTH=$(echo ${SEQ_INFO} | awk '{print $4}')
#--------------------------------------------------------------------------------------------------#
RGLB=${WES_LIB_KIT}
RUN_SCRIPT=/storage/home/kangsm/runScripts
BED_PREFIX=/storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
#--------------------------------------------------------------------------------------------------#

#---| PIPELINE RUN LOG |---------------------------------------------------------------------------#
    PIPELINE_LOG_DIR=${BASE_DIR}/${SEQ_FOLDER}/meta/analysis_log
    if [ ! -d ${PIPELINE_LOG_DIR} ]; then mkdir -p ${PIPELINE_LOG_DIR}; fi
    #--------------------------------------------------------------------------#
    PIPELINE_LOG=${PIPELINE_LOG_DIR}/${SEQ_ID}.sample.level.analysis.$(date '+%Y-%m-%d').log
    if [ ! -f ${PIPELINE_LOG} ]; then touch ${PIPELINE_LOG}; fi
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
    echo "   " >> ${PIPELINE_LOG}
    echo "--------------------------------------------------------------------------" >> ${PIPELINE_LOG}
    echo "$(date '+%Y-%m-%d %H:%M:%S') | SINGLE SAMPLE LEVEL ANALYSIS PIPELINE" >> ${PIPELINE_LOG}
#--------------------------------------------------------------------------------------------------#

#---| 1. FASTQ-SCREEN |----------------------------------------------------------------------------#
    if [ "${RUN_FASTQ_SCREEN}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FASTQ-SCREEN ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.FastqScreen.sh \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --threads ${THREADS}
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FASTQ-SCREEN ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 2. FASTQ-QC |--------------------------------------------------------------------------------#
    if [ "${RUN_FASTQ_QC}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FASTQ-QC ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.FastqQC.sh \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --threads ${THREADS} \
            --readLength ${READ_LENGTH} \
            --runMultiqc "false"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FASTQ-QC ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 3. FASTQ-TO-BAM |----------------------------------------------------------------------------#
    if [ "${RUN_FASTQ_TO_BAM}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FASTQ-TO-BAM ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.FastqToBAM.sh \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --threads ${THREADS} \
            --targetBed ${BED_PREFIX}.target.bed \
            --assembly ${GENOME_ASSEMBLY} \
            --RGID ${FLOWCELL_ID} \
            --RGLB ${RGLB} \
            --RGPL "ILLUMINA" \
            --RGCN "GCX"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ FASTQ-TO-BAM ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 4. BAM & ALIGNMENT QC |----------------------------------------------------------------------#
    if [ "${RUN_BAM_QC}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ BAM & ALIGNMENT QC ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.AnalysisReadyBAM_QC.sh \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --assembly ${GENOME_ASSEMBLY} \
            --targetBed ${BED_PREFIX}.target.bed \
            --intervalBait ${BED_PREFIX}.probe.interval_list \
            --intervalTarget ${BED_PREFIX}.target.interval_list \
            --removeTempBAM "true" \
            --runMultiqc "true"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ BAM & ALIGNMENT QC ] FINISHED." >> ${PIPELINE_LOG}    
    fi
#--------------------------------------------------------------------------------------------------#
#---| 5. VARIANT CALLING : MUTECT2 TUMOR ONLY |----------------------------------------------------#
    if [ "${RUN_MUTECT2_TONLY}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VARIANT CALLING : MUTECT2 TUMOR ONLY ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.Mutect2_TumorOnlyMode.sh \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --interval ${BED_PREFIX}.target.interval_list \
            --assembly ${GENOME_ASSEMBLY} \
            --seqDepth ${SEQ_DEPTH} \
            --normalTonlyCall "true" \
            --keepGermlineCall "true"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VARIANT CALLING : MUTECT2 TUMOR ONLY ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 6. VARIANT ANNOTATION : MUTECT2 TUMOR ONLY |-------------------------------------------------#
    if [ "${RUN_MUTECT2_TONLY_ANNOTATION}" = "ON" ]; then      
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VARIANT ANNOTATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.VariantAnnotation.sh \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --nType "TS" \
            --assembly ${GENOME_ASSEMBLY} \
            --variantType "somatic" \
            --variantCallMode "tonly"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VARIANT ANNOTATION ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 7. MAF POST-PROCESSING |---------------------------------------------------------------------#
    if [ "${RUN_MUTECT2_TONLY_MAF}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.MAF_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --nType "TS" \
            --writeSomaticVarinatsOnlyMAF TRUE \
            --DP ${MAF_DP} \
            --ALT ${MAF_ALT} \
            --POPAF ${MAF_AF} \
            --writeDB TRUE \
            --variantCallMode "tonly" \
            --GERMLINE FALSE      
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 8. QC SUMMARY |-----------------------------------------------------------------------------#
    if [ "${RUN_QC_SUMMARY}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ QC SUMMARY ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------# 
        Rscript ${RUN_SCRIPT}/WES_run.MultiQC_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --statusQC BAM \
            --summaryFastqScreen TRUE \
            --importDB TRUE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ QC SUMMARY ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 9. TMB CALCULATION |-------------------------------------------------------------------------#
    if [ "${RUN_TMB_TONLY}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.TMB_Calculation.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --nType "TS" \
            --useReferenceMethod TRUE \
            --useReferenceSize TRUE \
            --variantCallMode "tonly" \
            --importDB TRUE
        #----------------------------------------------------------------------#
        TMB_RES_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/tmb
        mv ${TMB_RES_DIR}/${SEQ_ID}.TMB.calculation.result.txt ${TMB_RES_DIR}/${SEQ_ID}.T-only.TMB.calculation.result.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] FINISHED." >> ${PIPELINE_LOG}    
    fi
#--------------------------------------------------------------------------------------------------#
#---| 10. MSI CALCULATION : TUMOR-ONLY MODE |-------------------------------------------------------#
    if [ "${RUN_MSI_TONLY}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSI CALCULATION : TUMOR-ONLY ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.MSI_Calculation.sh \
            --seqFolder ${SEQ_FOLDER} \
            --tumorSeqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --targetBed ${BED_PREFIX}.target.bed \
            --runMode 'tonly' \
            --assembly ${GENOME_ASSEMBLY}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.MSI_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --runMode 'tonly' \
            --nType "TS" \
            --importDB TRUE
        #----------------------------------------------------------------------#
        MSI_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/msi
        #----------------------------------------------------------------------#
        mv ${MSI_DIR}/${SEQ_ID}.MSI.status.calculation.results.txt ${MSI_DIR}/${SEQ_ID}.T-only.MSI.status.calculation.results.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSI CALCULATION : TUMOR-ONLY ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 11. HLA-TYPING |-----------------------------------------------------------------------------#
    if [ "${RUN_HLA}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ HLA-TYPING ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.HLA_Typing.sh \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --threads ${THREADS} \
            --inputType DNA
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.HLA_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --importDB TRUE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ HLA-TYPING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 12. DEEPVARIANT GERMLINE VARIANT CALLING |---------------------------------------------------#
    if [ "${RUN_DEEPVARIANT}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ DEEPVARIANT GERMLINE VARIANT CALLING ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.Deepvariant_GermlineVariantCall.sh \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --bed ${BED_PREFIX}.target.bed \
            --modelType WES \
            --threads 4 \
            --assembly ${GENOME_ASSEMBLY}
        #----------------------------------------------------------------------#    
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ DEEPVARIANT GERMLINE VARIANT CALLING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 13. VARIANT ANNOTATION : GERMLINE |----------------------------------------------------------#
    if [ "${RUN_DEEPVARIANT_ANNOTATION}" = "ON" ]; then      
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VARIANT ANNOTATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.VariantAnnotation.sh \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --nType "TS" \
            --threads ${THREADS} \
            --variantType "germline" \
            --assembly ${GENOME_ASSEMBLY}
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ VARIANT ANNOTATION ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 14. MAF POST-PROCESSING : GERMLINE |---------------------------------------------------------#
    if [ "${RUN_DEEPVARIANT_MAF}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.MAF_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --nType "TS" \
            --writeSomaticVarinatsOnlyMAF FALSE \
            --DP ${MAF_DP} \
            --ALT ${MAF_ALT} \
            --POPAF ${MAF_AF} \
            --variantCallMode "tonly" \
            --writeDB FALSE \
            --GERMLINE TRUE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 15. RUN ASCAT |------------------------------------------------------------------------------#
    if [ "${RUN_CNV_ASCAT}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ ASCAT CNV ANALYSIS ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_ASCAT_TumorOnly.R \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --tumorSeqID ${SEQ_ID} \
            --wesLibKit ${WES_LIB_KIT} \
            --threads ${THREADS} \
            --nType "TS" 
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ ASCAT CNV ANALYSIS ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 16. RUN PURECN |-----------------------------------------------------------------------------#
    if [ "${RUN_CNV_PURECN}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PURECN CNV ANALYSIS ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_PureCN.sh \
            --seqFolder ${SEQ_FOLDER} \
            --baseDir ${BASE_DIR} \
            --tumorSeqID ${SEQ_ID} \
            --runMode 'tonly' \
            --assembly ${GENOME_ASSEMBLY} \
            --wesLibKit ${WES_LIB_KIT} \
            --callMinDepth 20 \
            --nType "TS" \
            --forceCoverage "false"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ PURECN CNV ANALYSIS ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 17. RUN SUPERFREQ |--------------------------------------------------------------------------#
    if [ "${RUN_CNV_SUPERFREQ}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SUPERFREQ CNV ANALYSIS ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
            Rscript ${RUN_SCRIPT}/WES_run.CNV_superFreq.R \
                --baseDir ${BASE_DIR} \
                --seqFolder ${SEQ_FOLDER} \
                --tumorSeqID ${SEQ_ID} \
                --assembly ${GENOME_ASSEMBLY} \
                --wesLibKit ${WES_LIB_KIT} \
                --mode "exome" \
                --threads ${THREADS} \
                --nType "TS" \
                --setID ${SEQ_ID}.TumorOnly
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ SUPERFREQ CNV ANALYSIS ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 18. CNVkit : CNV Tonly Mode |----------------------------------------------------------------#
    if [ "${RUN_CNV_CNVKIT_TONLY}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNVKIT : CNVkit T-only Mode ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        ${RUN_SCRIPT}/WES_run.CNV_CNVkit.sh \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --tumorSeqID ${SEQ_ID} \
            --runMode "tonly" \
            --wesLibKit ${WES_LIB_KIT} \
            --assembly  ${GENOME_ASSEMBLY} \
            --threads ${THREADS} \
            --nType "TS" \
            --superfreqSetID "TumorOnly"
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNVKIT : CNVkit T-only Mode ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#

echo "Running Pipeline FINISHED." >> ${PIPELINE_LOG}
echo "--------------------------------------------------------------------------" >> ${PIPELINE_LOG}
echo "   " >> ${PIPELINE_LOG}
