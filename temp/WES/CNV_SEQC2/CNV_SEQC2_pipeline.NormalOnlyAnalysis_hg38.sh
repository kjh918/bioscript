
#! /bin/bash

#$ -N seqc2_normal_only_hg38
#$ -S /bin/bash
#$ -pe mpi 20
#$ -t 1-6
#$ -tc 6
#$ -e /storage/home/kangsm/analysis/cnv_seqc2/sge_log/normal.only.analysis.20240406.err
#$ -o /storage/home/kangsm/analysis/cnv_seqc2/sge_log/normal.only.analysis.20240406.out

#--------------------------------------------------------------------------------------------------#
BASE_DIR=/storage/home/kangsm/analysis/cnv_seqc2
SEQ_FOLDER=hg38
INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/normal.list
#--------------------------------------------------------------------------------------------------#
SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
#SEQ_DEPTH=$(echo ${SEQ_INFO} | awk '{print $2}')
#FLOWCELL_ID=$(echo ${SEQ_INFO} | awk '{print $3}')
#READ_LENGTH=$(echo ${SEQ_INFO} | awk '{print $4}')
#--------------------------------------------------------------------------------------------------#
THREADS=20
GENOME_ASSEMBLY="hg38"
WES_LIB_KIT="sureselect.v6"
#--------------------------------------------------------------------------------------------------#
RGLB=${WES_LIB_KIT}
RUN_SCRIPT=/storage/home/kangsm/runScripts
BED_PREFIX=/storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
MAF_DP=25
MAF_ALT=3
MAF_AF=0.01
#--------------------------------------------------------------------------------------------------#
RUN_FASTQ_SCREEN="OFF"
RUN_FASTQ_QC="OFF"
RUN_FASTQ_TO_BAM="OFF"
RUN_BAM_QC="OFF"
RUN_MUTECT2_TONLY="ON"
RUN_MUTECT2_TONLY_ANNOTATION="OFF"
RUN_MUTECT2_TONLY_MAF="OFF"
RUN_TMB_TONLY="OFF"
RUN_MSI_TONLY="OFF"
RUN_HLA="OFF"
RUN_DEEPVARIANT="OFF"
RUN_QC_SUMMARY="OFF"
RUN_CNV_SUPERFREQ="ON"
#--------------------------------------------------------------------------------------------------#

#---| PIPELINE RUN LOG |---------------------------------------------------------------------------#
    PIPELINE_LOG_DIR=${BASE_DIR}/${SEQ_FOLDER}/meta/analysis_run_log
    if [ ! -d ${PIPELINE_LOG_DIR} ]; then mkdir -p ${PIPELINE_LOG_DIR}; fi
    #--------------------------------------------------------------------------#
    PIPELINE_LOG=${PIPELINE_LOG_DIR}/${SEQ_ID}.normal.only.analysis.hg38.$(date '+%Y-%m-%d').log
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
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
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
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
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
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
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
            --intervalTarget ${BED_PREFIX}target.interval_list \
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
            --keepGermlineCall "false"
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
            --varCallMode "tonly" \
            --assembly ${GENOME_ASSEMBLY}
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
            --varCallMode "tonly" \
            --writeSomaticVarinatsOnlyMAF TRUE \
            --DP ${MAF_DP} \
            --ALT ${MAF_ALT} \
            --POPAF ${MAF_AF} \
            --writeDB FALSE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MAF POST-PROCESSING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 8. TMB CALCULATION |-------------------------------------------------------------------------#
    if [ "${RUN_TMB_TONLY}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
        Rscript ${RUN_SCRIPT}/WES_run.TMB_Calculation.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --useReference TRUE \
            --useReferenceSize TRUE \
            --variantCallMode "tonly" \
            --maf ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/vcf/${SEQ_ID}.mutect2.filtered.somatic.variants.only.maf             
        #----------------------------------------------------------------------#
        mv ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/tmb/${SEQ_ID}.TMB.calculation.result.txt ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/tmb/${SEQ_ID}.T-only.TMB.calculation.result.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ TMB CALCULATION ] FINISHED." >> ${PIPELINE_LOG}    
    fi
#--------------------------------------------------------------------------------------------------#
#---| 9. MSI CALCULATION : TUMOR-ONLY MODE |-------------------------------------------------------#
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
            --importDB FALSE
        #----------------------------------------------------------------------#
        MSI_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/msi
        #----------------------------------------------------------------------#
        mv ${MSI_DIR}/${SEQ_ID}.MSI.status.calculation.results.txt ${MSI_DIR}/${SEQ_ID}.T-only.MSI.status.calculation.results.txt
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ MSI CALCULATION : TUMOR-ONLY ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 10. HLA-TYPING |-----------------------------------------------------------------------------#
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
            --importDB FALSE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ HLA-TYPING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 11. DEEPVARIANT GERMLINE VARIANT CALLING |---------------------------------------------------#
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
            --assembly hg19
        #----------------------------------------------------------------------#    
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ DEEPVARIANT GERMLINE VARIANT CALLING ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#
#---| 12. QC SUMMARY |-----------------------------------------------------------------------------#
    if [ "${RUN_QC_SUMMARY}" = "ON" ]; then        
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ QC SUMMARY ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------# 
        Rscript ${RUN_SCRIPT}/WES_run.MultiQC_PostProcessing.R \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --baseDir ${BASE_DIR} \
            --statusQC BAM \
            --summaryFastqScreen TRUE \
            --importDB FALSE
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ QC SUMMARY ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#

#---| 13. CNV : SUPERFREQ |----------------------------------------------------------------------------#
    if [ "${RUN_CNV_SUPERFREQ}" = "ON" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : SUPERFREQ ] START." >> ${PIPELINE_LOG}
        #----------------------------------------------------------------------#
            Rscript ${RUN_SCRIPT}/WES_run.CNV_superFreq_rev.R \
                --baseDir ${BASE_DIR} \
                --seqFolder ${SEQ_FOLDER} \
                --tumorSeqID ${SEQ_ID} \
                --setID ${SEQ_ID} \
                --genomeAssembly ${GENOME_ASSEMBLY} \
                --wesLibKit ${WES_LIB_KIT} \
                --mode "exome" \
                --threads ${THREADS}
        #----------------------------------------------------------------------#
        echo "$(date '+%Y-%m-%d %H:%M:%S') | [ CNV : SUPERFREQ ] FINISHED." >> ${PIPELINE_LOG}
    fi
#--------------------------------------------------------------------------------------------------#

echo "Running Pipeline FINISHED." >> ${PIPELINE_LOG}
echo "--------------------------------------------------------------------------" >> ${PIPELINE_LOG}
echo "   " >> ${PIPELINE_LOG}
