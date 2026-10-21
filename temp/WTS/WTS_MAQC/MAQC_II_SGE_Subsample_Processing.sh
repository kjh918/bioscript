
#! /bin/bash

#$ -N MAQC_II
#$ -S /bin/bash
#$ -pe mpi 15
#$ -t 1-4
#$ -tc 4
#$ -e /storage/home/kangsm/analysis/gcx.rnaseq.maqc/no_rep_subsampling_test/meta/MAQC_II_Processing_20251222.err
#$ -o /storage/home/kangsm/analysis/gcx.rnaseq.maqc/no_rep_subsampling_test/meta/MAQC_II_Processing_20251222.out

# >>>>> USER IMPUT START
#--------------------------------------------------------------------------------------------------#
BASE_DIR="/storage/home/kangsm/analysis/gcx.rnaseq.maqc"
SEQ_FOLDER="no_rep_subsampling_test"
INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/sample.list
# SEQ_ID SEQ_TYPE PAIRED THREADS STRAND GENOME_ASSEMBLY #
ANALYSIS_CONFIG_YAML=
READ_LENGTH=100
#------------------------------------------------------------------------------#
RUN_FASTQC="OFF"
RUN_FASTQ_TO_BAM="ON"
RUN_RSEM_QUANT="ON"
RUN_BAM_QC="ON"
RUN_FUSION_GENES="OFF"
CREATE_EXPR_MATRIX="OFF"
RUN_WTS_ANALYSIS="OFF"
#--------------------------------------------------------------------------------------------------#
# <<<<< USER IMPUT END

#--------------------------------------------------------------------------------------------------#
SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
SEQ_TYPE=$(echo ${SEQ_INFO} | awk '{print $2}')
PAIRED=$(echo ${SEQ_INFO} | awk '{print $3}')
THREADS=$(echo ${SEQ_INFO} | awk '{print $4}')
STRAND=$(echo ${SEQ_INFO} | awk '{print $5}')
GENOME_ASSEMBLY=$(echo ${SEQ_INFO} | awk '{print $6}')
#--------------------------------------------------------------------------------------------------#
RUN_SCRIPT=/storage/home/kangsm/runScripts
#--------------------------------------------------------------------------------------------------#

#---| FASTQC |-------------------------------------------------------------------------------------#
    if [ "${RUN_FASTQC}" = "ON" ]; then 
        ${RUN_SCRIPT}/WTS_run.FastqQC.sh \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --readLength ${READ_LENGTH} \
            --threads ${THREADS} \
            --runMultiqc "false" \
            --seqType ${SEQ_TYPE}
    fi
#--------------------------------------------------------------------------------------------------#

#---| FASTQ_TO_BAM |-------------------------------------------------------------------------------#
    if [ "${RUN_FASTQ_TO_BAM}" = "ON" ]; then 
        ${RUN_SCRIPT}/WTS_run.FastqToBAM.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS ${THREADS} \
            --READ_LENGTH ${READ_LENGTH} \
            --STRANDED ${STRAND} \
            --GENOME_ASSEMBLY ${GENOME_ASSEMBLY} \
            --PAIRED ${PAIRED}
    fi
#--------------------------------------------------------------------------------------------------#

#---| RSEM_QUANTIFICATION |------------------------------------------------------------------------#
    if [ "${RUN_RSEM_QUANT}" = "ON" ]; then 
        ${RUN_SCRIPT}/WTS_run.Quantification_RSEM.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS ${THREADS} \
            --READ_LENGTH ${READ_LENGTH} \
            --STRANDED ${STRAND} \
            --GENOME_ASSEMBLY ${GENOME_ASSEMBLY} \
            --PAIRED ${PAIRED}
    fi
#--------------------------------------------------------------------------------------------------#

#---| BAM QC |-------------------------------------------------------------------------------------#
    if [ "${RUN_BAM_QC}" = "ON" ]; then 
        ${RUN_SCRIPT}/WTS_run.BAM_QC.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS ${THREADS} \
            --STRANDED ${STRAND} \
            --GENOME_ASSEMBLY ${GENOME_ASSEMBLY} \
            --PAIRED ${PAIRED} \
            --RUN_MULTIQC "true" \
            --QC_RES_DB_IMPORT "false"
    fi
#--------------------------------------------------------------------------------------------------#

#---| FUSION_GENES |-------------------------------------------------------------------------------#
    if [ "${RUN_FUSION_GENES}" = "ON" ]; then 
        ${RUN_SCRIPT}/WTS_run.FusionGenes_and_SpliceVariants.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS ${THREADS} \
            --GENOME_ASSEMBLY ${GENOME_ASSEMBLY} \
            --RUN_ARRIBA "true" \
            --RUN_STARFUSION "true" \
            --RUN_SPLICEVARIANT "true" \
            --DRAW_ARRIBA_PLOTS "true" \
            --PAIRED_END ${PAIRED}
    fi
#--------------------------------------------------------------------------------------------------#

#---| CREATE_EXPR_MATRIX |-------------------------------------------------------------------------#
    if [ "${CREATE_EXPR_MATRIX}" = "ON" ]; then 
        Rscript ${RUN_SCRIPT}/WTS_run.CreateExprValuesObjects.R \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --ANALYSIS_CONFIG ${ANALYSIS_CONFIG_YAML} \
            --DB_IMPORT TRUE
    fi
#--------------------------------------------------------------------------------------------------#
 
#---| WTS ANALYSIS RUN |---------------------------------------------------------------------------#
    if [ "${RUN_WTS_ANALYSIS}" = "ON" ]; then 
        Rscript ${RUN_SCRIPT}/WTS_run.WTS_AnalysisRun.R \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --ANALYSIS_CONFIG ${ANALYSIS_CONFIG_YAML} \
            --GENOME_ASSEMBLY ${GENOME_ASSEMBLY}
    fi
#--------------------------------------------------------------------------------------------------#
