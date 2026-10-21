
#! /bin/bash

#$ -N cnv_test_run
#$ -S /bin/bash
#$ -pe mpi 1
#$ -t 1-16
#$ -tc 7
#$ -e /storage/home/kangsm/sge_log/cnv_test_run/cnv_test_run.err
#$ -o /storage/home/kangsm/sge_log/cnv_test_run/cnv_test_run.out

RUN_CONFIG=/storage/home/kangsm/analysis/cnv_test/wes_cnv_results_02_new/meta/cnv_test_run_list2.txt
CNV_RUN=$(sed -n -e "${SGE_TASK_ID} p" ${RUN_CONFIG})
#------------------------------------------------------------------------------#
tumorBamDir=/storage/home/kangsm/analysis/cnv_test/tumor_bams
normalBamDir=/storage/home/kangsm/analysis/cnv_test/normal_bams
#------------------------------------------------------------------------------#
TUMOR_BAM_ID=$(echo $CNV_RUN | awk '{print $1}')
TUMOR_BAM=${tumorBamDir}/${TUMOR_BAM_ID}.analysis.ready.bam
#------------------------------------------------------------------------------#
NORMAL_BAM_ID=$(echo $CNV_RUN | awk '{print $2}')
NORMAL_BAM=${normalBamDir}/${NORMAL_BAM_ID}.analysis.ready.bam
#------------------------------------------------------------------------------#
PON_BAM_DIR=$(echo $CNV_RUN | awk '{print $3}')
#------------------------------------------------------------------------------#
SAMPLE_ID=$(echo $CNV_RUN | awk '{print $4}')
#------------------------------------------------------------------------------#
RUN_CALLER=$(echo $CNV_RUN | awk '{print $5}')
#------------------------------------------------------------------------------#
RUN_MODE=$(echo $CNV_RUN | awk '{print $6}')
#------------------------------------------------------------------------------#
CNV_RES_DIR=$(echo $CNV_RUN | awk '{print $7}')
#------------------------------------------------------------------------------#
WES_LIB=$(echo $CNV_RUN | awk '{print $8}')
#------------------------------------------------------------------------------#
REF_VERSION=$(echo $CNV_RUN | awk '{print $9}')
#------------------------------------------------------------------------------#
SCRIPT_DIR=/storage/home/kangsm/runScripts
#------------------------------------------------------------------------------#
####
    source /storage/apps/miniconda3/bin/activate cnv

    ${SCRIPT_DIR}/WES.v1.2_runCNV_${RUN_CALLER}.sh \
        --tumor ${TUMOR_BAM} \
        --normal ${NORMAL_BAM} \
        --ponDir ${PON_BAM_DIR} \
        --sampleId ${SAMPLE_ID} \
        --runMode ${RUN_MODE} \
        --runDir ${CNV_RES_DIR} \
        --resultDir ${CNV_RES_DIR} \
        --bed ${WES_LIB} \
        --genomeVersion ${REF_VERSION}
####




