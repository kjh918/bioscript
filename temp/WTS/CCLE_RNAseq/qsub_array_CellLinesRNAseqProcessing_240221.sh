
#! /bin/bash

#$ -N CCLE_240221
#$ -S /bin/bash
#$ -pe mpi 10
#$ -t 1-22
#$ -tc 5
#$ -e /storage/home/kangsm/sge_log/CCLE/CCLE_240221.err
#$ -o /storage/home/kangsm/sge_log/CCLE/CCLE_240221.out

#--------------------------------------------------------------------------------------------------#
BASE_DIR="/data/wts"
#--------------------------------------------------------------------------------------------------#
INPUT_LIST_FILE="/storage/home/kangsm/analysis/ccle_rnaseq_processing/sample_240221.list"
SAMPLE_ID=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
#--------------------------------------------------------------------------------------------------#
SEQ_FOLDER=$(echo ${SAMPLE_ID} | awk '{print $1}')
SEQ_ID=$(echo ${SAMPLE_ID} | awk '{print $2}')
READ_LENGTH=$(echo ${SAMPLE_ID} | awk '{print $3}')
#--------------------------------------------------------------------------------------------------#
SCRIPT_DIR="/storage/home/kangsm/runScripts"
NCPU=10
#--------------------------------------------------------------------------------------------------#

# [ RUN STAR-RSEM PIPELINE ]
#--------------------------------------------------------------------------------------------------#
    ${SCRIPT_DIR}/WTS.v1.0_FastqToStarRsem_v0.1.sh \
        --seqFolder ${SEQ_FOLDER} \
        --seqID ${SEQ_ID} \
        --logFile ${SEQ_ID}.STAR-RSEM.pipeline.processing.log \
        --threads ${NCPU} \
        --baseDir ${BASE_DIR} \
        --readLength ${READ_LENGTH} \
        --stranded "no"
#--------------------------------------------------------------------------------------------------#

