
#! /bin/bash

#$ -N WTSTEST01
#$ -S /bin/bash
#$ -pe mpi 25
#$ -t 1-10
#$ -tc 4
#$ -e /storage/home/kangsm/sge_log/WTSTEST01/WTSTEST01_rnaseq.err
#$ -o /storage/home/kangsm/sge_log/WTSTEST01/WTSTEST01_rnaseq.out

#---| Run Script
##
SEQ_FOLDER="WTSTEST01"  # <-- This is what you input manually
SEQ_TYPE="wts"
##
INPUT_LIST_FILE=/data/${SEQ_TYPE}/${SEQ_FOLDER}/meta/sample.list
SEQ_ID=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
##
export LANG="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
/storage/home/kangsm/runScripts/gcx_wts_v1.0_star_rsem_pipeline.sh \
    --seqFolder ${SEQ_FOLDER} \
    --seqID ${SEQ_ID} \
    --logFile "${SEQ_ID}.processing.log" \
    --readLength 100 \
    --threads 25 \
    --stranded "no"



