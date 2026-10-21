
#! /bin/bash

#$ -N wts_download
#$ -S /bin/bash
#$ -pe mpi 1
#$ -t 1-38
#$ -tc 8
#$ -e /storage/home/kangsm/analysis/geo_organoid_rnaseq/download_log/download.err
#$ -o /storage/home/kangsm/analysis/geo_organoid_rnaseq/download_log/download.out

#--------------------------------------------------------------------------------------------------#
DOWNLOAD_LIST=/storage/home/kangsm/analysis/geo_organoid_rnaseq/gse_info/download.list
SAMPLE_LIST=$(sed -n -e "${SGE_TASK_ID} p" ${DOWNLOAD_LIST})
SAMPLE_NAME=$(echo ${SAMPLE_LIST} | awk '{print $2}')
SRR_ID=$(echo ${SAMPLE_LIST} | awk '{print $1}')
#--------------------------------------------------------------------------------------------------#
DOWNLOAD_DIR=/storage/home/kangsm/analysis/geo_organoid_rnaseq/fastq
LOG_FILE=/storage/home/kangsm/analysis/geo_organoid_rnaseq/download_log/download.log
#--------------------------------------------------------------------------------------------------#
echo "$(date '+%Y-%m-%d %H:%M:%S') | ${SRR_ID} Download Start." >> ${LOG_FILE}
#--------------------------------------------------------------------------------------------------#
/storage/apps/sratoolkit.3.0.10/bin/fasterq-dump ${SRR_ID} \
    --outdir ${DOWNLOAD_DIR} \
    --threads 12 \
    --split-3 \
    --temp ${DOWNLOAD_DIR}/tmp \
    --force
#--------------------------------------------------------------------------------------------------#
echo "$(date '+%Y-%m-%d %H:%M:%S') | ${SRR_ID} Download Finished." >> ${LOG_FILE}
#--------------------------------------------------------------------------------------------------#

# RENAME AND COMPRESSION
# for SGE_TASK_ID in $(seq 2 12)
# do
#     SAMPLE_LIST=$(sed -n -e "${SGE_TASK_ID} p" ${DOWNLOAD_LIST})
#     SAMPLE_NAME=$(echo ${SAMPLE_LIST} | awk '{print $1}')
#     SRR_ID=$(echo ${SAMPLE_LIST} | awk '{print $2}')
#     bgzip ${DOWNLOAD_DIR}/${SRR_ID}_1.fastq -c -@ 30 > ${DOWNLOAD_DIR}/${SAMPLE_NAME}_R1.fastq.gz
#     bgzip ${DOWNLOAD_DIR}/${SRR_ID}_2.fastq -c -@ 30 > ${DOWNLOAD_DIR}/${SAMPLE_NAME}_R2.fastq.gz
#     echo "${SRR_ID} > ${SAMPLE_NAME} DONE"
# done
#--------------------------------------------------------------------------------------------------#
# rename 's/WES_FD_N_1/FD1_N/' WES_FD_N_1*
# rename 's/WES_FD_N_2/FD2_N/' WES_FD_N_2*
# rename 's/WES_FD_N_3/FD3_N/' WES_FD_N_3*
# rename 's/WES_FD_T_1/FD1_T/' WES_FD_T_1*
# rename 's/WES_FD_T_2/FD2_T/' WES_FD_T_2*
# rename 's/WES_FD_T_3/FD3_T/' WES_FD_T_3*
# rename 's/WES_NV_N_1/NV1_N/' WES_NV_N_1*
# rename 's/WES_NV_N_2/NV2_N/' WES_NV_N_2*
# rename 's/WES_NV_N_3/NV3_N/' WES_NV_N_3*
# rename 's/WES_NV_T_1/NV1_T/' WES_NV_T_1*
# rename 's/WES_NV_T_2/NV2_T/' WES_NV_T_2*
# rename 's/WES_NV_T_3/NV3_T/' WES_NV_T_3*
#--------------------------------------------------------------------------------------------------#
