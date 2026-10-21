
#! /bin/bash

#$ -N hg38_superfreq
#$ -S /bin/bash
#$ -pe mpi 15
#$ -t 2-7
#$ -tc 7
#$ -e /storage/home/kangsm/sge_log/cnv_seqc2/paired_cnv.01.err
#$ -o /storage/home/kangsm/sge_log/cnv_seqc2/paired_cnv.01.out

#--------------------------------------------------------------------------------------------------#
BASE_DIR=/storage/home/kangsm/analysis/cnv_seqc2
SEQ_FOLDER=hg38
INPUT_LIST_FILE=${BASE_DIR}/${SEQ_FOLDER}/meta/NT_match.list
#--------------------------------------------------------------------------------------------------#
SEQ_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${INPUT_LIST_FILE})
TUMOR_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $2}')
NORMAL_SEQ_ID=$(echo ${SEQ_INFO} | awk '{print $1}')
SEQ_DEPTH=$(echo ${SEQ_INFO} | awk '{print $3}')
#--------------------------------------------------------------------------------------------------#
THREADS=15
GENOME_ASSEMBLY="hg38"
WES_LIB_KIT="sureselect.v6"
#--------------------------------------------------------------------------------------------------#
RUN_SCRIPT=/storage/home/kangsm/runScripts
BED_PREFIX=/storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
#--------------------------------------------------------------------------------------------------#
    # ln -s ${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/vcf/${TUMOR_SEQ_ID}.mutect2.NT.purecn.filtered.vcf ${BASE_DIR}/${SEQ_FOLDER}/${TUMOR_SEQ_ID}/vcf/${TUMOR_SEQ_ID}.mutect2.NT.keep.germline.filtered.vcf
#--------------------------------------------------------------------------------------------------#

    # ${RUN_SCRIPT}/WES_run.CNV_PureCN.sh \
    #     --seqFolder ${SEQ_FOLDER} \
    #     --baseDir ${BASE_DIR} \
    #     --tumorSeqID ${TUMOR_SEQ_ID} \
    #     --normalSeqID ${NORMAL_SEQ_ID} \
    #     --runMode 'nt' \
    #     --assembly ${GENOME_ASSEMBLY} \
    #     --wesLibKit ${WES_LIB_KIT} \
    #     --callMinDepth 20
#--------------------------------------------------------------------------------------------------#
#     Rscript ${RUN_SCRIPT}/WES_run.CNV_superFreq_rev.R \
#         --baseDir ${BASE_DIR} \
#         --seqFolder ${SEQ_FOLDER} \
#         --tumorSeqID ${TUMOR_SEQ_ID} \
#         --normalSeqID ${NORMAL_SEQ_ID} \
#         --setID ${SEQ_ID}.NT \
#         --genomeAssembly ${GENOME_ASSEMBLY} \
#         --wesLibKit ${WES_LIB_KIT} \
#         --mode "exome" \
#         --threads ${THREADS} \
#         --vcfTag "mutect2.keep.germline.filtered"
# #--------------------------------------------------------------------------------------------------#
#     Rscript /storage/home/kangsm/runScripts/WES_run.CNV_ASCAT.R \
#         --baseDir ${BASE_DIR} \
#         --seqFolder ${SEQ_FOLDER} \
#         --tumorSeqID ${TUMOR_SEQ_ID} \
#         --normalSeqID ${NORMAL_SEQ_ID}  \
#         --bed ${BED_PREFIX}.target.bed \
#         --assembly ${GENOME_ASSEMBLY} \
#         --threads ${THREADS}
#--------------------------------------------------------------------------------------------------#

    ${RUN_SCRIPT}/WES_run.CNV_CNVkit.sh \
        --tumorSeqID ${TUMOR_SEQ_ID} \
        --normalSeqID ${NORMAL_SEQ_ID} \
        --baseDir ${BASE_DIR} \
        --seqFolder ${SEQ_FOLDER} \
        --runMode 'nt' \
        --wesLibKit ${WES_LIB_KIT} \
        --assembly ${GENOME_ASSEMBLY} \
        --threads ${THREADS}
