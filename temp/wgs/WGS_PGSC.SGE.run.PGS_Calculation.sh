
#! /bin/bash

#$ -N PGSC_JOB_NAME
#$ -S /bin/bash
#$ -pe mpi 8
#$ -t 1-31
#$ -tc 1
#$ -e ERROR_LOG
#$ -o OUT_LOG

#--------------------------------------------------------------------------------------------------#
# >>
    SEQ_ID=PGSC_SEQ_ID
    WORK_DIR=PGSC_RUN_DIR
    INPUT_VCF=${WORK_DIR}/${SEQ_ID}.PGSC.READY.vcf.gz
    GENOME_BUILD=PGSC_GENOME
# <<
#--------------------------------------------------------------------------------------------------#
    GCX_PGC_LIST=/storage/home/kangsm/PGS_calculation/GCX_PGS_v2.list
#--------------------------------------------------------------------------------------------------#

    RUN_LOG=${WORK_DIR}/pgsc_run.log
    if [ ! -f $RUN_LOG ]; then touch $RUN_LOG; fi
    echo " " >> $RUN_LOG

    PGS_INFO=$(sed -n -e "${SGE_TASK_ID} p" ${GCX_PGC_LIST})
    PGS_ID=$(echo $PGS_INFO | awk '{print $1}')
    PHENOTYPE=$(echo $PGS_INFO | awk '{print $2}')

#---| shamplesheet |
    echo -e "sampleset,path_prefix,chrom,format,vcf_genotype_field" > ${WORK_DIR}/${PGS_ID}_samplesheet.csv
    echo -e "PGS-${PHENOTYPE},${WORK_DIR}/${INPUT_VCF},,vcf,GT" >> ${WORK_DIR}/${PGS_ID}_samplesheet.csv

#---| RUN PGSC |
    cd ${WORK_DIR}

    /storage/home/kangsm/apps/bin/nextflow run pgscatalog/pgsc_calc \
        -profile singularity \
        --input ${WORK_DIR}/${PGS_ID}_samplesheet.csv \
        --target_build ${GENOME_BUILD} \
        --outdir ${WORK_DIR}/results \
        --pgs_id ${PGS_ID} \
        --min_overlap 0.001 \
        --run_ancestry /storage/home/kangsm/myDB/PGS_Catalog/pgsc_HGDP+1kGP_v1.tar.zst \
        --max_cpus 8
        
#---| LOG 
    echo -e "$(date) PGSC : ${PGS_ID} - ${PHENOTYPE} >> DONE." >> $RUN_LOG

#--------------------------------------------------------------------------------------------------#

