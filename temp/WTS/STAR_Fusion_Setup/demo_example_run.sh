
#---| RUN SETUP
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_STAR_FUSION="${RUN_SINGULARITY}/star-fusion-0.13.0.sif STAR-Fusion"

    CTAT_RESOURCE_LIB=/storage/references_and_index/hg19/star_fusion/CTAT_resource_lib
    DEMO_DIR=/storage/home/kangsm/analysis/star-fusion/internal_demo_run

#---| FASTQ INPUT RUN TEST   
    ${RUN_STAR_FUSION} \
        --genome_lib_dir ${CTAT_RESOURCE_LIB} \
        --left_fq ${DEMO_DIR}/reads_1.fq.gz \
        --right_fq ${DEMO_DIR}/reads_2.fq.gz \
        --output_dir ${DEMO_DIR}/fastq_input_result \
        --CPU 8

#---| JUNCTION INPUT RUN TEST
    ${RUN_STAR_FUSION} \
        --genome_lib_dir ${CTAT_RESOURCE_LIB} \
        --chimeric_junction ${DEMO_DIR}/hg19.Chimeric.out.junction \
        --output_dir ${DEMO_DIR}/junction_input_result \
        --CPU 8

#--------------------------------------------------------------------------------------------------#

    ${RUN_STAR_FUSION} \
        --genome_lib_dir ${CTAT_RESOURCE_LIB} \
        --left_fq ${DEMO_DIR}/reads_1.fq.gz \
        --right_fq ${DEMO_DIR}/reads_2.fq.gz \
        --output_dir ${DEMO_DIR}/fastq_input_result \
        --CPU 8 \
        --FusionInspector validate



