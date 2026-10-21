

    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_CTAT_SPLICING="${RUN_SINGULARITY}/ctat-splicing-0.0.2.sif /usr/local/src/CTAT-SPLICING"

    CTAT_RESOURCE_LIB=/storage/references_and_index/hg19/star_fusion/CTAT_resource_lib


#---| INSTALL RESOURCE LIB |---------------------------------------------------#
    ${RUN_CTAT_SPLICING}/prep_genome_lib/ctat-splicing-lib-integration.py \
        --cancer_introns_tsv /storage/references_and_index/hg19/star_fusion/cancer_introns.GRCh37.Jun232020.tsv \
        --genome_lib_dir ${CTAT_RESOURCE_LIB}
#------------------------------------------------------------------------------#
    