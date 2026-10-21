

#---|
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK4="${RUN_SINGULARITY}/gatk-4.4.0.0.sif gatk"
#--------------------------------------------------------------------------------------------------#

    GENOME_ASSEMBLY="hg38"
    WES_LIB_KIT="sureselect.v6"


    if [ "${GENOME_ASSEMBLY}" = "hg38" ]; then
        GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    else
        GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    fi
    BED_PREFIX=/storage/references_and_index/${GENOME_ASSEMBLY}/bed
    CNV_GATK_INTERVALS=/storage/home/kangsm/myDB/cnv_gatk_intervals

    ${GATK4} PreprocessIntervals \
        -L ${BED_PREFIX}/Agilent_SureSelect_v6/hg38_sureselect.v6.target.bed \
        -R ${GENOME_FASTA} \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${CNV_GATK_INTERVALS}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.target_gatk4.wes.cnv.interval_list

    ${GATK4} AnnotateIntervals \
          -R ${GENOME_FASTA} \
          -L /storage/references_and_index/${GENOME_ASSEMBLY}/bed/${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.target.interval_list \
          --interval-merging-rule OVERLAPPING_ONLY \
          -O ${GATK_CNV_INTERVALS}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_wes.cnv.annotated.intervals.tsv


#---|
    BASE_DIR=/storage/home/kangsm/analysis
    SEQ_FOLDER=crc.nto.5samples
    SEQ_ID=
    BIN_SIZE=200
    GENOME_ASSEMBLY="hg38"
    WES_LIB_KIT="sureselect.v6"
    #--------------------------------------------------------------------------#
    GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
    GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
    #--------------------------------------------------------------------------#
    GATK_CNV_INTERVALS=/storage/home/kangsm/myDB/cnv_gatk_intervals
    BIN_INTERVALS=${GATK_CNV_INTERVALS}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.target_gatk4.wes.cnv.interval_list
    #--------------------------------------------------------------------------#
    PON_HDF5_DIR=/data/pon/pon_cnv.gatk/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
    REV="0.1"
    PON_HDF5=${PON_HDF5_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_gatk.wes.cnv.pon.r${REV}.hdf5

    if [ ! -d ${PON_HDF5} ]; then mkdir -p ${PON_HDF5}; fi
    #--------------------------------------------------------------------------#

    for SEQ_ID in FD1_N FD2_N FD3_N NV1_N NV2_N NV3_N
    do
        ${GATK4} CollectReadCounts \
            --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=14 -Xmx16384m" \
            --input ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam/${SEQ_ID}.analysisReady.bam \
            --intervals ${BIN_INTERVALS} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${PON_HDF5_DIR}/${SEQ_ID}.pon.read.counts.hdf5
    done
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#

    ${GATK4} CreateReadCountPanelOfNormals \
        --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=14 -Xmx32768m" \
        --minimum-interval-median-percentile 5.0 \
        --maximum-zeros-in-sample-percentage 5.0 \
        --number-of-eigensamples 20 \
        -I ${PON_HDF5_DIR}/FD1_N.pon.read.counts.hdf5 \
        -I ${PON_HDF5_DIR}/FD2_N.pon.read.counts.hdf5 \
        -I ${PON_HDF5_DIR}/FD3_N.pon.read.counts.hdf5 \
        -I ${PON_HDF5_DIR}/NV1_N.pon.read.counts.hdf5 \
        -I ${PON_HDF5_DIR}/NV2_N.pon.read.counts.hdf5 \
        -I ${PON_HDF5_DIR}/NV3_N.pon.read.counts.hdf5 \
        -O ${PON_HDF5}

        --annotated-intervals ${GATK_CNV_INTERVALS}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_wes.cnv.annotated.intervals.tsv \

GERMLINE_VCF=/storage/references_and_index/hg19/vcf/af-only-gnomad.raw.sites.grch37.vcf.gz


BED=/storage/references_and_index/hg19/bed/Twist_WES_2.0/GRCh37_twist.exome.2.0.target.bed
BED_VCF=/storage/home/kangsm/myDB/cnv_gatk_intervals/GRCh37_twist.exome.2.0_gnomad.raw.sites.vcf
bcftools view -R ${BED} ${GERMLINE_VCF} > ${BED_VCF}

BED=/storage/references_and_index/hg19/bed/Agilent_SureSelect_v6/GRCh37_sureselect.v6.target.bed
BED_VCF=/storage/home/kangsm/myDB/cnv_gatk_intervals/GRCh37_sureselect.v6_gnomad.raw.sites.vcf
bcftools view -R ${BED} ${GERMLINE_VCF} > ${BED_VCF}


BED=/storage/references_and_index/hg19/bed/Agilent_SureSelect_v6_UTR/GRCh37_sureselect.v6.utr.target.bed
BED_VCF=/storage/home/kangsm/myDB/cnv_gatk_intervals/GRCh37_sureselect.v6.utr_gnomad.raw.sites.vcf
bcftools view -R ${BED} ${GERMLINE_VCF} > ${BED_VCF}


BED=/storage/references_and_index/hg19/bed/Agilent_SureSelect_v8/GRCh37_sureselect.v8.target.bed
BED_VCF=/storage/home/kangsm/myDB/cnv_gatk_intervals/GRCh37_sureselect.v8_gnomad.raw.sites.vcf
bcftools view -R ${BED} ${GERMLINE_VCF} > ${BED_VCF}

GERMLINE_VCF=/storage/references_and_index/hg38/vcf/af-only-gnomad.hg38.vcf.gz

BED=/storage/references_and_index/hg38/bed/Agilent_SureSelect_v6/hg38_sureselect.v6.target.bed
BED_VCF=/storage/home/kangsm/myDB/cnv_gatk_intervals/hg38_sureselect.v6_gnomad.raw.sites.vcf
bcftools view -R ${BED} ${GERMLINE_VCF} > ${BED_VCF}

BED=/storage/references_and_index/hg38/bed/Agilent_SureSelect_v6/hg38_sureselect.v6.target.bed
BED_VCF=/storage/home/kangsm/myDB/cnv_gatk_intervals/hg38_sureselect.v6_gnomad.raw.sites.vcf
bcftools view -R ${BED} ${GERMLINE_VCF} > ${BED_VCF}










