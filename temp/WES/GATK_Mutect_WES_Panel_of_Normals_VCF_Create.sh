
#! /bin/bash

#---/ Run CMD 
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK4="${RUN_SINGULARITY}/gatk-4.4.0.0.sif gatk"
#
#---/ Default hg19 Params
    GENOME_VERSION="hg19"
    GENOME_FASTA="/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta"
    GERMLINE_VCF="/storage/references_and_index/hg19/vcf/af-only-gnomad.raw.sites.grch37.vcf.gz"
    NGS_LIB=
    TARGET_INTERVAL=/storage/references_and_index/hg19/bed/${NGS_LIB}/hg19_${NGS_LIB}.target.interval_list
#
#
#---/ Manual Input Params
    GENOME_VERSION="canFam3.1"
    GENOME_FASTA="/storage/references_and_index/canFam3.1/fasta/canFam3.1.fasta"
    NGS_LIB="canfam3.sureselect.all.exon.v2"
    GERMLINE_VCF="/storage/references_and_index/canFam3.1/vcf/canFam3.dbsnp151.vcf.gz"
    TARGET_INTERVAL="/storage/references_and_index/canFam3.1/bed/canfam3.sureselect.all.exon.v2.target.interval_list"  
#
#
#---/ Panel of Normals VCF Processing
    BASE_DIR=/storage/home/kangsm/analysis/gcx.podo.canis.cancer.wes/03.WES_DATA_RESULTS
    BATCH_ID=CANINE_PON
    SAMPLE_LIST_FILE=/storage/home/kangsm/analysis/gcx.podo.canis.cancer.wes/03.WES_DATA_RESULTS/CANINE_PON/meta/sample_id_only.list
    echo -n "VCF processing  " 
    for SEQ_ID in $(cat $SAMPLE_LIST_FILE)
    do
        VCF_DIR=${BASE_DIR}/${BATCH_ID}/${SEQ_ID}/vcf
        grep ^# ${VCF_DIR}/${SEQ_ID}.mutect2.bias.filtered.vcf > ${VCF_DIR}/${SEQ_ID}.mutect2.panel.of.normals.filtered.pass.vcf
        bcftools query -i "FILTER='PASS' && FMT/DP>25 && FMT/AD[0:1]>3 && INFO/artiStatus='snv'" -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT" \
            ${VCF_DIR}/${SEQ_ID}.mutect2.bias.filtered.vcf >> ${VCF_DIR}/${SEQ_ID}.mutect2.panel.of.normals.filtered.pass.vcf
        bgzip -c ${VCF_DIR}/${SEQ_ID}.mutect2.panel.of.normals.filtered.pass.vcf > ${VCF_DIR}/${SEQ_ID}.mutect2.panel.of.normals.filtered.pass.vcf.gz
        tabix -p vcf ${VCF_DIR}/${SEQ_ID}.mutect2.panel.of.normals.filtered.pass.vcf.gz
        echo -n -e "...${SEQ_ID}"
    done
    echo "DONE."
#
#
#---/ Sample Map File
    SAMPLE_MAP_FILE=/data/pon/pon_meta/canfam3.sureselect.all.exon.v2.pon.vcf.sample.list
    if [ -f "${SAMPLE_MAP_FILE}" ]; then 
        rm $SAMPLE_MAP_FILE 
        touch ${SAMPLE_MAP_FILE}
    else
        touch ${SAMPLE_MAP_FILE}
    fi
    ###
    for  SEQ_ID in $(cat $SAMPLE_LIST_FILE)
    do
        echo -e "${SEQ_ID}\t${BASE_DIR}/${BATCH_ID}/${SEQ_ID}/vcf/${SEQ_ID}.mutect2.panel.of.normals.filtered.pass.vcf.gz" >> ${SAMPLE_MAP_FILE}
    done
#
#
#---/ Genomic DB Import 
    cd /data/pon
    GENOMIC_DB_NAME=pondb_canfam3.sureselect.all.exon.v2

    ${GATK4} GenomicsDBImport --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=14 -Xmx32768m" \
        --reference ${GENOME_FASTA} --intervals ${TARGET_INTERVAL} \
        --genomicsdb-workspace-path ${GENOMIC_DB_NAME} \
        --sample-name-map ${SAMPLE_MAP_FILE} \
        --merge-input-intervals true \
        --batch-size 10 --tmp-dir /data/pon/tmp
#
#---/ Panel of Normals VCF Crreation
    OUTPUT_PON_VCF=/data/pon/pon_vcf/canfam3.sureselect.all.exon.v2.gatk.mutect2.panel.of.normals.vcf.gz

    ${GATK4} CreateSomaticPanelOfNormals \
        --reference ${GENOME_FASTA} \
        --germline-resource ${GERMLINE_VCF} \
        -V gendb://${GENOMIC_DB_NAME} \
        -O ${OUTPUT_PON_VCF}
