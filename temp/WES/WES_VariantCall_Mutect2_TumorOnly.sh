
#! /usr/bin/bash

    DATA_DIR=/BiO/DATA/WES/23WES02
    SAMPLE_LIST=$(cat ${DATA_DIR}/sample.list)
    SAMPLE_NUMBERS=$(echo ${SAMPLE_LIST} | wc -w)
    if [[ "${JOBS}" -eq "1"  &&  "${SAMPLE_NUMBERS}" -eq "1" ]]; then
        JOBS=${JOBS}
    else
        JOBS=${SAMPLE_NUMBERS}
    fi


#------------------------------------------------------------------------------#

    BAM_DIR=${DATA_DIR}/bam_files
    QC_RES_DIR=${DATA_DIR}/qc_files
    VCF_RES_DIR=${DATA_DIR}/vcf_files
    FINAL_RES_DIR=${DATA_DIR}/Results
    TMP_DIR=${DATA_DIR}/tmp
#------------------------------------------------------------------------------#
    REF_DIR="/storage/references_and_index/gcx.hg19"
    REF_GENOME_FASTA=${REF_DIR}/gcx.hg19.fasta
    VCF_GERMLINE=${REF_DIR}/vcf/gcx.hg19.af-only-gnomad.raw.sites.vcf.gz
    VCF_PON=${REF_DIR}/vcf/gcx.hg19.mutect2-exome-panel.vcf.gz
#

#---| Generate pileup summaries on tumor sample
java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx3G -jar /storage/apps/bin/gatk-local.jar GetPileupSummaries \
    -I ${BAM_DIR}/17-100X.sorted.dup.marked.recal.bam \
    -O ${BAM_DIR}/17-100X.targeted_sequencing.table \
    -V /storage/references_and_index/gcx.hg19/vcf/gcx.hg19.af-only-gnomad.raw.sites.vcf.gz \
    -L /storage/DB/gnomAD/v2.1/exome_calling_regions.v1.interval_list \
    -R ${REF_GENOME_FASTA}

#---| Calculate contamination on tumor sample
java -XX:+UseParallelGC -XX:ParallelGCThreads=10 -Xmx3G -jar /storage/apps/bin/gatk-local.jar CalculateContamination \
    -I ${BAM_DIR}/17-100X.targeted_sequencing.table \
    -O ${BAM_DIR}/17-100X.targeted_sequencing.contamination.table

#---| Run MuTect2 using only tumor sample
/storage/apps/bin/gatk Mutect2 \
    -R ${REF_GENOME_FASTA} \
    -I ${BAM_DIR}/17-100X.sorted.dup.marked.recal.bam \
    -O ${VCF_RES_DIR}/17-100X.mutect2.T-only.wPoN.rev.vcf \
    --af-of-alleles-not-in-resource 2.5e-06 \
    --germline-resource ${VCF_GERMLINE} \
    --panel-of-normals ${VCF_PON}

#---| Mutect2 Run with PoN |
    parallel -j ${JOBS} -k "/BiO/SWTools/gatk Mutect2 \
        -I ${BAM_DIR}/{}.sorted.dup.marked.recal.bam \
        -R ${REF_GENOME_FASTA} \
        --germline-resource ${VCF_GERMLINE} \
        --panel-of-normals ${VCF_PON} \
        --tmp-dir ${TMP_DIR} \
        -O ${VCF_RES_DIR}/{}.mutect2.T-only.wPoN.vcf" ::: ${SAMPLE_LIST}
#---|

#---| VCF FILTERING


#---| Mutect2 VCF Filtering |
    parallel -j ${JOBS} -k "/BiO/SWTools/gatk FilterMutectCalls \
        -V ${VCF_RES_DIR}/{}.mutect2.T-only.wPoN.vcf \
        -R ${REF_GENOME_FASTA} \
        --QUIET true \
        --tmp-dir ${TMP_DIR} \ 
        -O ${VCF_RES_DIR}/{}.mutect2.T-only.wPoN.filtered.vcf" ::: ${SAMPLE_LIST}
#---|
    parallel -j ${JOBS} -k "/BiO/SWTools/gatk FilterMutectCalls \
        -V ${VCF_RES_DIR}/{}.mutect2.T-only.vcf \
        -R ${REF_GENOME_FASTA} \
        --QUIET true \
        --tmp-dir ${TMP_DIR} \ 
        -O ${VCF_RES_DIR}/{}.mutect2.T-only.filtered.vcf" ::: ${SAMPLE_LIST}

SINGULARITY_SIF=/storage/home/kangsm/Downloads/vep.sif

singularity run /storage/home/kangsm/Downloads/vep.sif vep --offline --vcf --dir_cache "/storage/references_and_index/gcx.hg19/vep-cache" --species homo_sapiens --assembly GRCh37 -o "/data/wes/23WES02/vcf_files/17-100X.mutect2.T-only.wPoN.vep3.vcf" -i "/data/wes/23WES02/vcf_files/17-100X.mutect2.T-only.wPoN.vcf"

singularity run "${SINGULARITY_SIF}"/ensembl-vep_v99.2.sif vep --offline --vcf --dir_cache "${RVANNOT}/vep_core/v99" --species homo_sapiens --assembly GRCh37 -o ffe4bb51-e98a-41a7-a4e1-c3970386889c.consensus.20161006.somatic.indel.vep.vcf -offline -i /projects/verhaak-lab/ecdna/datasets/pcawg/dump/final_consensus_12oct/icgc/indel/ffe4bb51-e98a-41a7-a4e1-c3970386889c.consensus.20161006.somatic.indel.vcf.gz

singularity exec /storage/home/kangsm/Downloads/vep.sif INSTALL.pl -c /storage/references_and_index/gcx.hg19/vep-cache -a cf -s homo_sapiens -y GRCh37

CACHE_DIR=/storage/references_and_index/gcx.hg19/vep-cache
INPUT_DIR=/data/wes/23WES02/vcf_files
CUSTOM_DIR=/storage/DB/clinvar/GRCh37
REF_DIR=/storage/references_and_index/gcx.hg19
singularity exec -B ${CACHE_DIR}:/vep-cache,${INPUT_DIR}:/vep-input,${CUSTOM_DIR}:/vep-custom,${REF_DIR}:/vep-ref /storage/home/kangsm/Downloads/vep.sif \
  vep --cache --merged --vcf --force_overwrite --af_gnomad \
      --dir /vep-cache \
      --assembly GRCh37 \
      --input_file /vep-input/17-100X.mutect2.T-only.wPoN.vcf \
      --output_file /vep-input/17-100X.mutect2.T-only.wPoN.vep3.vcf \
      --fasta /vep-ref/gcx.hg19.fasta \
      --custom /vep-custom/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN



#---| ENSEMBL-VEP RUN

    parallel -j ${JOBS} -k "/storage/apps/bin/vep --cache --merged --force_overwrite \
        --assembly GRCh37 \
        --stats_text \
        --dir /storage/DB/ensembl-vep-cache/release-110 \
        --fasta /storage/references_and_index/gcx.hg19/gcx.hg19.fasta \
        --custom /storage/DB/clinvar/GRCh37/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN \
        --custom /storage/DB/gnomAD/v2.1/gnomad.exomes.r2.1.sites.INFO_ANNOTATIONS_FIXED.vcf.bgz,gnomAD,vcf,exact,0,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth \
        --input_file ${VCF_RES_DIR}/{}.mutect2.T-only.wPoN.vcf \
        --vcf --output_file ${VCF_RES_DIR}/{}.mutect2.T-only.wPoN.vep.vcf" ::: ${SAMPLE_LIST}


    parallel -j ${JOBS} -k "perl /storage/apps/vcf2maf-1.6.21/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf /data/wes/23WES02/vcf_files/17-100X.mutect2.T-only.wPoN.vep3.vcf \
        --output-maf /data/wes/23WES02/vcf_files/17-100X.mutect2.T-only.wPoN.vep3.maf \
        --tumor-id 17-100X \
        --ncbi-build GRCh37 \
        --ref-fasta /storage/references_and_index/gcx.hg19/gcx.hg19.fasta \
        --verbose false" ::: ${SAMPLE_LIST}

perl /storage/apps/vcf2maf-1.6.21/vcf2maf.pl \
        --inhibit-vep \
        --input-vcf /data/wes/23WES02/vcf_files/17-100X.mutect2.T-only.wPoN.vep3.vcf \
        --output-maf /data/wes/23WES02/vcf_files/17-100X.mutect2.T-only.wPoN.vep3.maf \
        --tumor-id 17-100X \
        --ncbi-build GRCh37 \
        --ref-fasta /storage/references_and_index/gcx.hg19/gcx.hg19.fasta \
        --verbose false

#---| GATK FUNCOTATOR

    source /storage/apps/miniconda3/bin/activate gatk 
    parallel -j ${JOBS} -k "/storage/apps/bin/gatk Funcotator \
        --variant ${VCF_RES_DIR}/{}.mutect2.T-only.wPoN.vcf \
        --reference /storage/references_and_index/gcx.hg19/gcx.hg19.fasta \
        --ref-version hg19 \
        --data-sources-path /storage/DB/gatk-funcotator/funcotator_dataSources.v1.7.20200521s \
        --output ${VCF_RES_DIR}/{}.mutect2.T-only.wPoN.gatk.funcotator.maf \
        --output-file-format MAF" ::: ${SAMPLE_LIST}


--vep-path /opt/vep --vep-data /srv/vep

#---/ Compression and Indexing
    parallel -j 7 -k "bgzip -c ${VCF_RES_DIR}/{}.mutect2.single.vcf > ${VCF_RES_DIR}/{}.mutect2.single.vcf.gz" ::: ${SAMPLE_LIST}
    parallel -j 7 -k "tabix -p vcf ${VCF_RES_DIR}/{}.mutect2.single.vcf.gz" ::: ${SAMPLE_LIST}

#---/ 




#---/ VCF FILTERING WITH BED-FILE
    parallel -j 6 -k "/BiO/SWTools/bedtools intersect \
        -a ${VCF_RES_DIR}/{}.mutect2.single.vcf \
        -b /BiO/DB/Twist-Exome-2.0/b37_exome_v2.0.2_targets_sorted_validated.re_annotated.bed > ${VCF_RES_DIR}/{}.mutect2.single.bed.filtered.vcf" ::: ${SAMPLE_LIST}

# VCF FILTERING



VCF=/BiO/DATA/WES/23WES02/vcf_out/filter_test/3-100X.mutect2.single.vcf.gz
OUT=/BiO/DATA/WES/23WES02/vcf_out/filter_test/3-100X.mutect2.single.stat.test

# allele frequency
/BiO/SWTools/vcftools --gzvcf ${VCF} --freq2 --out ${OUT} --max-alleles 2 
# mean depth per individual
/BiO/SWTools/vcftools --gzvcf ${VCF} --depth --out ${OUT} 
# mean depth per site
/BiO/SWTools/vcftools --gzvcf ${VCF} --site-mean-depth --out ${OUT} 
# site quality
/BiO/SWTools/vcftools --gzvcf ${VCF} --site-quality --out ${OUT} 
# proportion of missing data per individual
/BiO/SWTools/vcftools --gzvcf ${VCF} --missing-indv --out ${OUT} 
# proportion of missing data per site
/BiO/SWTools/vcftools --gzvcf ${VCF} --missing-site --out ${OUT} 
# heterozygosity
/BiO/SWTools/vcftools --gzvcf ${VCF} --het --out ${OUT} 
