
## CNVKIT FINE-TUNE

## HCC1395_FFPE_N = SRR7890963 = FFX_GZ_N_1h_1
## HCC1395_FFPE_T = SRR7890945 = FFX_GZ_T_1h_1
## HCC1395_N = SRR7890845 = WES_NC_N_1
## HCC1395_T = SRR7890844 = WES_NC_T_1


# 1. VARIANT CALLING : MUTECT2 Paired-Mode
#------------------------------------------------------------------------------#
    TARGET_INTERVAL=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.extended.target.interval_list
    REF_DIR=/storage/references_and_index/hg38
    GENOME_FASTA=${REF_DIR}/fasta/Homo_sapiens_assembly38.fasta
    VCF_GNOMAD=${REF_DIR}/vcf/af-only-gnomad.hg38.vcf.gz
    VCF_PON=/storage/home/kangsm/myDB/vcf/gatk_1kg_panel_of_normal_hg38.vcf.gz

    VEP_CACHE_DIR=/storage/references_and_index/hg38/vep
    VEP_PLUGIN_DIR=${VEP_CACHE_DIR}/Plugins/110
    VEP_BAM=${VEP_CACHE_DIR}/bam/GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bam

#------------------------------------------------------------------------------#
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    GATK4="${RUN_SINGULARITY}/gatk-4.4.0.0.sif gatk"
    RUN_VEP="singularity exec -B /storage,/data /storage/images/vep-110.1.sif /opt/vep/src/ensembl-vep/vep"
    GATK3="${RUN_SINGULARITY}/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar"

#------------------------------------------------------------------------------#
    NORMAL_BAM_FILE=/storage/home/kangsm/analysis/cnv_test2/normal_bam/HCC1395_N.analysis.ready.bam
    TUMOR_BAM_FILE=/storage/home/kangsm/analysis/cnv_test2/tumor_bam/HCC1395_T.analysis.ready.bam
    QC_DIR=/storage/home/kangsm/analysis/cnv_test2/qc_tmp
    VCF_DIR=/storage/home/kangsm/analysis/cnv_test2/vcf
#------------------------------------------------------------------------------#
    ${GATK4} Mutect2 \
        --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=14 -Xmx32768m" \
        -pairHMM AVX_LOGLESS_CACHING_OMP \
        --native-pair-hmm-threads 15 \
        --reference ${GENOME_FASTA} \
        --intervals ${TARGET_INTERVAL} \
        --input ${TUMOR_BAM_FILE} \
        --input ${NORMAL_BAM_FILE} \
        -normal "SRR7890845" \
        --germline-resource ${VCF_GNOMAD} \
        --panel-of-normals ${VCF_PON} \
        --f1r2-tar-gz ${QC_DIR}/HCC1395_T_rev.f1r2.tar.gz \
        --min-base-quality-score 20 \
        --base-quality-score-threshold 20 \
        --interval-padding 50 \
        --genotype-germline-sites true \
        --genotype-pon-sites true \
        --output ${VCF_DIR}/HCC1395_T.mutect2_rev.vcf

    ${GATK4} LearnReadOrientationModel \
        --java-options "-XX:ParallelGCThreads=14 -Xmx8192m" \
        -I ${QC_DIR}/HCC1395_T_rev.f1r2.tar.gz \
        -O ${QC_DIR}/HCC1395_T_rev.read-orientation-model.tar.gz

    ${GATK4} GetPileupSummaries \
        --java-options "-XX:ParallelGCThreads=14 -Xmx8192m" \
        -I ${TUMOR_BAM_FILE} \
        -O ${QC_DIR}/HCC1395_T.targeted_sequencing.table \
        -V ${VCF_GNOMAD} \
        -L ${TARGET_INTERVAL} \
        -R ${GENOME_FASTA}

    ${GATK4} GetPileupSummaries \
        --java-options "-XX:ParallelGCThreads=14 -Xmx8192m" \
        -I ${NORMAL_BAM_FILE} \
        -O ${QC_DIR}/HCC1395_N.targeted_sequencing.table \
        -V ${VCF_GNOMAD} \
        -L ${TARGET_INTERVAL} \
        -R ${GENOME_FASTA}

    ${GATK4} CalculateContamination \
        --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        -I ${QC_DIR}/HCC1395_T.targeted_sequencing.table \
        -matched ${QC_DIR}/HCC1395_N.targeted_sequencing.table \
        -O ${QC_DIR}/HCC1395_T.targeted_sequencing.contamination.table

    ${GATK4} FilterMutectCalls \
        --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        -V ${VCF_DIR}/HCC1395_T.mutect2_rev.vcf \
        -L ${TARGET_INTERVAL} \
        --reference ${GENOME_FASTA} \
        --contamination-table ${QC_DIR}/HCC1395_T.targeted_sequencing.contamination.table \
        --ob-priors ${QC_DIR}/HCC1395_T_rev.read-orientation-model.tar.gz \
        -O ${VCF_DIR}/HCC1395_T.mutect2_rev.filtered.vcf

    ${GATK3} -T DepthOfCoverage \
        -I ${NORMAL_BAM_FILE} \
        -L ${TARGET_INTERVAL} \
        -R ${GENOME_FASTA} \
        -dt BY_SAMPLE \
        -dcov 5000 \
        -l INFO \
        --omitDepthOutputAtEachBase \
        --omitLocusTable \
        --minBaseQuality 0 \
        --minMappingQuality 20 \
        --start 1 \
        --stop 5000 \
        --nBins 200 \
        --includeRefNSites \
        --countType COUNT_FRAGMENTS \
        -o /storage/home/kangsm/analysis/cnv_test2/purecn/HCC1395_N.depth.cov

    ${GATK3} -T DepthOfCoverage \
        -I ${TUMOR_BAM_FILE} \
        -L ${TARGET_INTERVAL} \
        -R ${GENOME_FASTA} \
        -dt BY_SAMPLE \
        -dcov 5000 \
        -l INFO \
        --omitDepthOutputAtEachBase \
        --omitLocusTable \
        --minBaseQuality 0 \
        --minMappingQuality 20 \
        --start 1 \
        --stop 5000 \
        --nBins 200 \
        --includeRefNSites \
        --countType COUNT_FRAGMENTS \
        -o /storage/home/kangsm/analysis/cnv_test2/purecn/HCC1395_T.depth.cov





#--------------------------------------------------------------------------------------------------#
    # hg38 mappability
    wget https://s3.amazonaws.com/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_76.bw
    wget https://s3.amazonaws.com/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw

    PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata

    BAIT_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.probe.bed
    REF_DIR=/storage/references_and_index/hg38
    GENOME_FASTA=${REF_DIR}/fasta/Homo_sapiens_assembly38.fasta

    OUT_INTERVAL_DIR=/storage/home/kangsm/myDB/PureCN/interval_files
    MAP_FILE_DIR=/storage/home/kangsm/myDB/PureCN/mappability_files
    MAPPABILITY_FILE=${MAP_FILE_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw
    
    Rscript $PURECN/IntervalFile.R \
        --in-file ${BAIT_BED} \
        --fasta ${GENOME_FASTA} \
        --out-file ${OUT_INTERVAL_DIR}/hg38_sureselect.v6_baits_intervals.txt \
        --off-target --genome hg38 \
        --export ${OUT_INTERVAL_DIR}/hg38_sureselect.v6_baits_optimized.bed \
        --mappability ${MAPPABILITY_FILE}

    # internal segmentation
    # coverage
    COV_N=/storage/home/kangsm/analysis/cnv_test2/purecn/coverage_normal
    COV_T=/storage/home/kangsm/analysis/cnv_test2/purecn/coverage_tumor

    Rscript $PURECN/Coverage.R \
        --out-dir ${COV_N} \
        --bam ${NORMAL_BAM_FILE} \
        --intervals ${OUT_INTERVAL_DIR}/hg38_sureselect.v6.utr_baits_intervals.txt

    Rscript $PURECN/Coverage.R \
        --out-dir ${COV_N} \
        --bam /storage/home/kangsm/analysis/cnv_test2/normal_bam/HCC1395_FFPE_N.analysis.ready.bam \
        --intervals ${OUT_INTERVAL_DIR}/hg38_sureselect.v6.utr_baits_intervals.txt


    Rscript $PURECN/Coverage.R \
        --out-dir ${COV_T} \
        --bam ${TUMOR_BAM_FILE} \
        --intervals ${OUT_INTERVAL_DIR}/hg38_sureselect.v6.utr_baits_intervals.txt

    # normaldb
$GATK Mutect2 \
 -R $REF_GENOME_b37 \
 -I ${GERMLINE_PFX}_G.bam \
 --max-mnp-distance 0 \
 -O ${GERMLINE_PFX}.vcf.gz


$GATK GenomicsDBImport \
 -R $REF_GENOME_b37 \
-L $BEDFILE \
 --genomicsdb-workspace-path pon_db \
 -V ${GERMLINE_PFX}.vcf.gz \ #sample 1 germline
.
.
.
 -V${GERMLINE_PFX}.vcf.gz \ #sample 9 germline


$GATK CreateSomaticPanelOfNormals \
 -R $REF_GENOME_b37 \
 --germline-resource af-only-gnomad.raw.sites.vcf \
 -V gendb://pon_db \
 -O pon.vcf.gz

#---| PureCN NORMAL DB |---------------------------------------------------------------------------#
    BASE_DIR=/storage/home/kangsm/analysis
    SEQ_FOLDER="crc.nto.5samples"
    GENOME_ASSEMBLY="hg19"
    WES_LIB_KIT="sureselect.v6"
    NORMAL_SAMPLE_LIST="CRC01_N CRC02_N CRC03_N CRC04_N CRC05_N"

    PURECN=/storage/apps/R/4.3.1/lib/R/library/PureCN/extdata
    INTERVALS=/storage/home/kangsm/myDB/cnv_purecn/interval_files/${GENOME_ASSEMBLY}_${WES_LIB_KIT}.baits_intervals.txt

    # 1) COVERAGES of NORMALS
    for SEQ_ID in ${NORMAL_SAMPLE_LIST}
    do
        PURECN_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/cnv/purecn
        COVERAGE_RES=${PURECN_DIR}/${SEQ_ID}.analysisReady_coverage_loess.txt.gz

        if [ ! -f ${COVERAGE_RES} ]; then
            if [ ! -d ${PURECN_DIR} ]; then mkdir -p ${PURECN_DIR}; fi
            #------------------------------------------------------------------#
            Rscript $PURECN/Coverage.R \
                --out-dir ${PURECN_DIR} \
                --bam ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam/${SEQ_ID}.analysisReady.bam \
                --intervals ${INTERVALS}
        fi
    done

    # 2) CREATE NORMAL COVERAGE LISTS
    NORMALDB_DIR=/storage/home/kangsm/myDB/cnv_purecn/normalDB/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
    if [ ! -d ${NORMALDB_DIR} ]; then mkdir -p ${NORMALDB_DIR}; fi
    #--------------------------------------------------------------------------#
    COV_LIST_FILE=/storage/home/kangsm/myDB/cnv_purecn/normal_coverage_list/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_normal_coverages.list
    if [ ! -f ${COV_LIST_FILE} ]; then touch ${COV_LIST_FILE}; fi    
    #--------------------------------------------------------------------------#
    for SEQ_ID in ${NORMAL_SAMPLE_LIST}
    do
        PURECN_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/cnv/purecn
        echo "${PURECN_DIR}/${SEQ_ID}.analysisReady_coverage_loess.txt.gz" >> ${COV_LIST_FILE}
    done

    # 3) CREATE PoN VCF
    RUN_SCRIPT=/storage/home/kangsm/runScripts
    for SEQ_ID in ${NORMAL_SAMPLE_LIST}
    do
        VCF_DIR=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/vcf
        Rscript ${RUN_SCRIPT}/WES_run.CNV_CreatePureCNPoNVCF.R --SEQ_ID ${SEQ_ID} --VCF_DIR ${VCF_DIR}
    done
    
    PURECN_PON_VCF_DIR=/storage/home/kangsm/myDB/cnv_purecn/pon_vcf/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
    for SEQ_ID in ${NORMAL_SAMPLE_LIST}
    do
        cp ${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/vcf/${SEQ_ID}.purecn.pon.vcf ${PURECN_PON_VCF_DIR}/
    cd 

    for SEQ_ID in ${NORMAL_SAMPLE_LIST}
    do
        bgzip -c ${PURECN_PON_VCF_DIR}/${SEQ_ID}.purecn.pon.vcf > ${PURECN_PON_VCF_DIR}/${SEQ_ID}.purecn.pon.vcf.gz
        tabix -p vcf ${PURECN_PON_VCF_DIR}/${SEQ_ID}.purecn.pon.vcf.gz
    done

    bcftools merge ${PURECN_PON_VCF_DIR}/*.purecn.pon.vcf.gz -o ${PURECN_PON_VCF_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_PureCN.Panel.of.Normal.vcf

    bgzip -c ${PURECN_PON_VCF_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_PureCN.Panel.of.Normal.vcf > ${PURECN_PON_VCF_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_PureCN.Panel.of.Normal.vcf.gz
    tabix -p vcf ${PURECN_PON_VCF_DIR}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_PureCN.Panel.of.Normal.vcf.gz

    # 3) CREATE NORMAL-DB RDS
    NORMALDB_DIR=/storage/home/kangsm/myDB/cnv_purecn/normalDB/${GENOME_ASSEMBLY}_${WES_LIB_KIT}
    PON_VCF=/storage/home/kangsm/myDB/cnv_purecn/pon_vcf/${GENOME_ASSEMBLY}_${WES_LIB_KIT}/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_PureCN.Panel.of.Normal.vcf.gz
    COV_LIST_FILE=/storage/home/kangsm/myDB/cnv_purecn/normal_coverage_list/${GENOME_ASSEMBLY}_${WES_LIB_KIT}_normal_coverages.list
    
    Rscript $PURECN/NormalDB.R \
        --out-dir ${NORMALDB_DIR} \
        --coverage-files ${COV_LIST_FILE} \
        --normal-panel ${PON_VCF} \
        --genome ${GENOME_ASSEMBLY} \
        --assay ${GENOME_ASSEMBLY}_${WES_LIB_KIT} \
        --force
#--------------------------------------------------------------------------------------------------#




    HG38_PON_VCF=/storage/home/kangsm/myDB/vcf/gatk_1kg_panel_of_normal_hg38.vcf.gz
    
    COV_LIST_FILE=/storage/home/kangsm/myDB/PureCN/normal_coverage_list/hg38_sureselect.v6.utr_normal_coverages.list

    ls -a ${COV_N}/*_loess.txt.gz | cat > ${COV_LIST_FILE}

    ${WES_LIB_KIT} \
        --out-dir ${NDB_OUT_DIR} \
        --coverage-files ${COV_LIST_FILE} \
        --genome hg38 \
        --assay hg38_sureselect.v6.utr \
        --force

    Rscript $PURECN/NormalDB.R \
        --outdir ${NDB_OUT_DIR} \
        --normal_panel ${HG38_PON_VCF} \
        --assay hg38.GATK.PoN \
        --genome hg38 \
        --force

    # PureCN
    RES_INTERNAL=/storage/home/kangsm/analysis/cnv_test2/purecn/results_internal_segmentation_02

    Rscript $PURECN/PureCN.R \
        --out ${RES_INTERNAL} \
        --tumor ${COV_T}/HCC1395_T.analysis.ready_coverage_loess.txt.gz \
        --normal ${COV_N}/HCC1395_N.analysis.ready_coverage_loess.txt.gz \
        --sampleid SRR7890844 \
        --vcf ${VCF_DIR}/HCC1395_T.mutect2_rev.filtered.vcf \
        --intervals ${OUT_INTERVAL_DIR}/hg38_sureselect.v6.utr_baits_intervals.txt \
        --genome hg38

    --normaldb ${NDB_OUT_DIR}/normalDB_hg38_sureselect.v6.utr_hg38.rds \

https://bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/Quick.html#41_Coverage

https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/PureCN/inst/doc/Quick.pdf

https://rstudio-pubs-static.s3.amazonaws.com/444222_0505f4eecb74455ca9fda115989dad5f.html

    ### CNVkit
    RUN_CNVKIT="singularity exec -B /storage,/data /storage/images/cnvkit-0.9.10.sif cnvkit.py"
    CNVKIT_RES=/storage/home/kangsm/analysis/cnv_test2/purecn/cnvkit_results
    
    ${RUN_CNVKIT} export seg \
        ${CNVKIT_RES}/HCC1395_T.analysis.ready.cns \
        --enumerate-chroms \
        -o ${CNVKIT_RES}/HCC1395_T.analysis.ready.seg

    PURECN_3RD_RES=/storage/home/kangsm/analysis/cnv_test2/purecn/results_3rd_party_sw
    SIMPLE_REPEAT=/storage/home/kangsm/myDB/PureCN/snp_blacklist/hg38_ucsc.simple.repeat.bed

    Rscript $PURECN/PureCN.R \
        --out ${PURECN_3RD_RES}/HCC1395_T \
        --sampleid HCC1395_T \
        --tumor ${CNVKIT_RES}/HCC1395_T.analysis.ready.cnr \
        --seg-file ${CNVKIT_RES}/HCC1395_T.analysis.ready.seg \
        --vcf ${VCF_DIR}/HCC1395_T.mutect2_rev.filtered.vcf \
        --snp-blacklist ${SIMPLE_REPEAT} \
        --genome hg38 \
        --fun-segmentation Hclust \
        --force --post-optimize --seed 123

#--------------------------------------------------------------------------------------------------#

    ${GATK3} -T CallableLoci \
        -R ${GENOME_FASTA} \
        -I:tumor ${TUMOR_BAM_FILE}  \
        --summary ${COV_T}/HCC1395_T_minDP10_callable_table.txt \
        -o ${COV_T}/HCC1395_T_minDP10_callable_status.bed \
        --minDepth 10

    grep CALLABLE ${COV_T}/HCC1395_T_minDP10_callable_status.bed > ${COV_T}/HCC1395_T_minDP10_callable_status_filtered.bed

    # mutational burden 
    Rscript $PURECN/FilterCallableLoci.R \
    --genome hg38 \
    --in-file ${COV_T}/HCC1395_T_minDP10_callable_status_filtered.bed \
    --out-file ${COV_T}/HCC1395_T_minDP10_callable_status_CDS.bed


    # CNV (INTERNAL)
    Rscript $PURECN/Dx.R \
        --out ${RES_INTERNAL}/SRR7890844 \
        --rds ${RES_INTERNAL}/SRR7890844.rds \
        --callable ${COV_T}/HCC1395_T_minDP30_callable_status_CDS.bed \
        --exclude ${SIMPLE_REPEAT} \
        --signatures \
        --force

    # CNV (CNVkit)
    Rscript $PURECN/Dx.R \
        --out ${PURECN_3RD_RES}/HCC1395_T \
        --rds ${PURECN_3RD_RES}/HCC1395_T.rds \
        --callable ${COV_T}/HCC1395_T_minDP30_callable_status_CDS.bed \
        --exclude ${SIMPLE_REPEAT} \
        --signatures \
        --force


    











