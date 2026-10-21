
#! /bin/bash

## https://github.com/cphgeno/CNVbench

### DATASET2 ###
  # WES_NC_N_1 = SRR7890845
  # WES_NC_T_1 = SRR7890844
  # FFX_GZ_N_1h_1 = SRR7890963
  # FFX_GZ_T_1h_1 = SRR7890945

  cp /data/wes/hg38/ncbi/SRR7890844/bam/SRR7890844.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_T.analysis.ready.bam
  cp /data/wes/hg38/ncbi/SRR7890845/bam/SRR7890845.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_N.analysis.ready.bam
  cp /data/wes/hg38/ncbi/SRR7890844/bam/SRR7890844.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_T.analysis.ready.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890845/bam/SRR7890845.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_N.analysis.ready.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890963/bam/SRR7890963.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_FFPE_N.analysis.ready.bam
  cp /data/wes/hg38/ncbi/SRR7890945/bam/SRR7890945.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_FFPE_T.analysis.ready.bam
  cp /data/wes/hg38/ncbi/SRR7890963/bam/SRR7890963.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_FFPE_N.analysis.ready.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890945/bam/SRR7890945.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2/HCC1395_FFPE_T.analysis.ready.bam.bai

  cp /data/wes/hg38/ncbi/SRR7890963/bam/SRR7890963.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N01.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890966/bam/SRR7890966.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N02.analysis.bamhead
  cp /data/wes/hg38/ncbi/SRR7890952/bam/SRR7890952.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N03.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890951/bam/SRR7890951.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N04.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890954/bam/SRR7890954.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N05.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890954/bam/SRR7890954.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N06.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890963/bam/SRR7890963.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N01.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890966/bam/SRR7890966.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N02.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890952/bam/SRR7890952.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N03.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890951/bam/SRR7890951.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N04.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890954/bam/SRR7890954.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N05.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890954/bam/SRR7890954.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe/HCC1395_FFPE_N06.analysis.bam.bai

  cp /data/wes/hg38/ncbi/SRR7890919/bam/SRR7890919.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N01.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890880/bam/SRR7890880.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N02.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890881/bam/SRR7890881.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N03.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890878/bam/SRR7890878.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N04.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890874/bam/SRR7890874.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N05.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890875/bam/SRR7890875.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N06.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890882/bam/SRR7890882.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N07.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890851/bam/SRR7890851.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N08.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890845/bam/SRR7890845.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N09.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890847/bam/SRR7890847.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N10.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890846/bam/SRR7890846.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N11.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890853/bam/SRR7890853.recaled.bam /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N12.analysis.bam
  cp /data/wes/hg38/ncbi/SRR7890919/bam/SRR7890919.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N01.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890880/bam/SRR7890880.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N02.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890881/bam/SRR7890881.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N03.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890878/bam/SRR7890878.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N04.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890874/bam/SRR7890874.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N05.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890875/bam/SRR7890875.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N06.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890882/bam/SRR7890882.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N07.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890851/bam/SRR7890851.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N08.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890845/bam/SRR7890845.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N09.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890847/bam/SRR7890847.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N10.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890846/bam/SRR7890846.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N11.analysis.bam.bai
  cp /data/wes/hg38/ncbi/SRR7890853/bam/SRR7890853.recaled.bai /storage/home/kangsm/analysis/cnv_test/testset2_normal_cell/HCC1395_N12.analysis.bam.bai

###
### DATASET3 ###
  # /storage/home/kangsm/analysis/cnv_test/testset3
  cp /data/pon/PON02/PON02_01/bam/PON02_01.merged.dup.marked.realign.recal.bam /storage/home/kangsm/analysis/cnv_test/testset3/PON02_01.analysis.ready.bam
  cp /data/pon/PON02/PON02_02/bam/PON02_02.merged.dup.marked.realign.recal.bam /storage/home/kangsm/analysis/cnv_test/testset3/PON02_02.analysis.ready.bam
  cp /data/pon/PON02/PON02_03/bam/PON02_03.merged.dup.marked.realign.recal.bam /storage/home/kangsm/analysis/cnv_test/testset3/PON02_03.analysis.ready.bam
  cp /data/pon/PON02/PON02_04/bam/PON02_04.merged.dup.marked.realign.recal.bam /storage/home/kangsm/analysis/cnv_test/testset3/PON02_04.analysis.ready.bam
  cp /data/pon/PON02/PON02_05/bam/PON02_05.merged.dup.marked.realign.recal.bam /storage/home/kangsm/analysis/cnv_test/testset3/PON02_05.analysis.ready.bam

  cp /data/pon/PON02/PON02_01/bam/PON02_01.merged.dup.marked.realign.recal.bam.bai /storage/home/kangsm/analysis/cnv_test/testset3/PON02_01.analysis.ready.bam.bai
  cp /data/pon/PON02/PON02_02/bam/PON02_02.merged.dup.marked.realign.recal.bai /storage/home/kangsm/analysis/cnv_test/testset3/PON02_02.analysis.ready.bam.bai
  cp /data/pon/PON02/PON02_03/bam/PON02_03.merged.dup.marked.realign.recal.bai /storage/home/kangsm/analysis/cnv_test/testset3/PON02_03.analysis.ready.bam.bai
  cp /data/pon/PON02/PON02_04/bam/PON02_04.merged.dup.marked.realign.recal.bai /storage/home/kangsm/analysis/cnv_test/testset3/PON02_04.analysis.ready.bam.bai
  cp /data/pon/PON02/PON02_05/bam/PON02_05.merged.dup.marked.realign.recal.bai /storage/home/kangsm/analysis/cnv_test/testset3/PON02_05.analysis.ready.bam.bai
###

### TESTSET3 RUN
GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
REF_FLAT=/storage/references_and_index/hg19/extra/human_g1k_v37.refFlat.txt
TWIST_BAIT_BED=/storage/references_and_index/hg19/bed/Twist_WES_2.0/probe.bed
TWIST_TARGET_BED=/storage/references_and_index/hg19/bed/Twist_WES_2.0/target.bed
TWIST_TARGET_INTERVAL=/storage/references_and_index/hg19/bed/Twist_WES_2.0/target.interval_list

NORMAL_BAM_DIR=/storage/home/kangsm/analysis/cnv_test/testset3
TUMOR_BAM_DIR=/storage/home/kangsm/analysis/cnv_test/tumor_bams

NA12878_LIST="NA12878_01 NA12878_02"
GCX_PON_LIST="PON02_01 PON02_02 PON02_03 PON02_04 PON02_05"

### CNVkit ###-------------------------------------------------------------------------------------#
  RUN_CNVKIT="singularity exec -B /storage,/data /storage/images/cnvkit-0.9.10.sif cnvkit.py"
  ACCESS_BED=/storage/home/kangsm/myDB/bed/CNVkit/access-5k-mappable.grch37.bed

  ### CNVkit pooled noemal : PON
  POOLED_NORMAL=/storage/home/kangsm/analysis/cnv_test/cnvkit/pooled_normal
  parallel -j 2 -k "${RUN_CNVKIT} batch ${TUMOR_BAM_DIR}/{}.analysis.ready.bam \
      --normal ${NORMAL_BAM_DIR}/*.analysis.ready.bam \
      --targets ${TWIST_BAIT_BED} \
      --annotate ${REF_FLAT} \
      --fasta ${GENOME_FASTA} \
      --access ${ACCESS_BED} \
      --output-reference ${POOLED_NORMAL}/{}/{}.reference.cnn \
      --output-dir ${POOLED_NORMAL}/{} \
      -p 6 \
      --diagram --scatter" ::: $(echo ${NA12878_LIST} | awk '{print $1,$2}')

###------------------------------------------------------------------------------------------------#

### XHMM ###---------------------------------------------------------------------------------------#
  XHMM_DIR=/storage/home/kangsm/analysis/cnv_test/xhmm
  source /storage/apps/miniconda3/bin/activate cnv
  RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
  GATK3="${RUN_SINGULARITY}/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar"

  TWIST_TARGET_INTERVAL=/storage/references_and_index/hg19/bed/Twist_WES_2.0/target.interval_list
  XHMM_COVERAGE_RES=${XHMM_DIR}/bam_coverages
  # PON-N coverage
    parallel -j 5 -k "${GATK3} -T DepthOfCoverage \
      -I ${NORMAL_BAM_DIR}/{}.analysis.ready.bam \
      -L ${TWIST_TARGET_INTERVAL} \
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
      -o ${XHMM_COVERAGE_RES}/{}.depth.cov" ::: $(echo ${GCX_PON_LIST} | awk '{print $1,$2,$3,$4,$5}')
  
  # extreme GC : TWIST
    ${GATK3} -T GCContentByInterval \
      -L ${TWIST_TARGET_INTERVAL} \
      -R ${GENOME_FASTA} \
      -o ${XHMM_DIR}/grch37.twist.locus_gc.txt
    cat ${XHMM_DIR}/grch37.twist.locus_gc.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > ${XHMM_DIR}/grch37.twist.extreme_gc.targets.txt

  # pooled_normal
    for SEQ_ID in NA12878_01 NA12878_02
    do
      XHMM_POOLED_NORMAL=${XHMM_DIR}/${SEQ_ID}_pooled.normal
      mkdir -p ${XHMM_POOLED_NORMAL}
      xhmm --mergeGATKdepths -o ${XHMM_POOLED_NORMAL}/analysisStart.RD.txt \
        --GATKdepths ${XHMM_COVERAGE_RES}/${SEQ_ID}.depth.cov.sample_interval_summary \
        --GATKdepths ${XHMM_COVERAGE_RES}/PON02_01.depth.cov.sample_interval_summary \
        --GATKdepths ${XHMM_COVERAGE_RES}/PON02_02.depth.cov.sample_interval_summary \
        --GATKdepths ${XHMM_COVERAGE_RES}/PON02_03.depth.cov.sample_interval_summary \
        --GATKdepths ${XHMM_COVERAGE_RES}/PON02_04.depth.cov.sample_interval_summary \
        --GATKdepths ${XHMM_COVERAGE_RES}/PON02_05.depth.cov.sample_interval_summary
    done
  # Run Analysis
    for Nr in "01" "02" 
    do
      XHMM_RUN_DIR=${XHMM_DIR}/NA12878_${Nr}_pooled.normal
      RUN_GC_TARGETS=${XHMM_DIR}/grch37.sure.select.extreme_gc.targets.txt

      {
      # filtering & mean-centering
        xhmm --matrix -r ${XHMM_RUN_DIR}/analysisStart.RD.txt \
          --centerData --centerType target \
          --excludeTargets ${RUN_GC_TARGETS} \
          --minTargetSize 10 --maxTargetSize 10000 --minMeanTargetRD 10 --maxMeanTargetRD 500 \
          --minMeanSampleRD 25 --maxMeanSampleRD 200 --maxSdSampleRD 150 \
          --outputExcludedTargets ${XHMM_RUN_DIR}/meanCentering.excludedTargets.txt \
          --outputExcludedSamples ${XHMM_RUN_DIR}/meanCentering.excludedSamples.txt \
          -o ${XHMM_RUN_DIR}/meanCentering.RD.txt
      # PCA of mean-centered data
        xhmm --PCA -r ${XHMM_RUN_DIR}/meanCentering.RD.txt \
          --PCAfiles ${XHMM_RUN_DIR}/RD_PCA
      # normalize based on PCA
        xhmm --normalize -r ${XHMM_RUN_DIR}/meanCentering.RD.txt \
          --PCAfiles ${XHMM_RUN_DIR}/RD_PCA \
          --normalizeOutput ${XHMM_RUN_DIR}/pcaNormalized.RD.txt \
          --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
      # filter based on normalized z-score
        xhmm --matrix -r ${XHMM_RUN_DIR}/pcaNormalized.RD.txt \
          --centerData --centerType sample --zScoreData --maxSdTargetRD 30 \
          -o ${XHMM_RUN_DIR}/pcaNormalized.zScores.RD.txt \
          --outputExcludedTargets ${XHMM_RUN_DIR}/pcaNormalized.zScore.excludedTargets.txt \
          --outputExcludedSamples ${XHMM_RUN_DIR}/pcaNormalized.zScore.excludedSamples.txt 
      # filter original data with same filtering cutoff
        xhmm --matrix -r ${XHMM_RUN_DIR}/analysisStart.RD.txt \
        --excludeTargets ${XHMM_RUN_DIR}/meanCentering.excludedTargets.txt \
        --excludeTargets ${XHMM_RUN_DIR}/pcaNormalized.zScore.excludedTargets.txt \
        --excludeSamples ${XHMM_RUN_DIR}/meanCentering.excludedSamples.txt \
        --excludeSamples ${XHMM_RUN_DIR}/pcaNormalized.zScore.excludedSamples.txt \
        -o ${XHMM_RUN_DIR}/finalFiltered.RD.txt
      # find CNV
        xhmm --discover -p ${XHMM_DIR}/params.txt \
          -r ${XHMM_RUN_DIR}/pcaNormalized.zScores.RD.txt \
          -R ${XHMM_RUN_DIR}/finalFiltered.RD.txt \
          -c ${XHMM_RUN_DIR}/CNV.analysisResult.xcnv \
          -a ${XHMM_RUN_DIR}/CNV.analysisResult.aux_xcnv \
          -s ${XHMM_RUN_DIR}/CNV.analysisResult
      # genotype
        xhmm --genotype -p ${XHMM_DIR}/params.txt \
          -r ${XHMM_RUN_DIR}/pcaNormalized.zScores.RD.txt \
          -R ${XHMM_RUN_DIR}/finalFiltered.RD.txt \
          -g ${XHMM_RUN_DIR}/CNV.analysisResult.xcnv \
          -F ${GENOME_FASTA} \
          -v ${XHMM_RUN_DIR}/CNV.analysisResult.vcf
      }

    done
    ##
#--------------------------------------------------------------------------------------------------#

### EXCAVATOR2 ###---------------------------------------------------------------------------------#

  ** libpng12.so.0 error
  sudo add-apt-repository ppa:linuxuprising/libpng12
  sudo apt update
  sudo apt install libpng12-0

  ## initialize
  #mv SourceTarget.txt SwtoolsDefaultSourceTarget.txt
  cd /storage/apps/EXCAVATOR2-1.1.2
  mv SourceTarget.txt SwtoolsDefaultSourceTarget.txt
  echo "/storage/apps/EXCAVATOR2-1.1.2/data/gcx.hg19.bw /storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta" > /storage/apps/EXCAVATOR2-1.1.2/SourceTarget.txt

  SOURCE_TARGET=/storage/apps/EXCAVATOR2-1.1.2/SourceTarget.txt
  TWIST_TARGET_BED=/storage/home/kangsm/myDB/bed/hg19/GRCh37_twist.exome.2.0.target.bed
  TWIST_TARGET_NOL_BED=/storage/home/kangsm/myDB/bed/TwistExome2.0/twist_exome2.0.target.non.overlap.bed

  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_BED} twistExome2_w1k 1000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_BED} twistExome2_w5k 5000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_BED} twistExome2_w10k 10000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_BED} twistExome2_w50k 50000 hg19

  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_NOL_BED} twistExome2_NOL_w1k 1000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_NOL_BED} twistExome2_NOL_w5k 5000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_NOL_BED} twistExome2_NOL_w10k 10000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${TWIST_TARGET_NOL_BED} twistExome2_NOL_w50k 50000 hg19

  perl TargetPerla.pl ${SOURCE_TARGET} ${SURESELECT_TARGET_BED} sureSelect.V6_w50K 50000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${SURESELECT_TARGET_BED} sureSelect.V6_w10K 10000 hg19
  perl TargetPerla.pl ${SOURCE_TARGET} ${SURESELECT_TARGET_BED} sureSelect.V6_w100K 100000 hg19


  ## data prepare
    WORK_DIR=/storage/home/kangsm/analysis/cnv_test/excavator2/NA12878_w50k
    TARGET_WINDOW="twistExome2_w50K"
    mkdir -p ${WORK_DIR}
    CONFIG_DIR=${WORK_DIR}/config && mkdir -p ${CONFIG_DIR}
    PREPARE_DATA=${WORK_DIR}/prepareData && mkdir -p ${PREPARE_DATA}
    
  # 1) prepare config
    echo "${TUMOR_BAM_DIR}/NA12878_01.analysis.ready.bam ${PREPARE_DATA}/NA12878_01 NA12878_01" > ${CONFIG_DIR}/prepareData.txt
    echo "${TUMOR_BAM_DIR}/NA12878_02.analysis.ready.bam ${PREPARE_DATA}/NA12878_02 NA12878_02" >> ${CONFIG_DIR}/prepareData.txt
    for N in "01" "02" "03" "04" "05"
    do
      echo "${NORMAL_BAM_DIR}/PON02_${N}.analysis.ready.bam ${PREPARE_DATA}/PON02_${N} PON02_${N}" >> ${CONFIG_DIR}/prepareData.txt
    done
  # 2) prepare data
    NCPU=$(cat ${CONFIG_DIR}/prepareData.txt | wc -l)
    perl EXCAVATORDataPrepare.pl ${CONFIG_DIR}/prepareData.txt --processors ${NCPU} --target ${TARGET_WINDOW} --assembly hg19 --mapq 20

  # 4-1) analysis.2 Pooling-Normal + 1 Sample : config
    for seqID in NA12878_01 NA12878_02
    do
        echo "T1 ${PREPARE_DATA}/${seqID} ${seqID}" > ${CONFIG_DIR}/analysisPooling.${seqID}.txt
        for N in 1 2 3 4 5
        do
          echo "C${N} ${PREPARE_DATA}/PON02_0${N} PON02_0${N}" >> ${CONFIG_DIR}/analysisPooling.${seqID}.txt
        done
    done
  # 4-2) analysis.2 Pooling-Normal + 1 Sample : analysis run
    for seqID in NA12878_01 NA12878_02
    do
      perl EXCAVATORDataAnalysis.pl ${CONFIG_DIR}/analysisPooling.${seqID}.txt \
        --processors ${NCPU} --target ${TARGET_WINDOW} --assembly hg19 --output ${WORK_DIR}/result.pooling.${seqID} --mode pooling --verbose
    done

#--------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------------------------------------#

### TESTSET2 RUN

GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
REF_FLAT=/storage/home/kangsm/myDB/annotation/hg38/hg38.refFlat.txt
AGILENT_BAIT_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.probe.bed
AGILENT_TARGET_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.extended.target.bed
AGILENT_TARGET_INTERVAL=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.extended.target.interval_list

NORMAL_BAM_CELL=/storage/home/kangsm/analysis/cnv_test/testset2_normal_cell
NORMAL_BAM_FFPE=/storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe
TUMOR_BAM_DIR=/storage/home/kangsm/analysis/cnv_test/testset2

SEQC2_LIST="HCC1395 HCC1395_FFPE"


### CNVkit ###-------------------------------------------------------------------------------------#
  RUN_CNVKIT="singularity exec -B /storage,/data /storage/images/cnvkit-0.9.10.sif cnvkit.py"
  CNVKIT_HG38_ACCESS_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_access.10k_CNVkit.bed


  ### CNVkit N-T paired mode run
  PAIRED_RES_DIR=/storage/home/kangsm/analysis/cnv_test/cnvkit/paired
    parallel -j 2 -k "${RUN_CNVKIT} batch ${TUMOR_BAM_DIR}/{}_T.analysis.ready.bam \
      --normal ${TUMOR_BAM_DIR}/{}_N.analysis.ready.bam \
      --targets ${BAIT_BED} \
      --annotate ${REF_FLAT} \
      --fasta ${GENOME_FASTA} \
      --access ${ACCESS_BED} \
      --output-reference ${PAIRED_RES_DIR}/{}.paired_mode_reference.cnn \
      --output-dir ${PAIRED_RES_DIR}/ \
      -p 2 \
      --diagram --scatter" ::: $(echo ${SEQC2_LIST} | awk '{print $1,$2}')

  ### CNVkit pooled noemal : PON
  POOLED_NORMAL=/storage/home/kangsm/analysis/cnv_test/cnvkit/pooled_normal

  ${RUN_CNVKIT} batch ${TUMOR_BAM_DIR}/HCC1395_T.analysis.ready.bam \
      --normal ${NORMAL_CELL_BAM_DIR}/*.analysis.ready.bam \
      --targets ${BAIT_BED} \
      --annotate ${REF_FLAT} \
      --fasta ${GENOME_FASTA} \
      --access ${ACCESS_BED} \
      --output-reference ${POOLED_NORMAL}/HCC1395/HCC1395.reference.cnn \
      --output-dir ${POOLED_NORMAL}/HCC1395 \
      -p 6 \
      --diagram --scatter

  ${RUN_CNVKIT} batch ${TUMOR_BAM_DIR}/HCC1395_FFPE_T.analysis.ready.bam \
      --normal ${NORMAL_FFPE_BAM_DIR}/*.analysis.ready.bam \
      --targets ${BAIT_BED} \
      --annotate ${REF_FLAT} \
      --fasta ${GENOME_FASTA} \
      --access ${ACCESS_BED} \
      --output-reference ${POOLED_NORMAL}/HCC1395_FFPE/HCC1395_FFPE.reference.cnn \
      --output-dir ${POOLED_NORMAL}/HCC1395_FFPE \
      -p 6 \
      --diagram --scatter

###------------------------------------------------------------------------------------------------#

### ADTEx ###--------------------------------------------------------------------------------------#
  for SEQ_ID in PON01_01 PON01_02 PON01_03 PON01_04 PON01_05
  do
    samtools view -@ 50 -b ${NORMAL_BAM_DIR}/${SEQ_ID}_N.analysis.ready.bam $(seq 1 22) X Y > ${CHR_BAM_DIR}/${SEQ_ID}_N.chr.bam
    samtools view -@ 50 -b ${TUMOR_BAM_DIR}/${SEQ_ID}_T.analysis.ready.bam $(seq 1 22) X Y > ${CHR_BAM_DIR}/${SEQ_ID}_T.chr.bam
    samtools view -H ${CHR_BAM_DIR}/${SEQ_ID}_N.chr.bam | grep -v '^@SQ\sSN:GL' | grep -v '^@SQ\sSN:NC' | grep -v '^@SQ\sSN:hs' \
      | grep -v '^@SQ\sSN:MT' | samtools reheader - ${CHR_BAM_DIR}/${SEQ_ID}_N.chr.bam > ${CHR_BAM_DIR}/${SEQ_ID}_N.chr.rev.bam
    samtools view -H ${CHR_BAM_DIR}/${SEQ_ID}_T.chr.bam | grep -v '^@SQ\sSN:GL' | grep -v '^@SQ\sSN:NC' | grep -v '^@SQ\sSN:hs' \
      | grep -v '^@SQ\sSN:MT' | samtools reheader - ${CHR_BAM_DIR}/${SEQ_ID}_T.chr.bam > ${CHR_BAM_DIR}/${SEQ_ID}_T.chr.rev.bam
  done

  RUN_ADTEX="singularity exec -B /storage,/data /storage/images/adtex.sif ADTEx.py"
  CENTROMERE=/storage/home/kangsm/myDB/bed/ADTEx_hg19/GRCh37.centromeres.bed
  AGILENT_BED=/storage/home/kangsm/myDB/bed/ADTEx_hg19/ADTEx_SureSelect.HumanAllExon.v6.extended.target.bed
  CHR_BAM_DIR=/storage/home/kangsm/analysis/cnv_test/chr_bams

  ADTEX_RES_DIR=/storage/home/kangsm/analysis/cnv_test/adtex/

  
  samtools view -@ 50 -b ${TUMOR_BAM_DIR}/${SEQ_ID}_N.analysis.ready.bam $(seq 1 22) X Y > ${TUMOR_BAM_DIR}/${SEQ_ID}_N.chr.bam

    ${RUN_ADTEX} \
        -n ${TUMOR_BAM_DIR}/HCC1395_N.chr.only.bam \
        -t ${TUMOR_BAM_DIR}/HCC1395_T.chr.only.bam \
        -s ${SEQ_ID} \
        -b ${AGILENT_BED} \
        -c ${CENTROMERE} \
        -o ${ADTEX_RES_DIR}/${SEQ_ID}

    ${RUN_ADTEX} \
        -n ${TUMOR_BAM_DIR}/${SEQ_ID}_N.chr.rev.bam \
        -t ${TUMOR_BAM_DIR}/${SEQ_ID}_T.chr.rev.bam \
        -s ${SEQ_ID} \
        -b ${AGILENT_BED} \
        -c ${CENTROMERE} \
        -o ${ADTEX_RES_DIR}/${SEQ_ID}


  for SEQ_ID in $(echo ${PON_SAMPLE_LIST} | awk '{print $1,$2,$3,$4,$5}')
  do
    ${RUN_ADTEX} \
        -n ${CHR_BAM_DIR}/${SEQ_ID}_N.chr.rev.bam \
        -t ${CHR_BAM_DIR}/${SEQ_ID}_T.chr.rev.bam \
        -s ${SEQ_ID} \
        -b ${AGILENT_BED} \
        -c ${CENTROMERE} \
        -o ${ADTEX_RES_DIR}/${SEQ_ID}
  done

  parallel -j 2 -k "${RUN_ADTEX} \
        -n ${CHR_BAM_DIR}/{}_N.chr.rev.bam \
        -t ${CHR_BAM_DIR}/{}_T.chr.rev.bam \
        -s {} \
        -b ${AGILENT_BED} \
        -c ${CENTROMERE} \
        -o ${ADTEX_RES_DIR}/{}" ::: $(echo ${PON_SAMPLE_LIST} | awk '{print $1,$2,$3,$4,$5}')

###------------------------------------------------------------------------------------------------#



