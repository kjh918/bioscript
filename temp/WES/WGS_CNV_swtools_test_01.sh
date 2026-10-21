
#! /bin/bash

GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta
REF_FLAT=/storage/references_and_index/hg19/extra/human_g1k_v37.refFlat.txt

TWIST_BAIT_BED=/storage/references_and_index/hg19/bed/Twist_WES_2.0/probe.bed
TWIST_TARGET_BED=/storage/references_and_index/hg19/bed/Twist_WES_2.0/target.bed
TWIST_TARGET_INTERVAL=/storage/references_and_index/hg19/bed/Twist_WES_2.0/target.interval_list

SURESELECT_PROBE_BED=/storage/home/kangsm/myDB/bed/Agilent_SureSelect_v6_hg19/agilent_SureSelect_HumanAllExon_v6.probe.bed
SURESELECT_TARGET_BED=/storage/home/kangsm/myDB/bed/Agilent_SureSelect_v6_hg19/agilent_SureSelect_HumanAllExon_v6.target.extended.bed
SURESELECT_TARGET_INTERVAL=/storage/home/kangsm/myDB/bed/Agilent_SureSelect_v6_hg19/agilent_SureSelect_HumanAllExon_v6.target.extended.interval_list

NORMAL_BAM_DIR=/storage/home/kangsm/analysis/cnv_test/normal_bams
TUMOR_BAM_DIR=/storage/home/kangsm/analysis/cnv_test/tumor_bams
VCF_DIR=/storage/home/kangsm/analysis/cnv_test/vcfs

PON_SAMPLE_LIST="PON01_01 PON01_02 PON01_03 PON01_04 PON01_05"
NA12878_LIST="NA12878_01 NA12878_02"

### CNVpytor ###-----------------------------------------------------------------------------------#
  source /storage/apps/miniconda3/bin/activate cnv

  RES_DIR=/storage/home/kangsm/analysis/cnv_test/cnvpytor
  cd ${RES_DIR} # WORK_DIR

  # TUMOR BAM CNV
  for SEQ_ID in NA12878_01 NA12878_02
  do
    OUT_PYTOR=${RES_DIR}/${SEQ_ID}.pytor
    cnvpytor -root ${OUT_PYTOR} -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y \
      -rd ${TUMOR_BAM_DIR}/${SEQ_ID}.analysis.ready.bam

    cnvpytor -root ${OUT_PYTOR} -his 1000 10000 100000
    cnvpytor -root ${OUT_PYTOR} -partition 1000 10000 100000
    cnvpytor -root ${OUT_PYTOR} -call 1000 > ${SEQ_ID}.calls.1000.tsv
    cnvpytor -root ${OUT_PYTOR} -call 10000 > ${SEQ_ID}.calls.10000.tsv
    cnvpytor -root ${OUT_PYTOR} -call 100000 > ${SEQ_ID}.calls.100000.tsv
  done

  for SEQ_ID in PON01_01 PON01_02 PON01_03 PON01_04 PON01_05
  do
    OUT_PYTOR=${RES_DIR}/${SEQ_ID}_T.pytor
    cnvpytor -root ${OUT_PYTOR} -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y \
      -rd ${TUMOR_BAM_DIR}/${SEQ_ID}_T.analysis.ready.bam

    cnvpytor -root ${OUT_PYTOR} -his 1000 10000 100000
    cnvpytor -root ${OUT_PYTOR} -partition 1000 10000 100000
    cnvpytor -root ${OUT_PYTOR} -call 1000 > ${SEQ_ID}.calls.1000.tsv
    cnvpytor -root ${OUT_PYTOR} -call 10000 > ${SEQ_ID}_T.calls.10000.tsv
    cnvpytor -root ${OUT_PYTOR} -call 100000 > ${SEQ_ID}_T.calls.100000.tsv
  done

###------------------------------------------------------------------------------------------------#

### DELLY ###--------------------------------------------------------------------------------------#
  export OMP_NUM_THREADS=2
  RUN_DELLY="singularity exec -B /storage,/data /storage/images/delly_v1.2.6.sif delly"
  EXCLURE_REGION_FILE=/storage/home/kangsm/myDB/bed/Delly_excludeRegions/human.grch37.excl.tsv

  DELLY_RES_DIR=/storage/home/kangsm/analysis/cnv_test/delly
  INPUT_LIST=$(echo "NA12878_01 NA12878_02 PON01_01_T PON01_02_T PON01_03_T PON01_04_T PON01_05_T" | awk '{print $1,$2,$3,$4,$5,$6,$7}')

  parallel -j 7 -k "${RUN_DELLY} call \
      --exclude ${EXCLURE_REGION_FILE} \
      --outfile ${DELLY_RES_DIR}/{}.delly.default.run.sv.bcf \
      --genome ${GENOME_FASTA} \
      ${TUMOR_BAM_DIR}/{}.analysis.ready.bam" ::: $INPUT_LIST


  SEQ_ID="PON01_01"

  for SEQ_ID in PON01_02 PON01_03 PON01_04 PON01_05
  do
    TUMOR_BAM_FILE=${TUMOR_BAM_DIR}/${SEQ_ID}_T.analysis.ready.bam
    NORMAL_BAM_FILE=${NORMAL_BAM_DIR}/${SEQ_ID}_N.analysis.ready.bam

      ${RUN_DELLY} call \
          --exclude ${EXCLURE_REGION_FILE} \
          --outfile ${DELLY_RES_DIR}/${SEQ_ID}.NT.paired.sv.bcf \
          --genome ${GENOME_FASTA} \
          ${TUMOR_BAM_FILE} ${NORMAL_BAM_FILE}

      PAT_NO=$(echo $SEQ_ID | awk -F "_" '{print $2}')
      echo -e "coad.tumor.$PAT_NO\ttumor" > ${DELLY_RES_DIR}/sample_desc.tsv
      echo -e "coad.normal.$PAT_NO\tcontrol" >> ${DELLY_RES_DIR}/sample_desc.tsv

      ${RUN_DELLY} filter \
          -f somatic \
          -o ${DELLY_RES_DIR}/${SEQ_ID}.NT.paired.sv.filtered.bcf \
          -s ${DELLY_RES_DIR}/sample_desc.tsv \
          ${DELLY_RES_DIR}/${SEQ_ID}.NT.paired.sv.bcf

      ${RUN_DELLY} call \
          -g ${GENOME_FASTA} \
          -v ${DELLY_RES_DIR}/${SEQ_ID}.NT.paired.sv.filtered.bcf \
          -o ${DELLY_RES_DIR}/${SEQ_ID}.NT.paired.sv.genotype.bcf \
          -x ${EXCLURE_REGION_FILE} \
          ${TUMOR_BAM_FILE} ${NORMAL_BAM_FILE}
      
      ${RUN_DELLY} filter \
          -f somatic \
          -o ${DELLY_RES_DIR}/${SEQ_ID}.NT.paired.somatic.bcf \
          -s ${DELLY_RES_DIR}/sample_desc.tsv \
          ${DELLY_RES_DIR}/${SEQ_ID}.NT.paired.sv.genotype.bcf
  done

  # DELLY SCNA
    GRCH37_MAP=/storage/home/kangsm/myDB/fasta/delly_mappability_maps/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz
    SCNA_RES_DIR=/storage/home/kangsm/analysis/cnv_test/delly/SCNA
    cd ${SCNA_RES_DIR} 

  SEQ_ID="PON01_01"

  TUMOR_BAM_FILE=${TUMOR_BAM_DIR}/${SEQ_ID}_T.analysis.ready.bam
  NORMAL_BAM_FILE=${NORMAL_BAM_DIR}/${SEQ_ID}_N.analysis.ready.bam
  OUTPUT_DIR=${SCNA_RES_DIR}/${SEQ_ID}
  mkdir -p ${OUTPUT_DIR}
  cd ${OUTPUT_DIR}

  ${RUN_DELLY} cnv -i 5000000 -r ${TARGET_BED} -b ${TARGET_BED} \
      -g ${GENOME_FASTA} -m ${GRCH37_MAP} ${TUMOR_BAM_FILE}
###
