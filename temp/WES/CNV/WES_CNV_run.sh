
source /storage/apps/miniconda3/bin/activate cnv

#---| REQUIRED FILES based on GENOME VER. |----------------------------------------------------------------------------------------------------------#
CNV_FILES_DIR=/storage/home/kangsm/myDB/cnv
REFERENCE_DIR=/storage/references_and_index

if [ "${REF_VERSION}" = "hg38" ]; then
    GENOME_FASTA=${REFERENCE_DIR}/hg38/fasta/Homo_sapiens_assembly38.fasta
    REF_FLAT=/storage/home/kangsm/myDB/annotation/hg38/hg38.refFlat.txt
    CNVKIT_ACCESS_BED=${CNV_FILES_DIR}/hg38_access.10k_CNVkit.bed
    CENTROMERE_BED=${CNV_FILES_DIR}/hg38_centromere_ADTEx.bed
else
    GENOME_FASTA=${REFERENCE_DIR}/hg19/fasta/human_g1k_v37_decoy.fasta
    REF_FLAT=${REFERENCE_DIR}/hg19/extra/human_g1k_v37.refFlat.txt
    CNVKIT_ACCESS_BED=${CNV_FILES_DIR}/GRCh37_access.10k_CNVkit.bed
    CENTROMERE_BED=${CNV_FILES_DIR}/GRCh37_centromeres_ADTEx.bed
fi

RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
CNVKIT_RUN_CMD="${RUN_SINGULARITY}/cnvkit-0.9.10.sif cnvkit.py"
ADTEX_RUN_CMD="${RUN_SINGULARITY}/adtex.sif ADTEx.py"
GATK3="${RUN_SINGULARITY}/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar"

#--------------------------------------------------------------------------------------------------#
INPUT_TUMOR_BAM=
INPUT_NORMAL_BAM=
SAMPLE_ID=
TARGET_BED=
PROBE_BED=
#

if [ "${RUN_MODE}" = "paired" ]; then

    if [ "${RUN_CNVKIT}" == "true" ]; then

        CNVKIT_RES_DIR=${CNV_RES_DIR}/CNVkit_paired
        if [ ! -d ${CNVKIT_RES_DIR} ]; then mkdir -p ${CNVKIT_RES_DIR}; fi

        ${CNVKIT_RUN_CMD} batch ${INPUT_TUMOR_BAM} \
            --normal ${INPUT_NORMAL_BAM} \
            --targets ${PROBE_BED} \
            --fasta ${GENOME_FASTA} \
            --access ${CNVKIT_ACCESS_BED} \
            --output-reference ${CNVKIT_RES_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_reference.cnn \
            --output-dir ${CNVKIT_RES_DIR}/${SAMPLE_ID} \
            -p 2 \
            --diagram --scatter

    fi


    if [ "${RUN_ADTEX}" = "true" ]; then
        
            





else





fi








#---| Required Files < hg38 > |--------------------------------------------------------------------#
HG38_GENOME_FASTA=/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta
HG38_REF_FLAT=/storage/home/kangsm/myDB/annotation/hg38/hg38.refFlat.txt
HG38_AGILENT_BAIT_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.probe.bed
HG38_AGILENT_TARGET_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.target.bed
HG38_AGILENT_TARGET_INTERVAL=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.extended.target.interval_list

HG38_U_ACCESS_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_access.10k_CNVkit.bed
HG38_CENTROMERE_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_centromere_ADTEx.bed

#---| Run SWTOOLS |--------------------------------------------------------------------------------#
RUN_CNVKIT="singularity exec -B /storage,/data /storage/images/cnvkit-0.9.10.sif cnvkit.py"
RUN_ADTEX="singularity exec -B /storage,/data /storage/images/adtex.sif ADTEx.py"


#---| DATA_FOLDER |--------------------------------------------------------------------------------#
SEQC_PAIRED_TEST_SAMPLE=/storage/home/kangsm/analysis/cnv_test/testset2
SEQC_PON_TEST_CELL=/storage/home/kangsm/analysis/cnv_test/testset2_normal_cell
SEQC_PON_TEST_FFPE=/storage/home/kangsm/analysis/cnv_test/testset2_normal_ffpe

CNVKIT_PAIRED_RESULT=${CNVKIT_WORK_DIR}/paired
CNVKIT_PON_RESULT=${CNVKIT_WORK_DIR}/pooled_normal

###| CNVkit |###

CNVKIT_WORK_DIR=/storage/home/kangsm/analysis/cnv_test/cnvkit
CNVKIT_PAIRED_RESULT=${CNVKIT_WORK_DIR}/paired
CNVKIT_PON_RESULT=${CNVKIT_WORK_DIR}/pooled_normal

CNVKIT_RUN_MODE=
INPUT_TUMOR=
INPUT_NORMAL=
PON_BAM_DIR=

if [ "${CNVKIT_RUN_MODE}" = "paired" ]; then
    CNVkitResultFolder=${CNVKIT_PAIRED_RESULT}
    
    ${RUN_CNVKIT} batch ${INPUT_TUMOR} \
        --normal ${INPUT_NORMAL} \
        --targets ${CNVKIT_BED} \



else
    CNVkitResultFolder=${CNVKIT_PON_RESULT}


fi


parallel -j 2 -k "${RUN_CNVKIT} batch ${SEQC_PAIRED_TEST_SAMPLE}/{}_T.analysis.ready.bam \
    --normal ${SEQC_PAIRED_TEST_SAMPLE}/{}_N.analysis.ready.bam \
    --targets ${HG38_AGILENT_BAIT_BED} \
    --annotate ${HG38_REF_FLAT} \
    --fasta ${HG38_GENOME_FASTA} \
    --access ${HG38_CNVKIT_ACCESS_BED} \
    --output-reference ${CNVKIT_PAIRED_RESULT}/{}.paired_mode_reference.cnn \
    --output-dir ${CNVKIT_PAIRED_RESULT}/ \
    -p 2 \
    --diagram --scatter" ::: $(echo ${SEQC2_LIST} | awk '{print $1,$2}')

 ${RUN_CNVKIT} batch ${SEQC_PAIRED_TEST_SAMPLE}/HCC1395_T.analysis.ready.bam \
    --normal ${SEQC_PON_TEST_CELL}/HCC1395_N*.analysis.bam \
    --targets ${HG38_AGILENT_BAIT_BED} \
    --annotate ${HG38_REF_FLAT} \
    --fasta ${HG38_GENOME_FASTA} \
    --access ${HG38_CNVKIT_ACCESS_BED} \
    --output-reference ${CNVKIT_PON_RESULT}/HCC1395/HCC1395.reference.cnn \
    --output-dir ${CNVKIT_PON_RESULT}/HCC1395 \
    -p 6 \
    --diagram --scatter

 ${RUN_CNVKIT} batch ${SEQC_PAIRED_TEST_SAMPLE}/HCC1395_FFPE_T.analysis.ready.bam \
    --normal ${SEQC_PON_TEST_FFPE}/HCC1395_FFPE_N*.analysis.bam \
    --targets ${HG38_AGILENT_BAIT_BED} \
    --annotate ${HG38_REF_FLAT} \
    --fasta ${HG38_GENOME_FASTA} \
    --access ${HG38_CNVKIT_ACCESS_BED} \
    --output-reference ${CNVKIT_PON_RESULT}/HCC1395_FFPE/HCC1395_FFPE.reference.cnn \
    --output-dir ${CNVKIT_PON_RESULT}/HCC1395_FFPE \
    -p 6 \
    --diagram --scatter



### ADTEx
CHR_CONTG_FILE=/storage/home/kangsm/myDB/annotation/hg38/hg38_chr.contigs
HG38_CENTROMERE=/storage/home/kangsm/myDB/bed/hg38/hg38_centromere_ADTEx.bed
HG38_AGILENT_TARGET_BED=/storage/home/kangsm/myDB/bed/hg38/hg38_sureselect.v6.utr.extended.target.bed

RUN_ADTEX="singularity exec -B /storage,/data /storage/images/adtex.sif ADTEx.py"

for INPUT_ID in HCC1395_T HCC1395_N HCC1395_FFPE_T HCC1395_FFPE_N
do
    #INPUT_BAMFILE=${TUMOR_BAM_DIR}/${INPUT_ID}.analysis.ready.bam
    OUTPUT_CHR_BAM=${TUMOR_BAM_DIR}/${INPUT_ID}.chr.only.bam
    OUTPUT_REHEAD_BAM=${TUMOR_BAM_DIR}/${INPUT_ID}.chr.rehead.bam
    #samtools view -@ 50 -b ${INPUT_BAMFILE} $(cat ${CHR_CONTG_FILE}) > ${OUTPUT_CHR_BAM}
    samtools view -H ${OUTPUT_CHR_BAM} | grep -v '^@SQ\sSN:chrUn' | grep -v '^@SQ\sSN:HLA' | grep -v '^@SQ\sSN:chr[1-9,X,Y]_' \
      | grep -v '^@SQ\sSN:chr[1-9][0-9]_' | grep -v '^@SQ\sSN:chrEBV' | samtools reheader - ${OUTPUT_CHR_BAM} > ${OUTPUT_REHEAD_BAM}
done

    ${RUN_ADTEX} --keeptemp \
        -n ${TUMOR_BAM_DIR}/HCC1395_FFPE_N.chr.rehead.bam \
        -t ${TUMOR_BAM_DIR}/HCC1395_FFPE_T.chr.rehead.bam \
        -s "HCC1395_FFPE" \
        -b ${HG38_AGILENT_TARGET_BED} \
        -c ${HG38_CENTROMERE} \
        -o /storage/home/kangsm/analysis/cnv_test/adtex/HCC1395_FFPE


    ${RUN_ADTEX} --keeptemp \
        -n ${TUMOR_BAM_DIR}/HCC1395_N.chr.rehead.bam \
        -t ${TUMOR_BAM_DIR}/HCC1395_T.chr.rehead.bam \
        -s "HCC1395" \
        -b ${HG38_AGILENT_TARGET_BED} \
        -c ${HG38_CENTROMERE} \
        -o /storage/home/kangsm/analysis/cnv_test/adtex/HCC1395


