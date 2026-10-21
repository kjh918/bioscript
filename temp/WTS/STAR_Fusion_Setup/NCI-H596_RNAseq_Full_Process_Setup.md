

# FusionGene Detection Full Process Setup with NCI-H596 RNAseq Data

## 1. DATA PRPARE --------------------------------------------------------------
- RNAseq Fastq Symbolic-link 생성

<pre  style="font-family: Consolas">
FASTQ_DIR=/storage/home/kangsm/analysis/gcx.wts.pipeline.fusionGenesDetection/NCI-H596_RNAseq_test_run/fastq

ln -s /storage/DB/ccle/wts_fastq/NCI-H596_C1EHH_R1.fastq.gz NCI-H596_rep1.trimmed_R1.fastq.gz
ln -s /storage/DB/ccle/wts_fastq/NCI-H596_C1EHH_R2.fastq.gz NCI-H596_rep1.trimmed_R2.fastq.gz
ln -s /storage/DB/ccle/wts_fastq/NCI-H596_D1KNF_R1.fastq.gz NCI-H596_rep2.trimmed_R1.fastq.gz
ln -s /storage/DB/ccle/wts_fastq/NCI-H596_D1KNF_R2.fastq.gz NCI-H596_rep2.trimmed_R2.fastq.gz
ln -s /storage/DB/ccle/wts_fastq/NCI-H596_D1KYP_R1.fastq.gz NCI-H596_rep3.trimmed_R1.fastq.gz
ln -s /storage/DB/ccle/wts_fastq/NCI-H596_D1KYP_R2.fastq.gz NCI-H596_rep3.trimmed_R2.fastq.gz
</pre>

## 2. Default Parameters -------------------------------------------------------

<pre style="font-family: Consolas">
    BASE_DIR=/storage/home/kangsm/analysis/gcx.wts.pipeline.fusionGenesDetection
    SEQ_FOLDER=NCI-H596_RNAseq_test_run
    THREADS=15
    READ_LENGTH="100"
    STRANDED="false"
    GENOME_ASSEMBLY="hg19
</pre>


## 3. STAR Alignment -----------------------------------------------------------
- STAR-Alignment Run 
- STAR option = <b>Arriba compartible</b>

<pre  style="font-family: Consolas">
    SEQ_ID=NCI-H596_rep1
    ##
    STAR_INDEX=/storage/references_and_index/${GENOME_ASSEMBLY}/rsem_index/pe${READ_LENGTH}
    GTF_FILE=/storage/references_and_index/${GENOME_ASSEMBLY}/gtf/${GENOME_ASSEMBLY}.gtf
    ##
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    RUN_STAR="${RUN_SINGULARITY}/star-2.7.11.sif STAR"
    RUN_SAMTOOLS="${RUN_SINGULARITY}/samtools-1.18.sif samtools"
    ##
    /storage/home/kangsm/runScripts/WTS_run.FastqToBAM.sh \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --THREADS 15 \
        --READ_LENGTH 100 \
        --STRANDED "false" \
        --GENOME_ASSEMBLY "hg19"
</pre>

## Individual Tool Test --------------------------------------------------------

1) Arriba
<pre style="font-family: Consolas">
    GENOME_FASTA=/storage/references_and_index/${GENOME_ASSEMBLY}/fasta/human_g1k_v37_decoy.fasta
    GTF_FILE=/storage/references_and_index/${GENOME_ASSEMBLY}/gtf/${GENOME_ASSEMBLY}.gtf

    ARRIBA_REFERENCE=/storage/references_and_index/${GENOME_ASSEMBLY}/arriba
    ARRIBA_BLACKLIST=${ARRIBA_REFERENCE}/blacklist_${GENOME_ASSEMBLY}.tsv.gz
    ARRIBA_KNOWN_FUSIONS=${ARRIBA_REFERENCE}/known_fusions_${GENOME_ASSEMBLY}.tsv.gz
    ARRIBA_PROTEIN_FEATURES=${ARRIBA_REFERENCE}/protein_domains_${GENOME_ASSEMBLY}.gff3
    ARRIBA_CYTOBAND=${ARRIBA_REFERENCE}/cytobands_${GENOME_ASSEMBLY}.tsv

    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    ARRIBA_RUN="${RUN_SINGULARITY}/arriba-2.4.0.sif /arriba_v2.4.0"

    ${ARRIBA_RUN}/arriba \
        -x ${BAM_DIR}/${SEQ_ID}.Aligned.out.Sorted.bam \
        -o ${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.fusion.genes.tsv \
        -O ${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.filtered.out.fusion.genes.tsv \
        -a ${GENOME_FASTA} \
        -g ${GTF_FILE} \
        -b ${ARRIBA_BLACKLIST} \
        -k ${ARRIBA_KNOWN_FUSIONS} \
        -p ${ARRIBA_PROTEIN_FEATURES} \
        -i "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

    ${RUN_ARRIBA}/draw_fusions.R \
        --fusions=${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.fusion.genes.gencode.tsv \
        --alignments=${BASE_DIR}/${SEQ_FOLDER}/${SEQ_ID}/bam/${SEQ_ID}.Aligned.out.Sorted.bam \
        --output=${ARRIBA_RES_DIR}/${SEQ_ID}.arriba.fusion.genes.gencode.plots.pdf \
        --annotation=${ARRIBA_GTF_GENCODE} \
        --cytobands=${ARRIBA_CYTOBAND} \
        --proteinDomains=${ARRIBA_PROTEIN_FEATURES}
</pre>

2) STAR-Fusion
<pre style="font-family: Consolas">
    RUN_SINGULARITY="singularity exec -B /storage,/data /storage/images"
    STAR_FUSION_RUN="${RUN_SINGULARITY}/star-fusion-0.13.0.sif STAR-Fusion"
    CTAT_RESOURCE_LIB=/storage/references_and_index/hg19/star_fusion/CTAT_resource_lib

    ${STAR_FUSION_RUN} \
        --genome_lib_dir ${CTAT_RESOURCE_LIB} \
        --left_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed_R1.fastq.gz \
        --right_fq ${FASTQ_DIR}/${SEQ_ID}.trimmed_R2.fastq.gz \
        --output_dir ${STARFUSION_RES_DIR} \
        --CPU ${THREADS} \
        --STAR_twopass  
</pre>

3) CTAT_SPLICING
<pre>
    CTAT_SPLICING_RUN="${RUN_SINGULARITY}/ctat-splicing-0.0.2.sif /usr/local/src/CTAT-SPLICING"
    RUN_SAMTOOLS="${RUN_SINGULARITY}/samtools-1.18.sif samtools"

    ${RUN_SAMTOOLS} sort --threads ${THREADS} --output-fmt BAM -o ${STARFUSION_RES_DIR}/Aligned.out.Sorted.bam ${STARFUSION_RES_DIR}/Aligned.out.bam  
    
    ${RUN_SAMTOOLS} index --threads ${THREADS} -b ${STARFUSION_RES_DIR}/Aligned.out.Sorted.bam

    ${CTAT_SPLICING_RUN}/STAR_to_cancer_introns.py \
         --ctat_genome_lib ${CTAT_RESOURCE_LIB} \
         --SJ_tab_file  ${STARFUSION_RES_DIR}/SJ.out.tab \
         --chimJ_file ${STARFUSION_RES_DIR}/Chimeric.out.junction \
         --vis \
         --bam_file ${STARFUSION_RES_DIR}/Aligned.out.Sorted.bam \
         --output_prefix ${SPLICE_VAR_RES_DIR}/${SEQ_ID} \
         --sample_name ${SEQ_ID} 
</pre>

## SCRIPT PIPELINE TEST --------------------------------------------------------
<pre>

    BASE_DIR="/storage/home/kangsm/analysis/gcx.wts.pipeline.fusionGenesDetection"
    SEQ_FOLDER="NCI-H596_RNAseq_test_run"
    SEQ_ID="NCI-H596_rep3"

    /storage/home/kangsm/runScripts/WTS_run.FastqToBAM.sh \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --THREADS 15 \
        --READ_LENGTH 100 \
        --STRANDED "false" \
        --GENOME_ASSEMBLY "hg19"

    /storage/home/kangsm/runScripts/WTS_run.FusionGenes_and_SpliceVariants.sh \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --THREADS 15 \
        --GENOME_ASSEMBLY "hg19" \
        --RUN_ARRIBA "false" \
        --RUN_STARFUSION "true" \
        --RUN_SPLICEVARIANT "false" \
        --DRAW_ARRIBA_PLOTS "true"

    /storage/home/kangsm/runScripts/WTS_run.AlignedBAM_QC.sh \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --THREADS 15 \
        --STRANDED "false" \
        --GENOME_ASSEMBLY "hg19"
</pre>


