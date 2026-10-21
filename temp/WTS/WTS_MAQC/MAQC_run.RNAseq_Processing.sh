
BASE_DIR="/data/wts"
SEQ_FOLDER="WTS_MAQC"

for SEQ_ID in "MAQC_C_F1_rep2" 
do 
    #---| FASTQC 
        /storage/home/kangsm/runScripts/WTS_run.FastqQC.sh \
            --baseDir ${BASE_DIR} \
            --seqFolder ${SEQ_FOLDER} \
            --seqID ${SEQ_ID} \
            --readLength 100 \
            --threads 15 \
            --runMultiqc "false" \
            --seqType "PE"

    #---| FASTQ to BAM 
        /storage/home/kangsm/runScripts/WTS_run.FastqToBAM.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS 15 \
            --READ_LENGTH 100 \
            --STRANDED "false" \
            --GENOME_ASSEMBLY "hg19" \
            --PAIRED "true"

    #---| RSEM QUANTIFICATION
        /storage/home/kangsm/runScripts/WTS_run.Quantification_RSEM.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS 15 \
            --READ_LENGTH 100 \
            --STRANDED "false" \
            --GENOME_ASSEMBLY "hg19" \
            --PAIRED "true"

    #---| BAM QC 
        /storage/home/kangsm/runScripts/WTS_run.BAM_QC.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS 15 \
            --STRANDED "false" \
            --GENOME_ASSEMBLY "hg19" \
            --PAIRED "true" \
            --RUN_MULTIQC "true" \
            --QC_RES_DB_IMPORT "true"
        
    #---| FUSION-GENES
        /storage/home/kangsm/runScripts/WTS_run.FusionGenes_and_SpliceVariants.sh \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --THREADS 15 \
            --GENOME_ASSEMBLY "hg19" \
            --RUN_ARRIBA "false" \
            --RUN_STARFUSION "false" \
            --RUN_SPLICEVARIANT "false" \
            --DRAW_ARRIBA_PLOTS "true" \
            --PAIRED_END "false"
    
        Rscript /storage/home/kangsm/runScripts/WTS_run.FusionGenes_SpliceVariants_Summary.R \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --DRAW_STAR_FUSION_CIRCOS TRUE \
            --GENOME_ASSEMBLY "hg19" \
            --DB_IMPORT TRUE 

    #---| SAVE EXPRESSION MATRIX
        Rscript /storage/home/kangsm/runScripts/WTS_run.CreateExprValuesObjects.R \
            --BASE_DIR ${BASE_DIR} \
            --SEQ_FOLDER ${SEQ_FOLDER} \
            --SEQ_ID ${SEQ_ID} \
            --ANALYSIS_CONFIG /data/wts/params/default_wts_analysis_config.yaml \
            --DB_IMPORT TRUE
done

for SEQ_ID in "MAQC_A_rep1" "MAQC_A_rep2" "MAQC_B_rep1" "MAQC_B_rep2" "MAQC_C_F1_rep1" "MAQC_C_F1_rep2" "MAQC_C_F2_rep1" "MAQC_C_F2_rep2"
do 
    Rscript /storage/home/kangsm/runScripts/WTS_run.CreateExprValuesObjects.R \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --ANALYSIS_CONFIG /data/wts/params/default_wts_analysis_config.yaml \
        --DB_IMPORT TRUE
done

for SEQ_ID in "MAQC_C_F1_rep1" "MAQC_C_F1_rep2" "MAQC_C_F2_rep1" "MAQC_C_F2_rep2"
do 
    /storage/home/kangsm/runScripts/WTS_run.FusionGenes_and_SpliceVariants.sh \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --THREADS 15 \
        --GENOME_ASSEMBLY "hg19" \
        --RUN_ARRIBA "true" \
        --RUN_STARFUSION "true" \
        --RUN_SPLICEVARIANT "true" \
        --DRAW_ARRIBA_PLOTS "true" \
        --PAIRED_END "true"

    Rscript /storage/home/kangsm/runScripts/WTS_run.FusionGenes_SpliceVariants_Summary.R \
        --BASE_DIR ${BASE_DIR} \
        --SEQ_FOLDER ${SEQ_FOLDER} \
        --SEQ_ID ${SEQ_ID} \
        --DRAW_STAR_FUSION_CIRCOS TRUE \
        --GENOME_ASSEMBLY "hg19" \
        --DB_IMPORT TRUE 
done

"MAQC_C_F1_rep1" "MAQC_C_F1_rep2" "MAQC_C_F2_rep1" "MAQC_C_F2_rep2"


#--------------------------------------------------------------------------------------------------#

    Rscript /storage/home/kangsm/runScripts/WTS_run.WTS_AnalysisRun.R \
        --BASE_DIR /data/wts \
        --SEQ_FOLDER WTS_MAQC \
        --ANALYSIS_CONFIG /data/wts/WTS_MAQC/meta/WTS_MAQC_wts_analysis_config.yaml \
        --GENOME_ASSEMBLY hg19

    


