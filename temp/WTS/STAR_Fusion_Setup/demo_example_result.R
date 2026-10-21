

# https://github.com/STAR-Fusion/STAR-Fusion_benchmarking_data/blob/master/cancer_cell_lines/SuppTable-cancer_cell_lines.csv
# https://data.broadinstitute.org/Trinity/STAR_FUSION_PAPER/SupplementaryData/cancer_cell_lines/FASTQ/


    library(rtracklayer)
    library(dplyr)
    library(chimeraviz)

    source("/storage/home/kangsm/myScripts/R_FUNCTIONAL_MODULES/FunMod_StarFusionResults.R")

    DEMO_DIR     = "/storage/home/kangsm/analysis/star-fusion/tso500_fastq_run"
    GENOME       = "hg19"

    FASTQ_INPUT_RES    = sprintf("%s/23Cancer01_09/star-fusion.fusion_predictions.abridged.tsv", DEMO_DIR)
    #JUNCTION_INPUT_RES = sprintf("%s/junction_input_result/star-fusion.fusion_predictions.abridged.tsv", DEMO_DIR)

    fusionFASTQ    = import_starfusion(filename=FASTQ_INPUT_RES, genome_version=GENOME)
    #fusionJUNCTION = import_starfusion(filename=JUNCTION_INPUT_RES, genome_version=GENOME)

    resFASTQ = StarFusionList.Summary(fusionFASTQ, COMPACT=T)
    #resJUNCTION = StarFusionList.Summary(fusionJUNCTION, COMPACT=T)

    openxlsx::write.xlsx(resFASTQ, "23Cancer01_09.star_fusion.fastq.input.test.run.result_raw.xlsx", rowNames=F, overwrite=T)

    fusionFASTQ_filtered = Filter.STAR.Fusion.Result( STAR_FUSION_LIST=fusionFASTQ, FFPM_CUTOFF=1)

    plot_circle(fusionFASTQ_filtered)

    plot_circle(fusionFASTQ)

    
    for( ID in c(5,8,9,11))
    {
        FUSION = get_fusion_by_id(fusionFASTQ_filtered, 5)
        png(sprintf("23Cancer01_10/FUSION_ID_%s_plot_fusion.png", ID), width=24, height=2, units="in", res=200, type="cairo")
        plot_fusion( fusion = FUSION, edb = ENSDB_hg19, reduce_transcripts=T, non_ucsc = TRUE)
        dev.off()
        png(sprintf("23Cancer01_10/FUSION_ID_%s_plot_transcript.png", ID), width=24, height=2, units="in", res=200, type="cairo")
        plot_transcripts( fusion = FUSION, edb = ENSDB_hg19, reduce_transcripts=T )
        dev.off()
    }



    FUSION = get_fusion_by_id(fusionFASTQ_filtered, 5)
    #BAM_FILE ="fastq_input_result/Aligned.out.sorted.bam"
    #FUSION = add_fusion_reads_alignment(FUSION, "fastq_input_result/Aligned.out.sorted.bam")
    

    # 1) plot_fusion
    plot_fusion( fusion = FUSION, edb = ENSDB_hg19, reduce_transcripts=T, non_ucsc = TRUE)
    plot_fusion( fusion = FUSION, edb = ENSDB_hg19, non_ucsc = TRUE)

    # 2) plot_fusion_reads
    # plot_fusion_reads(FUSION)
    
    # 3) plot_fusion_transcripts
    plot_fusion_transcript( FUSION, ENSDB_hg19 )

    # 4) plot_fusion_transcripts_graph
    plot_fusion_transcripts_graph( fusion = FUSION, edb = ENSDB_hg19 )

    # 5) plot_transcripts
    plot_transcripts( fusion = FUSION, edb = ENSDB_hg19, reduce_transcripts=T )
    plot_transcripts( fusion = FUSION, edb = ENSDB_hg19 )


    # 
    Calculate.Precision.Recall( TP = 15, TN = 0, FP = 56, FN = 0 )

    Calculate.Precision.Recall( TP = 15, TN = 0, FP = 0, FN = 0 )

    #
    splicing_res_dir="/storage/home/kangsm/analysis/star-fusion/tso500_fastq_run/23Cancer01_09/ctat_splicing"

    cancerIntron = read.delim(sprintf("%s/23Cancer01_09.cancer.introns",splicing_res_dir))

    CI = cancerIntron %>% dplyr::select(c(1,2,3,4,5,8)) 


