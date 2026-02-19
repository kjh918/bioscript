
# RNAseq BAM-QC RESULT INTEGRATION AND SAVE SCRIPT

#---| PACKAGES |-------------------------------------------------------------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("reshape2"))
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| PARSE ARGS |-----------------------------------------------------------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--baseDir"), action="store", default=NULL, type="character", help="BaseDir"),   
        make_option(c("--batchID"), action="store", default=NULL, type="character", help="BatchID"),
        make_option(c("--seqID"), action="store", default=NULL, type="character", help="SeqID"),
        make_option(c("--databaseImport"), action="store", default=NULL, type="character", help="DatabaseImport")
    )
    #---------------------------------------------------------------------------  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #---------------------------------------------------------------------------
    BaseDir        <- ARGS$baseDir
    BatchID        <- ARGS$batchID
    SeqID          <- ARGS$seqID
    DatabaseImport <- ARGS$databaseImport
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| READ MULTI-QC RESULTS |------------------------------------------------------------------------------------------------------------------------#
    message("|--->> Read QC results...")
    multiqcResDir <- sprintf("%s/%s/%s/qcfiles/qc_res/%s_qc_data", BaseDir, BatchID, SeqID, SeqID)
    #---------------------------------------------------------------------------
    stats_general <- read.delim(sprintf("%s/multiqc_general_stats.txt",        multiqcResDir))
    stats_rnaseqc <- read.delim(sprintf("%s/multiqc_rna_seqc.txt",             multiqcResDir))
    stats_picard  <- read.delim(sprintf("%s/multiqc_picard_RnaSeqMetrics.txt", multiqcResDir))
    stats_star    <- read.delim(sprintf("%s/multiqc_star.txt",                 multiqcResDir))
    #---------------------------------------------------------------------------
    picard_hist_file <- sprintf("%s/%s/%s/qcfiles/%s.gatk.rnaseq.metrics", BaseDir, BatchID, SeqID, SeqID)
    skipRows         <- system(sprintf("grep -n '^## HISTOGRAM' %s", picard_hist_file), intern=T)
    skipRows         <- as.numeric(unlist(strsplit(skipRows, ":"))[1])
    histo_data       <- read.delim(sprintf("%s", picard_hist_file), header=TRUE, skip=skipRows)
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| QC RESULT TABLE |------------------------------------------------------------------------------------------------------------------------------#
    message("|--->> QC stats integration and summarization...")
    QC_RES <- data.frame(
        seq_folder = BatchID,
        seq_id     = SeqID,
        #---| STATS : GENERAL --------------------------------------------------
        pf_rds            = stats_general %>% filter( Sample == SeqID ) %>% .$fastp_mqc.generalstats.fastp.filtering_result_passed_filter_reads,
        pf_q30_base       = stats_general %>% filter( Sample == SeqID ) %>% .$fastp_mqc.generalstats.fastp.after_filtering_q30_bases,
        uniq_mapped_rds   = stats_general %>% filter( Sample == SeqID ) %>% .$STAR_mqc.generalstats.star.uniquely_mapped,
        pct_pf_rds        = stats_general %>% filter( Sample == SeqID ) %>% .$fastp_mqc.generalstats.fastp.pct_surviving %>% round(1),
        pct_pf_q30_rate   = stats_general %>% filter( Sample == SeqID ) %>% .$fastp_mqc.generalstats.fastp.after_filtering_q30_rate,
        pct_uniq_mapped   = stats_general %>% filter( Sample == SeqID ) %>% .$STAR_mqc.generalstats.star.uniquely_mapped_percent,
        pct_duplicates    = stats_general %>% filter( Sample == SeqID ) %>% .$fastp_mqc.generalstats.fastp.pct_duplication %>% round(1),
        pct_gc_content    = stats_general %>% filter( Sample == SeqID ) %>% .$fastp_mqc.generalstats.fastp.after_filtering_gc_content %>% round(3),
        rate_rRNA         = stats_general %>% filter( Sample == SeqID ) %>% .$RNA.SeQC_mqc.generalstats.rna_seqc.rRNA_rate,
        avg_seq_length    = stats_general %>% filter( Sample == SeqID ) %>% .$FastQC_mqc.generalstats.fastqc.avg_sequence_length,
        median_seq_length = stats_general %>% filter( Sample == SeqID ) %>% .$FastQC_mqc.generalstats.fastqc.median_sequence_length,
        #---| STATS : PICARD ---------------------------------------------------
        pf_base              = stats_picard$PF_BASES,
        pf_align_base        = stats_picard$PF_ALIGNED_BASES,
        pf_unaligned_base    = stats_picard$PF_NOT_ALIGNED_BASES,
        coding_base          = stats_picard$CODING_BASES,
        ribosomal_base       = stats_picard$RIBOSOMAL_BASES,
        intronic_base        = stats_picard$INTRONIC_BASES,
        utr_base             = stats_picard$UTR_BASES,
        intergenic_base      = stats_picard$INTERGENIC_BASES,
        pct_usable_base      = stats_picard$PCT_USABLE_BASES,
        pct_mrna_base        = stats_picard$PCT_MRNA_BASES,
        pct_coding_base      = stats_picard$PCT_CODING_BASES,
        pct_ribosomal_base   = stats_picard$PCT_RIBOSOMAL_BASES,
        pct_intronic_base    = stats_picard$PCT_INTRONIC_BASES,
        pct_utr_base         = stats_picard$PCT_UTR_BASES,
        pct_intergenic_base  = stats_picard$PCT_INTERGENIC_BASES,
        pct_r1_tx_strand_rds = stats_picard$PCT_R1_TRANSCRIPT_STRAND_READS,
        pct_r2_tx_strand_rds = stats_picard$PCT_R2_TRANSCRIPT_STRAND_READS,
        r1_tx_strand_rds     = stats_picard$NUM_R1_TRANSCRIPT_STRAND_READS,
        r2_tx_strand_rds     = stats_picard$NUM_R2_TRANSCRIPT_STRAND_READS,
        correct_strand_rds   = stats_picard$CORRECT_STRAND_READS,
        incorrect_strand_rds = stats_picard$INCORRECT_STRAND_READS,
        ignored_rds          = stats_picard$IGNORED_READS,
        unexplained_rds      = stats_picard$NUM_UNEXPLAINED_READS,
        median_cv_cov        = stats_picard$MEDIAN_CV_COVERAGE,
        median_3bias         = stats_picard$MEDIAN_3PRIME_BIAS,
        median_5bias         = stats_picard$MEDIAN_5PRIME_BIAS,
        median_5to3bias      = stats_picard$MEDIAN_5PRIME_TO_3PRIME_BIAS,
        #---| STATS : RNASEQC --------------------------------------------------
        seqc_total_rds           = stats_rnaseqc$Total.Read.Number,
        seqc_mapped_rds          = stats_rnaseqc$Mapped.Reads,
        seqc_pf_uniq_mapped_rds  = stats_rnaseqc$Unique.Mapping..Vendor.QC.Passed.Reads,
        seqc_uniq_mapped_rds     = stats_rnaseqc$Mapped.Unique.Reads,
        seqc_hq_rds              = stats_rnaseqc$High.Quality.Reads,
        seqc_lowqual_rds         = stats_rnaseqc$Low.Quality.Reads,
        seqc_low_mapq_rds        = stats_rnaseqc$Low.Mapping.Quality,
        seqc_rrna_rds            = stats_rnaseqc$rRNA.Reads,
        seqc_intragenic_rds      = stats_rnaseqc$Intragenic.Reads,
        seqc_intronic_rds        = stats_rnaseqc$Intronic.Reads,
        seqc_non_globin_rds      = stats_rnaseqc$Non.Globin.Reads,
        seqc_paired_rds          = stats_rnaseqc$Total.Mapped.Pairs,
        seqc_unpaired_rds        = stats_rnaseqc$Unpaired.Reads,
        seqc_ambiguous_rds       = stats_rnaseqc$Ambiguous.Reads,
        seqc_alternative_align   = stats_rnaseqc$Alternative.Alignments,
        seqc_rds_intron_exon     = stats_rnaseqc$Reads.used.for.Intron.Exon.counts,
        seqc_total_base          = stats_rnaseqc$Total.Bases,
        seqc_pct_hq_ambigous     = stats_rnaseqc$High.Quality.Ambiguous.Alignment.Rate,
        seqc_read_length         = stats_rnaseqc$Read.Length,
        seqc_avg_split_rds       = stats_rnaseqc$Avg..Splits.per.Read,
        seqc_exonic_rds          = stats_rnaseqc$Exonic.Reads,
        seqc_intergenic_rds      = stats_rnaseqc$Intergenic.Reads,
        seqc_genes_detect        = stats_rnaseqc$Genes.Detected,
        seqc_mismatch_base       = stats_rnaseqc$Base.Mismatch,
        seqc_pct_intronic        = stats_rnaseqc$Intronic.Rate,
        seqc_pct_hq_intronic     = stats_rnaseqc$High.Quality.Intronic.Rate,
        seqc_pct_intragenic      = stats_rnaseqc$Intragenic.Rate,
        seqc_pct_intergenic      = stats_rnaseqc$Intergenic.Rate,
        seqc_pct_exonic          = stats_rnaseqc$Exonic.Rate,
        seqc_pct_hq_exonic       = stats_rnaseqc$High.Quality.Exonic.Rate,
        seqc_pct_hq_intragenic   = stats_rnaseqc$High.Quality.Intragenic.Rate,
        seqc_pct_hq_intergenic   = stats_rnaseqc$High.Quality.Intergenic.Rate,
        seqc_pct_ambiguous_align = stats_rnaseqc$Ambiguous.Alignment.Rate,
        seqc_median_tx_cov       = stats_rnaseqc$Median.of.Avg.Transcript.Coverage,
        seqc_median_exon_cv      = stats_rnaseqc$Median.Exon.CV,
        seqc_median_tx_cov_cv    = stats_rnaseqc$Median.of.Transcript.Coverage.CV,
        #---| STATS : STAR -----------------------------------------------------
        star_total_rds            = stats_star$total_reads,
        star_uniq_mapped_rds      = stats_star$uniquely_mapped,
        star_multi_mapped_rds     = stats_star$multimapped,
        star_unmapped             = stats_star$unmapped_tooshort + stats_star$unmapped_other,
        star_total_splices        = stats_star$num_splices,
        star_annot_splices        = stats_star$num_annotated_splices,
        star_GTAG_splices         = stats_star$num_GTAG_splices,
        star_GCAG_splices         = stats_star$num_GCAG_splices,
        star_ATAC_splices         = stats_star$num_ATAC_splices,
        star_noncanonical_splices = stats_star$num_noncanonical_splices,
        star_pct_uniq_mapped      = stats_star$uniquely_mapped_percent,
        star_pct_multi_mapped     = stats_star$multimapped_percent,
        star_pct_mismatch         = stats_star$mismatch_rate,
        #---| Picard Coverage Bias Histogram Values |---------------------------
        cov_bias_hist_value       = paste(histo_data$All_Reads.normalized_coverage, collapse=";")
    )
#----------------------------------------------------------------------------------------------------------------------------------------------------#        

#---| REPORT VALUES TABLE |--------------------------------------------------------------------------------------------------------------------------#
    message("|--->> Prepare table for Report...")
    QC_REPORT <- QC_RES[, c(
        "seq_folder","seq_id",
        "seqc_total_rds","star_total_rds","star_uniq_mapped_rds","star_multi_mapped_rds","star_pct_uniq_mapped","star_pct_multi_mapped",
        "pct_pf_q30_rate","pct_duplicates","pct_gc_content","pct_coding_base","pct_utr_base","pct_intronic_base","pct_intergenic_base","pct_ribosomal_base",
        "pct_r1_tx_strand_rds","pct_r2_tx_strand_rds","median_5to3bias","cov_bias_hist_value"
    )] %>% mutate(
        star_pct_uniq_mapped  = round(star_pct_uniq_mapped, 1),
        star_pct_multi_mapped = round(star_pct_multi_mapped, 1),
        pct_pf_q30_rate       = round(pct_pf_q30_rate*100, 1),
        pct_duplicates        = round(pct_duplicates, 1),
        pct_gc_content        = round(pct_gc_content*100, 1),
        pct_coding_base       = round(pct_coding_base, 1),
        pct_utr_base          = round(pct_utr_base, 1),
        pct_intronic_base     = round(pct_intronic_base, 1),
        pct_intergenic_base   = round(pct_intergenic_base, 1),
        pct_r1_tx_strand_rds  = round(pct_r1_tx_strand_rds, 2),
        pct_r2_tx_strand_rds  = round(pct_r2_tx_strand_rds, 2),
        median_5to3bias       = round(median_5to3bias, 2)
    )
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| SAVE RESULTS AS FILE |-------------------------------------------------------------------------------------------------------------------------#
    message("|--->> Write stats as TSV...")
    write.table(
        QC_RES,
        sprintf("%s/%s/%s/qcfiles/%s_QC_RESULT_TABLE.tsv", BaseDir, BatchID, SeqID, SeqID),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
#----------------------------------------------------------------------------------------------------------------------------------------------------#        

#---| IMPORT INTO DATABASE |-------------------------------------------------------------------------------------------------------------------------#
    if( DatabaseImport == "true" )
    {
        message("|--->> Import into database...")
        source("/data/wts/params/ruo_wts_db.R")
        options(echo=FALSE)
        #-----------------------------------------------------------------------
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
        clearPreExistData1 <- dbGetQuery(dbCon, sprintf("DELETE FROM qc_bam WHERE seq_folder = '%s' AND seq_id = '%s'", BatchID, SeqID))
        writeNewData1      <- dbWriteTable(dbCon, name='qc_bam', value=QC_RES, row.names=FALSE, append=TRUE)
        clearPreExistData2 <- dbGetQuery(dbCon, sprintf("DELETE FROM qc_report WHERE seq_folder = '%s' AND seq_id = '%s'", BatchID, SeqID))
        writeNewData2      <- dbWriteTable(dbCon, name='qc_report', value=QC_REPORT, row.names=FALSE, append=TRUE)
        dbDisconnect(dbCon)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
    message("|--->> DONE.")
#----------------------------------------------------------------------------------------------------------------------------------------------------#
