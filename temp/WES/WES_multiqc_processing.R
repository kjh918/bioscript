
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("reshape2"))

#---| 
    option_list = list(
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="$SEQ_FOLDER"),
        make_option(c("--seqID")    , action="store", default=NA, type="character", help="$SEQ_ID")
    )
    ARGS = parse_args(OptionParser(option_list=option_list))
    seqFolder = ARGS$seqFolder
    seqID     = ARGS$seqID
#---/
    qcDataDir = sprintf("/data/wes/%s/%s/qcfiles/qc_res/%s_qc_data", seqFolder, seqID, seqID)

# pf ==> pass filter

#---| multiqc
    mqc.general  = read.delim(sprintf("%s/multiqc_general_stats.txt", qcDataDir))
    mqc.dup      = read.delim(sprintf("%s/multiqc_picard_dups.txt", qcDataDir))
    mqc.align    = read.delim(sprintf("%s/multiqc_picard_AlignmentSummaryMetrics.txt", qcDataDir))
    mqc.hs       = read.delim(sprintf("%s/multiqc_picard_HsMetrics.txt", qcDataDir))

    qc_raw = data.frame(
        seq_folder = seqFolder,
        seq_id     = seqID,
        date       = today(), 
        # from general stats
        pct_q30_pf        = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.after_filtering_q30_rate *100, 3),
        q30_base_pf       = mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.after_filtering_q30_bases ,
        q30_rds_pf        = mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.filtering_result_passed_filter_reads,
        pct_fq_pf         = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.pct_surviving, 3),
        pct_fq_adapter    = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.pct_adapter, 3),
        pct_gc_content    = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.after_filtering_gc_content *100, 3),
        fq_r1_seqs        = mqc.general %>% filter(Sample == paste0(seqID, "_R1")) %>% .$FastQC_mqc.generalstats.fastqc.total_sequences,
        fq_r2_seqs        = mqc.general %>% filter(Sample == paste0(seqID, "_R2")) %>% .$FastQC_mqc.generalstats.fastqc.total_sequences,
        # from picard duplicates
        rds_pairs_dup_chk = mqc.dup$READ_PAIRS_EXAMINED,
        rds_dup           = mqc.dup$READ_PAIR_DUPLICATES,
        rds_opt_dup       = mqc.dup$READ_PAIR_OPTICAL_DUPLICATES,
        rds_2nd_or_suppl  =  mqc.dup$SECONDARY_OR_SUPPLEMENTARY_RDS,
        pct_rds_dup       = round(mqc.dup$PERCENT_DUPLICATION*100, 3)
    )

    qc_bam_raw = data.frame(
        seq_folder = seqFolder,
        seq_id     = seqID,
        date       = today(), 
        bam_status = c('m','m.dm','m.dm.br', 'm.dm.br.lr'),
        # picard alignment summary
        rds_pf                = mqc.align$TOTAL_READS,                               # mqc.hs$PF_READS
        rds_pf_align          = mqc.align$PF_READS_ALIGNED,
        rds_pf_align_hq       = mqc.align$PF_HQ_ALIGNED_READS,
        rds_align_pairs       = mqc.align$READS_ALIGNED_IN_PAIRS,
        rds_pf_improper_pairs = mqc.align$PF_READS_IMPROPER_PAIRS,
        rds_pf_noise          = mqc.align$PF_NOISE_READS,
        pct_rds_pf            = round(mqc.align$PCT_PF_READS*100,3),
        pct_rds_pf_align      = round(mqc.align$PCT_PF_READS_ALIGNED*100,3),
        pct_rds_align_pairs   = round(mqc.align$PCT_READS_ALIGNED_IN_PAIRS*100,3),
        pct_mismatch          = round(mqc.align$PF_MISMATCH_RATE*100,3),
        pct_adapter           = round(mqc.align$PCT_ADAPTER*100,3),
        pct_softclip          = round(mqc.align$PCT_SOFTCLIP*100,3),
        pct_hardclip          = round(mqc.align$PCT_HARDCLIP*100,3),
        pct_pf_hq_error       = round(mqc.align$PF_HQ_ERROR_RATE*100,3),
        pct_pf_indel          = round(mqc.align$PF_INDEL_RATE*100,3),
        pct_chimeras          = round(mqc.align$PCT_CHIMERAS*100,3),
        strand_balance        = round(mqc.align$STRAND_BALANCE, 4),
        avg_softclip_length   = round(mqc.align$AVG_POS_3PRIME_SOFTCLIP_LENGTH,3),
        pf_hq_mismatch_median = mqc.align$PF_HQ_MEDIAN_MISMATCHES,
        rds_length_mean       = mqc.align$MEAN_READ_LENGTH,
        rds_length_align_mean = round(mqc.align$MEAN_ALIGNED_READ_LENGTH, 3),
        bad_cycles            = mqc.align$BAD_CYCLES,
        base_pf_align         = mqc.align$PF_ALIGNED_BASES,                          # mqc.hs$PF_BASES_ALIGNED
        base_pf_align_hq      = mqc.align$PF_HQ_ALIGNED_BASES,
        base_pf_align_hq_q20  = mqc.align$PF_HQ_ALIGNED_Q20_BASES,
        pct_base_align_q20    = round(mqc.align$PF_HQ_ALIGNED_Q20_BASES/mqc.align$PF_ALIGNED_BASES*100, 3),
        # from picard hsmetric
        rds_total             = mqc.hs$TOTAL_READS,
        rds_pf_uniq           = mqc.hs$PF_UNIQUE_READS,
        rds_pf_uniq_align     = mqc.hs$PF_UQ_READS_ALIGNED,
        base_pf               = mqc.hs$PF_BASES,
        base_pf_uniq_align    = mqc.hs$PF_UQ_BASES_ALIGNED,
        base_on_target        = mqc.hs$ON_TARGET_BASES,
        base_on_bait          = mqc.hs$ON_BAIT_BASES,
        base_off_bait         = mqc.hs$OFF_BAIT_BASES,
        base_near_bait        = mqc.hs$NEAR_BAIT_BASES,
        base_genome_ref       = mqc.hs$GENOME_SIZE,
        pct_base_on_bait      = round(mqc.hs$PCT_USABLE_BASES_ON_BAIT*100, 3),
        pct_base_on_target    = round(mqc.hs$PCT_USABLE_BASES_ON_TARGET*100, 3),
        pct_off_bait          = round(mqc.hs$PCT_OFF_BAIT*100, 3),
        pct_on_bait_vs_select = round(mqc.hs$ON_BAIT_VS_SELECTED*100, 3),
        pct_zero_cov_target   = round(mqc.hs$ZERO_CVG_TARGETS_PCT*100, 3),
        pct_exc_mapq          = round(mqc.hs$PCT_EXC_MAPQ*100, 3),
        pct_exc_baseq         = round(mqc.hs$PCT_EXC_BASEQ*100, 3),
        pct_exc_adapter       = round(mqc.hs$PCT_EXC_ADAPTER*100, 3),
        pct_exc_off_target    = round(mqc.hs$PCT_EXC_OFF_TARGET*100, 3),
        pct_exc_dup           = round(mqc.hs$PCT_EXC_DUPE*100, 3),
        pct_exc_overlap       = round(mqc.hs$PCT_EXC_OVERLAP*100, 3),
        pct_base_target_1x    = round(mqc.hs$PCT_TARGET_BASES_1X*100, 3),
        pct_base_target_10x   = round(mqc.hs$PCT_TARGET_BASES_10X*100, 3),
        pct_base_target_30x   = round(mqc.hs$PCT_TARGET_BASES_30X*100, 3),
        pct_base_target_50x   = round(mqc.hs$PCT_TARGET_BASES_50X*100, 3),
        pct_base_target_100x  = round(mqc.hs$PCT_TARGET_BASES_100X*100, 3),
        pct_base_target_250x  = round(mqc.hs$PCT_TARGET_BASES_250X*100, 3),
        pct_base_target_500x  = round(mqc.hs$PCT_TARGET_BASES_500X*100, 3),
        pct_base_target_1000x = round(mqc.hs$PCT_TARGET_BASES_1000X*100, 3),
        territory_base        = mqc.hs$BAIT_TERRITORY,
        territory_target      = mqc.hs$TARGET_TERRITORY,
        depth_bait_mean       = round(mqc.hs$MEAN_BAIT_COVERAGE, 3),
        depth_target_mean     = round(mqc.hs$MEAN_TARGET_COVERAGE, 3),
        depth_target_median   = mqc.hs$MEDIAN_TARGET_COVERAGE,
        bait_design_eff_rate  = round(mqc.hs$BAIT_DESIGN_EFFICIENCY*100, 3),
        fold_enrich           = round(mqc.hs$FOLD_ENRICHMENT, 3),
        het_snp_sensitivity   = round(mqc.hs$HET_SNP_SENSITIVITY, 3),
        het_snp_q             = mqc.hs$HET_SNP_Q
    )
#---/
#---| mosdepth
    mosd  = read.delim(sprintf("%s/mosdepth_cumcov_regions_dist.txt", qcDataDir))   
    mosd2 = melt(mosd) %>% select(2:3) %>% rename( depth=1, coverage=2) %>% mutate(depth = as.numeric(as.character(gsub("^X","", depth)))) %>% arrange((depth))
    mosd.lbl = mosd2 %>% filter(depth %in% c(30,100)) 
    mosd.lbl$lbl = apply(mosd.lbl, 1, function(y) paste0(y[1], "X=", y[2],"%"))

    mosch = read.delim(sprintf("%s/mosdepth_perchrom_regions.txt", qcDataDir), check.names=F) %>% melt() %>%
        select(2:3) %>% rename(chromosome=1, coverage=2) %>% mutate( chromosome = factor(chromosome, levels=(c(1:22,'X','Y'))))
    #
        #ggplot( mosd2, aes(x=depth, y=coverage)) + geom_line( size = 1.2, color = '#2261e6' ) +
        #    labs( x = 'Depth(X)', y = 'Coverage(%)' ) + 
        #    geom_vline(xintercept=c(30,100), linetype=2, color='#e64522') +
        #    coord_cartesian( xlim = c(0,200), ylim = c(0, 100) ) +
        #    theme(
        #        panel.background = element_rect(fill='white'), panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'),
        #        axis.text.x  = element_text(colour="black", size=10, vjust=1, hjust=1),
        #        axis.text.y  = element_text(colour="black", size=10, hjust=0.25),
        #        axis.title.x = element_text(colour="black",size=10*1.2 ,vjust=1),
        #        axis.title.y = element_text(colour="black",size=10*1.2 ,vjust=1),
        #        axis.line.x  = element_line(colour = "#6f6d6d"), 
        #        axis.line.y  = element_line(colour = "#6f6d6d"),
        #    ) +
        #   geom_point( x = 30, y = mosd2[mosd2$depth==30, 'coverage'], size = 3, color='#e64522') +
        #    geom_point( x = 100, y = mosd2[mosd2$depth==100, 'coverage'], size = 3, color='#e64522') +  
        #    geom_text( data=mosd.lbl, aes(x = depth, y =coverage, label=lbl), hjust=0, nudge_x = 1, nudge_y =3  ) +
        #    ggtitle("Target Coverage (Base)")

        # ggplot(mosch, aes(x=chromosome, y=coverage)) + geom_bar(position="dodge", stat="identity", width =0.7)
    #
#---| alfred
    alf.chr  = read.delim(sprintf("/data/wes/%s/%s/qcfiles/%s.realign.alfred.chr.map.stats.txt", seqFolder, seqID, seqID), check.names=F)
    alf.chr2 = alf.chr %>% filter( Chrom %in% c(1:22, 'X','Y','MT') ) 
    alf.tcov = read.delim(sprintf("/data/wes/%s/%s/qcfiles/%s.realign.alfred.target.coverage.txt", seqFolder, seqID, seqID), check.names=F)

    qc_extra = data.frame(
        seq_folder = seqFolder,
        seq_id     = seqID,
        date       = today(), 
        # mosdepth
        mosdepth_cov_x_depth = paste(mosd2$depth, collapse=";"),
        mosdepth_cov_y_cov   = paste(mosd2$coverage, collapse=";"),
        mosdepth_chr_cov_x_chr = paste(mosch$chromosome, collapse=";"),
        mosdepth_chr_cov_y_cov = paste(round(mosch$coverage, 3), collapse=";"),
        # alfred
        alf_chr_chr = paste(alf.chr2$Chrom, collapse=";"),
        alf_chr_map_rate = paste(round(alf.chr2$MappedFraction,3), collapse=";"),
        alf_chr_obsexp_ratio = paste(round(alf.chr2$ObsExpRatio,3), collapse=";"),
        alf_tg_cov_dp = "1x;30x;50x;100x;150x;200x;250x;300x;500x;1000x",
        alf_tg_cov = paste(
            round(alf.tcov %>% filter(AvgCov >= 1) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 30) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 50) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 100) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 150) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 200) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 250) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 300) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 400) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 500) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
            round(alf.tcov %>% filter(AvgCov >= 1000) %>% summarise(sum(AvgCov)/nrow(alf.tcov)),3),
        sep=";")
    )
#---/

#---| qc_raw re-format
    qc_raw2 = cbind(
        qc_raw,
        qc_bam_raw %>% filter( bam_status == 'm.dm.br.lr' ) %>% select( 5:ncol(qc_bam_raw) ),
        qc_extra %>% select( 4:ncol(qc_extra) )
    )

#---| Import to DB
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='gcx_wes' )
    # qc-report
    clearPreExistData.report = dbGetQuery(dbCon, sprintf("DELETE FROM qc_report WHERE seq_folder = '%s' AND seq_id = '%s'", seqFolder, seqID) )
    writeNewData.report      = dbWriteTable(dbCon, name = 'qc_report', value = qc_raw2, row.names = FALSE, append = TRUE)
    # qc-bam
    clearPreExistData.bam = dbGetQuery(dbCon, sprintf("DELETE FROM qc_bam WHERE seq_folder = '%s' AND seq_id = '%s'", seqFolder, seqID) )
    writeNewData.bam      = dbWriteTable(dbCon, name = 'qc_bam', value = qc_bam_raw, row.names = FALSE, append = TRUE)
    dbDisconnect(dbCon)

#----| message
    message(sprintf("SeqID : %s , DB imported.", seqID))