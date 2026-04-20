
#   script version v1.1
#   update 2023-12-08
#   summary of change : merge default-vep, refseq-vep (remove gatk-funcotator)

#---| PACKAGES |-----------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("reshape2"))
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="$SEQ_FOLDER"),
        make_option(c("--seqID")    , action="store", default=NA, type="character", help="$SEQ_ID"),
        make_option(c("--baseDir")  , action="store", default=NA, type="character", help="$BASE_DIR"),
        make_option(c("--statusQC") , action="store", default=NA, type="character", help="'FASTQ' or 'BAM'"),
        make_option(c("--summaryFastqScreen"), action="store", default=TRUE, type="logical", help="Run QC Summary for FastqScreen"),
        make_option(c("--importDB") , action="store", default=FALSE, type="logical", 
            help="Import Into DB. default = FALSE. only activated when statusQC = 'BAM' ")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    seqFolder        <- ARGS$seqFolder
    seqID            <- ARGS$seqID
    baseDir          <- ARGS$baseDir
    qcLevel          <- ARGS$statusQC
    includeFQScreen  <- ARGS$summaryFastqScreen
    importDB         <- ARGS$importDB
#--------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTIONS |-----------------------------------------------------------------------#
    source("/data/wes/params/ruo_wes_db.R")
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    qcDataDir <- sprintf("%s/%s/%s/qcfiles/qc_res/%s_qc_data", baseDir, seqFolder, seqID, seqID)
    LOG_DIR   <- sprintf("%s/%s/%s/log", baseDir, seqFolder, seqID)
    if( !dir.exists(LOG_DIR) ){ system(sprintf("mkdir -p %s", LOG_DIR)) }
#--------------------------------------------------------------------------------------------------#

#---| LOG FILE |-----------------------------------------------------------------------------------#
    LOG_FILE <- sprintf("%s/%s.multiqc.summary.%s.log", LOG_DIR, seqID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES MultiQC QC-Result Summary log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE             : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER           : %s", seqFolder), LOG_FILE, append=T) 
    write(sprintf(" SEQ ID               : %s", seqID), LOG_FILE, append=T) 
    write(sprintf(" QC LEVEL             : %s", qcLevel), LOG_FILE, append=T) 
    write(sprintf(" INCLUDE FASTQ-SCREEN : %s", includeFQScreen), LOG_FILE, append=T) 
    write(sprintf(" DB IMPORT            : %s", importDB), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

#==================================================================================================# 
    write(sprintf("%s | MULTIQC RESULT INTEGRATION & SUMMARY START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("MULTIQC RESULT INTEGRATION & SUMMARY START.")
#==================================================================================================# 

#---| READ GENERAL STATS |-------------------------------------------------------------------------#
    MQC_GENEREAL_RES <- sprintf("%s/multiqc_general_stats.txt", qcDataDir)
    if( file.exists(MQC_GENEREAL_RES) ) { mqc.general <- read.delim(MQC_GENEREAL_RES) 
    }else{ message("\n ** NO General stats result file. Stopped.") ; quit( save="no", status=0 ) }
#--------------------------------------------------------------------------------------------------#

#---| QC TABLE GENERATION |------------------------------------------------------------------------#
    if( qcLevel == "FASTQ" )
    {
        #---| FASTQ LEVEL QC : RESULT |------------------------------------------------------------#
        fq.level.qc_res <- data.frame(
            seq_folder     = seqFolder,
            seq_id         = seqID,
            date           = today(),
            pct_q30_pf     = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.after_filtering_q30_rate *100, 3),
            q30_base_pf    = mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.after_filtering_q30_bases ,
            q30_rds_pf     = mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.filtering_result_passed_filter_reads,
            pct_fq_pf      = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.pct_surviving, 3),
            pct_fq_adapter = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.pct_adapter, 3),
            pct_gc_content = round( mqc.general %>% filter(Sample==seqID) %>% .$fastp_mqc.generalstats.fastp.after_filtering_gc_content *100, 3),
            fq_r1_seqs     = mqc.general %>% filter(Sample == paste0(seqID, "_R1")) %>% .$FastQC_mqc.generalstats.fastqc.total_sequences,
            fq_r2_seqs     = mqc.general %>% filter(Sample == paste0(seqID, "_R2")) %>% .$FastQC_mqc.generalstats.fastqc.total_sequences
        )
        #------------------------------------------------------------------------------------------#
        #---| FASTQ LEVEL QC : SAVE AS FILE |------------------------------------------------------#
            fq.level.res.filename <- sprintf("%s/%s/%s/qcfiles/%s.FASTQ.ONLY.QC.RESULT.txt", baseDir, seqFolder, seqID, seqID)
            write.table( fq.level.qc_res, file=fq.level.res.filename, quote=F, col.names=T, row.names=F, sep="\t" )
        #------------------------------------------------------------------------------------------#

        #==========================================================================================#
            write(sprintf("%s | %s LEVEL QC RESULT CREATED AND SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S"), qcLevel), LOG_FILE, append=T)
            message(sprintf("%s | %s LEVEL QC RESULT CREATED AND SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S"), qcLevel))
        #==========================================================================================#
    }else{
        #---| LOAD BAM QC DATA  |------------------------------------------------------------------#
        MQC_DUP_RES <- sprintf("%s/multiqc_picard_dups.txt", qcDataDir)
        if( file.exists(MQC_DUP_RES) ){ mqc.dup <- read.delim(MQC_DUP_RES) 
        }else{ message("\n ** NO MarkDuplicates stats result file. Stopped.") ; quit( save="no", status=0 )}
        #----------------------------------------------------------------------#
        MQC_ALIGN_RES <- sprintf("%s/multiqc_picard_AlignmentSummaryMetrics.txt", qcDataDir)
        if( file.exists(MQC_ALIGN_RES) ){ mqc.align <- read.delim(MQC_ALIGN_RES) 
        }else{ message("\n ** NO BAM Alignment stats result file. Stopped.") ; quit( save="no", status=0 )}
        #----------------------------------------------------------------------#
        MQC_HSMTERIC_RES <- sprintf("%s/multiqc_picard_HsMetrics.txt", qcDataDir)
        if( file.exists(MQC_HSMTERIC_RES) ){ mqc.hs <- read.delim(MQC_HSMTERIC_RES) 
        }else{ message("\n ** NO HS-Metrics stats result file. Stopped.") ; quit( save="no", status=0 )}
        #----------------------------------------------------------------------#
        MOSDEPTH_COV_RES <- sprintf("%s/mosdepth_cumcov_regions_dist.txt", qcDataDir)
        if( file.exists(MOSDEPTH_COV_RES) ){  mosd <- read.delim(MOSDEPTH_COV_RES)   
        }else{ message("\n ** NO Mosdepth Coverage result file. Stopped.") ; quit( save="no", status=0 )}
        #----------------------------------------------------------------------#
        mosd2 <- melt(mosd) %>% 
            select(2:3) %>% 
            dplyr::rename( depth=1, coverage=2) %>% 
            mutate(depth = as.numeric(as.character(gsub("^X","", depth)))) %>% 
            arrange((depth))
        #----------------------------------------------------------------------#
        mosd.lbl     <- mosd2 %>% filter(depth %in% c(30,100)) 
        mosd.lbl$lbl <- apply(mosd.lbl, 1, function(y) paste0(y[1], "X=", y[2],"%"))
        #----------------------------------------------------------------------#
        MOSDEPTH_CHR_RES <- sprintf("%s/mosdepth_perchrom_regions.txt", qcDataDir)
        if( file.exists(MOSDEPTH_CHR_RES) ){  mosch <- read.delim(MOSDEPTH_CHR_RES)   
        }else{ message("\n ** NO Mosdepth Chromosome stats result file. Stopped.") ; quit( save="no", status=0 )}
        #----------------------------------------------------------------------#
        mosch <- mosch %>%
            select(2:3) %>% 
            dplyr::rename(chromosome=1, coverage=2) %>% 
            mutate( chromosome = factor(chromosome, levels=(c(1:22,'X','Y'))))
        #----------------------------------------------------------------------#
        ALFRED_COV_RES <- sprintf("%s/%s/%s/qcfiles/%s.alfred.target.coverage.txt", baseDir, seqFolder, seqID, seqID) 
        if( file.exists(ALFRED_COV_RES) ){  alf.tcov <- read.delim(ALFRED_COV_RES, check.names=F)   
        }else{ message("\n ** NO Alfred Coverage stats result file. Stopped.") ; quit( save="no", status=0 )} 
        #----------------------------------------------------------------------#
        ALFRED_CHR_RES <- sprintf("%s/%s/%s/qcfiles/%s.alfred.chr.map.stats.txt", baseDir, seqFolder, seqID, seqID)
        if( file.exists(ALFRED_CHR_RES) ){  alf.chr <- read.delim(ALFRED_CHR_RES, check.names=F)   
        }else{ message("\n ** NO Alfred Chromosome stats result file. Stopped.") ; quit( save="no", status=0 )}
        #----------------------------------------------------------------------#
        alf.chr2 <- alf.chr %>% 
            filter( Chrom %in% c(1:22, 'X','Y','MT') ) 
        #----------------------------------------------------------------------#
    
        #---| QC-TABLES  |-------------------------------------------------------------------------#
        qc.FASTQ <- data.frame(
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
        #----------------------------------------------------------------------#
        qc.BAM <- data.frame(
            seq_folder = seqFolder,
            seq_id     = seqID,
            date       = today(), 
            bam_status = c('m','m.dm','m.dm.lr', 'm.dm.lr.br'),
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
        #----------------------------------------------------------------------#
        qc.COV <- data.frame(
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
        #------------------------------------------------------------------------------------------#

        #---| FINAL QC RESULTS |-------------------------------------------------------------------#
        QC_REPORT_RES <- cbind(
            qc.FASTQ,
            qc.BAM %>% filter( bam_status == 'm.dm.lr.br' ) %>% select( 5:ncol(qc.BAM) ),
            qc.COV %>% select( 4:ncol(qc.COV) )
        )
        QC_BAM_RES <- qc.BAM
        #------------------------------------------------------------------------------------------#

        #---| WRITE AS FILES |---------------------------------------------------------------------#
        BAM_QC_RES_NAME    <- sprintf("%s/%s/%s/qcfiles/%s.BAM.QC.RESULT.txt",     baseDir, seqFolder, seqID, seqID)
        QC_REPORT_RES_NAME <- sprintf("%s/%s/%s/qcfiles/%s.QC.REPORT.RESULT.txt", baseDir, seqFolder, seqID, seqID)
        #--------------------------------------------------------------------------#
        write.table( QC_BAM_RES,    file=BAM_QC_RES_NAME,    quote=F, col.names=T, row.names=F, sep="\t" )
        write.table( QC_REPORT_RES, file=QC_REPORT_RES_NAME, quote=F, col.names=T, row.names=F, sep="\t" )
        #------------------------------------------------------------------------------------------#

        #==========================================================================================#
            write(sprintf("%s | %s LEVEL QC RESULT CREATED AND SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S"), qcLevel), LOG_FILE, append=T)
            message(sprintf("%s | %s LEVEL QC RESULT CREATED AND SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S"), qcLevel))
        #==========================================================================================#
    }

#---| FASTQ-SCREEN |-------------------------------------------------------------------------------#
    if( includeFQScreen )
    {
        #---| FASTQ SCREEN QC TABLE |--------------------------------------------------------------#
        FASTQ_SCREEN_RES <- sprintf("%s/%s/%s/qcfiles/%s_R1_screen.txt", baseDir, seqFolder, seqID, seqID)
        if( file.exists(FASTQ_SCREEN_RES) ){ fqscnRes <- read.delim(FASTQ_SCREEN_RES, skip=1, check.names=FALSE )
        }else{ message("\n ** NO FastqScreen stats result file. Stopped.") ; quit( save="no", status=0 )}
        #----------------------------------------------------------------------#
        fixCols    <- c("genome","tested_reads","unmap","pct_unmap","H1_G1","pct_H1_G1","Hm_G1","pct_Hm_G1",
                        "H1_Gm","pct_H1_Gm","Hm_Gm","pct_Hm_Gm")
        pickCols   <- c("pct_H1_G1","pct_Hm_G1","pct_H1_Gm","pct_Hm_Gm","pct_unmap")
        reportCols <- c("seq_folder","seq_id","fastq","species","human.pct_H1_G1","mouse.pct_H1_G1","rat.pct_H1_G1",
                        "ecoli.pct_H1_G1","phix.pct_H1_G1","pseudomonas.pct_H1_G1", "mycoplasma.pct_H1_G1",
                        "adapters.pct_H1_G1","vectors.pct_H1_G1")
        #----------------------------------------------------------------------#
        colnames(fqscnRes) <- fixCols
        qc.FASTQSCREEN <- data.frame(
            human       = fqscnRes %>% filter (genome == "Human" )       %>% select(any_of(pickCols)),
            mouse       = fqscnRes %>% filter (genome == "Mouse" )       %>% select(any_of(pickCols)),
            rat         = fqscnRes %>% filter (genome == "Rat" )         %>% select(any_of(pickCols)),
            ecoli       = fqscnRes %>% filter (genome == "Ecoli" )       %>% select(any_of(pickCols)),
            phix        = fqscnRes %>% filter (genome == "PhiX" )        %>% select(any_of(pickCols)),
            pseudomonas = fqscnRes %>% filter (genome == "Pseudomonas" ) %>% select(any_of(pickCols)),
            mycoplasma  = fqscnRes %>% filter (genome == "Mycoplasma" )  %>% select(any_of(pickCols)),
            adapters    = fqscnRes %>% filter (genome == "Adapters" )    %>% select(any_of(pickCols)),
            vectors     = fqscnRes %>% filter (genome == "Vectors" )     %>% select(any_of(pickCols))
        ) %>% mutate( seq_folder=seqFolder, seq_id=seqID )
        #----------------------------------------------------------------------#
        qc.FASTQSCREEN$species = ""
        #----------------------------------------------------------------------#
        qc.summary.FASTQSCREEN <- qc.FASTQSCREEN %>% select(any_of(reportCols))
        colnames(qc.summary.FASTQSCREEN) <- gsub("\\.pct_H1_G1", "_H1G1", colnames(qc.summary.FASTQSCREEN))
        qc.summary.FASTQSCREEN$species   <- apply(qc.summary.FASTQSCREEN[,4:12], 1, function(z) 
        { paste(gsub("_H1G1","", colnames(qc.summary.FASTQSCREEN)[4:12][which(z >= 10)]),collapse=",") })
        #----------------------------------------------------------------------#
        qc.FASTQSCREEN$species <- qc.summary.FASTQSCREEN$species
        colnames(qc.FASTQSCREEN) <- gsub( "pct_", "" , colnames(qc.FASTQSCREEN) )
        colnames(qc.FASTQSCREEN) <- gsub( "_"   , "" , colnames(qc.FASTQSCREEN) )
        colnames(qc.FASTQSCREEN) <- gsub( "\\." , "_", colnames(qc.FASTQSCREEN) )
        colnames(qc.FASTQSCREEN) <- gsub( "seqfolder" , "seq_folder", colnames(qc.FASTQSCREEN) )
        colnames(qc.FASTQSCREEN) <- gsub( "seqid" , "seq_id", colnames(qc.FASTQSCREEN) )
        #------------------------------------------------------------------------------------------#
    
        #---| FASTQ-SCREEN RESULT SAVE AS FILE |---------------------------------------------------#
        FQSCREEN_SUMMARY_NAME <- sprintf("%s/%s/%s/qcfiles/%s.FASTQ-SCREEN.SUMMARY.txt",   baseDir, seqFolder, seqID, seqID)
        FQSCREEN_QC_RES_NAME  <- sprintf("%s/%s/%s/qcfiles/%s.FASTQ-SCREEN.QC.RESULT.txt", baseDir, seqFolder, seqID, seqID)
        #----------------------------------------------------------------------#
        write.table( qc.summary.FASTQSCREEN, file=FQSCREEN_SUMMARY_NAME, quote=F, col.names=T, row.names=F, sep="\t" )
        write.table( qc.FASTQSCREEN,         file=FQSCREEN_QC_RES_NAME,  quote=F, col.names=T, row.names=F, sep="\t" )
        #------------------------------------------------------------------------------------------#

        #==========================================================================================#
            write(sprintf("%s | FASTQ-SCREEN QC SUMMARY CREATED AND SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message(sprintf("%s | FASTQ-SCREEN QC SUMMARY CREATED AND SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")))
        #==========================================================================================#
    }else{
        #==========================================================================================#
            write(sprintf("%s | FASTQ-SCREEN QC SUMMARY SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message(sprintf("%s | FASTQ-SCREEN QC SUMMARY SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")))
        #==========================================================================================#
    }

#---| IMPORT RESULTS INTO DB |---------------------------------------------------------------------#
    if( importDB )
    {
        dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
        #----------------------------------------------------------------------#
        # FASTQ-SCREEN 
        if( includeFQScreen )
        {
            preClearData <- dbGetQuery(dbCon, sprintf("DELETE FROM qc_fastqscreen WHERE seq_id = '%s'", seqID))
            writeNewData <- dbWriteTable(dbCon, name="qc_fastqscreen", value=qc.FASTQSCREEN, row.names=F, append=T)
        }
        #----------------------------------------------------------------------#
        # QC-BAM
        preClearData <- dbGetQuery(dbCon, sprintf("DELETE FROM qc_bam WHERE seq_folder = '%s' AND seq_id = '%s'", seqFolder, seqID) )    
        writeNewData <- dbWriteTable(dbCon, name = 'qc_bam', value = QC_BAM_RES, row.names = FALSE, append = TRUE)
        #----------------------------------------------------------------------#
        # QC-REPORT
        preClearData <- dbGetQuery(dbCon, sprintf("DELETE FROM qc_report WHERE seq_folder = '%s' AND seq_id = '%s'", seqFolder, seqID) )
        writeNewData <- dbWriteTable(dbCon, name = 'qc_report', value = QC_REPORT_RES, row.names = FALSE, append = TRUE)
        #----------------------------------------------------------------------#
        dbDisconnect(dbCon)
        #----------------------------------------------------------------------#

        #==========================================================================================#
            write(sprintf("%s |  RESULTS DATABASE IMPORTED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message(sprintf("%s |  RESULTS DATABASE IMPORTED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")))
        #==========================================================================================#
    }else{
        #==========================================================================================#
            write(sprintf("%s | RESULTS DATABASE IMPORT SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message(sprintf("%s | RESULTS DATABASE IMPORT SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")))
        #==========================================================================================#
    }   
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | MULTIQC RESULT INTEGRATION & SUMMARY FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("MULTIQC RESULT INTEGRATION & SUMMARY FINISHED.")
#==================================================================================================# 

#==================================================================================================#
    write("   ", LOG_FILE, append=T)
    write(sprintf("%s | WES MultiQC QC-Result Summary DONE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
    write("   ", LOG_FILE, append=T)
#==================================================================================================# 


#---| COVERAGE CURVE PLOT in R --------------------------------------------------------------------#
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
#--------------------------------------------------------------------------------------------------#
