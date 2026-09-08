
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
suppressPackageStartupMessages(library("reshape2"))
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--BASE_DIR"), action="store", default=NA, type="character", help="[ Required ] analysis basse folder"),
        make_option(c("--BATCH_ID"), action="store", default=NA, type="character", help="[ Required ] analysis batch-id"),
        make_option(c("--SEQ_ID"), action="store", default=NA, type="character", help="[ Required ] sample-id"),
        make_option(c("--BAM_TYPE"), action="store", default=NA, type="character", help="bam file type. 'hmftools' or 'gatk'. default = 'hmftools'"),
        make_option(c("--INCLUDE_FASTQ_SCREEN"), action="store", default=TRUE, type="logical", help="include fastq-screen result"),
        make_option(c("--DATABASE_IMPORT") , action="store", default=FALSE, type="logical", help="database import")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BaseDir            <- ARGS$BASE_DIR
    BatchID            <- ARGS$BATCH_ID
    SeqID              <- ARGS$SEQ_ID
    BamType            <- ARGS$BAM_TYPE
    IncludeFastqScreen <- as.logical(ARGS$INCLUDE_FASTQ_SCREEN)
    DatabaseImport     <- as.logical(ARGS$DATABASE_IMPORT)
#--------------------------------------------------------------------------------------------------#

### CHECK PARAMS
    if( is.null(IncludeFastqScreen) ){ IncludeFastqScreen <- FALSE      }
    if( is.null(DatabaseImport)     ){ DatabaseImport     <- FALSE      }
    if( is.null(BamType)            ){ BamType            <- "hmftools" }
###
### manual params
    # BaseDir            = "/data/wgs"
    # BatchID            = "WGS-26-01"
    # SeqID              = "WGS-26-01_01"
    # BamType            = "hmftools"
    # IncludeFastqScreen = FALSE
    # DatabaseImport     = TRUE
###
### FUNCTIONS
    rounding <- function(x, digits = 0) {
        posneg <- sign(x)
        z <- abs(x) * 10^digits
        z <- floor(z + 0.5)
        v <- z / 10^digits * posneg
        return(v)
    }
    Calculate.Uniformity <- function( GATK_WGS_MATRIC ){
        unif_target_dp <- GATK_WGS_MATRIC$MEAN_COVERAGE * 0.2
        cov_tb <- GATK_WGS_MATRIC %>% dplyr::select( 
            all_of(colnames(GATK_WGS_MATRIC)[colnames(GATK_WGS_MATRIC) %in% paste0("PCT_", 1:100, "X")])
        ) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("dp") %>% mutate(
            dp  =gsub("^PCT_", "", dp),
            dp = gsub("X$", "", dp),
            dp = as.numeric(dp)
        ) %>% arrange(dplyr::desc(dp)) %>% dplyr::rename( cov=2 )
        bin_min <- cov_tb %>% filter( dp < unif_target_dp ) %>% .$dp %>% max
        bin_max <- cov_tb %>% filter( dp > unif_target_dp ) %>% .$dp %>% min
        unif_cov <- approx( 
            x = c( bin_min, bin_max ),
            y = c( cov_tb %>% filter( dp == bin_min ) %>% .$cov, cov_tb %>% filter( dp == bin_max ) %>% .$cov ),
            xout = unif_target_dp
        )$y
        return(unif_cov)
    }
###
### QC-DATA PREFIX
    if( BamType == "hmftools" ){
        QcDataPrefix <- sprintf("%s/%s/%s/qcfiles/%s.redux", BaseDir, BatchID, SeqID, SeqID)
    }else{
        QcDataPrefix <- sprintf("%s/%s/%s/qcfiles/%s.gatk", BaseDir, BatchID, SeqID, SeqID)
    }
###
### READ QC-DATA    
    QcDataDir <- sprintf("%s.qc_res/%s_qc_data", QcDataPrefix, SeqID)
    
    mqc_stats_file  <- sprintf("%s/multiqc_general_stats.txt", QcDataDir)
    mqc_dups_file   <- sprintf("%s/multiqc_picard_dups.txt", QcDataDir)
    mqc_align_file  <- sprintf("%s/multiqc_picard_AlignmentSummaryMetrics.txt", QcDataDir)
    mqc_insize_file <- sprintf("%s/multiqc_picard_insertSize.txt", QcDataDir)
    mqc_wgsmetric_file <- sprintf("%s/multiqc_picard_wgsmetrics.txt", QcDataDir)

    if( !file.exists(mqc_stats_file) ){
        message("No MultiQC result files. Please Run MultiQC first."); quit( save="no", status=0 ) 
    }else{
        mqc_stats <- read.delim(mqc_stats_file)
        mqc_dups  <- read.delim(mqc_dups_file)
        mqc_align <- read.delim(mqc_align_file)
        mqc_insertsize <- read.delim(mqc_insize_file)
        mqc_wgsmetric  <- read.delim(mqc_wgsmetric_file)
    }

    cov_mosdp_file <- sprintf("%s/mosdepth_perchrom.txt", QcDataDir)
    if( !file.exists(cov_mosdp_file) ){ message("No MultiQC result files. Please Run MultiQC first."); quit( save="no", status=0 ) }
    cov_mosdp <- read.delim(cov_mosdp_file, check.names=F) 
###
### QC-DATA ARRANGMENT
    # fastq 
    QC_fastq <- mqc_stats %>% filter( Sample == SeqID ) %>% dplyr::select(c(
        "Sample",
        "fastp_mqc.generalstats.fastp.after_filtering_q30_bases",
        "fastp_mqc.generalstats.fastp.after_filtering_q30_rate",
        "fastp_mqc.generalstats.fastp.after_filtering_gc_content"
    )) %>% dplyr::rename(SeqID=1, bases_q30=2, pct_q30_bases=3, gc_ratio=4 ) %>% mutate(
        pct_q30_bases = rounding((pct_q30_bases*100), digits=1),
        pct_gc        = rounding((gc_ratio*100), digits=1)
    ) %>% dplyr::select(c("SeqID","bases_q30","pct_q30_bases","pct_gc"))
    # duplicates
    QC_dups <- mqc_dups %>% mutate( SeqID = SeqID, pct_opt_dup = READ_PAIR_OPTICAL_DUPLICATES/READ_PAIRS_EXAMINED ) %>% 
        dplyr::rename( pct_dup=PERCENT_DUPLICATION ) %>% 
        mutate(
            pct_dup     = rounding((pct_dup*100), digits=1),
            pct_opt_dup = rounding((pct_opt_dup*100), digits=1)
        ) %>% dplyr::select(c( "SeqID","pct_dup","pct_opt_dup" ))
    # align
    QC_align <- mqc_align %>% mutate( SeqID = SeqID ) %>% dplyr::select(c(
        "SeqID","PF_READS","PF_READS_ALIGNED","PF_ALIGNED_BASES","PCT_PF_READS_ALIGNED",
        "PCT_READS_ALIGNED_IN_PAIRS","MEAN_ALIGNED_READ_LENGTH","PCT_HARDCLIP","PCT_SOFTCLIP","PCT_ADAPTER"
    )) %>% dplyr::rename(
        rds_pf=2, rds_align=3, bases_align=4, pct_rds_align=5,
        pct_rds_align_pair=6, mean_align_read_length=7, pct_hardclip=8, pct_soft_clip=9, pct_adapter=10
    ) %>% mutate(
        pct_rds_align = rounding((pct_rds_align*100), digits=1),
        pct_rds_align_pair = rounding((pct_rds_align_pair*100), digits=1),
        mean_align_read_length = rounding(mean_align_read_length, digits=1),
        pct_hardclip = rounding((pct_hardclip*100), digits=1),
        pct_soft_clip = rounding((pct_soft_clip*100), digits=1),
        pct_adapter = rounding((pct_adapter*100), digits=1)
    )
    # insert-size
    QC_insert <- mqc_insertsize %>% mutate( SeqID = SeqID ) %>% dplyr::select(c(
        "SeqID","MEAN_INSERT_SIZE","MEDIAN_INSERT_SIZE","MIN_INSERT_SIZE","MAX_INSERT_SIZE"
    )) %>% dplyr::rename(
        mean_insert_size=2, median_insert_size=3, min_insert_size=4, max_insert_size=5
    ) %>% mutate( mean_insert_size = rounding(as.numeric(mean_insert_size), digits=1) )
    # wgs-metric
    QC_wgs <- mqc_wgsmetric %>% mutate( SeqID = SeqID ) %>% mutate(
        uniformity = Calculate.Uniformity(.)
    ) %>% dplyr::select(c(
        "SeqID","MEAN_COVERAGE","MEDIAN_COVERAGE","uniformity","PCT_1X","PCT_5X","PCT_15X","PCT_30X","PCT_60X","PCT_80X",'PCT_100X'
    )) %>% dplyr::rename(
        cov_gatk_mean=2, cov_gatk_median=3, pct_1x=5, pct_5x=6, pct_15x=7, pct_30x=8, pct_60x=9, pct_80x=10, pct_100x=11
    ) %>% mutate(
        cov_gatk_mean = rounding(cov_gatk_mean, digits=1),
        uniformity = rounding((uniformity*100), digits=1),
        pct_1x = rounding((pct_1x*100), digits=1),
        pct_5x = rounding((pct_5x*100), digits=1),
        pct_15x = rounding((pct_15x*100), digits=1),
        pct_30x = rounding((pct_30x*100), digits=1),
        pct_60x = rounding((pct_60x*100), digits=1),
        pct_80x = rounding((pct_80x*100), digits=1),
        pct_100x = rounding((pct_100x*100), digits=1)
    )
    # chrom-level-cov
    QC_CHROM_COV <- cov_mosdp %>% dplyr::select(!c("Sample")) %>% t() %>% as.data.frame() %>% 
        tibble::rownames_to_column("chrom") %>% dplyr::rename(cov=2) %>% 
        filter( chrom != "MT" ) %>% 
        mutate( chrom = factor(chrom, levels=c(1:22, "X","Y")) ) %>% 
        arrange(chrom) %>% mutate( 
            chrom = gsub("^", "chr", chrom ),
            cov = rounding(cov,1), 
            db_value = paste(chrom,cov,sep=",")
        ) %>% .$db_value %>% paste0(collapse=";")
###
### SUMMARIZAITON
    QC_RES <- QC_fastq %>% 
        left_join(., QC_align, by="SeqID") %>% 
        left_join(., QC_dups, by="SeqID") %>% 
        left_join(., QC_insert, by="SeqID") %>% 
        left_join(., QC_wgs, by="SeqID") %>% 
        mutate( chrom_cov=QC_CHROM_COV )

    write.table( QC_RES, sprintf("%s.QC.SUMMARY.tsv", QcDataPrefix), quote=F, col.names=T, row.names=F, sep="\t")
      
###
### FASTQ-SCREEN QC
    if( IncludeFastqScreen )
    {
        # read fq-screen data
        res_fqscreen_file <- sprintf("%s/%s/%s/qcfiles/%s_R1_screen.txt", BaseDir, BatchID, SeqID, SeqID)
        if( file.exists(res_fqscreen_file) ){
            res_fqscreen <- read.delim(res_fqscreen_file, skip=1, check.names=FALSE )
            fixCols    <- c("genome","tested_reads","unmap","pct_unmap","H1_G1","pct_H1_G1","Hm_G1","pct_Hm_G1",
                            "H1_Gm","pct_H1_Gm","Hm_Gm","pct_Hm_Gm")
            pickCols   <- c("pct_H1_G1","pct_Hm_G1","pct_H1_Gm","pct_Hm_Gm","pct_unmap")
            reportCols <- c("seq_folder","seq_id","fastq","species","human.pct_H1_G1","mouse.pct_H1_G1","rat.pct_H1_G1",
                            "ecoli.pct_H1_G1","phix.pct_H1_G1","pseudomonas.pct_H1_G1", "mycoplasma.pct_H1_G1",
                            "adapters.pct_H1_G1","vectors.pct_H1_G1")
            colnames(res_fqscreen) <- fixCols
            
            QC_FQSCREEN <- data.frame(
                human       = res_fqscreen %>% filter (genome == "Human" )       %>% select(any_of(pickCols)),
                mouse       = res_fqscreen %>% filter (genome == "Mouse" )       %>% select(any_of(pickCols)),
                rat         = res_fqscreen %>% filter (genome == "Rat" )         %>% select(any_of(pickCols)),
                ecoli       = res_fqscreen %>% filter (genome == "Ecoli" )       %>% select(any_of(pickCols)),
                phix        = res_fqscreen %>% filter (genome == "PhiX" )        %>% select(any_of(pickCols)),
                pseudomonas = res_fqscreen %>% filter (genome == "Pseudomonas" ) %>% select(any_of(pickCols)),
                mycoplasma  = res_fqscreen %>% filter (genome == "Mycoplasma" )  %>% select(any_of(pickCols)),
                adapters    = res_fqscreen %>% filter (genome == "Adapters" )    %>% select(any_of(pickCols)),
                vectors     = res_fqscreen %>% filter (genome == "Vectors" )     %>% select(any_of(pickCols))
            ) %>% mutate( seq_folder=BatchID, seq_id=SeqID, species="" )
            
            sum_fqscreen <- QC_FQSCREEN %>% select(any_of(reportCols))
            colnames(sum_fqscreen) <- gsub("\\.pct_H1_G1", "_H1G1", colnames(sum_fqscreen))
            sum_fqscreen$species   <- apply(sum_fqscreen[,4:12], 1, function(z) { 
                paste(gsub("_H1G1","", colnames(sum_fqscreen)[4:12][which(z >= 10)]),collapse=",") 
            })
            QC_FQSCREEN$species <- sum_fqscreen$species

            colnames(QC_FQSCREEN) <- gsub( "pct_", "" , colnames(QC_FQSCREEN) )
            colnames(QC_FQSCREEN) <- gsub( "_"   , "" , colnames(QC_FQSCREEN) )
            colnames(QC_FQSCREEN) <- gsub( "\\." , "_", colnames(QC_FQSCREEN) )
            colnames(QC_FQSCREEN) <- gsub( "seqfolder" , "seq_folder", colnames(QC_FQSCREEN) )
            colnames(QC_FQSCREEN) <- gsub( "seqid" , "seq_id", colnames(QC_FQSCREEN) )
        
            fqscreen_output_file <- sprintf("%s/%s/%s/qcfiles/%s.FASTQ-SCREEN.QC.RESULT.txt", BaseDir, BatchID, SeqID, SeqID)
            write.table( QC_FQSCREEN, file=fqscreen_output_file,  quote=F, col.names=T, row.names=F, sep="\t" )
        }else{
            QC_FQSCREEN <- data.frame()
        }
    }
###
### IMPORT DATABASE 
    if( DatabaseImport )
    {
        # bam-qc
        qc_db_data <- QC_RES %>% mutate( BatchID = BatchID )

        dbCon <- dbConnect(dbDriver("MySQL"), host="192.168.0.34", user="gcx", port=3306, password="gencurix!!", db="gcx_wgs" )
        rm_prev_data <- dbGetQuery(dbCon, sprintf("DELETE FROM qc_report WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, SeqID) )    
        write_data   <- dbWriteTable(dbCon, name = 'qc_report', value = qc_db_data, row.names = FALSE, append = TRUE)
        dbDisconnect(dbCon)

        # fastq-screen
        if( IncludeFastqScreen ){
            if( nrow(QC_FQSCREEN) > 0 ){
                dbCon <- dbConnect(dbDriver("MySQL"), host="192.168.0.34", user="gcx", port=3306, password="gencurix!!", db="gcx_wgs" )
                rm_prev_data <- dbGetQuery(dbCon, sprintf("DELETE FROM qc_fastqscreen WHERE BatchID = '%s' AND SeqID = '%s'", BatchID, SeqID))  
                write_data   <- dbWriteTable(dbCon, name="qc_fastqscreen", value=QC_FQSCREEN, row.names=F, append=T)
                dbDisconnect(dbCon)
            }
        }
    }
###
   
#---| COVERAGE CURVE PLOT in R --------------------------------------------------------------------#
    #ggplot( cov_mosdp2, aes(x=depth, y=coverage)) + geom_line( size = 1.2, color = '#2261e6' ) +
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
    #   geom_point( x = 30, y = cov_mosdp2[cov_mosdp2$depth==30, 'coverage'], size = 3, color='#e64522') +
    #    geom_point( x = 100, y = cov_mosdp2[cov_mosdp2$depth==100, 'coverage'], size = 3, color='#e64522') +  
    #    geom_text( data=cov_mosdpepth_labels, aes(x = depth, y =coverage, label=lbl), hjust=0, nudge_x = 1, nudge_y =3  ) +
    #    ggtitle("Target Coverage (Base)")

    # ggplot(mosch, aes(x=chromosome, y=coverage)) + geom_bar(position="dodge", stat="identity", width =0.7)
#--------------------------------------------------------------------------------------------------#
