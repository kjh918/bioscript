
# Prepare Personalized NGS Report (PNR)

#---| PACKAGES |-----------------------------------------------------------------------------------#
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("data.table"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("Cairo"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("superFreq"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(stringsAsFactors=FALSE) 
    options(bitmapType='cairo')
#--------------------------------------------------------------------------------------------------#

#---| PARSE ARGUMENTS |----------------------------------------------------------------------------#
    option_list = list( 
      make_option(c("--BASE_DIR"),  action="store", default=NA, type="character", help="data base dir"),   
      make_option(c("--SEQ_ID")  ,  action="store", default=NA, type="character", help="seq id"),
      make_option(c("--PGSC_DIR"),  action="store", default=NA, type="character", help="pgsc dir (pgsc run dir)"),
      make_option(c("--IMPORT_DB"), action="store", default=TRUE, type="logical", help="Rimport results into db")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR  <- ARGS$BASE_DIR
    SEQ_ID    <- ARGS$SEQ_ID
    PGSC_DIR  <- ARGS$PGSC_DIR
    IMPORT_DB <- ARGS$IMPORT_DB

    # BASE_DIR = "/data/pon/pon_wes_gcx.03"
    # SEQ_ID   = "PON03_02"
    # PGSC_DIR = "/storage/home/kangsm/PGS_calculation/20240920_han"
#--------------------------------------------------------------------------------------------------#

#---| SET-UP |-------------------------------------------------------------------------------------#
    pnr_dir = sprintf("%s/%s/pnr", BASE_DIR, SEQ_ID)
    if( !dir.exists(pnr_dir) ){ system(sprintf("mkdir -p %s", pnr_dir)) }
    #--------------------------------------------------------------------------#
    vcf_dir      = sprintf("%s/%s/vcf", BASE_DIR, SEQ_ID)
    sf_cnv_dir   = sprintf("%s/%s/cnv/superfreq", BASE_DIR, SEQ_ID)
    pgsc_res_dir = sprintf("%s/results", PGSC_DIR)
#--------------------------------------------------------------------------------------------------#


#---| 1. PGSC Results |----------------------------------------------------------------------------#
    pgs_list_file = "/storage/home/kangsm/PGS_calculation/GCX_PGS_v2.list"
    pgs_list_rev = as.data.frame(fread(pgs_list_file, header=F)) %>% 
        dplyr::rename(pgs_id=1, setname=2, setname_kor=3)
    ExcludePhenotypes = c("AlcoholDrinkingStatus","NicotineDependence","AnginaPectoris")
    pgs_list_rev = pgs_list_rev %>% filter( setname %nin% ExcludePhenotypes )
    #--------------------------------------------------------------------------#
    IMPUTE_PGS_PLOT_LIST = IMPUTE_PGS_PLOT_LIST2 = list()
    IMPUTE_PGS_RESULT    = data.frame()
    # pgsc result integration and summary -------------------------------------#
    message(" >> PGSC result integration and summary.")
    #
    for( n in 1:nrow(pgs_list_rev) )
    {
        pgsc_run_id = sprintf("PGS-%s", pgs_list_rev[n,2])
        pgsc_res_file = sprintf("%s/%s/score/%s_pgs.txt.gz", pgsc_res_dir, pgsc_run_id, pgsc_run_id)
        if( !file.exists(pgsc_res_file) )
        {
            next
        }else{
            pgsc_res    = fread(sprintf("%s/%s/score/%s_pgs.txt.gz", pgsc_res_dir, pgsc_run_id, pgsc_run_id))
            pgsc_res    = pgsc_res %>% arrange(Z_norm2) %>% mutate( pRank = round((1:nrow(pgsc_res)) / nrow(pgsc_res) *100, 1) )
            patient_res = pgsc_res %>% filter( IID == SEQ_ID ) %>% dplyr::select(c(1,3,4,8,10)) %>% 
                dplyr::rename(sampleset=1, seq_id=2, score_file=3, rank_score=4) %>%
                mutate( trait = pgs_list_rev[n,2], pgs_catalog_id = pgs_list_rev[n,1] )
            ###
            if( patient_res$pRank > 80 ) { 
                risk_status = "High Risk"
            }else if( patient_res$pRank < 40  ) {
                risk_status = "Low Risk"
            }else{ risk_status = "Average Risk" }
            ###
            patient_res$risk_level = risk_status
            ###
            patient_res = patient_res[,c(2,6,7,3,4,5,8)]
            TAIL = ifelse(patient_res$rank_score > 0, FALSE, TRUE)
            patient_res$pvalue=pnorm( patient_res$rank_score, mean = mean(pgsc_res$Z_norm2), sd = sd(pgsc_res$Z_norm2), lower.tail=TAIL )
            ###
            IMPUTE_PGS_RESULT = rbind(IMPUTE_PGS_RESULT, patient_res)
            # Plots ---------------------------------------------------------------#
            maxScore = max(abs(pgsc_res$Z_norm2))
            patient_location = pgsc_res %>% filter( IID == SEQ_ID ) %>% .$Z_norm2
            pat_locate_col = ifelse( patient_res$risk_level == "High Risk", "#B85450", "#82B366" )
            pgs_plot = ggplot(pgsc_res) + geom_density(aes(x=Z_norm2), col="#6C8EBF", linewidth=0.95, alpha=0.8 ) +
                geom_vline( xintercept=patient_location, col=pat_locate_col, linewidth=2.5, alpha=0.7) + 
                coord_cartesian(xlim=c(-maxScore, maxScore)) +
                ggtitle( pgs_list_rev[n,2] ) +
                labs(x="", y="") +
                theme(
                    panel.background = element_rect(fill='white'), 
                    panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'), 
                    legend.position  = "none",
                    axis.line.x  = element_line(colour = "#6f6d6d"), 
                    axis.line.y  = element_line(colour = "#6f6d6d"),
                    axis.text.x  = element_text(face = "bold", size = 12),
                    axis.text.y  = element_text(face = "bold", size = 10),
                    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 8)),
                    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 8)),
                    plot.title   = element_text(face = "bold", size=15),
                    plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
                )
            IMPUTE_PGS_PLOT_LIST = c(IMPUTE_PGS_PLOT_LIST, list(pgs_plot))
            names(IMPUTE_PGS_PLOT_LIST)[length(IMPUTE_PGS_PLOT_LIST)] = pgs_list_rev[n,2]
            ###
            pgs_plot2 = ggplot(pgsc_res) + geom_density(aes(x=Z_norm2), col="#6C8EBF", linewidth=0.95, alpha=0.8 ) +
                geom_vline( xintercept=patient_location, col=pat_locate_col, linewidth=2.5, alpha=0.7) + 
                coord_cartesian(xlim=c(-maxScore, maxScore)) +
                ggtitle( pgs_list_rev[n,2] ) +
                labs(x="", y="") +
                theme(
                    panel.background = element_rect(fill='white'), 
                    panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'), 
                    legend.position  = "none",
                    axis.line.x  = element_line(colour = "#6f6d6d"), 
                    axis.line.y  = element_line(colour = "#6f6d6d"),
                    axis.text.x  = element_text(face = "bold", size = 9),
                    axis.text.y  = element_text(face = "bold", size = 8),
                    axis.title.x = element_text(face = "bold", size = 9, margin = margin(t = 8)),
                    axis.title.y = element_text(face = "bold", size = 9, margin = margin(r = 8)),
                    plot.title   = element_text(face = "bold", size=10),
                    plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
                )
            IMPUTE_PGS_PLOT_LIST2 = c(IMPUTE_PGS_PLOT_LIST2, list(pgs_plot2))
            names(IMPUTE_PGS_PLOT_LIST2)[length(IMPUTE_PGS_PLOT_LIST2)] = pgs_list_rev[n,2]
            ###
            cat(sprintf("%s...", n))
        }
    }
    #
    cat("\n")
    message(" >> create PGSC result tables.")
    # save pgsc results
    IMPUTE_PGS_RESULT = IMPUTE_PGS_RESULT %>% arrange(pRank) %>% 
        mutate( trait_kor = pgs_list_rev[match(trait, pgs_list_rev$setname), "setname_kor"])
    write.table(IMPUTE_PGS_RESULT, sprintf("%s/%s.PGSC.results.tsv", pnr_dir, SEQ_ID), quote=F, col.names=T, row.names=F ,sep="\t")
    openxlsx::write.xlsx(IMPUTE_PGS_RESULT, sprintf("%s/%s.PGSC.results.xlsx", pnr_dir, SEQ_ID))
    #
    IMPUTE_PGS_PLOT_LIST  = IMPUTE_PGS_PLOT_LIST[IMPUTE_PGS_RESULT$trait]
    IMPUTE_PGS_PLOT_LIST2 = IMPUTE_PGS_PLOT_LIST2[IMPUTE_PGS_RESULT$trait]
    #
    message(" >> create PGSC result plots.")
    #
    png(sprintf("%s/%s.PGSC.result.plots.png", pnr_dir, SEQ_ID), width=25,height=30,units="in",res=150)
    cowplot::plot_grid(plotlist=IMPUTE_PGS_PLOT_LIST, ncol=4)
    dev.off()
    #
    save(IMPUTE_PGS_PLOT_LIST2, file=sprintf("%s/%s.PGSC.plot.list.Rdata", pnr_dir, SEQ_ID))
    #
    if( IMPORT_DB )
    {
        message(" >> PGSC results import into DB.")
        require(RMySQL)
        dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db="gcx_pon" )
        clearing = dbGetQuery(dbCon, sprintf("DELETE FROM %s WHERE seq_id= '%s'", "pgsc", SEQ_ID))
        dbWriteTable(dbCon, name="pgsc", value=IMPUTE_PGS_RESULT, row.names=F, append=T)
        dbDisconnect(dbCon)
    }
#--------------------------------------------------------------------------------------------------#

#---| 2. superFreq CNV |---------------------------------------------------------------------------#
    #
    message(" >> CNV results summary.")
    source("/storage/home/kangsm/myScripts/R_FUNCTIONAL_MODULES/FunMod_Customized_superFreq_PlotCR.R")
    #
    sf_res_rdata = sprintf("%s/%s.TumorOnly/clusters.Rdata", sf_cnv_dir, SEQ_ID)
    load( sf_res_rdata ) # clusters
    # cnv plot 
    message(" >> create CNV plots.")
    png(sprintf("%s/%s.superfreq.cnv.plot.png", pnr_dir, SEQ_ID), width=18, height=4, unit='in',,res=150)
    plotCR.custom( cR=clusters[[1]]$CR, chr='all', genome='hg19', alpha=0.1, pt.cex=0.8, colourDeviation=T, lwd=1, sideSpace=0.05, noSex=T)
    dev.off()
    # cnv genes table
    message(" >> create CNV genes table.")
    #
    sf_cnv_res_file = sprintf("%s/%s.TumorOnly/%s.TumorOnly/data/CNAsegments_%s.tsv", sf_cnv_dir, SEQ_ID, SEQ_ID, SEQ_ID)
    sf_cnv_res      = read.delim(sf_cnv_res_file)
    #
    RES_GAIN = RES_LOSS = data.frame()
    # gain genes
    cn_gain = sf_cnv_res %>% filter( M > log2(1.5) )
    if( nrow(cn_gain) > 0 )
    {
        RES_GAIN = data.frame()
        for( k in 1:nrow(cn_gain))
        {
            genes_gain = unlist(strsplit(cn_gain$genes[k] ,","))
            genes_gain = genes_gain[grep("^ENST", genes_gain, invert=T)]
            genes_gain_N = length(genes_gain)
            genes_gain_list = paste(genes_gain, collapse=",")
            res_gain = data.frame( cna = "gain", chr = cn_gain$chr[k], start=cn_gain$start[k], end=cn_gain$end[k], logCR=cn_gain$M[k], genes=genes_gain_N, gene_list=genes_gain_list, cosmic_genes = cn_gain$COSMIC_genes[k])
            RES_GAIN = rbind(RES_GAIN, res_gain)
        }
    }else{
        res_gain = data.frame( cna = "gain", chr = "", start=NA, end=NA, logCR=NA, genes=NA, gene_list="", cosmic_genes="")
        RES_GAIN = rbind(RES_GAIN, res_gain)
    }
    # loss genes
    cn_loss = sf_cnv_res %>% filter( M < -1 )
    if( nrow(cn_loss) > 0 )
    {
        RES_LOSS = data.frame()
        for( k in 1:nrow(cn_loss))
        {
            genes_loss = unlist(strsplit(cn_loss$genes[k] ,","))
            genes_loss = genes_loss[grep("^ENST", genes_loss, invert=T)]
            genes_loss_N = length(genes_loss)
            genes_loss_list = paste(genes_loss, collapse=",")
            res_loss = data.frame( cna = "loss", chr = cn_loss$chr[k], start=cn_loss$start[k], end=cn_loss$end[k], logCR=cn_loss$M[k], genes=genes_loss_N, gene_list=genes_loss_list, cosmic_genes = cn_gain$COSMIC_genes[k])
            RES_LOSS = rbind(RES_LOSS, res_loss)
        }
    }else{
        res_loss = data.frame( cna = "loss", chr = "", start=NA, end=NA, logCR=NA, genes=NA, gene_list="", cosmic_genes="")
        RES_LOSS = rbind(RES_LOSS, res_loss)
    }   
    # integration gain and loss
    RES_CNA = data.frame( seq_id = SEQ_ID, swtool = "superFreq", rbind( RES_GAIN, RES_LOSS) )
    write.table(RES_CNA, sprintf("%s/%s.CNV.genes.results.tsv", pnr_dir, SEQ_ID), quote=F, col.names=T, row.names=F ,sep="\t")
    openxlsx::write.xlsx(RES_CNA, sprintf("%s/%s.CNV.genes.results.xlsx", pnr_dir, SEQ_ID))
    # import DB
    if( IMPORT_DB )
    {
        message(" >> CNV results import into DB.")
        require(RMySQL)
        dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db="gcx_pon" )
        clearing = dbGetQuery(dbCon, sprintf("DELETE FROM %s WHERE seq_id= '%s'", "cnv", SEQ_ID))
        dbWriteTable(dbCon, name="cnv", value=RES_CNA, row.names=F, append=T)
        dbDisconnect(dbCon)        
    }
#--------------------------------------------------------------------------------------------------#

#---| 3. Germline Variants Stats |-----------------------------------------------------------------#
    #
    message(" >> Germline Variants Stats summary.")
    #
    gvcf_stat_file = sprintf("%s/%s.deepvariant.germline.variant.filtered.maf.variant.class.summary.txt", vcf_dir, SEQ_ID)
    if( !file.exists(gvcf_stat_file) )
    {
        write("No Germline Variants Stats File.", sprintf("%s/%s.no.germmline.variants.stats.files", pnr_dir, SEQ_ID))
    }else{
        GV_stats = read.delim(gvcf_stat_file)
        GV_stats$variant_call_mode = "germline"
        require(RMySQL)
        dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db="gcx_pon" )
        clearing = dbGetQuery(dbCon, sprintf("DELETE FROM %s WHERE seq_id= '%s' AND variant_call_mode = '%s'", "variants_summary", SEQ_ID, "germline"))
        dbWriteTable(dbCon, name="variants_summary", value=GV_stats, row.names=F, append=T)
        dbDisconnect(dbCon)
        write("Germline Variants Stats Imported Into DB.", sprintf("%s/%s.germmline.variants.db.imported", pnr_dir, SEQ_ID))
    }
#--------------------------------------------------------------------------------------------------#

#---| 4. PNR Report RMD Prepare |------------------------------------------------------------------#
    #
    message(" >> PREPARE PNR-REPORT DIR and DATA.")
    #
    pnr_report_base_dir = "/storage/home/kangsm/shinyWeb/GCX_NGS_REPORT"
    pnr_report_dir = sprintf("%s/PNR_%s", pnr_report_base_dir, SEQ_ID)
    if(!dir.exists(pnr_report_dir))
    { 
        # make dir
        system(sprintf("mkdir -p %s", pnr_report_dir)) 
        # symbolic link resources 
        system(sprintf("ln -s /storage/home/kangsm/shinyWeb/resources %s/resources", pnr_report_dir))
    }
    # copy necessary files into pnr dir ---------------------------------------#
    # CNV PLOTS
    system(sprintf("cp %s/%s.superfreq.cnv.plot.png %s", pnr_dir, SEQ_ID, pnr_report_dir))
    # PGSC PLOT LIST
    system(sprintf("cp %s/%s.PGSC.plot.list.Rdata %s", pnr_dir, SEQ_ID, pnr_report_dir))
    #
    message("DONE.")
#--------------------------------------------------------------------------------------------------#

