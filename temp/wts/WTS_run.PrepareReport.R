
#---/ PAKCAGES /-------------------------------------------------------------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("dplyr"))
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ ARGUMENTS /------------------------------------------------------------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"),     action="store", default="/data/wts", type="character", help="wts analysis base folder. REQUIRED."),   
        make_option(c("--SEQ_FOLDER"),   action="store", default=NA,          type="character", help="analysis batch id (seq_folder). REQUIRED."),
        make_option(c("--ADVANCED"),     action="store", default=FALSE,       type="logical",   help="include advanced report results"),
        make_option(c("--FUSION_GENES"), action="store", default=FALSE,       type="logical",   help="include fusion-genes and splice-variants results")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR        <- ARGS$BASE_DIR
    SEQ_FOLDER      <- ARGS$SEQ_FOLDER
    ADVANCED_REPORT <- ARGS$ANALYSIS_CONFIG
    FUSION_GENES    <- ARGS$FUSION_GENES
    # config file check -------------------------------------------------------#
    if( is.null(SEQ_FOLDER) ){ stop("|---!!! No SeqFolderID. REQUIRD. please check again. STOPPED.") } 
#----------------------------------------------------------------------------------------------------------------------------------------------------#

# manual params ................................................................
    BASE_DIR        = "/data/wts"
    SEQ_FOLDER      = "WTS_25_02"
    ADVANCED_REPORT = TRUE
    FUSION_GENES    = FALSE
#...............................................................................

#---/ GET INFO FROM DATABASE /-----------------------------------------------------------------------------------------------------------------------#
    source("/data/wts/params/ruo_wts_db.R")
    dbCon      <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_INFO )
    OrderInfo  <- dbGetQuery( dbCon, sprintf("SELECT * FROM seqfolder_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    SampleInfo <- dbGetQuery( dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    dbDisconnect(dbCon)
    #---------------------------------------------------------------------------
    dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_DATA )
    ainfo <- dbGetQuery( dbCon, sprintf("SELECT seq_folder,analysis_id FROM deg_analysis_info WHERE seq_folder = '%s'", SEQ_FOLDER)) %>% unique()
    dbDisconnect(dbCon)
    #---------------------------------------------------------------------------
    if( nrow(OrderInfo) < 1 ){ stop("|---!!! No OrderInfo of SeqFolderID found. check again. STOPPED.") }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ DATA SOURCE FOLDER & TARGET FOLDERS /----------------------------------------------------------------------------------------------------------#
    # data source folders ------------------------------------------------------
    base_data_folder <- sprintf("%s/%s/analysis", BASE_DIR, SEQ_FOLDER)
    meta_data_folder <- sprintf("%s/%s/meta",     BASE_DIR, SEQ_FOLDER)
    data_dir_clustering <- sprintf("%s/sample_clustering", base_data_folder)
    data_dir_deg        <- sprintf("%s/deg_analysis",      base_data_folder)
    if( ADVANCED_REPORT )
    {
        data_dir_pathway <- sprintf("%s/pathway_analysis",   base_data_folder)
        data_dir_gsea    <- sprintf("%s/GSEA",               base_data_folder)
        data_dir_go      <- sprintf("%s/go_enrich_analysis", base_data_folder)
    }
    # target folders -----------------------------------------------------------
    report_base_dir <- "/storage/home/kangsm/shinyWeb/REPORT_WTS"
    report_resource <- "/storage/home/kangsm/shinyWeb/resources"
    
    ReportFolder  <- sprintf("%s/%s.%s", report_base_dir, OrderInfo$ngs_order_id, SEQ_FOLDER)
    ReportDataDir <- sprintf("%s/data_figures", ReportFolder)
    if( !dir.exists(ReportFolder) )
    { 
        system(sprintf("mkdir -p %s", ReportFolder) )     
        system(sprintf("mkdir -p %s", ReportDataDir))     
        system(sprintf("ln -s %s %s/resources", report_resource, ReportFolder))
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ PREPARE RMD AND CONFIG /-----------------------------------------------------------------------------------------------------------------------#
    # main Rmd file ------------------------------------------------------------
    source_rmd_file <- sprintf("%s/DEFAULT_wts_analysis_report.Rmd", report_resource)
    report_rmd_file <- sprintf("%s/%s.WTS_Analysis_Report.Rmd",      ReportFolder, SEQ_FOLDER)
    # copy
    system(sprintf("cp %s %s", source_rmd_file, report_rmd_file))
    #system(sprintf("sed -i 's/DEFAULT_ORDER_ID/%s/g' %s", "WTS_MAQC", report_rmd_file))
    # modification
    system(sprintf("sed -i 's/DEFAULT_ORDER_ID/%s/g' %s", OrderInfo$ngs_order_id, report_rmd_file))

    # config yaml file ---------------------------------------------------------
    report_config_yml <- sprintf("%s/%s_wts_report_config.yaml", meta_data_folder, SEQ_FOLDER)
    if( !file.exists(report_config_yml) )
    { 
        stop("|---!!! No report config yaml file found. REQUIRED. check again. STOPPED.") 
    }else{
        system(sprintf("cp %s %s/report_config.yaml", report_config_yml, ReportFolder))
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ COPY FIGURES TO REPORT DATA FOLDER /-----------------------------------------------------------------------------------------------------------#

    #---/ PCA /----------------------------------------------------------------#
    source_pca_fig <- sprintf( "%s/%s.PCA.plot.png", data_dir_clustering, SEQ_FOLDER )
    target_pca_fig <- sprintf( "%s/pca_result.png", ReportDataDir )
    system(sprintf("cp %s %s", source_pca_fig, target_pca_fig))
    #--------------------------------------------------------------------------#

    #---/ UMAP /---------------------------------------------------------------#
    source_umap_3d_fig <- sprintf( "%s/%s.UMAP.3D.plot.png", data_dir_clustering, SEQ_FOLDER )
    target_umap_3d_fig <- sprintf( "%s/umap_3d_result.png", ReportDataDir )
    source_umap_2d_fig <- sprintf( "%s/%s.UMAP.2D.plot.png", data_dir_clustering, SEQ_FOLDER )
    target_umap_2d_fig <- sprintf( "%s/umap_2d_result.png", ReportDataDir )
    system(sprintf("cp %s %s", source_umap_3d_fig, target_umap_3d_fig))
    system(sprintf("cp %s %s", source_umap_2d_fig, target_umap_2d_fig))
    #--------------------------------------------------------------------------#

    #---/ SAMPLE DISTANCE /----------------------------------------------------#
    source_sample_dist_fig <- sprintf( "%s/%s.Sample.Distance.Hetamap.png", data_dir_clustering, SEQ_FOLDER )
    target_sample_dist_fig <- sprintf( "%s/sample_relative_dist_heatmap.png", ReportDataDir )
    system(sprintf("cp %s %s", source_sample_dist_fig, target_sample_dist_fig))
    #--------------------------------------------------------------------------#

    #---/ DEG RESULT /---------------------------------------------------------#
    # deg heatmap
    source_deg_heatmap <- sprintf("%s/%s.DEG.only.scales.row.heatmap.png", data_dir_deg, SEQ_FOLDER )
    target_deg_heatmap <- sprintf("%s/deg.only.heatmap.png", ReportDataDir )
    system(sprintf("cp %s %s", source_deg_heatmap, target_deg_heatmap)) 
   
    # volcano and corr.plot
    for( k in 1:nrow(ainfo) )
    {
        # volcano
        source_deg_volcano <- sprintf("%s/%s.%s.DEG.analysis.volcano.plot.png", data_dir_deg, SEQ_FOLDER, ainfo[k, "analysis_id"] )
        target_deg_volcano <- sprintf("%s/%s.deg.volcano.png", ReportDataDir, ainfo[k, "analysis_id"] )
        system(sprintf("cp %s %s", source_deg_volcano, target_deg_volcano)) 
        # corr 
        source_deg_corr     <- sprintf("%s/%s.%s.CPM.correlation.scatter.pair.plot.png", data_dir_deg, SEQ_FOLDER, ainfo[k, "analysis_id"] )
        target_deg_corr_tmp <- sprintf("%s/%s.deg.corr.tmp.png",  ReportDataDir, ainfo[k, "analysis_id"] )
        target_deg_corr     <- sprintf("%s/%s.deg.corr.plot.png", ReportDataDir, ainfo[k, "analysis_id"] )
        system(sprintf("cp %s %s", source_deg_corr, target_deg_corr_tmp)) 
        # deg plots (merge)
        deg_plot_name <- sprintf("%s/%s.deg.plot.png", ReportDataDir, ainfo[k, "analysis_id"] )
        system(sprintf("convert %s -resize %s %s", target_deg_corr_tmp, "75%", target_deg_corr))
        system(sprintf("convert +append %s %s %s", target_deg_volcano, target_deg_corr, deg_plot_name))
    }
    #--------------------------------------------------------------------------#

    if( ADVANCED_REPORT )
    {
        analysis_config_file <- sprintf("%s/%s_wts_analysis_config.yaml", meta_data_folder, SEQ_FOLDER)
        #----------------------------------------------------------------------#
        if( !file.exists(analysis_config_file) )
        {
            stop("|---!!! No WTS analysis config file found. REQUIRED for ADVANCED REPORT data preparation. STOPPED.")
        }else{
            analysis_config <- yaml::read_yaml(analysis_config_file)
            pathway_id_list <- unlist(strsplit(analysis_config$ORA$SELECTED_PATHWAY_IDS, ","))
        }

        #---/ PATHWAY ANALYSIS /-----------------------------------------------#
        for( pid in  pathway_id_list )
        {
            for( k in 1:nrow(ainfo) )
            {
                # top rank dotplot
                source_ora_rank <- sprintf("%s/%s.%s.pathway.%s.top.rank.10.merged.dotplot.png", data_dir_pathway, SEQ_FOLDER, ainfo[k, "analysis_id"], pid )
                target_ora_rank <- sprintf("%s/%s.top10.ORA.%s.dotplot.png", ReportDataDir, ainfo[k, "analysis_id"], pid )
                system(sprintf("cp %s %s", source_ora_rank, target_ora_rank))
                # top rank network
                source_ora_net <- sprintf("%s/%s.%s.ORA.%s.network.png", data_dir_pathway, SEQ_FOLDER, ainfo[k, "analysis_id"], pid )
                target_ora_net <- sprintf("%s/%s.ORA.%s.network.png", ReportDataDir, ainfo[k, "analysis_id"], pid )
                system(sprintf("cp %s %s", source_ora_net, target_ora_net))
            }
        }
        #----------------------------------------------------------------------#

        #---/ GSEA /-----------------------------------------------------------#
        for( pid in pathway_id_list )
        {
            for( k in 1:nrow(ainfo) )
            {
                target_gsea_prefix <- sprintf("%s/%s.GSEA.%s.rank", ReportDataDir, ainfo[k, "analysis_id"], pid)

                for( i in 1:5 )
                {
                    source_gsea_plot_plus <- sprintf("%s/%s.%s.GSEA.%s.NES.plus.rank.%s.plot.png", data_dir_gsea, SEQ_FOLDER, ainfo[k, "analysis_id"], pid, i )
                    target_gsea_plot_plus <- sprintf("%s.NES.plus.%s.png", target_gsea_prefix, i)
                    system(sprintf("cp %s %s", source_gsea_plot_plus, target_gsea_plot_plus))
                    system(sprintf("convert %s -bordercolor white -border 20x10 %s", target_gsea_plot_plus, target_gsea_plot_plus))

                    source_gsea_plot_minus <- sprintf("%s/%s.%s.GSEA.%s.NES.minus.rank.%s.plot.png", data_dir_gsea, SEQ_FOLDER, ainfo[k, "analysis_id"], pid, i )
                    target_gsea_plot_minus <- sprintf("%s.NES.minus.%s.png", target_gsea_prefix, i)
                    system(sprintf("cp %s %s", source_gsea_plot_minus, target_gsea_plot_minus))
                    system(sprintf("convert %s -bordercolor white -border 20x10 %s", target_gsea_plot_minus, target_gsea_plot_minus))
                }
                system(sprintf("convert +append %s.NES.plus.1.png %s.NES.plus.2.png %s.NES.plus.3.png %s.NES.plus.4.png %s.NES.plus.5.png %s/%s.top5.NES.plus.GSEA.%s.plots.png", 
                    target_gsea_prefix, target_gsea_prefix, target_gsea_prefix, target_gsea_prefix, target_gsea_prefix, ReportDataDir, ainfo[k, "analysis_id"], pid
                ))
                system(sprintf("convert +append %s.NES.minus.1.png %s.NES.minus.2.png %s.NES.minus.3.png %s.NES.minus.4.png %s.NES.minus.5.png %s/%s.top5.NES.minus.GSEA.%s.plots.png", 
                    target_gsea_prefix, target_gsea_prefix, target_gsea_prefix, target_gsea_prefix, target_gsea_prefix, ReportDataDir, ainfo[k, "analysis_id"], pid
                ))
            }
        }
        #----------------------------------------------------------------------#

        #---/ GO ENRICH /------------------------------------------------------#
        for( pid in pathway_id_list )
        {
            for( k in 1:nrow(ainfo) )
            {
                # dotplot
                source_go_dotplot <- sprintf("%s/%s.%s.Go.enrich.top.rank.10.merged.dotplot.png", data_dir_go, SEQ_FOLDER, ainfo[k, "analysis_id"] )
                target_go_dotplot <- sprintf("%s/%s.top10.GO.dotplot.png", ReportDataDir, ainfo[k, "analysis_id"] )
                system(sprintf("cp %s %s", source_go_dotplot, target_go_dotplot))
                # top net
                source_go_dotplot <- sprintf("%s/%s.%s.GO.enrich.gobp.top.rank.10.up.down.degs.network.png", data_dir_go, SEQ_FOLDER, ainfo[k, "analysis_id"] )
                target_go_dotplot <- sprintf("%s/%s.top10.GO.network.png", ReportDataDir, ainfo[k, "analysis_id"] )
                system(sprintf("cp %s %s", source_go_dotplot, target_go_dotplot))
            }
        }
        #----------------------------------------------------------------------#
    }
    #--------------------------------------------------------------------------#

    if( FUSION_GENES )
    {
        # load fusion-gene and splice-variant data from database ---------------
        source("/storage/home/kangsm/runScripts/WTS_Fun.ReportModules.R")
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_DATA )
        fg_data <- dbGetQuery(dbCon, sprintf("SELECT * FROM fusion_genes WHERE seq_folder = '%s'", SEQ_FOLDER))
        sv_data <- dbGetQuery(dbCon, sprintf("SELECT * FROM splice_variants WHERE seq_folder = '%s'", SEQ_FOLDER))
        dbDisconnect(dbCon)
        # load fusion gene related disease info -------------------------------- 
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_RESOURCE )
        fusion_gene_disease <- dbGetQuery( dbCon, "SELECT * FROM fusion_gene_related_disease WHERE validated = '1'")
        dbDisconnect(dbCon)
        # filter fusion-gene result --------------------------------------------
        FusionGeneResTable <- WTS_report.FusionGenes.Table(
            SeqFolderId      = SEQ_FOLDER, 
            SampleInfo       = SampleInfo, 
            FusionGeneResult = fg_data
        )
        FusionGeneResTable <- FusionGeneResTable %>% mutate( fgindex = paste0( geneA, "--", geneB ) ) 
        # add disease info and modification for report -------------------------
        FusionGeneResTable <- data.frame(
            seq_folder = SEQ_FOLDER,
            FusionGeneResTable,
            fusion_gene_disease[match(FusionGeneResTable$fgindex, fusion_gene_disease$fusionGene), c(5,7,8,9,10,11)]
        )
        {
            FusionGeneResTable$geneA_breakpoint <- gsub( "\\:", "_", FusionGeneResTable$geneA_breakpoint )
            FusionGeneResTable$geneA_breakpoint <- gsub( "^", "chr", FusionGeneResTable$geneA_breakpoint )
            FusionGeneResTable$geneB_breakpoint <- gsub( "\\:", "_", FusionGeneResTable$geneB_breakpoint )
            FusionGeneResTable$geneB_breakpoint <- gsub( "^", "chr", FusionGeneResTable$geneB_breakpoint )
            FusionGeneResTable$geneA_strand <- sapply(FusionGeneResTable$geneA_strand, function(st) unlist(strsplit(st, "\\/"))[1] )
            FusionGeneResTable$geneB_strand <- sapply(FusionGeneResTable$geneB_strand, function(st) unlist(strsplit(st, "\\/"))[1] )
            FusionGeneResTable$non_cancer_disease <- gsub(" ", "_", FusionGeneResTable$non_cancer_disease)
            FusionGeneResTable$non_cancer_disease <- gsub(";", " ", FusionGeneResTable$non_cancer_disease)
            FusionGeneResTable[which(is.na(FusionGeneResTable$non_cancer_disease)), "non_cancer_disease"] <- ""
            FusionGeneResTable$cancer_related     <- sapply(FusionGeneResTable$cancer_related,     function(w) ifelse( is.na(w) , "", ifelse( w == 1, "O", "" ) ))
            FusionGeneResTable$kinase             <- sapply(FusionGeneResTable$kinase,             function(w) ifelse( is.na(w) , "", ifelse( w == 1, "O", "" ) ))
            FusionGeneResTable$oncogene           <- sapply(FusionGeneResTable$oncogene,           function(w) ifelse( is.na(w) , "", ifelse( w == 1, "O", "" ) ))
            FusionGeneResTable$tumorSuppressor    <- sapply(FusionGeneResTable$tumorSuppressor,    function(w) ifelse( is.na(w) , "", ifelse( w == 1, "O", "" ) ))
            FusionGeneResTable$transcrptionFactor <- sapply(FusionGeneResTable$transcrptionFactor, function(w) ifelse( is.na(w) , "", ifelse( w == 1, "O", "" ) ))
        }
        # import into database -------------------------------------------------
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_DATA )
        DeletePreData <- dbGetQuery(dbCon, sprintf("DELETE FROM fusion_genes_report WHERE seq_folder = '%s'", SEQ_FOLDER))
        WriteNewData  <- dbWriteTable(dbCon, name="fusion_genes_report", value=FusionGeneResTable, row.names=FALSE, append=TRUE)
        dbDisconnect(dbCon)
        # copy plot into report data folder ------------------------------------
        for( k in 1:nrow(FusionGeneResTable) )
        {
            fg_data_dir <- sprintf("%s/%s/%s_rep1/fusionGenes_spliceVariants/FusionGenesPlots", BASE_DIR, SEQ_FOLDER, FusionGeneResTable[k, "sample_group"])
            file_prefix <- ifelse( FusionGeneResTable[k, "analysis_method"] == "arriba", "arriba.only.FusionGene", "FusionGene" )

            source_fusion_img <- sprintf("%s/%s_rep1.%s.%s-%s.plot.png", fg_data_dir, FusionGeneResTable[k, "sample_group"], file_prefix, FusionGeneResTable[k, "geneA"], FusionGeneResTable[k, "geneB"])
            target_fusion_img <- sprintf("%s/%s_fusion_%s_%s.png", ReportDataDir, SEQ_FOLDER, FusionGeneResTable[k, "sample_group"], FusionGeneResTable[k, "index"] )

            system(sprintf("convert '%s' -resize %s %s", source_fusion_img, "80%", target_fusion_img))
        }
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
  
    


    
