

#===| PACKAGES |=====================================================================================================================================#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("parallel"))
    suppressPackageStartupMessages(library("yaml"))
    suppressPackageStartupMessages(library("reshape2"))
    suppressPackageStartupMessages(library("openxlsx"))
    suppressPackageStartupMessages(library("dplyr"))
#====================================================================================================================================================#

#===| ARGUMENTS |====================================================================================================================================#
    option_list = list( 
        make_option(c("--BASE_DIR"),        action="store", default="/data/wts", type="character", help="wts analysis base folder. REQUIRED."),   
        make_option(c("--SEQ_FOLDER"),      action="store", default=NA,          type="character", help="analysis batch id (seq_folder). REQUIRED."),
        make_option(c("--ANALYSIS_CONFIG"), action="store", default=NA,          type="character", help="analysis config YAML file. REQUIRED."),
        make_option(c("--GENOME_ASSEMBLY"), action="store", default="hg19",      type="character", help="genome assembly version. default = hg19")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR             <- ARGS$BASE_DIR
    SEQ_FOLDER           <- ARGS$SEQ_FOLDER
    ANALYSIS_CONFIG_FILE <- ARGS$ANALYSIS_CONFIG
    GENOME_ASSEMBLY      <- ARGS$GENOME_ASSEMBLY
    # config file check -------------------------------------------------------#
    if( is.null(ANALYSIS_CONFIG_FILE) ){ stop("|---!!! wts analysis config yaml file is not found. REQUIRD. please check again. STOPPED.") } 
#====================================================================================================================================================#

#---| MANUAL PARAMS |---------------------------------------------------------------------------------------------------------------------------------
    # BASE_DIR        = "/data/wts"
    # SEQ_FOLDER      = "WTS_25_02"
    # ANALYSIS_CONFIG_FILE = "/data/wts/WTS_25_02/meta/WTS_25_02_wts_analysis_config.yaml"
    # GENOME_ASSEMBLY = "hg19"
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#===| READ CONFIG and LOAD MODULES, FUNCTIONS, DATA |================================================================================================#    
    # read config -------------------------------------------------------------#
    if( tail(unlist(strsplit(ANALYSIS_CONFIG_FILE, "\\.")),1) %in% c("yaml","yml") )
    {
        CONFIG <- read_yaml( ANALYSIS_CONFIG_FILE )
    }else{
        stop("|---!!! wts analysis config file is not yaml format. please check again. STOPPED.")
    }
    # load modules and functions ----------------------------------------------#
    source("/storage/home/kangsm/runScripts/WTS_Fun.AnalysisModules.R")      # modules
    source("/storage/home/kangsm/runScripts/WTS_Fun.PlotModules.R")          # functions
    source("/data/wts/params/ruo_wts_db.R")  # database connection config
    # define function ---------------------------------------------------------#
    indexing = function(x){ sapply(unique(x), function(z) list(z)) }
    # load gene info table ----------------------------------------------------#
    dbCon     <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db="public_database" )
    GENE_INFO <- dbGetQuery(dbCon, "SELECT symbol,name,locus_group,locus_type,entrez,ens_geneid, hgnc_id FROM gene_info_hgnc") %>% unique()
    dbDisconnect(dbCon)
#====================================================================================================================================================#

#===| ANALYSIS BASE FOLDERS |========================================================================================================================#
    BaseAnalysisDir <- sprintf("%s/%s/analysis", BASE_DIR, SEQ_FOLDER)
    if( !dir.exists(BaseAnalysisDir) ){ system(sprintf("mkdir -p %s", BaseAnalysisDir)) }
#====================================================================================================================================================#

#===[ SAMPLE INFORMATION ]===========================================================================================================================#
    if( CONFIG$SAMPLE_INFO_DATA_SOURCE$LOAD_SAMPLE_INFO_FROM_DB )
    {
        # sample info load from database --------------------------------------#
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_INFO )
        SINFO <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SEQ_FOLDER) )    
        dbDisconnect(dbCon)
        # sample id and sample groups -----------------------------------------#
        SeqIdList        <- SINFO$seq_id
        SeqIdGroup       <- SINFO$sample_group
        names(SeqIdList) <- SeqIdGroup   
        message("|--->> sample info loaded from database.")
    }else{
        if( is.null(CONFIG$SAMPLE_INFO_DATA_SOURCE$SEQ_ID_LIST) )
        { 
            stop("|---!!! no sample id (seq_id) list. REQUIRED. please check again. STOPPED.") 
        }else{
            # sample id manual input ------------------------------------------#
            SeqIdList <- CONFIG$SAMPLE_INFO_DATA_SOURCE$SEQ_ID_LIST
            # sample group manual input ---------------------------------------#
            if( !is.null(CONFIG$SAMPLE_INFO_DATA_SOURCE$SEQ_GROUP_LIST) )
            {   
                names(SeqIdList) <- CONFIG$SAMPLE_INFO_DATA_SOURCE$SEQ_GROUP_LIST 
            }else{ 
                names(SeqIdList) <- CONFIG$SAMPLE_INFO_DATA_SOURCE$SEQ_ID_LIST 
                message("|--->> no sample groups found. sample id will be used as sample group.")
            }
            # sample info table creation --------------------------------------#
            SINFO <- data.frame(
                seq_folder   = SEQ_FOLDER,
                seq_id       = SeqIdList,
                sample_name  = SeqIdList,
                sample_group = names(SeqIdList)
            )
        }
    }
#====================================================================================================================================================#

#===[ EXPRESSION VALUE XLSX OUTPUT ]=================================================================================================================#
    if( CONFIG$GLOBAL$CREATE_EXPR_VALUES_XLSX )
    {
        message("|--->> create expression values matrix and export as file.")
        # count-matrix ---------------------------------------------------------
        count_matrix <- Create.Expression.Values.Matrix( SeqIdList=SeqIdList, ValueType="count", DataSource="rds", CountAsInteger=TRUE, Convert2Matrix=TRUE )
        colnames(count_matrix) <- SINFO[match(colnames(count_matrix), SINFO$seq_id), "sample_name"]
        # TPM matrix -----------------------------------------------------------
        tpm_matrix   <- Create.Expression.Values.Matrix( SeqIdList=SeqIdList, ValueType="tpm",   DataSource="rds", Convert2Matrix=TRUE )
        colnames(tpm_matrix) <- SINFO[match(colnames(tpm_matrix), SINFO$seq_id), "sample_name"]
        # log2 TPM matrix ------------------------------------------------------
        log_tpm_matrix <- round(log2(tpm_matrix + 1), 4)
        # CPM matrix -----------------------------------------------------------
        cpm_matrix <- edgeR::cpm( count_matrix, normalized.lib.sizes=FALSE, log=FALSE )
        # log2 CPM matrix -----------------------------------------------------
        log_cpm_matrix <- log2( cpm_matrix + 1 )
        # genes and gene information -------------------------------------------
        ens_geneid <- sapply( rownames(count_matrix), function(y) unlist(strsplit(y, "\\."))[1] )
        geneinfo   <- GENE_INFO[match(ens_geneid, GENE_INFO$ens_geneid), c("symbol","entrez","name","locus_type")] %>% dplyr::rename( gene=symbol, gene_type=locus_type )
        # create export list ---------------------------------------------------
        READ_COUNT = data.frame( ensembl = ens_geneid, count_matrix,   geneinfo )
        TPM        = data.frame( ensembl = ens_geneid, tpm_matrix,     geneinfo )
        log2TPM    = data.frame( ensembl = ens_geneid, log_tpm_matrix, geneinfo )
        CPM        = data.frame( ensembl = ens_geneid, cpm_matrix,     geneinfo )
        log2CPM    = data.frame( ensembl = ens_geneid, log_cpm_matrix, geneinfo )
        rownames(READ_COUNT) = rownames(TPM) = rownames(log2TPM) = rownames(CPM) = rownames(log2CPM) = NULL
        # export as XLSX file --------------------------------------------------
        openxlsx::write.xlsx( READ_COUNT, sprintf("%s/%s.expression.values.Read.Count.xlsx", BaseAnalysisDir, SEQ_FOLDER ), rowNames=FALSE, overwrite=TRUE )
        openxlsx::write.xlsx( list( TPM = TPM, log2TPM = log2TPM ),
            sprintf("%s/%s.expression.values.TPM.xlsx", BaseAnalysisDir, SEQ_FOLDER ), rowNames=FALSE, overwrite=TRUE 
        )
        openxlsx::write.xlsx( list( CPM = CPM, log2CPM = log2CPM ),
            sprintf("%s/%s.expression.values.CPM.xlsx", BaseAnalysisDir, SEQ_FOLDER ), rowNames=FALSE, overwrite=TRUE 
        )
        openxlsx::write.xlsx(
            list( READ_COUNT = READ_COUNT, TPM = TPM, log2TPM = log2TPM, CPM = CPM, log2CPM = log2CPM ),
            sprintf("%s/%s.expression.values.all.merged.xlsx", BaseAnalysisDir, SEQ_FOLDER ), rowNames=FALSE, overwrite=TRUE
        )
    }
#====================================================================================================================================================#

#===[ SAMPLE CLUSTERING ANALYSIS ]===================================================================================================================#
    # sample clustering analysis result folder --------------------------------#
    clustering_res_dir <- sprintf("%s/sample_clustering", BaseAnalysisDir)
    if( !dir.exists(clustering_res_dir) ){ system(sprintf("mkdir -p %s", clustering_res_dir)) }
    message(sprintf("|--->> sample clustering analysis result folder = %s", clustering_res_dir))
    #--------------------------------------------------------------------------#
#---/ UMAP CLUSTERING /----------------------------------------------------------------------------#    
    if( CONFIG$CLUSTERING$UMAP$RUN_UMAP_ANALYSIS )
    {
        # params ---------------------------------------------------------------
        #if( is.null(CONFIG$CLUSTERING$UMAP$UMAP_COMPONENTS) ){ UMAP_COMPONENTS   <- 3     }else{ UMAP_COMPONENTS  <- CONFIG$CLUSTERING$UMAP$UMAP_COMPONENTS  }
        if( is.null(CONFIG$CLUSTERING$UMAP$UMAP_EXPR_VALUE) ){ UMAP_EXPR_VALUE   <- "tpm" }else{ UMAP_EXPR_VALUE  <- CONFIG$CLUSTERING$UMAP$UMAP_EXPR_VALUE  }
        if( is.null(CONFIG$CLUSTERING$UMAP$UMAP_GENE_FILTER)){ UMAP_GENE_FILTER  <- 0.2   }else{ UMAP_GENE_FILTER <- CONFIG$CLUSTERING$UMAP$UMAP_GENE_FILTER }      
        # check gene filter ----------------------------------------------------
        if( !is.null(UMAP_GENE_FILTER) ){ 
            if( UMAP_GENE_FILTER == "deg" ){ 
                if( is.null(UMAP_DEG_LIST)) 
                { stop("|---!!! DEG is selected as umap gene filtering method, however DEG list NOT FOUND. REQUIRED. please check again. STOPPED.") } 
            }else{
                UMAP_GENE_FILTER <- as.numeric(UMAP_GENE_FILTER) 
            }           
        }
        # run UMAP analysis =======================================================================#
            UmapAnalysisResult_3D <- run.UMAP.Clustering(
                ExprValueMatrix     = Create.Expression.Values.Matrix( SeqIdList=SeqIdList, ValueType="tpm", DataSource="rds", Convert2Matrix=TRUE ),
                SeqIdList           = SeqIdList,
                SeqIdGroup          = names(SeqIdList), 
                LowExpressionFilter = TRUE, 
                ValueType           = UMAP_EXPR_VALUE, 
                GeneFilter          = UMAP_GENE_FILTER,  
                ComponentN          = 3
            )
            UmapAnalysisResult_2D <- run.UMAP.Clustering(
                ExprValueMatrix     = Create.Expression.Values.Matrix( SeqIdList=SeqIdList, ValueType="tpm", DataSource="rds", Convert2Matrix=TRUE ),
                SeqIdList           = SeqIdList,
                SeqIdGroup          = names(SeqIdList), 
                LowExpressionFilter = TRUE, 
                ValueType           = UMAP_EXPR_VALUE, 
                GeneFilter          = UMAP_GENE_FILTER,  
                ComponentN          = 2
            )
        #==========================================================================================#

        # add seq_folder id ---------------------------------------------------#
        UmapAnalysisResult_3D$seq_folder <- SEQ_FOLDER
        UmapAnalysisResult_2D$seq_folder <- SEQ_FOLDER
        
        # save UMAP result as file --------------------------------------------#
        if( CONFIG$CLUSTERING$UMAP$UMAP_RESULT_SAVE_AS_FILE )
        {
            write.table( 
                UmapAnalysisResult_3D,
                sprintf("%s/%s_umap.clustering.3D.%s.gene.filter.%s.result.table.tsv", clustering_res_dir, SEQ_FOLDER, UMAP_EXPR_VALUE, UMAP_GENE_FILTER ),
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
            )
            write.table( 
                UmapAnalysisResult_2D,
                sprintf("%s/%s_umap.clustering.2D.%s.gene.filter.%s.result.table.tsv", clustering_res_dir, SEQ_FOLDER, UMAP_EXPR_VALUE, UMAP_GENE_FILTER ),
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
            )
            message("|--->> umap analysis result table saved as file.")
        }
        # import UMAP result into database -------------------------------------
        if( CONFIG$CLUSTERING$UMAP$UMAP_RESULT_DB_IMPORT )
        {
            dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData_3D <- dbGetQuery(dbCon, sprintf("DELETE FROM sample_clustering_result WHERE seq_folder = '%s' AND analysis_method = 'UMAP' AND analysis_components = '3'", SEQ_FOLDER) )
            UpdateData_3D      <- dbWriteTable(dbCon, name="sample_clustering_result", value=UmapAnalysisResult_3D, row.names=FALSE, append=TRUE)
            DeleteExtstData_2D <- dbGetQuery(dbCon, sprintf("DELETE FROM sample_clustering_result WHERE seq_folder = '%s' AND analysis_method = 'UMAP' AND analysis_components = '2'", SEQ_FOLDER) )
            UpdateData_2D      <- dbWriteTable(dbCon, name="sample_clustering_result", value=UmapAnalysisResult_2D, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
            message("|--->> umap analysis result table imported into database.")
        }
    }
#--------------------------------------------------------------------------------------------------#
#---/ PCA CLUSTERING /-----------------------------------------------------------------------------#
    if( CONFIG$CLUSTERING$PCA$RUN_PCA_ANALYSIS )
    {
        # run PCA analysis ------------------------------------------------------------------------#
        PcaAnalysisResult <- run.PCA.Clustering( 
            ExprValueMatrix     = Create.Expression.Values.Matrix( SeqIdList=SeqIdList, ValueType="tpm", DataSource="rds", Convert2Matrix=TRUE ),
            SeqIdList           = SeqIdList, 
            SeqIdGroup          = names(SeqIdList), 
            LowExpressionFilter = TRUE, 
            Log2                = CONFIG$CLUSTERING$PCA$LOG2_PCA_INPUT, 
            PcaComponenetsN     = CONFIG$CLUSTERING$PCA$PCA_COMPONENTS, 
            GroupColors         = NULL
        )
        PcaAnalysisResult$sample_stats$seq_folder   <- SEQ_FOLDER
        PcaAnalysisResult$feature_stats$seq_folder  <- SEQ_FOLDER
        PcaAnalysisResult$feature_stats$var_contrib <- sapply( PcaAnalysisResult$feature_stats$var_contrib, function(z) unlist(strsplit(z, "\\."))[1] )
        #------------------------------------------------------------------------------------------#

        # add gene info to feature stats --------------------------------------#
        PcaAnalysisResult$feature_stats$gene_symbol <- GENE_INFO[match(PcaAnalysisResult$feature_stats$var_contrib, GENE_INFO$ens_geneid), "symbol" ]
        PcaAnalysisResult$feature_stats$gene_type   <- GENE_INFO[match(PcaAnalysisResult$feature_stats$var_contrib, GENE_INFO$ens_geneid), "locus_type" ]

        # PCA result save as file ---------------------------------------------#
        if( CONFIG$CLUSTERING$PCA$PCA_RESULT_SAVE_AS_FILE )
        {
            openxlsx::write.xlsx(
                PcaAnalysisResult[1:2],
                sprintf("%s/%s.PCA.analysis.result.xlsx", clustering_res_dir, SEQ_FOLDER),
                rowNames=FALSE, roverwrite=TRUE
            )
            message("|--->> PCA analysis result saved as XLSX file.")
        }
        # PCA result object save as RDS ---------------------------------------#
        if( CONFIG$CLUSTERING$PCA$PCA_RESULT_SAVE_AS_RDS )
        { 
            saveRDS( PcaAnalysisResult$analysis_object, file=sprintf("%s/%s.PCA.analysis.result.object.rds", clustering_res_dir, SEQ_FOLDER) )
            message("|--->> PCA analysis result object saved as RDS.")
        }
        # PCA result import into database -------------------------------------#
        if( CONFIG$CLUSTERING$PCA$PCA_RESULT_DB_IMPORT )
        {
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            # sample stats -----------------------------------------------------
            DeleteExtstData1 <- dbGetQuery(dbCon, sprintf("DELETE FROM sample_clustering_result WHERE seq_folder = '%s' AND analysis_method = 'PCA'", SEQ_FOLDER) )
            UpdateData1      <- dbWriteTable(dbCon, name="sample_clustering_result", value=PcaAnalysisResult$sample_stats, row.names=FALSE, append=TRUE)
            # feature stats ----------------------------------------------------
            DeleteExtstData2 <- dbGetQuery(dbCon, sprintf("DELETE FROM pca_analysis_feature_stats WHERE seq_folder = '%s'", SEQ_FOLDER))
            UpdateData2      <- dbWriteTable(dbCon, name="pca_analysis_feature_stats", value=PcaAnalysisResult$feature_stats, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
            message("|--->> PCA analysis result table imported into database.")
        }
    }
#--------------------------------------------------------------------------------------------------#
#---/ SAMPLE DISTANCE /----------------------------------------------------------------------------#
    if( CONFIG$CLUSTERING$DISTANCE$DISTANCE_CALCULATION )  
    {
        SampleDistanceResult <- Expression.Profile.Distance(
            ExprMatrix     = Create.Expression.Values.Matrix( 
                                SeqIdList      = SeqIdList, 
                                ValueType      = CONFIG$CLUSTERING$DISTANCE$EXPR_VALUE_TYPE, 
                                DataSource     = "rds", 
                                Convert2Matrix = TRUE ),
            SeqFolderID    = SEQ_FOLDER,
            DistanceMetric = CONFIG$CLUSTERING$DISTANCE$DISTANCE_METRIC,
            DbImport       = CONFIG$CLUSTERING$DISTANCE$RESULT_DB_IMPORT
        )
        if( CONFIG$CLUSTERING$DISTANCE$SAVE_AS_FILE )
        {
            write.table(
                SampleDistanceResult,
                sprintf("%s/%s.%s.metric.%s.sample.distance.calcuation.result.tsv", 
                    clustering_res_dir, SEQ_FOLDER, CONFIG$CLUSTERING$DISTANCE$DISTANCE_METRIC, CONFIG$CLUSTERING$DISTANCE$EXPR_VALUE_TYPE 
                ),
                quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
            )
        }
        if( CONFIG$CLUSTERING$DISTANCE$SAVE_AS_RDS )
        {
            saveRDS( SampleDistanceResult, file=sprintf("%s/%s.%s.metric.%s.sample.distance.calcuation.result.rds", 
                clustering_res_dir, SEQ_FOLDER, CONFIG$CLUSTERING$DISTANCE$DISTANCE_METRIC, CONFIG$CLUSTERING$DISTANCE$EXPR_VALUE_TYPE ) 
            )
        }
    }
#--------------------------------------------------------------------------------------------------#
#---/ SAMPLE CLUSTERING RESULT PLOTS /-------------------------------------------------------------#
    # UMAP result plot --------------------------------------------------------#
    if( CONFIG$CLUSTERING_RESULT_PLOTS$DRAW_UMAP_PLOT )
    {
        # data load from database ----------------------------------------------
        dbCon     <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
        umap_data <- dbGetQuery(dbCon, sprintf("SELECT * FROM sample_clustering_result WHERE analysis_method = 'UMAP' AND seq_folder = '%s'", SEQ_FOLDER))
        dbDisconnect(dbCon)
        # data check and draw plots --------------------------------------------
        if( nrow(umap_data) < 1 ){ message("|---!!! no umap analysis result found on database. please check again.") 
        }else{
            umap_plot_name <- sprintf("%s/%s.UMAP.3D.clustering.analysis.%s.gene.filter.%s.plot.png",
                clustering_res_dir, SEQ_FOLDER, CONFIG$CLUSTERING$UMAP$UMAP_EXPR_VALUE, CONFIG$CLUSTERING$UMAP$UMAP_GENE_FILTER            
            )
            ## draw plots 3D ===================================================
            umap_result_plot <- WTS_plot.UMAP.Clustering.Scatter2(
                UmapResult        = umap_data %>% filter( analysis_components == 3 ),
                PlotName          = umap_plot_name,
                PlotTitle         = "UMAP-3D",
                PlotDimension     = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_DIMENSION,
                Width             = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_WIDTH,
                Height            = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_HEIGHT,
                PointSize         = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_POINT_SIZE,
                PointPch          = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_POINT_PCH,
                AxisFontSize      = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_FONT_SIZE,
                AxisTitleFontSize = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_TITLE_FONT_SIZE,
                PlotSave          = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_SAVE,
                DisplayOnScreen   = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_SCREEN_DISPLAY
            )
            ##
            umap_plot_name2 <- sprintf("%s/%s.UMAP.2D.clustering.analysis.%s.gene.filter.%s.plot.png",
                clustering_res_dir, SEQ_FOLDER, CONFIG$CLUSTERING$UMAP$UMAP_EXPR_VALUE, CONFIG$CLUSTERING$UMAP$UMAP_GENE_FILTER            
            )
            ## draw plots 2D ===================================================
            umap_result_plot <- WTS_plot.UMAP.Clustering.Scatter2(
                UmapResult        = umap_data %>% filter( analysis_components == 2 ),
                PlotName          = umap_plot_name2,
                PlotTitle         = "UMAP-2D",
                PlotDimension     = 2,
                Width             = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_WIDTH,
                Height            = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_HEIGHT,
                PointSize         = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_POINT_SIZE,
                PointPch          = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_POINT_PCH,
                AxisFontSize      = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_FONT_SIZE,
                AxisTitleFontSize = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_TITLE_FONT_SIZE,
                PlotSave          = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_SAVE,
                DisplayOnScreen   = CONFIG$CLUSTERING_RESULT_PLOTS$UMAP_PLOT_SCREEN_DISPLAY
            )
            ## symbolic link for report preparation ============================
            system( sprintf("ln -Tsf %s %s/%s.UMAP.3D.plot.png", umap_plot_name, clustering_res_dir, SEQ_FOLDER ) )
            system( sprintf("ln -Tsf %s %s/%s.UMAP.2D.plot.png", umap_plot_name2, clustering_res_dir, SEQ_FOLDER ) )
            #-------------------------------------------------------------------
            message("|--->> umap analysis result plot created.")
        }
    }
    # PCA result plot ---------------------------------------------------------#
    if( CONFIG$CLUSTERING_RESULT_PLOTS$DRAW_PCA_PLOT )
    {
        # data load from database ----------------------------------------------
        dbCon             <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
        sample_stat_data  <- dbGetQuery(dbCon, sprintf("SELECT * FROM sample_clustering_result WHERE analysis_method = 'PCA' AND seq_folder = '%s'", SEQ_FOLDER))
        feature_stat_data <- dbGetQuery(dbCon, sprintf("SELECT * FROM pca_analysis_feature_stats WHERE seq_folder = '%s'", SEQ_FOLDER))
        dbDisconnect(dbCon)
        # data check -----------------------------------------------------------
        if( nrow(sample_stat_data) < 1 ){ stop("|---!!! no PCA analysis result sample stats data found on database. please check again.") 
        }else{
           if( CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_TYPE == "biplot" )
           {
                if( nrow(feature_stat_data) < 1 )
                {
                    stop("|---!!! no PCA analysis result feature stats data found on database. please check again. required for biplot type.") 
                }else{
                    pca_plot_name <- sprintf("%s/%s.PCA.clustering.analysis.biplot.png", clustering_res_dir, SEQ_FOLDER )
                }
            }else{
                pca_plot_name     <- sprintf("%s/%s.PCA.clustering.analysis.sampleOnly.scatter.plot.png", clustering_res_dir, SEQ_FOLDER )
                #feature_stat_data <- NULL
            }
        }
        sample_stat_data$sample_name <- SINFO[match(sample_stat_data$seq_id, SINFO$seq_id), "sample_name"]
        # draw plot ------------------------------------------------------------
        pca_result_plot <- WTS_plot.PCA.Clustering.Scatter(
            pca_sample_stats  = sample_stat_data,
            pca_feature_stats = feature_stat_data,
            PlotType          = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_TYPE,
            PlotName          = pca_plot_name,
            PlotTitle         = "PCA",
            Width             = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_WIDTH,
            Height            = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_HEIGHT,
            PointSize         = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_POINT_SIZE,
            PointPch          = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_POINT_PCH,
            AxisFontSize      = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_FONT_SIZE,
            AxisTitleFontSize = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_TITLE_FONT_SIZE,
            LabelsRepel       = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_LABEL_REPEL,
            PlotSave          = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_SAVE,
            DisplayOnScreen   = CONFIG$CLUSTERING_RESULT_PLOTS$PCA_PLOT_SCREEN_DISPLAY   
        )
        system(sprintf("ln -Tsf %s %s/%s.PCA.plot.png", pca_plot_name, clustering_res_dir, SEQ_FOLDER))
        message("|--->> PCA analysis result plot created.")
    }
    # Sample Distance Heatmap -------------------------------------------------#
    if( CONFIG$CLUSTERING_RESULT_PLOTS$DRAW_DISTANCE_HEATMAP )
    {
        # data load from database ----------------------------------------------
        dbCon     <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
        dist_data <- dbGetQuery(dbCon, sprintf("SELECT * FROM sample_dist_calculation WHERE metric = '%s' AND seq_folder = '%s'", CONFIG$CLUSTERING_RESULT_PLOTS$HEATMAP_DISTANCE_METRIC, SEQ_FOLDER))
        dbDisconnect(dbCon)

        #-----------------------------------------------------------------------
        dist_matrix <- Convert.Distance.Table.To.Matrix( DistanceResultTable = dist_data, SampleOrderList = SeqIdList )
        #-----------------------------------------------------------------------
        if( CONFIG$CLUSTERING_RESULT_PLOTS$DRAW_DISTANCE_HEATMAP )
        {
            dist_matrix_heatmap <- WTS_plot.Sample.Distance.Heatmap(
                DistanceMatrix = dist_matrix,
                FontSize       = 11,
                Clustering     = TRUE,
                CellSize       = 6
            )
            # heatmap save 
            if( CONFIG$CLUSTERING_RESULT_PLOTS$DIST_HEATMAP_SAVE )
            {
                dist_heatmap_name <- sprintf("%s/%s.metric.%s.sample.distance.heatmap.png", clustering_res_dir, SEQ_FOLDER, CONFIG$CLUSTERING_RESULT_PLOTS$HEATMAP_DISTANCE_METRIC)
                dist_heatmap_tmp  <- sprintf("%s/%s.metric.%s.sample.distance.heatmap_tmp.png", clustering_res_dir, SEQ_FOLDER, CONFIG$CLUSTERING_RESULT_PLOTS$HEATMAP_DISTANCE_METRIC)
                #---------------------------------------------------------------
                png( dist_heatmap_tmp, width=ncol(dist_matrix), height=ncol(dist_matrix), units="in", res=150 )
                draw( dist_matrix_heatmap, padding = unit(c(5, 5, 5, 5), "cm") )
                dev.off()
                # trimming png and remvoe temp file 
                system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", dist_heatmap_tmp, dist_heatmap_name))
                system(sprintf("rm %s", dist_heatmap_tmp))
                system(sprintf("ln -Tsf %s %s/%s.Sample.Distance.Hetamap.png", dist_heatmap_name, clustering_res_dir, SEQ_FOLDER))
            }    
            if( CONFIG$CLUSTERING_RESULT_PLOTS$DIST_HEATMAP_SCREEN_DESPLAY )
            {
                x11( width=ncol(dist_matrix), height=ncol(dist_matrix), bg="white")
                draw( dist_matrix_heatmap )
            }
        }
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===[ ANALYSIS INFO PREPARE ]========================================================================================================================#
    # load config file ---------------------------------------------------------
    ainfo_config_file_type <- tail(unlist(strsplit( CONFIG$ANALYSIS_INFO$ANALYSIS_SET_PAIR_CONFIG_FILE, "\\.")), 1)
    if( ainfo_config_file_type %in% c("yaml","yml") )
    {
        message(sprintf("|---> load DEG analysis set pair config file. [%s]", CONFIG$ANALYSIS_INFO$ANALYSIS_SET_PAIR_CONFIG_FILE))
        ControlTreatPairList <- read_yaml(CONFIG$ANALYSIS_INFO$ANALYSIS_SET_PAIR_CONFIG_FILE)
        ControlTreatPairList <- lapply(ControlTreatPairList, function(ctp) lapply(ctp, function(z) if(!is.null(z)){ unlist(strsplit(z, ","))} ))
    }else if( ainfo_config_file_type == "rds" )
    {
        message(sprintf("|---> load DEG analysis set pair config file. [%s]", CONFIG$ANALYSIS_INFO$ANALYSIS_SET_PAIR_CONFIG_FILE))
        ControlTreatPairList <- readRDS(CONFIG$ANALYSIS_INFO$ANALYSIS_SET_PAIR_CONFIG_FILE)
    }else if( ainfo_config_file_type %in% c("tsv","txt") )
    {
        message(sprintf("|---> load DEG analysis set pair config file. [%s]", CONFIG$ANALYSIS_INFO$ANALYSIS_SET_PAIR_CONFIG_FILE))
        CTpairs <- read.delim(CONFIG$ANALYSIS_INFO$ANALYSIS_SET_PAIR_CONFIG_FILE)
        ControlTreatPairList <- lapply( sapply( CTpairs$aid, function(id) list(id) ), function(z)
        {   
            list(
                ctrl  = CTpairs %>% filter( aid == z ) %>% .$ctrl  %>% unique(),
                treat = CTpairs %>% filter( aid == z ) %>% .$treat %>% unique()
            )
        })
    }else{
        stop("|---!!! no DEG analysis set pair config file. REQUIRED. please check again. config file should be the one of yaml, rds, and tab-deliminated txt or tsv.")
    }
    # create DEG analysis info table list --------------------------------------
    ainfo <- Prepapre.WTS.Analysis.Info.Table(
        SeqFolderId          = SEQ_FOLDER, 
        SampleInfoTable      = SINFO,
        ControlTreatPairList = ControlTreatPairList, 
        GetInfoFromDb        = CONFIG$ANALYSIS_INFO$ANLAYSIS_INFO_GET_DATA_FROM_DB, 
        CreateInfo           = CONFIG$ANALYSIS_INFO$ANALYSIS_INFO_CREATION, 
        DbImport             = CONFIG$ANALYSIS_INFO$ANALYSIS_INFO_DB_IMPROT,
        GenomeAssembly       = 'hg19'
    )
    # analysis id index ---------------------------------------------------------
    ANALYSIS_ID_INDEX <- sapply(unique(ainfo$analysis_id), function(z) list(z))
#====================================================================================================================================================#

#===[ DEG ANALYSIS ]=================================================================================================================================#
    # deg analysis result folder ----------------------------------------------#
    deg_res_dir <- sprintf("%s/deg_analysis", BaseAnalysisDir)
    if( !dir.exists(deg_res_dir) ){ system(sprintf("mkdir -p %s", deg_res_dir)) }
    message(sprintf("DEG analysis result folder = %s ", deg_res_dir))
    # global deg analysis rds folder ------------------------------------------#
    DEG_RDS_DIR = sprintf( "%s/RDS_DegAnalysis", CONFIG$GLOBAL$RDS_DIR )
    #--------------------------------------------------------------------------#
    
#---/ RUN DEG ANALYSIS /---------------------------------------------------------------------------#
    if( CONFIG$DEG_ANALYSIS$RUN_DEG_ANALYSIS )
    {
        # DEG analysis set index -----------------------------------------------
        deg_analysis_set_index <- Create.DEG.Analysis.Set.Table.List( SeqFolderId = SEQ_FOLDER, LoadFromDb = CONFIG$DEG_ANALYSIS$LOAD_ANALYSIS_INFO_FROM_DB )
        
        # run DESeq2 DEG analysis -----------------------------------------------------------------#
        deg_analysis_result_list <- mclapply( deg_analysis_set_index, function(aid) 
        {
            deseq2_analysis_result <- run.DESeq2.DEG.Analysis(
                ExprCountMatrix = Create.Expression.Values.Matrix( SeqIdList=aid$sample_id, ValueType="count", DataSource="rds", CountAsInteger=TRUE, Convert2Matrix=TRUE ),
                AnalysisDesign           = aid,
                FilterByCount            = TRUE,
                FilterByTPM              = TRUE,
                DegAssignMethod          = CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_CUTOFF_METHOD,
                PvalueCutOff             = CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_PVAL_CUTOFF,
                FDRCutOff                = CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_FDR_CUTOFF,
                logFoldChangeCutOff      = log2(CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_FC_CUTOFF),
                DegRatio                 = CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_RATIO_CUTOFF,
                MatrixCreationDataSource = "rds",
                WriteParamsFile          = TRUE,
                ParamsFileSaveDir        = deg_res_dir,
                ImportParamsIntoDb       = TRUE,
                GeneInfo                 = GENE_INFO
            )
            return(deseq2_analysis_result)
        }, mc.cores = length(deg_analysis_set_index) )
        #------------------------------------------------------------------------------------------#
      
        # export as file : individual TSV and RDS ------------------------------
        DegResExport <- list()
        DegResTable  <- data.frame()
        # integration ans summarization of all AID's results -------------------
        for( k in 1:length(deg_analysis_result_list) )
        {
            deg_res.k <- deg_analysis_result_list[[k]] %>% mutate( seq_folder = SEQ_FOLDER, analysis_id=names(deg_analysis_result_list)[k] ) %>% filter( keep )
            rownames(deg_res.k) <- NULL
            deg_res.k2 <- deg_res.k %>% mutate( 
                id_score = paste( ens_geneid, score, sep=","),
                id_fc    = paste( ens_geneid, log2_foldchange, sep=","),
                id_pval  = paste( ens_geneid, pvalue, sep=",")
            )
            deg_res.k2[which(deg_res.k2$deg ==""), "deg"] <- "noDEG"
            deg_res.k2 <- deg_res.k2 %>% mutate( id_degtag = paste( ens_geneid, deg, sep=",") )
            # save as tsv ------------------------------------------------------
            if( CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_RESULT_SAVE_AID_TEXT )
            {
                write.table(
                    deg_res.k %>% filter( keep ), 
                    sprintf("%s/%s.%s.DEG.analysis.result.table.tsv", deg_res_dir, SEQ_FOLDER, names(deg_analysis_result_list)[k] ),
                    quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
                )
                message(sprintf("|---> individual AID DEG result saved as text file. AID = %s", names(deg_analysis_result_list)[k]))
            }
            # integrate into deg-analysis-summary-list -------------------------
            DegResExport[[ names(deg_analysis_result_list)[k] ]] <- deg_res.k %>% filter( keep )
            # save RDS ---------------------------------------------------------
            if( CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_RESULT_SAVE_AID_RDS )
            {
                saveRDS( deg_res.k %>% filter( keep ), sprintf("%s/%s.%s.DEG.analysis.result.table.rds", deg_res_dir, SEQ_FOLDER, names(deg_analysis_result_list)[k] ) )
                saveRDS( deg_res.k %>% filter( keep ), sprintf("%s/%s.%s.DEG.analysis.result.table.rds", DEG_RDS_DIR, SEQ_FOLDER, names(deg_analysis_result_list)[k] ) )
                message(sprintf("|---> individual AID DEG result saved as RDS. AID = %s", names(deg_analysis_result_list)[k]))
            }
            # tables to import into database -----------------------------------
            DegResTable <- rbind(DegResTable,
                data.frame(
                    seq_folder  = SEQ_FOLDER,
                    analysis_id = names(deg_analysis_result_list)[k],
                    rds_path    = sprintf("%s/%s.%s.DEG.analysis.result.table.rds", DEG_RDS_DIR, SEQ_FOLDER, names(deg_analysis_result_list)[k] ),
                    deg_tag     = paste(deg_res.k2$id_degtag, collapse=";"),
                    deg_score   = paste(deg_res.k2$id_score,  collapse=";"),
                    deg_logfc   = paste(deg_res.k2$id_fc,     collapse=";"),
                    deg_pvalue  = paste(deg_res.k2$id_pval,   collapse=";")
                )
            )
        }
        # write summary-list as XLSX -------------------------------------------
        if( CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_RESULT_WRITE_XLSX )
        {
            aid_name_list = lapply(indexing(ainfo$analysis_id), function(idx){
                data.frame(aid_name = paste0(
                    ainfo %>% filter( analysis_id == idx, analysis_group == "TREAT" ) %>% .$sample_group %>% unique(),
                    "_vs_",
                    ainfo %>% filter( analysis_id == idx, analysis_group == "CTRL" ) %>% .$sample_group %>% unique()
                ))
            }) %>% ldply(.id="analysis_id") %>% mutate(analysis_id=as.character(analysis_id))
            DegResExport <- lapply(DegResExport, function(DRE) {
                DRE <- DRE %>% left_join(., aid_name_list, by="analysis_id") %>% dplyr::select(!c("seq_folder"))
                return(DRE)
            })
            openxlsx::write.xlsx( 
                DegResExport,
                sprintf("%s/%s.DEG.analysis.result.table.xlsx", deg_res_dir, SEQ_FOLDER ),
                rowNames=FALSE, overwite=TRUE
            )
            message(sprintf("|---> all AIDs DEG results saved as XLSX. seq_folder = %s", SEQ_FOLDER) )
        }
        # save seq-folder level RDS --------------------------------------------
        if( CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_RESULT_SAVE_SEQ_FOLDER_RDS )
        {
            saveRDS( DegResExport, file=sprintf("%s/%s.DEG.analysis.result.table.list.rds", DEG_RDS_DIR, SEQ_FOLDER) )
            message(sprintf("|---> all AIDs DEG results saved as RDS. seq_folder = %s", SEQ_FOLDER) )
        }
        # import seq-folder level result table into database -------------------
        if( CONFIG$DEG_ANALYSIS$DEG_ANALYSIS_RESULT_DB_IMPORT )
        {
            dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM deg_analysis_result WHERE seq_folder = '%s'", SEQ_FOLDER))
            UpdateData      <- dbWriteTable(dbCon, name="deg_analysis_result", value=DegResTable, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
            message(sprintf("|--->> DEG analysis result table imported into database. seq_folder = %s", SEQ_FOLDER) )
        }
    }
#--------------------------------------------------------------------------------------------------#
#---/ DEG ANALYSIS PLOTS /-------------------------------------------------------------------------#
    # load DEG analysis data --------------------------------------------------#
    if( CONFIG$DEG_PLOT$DEG_DATA_SOURCE == "rds" )
    {
        deg_data_list <- readRDS(sprintf("%s/%s.DEG.analysis.result.table.list.rds", DEG_RDS_DIR, SEQ_FOLDER))
    }else{
        dbCon    <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
        deg_data <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_result WHERE seq_folder = '%s'", SEQ_FOLDER))
        dbDisconnect(dbCon)
        # convert data foramt --------------------------------------------------
        deg_data_list <- mclapply( ANALYSIS_ID_INDEX, function(aid) 
            load.DEG.Analysis.Profile( SeqFolder=SEQ_FOLDER, AnalysisId=aid, GeneInfo=GENE_INFO ) %>% arrange(dplyr::desc(abs(score))), 
            mc.cores = length(ANALYSIS_ID_INDEX)
        )
    }
    #--------------------------------------------------------------------------#

    # deg set volcano plot ----------------------------------------------------#
    if( CONFIG$DEG_PLOT$VOLCANO_PLOT$DRAW_VOLCANO_PLOT )
    {
        # volcano plot list ----------------------------------------------------
        deg_analysis_volcano_plot_list <- lapply( ANALYSIS_ID_INDEX, function(aid) 
        {
            # sample group name ------------------------------------------------
            ctrl_group_name  <- ainfo %>% filter( analysis_id == aid, analysis_group == "CTRL" )   %>% .$sample_group %>% unique()
            treat_group_name <- ainfo %>% filter( analysis_id == aid, analysis_group == "TREAT" )  %>% .$sample_group %>% unique()
            # draw plots -------------------------------------------------------
            volcanoPlot <- WTS_plot.DEG.Analysis.Result.Volcano( 
                DegAnalysisResult = deg_data_list[[ aid ]], 
                ColorUp           = "#c73500", 
                ColorDown         = "#008000", 
                ScoreValue        = CONFIG$DEG_PLOT$VOLCANO_PLOT$SCORE_METRIC,
                PointSize         = CONFIG$DEG_PLOT$VOLCANO_PLOT$POINT_SIZE, 
                AxisFontSize      = CONFIG$DEG_PLOT$VOLCANO_PLOT$FONT_SZIE, 
                AxisTitleFontSize = CONFIG$DEG_PLOT$VOLCANO_PLOT$TITLE_FONT_SIZE,
                Width             = CONFIG$DEG_PLOT$VOLCANO_PLOT$PLOT_WIDTH, 
                Height            = CONFIG$DEG_PLOT$VOLCANO_PLOT$PLOT_HEIGHT, 
                PointPch          = CONFIG$DEG_PLOT$VOLCANO_PLOT$POINT_PCH,
                PlotSave          = CONFIG$DEG_PLOT$VOLCANO_PLOT$PLOT_SAVE, 
                DisplayOnScreen   = CONFIG$DEG_PLOT$VOLCANO_PLOT$PLOT_SCREEN_DISPLAY, 
                PlotName          = sprintf("%s/%s.%s.DEG.analysis.volcano.plot.png", deg_res_dir, SEQ_FOLDER, aid), 
                PlotTitle         = sprintf("%s : %s", treat_group_name, ctrl_group_name)
            )
            return(volcanoPlot)
        })
    }
    # deg only heatmap --------------------------------------------------------#
    if( CONFIG$DEG_PLOT$DEG_HEATMAP$DRAW_DEG_HEATMAP )
    {
        # sample list check ----------------------------------------------------
        deg_heatmap_sample_list <- SeqIdList
        deg_data_filter <- lapply( ANALYSIS_ID_INDEX, function(aid) 
        {
            ainfo_aid <- ainfo %>% filter( analysis_id == aid ) 
            if( all(ainfo_aid$sample_id %in% deg_heatmap_sample_list) ){ include_aid <- TRUE }else{ include_aid <- FALSE }
            return(include_aid) 
        })
        heatmap_deg_profiles <- names(which(deg_data_filter == TRUE))
        # heatmap deg data -----------------------------------------------------
        heatmap_deg_data_list <- deg_data_list[ heatmap_deg_profiles ]
        all_degs <- Reduce(union, lapply( heatmap_deg_data_list, function(HD) HD[which(HD$deg !=""), "ens_geneid"] ))
        # prepare input matrix 
        HeatmapInputMatrix <- Create.Expression.Values.Matrix( SeqIdList=deg_heatmap_sample_list, ValueType="tpm", DataSource="db", Convert2Matrix=TRUE )
        rownames(HeatmapInputMatrix) <- sapply( rownames(HeatmapInputMatrix), function(RN) unlist(strsplit(RN, "\\."))[1] )
        HeatmapInputMatrix <- HeatmapInputMatrix[all_degs, ]
        # show border option ---------------------------------------------------
        if( CONFIG$DEG_PLOT$DEG_HEATMAP$NO_BORDERS ){ HeatmapBorderColor <- NA }else{ HeatmapBorderColor <- NULL }
        # column label rotation angle and centering ----------------------------
        if( max(nchar(deg_heatmap_sample_list)) <= 5 ){ 
            ColumnLabelRotationAngle <- 0 
            ColmunNameCentering      <- TRUE
        }else{ 
            ColumnLabelRotationAngle <- 45 
            ColmunNameCentering      <- FALSE
        }

        # draw heatmap ----------------------------------------------------------------------------#
        deg_heatmap <- WTS_plot.Gene.Expression.Profiles.Heatmap(
            ExpressionValueMatrix = HeatmapInputMatrix,
            ColorUp               = "#b20000",
            ColorDown             = "#005900",
            Log2                  = CONFIG$DEG_PLOT$DEG_HEATMAP$INPUT_LOG2,
            Scales                = CONFIG$DEG_PLOT$DEG_HEATMAP$SCALES,
            RowClustering         = CONFIG$DEG_PLOT$DEG_HEATMAP$CLUSTERING_ROWS,
            ColumnClustering      = CONFIG$DEG_PLOT$DEG_HEATMAP$CLUSTERING_COLUMNS,
            ColorBorder           = ifelse( CONFIG$DEG_PLOT$DEG_HEATMAP$NO_BORDERS, NA, NULL ),
            RowLabelsSide         = CONFIG$DEG_PLOT$DEG_HEATMAP$ROW_LABEL_SIDE,
            ColumnLabelsSide      = CONFIG$DEG_PLOT$DEG_HEATMAP$COLUMN_LABEL_SIDE,          
            ColumnLabelsRotation  = ColumnLabelRotationAngle,
            RowLabelsSize         = 2,
            RowLabelsColors       = "#1e1e1e",
            ColumnNamesCentered   = ColmunNameCentering
        )
        #------------------------------------------------------------------------------------------#
        
        # heatmap save as png file --------------------------------------------#
        if( CONFIG$DEG_PLOT$DEG_HEATMAP$HEATMAP_SAVE_AS_PNG )
        {
            # heatmap save as png
            deg_heatmap_file_tmp <- sprintf( "%s/%s.DEG.only.heatmap_tmp.png", deg_res_dir, SEQ_FOLDER )
            deg_heatmap_filename <- sprintf( "%s/%s.DEG.only.scales.%s.heatmap.png", deg_res_dir, SEQ_FOLDER, CONFIG$DEG_PLOT$DEG_HEATMAP$SCALES )
            png(deg_heatmap_file_tmp, width=30, height=10, units="in", res=150 )
            draw( deg_heatmap, padding = unit(c(1, 1, 5, 1), "cm") )
            dev.off()
            # trimming png and remvoe temp file 
            system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", deg_heatmap_file_tmp, deg_heatmap_filename))
            system(sprintf("rm %s", deg_heatmap_file_tmp))
            message(sprintf("|--->> DEGs only heatmap (scales = %s) created.", CONFIG$DEG_PLOT$DEG_HEATMAP$SCALES))
        }
        # heatmap display on screen 
        if( CONFIG$DEG_PLOT$DEG_HEATMAP$HEATMAP_SCREEN_DISPLAY )
        {
            x11( width=15, height=10, bg="white")
            draw( deg_heatmap, padding = unit(c(5, 5, 15, 5), "mm") )
        }
    }
    # CPM correlation and scatter pair plot -----------------------------------#
    if( CONFIG$DEG_PLOT$DEG_SET_PAIR_PLOT$DRAW_DEG_SET_PAIR_PLOT )
    {
        # prepare CPM matrix ---------------------------------------------------
        CountMatrix <- Create.Expression.Values.Matrix( SeqIdList=SINFO$seq_id, ValueType="count", DataSource="rds", CountAsInteger=TRUE, Convert2Matrix=TRUE )
        CPM_MATRIX  <- edgeR::cpm(CountMatrix, log=TRUE, normalized.lib.sizes=TRUE, prior.count=1)
        # deg analysis set -----------------------------------------------------
        if( !"deg_analysis_set_index" %in% ls() )   
        { deg_analysis_set_index <- Create.DEG.Analysis.Set.Table.List( SeqFolderId = SEQ_FOLDER, LoadFromDb = TRUE) }
        
        # draw pair plots -------------------------------------------------------------------------#
        draw_aid_set_pair_plots <- lapply( deg_analysis_set_index, function(aidset)
        {
            aidset$sample_name           <- SINFO[match(aidset$sample_id, SINFO$seq_id), "sample_name"]
            pair_plot_sample_list        <- aidset$sample_id
            names(pair_plot_sample_list) <- aidset$sample_name
            pair_plot_title <- paste0( 
                unique(aidset[which(aidset$analysis_group == "TREAT"), "sample_group"]), 
                " vs. ", 
                unique(aidset[which(aidset$analysis_group == "CTRL"),  "sample_group"]) 
            )
            pair_plot_filename <- sprintf("%s/%s.%s.CPM.correlation.scatter.pair.plot.png", deg_res_dir, SEQ_FOLDER, unique(aidset$analysis_id))

            message(sprintf("|--->> %s.%s. creating pair plot...", SEQ_FOLDER, unique(aidset$analysis_id)))
            WTS_plot.CPM.Correlation.Pairs(
                CpmMatrix       = CPM_MATRIX,
                SampleList      = pair_plot_sample_list,
                PlotTitle       = pair_plot_title,
                PlotSave        = CONFIG$DEG_PLOT$DEG_SET_PAIR_PLOT$PLOT_SAVE,
                DisplayOnScreen = CONFIG$DEG_PLOT$DEG_SET_PAIR_PLOT$DISPLAY_ON_SCREEN,
                PlotName        = pair_plot_filename
            )
        })
    }
    #--------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===[ PATHWAY ANALYSIS AND PLOTS ]===================================================================================================================#
    # pathway analsyis result folder ------------------------------------------#
    pathway_res_dir <- sprintf("%s/pathway_analysis", BaseAnalysisDir)
    if( !dir.exists(pathway_res_dir) ){ system(sprintf("mkdir -p %s", pathway_res_dir)) }
    message(sprintf("pathway analysis result folder = %s ", pathway_res_dir))
    # global pathway analysis global rds folder -------------------------------#
    PATH_ANALYSIS_RDS_DIR <- sprintf("%s/RDS_PathwayAnalysis", CONFIG$GLOBAL$RDS_DIR)
    #--------------------------------------------------------------------------#
    
#---/ RUN PATHWAY ANALYSIS /-----------------------------------------------------------------------#
    message(rep("-", 150))
    message("    ### Pathway Analysis Running ###")
    message(" ")
    if( CONFIG$ORA$RUN_PATHWAY_ANALYSIS )
    {
        # load DEG analysis result ---------------------------------------------
        message("    [ Loading DEG Analysis Result Data ]")
        message(sprintf("|---> DEG analysis result data source : %s", CONFIG$ORA$DEG_DATA_SOURCE))
        deg_analysis_result_list <- lapply( ANALYSIS_ID_INDEX, function(aidx) 
            load.DEG.Analysis.Profile( SeqFolder = SEQ_FOLDER, AnalysisId = aidx, DataSource = CONFIG$ORA$DEG_DATA_SOURCE, GeneInfo = GENE_INFO )
        )
        # pathway ids for analysis ---------------------------------------------
        message("    [ Analysis Target PathwaySets and Subsets ]")
        message(sprintf("|---> analysis target PathwaySetID : %s", CONFIG$ORA$SELECTED_PATHWAY_IDS))
        PathIds              <- unlist(strsplit(CONFIG$ORA$SELECTED_PATHWAY_IDS, ","))
        AnalysisPathwayIndex <- sapply(unique(PathIds), function(z) list(z)) 
        # pathway subset ids for analysis --------------------------------------
        if( !is.null(CONFIG$ORA$SELECTED_PATHWAY_SUBSETS) )
        { 
            message(sprintf("|---> selected pathway subsets : %s", CONFIG$ORA$SELECTED_PATHWAY_SUBSETS))
            PathSubsets <- unlist(strsplit(CONFIG$ORA$SELECTED_PATHWAY_SUBSETS, ",")) 
        }else{ 
            message(sprintf("|---> no pathway subsets selected. "))
            PathSubsets <- NULL 
        }       
        # parallel cores -------------------------------------------------------
        if( length(deg_analysis_result_list) > 5 ){ ora_parallel_cores <- 5 }else{ ora_parallel_cores <- length(deg_analysis_result_list) }
        message(sprintf("|---> parallel processing cores set to %s ( total AnalysisIDs : %s )", ora_parallel_cores, length(deg_analysis_result_list) ))

        # run pathway analysis --------------------------------------------------------------------#
        message(" ")
        message("    [ Run Pathway Analysis ]")
        PathwayAnalysisResultList <- mclapply( deg_analysis_result_list, function(degaid) 
        {
            ora_res_list <- lapply( AnalysisPathwayIndex, function(pid)
            {
                message(rep(".", 50))
                message(" ")
                if( length(PathSubsets) > 0 )
                {
                    message(sprintf("|>>>> Pathway Analysis - AnalysisID : %s , PathwaySetID : %s , Subset : %s", unique(degaid$deg_analysis_id), pid, ps))
                    ora_res <- data.frame()
                    for( ps in PathSubsets ) 
                    {      
                        ora_res_subset <- run.DEG.Pathway.Analysis(
                            DegAnalysisResult = degaid, 
                            PathwaySetID      = pid, 
                            PathSubset        = ps, 
                            AddMergeDegRun    = CONFIG$ORA$ADD_TOTAL_DEG_RUN, 
                            GeneIdentifier    = CONFIG$ORA$GENE_IDENTIFIER, 
                            Threads           = CONFIG$ORA$THREADS
                        )
                        ora_res <- rbind(ora_res, ora_res_subset)
                    }
                }else{
                    message(sprintf("|>>>> Pathway Analysis - AnalysisID : %s , PathwaySetID : %s ", unique(degaid$deg_analysis_id), pid))
                    ora_res <- run.DEG.Pathway.Analysis(
                        DegAnalysisResult = degaid, 
                        PathwaySetID      = pid, 
                        PathwaySetName    = "cpdb",
                        PathSubset        = NULL, 
                        AddMergeDegRun    = CONFIG$ORA$ADD_TOTAL_DEG_RUN, 
                        GeneIdentifier    = CONFIG$ORA$GENE_IDENTIFIER, 
                        Threads           = CONFIG$ORA$THREADS
                    )
                }
                return(ora_res)
            })
            aid_ora_analysis_result           <- do.call( rbind, ora_res_list )
            rownames(aid_ora_analysis_result) <- NULL
            return(aid_ora_analysis_result)
        }, mc.cores = ora_parallel_cores )
        #------------------------------------------------------------------------------------------#
        message(" ")
        # top200 results -------------------------------------------------------
        message("|>>>> Extract TopRanked 200 Pathways ")
        PathResTop200List <- lapply(PathwayAnalysisResultList, function(PR) 
        {
            analysis_unit <- unique(PR[,c("input_geneset","analysis_pathset_id","analysis_subset_name")]) %>% filter( !is.na(input_geneset) )
            if( nrow(analysis_unit) > 0 ){    
                TOP200_RESULT <- data.frame()        
                for( i in 1:nrow(analysis_unit) )
                {
                    top200_res <- PR %>% filter( 
                        input_geneset        == analysis_unit[i, "input_geneset"       ],
                        analysis_pathset_id  == analysis_unit[i, "analysis_pathset_id" ],
                        analysis_subset_name == analysis_unit[i, "analysis_subset_name"]
                    ) %>% arrange(dplyr::desc(score))
                    if( nrow(top200_res) > 200 ){ top200_res <- top200_res %>% dplyr::slice(1:200) }
                    TOP200_RESULT <- rbind(TOP200_RESULT, top200_res)
                    top200_res    <- data.frame() 
                }
            }else{
                TOP200_RESULT <- PR
            }
            message(sprintf("|---> %s top200 pathway extraction is completed.", unique(PR$deg_analysis_id)))
            return(TOP200_RESULT)
        })
        # remove no-gene-enriched pathways
        # PathwayAnalysisResultList <- lapply(PathwayAnalysisResultList, function(res) res %>% filter( enrich_genes >= 3 ) )
        # PathResTop200List         <- lapply(PathResTop200List, function(res) res %>% filter( enrich_genes >= 3 ) )
        message(" ")
        message("    [ Save Pathway Analysis Results ]")
        # >> whole results save : RDS, TEXT | top200 results save : RDS, database
        # whole profiles save as RDS -------------------------------------------
        if( CONFIG$ORA$WHOLE_RESULT_SAVE_RDS )
        {
            message(">>>> Analysis Results Integration on AnalysisID-Level and Save as RDS ")
            message(" ")
            # analysis-id level whole result RDS save --------------------------
            for( k in 1:length(PathwayAnalysisResultList) )
            {
                saveRDS( PathwayAnalysisResultList[[k]], 
                    file=sprintf("%s/%s.%s.pathway.analysis.whole.result.table.rds", pathway_res_dir, SEQ_FOLDER, names(PathwayAnalysisResultList)[k] )
                )
                message(sprintf("|---> AnalysisID %s result RDS saved in pathway analysis result folder of AnalysisBatch.", unique(PathwayAnalysisResultList[[k]]$deg_analysis_id)))

                saveRDS( PathwayAnalysisResultList[[k]], 
                    file=sprintf("%s/%s.%s.pathway.analysis.whole.result.table.rds", PATH_ANALYSIS_RDS_DIR, SEQ_FOLDER, names(PathwayAnalysisResultList)[k] )
                )
                message(sprintf("|---> AnalysisID %s result RDS saved in global pathway analysis result RDS folder.", unique(PathwayAnalysisResultList[[k]]$deg_analysis_id)))
                message(rep(".", 50))
            }
            # seq_folder level result list RDS save ----------------------------
            message(" ")
            message(">>>> Analysis Results List Object Save as RDS ")
            saveRDS( PathwayAnalysisResultList, file=sprintf("%s/%s.pathway.analysis.whole.result.table.rds", pathway_res_dir, SEQ_FOLDER ) )
            message("|---> All results list RDS saved in pathway analysis result folder of AnalysisBatch.")
            saveRDS( PathwayAnalysisResultList, file=sprintf("%s/%s.pathway.analysis.whole.result.table.rds", PATH_ANALYSIS_RDS_DIR, SEQ_FOLDER ) )
            message("|---> All results list RDS saved in global pathway analysis result RDS folder.")
        }
        # seq_folder level whole result list XLSX save -------------------------
        if( CONFIG$ORA$WHOLE_RESULT_WRITE_XLSX )
        {
            message(">>>> All Analysis Results List Export as XLSX File ")
            AidSampleGroupName <- reshape(unique(ainfo[,c("analysis_id","analysis_group","sample_group")]), idvar='analysis_id', timevar='analysis_group', direction='wide') %>%
                mutate( labels = paste(sample_group.TREAT, sample_group.CTRL, sep="_"))
            aid_labels        <- AidSampleGroupName$labels
            names(aid_labels) <- AidSampleGroupName$analysis_id

            aid_name_list = lapply(indexing(ainfo$analysis_id), function(idx){
                data.frame(aid_name = paste0(
                    ainfo %>% filter( analysis_id == idx, analysis_group == "TREAT" ) %>% .$sample_group %>% unique(),
                    "_vs_",
                    ainfo %>% filter( analysis_id == idx, analysis_group == "CTRL" ) %>% .$sample_group %>% unique()
                ))
            }) %>% ldply(.id="analysis_id") %>% mutate(analysis_id=as.character(analysis_id))

            PathwayAnalysisResultList2 <- lapply(PathwayAnalysisResultList, function(PRE) {
                PRE_rev <- data.frame(
                    deg_analysis_name = aid_name_list[match(PRE$deg_analysis_id, aid_name_list$analysis_id), "aid_name"],
                    PRE %>% dplyr::select(!c("seq_folder"))
                )
                return(PRE_rev)
            })

            names(PathwayAnalysisResultList2) <- aid_labels[ names(PathwayAnalysisResultList2) ]
            openxlsx::write.xlsx( PathwayAnalysisResultList2,
                sprintf("%s/%s.pathway.analysis.whole.result.table.xlsx", pathway_res_dir, SEQ_FOLDER ),
                rowNames=FALSE, overwrite=TRUE
            )
            message(sprintf("|---> All analysis results list is saved as XLSX file in pathway analysis result folder of AnalysisBatch."))
            message(" ")
        }
        # top200 result save ---------------------------------------------------
        message(">>>> TopRank-200 Result Save as RDS and Import into Database")
        for( k in 1:length(PathResTop200List) )
        {
            message(rep(".", 50))
            message(" ")
            # save RDS ---------------------------------------------------------
            if( CONFIG$ORA$TOP200_RESULT_SAVE_RDS )
            {
                saveRDS( PathResTop200List[[k]], 
                    file=sprintf("%s/%s.%s.pathway.analysis.top200.result.table.rds", pathway_res_dir, SEQ_FOLDER, names(PathResTop200List)[k] )
                )
                message(sprintf("|---> AnalysisID %s TopRank-200 result RDS saved in pathway analysis result folder of AnalysisBatch.",names(PathResTop200List)[k]))
                saveRDS( PathResTop200List[[k]], 
                    file=sprintf("%s/%s.%s.pathway.analysis.top200.result.table.rds", PATH_ANALYSIS_RDS_DIR, SEQ_FOLDER, names(PathResTop200List)[k] )
                )
                message(sprintf("|---> AnalysisID %s TopRank-200 results list RDS saved in global pathway analysis result RDS folder.",names(PathResTop200List)[k]))
            }
            # import into database ---------------------------------------------
            if(  CONFIG$ORA$TOP200_RESULT_DB_IMPORT )
            {
                message(">>>> TopRank-200 Result Import into Databse ")
                dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                DeleteExtstData <- dbGetQuery(dbCon, 
                    sprintf("DELETE FROM pathway_analysis_top200_result WHERE seq_folder = '%s' AND deg_analysis_id = '%s'", 
                        unique(PathResTop200List[[k]]$seq_folder), unique(PathResTop200List[[k]]$deg_analysis_id) 
                    )
                )
                UpdateData      <- dbWriteTable(dbCon, name="pathway_analysis_top200_result", value=PathResTop200List[[k]], row.names=FALSE, append=TRUE)
                dbDisconnect(dbCon)
                message(sprintf("|---> AnalysisID %s TopRank-200 results imported into database.", names(PathResTop200List)[k]))
            }
        }
        # create pathway analysis params ---------------------------------------
        message(" ")
        message(">>>> Pathway Analysis Params Table Creation and Save")
        pathway_analysis_run_params <- do.call(rbind, 
            lapply( PathwayAnalysisResultList, function(y) unique(y[,c("seq_folder","deg_analysis_id","input_geneset","analysis_pathset_id","analysis_subset_name")]) )
        ) %>% group_by(seq_folder, deg_analysis_id) %>% reframe(
            input_deg = paste(unique(input_geneset), collapse=","), 
            run_pathset_id = paste(unique(analysis_pathset_id), collapse=","), 
            pathway_subsets = paste(unique(analysis_subset_name), collapse=",")
        ) %>% as.data.frame()
        # pathway analysis params save as text ---------------------------------
        if( CONFIG$ORA$ANALYSIS_PARAMS_SAVE_AS_TEXT )
        {
            write.table(
                pathway_analysis_run_params,
                sprintf("%s/%s.pathway.analysis.params.tsv", pathway_res_dir, SEQ_FOLDER),
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
            )
            message("|---> Analysis params table saved as a text file in pathway analysis result folder of AnalysisBatch.")
        }
        # pathway analysis params import into database -------------------------
        if( CONFIG$ORA$ANALYSIS_PARAMS_DB_IMPORT )
        {
            dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM pathway_analysis_params WHERE seq_folder = '%s'", SEQ_FOLDER))
            UpdateData      <- dbWriteTable(dbCon, name="pathway_analysis_params", value=pathway_analysis_run_params, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
            message("|---> Analysis params table imported into database." ) 
        }
        message(" ")
    }
    message(" ")
#--------------------------------------------------------------------------------------------------#
#---/ PATHWAY ANALYSIS RESULT PLOTS /--------------------------------------------------------------#
    message(rep("-", 150))
    message("    ### Pathway Analysis Result Plots ###")
    message(" ")
    # load pathway analysis data ----------------------------------------------#
    message("    [ Loading Pathway Analysis Result Data ]")
    if( CONFIG$ORA_PLOT$ORA_DATA_SOURCE == "rds" )
    {
        message(sprintf("|---> Load data from RDS, AnalysisBatch: %s", SEQ_FOLDER))
        path_res_data_list <- readRDS(sprintf("%s/%s.pathway.analysis.whole.result.table.rds", PATH_ANALYSIS_RDS_DIR, SEQ_FOLDER))
    }else{
        message(sprintf("|---> Load data from DATABASE, AnalysisBatch: %s", SEQ_FOLDER))
        path_res_data_list <- lapply( ANALYSIS_ID_INDEX, function(aid) load.Pathway.Analysis.Profile(SeqFolder=SEQ_FOLDER, AnalysisId=aid, DataSource='db', RdsDir=PATH_ANALYSIS_RDS_DIR ) )
    }
    # top-rank pathway result dot-plot & bar-plot -----------------------------#
    if( CONFIG$ORA_PLOT$DRAW_ORA_PLOTS )
    {
        message("    [ Create Plots for Top-Ranked Pathways ]")
        # TopRank Profiles Plots ---------------------------------------------------    
        drawPathAnalysisTopRankResPlots <- lapply( ANALYSIS_ID_INDEX, function(AnalysisID) 
        {    
            # pathwayids and subsets check -------------------------------------
            OraResult         <- path_res_data_list[[ AnalysisID ]]

            if( all(!is.na(OraResult$score)) ){
                plot_pathway_list <- unique(OraResult[,c("analysis_pathset_id","analysis_subset_name")])
                plot_pathway_list$pathway_name <- sapply( plot_pathway_list$analysis_pathset_id, function(pwid){
                    dbCon    <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db="public_database" )
                    pathname <- dbGetQuery(dbCon, sprintf("SELECT set_name FROM pathset_info WHERE set_id = '%s'", pwid))
                    dbDisconnect(dbCon)
                    return(pathname$set_name)
                })
                rownames(plot_pathway_list) <- NULL
                message(" ")
                message(sprintf("|>>>> PathwaySets and Subsets List : AnalysisBatch %s ", AnalysisID))
                message(" ")
                print(plot_pathway_list)
                message(" ")
                # dotplot
                if( CONFIG$ORA_PLOT$DOTPLOT$DRAW_TOP_RANK_DOTPLOT )
                {
                    message(rep(".", 50))
                    message(" ")
                    message(sprintf("|---> Create TopRank %s Pathways Dotplots for each AnalysisBatch.", CONFIG$ORA_PLOT$PLOT_TOP_RANK_N))
                    for( k in 1:nrow(plot_pathway_list) )
                    {
                        PathSetID.k <- plot_pathway_list[k, "analysis_pathset_id" ]
                        Subset.k    <- plot_pathway_list[k, "analysis_subset_name"]
                        PathName.K  <- plot_pathway_list[k, "pathway_name"        ]
                        if( Subset.k == "" )
                        {
                            OraRes.k       <- OraResult %>% filter( analysis_pathset_id == PathSetID.k )
                            PlotFilePrefix <- sprintf("%s/%s.%s.pathway.%s.top.rank.%s", pathway_res_dir, SEQ_FOLDER, AnalysisID, PathSetID.k, CONFIG$ORA_PLOT$PLOT_TOP_RANK_N) 
                            message(sprintf("|---> TopRank Pathway Dotplot : AnalysisBatch = %s, AnalysisID = %s, PathwaySet = %s", SEQ_FOLDER, AnalysisID, PathName.K))
                        }else{
                            OraRes.k       <- OraResult %>% filter( analysis_pathset_id == PathSetID.k, analysis_subset_name == Subset.k )
                            PlotFilePrefix <- sprintf("%s/%s.%s.pathway.%s.subset.%s.top.rank.%s", pathway_res_dir, SEQ_FOLDER, AnalysisID, PathSetID.k, Subset.k, CONFIG$ORA_PLOT$PLOT_TOP_RANK_N) 
                            message(sprintf("|---> TopRank Pathway Dotplot : AnalysisBatch = %s, AnalysisID = %s, PathwaySet = %s, Subset = %s", SEQ_FOLDER, AnalysisID, PathName.K, Subset.k))
                        }
                        toprankpathres_dotplot_up <- WTS_plot.Pathway.Analysis.TopRank.Result(
                            PathwayAnalysisResult = OraRes.k,
                            GeneType              = "UP",
                            PlotType              = "dot",
                            TopReults_N           = CONFIG$ORA_PLOT$PLOT_TOP_RANK_N,
                            AxisFontSize          = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_FONT_SIZE,
                            Width                 = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_WIDTH,
                            Height                = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_HEIGHT,
                            PlotSave              = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_SAVE,
                            DisplayOnScreen       = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_SCREEN_DISPLAY,
                            PlotName              = sprintf("%s.UP.genes.dotplot.png", PlotFilePrefix),
                            PlotTitle             = NULL
                        )
                        toprankpathres_dotplot_down <- WTS_plot.Pathway.Analysis.TopRank.Result(
                            PathwayAnalysisResult = OraRes.k,
                            GeneType              = "DOWN",
                            PlotType              = "dot",
                            TopReults_N           = CONFIG$ORA_PLOT$PLOT_TOP_RANK_N,
                            AxisFontSize          = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_FONT_SIZE,
                            Width                 = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_WIDTH,
                            Height                = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_HEIGHT,
                            PlotSave              = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_SAVE,
                            DisplayOnScreen       = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_SCREEN_DISPLAY,
                            PlotName              = sprintf("%s.DOWN.genes.dotplot.png", PlotFilePrefix),
                            PlotTitle             = NULL
                        )
                        if( CONFIG$ORA_PLOT$DOTPLOT$MERGED_DOTPLOT_SAVE )
                        {
                            if( CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_MERGE_DIRECTION == 1 )
                            {
                                toprankpathres_dotplot <- plot_grid(toprankpathres_dotplot_up, toprankpathres_dotplot_down, nrow=1, rel_widths=c(1,1))
                                suppressWarnings(ggsave( toprankpathres_dotplot, 
                                    file   = sprintf("%s.merged.dotplot.png", PlotFilePrefix),
                                    width  = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_WIDTH * 2, 
                                    height = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_HEIGHT, 
                                    unit = 'in', dpi = 150, type = 'cairo'   
                                ))
                                message("|---> Merged Dotplot also created. Merged direction : horizontal")
                            }else{
                                toprankpathres_dotplot <- plot_grid(toprankpathres_dotplot_up, toprankpathres_dotplot_down, ncol=1, rel_heights=c(1,1))
                                suppressWarnings(ggsave( toprankpathres_dotplot, 
                                    file   = sprintf("%s.merged.dotplot.png", PlotFilePrefix),
                                    width  = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_WIDTH, 
                                    height = CONFIG$ORA_PLOT$DOTPLOT$DOTPLOT_HEIGHT * 2, 
                                    unit ='in', dpi=150, type='cairo'   
                                ))
                                message("|---> Merged Dotplot also created. Merged direction : vertical")
                            }
                        }
                    }
                }
                # barplot               
                if( CONFIG$ORA_PLOT$BARPLOT$DRAW_TOP_RANK_BARPLOT )
                {
                    message(rep(".", 50))
                    message(" ")
                    message(sprintf("|---> Create TopRank %s Pathways Barplots for each AnalysisBatch.", CONFIG$ORA_PLOT$PLOT_TOP_RANK_N))
                    for( k in 1:nrow(plot_pathway_list) )
                    {
                        PathSetID.k <- plot_pathway_list[k, "analysis_pathset_id" ]
                        Subset.k    <- plot_pathway_list[k, "analysis_subset_name"]
                        PathName.K  <- plot_pathway_list[k, "pathway_name"        ]
                        if( Subset.k == "" )
                        {
                            OraRes.k       <- OraResult %>% filter( analysis_pathset_id == PathSetID.k )
                            PlotFilePrefix <- sprintf("%s/%s.%s.pathway.%s.top.rank.%s", pathway_res_dir, SEQ_FOLDER, AnalysisID, PathName.K, CONFIG$ORA_PLOT$PLOT_TOP_RANK_N) 
                            message(sprintf("|---> TopRank Pathway Barplot : AnalysisBatch = %s, AnalysisID = %s, PathwaySet = %s", SEQ_FOLDER, AnalysisID, PathName.K))
                        }else{
                            OraRes.k       <- OraResult %>% filter( analysis_pathset_id == PathSetID.k, analysis_subset_name == Subset.k )
                            PlotFilePrefix <- sprintf("%s/%s.%s.pathway.%s.subset.%s.top.rank.%s", pathway_res_dir, SEQ_FOLDER, AnalysisID, PathName.K, Subset.k, CONFIG$ORA_PLOT$PLOT_TOP_RANK_N) 
                            message(sprintf("|---> TopRank Pathway Barplot : AnalysisBatch = %s, AnalysisID = %s, PathwaySet = %s, Subset = %s", SEQ_FOLDER, AnalysisID, PathName.K, Subset.k))
                        }
                        toprankpathres_barplot_up <- WTS_plot.Pathway.Analysis.TopRank.Result(
                            PathwayAnalysisResult = OraRes.k,
                            GeneType              = "UP",
                            PlotType              = "bar",
                            TopReults_N           = CONFIG$ORA_PLOT$PLOT_TOP_RANK_N,
                            AxisFontSize          = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_FONT_SIZE,
                            Width                 = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_WIDTH,
                            Height                = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_HEIGHT,
                            PlotSave              = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_SAVE,
                            DisplayOnScreen       = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_SCREEN_DISPLAY,
                            PlotName              = sprintf("%s.UP.genes.barplot.png", PlotFilePrefix),
                            PlotTitle             = NULL
                        )
                        toprankpathres_barplot_down <- WTS_plot.Pathway.Analysis.TopRank.Result(
                            PathwayAnalysisResult = OraRes.k,
                            GeneType              = "DOWN",
                            PlotType              = "bar",
                            TopReults_N           = CONFIG$ORA_PLOT$PLOT_TOP_RANK_N,
                            AxisFontSize          = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_FONT_SIZE,
                            Width                 = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_WIDTH,
                            Height                = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_HEIGHT,
                            PlotSave              = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_SAVE,
                            DisplayOnScreen       = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_SCREEN_DISPLAY,
                            PlotName              = sprintf("%s.DOWN.genes.barplot.png", PlotFilePrefix),
                            PlotTitle             = NULL
                        )
                        if( CONFIG$ORA_PLOT$BARPLOT$MERGED_BARPLOT_SAVE )
                        {
                            if( CONFIG$ORA_PLOT$BARPLOT$BARPLOT_MERGE_DIRECTION == 1 )
                            {
                                toprankpathres_barplot <- plot_grid(toprankpathres_barplot_up, toprankpathres_barplot_down, nrow=1, rel_widths=c(1,1))
                                suppressWarnings(ggsave( toprankpathres_barplot, 
                                    file   =  sprintf("%s.merged.barplot.png", PlotFilePrefix),
                                    width  = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_WIDTH * 2, 
                                    height = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_HEIGHT, 
                                    unit='in', dpi=150, type='cairo'   
                                ))
                                message("|---> Merged Barplot also created. Merged direction : horizontal")
                            }else{
                                toprankpathres_barplot <- plot_grid(toprankpathres_barplot_up, toprankpathres_barplot_down, ncol=1, rel_heights=c(1,1))
                                suppressWarnings(ggsave( toprankpathres_barplot, 
                                    file   = sprintf("%s.merged.barplot.png", PlotFilePrefix),
                                    width  = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_WIDTH, 
                                    height = CONFIG$ORA_PLOT$BARPLOT$BARPLOT_HEIGHT * 2, 
                                    unit='in', dpi=150, type='cairo'   
                                ))
                                message("|---> Merged Barplot also created. Merged direction : vertical")
                            }
                        }

                    }
                }
            }
            return('done')
        })
    }
    # top-rank pathway result network 
    if( CONFIG$ORA_PLOT$TOP_RANK_NET$DRAW_ORA_NETWORK )
    {
        PathwayIDs <- unlist(strsplit(CONFIG$ORA$SELECTED_PATHWAY_IDS, ","))
        ora_net_data <- lapply( ANALYSIS_ID_INDEX, function(AnalysisID) 
        {
            for( pid in PathwayIDs )    
            {
                # create network data ----------------------------------------------
                ora_net_data <- create.TopRank.Pathways.and.Genes.NetData(
                    PathwayAnalysisResult     = path_res_data_list[[ AnalysisID ]],
                    MultiPathInvolveThreshold = CONFIG$ORA_PLOT$TOP_RANK_NET$MULTI_GENE_THRESHOLD,
                    PathwaySetId              = pid,
                    TopRank_N                 = CONFIG$ORA_PLOT$TOP_RANK_NET$TOP_RANK_N,
                    UseMergedAnalysisResult   = CONFIG$ORA_PLOT$TOP_RANK_NET$USE_MERGED_DEG_RESULT,
                    DegAnalysisResult         = deg_data_list[[ AnalysisID ]],
                    GeneInfo                  = GENE_INFO
                )
                # network plot -----------------------------------------------------
                ora_net_plotname <- sprintf( "%s/%s.%s.pathway.%s.analysis.top.rank.%s.multi.gene.threshold.%s.network.png", 
                    pathway_res_dir, SEQ_FOLDER, AnalysisID, pid, CONFIG$ORA_PLOT$TOP_RANK_NET$TOP_RANK_N, CONFIG$ORA_PLOT$TOP_RANK_NET$MULTI_GENE_THRESHOLD 
                )
                ora_net_plot_tmp <- sprintf("%s/%s.%s.pathway.%s.analysis.top.rank.network_tmp.png", pathway_res_dir, SEQ_FOLDER, AnalysisID, pid )
                # draw plot ----------------------------------------------------
                if( class(ora_net_data) != "data.frame" ){
                    png( ora_net_plot_tmp, width=CONFIG$ORA_PLOT$TOP_RANK_NET$PLOT_WIDTH, height=CONFIG$ORA_PLOT$TOP_RANK_NET$PLOT_HEIGHT, units="in",res=150 )
                    WTS_plot.Pathway.Analysis.TopRank.Network(
                        TopRankNetData        = ora_net_data,
                        NetworkLayout         = CONFIG$ORA_PLOT$TOP_RANK_NET$NETWORK_LAYOUT,
                        UseVisNetwork         = CONFIG$ORA_PLOT$TOP_RANK_NET$USE_VISNETWORK,
                        VisNetHighlightDegree = CONFIG$ORA_PLOT$TOP_RANK_NET$VISNET_HIGHLIGHT_DEGREE
                    )
                    dev.off()
                    #-----------------------------------------------------------
                    system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", ora_net_plot_tmp, ora_net_plotname))
                    system(sprintf("rm %s", ora_net_plot_tmp))
                    system(sprintf("ln -Tsf %s %s/%s.%s.ORA.%s.network.png", ora_net_plotname, pathway_res_dir, SEQ_FOLDER, AnalysisID, pid ))
                }
            }
        })
    }
    message(" ")
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===[ GSEA ANALYSIS ]================================================================================================================================#
    # GSEA analsyis result folder ---------------------------------------------#
    gsea_res_dir <- sprintf("%s/GSEA", BaseAnalysisDir)
    if( !dir.exists(gsea_res_dir) ){ system(sprintf("mkdir -p %s", gsea_res_dir)) }
    message(sprintf("pathway analysis result folder = %s ", gsea_res_dir))
    # global pathway analysis global rds folder -------------------------------#
    GSEA_RDS_DIR <- sprintf("%s/RDS_GSEA", CONFIG$GLOBAL$RDS_DIR)
    #--------------------------------------------------------------------------#
#---/ GSEA ANALYSIS /------------------------------------------------------------------------------#
    message(rep("-", 150))
    message("    ### GSEA Running ###")
    message(" ")
    if( CONFIG$GSEA$RUN_GSEA )
    {
        # load DEG analysis result ---------------------------------------------
        message("    [ Loading DEG Analysis Result Data ]")
        message(sprintf("|---> DEG analysis result data source : %s", CONFIG$GSEA$DEG_DATA_SOURCE))
        deg_analysis_result_list <- lapply( ANALYSIS_ID_INDEX, function(aidx) 
            load.DEG.Analysis.Profile( SeqFolder = SEQ_FOLDER, AnalysisId = aidx, DataSource = CONFIG$GSEA$DEG_DATA_SOURCE, GeneInfo = GENE_INFO )
        )
        # pathway ids for analysis ---------------------------------------------
        message("    [ Analysis Target PathwaySets and Subsets ]")
        message(sprintf("|---> analysis target PathwaySetID : %s", CONFIG$GSEA$SELECTED_PATHWAY_IDS))
        PathIds              <- unlist(strsplit(CONFIG$GSEA$SELECTED_PATHWAY_IDS, ","))
        AnalysisPathwayIndex <- sapply(unique(PathIds), function(z) list(z)) 
        # pathway subset ids for analysis --------------------------------------
        if( !is.null(CONFIG$GSEA$SELECTED_PATHWAY_SUBSETS) )
        { 
            message(sprintf("|---> selected pathway subsets : %s", CONFIG$GSEA$SELECTED_PATHWAY_SUBSETS))
            GseaSubsets <- unlist(strsplit(CONFIG$GSEA$SELECTED_PATHWAY_SUBSETS, ",")) 
        }else{ 
            message(sprintf("|---> no pathway subsets selected. "))
            GseaSubsets <- NULL 
        }       
        # parallel cores -------------------------------------------------------
        if( length(deg_analysis_result_list) > 5 ){ gsea_parallel_cores <- 5 }else{ gsea_parallel_cores <- length(deg_analysis_result_list) }
        message(sprintf("|---> parallel processing cores set to %s ( total AnalysisIDs : %s )", gsea_parallel_cores, length(deg_analysis_result_list) ))

        # run GSEA --------------------------------------------------------------------------------#
        message(" ")
        message("    [ Run Pathway Analysis ]")
        # pathwaysets of pathway analysis
        GseaResultList <- mclapply( deg_analysis_result_list, function(degaid) 
        {
            gsea_res_list <- lapply( AnalysisPathwayIndex, function(pid)
            {
                message(rep(".", 50))
                message(" ")
                if( length(GseaSubsets) > 0 )
                {
                    message(sprintf("|>>>> GSEA - AnalysisID : %s , PathwaySetID : %s , Subset : %s", unique(degaid$analysis_id), pid, ps))
                    gsea_res <- data.frame()
                    for( ps in GseaSubsets ) 
                    {      
                        gsea_res_subset <- run.GSEA.Analysis(
                            DegAnalysisResult = degaid,
                            AnalysisID        = unique(degaid$analysis_id),
                            PathwaySetID      = pid,
                            PathSubset        = ps, 
                            GeneIdentifier    = CONFIG$GSEA$GENE_IDENTIFIER,
                            RankMetric        = CONFIG$GSEA$RANK_METRIC,
                            Threads           = CONFIG$GSEA$THREADS,
                            MinPathSize       = CONFIG$GSEA$MIN_PATHWAY_SIZE,
                            MaxPathSize       = CONFIG$GSEA$MAX_PATHWAY_SIZE
                        )
                        gsea_res_subset$analysis_pathset_id  <- pid
                        gsea_res_subset$analysis_subset_name <- ps
                        #-------------------------------------------------------
                        gsea_res <- rbind(gsea_res, gsea_res_subset)
                    }
                }else{
                    message(sprintf("|>>>> GSEA - AnalysisID : %s , PathwaySetID : %s ", unique(degaid$analysis_id), pid))
                    gsea_res <- run.GSEA.Analysis(
                        DegAnalysisResult = degaid,
                        AnalysisID        = unique(degaid$analysis_id),
                        PathwaySetID      = pid,
                        PathSubset        = NULL, 
                        GeneIdentifier    = CONFIG$GSEA$GENE_IDENTIFIER,
                        RankMetric        = CONFIG$GSEA$RANK_METRIC,
                        Threads           = CONFIG$GSEA$THREADS,
                        MinPathSize       = CONFIG$GSEA$MIN_PATHWAY_SIZE,
                        MaxPathSize       = CONFIG$GSEA$MAX_PATHWAY_SIZE
                    )
                    gsea_res$analysis_pathset_id  <- pid
                    gsea_res$analysis_subset_name <- ""
                }
                return(gsea_res)
            })
            aid_gsea_result           <- do.call( rbind, gsea_res_list )
            rownames(aid_gsea_result) <- NULL
            return(aid_gsea_result)
        }, mc.cores = gsea_parallel_cores )
        #------------------------------------------------------------------------------------------#

        # top200 gsea results --------------------------------------------------
        message("|>>>> Extract TopRanked 200 GSEA Results ")
        GseaResTop200List <- lapply((GseaResultList), function(PR) 
        {
            TOP200_RESULT <- data.frame()
            analysis_unit <- unique(PR[,c("analysis_pathset_id","analysis_subset_name")])
            
            for( i in 1:nrow(analysis_unit) )
            {
                top100_plus <- PR %>% filter( 
                    analysis_pathset_id  == analysis_unit[i, "analysis_pathset_id" ],
                    analysis_subset_name == analysis_unit[i, "analysis_subset_name"],
                    NES > 0 
                ) %>% arrange(dplyr::desc(score))
                if( nrow(top100_plus) > 100 ){ top100_plus <- top100_plus %>% dplyr::slice(1:100) }
                top100_minus <- PR %>% filter( 
                    analysis_pathset_id  == analysis_unit[i, "analysis_pathset_id" ],
                    analysis_subset_name == analysis_unit[i, "analysis_subset_name"],
                    NES < 0 
                ) %>% arrange(dplyr::desc(abs(score)))
                if( nrow(top100_minus) > 100 ){ top100_minus <- top100_minus %>% dplyr::slice(1:100) }
                top200_res = rbind(top100_plus,top100_minus)
                TOP200_RESULT <- rbind(TOP200_RESULT, top200_res)
                top200_res    <- data.frame() 
            }
            message(sprintf("|---> %s top200 GSEA results extraction is completed.", unique(PR$analysis_id)))
            return(TOP200_RESULT)
        })
        message(" ")
        message("    [ Save GSEA Results ]")
        # >> whole results save : RDS, TEXT | top200 results save : RDS, database
        # whole profiles save as RDS -------------------------------------------
        if( CONFIG$GSEA$WHOLE_RESULT_SAVE_RDS )
        {
            message(">>>> Analysis Results Integration on AnalysisID-Level and Save as RDS ")
            message(" ")
            # analysis-id level whole result RDS save --------------------------
            for( k in 1:length(GseaResultList) )
            {
                saveRDS( GseaResultList[[k]], 
                    file=sprintf("%s/%s.%s.GSEA.whole.result.table.rds", gsea_res_dir, SEQ_FOLDER, names(GseaResultList)[k] )
                )
                message(sprintf("|---> AnalysisID %s result RDS saved in GSEA result folder of AnalysisBatch.", unique(GseaResultList[[k]]$analysis_id)))

                saveRDS( GseaResultList[[k]], 
                    file=sprintf("%s/%s.%s.GSEA.whole.result.table.rds", GSEA_RDS_DIR, SEQ_FOLDER, names(GseaResultList)[k] )
                )
                message(sprintf("|---> AnalysisID %s result RDS saved in global GSEA result RDS folder.", unique(GseaResultList[[k]]$analysis_id)))
                message(rep(".", 50))
            }
            # seq_folder level result list RDS save ----------------------------
            message(" ")
            message(">>>> Analysis Results List Object Save as RDS ")
            saveRDS( GseaResultList, file=sprintf("%s/%s.GSEA.whole.result.table.rds", gsea_res_dir, SEQ_FOLDER ) )
            message("|---> All results list RDS saved in GSEA result folder of AnalysisBatch.")
            saveRDS( GseaResultList, file=sprintf("%s/%s.GSEA.whole.result.table.rds", GSEA_RDS_DIR, SEQ_FOLDER ) )
            message("|---> All results list RDS saved in global GSEA result RDS folder.")
        }
        # seq_folder level whole result list XLSX save -------------------------
        if( CONFIG$GSEA$WHOLE_RESULT_WRITE_XLSX )
        {
            message(">>>> All Analysis Results List Export as XLSX File ")
            AidSampleGroupName <- reshape(unique(ainfo[,c("analysis_id","analysis_group","sample_group")]), idvar='analysis_id', timevar='analysis_group', direction='wide') %>%
                mutate( labels = paste(sample_group.TREAT, sample_group.CTRL, sep="_"))
            aid_labels        <- AidSampleGroupName$labels
            names(aid_labels) <- AidSampleGroupName$analysis_id

            aid_name_list = lapply(indexing(ainfo$analysis_id), function(idx){
                data.frame(aid_name = paste0(
                    ainfo %>% filter( analysis_id == idx, analysis_group == "TREAT" ) %>% .$sample_group %>% unique(),
                    "_vs_",
                    ainfo %>% filter( analysis_id == idx, analysis_group == "CTRL" ) %>% .$sample_group %>% unique()
                ))
            }) %>% ldply(.id="analysis_id") %>% mutate(analysis_id=as.character(analysis_id))

            GseaResultList2 <- lapply(GseaResultList, function(PRE) {
                PRE_rev <- data.frame(
                    deg_analysis_name = aid_name_list[match(PRE$analysis_id, aid_name_list$analysis_id), "aid_name"],
                    PRE %>% dplyr::select(!c("seq_folder"))
                )
                return(PRE_rev)
            })
            names(GseaResultList2) <- aid_labels[ names(GseaResultList2) ]
            openxlsx::write.xlsx( GseaResultList2,
                sprintf("%s/%s.GSEA.whole.result.table.xlsx", gsea_res_dir, SEQ_FOLDER ),
                rowNames=FALSE, overwrite=TRUE
            )
            message(sprintf("|---> All analysis results list is saved as XLSX file in GSEA result folder of AnalysisBatch."))
            message(" ")
        }
        # top200 result save ---------------------------------------------------
        message(">>>> TopRank-200 Result Save as RDS and Import into Database")
        for( k in 1:length(GseaResTop200List) )
        {
            message(rep(".", 50))
            message(" ")
            # save RDS ---------------------------------------------------------
            if( CONFIG$GSEA$TOP200_RESULT_SAVE_RDS )
            {
                saveRDS( GseaResTop200List[[k]], 
                    file=sprintf("%s/%s.%s.GSEA.top200.result.table.rds", gsea_res_dir, SEQ_FOLDER, names(GseaResTop200List)[k] )
                )
                message(sprintf("|---> AnalysisID %s TopRank-200 result RDS saved in GSEA result folder of AnalysisBatch.",names(GseaResTop200List)[k]))
                saveRDS( GseaResTop200List[[k]], 
                    file=sprintf("%s/%s.%s.GSEA.top200.result.table.rds", GSEA_RDS_DIR, SEQ_FOLDER, names(GseaResTop200List)[k] )
                )
                message(sprintf("|---> AnalysisID %s TopRank-200 results list RDS saved in global GSEA result RDS folder.",names(GseaResTop200List)[k]))
            }
            # import into database ---------------------------------------------
            if(  CONFIG$GSEA$TOP200_RESULT_DB_IMPORT )
            {
                message(">>>> TopRank-200 Result Import into Databse ")
                dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                DeleteExtstData <- dbGetQuery(dbCon, 
                    sprintf("DELETE FROM gsea_top200_result WHERE seq_folder = '%s' AND analysis_id = '%s' AND analysis_pathset_id = '%s'", 
                        unique(GseaResTop200List[[k]]$seq_folder), unique(GseaResTop200List[[k]]$analysis_id), unique(GseaResTop200List[[k]]$analysis_pathset_id) 
                    )
                )
                UpdateData      <- dbWriteTable(dbCon, name="gsea_top200_result", value=GseaResTop200List[[k]], row.names=FALSE, append=TRUE)
                dbDisconnect(dbCon)
                message(sprintf("|---> AnalysisID %s TopRank-200 results imported into database.", names(GseaResTop200List)[k]))
            }
        }
        # create GSEA run params -----------------------------------------------
        message(" ")
        message(">>>> GSEA Run Params Table Creation and Save")
        gsea_run_params <- do.call(rbind, 
            lapply( GseaResultList, function(y) unique(y[,c("seq_folder","analysis_id","analysis_pathset_id","analysis_subset_name")]) )
        ) %>% group_by(seq_folder, analysis_id) %>% reframe(
            run_pathset_id = paste(unique(analysis_pathset_id), collapse=","), 
            pathway_subsets = paste(unique(analysis_subset_name), collapse=",")
        ) %>% as.data.frame()
        # GSEA run params save as text -----------------------------------------
        if( CONFIG$GSEA$ANALYSIS_PARAMS_SAVE_AS_TEXT )
        {
            write.table(
                gsea_run_params,
                sprintf("%s/%s.GSEA.run.params.tsv", gsea_res_dir, SEQ_FOLDER),
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
            )
            message("|---> GSEA run params table saved as a text file in GSEA result folder of AnalysisBatch.")
        }
        # GSEA run params import into database ---------------------------------
        if( CONFIG$GSEA$ANALYSIS_PARAMS_DB_IMPORT )
        {
            dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM gsea_run_params WHERE seq_folder = '%s' AND run_pathset_id = '%s'", SEQ_FOLDER, unique(gsea_run_params$run_pathset_id)))
            UpdateData      <- dbWriteTable(dbCon, name="gsea_run_params", value=gsea_run_params, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
            message("|---> GSEA run params table imported into database." ) 
        }
        #----------------------------------------------------------------------#

        # GO-GSEA analysis -----------------------------------------------------
        if( CONFIG$GSEA$INCLUDE_GO_ANALYSIS )
        {
            # GO pathwayset id -------------------------------------------------
            GO_PathwayIDs <- c("gobp"="pw10", "gocc"="pw11", "gomf"="pw12")
            GO_PATH_ID    <- GO_PathwayIDs [ unlist(strsplit(CONFIG$GSEA$SELECTED_GENE_ONTOLOGY, ",")) ]
            #-------------------------------------------------------------------
            # GO set of GO enrich analysis
            GoGseaResultList <- mclapply( deg_analysis_result_list, function(degaid) 
            {
                gsea_go_res <- run.GSEA.Analysis(
                    DegAnalysisResult = degaid,
                    AnalysisID        = unique(degaid$analysis_id),
                    PathwaySetID      = GO_PATH_ID,
                    PathSubset        = NULL, 
                    GeneIdentifier    = CONFIG$GSEA$GENE_IDENTIFIER,
                    RankMetric        = CONFIG$GSEA$RANK_METRIC,
                    Threads           = CONFIG$GSEA$THREADS,
                    MinPathSize       = CONFIG$GSEA$MIN_PATHWAY_SIZE,
                    MaxPathSize       = CONFIG$GSEA$MAX_PATHWAY_SIZE
                )
                gsea_go_res$analysis_pathset_id  <- GO_PATH_ID
                gsea_go_res$analysis_subset_name <- CONFIG$GSEA$SELECTED_GENE_ONTOLOGY
                return(gsea_go_res)
            }, mc.cores = gsea_parallel_cores )

            # top200 gsea results --------------------------------------------------
            message("|>>>> Extract TopRanked 200 GO-GSEA Results ")
            GoGseaResTop200List <- lapply((GoGseaResultList), function(PR) 
            {
                go_top100_plus <- PR %>% filter( NES > 0 ) %>% arrange(dplyr::desc(score))
                if( nrow(go_top100_plus) > 100 ){ go_top100_plus <- go_top100_plus %>% dplyr::slice(1:100) }
                go_top100_minus <- PR %>% filter( NES < 0 ) %>% arrange(dplyr::desc(abs(score)))
                if( nrow(go_top100_minus) > 100 ){ go_top100_minus <- go_top100_minus %>% dplyr::slice(1:100) }
                go_top200_res <- rbind( go_top100_plus, go_top100_minus )
                message(sprintf("|---> %s top200 GO-GSEA results extraction is completed.", unique(PR$analysis_id)))
                return(go_top200_res)
            })
            message(" ")
            message("    [ Save GSEA Results ]")
            # >> whole results save : RDS, TEXT | top200 results save : RDS, database
            # whole profiles save as RDS -------------------------------------------
            if( CONFIG$GSEA$WHOLE_RESULT_SAVE_RDS )
            {
                message(">>>> Analysis Results Integration on AnalysisID-Level and Save as RDS ")
                message(" ")
                # analysis-id level whole result RDS save --------------------------
                for( k in 1:length(GoGseaResultList) )
                {
                    saveRDS( GoGseaResultList[[k]], 
                        file=sprintf("%s/%s.%s.GO.GSEA.whole.result.table.rds", gsea_res_dir, SEQ_FOLDER, names(GoGseaResultList)[k] )
                    )
                    message(sprintf("|---> AnalysisID %s result RDS saved in GSEA result folder of AnalysisBatch.", unique(GoGseaResultList[[k]]$analysis_id)))

                    saveRDS( GoGseaResultList[[k]], 
                        file=sprintf("%s/%s.%s.GO.GSEA.whole.result.table.rds", GSEA_RDS_DIR, SEQ_FOLDER, names(GoGseaResultList)[k] )
                    )
                    message(sprintf("|---> AnalysisID %s result RDS saved in global GSEA result RDS folder.", unique(GoGseaResultList[[k]]$analysis_id)))
                    message(rep(".", 50))
                }
                # seq_folder level result list RDS save ----------------------------
                message(" ")
                message(">>>> Analysis Results List Object Save as RDS ")
                saveRDS( GoGseaResultList, file=sprintf("%s/%s.GO.GSEA.whole.result.table.rds", gsea_res_dir, SEQ_FOLDER ) )
                message("|---> All results list RDS saved in GSEA result folder of AnalysisBatch.")
                saveRDS( GoGseaResultList, file=sprintf("%s/%s.GO.GSEA.whole.result.table.rds", GSEA_RDS_DIR, SEQ_FOLDER ) )
                message("|---> All results list RDS saved in global GSEA result RDS folder.")
            }
            # seq_folder level whole result list XLSX save -------------------------
            if( CONFIG$GSEA$WHOLE_RESULT_WRITE_XLSX )
            {
                message(">>>> All Analysis Results List Export as XLSX File ")
                AidSampleGroupName <- reshape(unique(ainfo[,c("analysis_id","analysis_group","sample_group")]), idvar='analysis_id', timevar='analysis_group', direction='wide') %>%
                    mutate( labels = paste(sample_group.TREAT, sample_group.CTRL, sep="_"))
                aid_labels        <- AidSampleGroupName$labels
                names(aid_labels) <- AidSampleGroupName$analysis_id

                aid_name_list = lapply(indexing(ainfo$analysis_id), function(idx){
                    data.frame(aid_name = paste0(
                        ainfo %>% filter( analysis_id == idx, analysis_group == "TREAT" ) %>% .$sample_group %>% unique(),
                        "_vs_",
                        ainfo %>% filter( analysis_id == idx, analysis_group == "CTRL" ) %>% .$sample_group %>% unique()
                    ))
                }) %>% ldply(.id="analysis_id") %>% mutate(analysis_id=as.character(analysis_id))

                GoGseaResultList2 <- lapply(GoGseaResultList, function(PRE) {
                    PRE_rev <- data.frame(
                        deg_analysis_name = aid_name_list[match(PRE$analysis_id, aid_name_list$analysis_id), "aid_name"],
                        PRE %>% dplyr::select(!c("seq_folder"))
                    )
                    return(PRE_rev)
                })

                names(GoGseaResultList2) <- aid_labels[ names(GoGseaResultList2) ]
                openxlsx::write.xlsx( GoGseaResultList2,
                    sprintf("%s/%s.GO.GSEA.whole.result.table.xlsx", gsea_res_dir, SEQ_FOLDER ),
                    rowNames=FALSE, overwrite=TRUE
                )
                message(sprintf("|---> All analysis results list is saved as XLSX file in pathway analysis result folder of AnalysisBatch."))
                message(" ")
            }
            # top200 result save ---------------------------------------------------
            message(">>>> TopRank-200 Result Save as RDS and Import into Database")
            for( k in 1:length(GoGseaResTop200List) )
            {
                message(rep(".", 50))
                message(" ")
                # save RDS ---------------------------------------------------------
                if( CONFIG$GSEA$TOP200_RESULT_SAVE_RDS )
                {
                    saveRDS( GoGseaResTop200List[[k]], 
                        file=sprintf("%s/%s.%s.GO.GSEA.top200.result.table.rds", gsea_res_dir, SEQ_FOLDER, names(GoGseaResTop200List)[k] )
                    )
                    message(sprintf("|---> AnalysisID %s TopRank-200 result RDS saved in GSEA result folder of AnalysisBatch.",names(GoGseaResTop200List)[k]))
                    saveRDS( GoGseaResTop200List[[k]], 
                        file=sprintf("%s/%s.%s.GO.GSEA.top200.result.table.rds", GSEA_RDS_DIR, SEQ_FOLDER, names(GoGseaResTop200List)[k] )
                    )
                    message(sprintf("|---> AnalysisID %s TopRank-200 results list RDS saved in global GSEA result RDS folder.",names(GoGseaResTop200List)[k]))
                }
                # import into database ---------------------------------------------
                if(  CONFIG$GSEA$TOP200_RESULT_DB_IMPORT )
                {
                    message(">>>> TopRank-200 Result Import into Databse ")
                    dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                    DeleteExtstData <- dbGetQuery(dbCon, 
                        sprintf("DELETE FROM gsea_top200_result WHERE seq_folder = '%s' AND analysis_id = '%s' AND analysis_pathset_id = '%s'", 
                            unique(GoGseaResTop200List[[k]]$seq_folder), unique(GoGseaResTop200List[[k]]$analysis_id), unique(GoGseaResTop200List[[k]]$analysis_pathset_id) 
                        )
                    )
                    UpdateData      <- dbWriteTable(dbCon, name="gsea_top200_result", value=GoGseaResTop200List[[k]], row.names=FALSE, append=TRUE)
                    dbDisconnect(dbCon)
                    message(sprintf("|---> AnalysisID %s TopRank-200 results imported into database.", names(GoGseaResTop200List)[k]))
                }
            }
            # create GSEA run params -------------------------------------------
            message(" ")
            message(">>>> GSEA Run Params Table Creation and Save")
            go_gsea_run_params <- do.call(rbind, 
                lapply( GoGseaResultList, function(y) unique(y[,c("seq_folder","analysis_id","analysis_pathset_id","analysis_subset_name")]) )
            ) %>% group_by(seq_folder, analysis_id) %>% reframe(
                run_pathset_id = paste(unique(analysis_pathset_id), collapse=","), 
                pathway_subsets = paste(unique(analysis_subset_name), collapse=",")
            ) %>% as.data.frame()
            # GSEA run params save as text -------------------------------------
            if( CONFIG$GSEA$ANALYSIS_PARAMS_SAVE_AS_TEXT )
            {
                write.table(
                    go_gsea_run_params,
                    sprintf("%s/%s.GO.GSEA.run.params.tsv", gsea_res_dir, SEQ_FOLDER),
                    quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
                )
                message("|---> GSEA run params table saved as a text file in GSEA run result folder of AnalysisBatch.")
            }
            # GSEA run params import into database -----------------------------
            if( CONFIG$GSEA$ANALYSIS_PARAMS_DB_IMPORT )
            {
                dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM gsea_run_params WHERE seq_folder = '%s' AND run_pathset_id = '%s'", SEQ_FOLDER, unique(go_gsea_run_params$run_pathset_id)))
                UpdateData      <- dbWriteTable(dbCon, name="gsea_run_params", value=go_gsea_run_params, row.names=FALSE, append=TRUE)
                dbDisconnect(dbCon)
                message("|---> GSEA run params table imported into database." ) 
            }
        }
    }
#--------------------------------------------------------------------------------------------------#
#---/ GSEA PLOTS /---------------------------------------------------------------------------------#
    message(rep("-", 150))
    message("    ### Top-Ranked Profiles GSEA Plots ###")
    message(" ")
    message("    [ Loading Pathway Analysis Result Data ]")
    # load GSEA data -----------------------------------------------------------
    if( CONFIG$GSEA_PLOT$GSEA_DATA_SOURCE == "rds" )
    {
        message(sprintf("|---> Load data from RDS, AnalysisBatch: %s", SEQ_FOLDER))
        gsea_data_list <- readRDS(sprintf("%s/%s.GSEA.whole.result.table.rds", PATH_ANALYSIS_RDS_DIR, SEQ_FOLDER))
    }else{
        message(sprintf("|---> Load data from DATABASE, AnalysisBatch: %s", SEQ_FOLDER))
        gsea_data_list <- lapply( ANALYSIS_ID_INDEX, function(aid) load.GSEA.Profile(SeqFolder=SEQ_FOLDER, AnalysisId=aid, DataSource='db', RdsDir=GSEA_RDS_DIR ) )
    }
    #---------------------------------------------------------------------------
    if( CONFIG$GSEA_PLOT$DRAW_GSEA_PLOT )
    {
        if( "deg_data_list" %nin% ls() )
        {
            if( CONFIG$DEG_PLOT$DEG_DATA_SOURCE == "rds" )
            {
                deg_data_list <- readRDS(sprintf("%s/%s.DEG.analysis.result.table.list.rds", DEG_RDS_DIR, SEQ_FOLDER))
            }else{
                dbCon    <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                deg_data <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_result WHERE seq_folder = '%s'", SEQ_FOLDER))
                dbDisconnect(dbCon)
                # convert data foramt --------------------------------------------------
                deg_data_list <- mclapply( ANALYSIS_ID_INDEX, function(aid) 
                    load.DEG.Analysis.Profile( SeqFolder=SEQ_FOLDER, AnalysisId=aid, GeneInfo=GENE_INFO ) %>% arrange(dplyr::desc(abs(score))), 
                    mc.cores = length(ANALYSIS_ID_INDEX)
                )
            }
        }
        #-----------------------------------------------------------------------
        toprank_gsea_plots <- lapply( ANALYSIS_ID_INDEX, function(AnalysisID) 
        {
            # pathawy set and subset list --------------------------------------
            gsea_path_sets <- unique(gsea_data_list[[ AnalysisID ]][,c("analysis_pathset_id","analysis_subset_name")])
            # draw gsea plots --------------------------------------------------
            for( k in 1:nrow(gsea_path_sets) )
            {
                # top rank N profiles ------------------------------------------
                gsea_plot_taraget_plus <- gsea_data_list[[ AnalysisID ]] %>% filter( 
                    analysis_pathset_id == gsea_path_sets[k, "analysis_pathset_id"],
                    analysis_subset_name == gsea_path_sets[k, "analysis_subset_name"]
                ) %>% filter( NES > 0 ) %>% arrange( dplyr::desc(score) ) %>% dplyr::slice(1:CONFIG$GSEA_PLOT$TOP_RANK_N)

                gsea_plot_taraget_minus <- gsea_data_list[[ AnalysisID ]] %>% filter( 
                    analysis_pathset_id == gsea_path_sets[k, "analysis_pathset_id"],
                    analysis_subset_name == gsea_path_sets[k, "analysis_subset_name"]
                ) %>% filter( NES < 0 ) %>% arrange( dplyr::desc(abs(score)) ) %>% dplyr::slice(1:CONFIG$GSEA_PLOT$TOP_RANK_N)
                
                # gene rank profiles -------------------------------------------
                DEG_PROFILE      <- deg_data_list[[ AnalysisID ]] %>% filter( !is.na(entrez) ) 
                GENE_RANK        <- DEG_PROFILE$score
                names(GENE_RANK) <- DEG_PROFILE$entrez
                # draw plots ---------------------------------------------------
                for( i in 1:nrow(gsea_plot_taraget_plus) )
                {
                    gsea_plot_name <- sprintf("%s/%s.%s.GSEA.%s.NES.plus.rank.%s.plot.png", 
                        gsea_res_dir, SEQ_FOLDER, AnalysisID, gsea_plot_taraget_plus[i,"analysis_pathset_id"], i
                    )
                    gsea_plot_draw <- WTS_plot.Top.Rank.Profile.GSEA.plot(
                        AnalysisID      = AnalysisID,
                        InputGeneRank   = GENE_RANK,
                        PathwayID       = gsea_plot_taraget_plus[i, "path_id"],
                        PathwayName     = gsea_plot_taraget_plus[i, "pathway"],
                        PlotSave        = CONFIG$GSEA_PLOT$PLOT_SAVE,
                        DisplayOnScreen = CONFIG$GSEA_PLOT$SCREEN_DISPLAY,
                        PlotName        = gsea_plot_name,
                        FontSize        = 11
                    )
                }
                for( i in 1:nrow(gsea_plot_taraget_minus) )
                {
                    gsea_plot_name <- sprintf("%s/%s.%s.GSEA.%s.NES.minus.rank.%s.plot.png", 
                        gsea_res_dir, SEQ_FOLDER, AnalysisID, gsea_plot_taraget_minus[i,"analysis_pathset_id"], i
                    )
                    gsea_plot_draw <- WTS_plot.Top.Rank.Profile.GSEA.plot(
                        AnalysisID      = AnalysisID,
                        InputGeneRank   = GENE_RANK,
                        PathwayID       = gsea_plot_taraget_minus[i, "path_id"],
                        PathwayName     = gsea_plot_taraget_minus[i, "pathway"],
                        PlotSave        = CONFIG$GSEA_PLOT$PLOT_SAVE,
                        DisplayOnScreen = CONFIG$GSEA_PLOT$SCREEN_DISPLAY,
                        PlotName        = gsea_plot_name,
                        FontSize        = 11
                    )
                }
            }
            res_message <- sprintf("%s_gsea_plot_done.", AnalysisID)
            return(res_message)
        })
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===[ GENE-ONTOLOGY ANALYSIS ]=======================================================================================================================#
    # go enrich analysis result folder ----------------------------------------#
    go_res_dir <- sprintf("%s/go_enrich_analysis", BaseAnalysisDir)
    if( !dir.exists(go_res_dir) ){ system(sprintf("mkdir -p %s", go_res_dir)) }
    message(sprintf("GO enrich analysis result folder = %s ", go_res_dir))
    # global pathway analysis global rds folder -------------------------------#
    GO_ANALYSIS_RDS_DIR <- sprintf("%s/RDS_GoEnrichAnalysis", CONFIG$GLOBAL$RDS_DIR)
    # GO pathwayset id ---------------------------------------------------------
    GO_PathwayIDs <- c("gobp"="pw10", "gocc"="pw11", "gomf"="pw12")
    #--------------------------------------------------------------------------#
#---/ GO ENRICH ANALYSIS /-------------------------------------------------------------------------#
    if( CONFIG$GO$RUN_GO_ANALYSIS )
    {
        # GO pathwayset id check -----------------------------------------------
        GO_PATH_ID <- GO_PathwayIDs [ unlist(strsplit(CONFIG$GO$SELECTED_GENE_ONTOLOGY, ",")) ]
        if( length(GO_PATH_ID) == 0 ){ stop("||---!!! no GeneOntology pathway-set-id found. please check again. STOPPED.")}
        # load DEG analysis result ---------------------------------------------
        deg_analysis_result_list <- lapply( ANALYSIS_ID_INDEX, function(aidx) 
            load.DEG.Analysis.Profile( SeqFolder = SEQ_FOLDER, AnalysisId = aidx, DataSource = CONFIG$GO$DEG_DATA_SOURCE, GeneInfo = GENE_INFO )
        )
        # set parallel cores ---------------------------------------------------
        if( length(deg_analysis_result_list) > 5 ){ ParallelCores <- 5 }else{ ParallelCores <- length(deg_analysis_result_list) }

        # run GO enrichment analysis --------------------------------------------------------------#
        GoEnrichResultList <- mclapply( deg_analysis_result_list, function(AD)
        {
            GoEnRes <- ldply(lapply(GO_PATH_ID, function(GOID)
            {
                run.DEG.Pathway.Analysis( 
                    DegAnalysisResult = AD, 
                    PathwaySetID      = GOID, 
                    AddMergeDegRun    = CONFIG$GO$ADD_TOTAL_DEG_RUN, 
                    GeneIdentifier    = CONFIG$GO$GENE_IDENTIFIER, 
                    Threads           = CONFIG$GO$THREADS 
                )
            }))
            GoEnRes <- GoEnRes %>% dplyr::select(!c(".id"))
            return(GoEnRes)
        }, mc.cores = ParallelCores)
        #------------------------------------------------------------------------------------------#

        # >> whole results save : RDS, TEXT | top200 results save : RDS, database
        # top200 results -------------------------------------------------------
        GoEnResTop200List <- lapply(GoEnrichResultList, function(GORES) 
        {
            analysis_unit <- unique(GORES[,c("input_geneset","analysis_pathset_id")])
            analysis_unit <- analysis_unit %>% filter( !is.na(input_geneset) )
            if( nrow(analysis_unit) > 0 ){ 
                TOP200_RESULT <- data.frame()
                for( i in 1:nrow(analysis_unit) )
                {
                    top200_res <- GORES %>% filter( 
                        input_geneset        == analysis_unit[i, "input_geneset"       ],
                        analysis_pathset_id  == analysis_unit[i, "analysis_pathset_id" ]
                    ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:200)
                    TOP200_RESULT <- rbind(TOP200_RESULT, top200_res)
                    top200_res    <- data.frame() 
                }
            }else{
                TOP200_RESULT <- GORES
            }
            return(TOP200_RESULT)
        })
        # whole profiles save as RDS -------------------------------------------
        if( CONFIG$GO$WHOLE_RESULT_SAVE_RDS )
        {
            # analysis-id level whole result RDS save --------------------------
            for( k in 1:length(GoEnrichResultList) )
            {
                saveRDS( GoEnrichResultList[[k]], file=sprintf("%s/%s.%s.GO.enrichment.analysis.whole.result.table.rds", go_res_dir,          SEQ_FOLDER, names(GoEnrichResultList)[k] ) )
                saveRDS( GoEnrichResultList[[k]], file=sprintf("%s/%s.%s.GO.enrichment.analysis.whole.result.table.rds", GO_ANALYSIS_RDS_DIR, SEQ_FOLDER, names(GoEnrichResultList)[k] ) )
                message(sprintf("|--->> GO enrich analysis-id level whole result save as RDS. SeqFolder = %s AnalysisID = %s.", SEQ_FOLDER, names(GoEnrichResultList)[k]))
            }
            # seq_folder level result list RDS save ----------------------------
            saveRDS( GoEnrichResultList, file=sprintf("%s/%s.GO.enrichment.analysis.whole.result.table.rds", GO_ANALYSIS_RDS_DIR, SEQ_FOLDER ) )
            message(sprintf("|--->> GO enrich analysis whole result save as RDS. SeqFolder = %s.", SEQ_FOLDER))
        }
        # seq_folder level whole result list XLSX save -------------------------
        if( CONFIG$GO$WHOLE_RESULT_WRITE_XLSX )
        {
            AidSampleGroupName <- reshape(unique(ainfo[,c("analysis_id","analysis_group","sample_group")]), idvar='analysis_id', timevar='analysis_group', direction='wide') %>%
                mutate( labels = paste(sample_group.TREAT, sample_group.CTRL, sep="_"))
            aid_labels        <- AidSampleGroupName$labels
            names(aid_labels) <- AidSampleGroupName$analysis_id

            aid_name_list = lapply(indexing(ainfo$analysis_id), function(idx){
                data.frame(aid_name = paste0(
                    ainfo %>% filter( analysis_id == idx, analysis_group == "TREAT" ) %>% .$sample_group %>% unique(),
                    "_vs_",
                    ainfo %>% filter( analysis_id == idx, analysis_group == "CTRL" ) %>% .$sample_group %>% unique()
                ))
            }) %>% ldply(.id="analysis_id") %>% mutate(analysis_id=as.character(analysis_id))

            GoEnrichResultList2 <- lapply(GoEnrichResultList, function(PRE) {
                PRE_rev <- data.frame(
                    deg_analysis_name = aid_name_list[match(PRE$deg_analysis_id, aid_name_list$analysis_id), "aid_name"],
                    PRE %>% dplyr::select(!c("seq_folder"))
                )
                return(PRE_rev)
            })

            names(GoEnrichResultList2) <- aid_labels[ names(GoEnrichResultList2) ]
            openxlsx::write.xlsx( GoEnrichResultList2,
                sprintf("%s/%s.pathway.analysis.whole.result.table.xlsx", go_res_dir, SEQ_FOLDER ),
                rowNames=FALSE, overwrite=TRUE
            )
            message(sprintf("|--->> pathway analysis whole result list save as XLSX. SeqFolder = %s.", SEQ_FOLDER))
        }
        # top200 result save ---------------------------------------------------
        for( k in 1:length(GoEnResTop200List) )
        {
            # save RDS ---------------------------------------------------------
            if( CONFIG$GO$TOP200_RESULT_SAVE_RDS )
            {
                saveRDS( GoEnResTop200List[[k]], file=sprintf("%s/%s.%s.Go.enrichment.analysis.top200.result.table.rds", go_res_dir,          SEQ_FOLDER, names(GoEnResTop200List)[k] ))
                saveRDS( GoEnResTop200List[[k]], file=sprintf("%s/%s.%s.Go.enrichment.analysis.top200.result.table.rds", GO_ANALYSIS_RDS_DIR, SEQ_FOLDER, names(GoEnResTop200List)[k] ))
                message(sprintf("|--->> GO enrich analysis AID level top200 result save as RDS. SeqFolder = %s AnalysisID = %s.", SEQ_FOLDER, names(GoEnResTop200List)[k]))
            }
            # import into database ---------------------------------------------
            if(  CONFIG$GO$TOP200_RESULT_DB_IMPORT )
            {
                dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                DeleteExtstData <- dbGetQuery(dbCon, 
                    sprintf("DELETE FROM go_enrich_analysis_top200_result WHERE seq_folder = '%s' AND deg_analysis_id = '%s'", 
                        unique(GoEnResTop200List[[k]]$seq_folder), unique(GoEnResTop200List[[k]]$deg_analysis_id) 
                    )
                )
                UpdateData      <- dbWriteTable(dbCon, name="go_enrich_analysis_top200_result", value=GoEnResTop200List[[k]], row.names=FALSE, append=TRUE)
                dbDisconnect(dbCon)
                message(sprintf("|--->> pathway analysis AID level top200 result table imported into database. seq_folder = %s, analysis_id = %s.", 
                    unique(GoEnResTop200List[[k]]$seq_folder), unique(GoEnResTop200List[[k]]$deg_analysis_id)) 
                )
            }
        }
        # create GO enrichment analysis params ---------------------------------
        go_analysis_run_params <- do.call(rbind, 
            lapply( GoEnrichResultList, function(y) unique(y[,c("seq_folder","deg_analysis_id","input_geneset","analysis_pathset_id","analysis_subset_name")]) )
        ) %>% group_by(seq_folder, deg_analysis_id) %>% reframe(
            input_deg = paste(unique(input_geneset), collapse=","), 
            run_pathset_id = paste(unique(analysis_pathset_id), collapse=","), 
            pathway_subsets = paste(unique(analysis_subset_name), collapse=",")
        ) %>% as.data.frame()
        # GO enrichment analysis params save as text ---------------------------
        if( CONFIG$GO$ANALYSIS_PARAMS_SAVE_AS_TEXT )
        {
            write.table(
                go_analysis_run_params,
                sprintf("%s/%s.GO.enrichment.analysis.params.tsv", go_res_dir, SEQ_FOLDER),
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
            )
        }
        # GO enrichment analysis params import into database -------------------
        if( CONFIG$GO$ANALYSIS_PARAMS_DB_IMPORT )
        {
            dbCon           <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM go_enrich_analysis_params WHERE seq_folder = '%s'", SEQ_FOLDER))
            UpdateData      <- dbWriteTable(dbCon, name="go_enrich_analysis_params", value=go_analysis_run_params, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
            message(sprintf("|--->> GO enrichment analysis params imported into database. seq_folder = %s", SEQ_FOLDER) ) 
        }
    }
#--------------------------------------------------------------------------------------------------#
#---/ GO ENRICH ANALYSIS RESULT PLOTS /------------------------------------------------------------#
    # load GO enrich result data ----------------------------------------------#
    if( CONFIG$GO_PLOT$GO_DATA_SOURCE == "rds" )
    {
        go_res_data_list <- readRDS(sprintf("%s/%s.GO.enrichment.analysis.whole.result.table.rds", GO_ANALYSIS_RDS_DIR, SEQ_FOLDER))
    }else{
        go_res_data_list <- lapply( ANALYSIS_ID_INDEX, function(aid) load.GO.Enrich.Analysis.Profile(SeqFolder=SEQ_FOLDER, AnalysisId=aid, DataSource='db', RdsDir=GO_ANALYSIS_RDS_DIR ) )
    }
    #--------------------------------------------------------------------------#
    GoRootNamespace <- c("gobp"="biological_process", "gocc"="cellular_component", "gomf"="molecular_function")
    # GO Terms ----------------------------------------------------------------#
    if( CONFIG$GO_PLOT$GO_TOP_NET$GO_TERM_DATA_SOURCE == "rds" )
    {
        GO_TERMS <- readRDS("/storage/DB/GeneOntology")
        GO_TERMS <- GO_TERMS %>% filter( GO_namespace == GoRootNamespace[ CONFIG$GO$SELECTED_GENE_ONTOLOGY ] )
    }else{
        dbCon   <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db="public_database" )
        GO_TERMS <- dbGetQuery(dbCon, sprintf("SELECT * FROM go_terms WHERE GO_namespace = '%s'", GoRootNamespace[ CONFIG$GO$SELECTED_GENE_ONTOLOGY ] ))
        dbDisconnect(dbCon)
    }
    #--------------------------------------------------------------------------#
    # GO enrich top-rank dotplot ----------------------------------------------#
    GoTopRankDotplotList <- lapply( go_res_data_list, function( GODATA ) 
    {
        #print(unique(GODATA$deg_analysis_id))
        InputGoTopNet  <- GODATA %>% filter( analysis_pathset_id == GO_PathwayIDs[ CONFIG$GO_PLOT$GO_TOP_DOTPLOT$GO_ROOT ] )
        PlotFilePrefix <- sprintf("%s/%s.%s.Go.enrich.top.rank.%s", go_res_dir, SEQ_FOLDER, unique(InputGoTopNet$deg_analysis_id), CONFIG$GO_PLOT$GO_TOP_DOTPLOT$PLOT_TOP_RANK_N) 
        #
        if( nrow(InputGoTopNet) > 0 ){
            toprankgores_dotplot_up <- WTS_plot.Pathway.Analysis.TopRank.Result(
                PathwayAnalysisResult = InputGoTopNet,
                GeneType              = "UP",
                PlotType              = "dot",
                TopReults_N           = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$PLOT_TOP_RANK_N,
                AxisFontSize          = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_FONT_SIZE,
                Width                 = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_WIDTH,
                Height                = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_HEIGHT,
                PathwayNameLength     = 30,
                PlotSave              = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_SAVE,
                DisplayOnScreen       = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_SCREEN_DISPLAY,
                PlotName              = sprintf("%s.UP.genes.dotplot.png", PlotFilePrefix),
                PlotTitle             = NULL #paste0(unique(InputGoTopNet$deg_analysis_id), "_", toupper(CONFIG$GO_PLOT$GO_TOP_DOTPLOT$GO_ROOT) )
            ) 
            toprankgores_dotplot_down <- WTS_plot.Pathway.Analysis.TopRank.Result(
                PathwayAnalysisResult = InputGoTopNet,
                GeneType              = "DOWN",
                PlotType              = "dot",
                TopReults_N           = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$PLOT_TOP_RANK_N,
                AxisFontSize          = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_FONT_SIZE,
                Width                 = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_WIDTH,
                Height                = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_HEIGHT,
                PathwayNameLength     = 30,
                PlotSave              = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_SAVE,
                DisplayOnScreen       = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_SCREEN_DISPLAY,
                PlotName              = sprintf("%s.DOWN.genes.dotplot.png", PlotFilePrefix),
                PlotTitle             = NULL #paste0(unique(InputGoTopNet$deg_analysis_id), "_", toupper(CONFIG$GO_PLOT$GO_TOP_DOTPLOT$GO_ROOT) )
            )
            if( CONFIG$GO_PLOT$GO_TOP_DOTPLOT$MERGED_DOTPLOT_SAVE )
            {
                if( CONFIG$GO_PLOT$GO_TOP_DOTPLOT$MERGE_DIRECTION == 1 )
                {
                    toprankpathres_dotplot <- cowplot::plot_grid(toprankgores_dotplot_up, toprankgores_dotplot_down, nrow=1, rel_widths=c(1,1))
                    suppressWarnings(ggsave( toprankpathres_dotplot, 
                        file   = sprintf("%s.merged.dotplot.png", PlotFilePrefix),
                        width  = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_WIDTH * 2, 
                        height = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_HEIGHT, 
                        unit = 'in', dpi = 150, type = 'cairo'   
                    ))
                    message("|---> Merged Dotplot also created. Merged direction : horizontal")
                }else{
                    toprankpathres_dotplot <- cowplot::plot_grid(toprankgores_dotplot_up, toprankgores_dotplot_down, ncol=1, rel_heights=c(1,1))
                    suppressWarnings(ggsave( toprankpathres_dotplot, 
                        file   = sprintf("%s.merged.dotplot.png", PlotFilePrefix),
                        width  = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_WIDTH, 
                        height = CONFIG$GO_PLOT$GO_TOP_DOTPLOT$DOTPLOT_HEIGHT * 2, 
                        unit ='in', dpi=150, type='cairo'   
                    ))
                    message("|---> Merged Dotplot also created. Merged direction : vertical")
                }  
            }
        }
    })
    # GO enrich top-rank result as network ------------------------------------#
    GoTopRankNetworkList <- lapply( go_res_data_list, function( GODATA ) 
    {
        print(unique(GODATA$deg_analysis_id))
        InputGoTopNet <- GODATA %>% filter( analysis_pathset_id == GO_PathwayIDs[ CONFIG$GO_PLOT$GO_TOP_NET$GO_ROOT ] )
        if( nrow(InputGoTopNet) > 0 ){
            # make network data ----------------------------------------------------
            TopRankGoAsNetwork <- create.TopRank.GO.Network.Data(
                GoPathwayAnalysisResult = InputGoTopNet ,
                GoTerms                 = GO_TERMS$GO_TermID,
                TopRank_N               = CONFIG$GO_PLOT$GO_TOP_NET$TOP_RANK_N,
                UseAllDegResult         = CONFIG$GO_PLOT$GO_TOP_NET$USE_ALL_DEG_RESULT,
                GoRoot                  = CONFIG$GO_PLOT$GO_TOP_NET$GO_ROOT,
                IncludeAllTerms         = CONFIG$GO_PLOT$GO_TOP_NET$INCLUDE_ALL_TERMS,
                Threads                 = CONFIG$GO_PLOT$GO_TOP_NET$THREADS
            )
            # draw as plot ---------------------------------------------------------
            # go_top_net_plot <- WTS_plot.GO.Top.Enrich.As.Network(
            #     GoNetworDataList      = TopRankGoAsNetwork,
            #     GoEnrichResult        = InputGoTopNet,
            #     GoTerms               = GO_TERMS,
            #     IsAllDegsNetwork      = CONFIG$GO_PLOT$GO_TOP_NET$USE_ALL_DEG_RESULT,
            #     ScoreCutOff           = CONFIG$GO_PLOT$GO_TOP_NET$SCORE_CUTOFF,
            #     TopRank_N             = CONFIG$GO_PLOT$GO_TOP_NET$TOP_RANK_N,
            #     NodeShape             = CONFIG$GO_PLOT$GO_TOP_NET$NODE_SHAPE,
            #     MarkRootNode          = CONFIG$GO_PLOT$GO_TOP_NET$MARK_ROOT_NODE,
            #     MarkTopRankNodes      = CONFIG$GO_PLOT$GO_TOP_NET$MARK_TOP_RANK_NODES,
            #     NetworkLayout         = CONFIG$GO_PLOT$GO_TOP_NET$NETOWRK_LAYOUT,
            #     UseVisNetwork         = CONFIG$GO_PLOT$GO_TOP_NET$USE_VISNETWORK,
            #     VisNetHighlightDegree = CONFIG$GO_PLOT$GO_TOP_NET$VISNET_HIGHLIGHT_DEGREE
            # )
            # save plots -----------------------------------------------------------
            if( class(TopRankGoAsNetwork) != "data.frame" ){
                if( CONFIG$GO_PLOT$GO_TOP_NET$NETOWRK_SAVE_AS_PNG )
                {
                    if( CONFIG$GO_PLOT$GO_TOP_NET$USE_VISNETWORK )
                    { 
                        message("|---!!> Save network as png file NOT SUPPORTED when use visNetwork.")
                    }else{
                        if( CONFIG$GO_PLOT$GO_TOP_NET$USE_ALL_DEG_RESULT ){ DegUseTag <- "all.degs" }else{ DegUseTag <- "up.down.degs" }
                        go_top_net_png_filename <- sprintf("%s/%s.%s.GO.enrich.%s.top.rank.%s.%s.network.png", 
                            go_res_dir, SEQ_FOLDER, unique(InputGoTopNet$deg_analysis_id), CONFIG$GO_PLOT$GO_TOP_NET$GO_ROOT, CONFIG$GO_PLOT$GO_TOP_NET$TOP_RANK_N, DegUseTag 
                        )
                        go_top_net_png_tmp <- sprintf("%s/%s.%s.go.top.network.tmp.png", go_res_dir, SEQ_FOLDER, unique(InputGoTopNet$deg_analysis_id) )
                        # make png -----------------------------------------------------
                        png( go_top_net_png_tmp, width=CONFIG$GO_PLOT$GO_TOP_NET$PLOT_WIDTH, height=CONFIG$GO_PLOT$GO_TOP_NET$PLOT_HEIGHT, units="in",res=150 )
                        go_top_net_plot <- WTS_plot.GO.Top.Enrich.As.Network(
                            GoNetworDataList      = TopRankGoAsNetwork,
                            GoEnrichResult        = InputGoTopNet,
                            GoTerms               = GO_TERMS,
                            IsAllDegsNetwork      = CONFIG$GO_PLOT$GO_TOP_NET$USE_ALL_DEG_RESULT,
                            ScoreCutOff           = CONFIG$GO_PLOT$GO_TOP_NET$SCORE_CUTOFF,
                            TopRank_N             = CONFIG$GO_PLOT$GO_TOP_NET$TOP_RANK_N,
                            NodeShape             = CONFIG$GO_PLOT$GO_TOP_NET$NODE_SHAPE,
                            MarkRootNode          = CONFIG$GO_PLOT$GO_TOP_NET$MARK_ROOT_NODE,
                            MarkTopRankNodes      = CONFIG$GO_PLOT$GO_TOP_NET$MARK_TOP_RANK_NODES,
                            NetworkLayout         = CONFIG$GO_PLOT$GO_TOP_NET$NETOWRK_LAYOUT,
                            UseVisNetwork         = CONFIG$GO_PLOT$GO_TOP_NET$USE_VISNETWORK,
                            VisNetHighlightDegree = CONFIG$GO_PLOT$GO_TOP_NET$VISNET_HIGHLIGHT_DEGREE
                        )
                        dev.off()
                        # trim png -----------------------------------------------------
                        system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", go_top_net_png_tmp, go_top_net_png_filename))
                        system(sprintf("rm %s", go_top_net_png_tmp))
                    }
                }
                # display on screen ----------------------------------------------------
                if( CONFIG$GO_PLOT$GO_TOP_NET$NETWORK_DISPLAY_ON_SCREEN )
                {
                    x11( width=20, height=12, bg="white")
                    go_top_net_plot <- WTS_plot.GO.Top.Enrich.As.Network(
                        GoNetworDataList      = TopRankGoAsNetwork,
                        GoEnrichResult        = InputGoTopNet,
                        GoTerms               = GO_TERMS,
                        IsAllDegsNetwork      = CONFIG$GO_PLOT$GO_TOP_NET$USE_ALL_DEG_RESULT,
                        ScoreCutOff           = CONFIG$GO_PLOT$GO_TOP_NET$SCORE_CUTOFF,
                        TopRank_N             = CONFIG$GO_PLOT$GO_TOP_NET$TOP_RANK_N,
                        NodeShape             = CONFIG$GO_PLOT$GO_TOP_NET$NODE_SHAPE,
                        MarkRootNode          = CONFIG$GO_PLOT$GO_TOP_NET$MARK_ROOT_NODE,
                        MarkTopRankNodes      = CONFIG$GO_PLOT$GO_TOP_NET$MARK_TOP_RANK_NODES,
                        NetworkLayout         = CONFIG$GO_PLOT$GO_TOP_NET$NETOWRK_LAYOUT,
                        UseVisNetwork         = CONFIG$GO_PLOT$GO_TOP_NET$USE_VISNETWORK,
                        VisNetHighlightDegree = CONFIG$GO_PLOT$GO_TOP_NET$VISNET_HIGHLIGHT_DEGREE
                    )
                }
            }
        }
    })
    # GO bubble plot ----------------------------------------------------------#
    # GoBubbleCategory <- c("gobp"="BP", "gocc"="CC", "gomf"="MF")

    # GoPlotObject <- create.GOplot.Object(
    #     GoEnrichProfile = go_res_data_list[[1]],
    #     DegProfile      = deg_data_list[[1]],
    #     GoTermTable     = GO_TERMS,
    #     GeneInfo        = GENE_INFO,
    #     GoCategory      = GoBubbleCategory[ CONFIG$GO_PLOT$GO_BUBBLE$GO_ROOT ]
    # )

    # label_cut_off <- quantile(-log10(GoPlotObject[which(GoPlotObject$adj_pval != 0 ), "adj_pval"]), 0.97)

    # GoPlotObject_Reduced <- reduce_overlap(GoPlotObject, overlap = 0.75)

    # GOBubble( GoPlotObject, labels = label_cut_off)


    # gobubbleplot <- GOBubble( GoPlotObject_Reduced, labels = label_cut_off) + coord_cartesian(ylim = c(0,50))




#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#















