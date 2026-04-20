
#---/ PACKAGES /-----------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("parallel"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("openxlsx"))
    suppressPackageStartupMessages(library("yaml"))
#--------------------------------------------------------------------------------------------------#

#---/ PRE-DEFINE VALUES /--------------------------------------------------------------------------#
    ExcludeSeqID           <- NULL
    IncludeSeqID           <- NULL
#--------------------------------------------------------------------------------------------------#

#---/ ARGUMENTS /----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--BASE_DIR"  ),      action="store", default=NULL, type="character", help="WES analysis base folder"),
        make_option(c("--SEQ_FOLDER"),      action="store", default=NULL, type="character", help="Analysis Batch ID (SeqFolderID)"),
        make_option(c("--PRESET_RDS"),      action="store", default=NULL, type="character", help="PresetGenes containing RDS file"),
        make_option(c("--CLIENT_ID"),       action="store", default=NULL, type="integer",   help="Client ID (as integer)"),
        make_option(c("--ORDER_ID"),        action="store", default=NULL, type="character", help="Order ID"),
        make_option(c("--EXCLUDE_SEQID"),   action="store", default=NULL, type="character", help="SampleID to remove in analysis"),
        make_option(c("--INCLUDE_SEQID"),   action="store", default=NULL, type="character", help="SampleID to add in analysis"),
        make_option(c("--ONCOPLOT1_LIMIT"), action="store", default=5,    type="integer",   help="Oncoplot-1 minimum sample limit"),
        make_option(c("--ONCOPLOT2_LIMIT"), action="store", default=10,   type="integer",   help="Oncoplot-2 minimum sampel limit"),
        make_option(c("--NGS_LIB"),         action="store", default="twist", type="character",  help="WES capture library"),
        make_option(c("--TONLY_CALL_MODE"), action="store", default="Tonly", type="character",  help="Tonly variant call mode. 'Tonly' or 'TonlyGermlineFiltered' ")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BaseDir                <- ARGS$BASE_DIR
    SeqFolderID            <- ARGS$SEQ_FOLDER
    ClientPreSetGenesRDS   <- ARGS$PRESET_RDS
    ClientID               <- ARGS$CLIENT_ID
    OrderID                <- ARGS$ORDER_ID
    ExcludeSeqID           <- ARGS$EXCLUDE_SEQID
    IncludeSeqID           <- ARGS$INCLUDE_SEQID
    OncoplotSampleLimit    <- ARGS$ONCOPLOT1_LIMIT
    CumOncoplotSmapleLimit <- ARGS$ONCOPLOT2_LIMIT
    NgsLib                 <- ARGS$NGS_LIB
    TonlyCallMode          <- ARGS$TONLY_CALL_MODE
#--------------------------------------------------------------------------------------------------#

#---/ MANUAL PARAMS (FOR DEV AND DEBUG) /----------------------------------------------------------#
    # BaseDir                = "/data/wes"
    # SeqFolderID            = "WES_25_09"
    # ClientPreSetGenesRDS   = "preset.genes.rds"
    # ClientID               = 15
    # OrderID                = "GCX-C15-251126"
    # ExcludeSeqID           = NULL
    # IncludeSeqID           = NULL
    # OncoplotSampleLimit    = 5
    # CumOncoplotSmapleLimit = 10
    # NgsLib                 = 'twist'
    # TonlyCallMode          = 'Tonly'
#--------------------------------------------------------------------------------------------------#

#---/ REQUIRED PARAMS CHECK /----------------------------------------------------------------------#
    if( is.null(BaseDir)               ){ stop("|---!!! No 'BASE_DIR' found. REQUIRED. STOP.")   }
    if( is.null(SeqFolderID)           ){ stop("|---!!! No 'SEQ_FOLDER' found. REQUIRED. STOP.") }
    if( is.null(ClientPreSetGenesRDS)  ){ stop("|---!!! No 'PRESET_RDS' found. REQUIRED. STOP.") }
    if( is.null(ClientID)              ){ stop("|---!!! No 'CLIENT_ID' found. REQUIRED. STOP.")  }
    if( is.null(OrderID)               ){ stop("|---!!! No 'ORDER_ID' found. REQUIRED. STOP.")   }
    if( is.null(NgsLib)                ){ NgsLib = "twist"            }
    if( is.null(OncoplotSampleLimit)   ){ OncoplotSampleLimit = 5     }
    if( is.null(CumOncoplotSmapleLimit)){ CumOncoplotSmapleLimit = 10 }
#--------------------------------------------------------------------------------------------------#

#---/ FUNCTIONS /----------------------------------------------------------------------------------#
    indexing = function(x){ sapply(unique(x), function(z) list(z)) }
    source("/storage/home/kangsm/runScripts/WES_Fun.ReportModules.R")
    source("/storage/home/kangsm/runScripts/NGS_Fun.ReportModules_v2.R")
    #cnf <- read_yaml("/storage/home/kangsm/myScripts/Default_Scripts/wes_report_config.yaml")
    cnf <- read_yaml(sprintf("%s/%s/meta/wes_report_config.yaml", BaseDir, SeqFolderID))
    #---------------------------------------------------------------------------
    db_host=cnf$db$host
    db_user=cnf$db$user
    db_pw=cnf$db$pw
    db_port=cnf$db$port
    db_name_info=cnf$db$db_info
    db_name_data=cnf$db$db_data
#--------------------------------------------------------------------------------------------------#

#---/ PRE-LOAD DATA /------------------------------------------------------------------------------#
    source("/storage/home/kangsm/myScripts/Default_Scripts/ruo_wes_db.R")
    dbCon <- dbConnect(dbDriver("MySQL"), host=cnf$db$host, user=cnf$db$user, port=cnf$db$port, password=cnf$db$pw, db=cnf$db$db_info )
    pinfo <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqfolder_info WHERE seq_folder = '%s'", SeqFolderID))
    sinfo <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info     WHERE seq_folder = '%s'", SeqFolderID))
    if( !is.null(ExcludeSeqID) )
    {
        exclude_samples <- unlist(strsplit(ExcludeSeqID, ","))
        sinfo           <- sinfo %>% filter( seq_id %nin% exclude_samples )
    }
    if( !is.null(IncludeSeqID) )
    {
        include_samples   <- unlist(strsplit(IncludeSeqID, ","))
        QueryIncludeSeqID <- paste(paste0("'", include_samples, "'"), collapse=",")
        sinfo_extra <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_id IN (%s)", QueryIncludeSeqID))
        sinfo <- rbind(sinfo, sinfo_extra)
    }    
    dbDisconnect(dbCon)
    #---/ ORDER ID /------------------------------------------------------------
    if( is.null(OrderID) )
    {
        if( nrow(pinfo) == 1 ){ OrderID <- pinfo$ngs_order_id }else{ OrderID <- "NO_ORDER_ID" }
    }else{
        OrderID <- OrderID
    }
    #---/ SAMPLE INFO ORDERING /------------------------------------------------
    sinfo[which(is.na(sinfo$sample_group)), "sample_group"] <- ""
    sinfo <- rbind(
        sinfo %>% filter( sample_group != "" ) %>% arrange( sample_info, sample_group, seq_id ),
        sinfo %>% filter( sample_group == "" ) %>% arrange( sample_info, seq_id )
    )
#--------------------------------------------------------------------------------------------------#

#---/ DATA FOLDERS AND RESULT FOLDERS /------------------------------------------------------------#
    # analysis result dir ------------------------------------------------------
    ResultDir <- ifelse(
        TonlyCallMode == "Tonly",
        sprintf("%s/%s/analysis", BaseDir, SeqFolderID),
        sprintf("%s/%s/analysis_GF", BaseDir, SeqFolderID)
    )
    if( ! dir.exists(ResultDir) ){ system(sprintf("mkdir -p %s", ResultDir)) }
    # client export data dir ---------------------------------------------------
    ExportDir <- sprintf("%s/export", ResultDir)
    if( ! dir.exists(ExportDir) ){ system(sprintf("mkdir -p %s", ExportDir)) }
    # report data dir ----------------------------------------------------------
    WES_REPORT_DIR <- ifelse(
        TonlyCallMode == "Tonly",
        sprintf("/storage/home/kangsm/shinyWeb/REPORT_WES/%s_%s", OrderID, SeqFolderID),
        sprintf("/storage/home/kangsm/shinyWeb/REPORT_WES/%s_%s_GermlineFiltered", OrderID, SeqFolderID)
    )    
    if( ! dir.exists(WES_REPORT_DIR) ){ system(sprintf("mkdir -p %s", WES_REPORT_DIR)) }
    # resource dir symbolic link -----------------------------------------------
    resource_dir <- "/storage/home/kangsm/shinyWeb/resources"
    system(sprintf("ln -Tsf %s %s/resources", resource_dir, WES_REPORT_DIR))
    # report data dir ----------------------------------------------------------
    report_data_dir <- sprintf("%s/data", WES_REPORT_DIR)
    if( ! dir.exists(report_data_dir) ){ system(sprintf("mkdir -p %s", report_data_dir)) }
#--------------------------------------------------------------------------------------------------#

#---/ ANALYSIS RUN PARAMS /------------------------------------------------------------------------#
    
    message("|--->> Prepare Analysis Parameters ... \n")
    
    SampleGroupList <- unique(sinfo$sample_group[!is.na(sinfo$sample_group)])
    #---/ Sample Groups /-------------------------------------------------------
    if( length(SampleGroupList) == 0 )
    { 
        WithinGroupAnalysis <- FALSE 
    }else{
        WithinGroupAnalysis <- TRUE
        AnalysisGroupIndex  <- indexing(SampleGroupList)
    }
    #---/ Within Group Sample Comparison /--------------------------------------
    if( WithinGroupAnalysis )
    {
        AnalysisSetInfoList <- lapply( AnalysisGroupIndex, function(idx)
        {
            GroupSamples <- sinfo %>% filter( sample_group == idx )
            SampleTypes  <- GroupSamples$matched_normal
            # tumor-only list --------------------------------------------------
            if( length(SampleTypes) >= 2 ){ 
                Do_Tonly_Compare      <- TRUE 
                Do_Tonly_Compare_List <- as.data.frame(t(combn(GroupSamples$seq_id, 2))) %>% dplyr::rename(s1=1, s2=2)
            }else{ 
                Do_Tonly_Compare      <- FALSE 
                Do_Tonly_Compare_List <- NULL
            }
            # if contain 'matched-normal tissue' -------------------------------
            if( '3' %in% SampleTypes ) 
            {
                if( nrow(GroupSamples[which(GroupSamples$matched_normal != '3'), ]) >= 2 ) 
                {
                    Do_MatchCall_N_TS_Compare      <- TRUE
                    MatchCall_N_TS_Samples         <- GroupSamples %>% filter( matched_normal != '3' )
                    Do_MatchCall_N_TS_Compare_List <- as.data.frame(t(combn(MatchCall_N_TS_Samples$seq_id, 2))) %>% dplyr::rename(s1=1, s2=2)
                }else{
                    Do_MatchCall_N_TS_Compare      <- FALSE
                    Do_MatchCall_N_TS_Compare_List <- NULL
                }
            }else{
                Do_MatchCall_N_TS_Compare      <- FALSE
                Do_MatchCall_N_TS_Compare_List <- NULL
            }
            # if contain 'matched normal organoid' -----------------------------
            if( '2' %in% SampleTypes )
            {
                if( nrow(GroupSamples[which(GroupSamples$matched_normal %nin% c('2','3')), ]) >= 2 )
                {
                    Do_MatchCall_N_ORG_Compare     <- TRUE
                    MatchCall_N_ORG_Samples        <- GroupSamples %>% filter( matched_normal %nin% c('2','3') )
                    Do_MatchCall_N_TS_Compare_List <- as.data.frame(t(combn(MatchCall_N_ORG_Samples$seq_id, 2))) %>% dplyr::rename(s1=1, s2=2)
                }else{
                    Do_MatchCall_N_ORG_Compare     <- FALSE
                    Do_MatchCall_N_ORG_Compare_List <- NULL
                }
            }else{
                Do_MatchCall_N_ORG_Compare     <- FALSE
                Do_MatchCall_N_ORG_Compare_List <- NULL
            }
            # comparison analysis list -----------------------------------------
            comparison_set_info <- list(
                'tumor_only'            = list( analysis_run = Do_Tonly_Compare,           analysis_run_list = Do_Tonly_Compare_List           ),
                'match_normal_tissue'   = list( analysis_run = Do_MatchCall_N_TS_Compare,  analysis_run_list = Do_MatchCall_N_TS_Compare_List  ),
                'match_normal_organoid' = list( analysis_run = Do_MatchCall_N_ORG_Compare, analysis_run_list = Do_MatchCall_N_ORG_Compare_List )
            )
            #-------------------------------------------------------------------
            return(comparison_set_info)
        })
        AnalysisSetInfoTable <- ldply(lapply(AnalysisSetInfoList, function(y)
        {
            data.frame( 
                TumorOnly = ifelse(y$tumor_only$analysis_run,            TRUE, FALSE),
                Match_TS  = ifelse(y$match_normal_tissue$analysis_run,   TRUE, FALSE),
                Match_ORG = ifelse(y$match_normal_organoid$analysis_run, TRUE, FALSE)
            )
        }), .id = 'SampleGroup')
        AnalysisSetInfoTable$cancer_code <- sinfo[match(AnalysisSetInfoTable$SampleGroup, sinfo$sample_group), "sample_info"]
        saveRDS(AnalysisSetInfoList,  file=sprintf("%s/wes_within_group_analysis_params_list.rds",  ResultDir))
        saveRDS(AnalysisSetInfoTable, file=sprintf("%s/wes_within_group_analysis_params_table.rds", ResultDir))
    }
    #---------------------------------------------------------------------------
    AnalysisSetInfoTable = AnalysisSetInfoTable %>% arrange( cancer_code, SampleGroup)
    AnalysisTonly        = AnalysisSetInfoTable[which(AnalysisSetInfoTable$TumorOnly), ]
    AnalysisMatchTS      = AnalysisSetInfoTable[which(AnalysisSetInfoTable$Match_TS),  ]  
    AnalysisMatchORG     = AnalysisSetInfoTable[which(AnalysisSetInfoTable$Match_ORG), ]  
    #---/ preset variants t-only samples /--------------------------------------
    sinfo_preset_vars_tonly = sinfo %>% filter( sample_group %nin% c(AnalysisMatchTS$SampleGroup, AnalysisMatchORG$SampleGroup) )
    #---/ Preset List By Cancer Type /------------------------------------------
    if( !is.null(ClientPreSetGenesRDS) )
    {
        # load gene presets ----------------------------------------------------
        ClientPresetGenes <- readRDS("/storage/home/kangsm/shinyWeb/resources/preset.genes.rds")
        # preset list to analysis ----------------------------------------------
        PresetList <- ldply(lapply(indexing(sinfo$sample_info), function(cc) 
        {
            preset_list <- c(ifelse( cc %in% names(ClientPresetGenes), ClientPresetGenes$cancer_code_table %>% filter( cancer_code == cc ) %>% .$cancer_tissue, "pan_cancers" ), "mutation_panel")
            samples_n <- nrow(sinfo %>% filter( sample_info == cc ))
            cancer_code_preset <- data.frame( 
                cancer_code = cc, 
                preset_name = preset_list, 
                samples_n   = samples_n,
                heatmap     = ifelse( samples_n >= 2, TRUE, FALSE )
            ) %>% mutate( index = paste( cancer_code, preset_name, sep="_"))
            cancer_code_preset$page_break <- c(FALSE, TRUE)
            return(cancer_code_preset)
        }))
        saveRDS(PresetList, file=sprintf("%s/preset_list_table.rds", ResultDir))
        system(sprintf("cp %s/preset_list_table.rds %s/preset_list_table.rds", ResultDir, report_data_dir))
        saveRDS(ClientPresetGenes, file=sprintf("%s/client_preset_gene_list.rds", ResultDir))
        system(sprintf("cp %s/client_preset_gene_list.rds %s/client_preset_gene_list.rds", ResultDir, report_data_dir))
        message("|----->> Analysis-Gene-Presets saved as RDS and copied for report.")
    }

    #---/ analysis params RDS /-------------------------------------------------
    analysis_run_info <- list(
        Tonly    = AnalysisTonly,
        MatchTS  = AnalysisMatchTS,
        MatchORG = AnalysisMatchORG
    )
    saveRDS( analysis_run_info, file=sprintf("%s/analysis_run_info.rds", ResultDir))
    system(sprintf("cp %s/analysis_run_info.rds %s/analysis_run_info.rds", ResultDir, report_data_dir))
    message("|----->> Analysis-Run-Info saved as RDS and copied for report.")

#--------------------------------------------------------------------------------------------------#

#---/ DATA LOAD FROM DATABASE /---------------------------------------------------------------------
    message("|--->> Load data from database ... \n")
    AnalysisSeqID        <- paste(paste0("'", sinfo$seq_id, "'"), collapse=",")
    dbCon                <- dbConnect(dbDriver("MySQL"), host=cnf$db$host, user=cnf$db$user, port=cnf$db$port, password=cnf$db$pw, db=cnf$db$db_data )
    data_QC              <- dbGetQuery(dbCon, sprintf("SELECT * FROM qc_report        WHERE seq_id IN (%s)", AnalysisSeqID))
    data_VarStatSummary  <- dbGetQuery(dbCon, sprintf("SELECT * FROM variants_summary WHERE seq_id IN (%s)", AnalysisSeqID))
    data_SomaticVariants <- dbGetQuery(dbCon, sprintf("SELECT * FROM variants_somatic WHERE seq_id IN (%s)", AnalysisSeqID))
    data_IndivMatch      <- dbGetQuery(dbCon, sprintf("SELECT * FROM indiv_match      WHERE seq_folder = '%s'", SeqFolderID))
    dbDisconnect(dbCon)
    if( TonlyCallMode == "Tonly" )
    {
        # data_VarStatSummary  = data_VarStatSummary  %>% filter( variant_call_mode == "Tonly" )
        # data_SomaticVariants = data_SomaticVariants %>% filter( variant_call_mode == "Tonly" )
        data_VarStatSummary  = data_VarStatSummary  %>% filter( variant_call_mode != "TonlyGermlineFiltered" )
        data_SomaticVariants = data_SomaticVariants %>% filter( variant_call_mode != "TonlyGermlineFiltered" )
    }else{
        # data_VarStatSummary  = data_VarStatSummary  %>% filter( variant_call_mode == "TonlyGermlineFiltered" )
        # data_SomaticVariants = data_SomaticVariants %>% filter( variant_call_mode == "TonlyGermlineFiltered" )
        data_VarStatSummary  = data_VarStatSummary  %>% filter( variant_call_mode != "Tonly" )
        data_SomaticVariants = data_SomaticVariants %>% filter( variant_call_mode != "Tonly" )
    }
#---------------------------------------------------------------------------------------------------

#---/ FIGURE LIST TABLE /---------------------------------------------------------------------------
    FIGUER_LIST <- data.frame()
#---------------------------------------------------------------------------------------------------

#---/ NGS STATS PLOTS /-----------------------------------------------------------------------------
    #ngs_qc_stat_plot_height <- as.integer(nrow(sinfo)/4)+1
    ngs_qc_stat_plot_height <- nrow(sinfo) * 0.5 + 0.5
    if( ngs_qc_stat_plot_height < 2 ){ ngs_qc_stat_plot_height <- 2 }
    if( ngs_qc_stat_plot_height > 10 ){ ngs_qc_stat_plot_height <- 10 }

    #---/ GC Contents Barplot /------------------------------------------------- 
        GcRatioPlot <- NGS_report.GC.Contents.Barplot( SeqFolderId=SeqFolderID, SampleInfo=sinfo )
        tmp_png <- sprintf("%s/gc_contents_plot_tmp.png", ResultDir )
        fig_png <- sprintf("%s/gc_contents_plot.png",     ResultDir )           
        ggsave( GcRatioPlot, file=tmp_png, width=16, height=ngs_qc_stat_plot_height, unit='in', dpi=200 )
        system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
        system(sprintf("convert %s -resize %s %s", fig_png, "50%", fig_png))
        # fig list 
        FIGUER_LIST <- rbind(
            FIGUER_LIST,
            data.frame( 
                report_category = "gc_contents",
                variant_call    = "",
                cancer_code     = "" ,
                preset_name     = "",
                sample_group    = "",
                png             = "data/gc_contents_plot.png"
            )
        )
        system(sprintf("cp %s %s/gc_contents_plot.png", fig_png, report_data_dir))
        message(sprintf("|----->> NGS QC ... GC-contents barplot created and PNG copied for report"))
        system(sprintf("rm %s", tmp_png))
    #---------------------------------------------------------------------------

    #---/ Duplicates Rate Barplot /---------------------------------------------
        DupRatePlot <- NGS_report.Duplicates.Ratio.Barplot( SeqFolderId=SeqFolderID, SampleInfo=sinfo, Application="WES" )
        tmp_png <- sprintf("%s/dup_ratio_plot_tmp.png", ResultDir )
        fig_png <- sprintf("%s/dup_ratio_plot.png",     ResultDir )           
        ggsave( DupRatePlot, file=tmp_png, width=16, height=ngs_qc_stat_plot_height, unit='in', dpi=200 )
        system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
        system(sprintf("convert %s -resize %s %s", fig_png, "50%", fig_png))
        # fig list 
        FIGUER_LIST <- rbind(
            FIGUER_LIST,
            data.frame( 
                report_category = "dup_ratio",
                variant_call    = "",
                cancer_code     = "" ,
                preset_name     = "",
                sample_group    = "",
                png             = "data/dup_ratio_plot.png"
            )
        )
        system(sprintf("cp %s %s/dup_ratio_plot.png", fig_png, report_data_dir))
        message(sprintf("|----->> NGS QC ... Duplicates ratio barplot created and PNG copied for report"))
        system(sprintf("rm %s", tmp_png))
    #---------------------------------------------------------------------------

    #---/ Alignment Stats Barplot /---------------------------------------------
        AlignStatPlot <- WES_analysis.Alignment.Stats.Barplot( SeqFolderId=SeqFolderID, SampleInfo=sinfo )
        tmp_png <- sprintf("%s/align_stats_plot_tmp.png", ResultDir )
        fig_png <- sprintf("%s/align_stats_plot.png",     ResultDir )           
        ggsave( AlignStatPlot, file=tmp_png, width=16, height=ngs_qc_stat_plot_height, unit='in', dpi=200 )
        system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
        system(sprintf("convert %s -resize %s %s", fig_png, "50%", fig_png))
        # fig list 
        FIGUER_LIST <- rbind(
            FIGUER_LIST,
            data.frame( 
                report_category = "align_stats",
                variant_call    = "",
                cancer_code     = "" ,
                preset_name     = "",
                sample_group    = "",
                png             = "data/align_stats_plot.png"
            )
        )
        system(sprintf("cp %s %s/align_stats_plot.png", fig_png, report_data_dir))
        message(sprintf("|----->> NGS QC ... Alignment stats barplot created and PNG copied for report"))
        system(sprintf("rm %s", tmp_png))
    #---------------------------------------------------------------------------

    #---/ Target Coverage /-----------------------------------------------------
        TargetCovPlot <- WES_analysis.Target.Coverage.Plot( SeqFolderId=SeqFolderID, SampleInfo=sinfo )
        tmp_png <- sprintf("%s/target_cov_plot_tmp.png", ResultDir )
        fig_png <- sprintf("%s/target_cov_plot.png",     ResultDir ) 
        ggsave( TargetCovPlot, file=tmp_png, width=15, height=10, unit='in', dpi=200 )
        system(sprintf("convert %s -resize %s %s", tmp_png, "50%", fig_png))
        # fig list 
        FIGUER_LIST <- rbind(
            FIGUER_LIST,
            data.frame( 
                report_category = "target_coverage",
                variant_call    = "",
                cancer_code     = "" ,
                preset_name     = "",
                sample_group    = "",
                png             = "data/target_cov_plot.png"
            )
        )
        system(sprintf("cp %s %s/target_cov_plot.png", fig_png, report_data_dir))
        message(sprintf("|----->> NGS QC ... Target coverage plot created and PNG copied for report"))
        system(sprintf("rm %s", tmp_png))
    #---------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------

#---/ ANALYSIS : PRESET GENES VARIANTS /------------------------------------------------------------
    message( sprintf('\n|--->> Run "PRESET GENES VARIANTS" analysis <<%s| \n', paste(rep("-", 100), collapse="") ))    

    #---/ TUMOR ONLY /----------------------------------------------------------
    if(  nrow(sinfo_preset_vars_tonly) > 0 )
    {
        message("\n|----->> Tumor-Only Variants ... ")

        PresetList_Tonly            <- PresetList[which(PresetList$cancer_code %in% sinfo_preset_vars_tonly$sample_info), ]
        PresetVarAnalysisIndexTonly <- indexing( PresetList_Tonly$index )
        VariantCallMethod           <- ifelse( TonlyCallMode == "Tonly", "tumor_only", "tumor_only_gerline_filtered" )
        # preset gene variants table 
        PresetGeneVarTableListTonly <- lapply( PresetVarAnalysisIndexTonly, function(idx) 
        {
            PSGV_TumorOnly <- WES_analysis.PresetGeneVariantsTable( 
                VariantData       = data_SomaticVariants, 
                CancerCode        = PresetList_Tonly[which(PresetList_Tonly$index == idx), "cancer_code"], 
                SampleInfo        = sinfo_preset_vars_tonly, 
                VariantCallMethod = VariantCallMethod, 
                PresetGeneList    = ClientPresetGenes[[ PresetList_Tonly[which(PresetList_Tonly$index == idx), "preset_name"] ]]
            )
            return(PSGV_TumorOnly)
        })
        message("\n|----->> Tumor-Only Variants ... preset genes analyzed.")
        
        # export result for report
        PresetGeneVarTableListTonlyReport <- lapply(PresetGeneVarTableListTonly, function(pgv)
        {
            if( nrow(pgv) > 0 )
            {
                res = pgv %>% dplyr::select(c("sample_group","sample_name","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) 
            }else{
                res = data.frame()
            }
            return(res)
        })
        saveRDS( PresetGeneVarTableListTonlyReport, file=sprintf("%s/preset_gene_variants_list_Tonly.rds", ResultDir))
        system(sprintf("cp %s/preset_gene_variants_list_Tonly.rds %s/preset_gene_variants_list_Tonly.rds", ResultDir, report_data_dir))
        message("\n|----->> Tumor-Only Variants ... result saved as RDS and copied for report.")
        
        # preset gene variants export as xlsx
        PresetGeneVarTableListTonlyExport <- lapply( PresetGeneVarTableListTonly, function(pgv)
        {
            if( nrow(pgv) > 0 )
            {
                res = pgv %>% dplyr::select(c("sample_group","sample_name","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) 
                res[which(res$sample_group == "No Grouped"), "sample_group"] = ""
                colnames(res) = c("Sample.Group","Sample.ID","Gene","DNA.Change","PROTEIN.Change","Depth","Tumor.Alt.Count","VAF")
            }else{
                res = data.frame()
            }
            return(res)
        })
        PresetVarResCheck = lapply(PresetGeneVarTableListTonlyExport, function(y) ifelse( nrow(y) == 0, FALSE, TRUE ))
        PresetGeneVarTableListTonlyExport = PresetGeneVarTableListTonlyExport[ names(which(unlist(PresetVarResCheck))) ]
        write.xlsx( PresetGeneVarTableListTonlyExport, sprintf("%s/%s.PresetGenes.TumorOnlyCall.Variants.xlsx", ResultDir, OrderID))
        system(sprintf("cp %s/%s.PresetGenes.TumorOnlyCall.Variants.xlsx %s/%s.PresetGenes.TumorOnlyCall.Variants.xlsx", ResultDir, OrderID, ExportDir, OrderID))
        message("\n|----->> Tumor-Only Variants ... result saved as XLSX for client.")
        
        # preset gene variants heatmap
        pvhs_chk = lapply(PresetGeneVarTableListTonly, function(y) ifelse(length(y$sample_name) >= 2, TRUE, FALSE) )
        heatmap_target_presets = names(which(unlist(pvhs_chk)))
        HeatmapPresetList_Tonly = PresetList_Tonly %>% filter( index %in% heatmap_target_presets )
        preset_n_chk = table(HeatmapPresetList_Tonly$cancer_code)
        HeatmapPresetList_Tonly[which(HeatmapPresetList_Tonly$cancer_code %in% names(which(preset_n_chk == 1))), "page_break"] = TRUE

        # preset gene variants heatmap
        for( k in 1:nrow(HeatmapPresetList_Tonly) )
        {
            PG_VARS      <- PresetGeneVarTableListTonly[[ HeatmapPresetList_Tonly[k, "index"] ]] 
            cancer_sinfo <- sinfo[which(sinfo$sample_info == HeatmapPresetList_Tonly[k, "cancer_code"]), ]
            if( all(nrow(PG_VARS) > 0 & nrow(cancer_sinfo) >= 2) )
            {
                PGV_HM  <- WES_analysis.Preset.Gene.Variants.Heatmap( PresetVariantsTable=PG_VARS, CancerSampleInfo=cancer_sinfo )
                tmp_png <- sprintf("%s/preset_variants_heatmap_Tonly_%s_tmp.png", ResultDir, HeatmapPresetList_Tonly[k, "index"] ) 
                fig_png <- sprintf("%s/preset_variants_heatmap_Tonly_%s.png",     ResultDir, HeatmapPresetList_Tonly[k, "index"] ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(PGV_HM)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "preset_variants_heatmap",
                        variant_call    = "Tumor Only",
                        cancer_code     = HeatmapPresetList_Tonly[k, "cancer_code"],
                        preset_name     = HeatmapPresetList_Tonly[k, "preset_name"],
                        sample_group    = "",
                        png             = sprintf("data/preset_variants_heatmap_Tonly_%s.png", HeatmapPresetList_Tonly[k, "index"])
                    )
                )
                system(sprintf("cp %s %s/preset_variants_heatmap_Tonly_%s.png", fig_png, report_data_dir, HeatmapPresetList_Tonly[k, "index"] ))
                message(sprintf("\n|----->> Tumor-Only Variants ... heatmap created and PNG copied for report : %s", HeatmapPresetList_Tonly[k, "index"] ))
                system(sprintf("rm %s", tmp_png))
            }else{
                next
            }
        }
    }else{
        message("|----->> Tumor-Only Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched TS /----------------------------------------------------------
    if(  nrow(AnalysisMatchTS) > 0 )
    {
        message("|----->> Matched-Normal TISSUE Variants ... ")

        PresetList_MatchTS            <- PresetList[which(PresetList$cancer_code %in% AnalysisMatchTS$cancer_code), ]
        PresetVarAnalysisIndexMatchTS <- indexing( PresetList_MatchTS$index )

        PresetGeneVarTableListMatchTS <- lapply( PresetVarAnalysisIndexMatchTS, function(idx) 
        {
            PSGV_MatchTS <- WES_analysis.PresetGeneVariantsTable( 
                VariantData       = data_SomaticVariants, 
                CancerCode        = PresetList_MatchTS[which(PresetList_MatchTS$index == idx), "cancer_code"], 
                SampleInfo        = sinfo, 
                VariantCallMethod = "match_normal_tissue", 
                PresetGeneList    = ClientPresetGenes[[ PresetList_MatchTS[which(PresetList_MatchTS$index == idx), "preset_name"] ]]
            )
            return(PSGV_MatchTS)
        })
        message("|----->> Matched-Normal TISSUE Variants ... preset genes analyzed.")
        
        # export result for report
        PresetGeneVarTableListMatchTSReport <- lapply(PresetGeneVarTableListMatchTS, function(pgv)
        {
            res <- pgv %>% dplyr::select(c("sample_group","sample_name","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) 
            return(res)
        })
        saveRDS( PresetGeneVarTableListMatchTSReport, file=sprintf("%s/preset_gene_variants_list_Matched_TS.rds", ResultDir))
        system(sprintf("cp %s/preset_gene_variants_list_Matched_TS.rds %s/preset_gene_variants_list_MatchTS.rds", ResultDir, report_data_dir))
        message("|----->> Matched-Normal TISSUE Variants ... result saved as RDS and copied for report.")

        # preset gene variants export as xlsx
        PresetGeneVarTableListMatchTSExport <- lapply( PresetGeneVarTableListMatchTS, function(pgv)
        {
            res <- pgv %>% dplyr::select(c("sample_group","sample_name","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) 
            res[which(res$sample_group == "No Grouped"), "sample_group"] <- ""
            colnames(res) <- c("Sample.Group","Sample.ID","Gene","DNA.Change","PROTEIN.Change","Depth","Tumor.Alt.Count","VAF")
            return(res)
        })
        write.xlsx( PresetGeneVarTableListMatchTSExport, sprintf("%s/%s.PresetGenes.MatchNormalTissueCall.Variants.xlsx", ResultDir, OrderID))
        system(sprintf("cp %s/%s.PresetGenes.MatchNormalTissueCall.Variants.xlsx %s/%s.PresetGenes.MatchNormalTissueCall.Variants.xlsx", ResultDir, OrderID, ExportDir, OrderID))
        message("|----->> Matched-Normal TISSUE Variants ... result saved as XLSX for client.")

        # preset gene variants heatmap
        for( k in 1:nrow(PresetList_MatchTS) )
        {
            PG_VARS      <- PresetGeneVarTableListMatchTS[[ PresetList_MatchTS[k, "index"] ]] 
            cancer_sinfo <- sinfo[which(sinfo$sample_info == PresetList_MatchTS[k, "cancer_code"]), ]
            if( all(nrow(PG_VARS) > 0 & nrow(cancer_sinfo) >= 2) )
            {
                PGV_HM  <- WES_analysis.Preset.Gene.Variants.Heatmap( PresetVariantsTable=PG_VARS, CancerSampleInfo=cancer_sinfo )
                tmp_png <- sprintf("%s/preset_variants_heatmap_MatchTS_%s_tmp.png", ResultDir, PresetList_MatchTS[k, "index"] ) 
                fig_png <- sprintf("%s/preset_variants_heatmap_MatchTS_%s.png", ResultDir, PresetList_MatchTS[k, "index"] ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(PGV_HM)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "preset_variants_heatmap",
                        variant_call    = "Matched Normal Tissue",
                        cancer_code     = PresetList_MatchTS[k, "cancer_code"],
                        preset_name     = PresetList_MatchTS[k, "preset_name"],
                        sample_group    = "",
                        png             = sprintf("data/preset_variants_heatmap_MatchTS_%s.png", PresetList_MatchTS[k, "index"])
                    )
                )
                system(sprintf("cp %s %s/preset_variants_heatmap_MatchTS_%s.png", fig_png, report_data_dir, PresetList_MatchTS[k, "index"] ))
                message(sprintf("|----->>  Matched-Normal TISSUE Variants ... heatmap created and PNG copied for report : %s", PresetList_MatchTS[k, "index"] ))
                system(sprintf("rm %s", tmp_png))
            }else{
                next
            }
        }
    }else{
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched ORG /---------------------------------------------------------
    if(  nrow(AnalysisMatchORG) > 0 )
    {
        message("|----->> Matched-Normal ORGANOID Variants ... ")

        PresetList_MatchORG            <- PresetList[which(PresetList$cancer_code %in% AnalysisMatchORG$cancer_code), ]
        PresetVarAnalysisIndexMatchORG <- indexing( PresetList_MatchORG$index )

        PresetGeneVarTableListMatchORG <- lapply( PresetVarAnalysisIndexMatchORG, function(idx) 
        {
            PSGV_MatchORG <- WES_analysis.PresetGeneVariantsTable( 
                VariantData       = data_SomaticVariants, 
                CancerCode        = PresetList_MatchORG[which(PresetList_MatchORG$index == idx), "cancer_code"], 
                SampleInfo        = sinfo, 
                VariantCallMethod = "match_normal_organoid", 
                PresetGeneList    = ClientPresetGenes[[ PresetList_MatchORG[which(PresetList_MatchORG$index == idx), "preset_name"] ]]
            )
            return(PSGV_MatchORG)
        })
        message("|----->> Matched-Normal ORGANOID Variants ... preset genes analyzed.")

        # export result for report
        PresetGeneVarTableListMatchORGReport <- lapply(PresetGeneVarTableListMatchORG, function(pgv)
        {
            res <- pgv %>% dplyr::select(c("sample_group","sample_name","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) 
            return(res)
        })
        saveRDS( PresetGeneVarTableListMatchORGReport, file=sprintf("%s/preset_gene_variants_list_Matched_ORG.rds", ResultDir))
        system(sprintf("cp %s/preset_gene_variants_list_Matched_ORG.rds %s/preset_gene_variants_list_MatchORG.rds", ResultDir, report_data_dir))
        message("|----->> Matched-Normal ORGANOID Variants ... result saved as RDS and copied for report.")
        # preset gene variants export as xlsx
        PresetGeneVarTableListMatchORGExport <- lapply( PresetGeneVarTableListMatchORG, function(pgv)
        {
            res <- pgv %>% dplyr::select(c("sample_group","sample_name","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) 
            res[which(res$sample_group == "No Grouped"), "sample_group"] <- ""
            colnames(res) <- c("Sample.Group","Sample.ID","Gene","DNA.Change","PROTEIN.Change","Depth","Tumor.Alt.Count","VAF")
            return(res)
        })
        write.xlsx( PresetGeneVarTableListMatchORGExport, sprintf("%s/%s.PresetGenes.MatchNormalOrganoidCall.Variants.xlsx", ResultDir, OrderID))
        system(sprintf("cp %s/%s.PresetGenes.MatchNormalOrganoidCall.Variants.xlsx %s/%s.PresetGenes.MatchNormalOrganoidCall.Variants.xlsx", ResultDir, OrderID, ExportDir, OrderID))
        message("|----->> Matched-Normal ORGANOID Variants ... result saved as XLSX for client.")

        # preset gene variants heatmap
        for( k in 1:nrow(PresetList_MatchORG) )
        {
            PG_VARS      <- PresetGeneVarTableListMatchORG[[ PresetList_MatchORG[k, "index"] ]] 
            cancer_sinfo <- sinfo[which(sinfo$sample_info == PresetList_MatchORG[k, "cancer_code"]), ]
            if( all(nrow(PG_VARS) > 0 & nrow(cancer_sinfo) >= 2) )
            {
                PGV_HM  <- WES_analysis.Preset.Gene.Variants.Heatmap( PresetVariantsTable=PG_VARS, CancerSampleInfo=cancer_sinfo )
                tmp_png <- sprintf("%s/preset_variants_heatmap_MatchORG_%s_tmp.png", ResultDir, PresetList_MatchORG[k, "index"] ) 
                fig_png <- sprintf("%s/preset_variants_heatmap_MatchORG_%s.png", ResultDir, PresetList_MatchORG[k, "index"] ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(PGV_HM)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "preset_variants_heatmap",
                        variant_call    = "Matched Normal Organoid",
                        cancer_code     = PresetList_MatchORG[k, "cancer_code"],
                        preset_name     = PresetList_MatchORG[k, "preset_name"],
                        sample_group    = "",
                        png             = sprintf("data/preset_variants_heatmap_MatchORG_%s.png", PresetList_MatchORG[k, "index"])
                    )
                )
                system(sprintf("cp %s %s/preset_variants_heatmap_MatchTS_%s.png", fig_png, report_data_dir, PresetList_MatchORG[k, "index"] ))
                message(sprintf("|----->>  Matched-Normal ORGANOID Variants ... heatmap created and PNG copied for report : %s", PresetList_MatchORG[k, "index"] ))
                system(sprintf("rm %s", tmp_png))
            }else{
                next
            }
        }
    }else{
        message("|----->> Matched-Normal ORGANOID Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------

#---/ ANALYSIS : CONCRDANCE /-----------------------------------------------------------------------
    message( sprintf('\n|--->> Run "VARIANT CONCORDANCE in SAMPLE GROUP" analysis <<%s| \n', paste(rep("-", 100), collapse="") ))

    #---/ TUMOR ONLY /----------------------------------------------------------
    AnalysisTonly_concordance <- AnalysisTonly %>% filter( SampleGroup != "" )
    if(  nrow(AnalysisTonly_concordance) > 0 )
    {
        ConcordanceResult_Tonly <- list()
        for( k in 1:nrow(AnalysisTonly_concordance) )
        {
            analysis_list <- AnalysisSetInfoList[[ AnalysisTonly_concordance[k, "SampleGroup"] ]]$tumor_only$analysis_run_list
            # concordacne table
            concordance_res <- data.frame()
            for( i in 1:nrow(analysis_list) )
            {
                v1 <- data_SomaticVariants %>% filter( variant_call_mode == TonlyCallMode, seq_id == analysis_list[i, "s1"] ) %>% .$vkey
                v2 <- data_SomaticVariants %>% filter( variant_call_mode == TonlyCallMode, seq_id == analysis_list[i, "s2"] ) %>% .$vkey
                tm    <- WES_analysis.Get.Tanimoto( SET1 = v1, SET2 = v2 )
                rate1 <- WES_analysis.Get.Common.Rate( SET1 = v1, SET2 = v2, p = 1 )
                rate2 <- WES_analysis.Get.Common.Rate( SET1 = v1, SET2 = v2, p = 2 )
                concordance_res <- rbind(
                    concordance_res,
                    data.frame(
                        Sample.ID1 = sinfo[which(sinfo$seq_id == analysis_list[i, "s1"]), "sample_name"], 
                        Sample.ID2 = sinfo[which(sinfo$seq_id == analysis_list[i, "s2"]), "sample_name"],
                        ID1.Variants = length(v1),
                        ID2.Variants = length(v2),
                        Common.Variants = length(intersect(v1, v2)),
                        MatchRate.ID1   = rate1,
                        MatchRate.ID2   = rate2,
                        Concordance.Score = tm
                    )
                )
            }
            SampleGroupTag <- paste0(gsub( "\\-", "_", AnalysisTonly_concordance[k, "SampleGroup"] ), "_")
            concordance_res$Sample.ID1 <- gsub( SampleGroupTag, "", concordance_res$Sample.ID1 )
            concordance_res$Sample.ID2 <- gsub( SampleGroupTag, "", concordance_res$Sample.ID2 )
            #
            ConcordanceResult_Tonly[[ as.character(AnalysisTonly_concordance[k, "SampleGroup"])  ]] <- concordance_res
            
            # venn-diagram
            variant_list <- lapply(indexing(unique(c(analysis_list$s1, analysis_list$s2))), function(idx) 
                data_SomaticVariants %>% filter( variant_call_mode == TonlyCallMode, seq_id == idx ) %>% .$vkey
            )
            names(variant_list) <- sinfo[match(names(variant_list), sinfo$seq_id), "sample_name"]
            names(variant_list) <- gsub(SampleGroupTag, "", names(variant_list))
            
            VarVenn <- WES_analysis.Variant.Concordance.VennDiagram( VariantList=variant_list )
            tmp_png <- sprintf("%s/concordance_venn_Tonly_%s_tmp.png", ResultDir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/concordance_venn_Tonly_%s.png",     ResultDir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) )           
            ggsave( VarVenn, file=tmp_png, width=6, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "50%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "concordance_venn",
                    variant_call    = "Tumor Only",
                    cancer_code     = as.character(AnalysisTonly_concordance[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisTonly_concordance[k, "SampleGroup"]),
                    png             = sprintf("data/concordance_venn_Tonly_%s.png", as.character(AnalysisTonly_concordance[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/concordance_venn_Tonly_%s.png", fig_png, report_data_dir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) ))
            message(sprintf("|----->> Tumor-Only Variants ... concordance venn-diagram created and PNG copied for report : %s", as.character(AnalysisTonly_concordance[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))

            # variant class compare barplot
            VarClassPlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    var_class_compare_data <- data.frame(rbind(
                        data_VarStatSummary %>% filter( variant_call_mode == TonlyCallMode, variant_group == "non_synonymous", seq_id == analysis_list[i, "s1"] ) %>%  
                            group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>% mutate( varN_PCT = varN/sum(varN)*100 , ID=analysis_list[i, "s1"] ),
                        data_VarStatSummary %>% filter( variant_call_mode == TonlyCallMode, variant_group == "non_synonymous", seq_id == analysis_list[i, "s2"] ) %>%  
                            group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>% mutate( varN_PCT = varN/sum(varN)*100 , ID=analysis_list[i, "s2"] )
                    )) 
                    var_class_compare_data$sample_name <- sinfo[match(var_class_compare_data$ID, sinfo$seq_id), "sample_name"]
                    var_class_compare_data$sample_name <- gsub( SampleGroupTag, "", var_class_compare_data$sample_name )
                    #-----------------------------------------------------------
                    vccp <- WES_analysis.Variant.Class.Compare.Barplot(VariantClassCompareData=var_class_compare_data)
                }else{
                    vccp <- ggplot() + theme_void()
                }
                return(vccp)
            })
            VarClassPlots <- cowplot::plot_grid( plotlist=VarClassPlotList, ncol=3 )
            tmp_png <- sprintf("%s/var_class_barplot_Tonly_%s_tmp.png", ResultDir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/var_class_barplot_Tonly_%s.png",     ResultDir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) )           
            ggsave( VarClassPlots, file=tmp_png, width=21, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "variant_class_barplot",
                    variant_call    = "Tumor Only",
                    cancer_code     = as.character(AnalysisTonly_concordance[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisTonly_concordance[k, "SampleGroup"]),
                    png             = sprintf("data/var_class_barplot_Tonly_%s.png", as.character(AnalysisTonly_concordance[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/var_class_barplot_Tonly_%s.png", fig_png, report_data_dir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) ))
            message(sprintf("|----->> Tumor-Only Variants ... variant class barplot created and PNG copied for report : %s", as.character(AnalysisTonly_concordance[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))

            # genomic change compare barplot
            GenomicChangePlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    genomic_change_data <- rbind(
                        data_SomaticVariants %>% filter( seq_id == analysis_list[i, "s1"], Variant_Type == "SNP", variant_call_mode == TonlyCallMode ) %>% mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2)), 
                        data_SomaticVariants %>% filter( seq_id == analysis_list[i, "s2"], Variant_Type == "SNP", variant_call_mode == TonlyCallMode ) %>% mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2))
                    ) %>% 
                        group_by( seq_id, VAC ) %>% summarise( GDC = length(VAC) ) %>%
                        group_by( seq_id ) %>% mutate( GDC_PCT = GDC/sum(GDC)*100 ) %>% 
                        mutate( sample_name = sinfo[match(seq_id, sinfo$seq_id ), "sample_name"] )
                    genomic_change_data$sample_name <- gsub( SampleGroupTag, "", genomic_change_data$sample_name )
                    #-----------------------------------------------------------
                    gcp <- WES_analysis.Genomic.Change.Barplot( GenomicChangeData=genomic_change_data )
                }else{
                    gcp <- ggplot() + theme_void()
                }
                return(gcp)
            })
            GenomicChangePlots <- cowplot::plot_grid( plotlist=GenomicChangePlotList, ncol=3 )
            tmp_png <- sprintf("%s/genomic_change_barplot_Tonly_%s_tmp.png", ResultDir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/genomic_change_barplot_Tonly_%s.png",     ResultDir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) )           
            ggsave( GenomicChangePlots, file=tmp_png, width=21, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "genomic_change_barplot",
                    variant_call    = "Tumor Only",
                    cancer_code     = as.character(AnalysisTonly_concordance[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisTonly_concordance[k, "SampleGroup"]),
                    png             = sprintf("data/genomic_change_barplot_Tonly_%s.png", as.character(AnalysisTonly_concordance[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/genomic_change_barplot_Tonly_%s.png", fig_png, report_data_dir, as.character(AnalysisTonly_concordance[k, "SampleGroup"]) ))
            message(sprintf("|----->> Tumor-Only Variants ... genomic change barplot created and PNG copied for report : %s", as.character(AnalysisTonly_concordance[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))
        }
        message("|----->> Tumor-Only Variants ... variant concordance analyzed.")

        # concordance table save -----------------------------------------------
        ConcordanceResult_Tonly <- ldply(ConcordanceResult_Tonly, .id="Sample.Group")
        saveRDS(ConcordanceResult_Tonly, file=sprintf("%s/concordance_table_Tonly.rds", ResultDir))
        #system(sprintf("cp %s/concordance_table_Tonly.rds %s/concordance_table_Tonly.rds", ResultDir, report_data_dir))
        message("|----->> Tumor-Only Variants ... concordance results RDS saved.")
    }else{
        ConcordanceResult_Tonly <- data.frame()
        message("|----->> Tumor-Only Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched TS /----------------------------------------------------------
    if(  nrow(AnalysisMatchTS) > 0 )
    {
        ConcordanceResult_MatchTS <- list()
        for( k in 1:nrow(AnalysisMatchTS) )
        {
            analysis_list <- AnalysisSetInfoList[[ as.character(AnalysisMatchTS[k, "SampleGroup"]) ]]$match_normal_tissue$analysis_run_list
            # concordacne table
            concordance_res <- data.frame()
            for( i in 1:nrow(analysis_list) )
            {
                v1 <- data_SomaticVariants %>% filter( variant_call_mode == "TS_N", seq_id == analysis_list[i, "s1"] ) %>% .$vkey
                v2 <- data_SomaticVariants %>% filter( variant_call_mode == "TS_N", seq_id == analysis_list[i, "s2"] ) %>% .$vkey
                tm    <- WES_analysis.Get.Tanimoto( SET1 = v1, SET2 = v2 )
                rate1 <- WES_analysis.Get.Common.Rate( SET1 = v1, SET2 = v2, p = 1 )
                rate2 <- WES_analysis.Get.Common.Rate( SET1 = v1, SET2 = v2, p = 2 )
                concordance_res <- rbind(
                    concordance_res,
                    data.frame(
                        Sample.ID1 = sinfo[which(sinfo$seq_id == analysis_list[i, "s1"]), "sample_name"], 
                        Sample.ID2 = sinfo[which(sinfo$seq_id == analysis_list[i, "s2"]), "sample_name"],
                        ID1.Variants = length(v1),
                        ID2.Variants = length(v2),
                        Common.Variants = length(intersect(v1, v2)),
                        MatchRate.ID1   = rate1,
                        MatchRate.ID2   = rate2,
                        Concordance.Score = tm
                    )
                )
            }
            SampleGroupTag <- paste0(gsub( "\\-", "_", AnalysisMatchTS[k, "SampleGroup"] ), "_")
            concordance_res$Sample.ID1 <- gsub( SampleGroupTag, "", concordance_res$Sample.ID1 )
            concordance_res$Sample.ID2 <- gsub( SampleGroupTag, "", concordance_res$Sample.ID2 )
            #
            ConcordanceResult_MatchTS[[ as.character(AnalysisMatchTS[k, "SampleGroup"])  ]] <- concordance_res
            
            # venn-diagram
            variant_list <- lapply(indexing(unique(c(analysis_list$s1, analysis_list$s2))), function(idx) 
                data_SomaticVariants %>% filter( variant_call_mode == "TS_N", seq_id == idx ) %>% .$vkey
            )
            names(variant_list) <- sinfo[match(names(variant_list), sinfo$seq_id), "sample_name"]
            names(variant_list) <- gsub(SampleGroupTag, "", names(variant_list))
            
            VarVenn <- WES_analysis.Variant.Concordance.VennDiagram( VariantList=variant_list )
            tmp_png <- sprintf("%s/concordance_venn_MatchTS_%s_tmp.png", ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/concordance_venn_MatchTS_%s.png",     ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )           
            ggsave( VarVenn, file=tmp_png, width=6, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "50%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "concordance_venn",
                    variant_call    = "Matched Normal Tissue",
                    cancer_code     = as.character(AnalysisMatchTS[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchTS[k, "SampleGroup"]),
                    png             = sprintf("data/concordance_venn_MatchTS_%s.png", as.character(AnalysisMatchTS[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/concordance_venn_MatchTS_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal TISSUE Variants ... concordance venn-diagram created and PNG copied for report : %s", as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))

            # variant class compare barplot
            VarClassPlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    var_class_compare_data <- data.frame(rbind(
                        data_VarStatSummary %>% filter( variant_call_mode == "TS_N", variant_group == "non_synonymous", seq_id == analysis_list[i, "s1"] ) %>%  
                            group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>% mutate( varN_PCT = varN/sum(varN)*100 , ID=analysis_list[i, "s1"] ),
                        data_VarStatSummary %>% filter( variant_call_mode == "TS_N", variant_group == "non_synonymous", seq_id == analysis_list[i, "s2"] ) %>%  
                            group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>% mutate( varN_PCT = varN/sum(varN)*100 , ID=analysis_list[i, "s2"] )
                    )) 
                    var_class_compare_data$sample_name <- sinfo[match(var_class_compare_data$ID, sinfo$seq_id), "sample_name"]
                    var_class_compare_data$sample_name <- gsub( SampleGroupTag, "", var_class_compare_data$sample_name )
                    #-----------------------------------------------------------
                    vccp <- WES_analysis.Variant.Class.Compare.Barplot(VariantClassCompareData=var_class_compare_data)
                }else{
                    vccp <- ggplot() + theme_void()
                }
                return(vccp)
            })
            VarClassPlots <- cowplot::plot_grid( plotlist=VarClassPlotList, ncol=3 )
            tmp_png <- sprintf("%s/var_class_barplot_MatchTS_%s_tmp.png", ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/var_class_barplot_MatchTS_%s.png",     ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )           
            ggsave( VarClassPlots, file=tmp_png, width=21, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "variant_class_barplot",
                    variant_call    = "Matched Normal Tissue",
                    cancer_code     = as.character(AnalysisMatchTS[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchTS[k, "SampleGroup"]),
                    png             = sprintf("data/var_class_barplot_MatchTS_%s.png", as.character(AnalysisMatchTS[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/var_class_barplot_MatchTS_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal TISSUE Variants ... variant class barplot created and PNG copied for report : %s", as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))

            # genomic change compare barplot
            GenomicChangePlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    genomic_change_data <- rbind(
                        data_SomaticVariants %>% filter( seq_id == analysis_list[i, "s1"], Variant_Type == "SNP", variant_call_mode == "TS_N" ) %>% mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2)), 
                        data_SomaticVariants %>% filter( seq_id == analysis_list[i, "s2"], Variant_Type == "SNP", variant_call_mode == "TS_N" ) %>% mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2))
                    ) %>% 
                        group_by( seq_id, VAC ) %>% summarise( GDC = length(VAC) ) %>%
                        group_by( seq_id ) %>% mutate( GDC_PCT = GDC/sum(GDC)*100 ) %>% 
                        mutate( sample_name = sinfo[match(seq_id, sinfo$seq_id ), "sample_name"] )
                    genomic_change_data$sample_name <- gsub( SampleGroupTag, "", genomic_change_data$sample_name )
                    #-----------------------------------------------------------
                    gcp <- WES_analysis.Genomic.Change.Barplot( GenomicChangeData=genomic_change_data )
                }else{
                    gcp <- ggplot() + theme_void()
                }
                return(gcp)
            })
            GenomicChangePlots <- cowplot::plot_grid( plotlist=GenomicChangePlotList, ncol=3 )
            tmp_png <- sprintf("%s/genomic_change_barplot_MatchTS_%s_tmp.png", ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/genomic_change_barplot_MatchTS_%s.png",     ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )           
            ggsave( GenomicChangePlots, file=tmp_png, width=21, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "genomic_change_barplot",
                    variant_call    = "Matched Normal Tissue",
                    cancer_code     = as.character(AnalysisMatchTS[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchTS[k, "SampleGroup"]),
                    png             = sprintf("data/genomic_change_barplot_MatchTS_%s.png", as.character(AnalysisMatchTS[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/genomic_change_barplot_MatchTS_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal TISSUE Variants ... genomic change barplot created and PNG copied for report : %s", as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))
        }
        message("|----->> Matched-Normal TISSUE Variants ... variant concordance analyzed.")

        # concordance table save -----------------------------------------------
        ConcordanceResult_MatchTS <- ldply(ConcordanceResult_MatchTS, .id="Sample.Group")
        saveRDS(ConcordanceResult_MatchTS, file=sprintf("%s/concordance_table_MatchTS.rds", ResultDir))
        #system(sprintf("cp %s/concordance_table_MatchTS.rds %s/concordance_table_MatchTS.rds", ResultDir, report_data_dir))
        message("|----->> Matched-Normal TISSUE Variants ... concordance results RDS saved.")
    }else{
        ConcordanceResult_MatchTS <- data.frame()
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched ORG /---------------------------------------------------------
    if(  nrow(AnalysisMatchORG) > 0 )
    {
        ConcordanceResult_MatchORG <- list()
        for( k in 1:nrow(AnalysisMatchORG) )
        {
            analysis_list <- AnalysisSetInfoList[[ as.character(AnalysisMatchORG[k, "SampleGroup"]) ]]$tumor_only$analysis_run_list
            # concordacne table
            concordance_res <- data.frame()
            for( i in 1:nrow(analysis_list) )
            {
                v1 <- data_SomaticVariants %>% filter( variant_call_mode == "ORG_N", seq_id == analysis_list[i, "s1"] ) %>% .$vkey
                v2 <- data_SomaticVariants %>% filter( variant_call_mode == "ORG_N", seq_id == analysis_list[i, "s2"] ) %>% .$vkey
                tm    <- WES_analysis.Get.Tanimoto( SET1 = v1, SET2 = v2 )
                rate1 <- WES_analysis.Get.Common.Rate( SET1 = v1, SET2 = v2, p = 1 )
                rate2 <- WES_analysis.Get.Common.Rate( SET1 = v1, SET2 = v2, p = 2 )
                concordance_res <- rbind(
                    concordance_res,
                    data.frame(
                        Sample.ID1 = sinfo[which(sinfo$seq_id == analysis_list[i, "s1"]), "sample_name"], 
                        Sample.ID2 = sinfo[which(sinfo$seq_id == analysis_list[i, "s2"]), "sample_name"],
                        ID1.Variants = length(v1),
                        ID2.Variants = length(v2),
                        Common.Variants = length(intersect(v1, v2)),
                        MatchRate.ID1   = rate1,
                        MatchRate.ID2   = rate2,
                        Concordance.Score = tm
                    )
                )
            }
            SampleGroupTag <- paste0(gsub( "\\-", "_", AnalysisMatchORG[k, "SampleGroup"] ), "_")
            concordance_res$Sample.ID1 <- gsub( SampleGroupTag, "", concordance_res$Sample.ID1 )
            concordance_res$Sample.ID2 <- gsub( SampleGroupTag, "", concordance_res$Sample.ID2 )
            #
            ConcordanceResult_MatchORG[[ as.character(AnalysisMatchORG[k, "SampleGroup"])  ]] <- concordance_res
            
            # venn-diagram
            variant_list <- lapply(indexing(unique(c(analysis_list$s1, analysis_list$s2))), function(idx) 
                data_SomaticVariants %>% filter( variant_call_mode == "Tonly", seq_id == idx ) %>% .$vkey
            )
            names(variant_list) <- sinfo[match(names(variant_list), sinfo$seq_id), "sample_name"]
            names(variant_list) <- gsub(SampleGroupTag, "", names(variant_list))
            
            VarVenn <- WES_analysis.Variant.Concordance.VennDiagram( VariantList=variant_list )
            tmp_png <- sprintf("%s/concordance_venn_MatchORG_%s_tmp.png", ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/concordance_venn_MatchORG_%s.png",     ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )           
            ggsave( VarVenn, file=tmp_png, width=6, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "50%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "concordance_venn",
                    variant_call    = "Matched Normal Organoid",
                    cancer_code     = as.character(AnalysisMatchORG[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchORG[k, "SampleGroup"]),
                    png             = sprintf("data/concordance_venn_MatchORG_%s.png", as.character(AnalysisMatchORG[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/concordance_venn_MatchORG_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchORG[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal ORGANOID Variants ... concordance venn-diagram created and PNG copied for report : %s", as.character(AnalysisMatchORG[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))

            # variant class compare barplot
            VarClassPlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    var_class_compare_data <- data.frame(rbind(
                        data_VarStatSummary %>% filter( variant_call_mode == "ORG_N", variant_group == "non_synonymous", seq_id == analysis_list[i, "s1"] ) %>%  
                            group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>% mutate( varN_PCT = varN/sum(varN)*100 , ID=analysis_list[i, "s1"] ),
                        data_VarStatSummary %>% filter( variant_call_mode == "ORG_N", variant_group == "non_synonymous", seq_id == analysis_list[i, "s2"] ) %>%  
                            group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>% mutate( varN_PCT = varN/sum(varN)*100 , ID=analysis_list[i, "s2"] )
                    )) 
                    var_class_compare_data$sample_name <- sinfo[match(var_class_compare_data$ID, sinfo$seq_id), "sample_name"]
                    var_class_compare_data$sample_name <- gsub( SampleGroupTag, "", var_class_compare_data$sample_name )
                    #-----------------------------------------------------------
                    vccp <- WES_analysis.Variant.Class.Compare.Barplot(VariantClassCompareData=var_class_compare_data)
                }else{
                    vccp <- ggplot() + theme_void()
                }
                return(vccp)
            })
            VarClassPlots <- cowplot::plot_grid( plotlist=VarClassPlotList, ncol=3 )
            tmp_png <- sprintf("%s/var_class_barplot_MatchORG_%s_tmp.png", ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/var_class_barplot_MatchORG_%s.png",     ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )           
            ggsave( VarClassPlots, file=tmp_png, width=21, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "variant_class_barplot",
                    variant_call    = "Matched Normal Organoid",
                    cancer_code     = as.character(AnalysisMatchORG[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchORG[k, "SampleGroup"]),
                    png             = sprintf("data/var_class_barplot_MatchORG_%s.png", as.character(AnalysisMatchORG[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/var_class_barplot_MatchORG_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchORG[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal ORGANOID Variants ... variant class barplot created and PNG copied for report : %s", as.character(AnalysisMatchORG[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))

            # genomic change compare barplot
            GenomicChangePlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    genomic_change_data <- rbind(
                        data_SomaticVariants %>% filter( seq_id == analysis_list[i, "s1"], Variant_Type == "SNP" ) %>% mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2)), 
                        data_SomaticVariants %>% filter( seq_id == analysis_list[i, "s2"], Variant_Type == "SNP" ) %>% mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2))
                    ) %>% 
                        group_by( seq_id, VAC ) %>% summarise( GDC = length(VAC) ) %>%
                        group_by( seq_id ) %>% mutate( GDC_PCT = GDC/sum(GDC)*100 ) %>% 
                        mutate( sample_name = sinfo[match(seq_id, sinfo$seq_id ), "sample_name"] )
                    genomic_change_data$sample_name <- gsub( SampleGroupTag, "", genomic_change_data$sample_name )
                    #-----------------------------------------------------------
                    gcp <- WES_analysis.Genomic.Change.Barplot( GenomicChangeData=genomic_change_data )
                }else{
                    gcp <- ggplot() + theme_void()
                }
                return(gcp)
            })
            GenomicChangePlots <- cowplot::plot_grid( plotlist=GenomicChangePlotList, ncol=3 )
            tmp_png <- sprintf("%s/genomic_change_barplot_MatchORG_%s_tmp.png", ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/genomic_change_barplot_MatchORG_%s.png",     ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )           
            ggsave( GenomicChangePlots, file=tmp_png, width=21, height=6, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "genomic_change_barplot",
                    variant_call    = "Matched Normal Organoid",
                    cancer_code     = as.character(AnalysisMatchORG[k, "cancer_code"]),
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchORG[k, "SampleGroup"]),
                    png             = sprintf("data/genomic_change_barplot_MatchORG_%s.png", as.character(AnalysisMatchORG[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/genomic_change_barplot_MatchORG_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchORG[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal ORGANOID Variants ... genomic change barplot created and PNG copied for report : %s", as.character(AnalysisMatchORG[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))
        }
        message("|----->> Matched-Normal ORGANOID Variants ... variant concordance analyzed.")

        # concordance table save -----------------------------------------------
        ConcordanceResult_MatchORG <- ldply(ConcordanceResult_MatchORG, .id="Sample.Group")
        saveRDS(ConcordanceResult_MatchORG, file=sprintf("%s/concordance_table_MatchORG.rds", ResultDir))
        #system(sprintf("cp %s/concordance_table_MatchORG.rds %s/concordance_table_MatchORG.rds", ResultDir, report_data_dir))
        message("|----->> Matched-Normal ORGANOID Variants ... concordance results RDS saved.")
    }else{
        ConcordanceResult_MatchORG <- data.frame()
        message("|----->> Matched-Normal ORGANOID Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    ConcordanceResultList <- list(
        Tonly    = ConcordanceResult_Tonly,
        MatchTS  = ConcordanceResult_MatchTS,
        MatchORG = ConcordanceResult_MatchORG
    )
    ConcordanceResultList <- lapply(ConcordanceResultList, function(y) y[which(y$Sample.Group != ""), ] )
    #---------------------------------------------------------------------------
    saveRDS(ConcordanceResultList, file=sprintf("%s/concordance_table_list.rds", ResultDir))
    system(sprintf("cp %s/concordance_table_list.rds %s/concordance_table_list.rds", ResultDir, report_data_dir))
    message("|----->> All Concordance results RDS saved and copied for report.")

#---------------------------------------------------------------------------------------------------

#---/ ANALYSIS : VAF CORRELATION /------------------------------------------------------------------
    message( sprintf('\n|--->> Run "VAF CORRELATION" analysis <<%s| \n', paste(rep("-", 100), collapse="") ))

    #---/ TUMOR ONLY /----------------------------------------------------------
    AnalysisTonly_vaf <- AnalysisTonly %>% filter( SampleGroup != "" )
    if(  nrow(AnalysisTonly_vaf) > 0 )
    {
        for( k in 1:nrow(AnalysisTonly_vaf) )
        {
            analysis_list <- AnalysisSetInfoList[[ as.character(AnalysisTonly_vaf[k, "SampleGroup"]) ]]$tumor_only$analysis_run_list
            # raw all variants 
            union(analysis_list$s1, analysis_list$s2)
            raw_variants_list <- lapply( indexing(union(analysis_list$s1, analysis_list$s2)), function(idx)
            {
                read.delim(sprintf("%s/%s/%s/vcf/%s.mutect2.bias.filtered.vep.annotated.maf", BaseDir, SeqFolderID, idx, idx)) %>%
                    mutate(VAF = round(TUMOR_ALT_COUNT/DEPTH, 3) ) %>% dplyr::select(c("vkey","seq_id","VAF"))
            })
            VafCorrPlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    VF <- merge( x=raw_variants_list[[ analysis_list[i, "s1"] ]], y=raw_variants_list[[ analysis_list[i, "s2"] ]], all.x=TRUE, all.y=TRUE, by="vkey") %>%
                        dplyr::rename(s1=2,s1_vaf=3,s2=4,s2_vaf=5)

                    v1 <- data_SomaticVariants %>% filter( variant_call_mode == TonlyCallMode, seq_id == analysis_list[i, "s1"] ) %>% .$vkey
                    v2 <- data_SomaticVariants %>% filter( variant_call_mode == TonlyCallMode, seq_id == analysis_list[i, "s2"] ) %>% .$vkey
                    VF[which(VF$vkey %in% v1), "s1_somatic"] = 1
                    VF[which(VF$vkey %in% v2), "s2_somatic"] = 1
                    
                    SampleGroupTag <- paste0(gsub( "\\-", "_", AnalysisTonly_vaf[k, "SampleGroup"] ), "_")
                    VF$s1          <- gsub( SampleGroupTag, "", sinfo[match(VF$s1, sinfo$seq_id), "sample_name"])
                    VF$s2          <- gsub( SampleGroupTag, "", sinfo[match(VF$s2, sinfo$seq_id), "sample_name"])
                    

                    if( unique(VF$s1[!is.na(VF$s1)]) == "TS"  ){
                            VF <- VF %>% dplyr::rename(X=s1,X.vaf=s1_vaf,X.mut=s1_somatic, Y=s2,Y.vaf=s2_vaf,Y.mut=s2_somatic)
                    }else if( unique(VF$s2[!is.na(VF$s2)]) == "TS"  ){
                            VF <- VF %>% dplyr::rename(Y=s1,Y.vaf=s1_vaf,Y.mut=s1_somatic, X=s2,X.vaf=s2_vaf,X.mut=s2_somatic)                   
                    }else{  VF <- VF %>% dplyr::rename(X=s1,X.vaf=s1_vaf,X.mut=s1_somatic, Y=s2,Y.vaf=s2_vaf,Y.mut=s2_somatic) }

                    VafCorrPlot <- WES_analysis.Vaf.Correlation.Plot( VafCorData=VF )
                }else{
                    VafCorrPlot <- ggplot() + theme_void()
                }
                return(VafCorrPlot)
            })
            VafCorrPlots <- cowplot::plot_grid( plotlist=VafCorrPlotList, ncol=3 )
            tmp_png <- sprintf("%s/vaf_corr_plot_Tonly_%s_tmp.png", ResultDir, as.character(AnalysisTonly_vaf[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/vaf_corr_plot_Tonly_%s.png",     ResultDir, as.character(AnalysisTonly_vaf[k, "SampleGroup"]) )           
            ggsave( VafCorrPlots, file=tmp_png, width=21, height=6.5, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "vaf_correlation",
                    variant_call    = "Tumor Only",
                    cancer_code     = as.character(AnalysisTonly_vaf[k, "cancer_code"]) ,
                    preset_name     = "",
                    sample_group    = as.character(AnalysisTonly_vaf[k, "SampleGroup"]),
                    png             = sprintf("data/vaf_corr_plot_Tonly_%s.png", as.character(AnalysisTonly_vaf[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/vaf_corr_plot_Tonly_%s.png", fig_png, report_data_dir, as.character(AnalysisTonly_vaf[k, "SampleGroup"]) ))
            message(sprintf("|----->> Tumor-Only Variants ... scatterplot created and PNG copied for report : %s", as.character(AnalysisTonly_vaf[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))
        }
    }else{
        message("|----->> Tumor-Only Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched TS /----------------------------------------------------------
    if(  nrow(AnalysisMatchTS) > 0 )
    {
        for( k in 1:nrow(AnalysisMatchTS) )
        {
            analysis_list <- AnalysisSetInfoList[[ as.character(AnalysisMatchTS[k, "SampleGroup"]) ]]$match_normal_tissue$analysis_run_list
            # raw all variants 
            union(analysis_list$s1, analysis_list$s2)
            raw_variants_list <- lapply( indexing(union(analysis_list$s1, analysis_list$s2)), function(idx)
            {
                read.delim(sprintf("%s/%s/%s/vcf/%s.mutect2.NT.bias.filtered.vep.annotated.maf", BaseDir, SeqFolderID, idx, idx)) %>%
                    mutate(VAF = round(TUMOR_ALT_COUNT/DEPTH, 3) ) %>% dplyr::select(c("vkey","seq_id","VAF"))
            })

            VafCorrPlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    VF <- merge( x=raw_variants_list[[ analysis_list[i, "s1"] ]], y=raw_variants_list[[ analysis_list[i, "s2"] ]], all.x=TRUE, all.y=TRUE, by="vkey") %>%
                        dplyr::rename(s1=2,s1_vaf=3,s2=4,s2_vaf=5)

                    v1 <- data_SomaticVariants %>% filter( variant_call_mode == "TS_N", seq_id == analysis_list[i, "s1"] ) %>% .$vkey
                    v2 <- data_SomaticVariants %>% filter( variant_call_mode == "TS_N", seq_id == analysis_list[i, "s2"] ) %>% .$vkey
                    VF[which(VF$vkey %in% v1), "s1_somatic"] = 1
                    VF[which(VF$vkey %in% v2), "s2_somatic"] = 1
                    
                    SampleGroupTag <- paste0(gsub( "\\-", "_", AnalysisMatchTS[k, "SampleGroup"] ), "_")
                    VF$s1          <- gsub( SampleGroupTag, "", sinfo[match(VF$s1, sinfo$seq_id), "sample_name"])
                    VF$s2          <- gsub( SampleGroupTag, "", sinfo[match(VF$s2, sinfo$seq_id), "sample_name"])
                    

                    if( unique(VF$s1[!is.na(VF$s1)]) == "TS"  ){
                            VF <- VF %>% dplyr::rename(X=s1,X.vaf=s1_vaf,X.mut=s1_somatic, Y=s2,Y.vaf=s2_vaf,Y.mut=s2_somatic)
                    }else if( unique(VF$s2[!is.na(VF$s2)]) == "TS"  ){
                            VF <- VF %>% dplyr::rename(Y=s1,Y.vaf=s1_vaf,Y.mut=s1_somatic, X=s2,X.vaf=s2_vaf,X.mut=s2_somatic)                   
                    }else{  VF <- VF %>% dplyr::rename(X=s1,X.vaf=s1_vaf,X.mut=s1_somatic, Y=s2,Y.vaf=s2_vaf,Y.mut=s2_somatic) }

                    VafCorrPlot <- WES_analysis.Vaf.Correlation.Plot( VafCorData=VF )
                }else{
                    VafCorrPlot <- ggplot() + theme_void()
                }
                return(VafCorrPlot)
            })
            VafCorrPlots <- cowplot::plot_grid( plotlist=VafCorrPlotList, ncol=3 )
            tmp_png <- sprintf("%s/vaf_corr_plot_MatchTS_%s_tmp.png", ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/vaf_corr_plot_MatchTS_%s.png",     ResultDir, as.character(AnalysisMatchTS[k, "SampleGroup"]) )           
            ggsave( VafCorrPlots, file=tmp_png, width=21, height=6.5, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "vaf_correlation",
                    variant_call    = "Matched Normal Tissue",
                    cancer_code     = as.character(AnalysisMatchTS[k, "cancer_code"]) ,
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchTS[k, "SampleGroup"]),
                    png             = sprintf("data/vaf_corr_plot_MatchTS_%s.png", as.character(AnalysisMatchTS[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/vaf_corr_plot_MatchTS_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal TISSUE Variants ... scatterplot created and PNG copied for report : %s", as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))
        }
    }else{
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched ORG /---------------------------------------------------------
    if(  nrow(AnalysisMatchORG) > 0 )
    {
        for( k in 1:nrow(AnalysisMatchORG) )
        {
            analysis_list <- AnalysisSetInfoList[[ as.character(AnalysisMatchORG[k, "SampleGroup"]) ]]$match_normal_organoid$analysis_run_list
            # raw all variants 
            union(analysis_list$s1, analysis_list$s2)
            raw_variants_list <- lapply( indexing(union(analysis_list$s1, analysis_list$s2)), function(idx)
            {
                read.delim(sprintf("%s/%s/%s/vcf/%s.mutect2.ORG.NT.bias.filtered.vep.annotated.maf", BaseDir, SeqFolderID, idx, idx)) %>%
                    mutate(VAF = round(TUMOR_ALT_COUNT/DEPTH, 3) ) %>% dplyr::select(c("vkey","seq_id","VAF"))
            })

            VafCorrPlotList <- lapply(1:3, function(i)
            {
                if( i <= nrow(analysis_list) )
                {
                    VF <- merge( x=raw_variants_list[[ analysis_list[i, "s1"] ]], y=raw_variants_list[[ analysis_list[i, "s2"] ]], all.x=TRUE, all.y=TRUE, by="vkey") %>%
                        dplyr::rename(s1=2,s1_vaf=3,s2=4,s2_vaf=5)

                    v1 <- data_SomaticVariants %>% filter( variant_call_mode == "ORG_N", seq_id == analysis_list[i, "s1"] ) %>% .$vkey
                    v2 <- data_SomaticVariants %>% filter( variant_call_mode == "ORG_N", seq_id == analysis_list[i, "s2"] ) %>% .$vkey
                    VF[which(VF$vkey %in% v1), "s1_somatic"] = 1
                    VF[which(VF$vkey %in% v2), "s2_somatic"] = 1
                    
                    SampleGroupTag <- paste0(gsub( "\\-", "_", AnalysisMatchORG[k, "SampleGroup"] ), "_")
                    VF$s1          <- gsub( SampleGroupTag, "", sinfo[match(VF$s1, sinfo$seq_id), "sample_name"])
                    VF$s2          <- gsub( SampleGroupTag, "", sinfo[match(VF$s2, sinfo$seq_id), "sample_name"])
                    

                    if( unique(VF$s1[!is.na(VF$s1)]) == "TS"  ){
                            VF <- VF %>% dplyr::rename(X=s1,X.vaf=s1_vaf,X.mut=s1_somatic, Y=s2,Y.vaf=s2_vaf,Y.mut=s2_somatic)
                    }else if( unique(VF$s2[!is.na(VF$s2)]) == "TS"  ){
                            VF <- VF %>% dplyr::rename(Y=s1,Y.vaf=s1_vaf,Y.mut=s1_somatic, X=s2,X.vaf=s2_vaf,X.mut=s2_somatic)                   
                    }else{  VF <- VF %>% dplyr::rename(X=s1,X.vaf=s1_vaf,X.mut=s1_somatic, Y=s2,Y.vaf=s2_vaf,Y.mut=s2_somatic) }

                    VafCorrPlot <- WES_analysis.Vaf.Correlation.Plot( VafCorData=VF )
                }else{
                    VafCorrPlot <- ggplot() + theme_void()
                }
                return(VafCorrPlot)
            })
            VafCorrPlots <- cowplot::plot_grid( plotlist=VafCorrPlotList, ncol=3 )
            tmp_png <- sprintf("%s/vaf_corr_plot_MatchORG_%s_tmp.png", ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )
            fig_png <- sprintf("%s/vaf_corr_plot_MatchORG_%s.png",     ResultDir, as.character(AnalysisMatchORG[k, "SampleGroup"]) )           
            ggsave( VafCorrPlots, file=tmp_png, width=21, height=6.5, unit='in', dpi=200 )
            #system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
            system(sprintf("convert %s -resize %s %s", tmp_png, "30%", fig_png))
            # fig list 
            FIGUER_LIST <- rbind(
                FIGUER_LIST,
                data.frame( 
                    report_category = "vaf_correlation",
                    variant_call    = "Matched Normal Organoid",
                    cancer_code     = as.character(AnalysisMatchTS[k, "cancer_code"]) ,
                    preset_name     = "",
                    sample_group    = as.character(AnalysisMatchTS[k, "SampleGroup"]),
                    png             = sprintf("data/vaf_corr_plot_MatchORG_%s.png", as.character(AnalysisMatchTS[k, "SampleGroup"]))
                )
            )
            system(sprintf("cp %s %s/vaf_corr_plot_MatchORG_%s.png", fig_png, report_data_dir, as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            message(sprintf("|----->> Matched-Normal ORGANOID Variants ... scatterplot created and PNG copied for report : %s", as.character(AnalysisMatchTS[k, "SampleGroup"]) ))
            system(sprintf("rm %s", tmp_png))
        }
    }else{
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------

#---/ ANALYSIS : ONCOPLOTS 1 /----------------------------------------------------------------------
    message( sprintf('\n|--->> Run "ONCOPLOT-1" analysis <<%s| \n', paste(rep("-", 100), collapse="") ))

    #---/ TUMOR ONLY /----------------------------------------------------------
    if(  nrow(sinfo_preset_vars_tonly) > 0 )
    {
        HeatmapTargetCancerCodeTonly <- sinfo_preset_vars_tonly %>% group_by(sample_info) %>% reframe( sample_counts = length(seq_id) ) %>% as.data.frame()
        HeatmapTargetCancerCodeTonly <- HeatmapTargetCancerCodeTonly %>% filter( sample_counts >= OncoplotSampleLimit )
        #-----------------------------------------------------------------------
        if( nrow(HeatmapTargetCancerCodeTonly) > 0 )
        {
            for( k in 1:nrow(HeatmapTargetCancerCodeTonly) )
            {
                CancerCodeHeatmapSampleInfo <- sinfo %>% filter( sample_info == HeatmapTargetCancerCodeTonly[k, "sample_info"] ) 
                # top30 gene oncoplot ------------------------------------------
                oncoplot_top30 <- WES_analysis.Draw.OncoPlots(
                    targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                selectGroupName  = "all", 
                                                selectCancerType = HeatmapTargetCancerCodeTonly[k, "sample_info"] , 
                                                VariantData      = data_SomaticVariants %>% filter( variant_call_mode == TonlyCallMode ), 
                                                topGenes         = 30, 
                                                geneSets         = NULL, 
                                                useBoth          = FALSE, 
                                                SampleInfo       = CancerCodeHeatmapSampleInfo
                                            ),
                    SampleInfo           = CancerCodeHeatmapSampleInfo,
                    SampleOrderByIdGroup = TRUE,
                    includeAllSamples    = TRUE
                )
                
                tmp_png <- sprintf("%s/oncoplot_top30_Tonly_%s_tmp.png", ResultDir, HeatmapTargetCancerCodeTonly[k, "sample_info"] ) 
                fig_png <- sprintf("%s/oncoplot_top30_Tonly_%s.png",     ResultDir, HeatmapTargetCancerCodeTonly[k, "sample_info"] ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(oncoplot_top30)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                # image size
                hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                if( hm_size[1] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                if( hm_size[2] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[2] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "oncoplot_1",
                        variant_call    = "Tumor Only",
                        cancer_code     = HeatmapTargetCancerCodeTonly[k, "sample_info"],
                        preset_name     = "top30",
                        sample_group    = "",
                        png             = sprintf("data/oncoplot_top30_Tonly_%s.png", HeatmapTargetCancerCodeTonly[k, "sample_info"])
                    )
                )
                system(sprintf("cp %s %s/oncoplot_top30_Tonly_%s.png", fig_png, report_data_dir, HeatmapTargetCancerCodeTonly[k, "sample_info"]) )
                message(sprintf("|----->> Tumor-Only Variants ... oncoplot(1) created and PNG copied for report : %s", HeatmapTargetCancerCodeTonly[k, "sample_info"]) )
                system(sprintf("rm %s", tmp_png))

                # preset gene oncoplot -----------------------------------------
                CancerCodePresetNameList <- PresetList[which(PresetList$cancer_code == HeatmapTargetCancerCodeTonly[k, "sample_info"]), "preset_name"]
                for( PSN in CancerCodePresetNameList )
                {
                    oncoplot_preset <- WES_analysis.Draw.OncoPlots(
                        targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                    selectGroupName  = "all", 
                                                    selectCancerType = HeatmapTargetCancerCodeTonly[k, "sample_info"] , 
                                                    VariantData      = data_SomaticVariants %>% filter( variant_call_mode == TonlyCallMode ), 
                                                    topGenes         = NULL, 
                                                    geneSets         = ClientPresetGenes[[ PSN ]], 
                                                    useBoth          = FALSE, 
                                                    SampleInfo       = CancerCodeHeatmapSampleInfo
                                                ),
                        SampleInfo           = CancerCodeHeatmapSampleInfo,
                        SampleOrderByIdGroup = TRUE,
                        includeAllSamples    = TRUE
                    )

                    tmp_png <- sprintf("%s/oncoplot_preset_%s_Tonly_%s_tmp.png", ResultDir, PSN, HeatmapTargetCancerCodeTonly[k, "sample_info"] ) 
                    fig_png <- sprintf("%s/oncoplot_preset_%s_Tonly_%s.png",     ResultDir, PSN, HeatmapTargetCancerCodeTonly[k, "sample_info"] ) 
                    png( tmp_png, width=20, height=20, units="in", res=200)
                    draw(oncoplot_preset)
                    dev.off()
                    system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                    # image size
                    hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                    hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                    if( hm_size[1] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    if( hm_size[2] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[2] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    # fig list 
                    FIGUER_LIST <- rbind(
                        FIGUER_LIST,
                        data.frame( 
                            report_category = "oncoplot_1",
                            variant_call    = "Tumor Only",
                            cancer_code     = HeatmapTargetCancerCodeTonly[k, "sample_info"],
                            preset_name     = PSN,
                            sample_group    = "",
                            png             = sprintf("data/oncoplot_preset_%s_Tonly_%s.png", PSN, HeatmapTargetCancerCodeTonly[k, "sample_info"])
                        )
                    )
                    system(sprintf("cp %s %s/oncoplot_preset_%s_Tonly_%s.png", fig_png, report_data_dir, PSN, HeatmapTargetCancerCodeTonly[k, "sample_info"]) )
                    message(sprintf("|----->> Tumor-Only Variants ... oncoplot(1) created and PNG copied for report : %s", HeatmapTargetCancerCodeTonly[k, "sample_info"]) )
                    system(sprintf("rm %s", tmp_png))
                    
                }
            }
        }
    }else{
        message("|----->> Tumor-Only Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched TS /----------------------------------------------------------
    if(  nrow(AnalysisMatchTS) > 0 )
    {
        HeatmapTargetCancerCodeMatchTS <- sinfo %>% filter( sample_group %in% AnalysisMatchTS[which(AnalysisMatchTS$Match_TS), "SampleGroup"], matched_normal %in% c(0,1) ) %>% 
            group_by(sample_info) %>% reframe( sample_counts = length(seq_id) ) %>% as.data.frame()
        HeatmapTargetCancerCodeMatchTS <- HeatmapTargetCancerCodeMatchTS %>% filter( sample_counts >= OncoplotSampleLimit )
        #-----------------------------------------------------------------------
        if( nrow(HeatmapTargetCancerCodeMatchTS) > 0 )
        {
            for( k in 1:nrow(HeatmapTargetCancerCodeMatchTS) )
            {
                CancerCodeHeatmapSampleInfo <- sinfo %>% filter( sample_info == HeatmapTargetCancerCodeMatchTS[k, "sample_info"], matched_normal %in% c(0,1) ) 
                # top30 gene oncoplot ----------------------------------------------
                oncoplot_top30 <- WES_analysis.Draw.OncoPlots(
                    targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                selectGroupName  = "all", 
                                                selectCancerType = HeatmapTargetCancerCodeMatchTS[k, "sample_info"] , 
                                                VariantData      = data_SomaticVariants %>% filter( variant_call_mode == "TS_N" ), 
                                                topGenes         = 30, 
                                                geneSets         = NULL, 
                                                useBoth          = FALSE, 
                                                SampleInfo       = CancerCodeHeatmapSampleInfo
                                            ),
                    SampleInfo           = CancerCodeHeatmapSampleInfo,
                    SampleOrderByIdGroup = TRUE,
                    includeAllSamples    = TRUE
                )
                
                tmp_png <- sprintf("%s/oncoplot_top30_MatchTS_%s_tmp.png", ResultDir, HeatmapTargetCancerCodeMatchTS[k, "sample_info"] ) 
                fig_png <- sprintf("%s/oncoplot_top30_MatchTS_%s.png",     ResultDir, HeatmapTargetCancerCodeMatchTS[k, "sample_info"] ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(oncoplot_top30)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                # image size
                hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                if( hm_size[1] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                if( hm_size[2] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[2] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "oncoplot_1",
                        variant_call    = "Matched Normal Tissue",
                        cancer_code     = HeatmapTargetCancerCodeMatchTS[k, "sample_info"],
                        preset_name     = "top30",
                        sample_group    = "",
                        png             = sprintf("data/oncoplot_top30_MatchTS_%s.png", HeatmapTargetCancerCodeMatchTS[k, "sample_info"])
                    )
                )
                system(sprintf("cp %s %s/oncoplot_top30_MatchTS_%s.png", fig_png, report_data_dir, HeatmapTargetCancerCodeMatchTS[k, "sample_info"]) )
                message(sprintf("|----->> Matched-Normal TISSUE Variants ... oncoplot(1) created and PNG copied for report : %s", HeatmapTargetCancerCodeMatchTS[k, "sample_info"]) )
                system(sprintf("rm %s", tmp_png))

                # preset gene oncoplot ---------------------------------------------
                CancerCodePresetNameList <- PresetList[which(PresetList$cancer_code == HeatmapTargetCancerCodeMatchTS[k, "sample_info"]), "preset_name"]
                for( PSN in CancerCodePresetNameList )
                {
                    oncoplot_preset <- WES_analysis.Draw.OncoPlots(
                        targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                    selectGroupName  = "all", 
                                                    selectCancerType = HeatmapTargetCancerCodeMatchTS[k, "sample_info"] , 
                                                    VariantData      = data_SomaticVariants %>% filter( variant_call_mode == "TS_N" ), 
                                                    topGenes         = NULL, 
                                                    geneSets         = ClientPresetGenes[[ PSN ]], 
                                                    useBoth          = FALSE, 
                                                    SampleInfo       = CancerCodeHeatmapSampleInfo
                                                ),
                        SampleInfo           = CancerCodeHeatmapSampleInfo,
                        SampleOrderByIdGroup = TRUE,
                        includeAllSamples    = TRUE
                    )

                    tmp_png <- sprintf("%s/oncoplot_preset_%s_MatchTS_%s_tmp.png", ResultDir, PSN, HeatmapTargetCancerCodeMatchTS[k, "sample_info"] ) 
                    fig_png <- sprintf("%s/oncoplot_preset_%s_MatchTS_%s.png",     ResultDir, PSN, HeatmapTargetCancerCodeMatchTS[k, "sample_info"] ) 
                    png( tmp_png, width=20, height=20, units="in", res=200)
                    draw(oncoplot_preset)
                    dev.off()
                    system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                    # image size
                    hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                    hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                    if( hm_size[1] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    if( hm_size[2] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[2] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    # fig list 
                    FIGUER_LIST <- rbind(
                        FIGUER_LIST,
                        data.frame( 
                            report_category = "oncoplot_1",
                            variant_call    = "Matched Normal Tissue",
                            cancer_code     = HeatmapTargetCancerCodeMatchTS[k, "sample_info"],
                            preset_name     = PSN,
                            sample_group    = "",
                            png             = sprintf("data/oncoplot_preset_%s_MatchTS_%s.png", PSN, HeatmapTargetCancerCodeMatchTS[k, "sample_info"])
                        )
                    )
                    system(sprintf("cp %s %s/oncoplot_preset_%s_MatchTS_%s.png", fig_png, report_data_dir, PSN, HeatmapTargetCancerCodeMatchTS[k, "sample_info"]) )
                    message(sprintf("|----->> Matched-Normal TISSUE Variants ... oncoplot(1) created and PNG copied for report : %s", HeatmapTargetCancerCodeMatchTS[k, "sample_info"]) )
                    system(sprintf("rm %s", tmp_png))
                    
                }
            }
        }
    }else{
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ Matched ORG /---------------------------------------------------------
    if(  nrow(AnalysisMatchORG) > 0 )
    {
        HeatmapTargetCancerCodeMatchORG <- sinfo %>% filter( sample_group %in% AnalysisMatchORG[which(AnalysisMatchORG$Match_ORG), "SampleGroup"], matched_normal %in% c(0,1) ) %>% 
            group_by(sample_info) %>% reframe( sample_counts = length(seq_id) ) %>% as.data.frame()
        HeatmapTargetCancerCodeMatchORG <- HeatmapTargetCancerCodeMatchORG %>% filter( sample_counts >= OncoplotSampleLimit )
        #-----------------------------------------------------------------------
        if( nrow(HeatmapTargetCancerCodeMatchORG) > 0 )
        {
            for( k in 1:nrow(HeatmapTargetCancerCodeMatchORG) )
            {
                CancerCodeHeatmapSampleInfo <- sinfo %>% filter( sample_info == HeatmapTargetCancerCodeMatchORG[k, "sample_info"], matched_normal %in% c(0,1) ) 
                # top30 gene oncoplot ------------------------------------------
                oncoplot_top30 <- WES_analysis.Draw.OncoPlots(
                    targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                selectGroupName  = "all", 
                                                selectCancerType = HeatmapTargetCancerCodeMatchORG[k, "sample_info"] , 
                                                VariantData      = data_SomaticVariants %>% filter( variant_call_mode == "ORG_N" ), 
                                                topGenes         = 30, 
                                                geneSets         = NULL, 
                                                useBoth          = FALSE, 
                                                SampleInfo       = CancerCodeHeatmapSampleInfo
                                            ),
                    SampleInfo           = CancerCodeHeatmapSampleInfo,
                    SampleOrderByIdGroup = TRUE,
                    includeAllSamples    = TRUE
                )
                
                tmp_png <- sprintf("%s/oncoplot_top30_MatchORG_%s_tmp.png", ResultDir, HeatmapTargetCancerCodeMatchORG[k, "sample_info"] ) 
                fig_png <- sprintf("%s/oncoplot_top30_MatchORG_%s.png",     ResultDir, HeatmapTargetCancerCodeMatchORG[k, "sample_info"] ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(oncoplot_top30)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                # image size
                hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                if( hm_size[1] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                if( hm_size[2] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[2] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "oncoplot_1",
                        variant_call    = "Matchced Normal Organoid",
                        cancer_code     = HeatmapTargetCancerCodeMatchORG[k, "sample_info"],
                        preset_name     = "top30",
                        sample_group    = "",
                        png             = sprintf("data/oncoplot_top30_MatchORG_%s.png", HeatmapTargetCancerCodeMatchORG[k, "sample_info"])
                    )
                )
                system(sprintf("cp %s %s/oncoplot_top30_MatchORG_%s.png", fig_png, report_data_dir, HeatmapTargetCancerCodeMatchORG[k, "sample_info"]) )
                message(sprintf("|----->> Matched-Normal ORGANOID Variants ... oncoplot(1) created and PNG copied for report : %s", HeatmapTargetCancerCodeMatchORG[k, "sample_info"]) )
                system(sprintf("rm %s", tmp_png))

                # preset gene oncoplot -----------------------------------------
                CancerCodePresetNameList <- PresetList[which(PresetList$cancer_code == HeatmapTargetCancerCodeMatchORG[k, "sample_info"]), "preset_name"]
                for( PSN in CancerCodePresetNameList )
                {
                    oncoplot_preset <- WES_analysis.Draw.OncoPlots(
                        targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                    selectGroupName  = "all", 
                                                    selectCancerType = HeatmapTargetCancerCodeMatchORG[k, "sample_info"] , 
                                                    VariantData      = data_SomaticVariants %>% filter( variant_call_mode == "ORG_N" ), 
                                                    topGenes         = NULL, 
                                                    geneSets         = ClientPresetGenes[[ PSN ]], 
                                                    useBoth          = FALSE, 
                                                    SampleInfo       = CancerCodeHeatmapSampleInfo
                                                ),
                        SampleInfo           = CancerCodeHeatmapSampleInfo,
                        SampleOrderByIdGroup = TRUE,
                        includeAllSamples    = TRUE
                    )

                    tmp_png <- sprintf("%s/oncoplot_preset_%s_MatchORG_%s_tmp.png", ResultDir, PSN, HeatmapTargetCancerCodeMatchORG[k, "sample_info"] ) 
                    fig_png <- sprintf("%s/oncoplot_preset_%s_MatchORG_%s.png",     ResultDir, PSN, HeatmapTargetCancerCodeMatchORG[k, "sample_info"] ) 
                    png( tmp_png, width=20, height=20, units="in", res=200)
                    draw(oncoplot_preset)
                    dev.off()
                    system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                    # image size
                    hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                    hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                    if( hm_size[1] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    if( hm_size[2] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[2] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    # fig list 
                    FIGUER_LIST <- rbind(
                        FIGUER_LIST,
                        data.frame( 
                            report_category = "oncoplot_1",
                            variant_call    = "Matched Normal Organoid",
                            cancer_code     = HeatmapTargetCancerCodeMatchORG[k, "sample_info"],
                            preset_name     = PSN,
                            sample_group    = "",
                            png             = sprintf("data/oncoplot_preset_%s_MatchORG_%s.png", PSN, HeatmapTargetCancerCodeMatchORG[k, "sample_info"])
                        )
                    )
                    system(sprintf("cp %s %s/oncoplot_preset_%s_MatchORG_%s.png", fig_png, report_data_dir, PSN, HeatmapTargetCancerCodeMatchORG[k, "sample_info"]) )
                    message(sprintf("|----->> Matched-Normal ORGANOID Variants ... oncoplot(1) created and PNG copied for report : %s", HeatmapTargetCancerCodeMatchORG[k, "sample_info"]) )
                    system(sprintf("rm %s", tmp_png))
                    
                }
            }
        }
    }else{
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------

#---/ ANALYSIS : ONCOPLOTS 2 /----------------------------------------------------------------------
    message( sprintf('\n|--->> Run "ONCOPLOT-2" analysis <<%s| \n', paste(rep("-", 100), collapse="") ))

    #---/ TUMOR ONLY /----------------------------------------------------------
    if(  nrow(sinfo_preset_vars_tonly) > 0 )
    {
        HeatmapTargetCancerCodeList <- sort(unique(sinfo_preset_vars_tonly$sample_info))
        
        for( CODE in HeatmapTargetCancerCodeList )
        {
            cum_oncoplot_data <- WES_analysis.Get.Cumulative.Oncoplot.Data( CancerCode=CODE, ClientID=ClientID, VariantCallMethod=TonlyCallMode, SampleLimit=CumOncoplotSmapleLimit )

            if( !is.null(cum_oncoplot_data) )
            {
                CumOncoplotSampleInfo <- cum_oncoplot_data$var_sinfo
                # top30 gene oncoplot ------------------------------------------
                oncoplot_top30 <- WES_analysis.Draw.OncoPlots(
                    targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                selectGroupName  = "all", 
                                                selectCancerType = CODE , 
                                                VariantData      = cum_oncoplot_data$var_data, 
                                                topGenes         = 30, 
                                                geneSets         = NULL, 
                                                useBoth          = FALSE, 
                                                SampleInfo       = CumOncoplotSampleInfo
                                            ),
                    SampleInfo           = CumOncoplotSampleInfo,
                    SampleOrderByIdGroup = TRUE,
                    includeAllSamples    = TRUE
                )

                tmp_png <- sprintf("%s/cum_oncoplot_top30_Tonly_%s_tmp.png", ResultDir, CODE ) 
                fig_png <- sprintf("%s/cum_oncoplot_top30_Tonly_%s.png",     ResultDir, CODE ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(oncoplot_top30)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 100x20 %s", tmp_png, fig_png))
                system(sprintf("cp %s %s/%s.%s.top30.genes.Tonly.cumulative.samples.oncoplots.png", fig_png, ExportDir, OrderID, CODE ) ) # export
                # image size
                hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                if( hm_size[1] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "oncoplot_2",
                        variant_call    = "Tumor Only",
                        cancer_code     = CODE,
                        preset_name     = "top30",
                        sample_group    = "",
                        png             = sprintf("data/cum_oncoplot_top30_Tonly_%s.png", CODE)
                    )
                )
                system(sprintf("cp %s %s/cum_oncoplot_top30_Tonly_%s.png", fig_png, report_data_dir, CODE) )
                message(sprintf("|----->> Tumor-Only Variants ... oncoplot(2) created and PNG copied for report : %s", CODE) )
                system(sprintf("rm %s", tmp_png))

                # preset gene oncoplot -----------------------------------------
                CancerCodePresetNameList <- PresetList[which(PresetList$cancer_code == CODE), "preset_name"]
                for( PSN in CancerCodePresetNameList )
                {
                    oncoplot_preset <- WES_analysis.Draw.OncoPlots(
                        targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                    selectGroupName  = "all", 
                                                    selectCancerType = CODE , 
                                                    VariantData      = cum_oncoplot_data$var_data, 
                                                    topGenes         = NULL, 
                                                    geneSets         = ClientPresetGenes[[ PSN ]], 
                                                    useBoth          = FALSE, 
                                                    SampleInfo       = CumOncoplotSampleInfo
                                                ),
                        SampleInfo           = CumOncoplotSampleInfo,
                        SampleOrderByIdGroup = TRUE,
                        includeAllSamples    = TRUE
                    )

                    tmp_png <- sprintf("%s/cum_oncoplot_preset_%s_Tonly_%s_tmp.png", ResultDir, PSN, CODE ) 
                    fig_png <- sprintf("%s/cum_oncoplot_preset_%s_Tonly_%s.png",     ResultDir, PSN, CODE ) 
                    png( tmp_png, width=20, height=20, units="in", res=200)
                    draw(oncoplot_preset)
                    dev.off()
                    system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                    system(sprintf("cp %s %s/%s.%s.preset.%s.Tonly.cumulative.samples.oncoplots.png", fig_png, ExportDir, OrderID, CODE, PSN) ) # export
                    
                    # image size
                    hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                    hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                    if( hm_size[1] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    # fig list 
                    FIGUER_LIST <- rbind(
                        FIGUER_LIST,
                        data.frame( 
                            report_category = "oncoplot_2",
                            variant_call    = "Tumor Only",
                            cancer_code     = CODE,
                            preset_name     = PSN,
                            sample_group    = "",
                            png             = sprintf("data/cum_oncoplot_preset_%s_Tonly_%s.png", PSN, CODE)
                        )
                    )
                    system(sprintf("cp %s %s/cum_oncoplot_preset_%s_Tonly_%s.png", fig_png, report_data_dir, PSN, CODE) )
                    message(sprintf("|----->> Tumor-Only Variants ... oncoplot(1) created and PNG copied for report : %s", CODE) )
                    system(sprintf("rm %s", tmp_png))
                }
            }
        }
    }else{
        message("|----->> Tumor-Only Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ MATCH TS /------------------------------------------------------------
    if(  nrow(AnalysisMatchTS) > 0 )
    {
        HeatmapTargetCancerCodeList <- sort(unique(AnalysisMatchTS[which(AnalysisMatchTS$Match_TS), "cancer_code"]))
        
        for( CODE in HeatmapTargetCancerCodeList )
        {
            cum_oncoplot_data <- WES_analysis.Get.Cumulative.Oncoplot.Data( CancerCode=CODE, ClientID=ClientID, VariantCallMethod="TS_N", SampleLimit=CumOncoplotSmapleLimit )

            if( !is.null(cum_oncoplot_data) )
            {
                CumOncoplotSampleInfo <- cum_oncoplot_data$var_sinfo
                # top30 gene oncoplot ------------------------------------------
                oncoplot_top30 <- WES_analysis.Draw.OncoPlots(
                    targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                selectGroupName  = "all", 
                                                selectCancerType = CODE , 
                                                VariantData      = cum_oncoplot_data$var_data, 
                                                topGenes         = 30, 
                                                geneSets         = NULL, 
                                                useBoth          = FALSE, 
                                                SampleInfo       = CumOncoplotSampleInfo
                                            ),
                    SampleInfo           = CumOncoplotSampleInfo,
                    SampleOrderByIdGroup = TRUE,
                    includeAllSamples    = TRUE
                )

                tmp_png <- sprintf("%s/cum_oncoplot_top30_MatchTS_%s_tmp.png", ResultDir, CODE ) 
                fig_png <- sprintf("%s/cum_oncoplot_top30_MatchTS_%s.png",     ResultDir, CODE ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(oncoplot_top30)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 100x20 %s", tmp_png, fig_png))
                system(sprintf("cp %s %s/%s.%s.top30.genes.MatchNormalTissue.cumulative.samples.oncoplots.png", fig_png, ExportDir, OrderID, CODE ) ) # export
                # image size
                hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                if( hm_size[1] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "oncoplot_2",
                        variant_call    = "Matched Normal Tissue",
                        cancer_code     = CODE,
                        preset_name     = "top30",
                        sample_group    = "",
                        png             = sprintf("data/cum_oncoplot_top30_MatchTS_%s.png", CODE)
                    )
                )
                system(sprintf("cp %s %s/cum_oncoplot_top30_MatchTS_%s.png", fig_png, report_data_dir, CODE) )
                message(sprintf("|----->> Matched-Normal TISSUE Variants ... oncoplot(1) created and PNG copied for report : %s", CODE) )
                system(sprintf("rm %s", tmp_png))

                # preset gene oncoplot -----------------------------------------
                CancerCodePresetNameList <- PresetList[which(PresetList$cancer_code == CODE), "preset_name"]
                for( PSN in CancerCodePresetNameList )
                {
                    oncoplot_preset <- WES_analysis.Draw.OncoPlots(
                        targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                    selectGroupName  = "all", 
                                                    selectCancerType = CODE , 
                                                    VariantData      = data_SomaticVariants %>% filter( variant_call_mode == "TS_N" ), 
                                                    topGenes         = NULL, 
                                                    geneSets         = ClientPresetGenes[[ PSN ]], 
                                                    useBoth          = FALSE, 
                                                    SampleInfo       = CumOncoplotSampleInfo
                                                ),
                        SampleInfo           = CumOncoplotSampleInfo,
                        SampleOrderByIdGroup = TRUE,
                        includeAllSamples    = TRUE
                    )

                    tmp_png <- sprintf("%s/cum_oncoplot_preset_%s_MatchTS_%s_tmp.png", ResultDir, PSN, CODE ) 
                    fig_png <- sprintf("%s/cum_oncoplot_preset_%s_MatchTS_%s.png",     ResultDir, PSN, CODE ) 
                    png( tmp_png, width=20, height=20, units="in", res=200)
                    draw(oncoplot_preset)
                    dev.off()
                    system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                    system(sprintf("cp %s %s/%s.%s.preset.%s.MatchNormalTissue.cumulative.samples.oncoplots.png", fig_png, ExportDir, OrderID, CODE, PSN) ) # export
                    
                    # image size
                    hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                    hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                    if( hm_size[1] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    # fig list 
                    FIGUER_LIST <- rbind(
                        FIGUER_LIST,
                        data.frame( 
                            report_category = "oncoplot_2",
                            variant_call    = "Matched Normal Tissue",
                            cancer_code     = CODE,
                            preset_name     = PSN,
                            sample_group    = "",
                            png             = sprintf("data/cum_oncoplot_preset_%s_MatchTS_%s.png", PSN, CODE)
                        )
                    )
                    system(sprintf("cp %s %s/cum_oncoplot_preset_%s_MatchTS_%s.png", fig_png, report_data_dir, PSN, CODE) )
                    message(sprintf("|----->> Matched-Normal TISSUE Variants ... oncoplot(1) created and PNG copied for report : %s", CODE) )
                    system(sprintf("rm %s", tmp_png))
                }
            }
        }
    }else{
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

    #---/ MATCH ORG /-----------------------------------------------------------
    if(  nrow(AnalysisMatchORG) > 0 )
    {
        HeatmapTargetCancerCodeList <- sort(unique(AnalysisMatchORG[which(AnalysisMatchORG$Match_ORG), "cancer_code"]))
        
        for( CODE in HeatmapTargetCancerCodeList )
        {
            cum_oncoplot_data <- WES_analysis.Get.Cumulative.Oncoplot.Data( CancerCode=CODE, ClientID=ClientID, VariantCallMethod="ORG_N", SampleLimit=CumOncoplotSmapleLimit )

            if( !is.null(cum_oncoplot_data) )
            {
                CumOncoplotSampleInfo <- cum_oncoplot_data$var_sinfo
                # top30 gene oncoplot ------------------------------------------
                oncoplot_top30 <- WES_analysis.Draw.OncoPlots(
                    targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                selectGroupName  = "all", 
                                                selectCancerType = CODE , 
                                                VariantData      = cum_oncoplot_data$var_data, 
                                                topGenes         = 30, 
                                                geneSets         = NULL, 
                                                useBoth          = FALSE, 
                                                SampleInfo       = CumOncoplotSampleInfo
                                            ),
                    SampleInfo           = CumOncoplotSampleInfo,
                    SampleOrderByIdGroup = TRUE,
                    includeAllSamples    = TRUE
                )

                tmp_png <- sprintf("%s/cum_oncoplot_top30_MatchORG_%s_tmp.png", ResultDir, CODE ) 
                fig_png <- sprintf("%s/cum_oncoplot_top30_MatchORG_%s.png",     ResultDir, CODE ) 
                png( tmp_png, width=20, height=20, units="in", res=200)
                draw(oncoplot_top30)
                dev.off()
                system(sprintf("convert %s -trim -bordercolor white -border 100x20 %s", tmp_png, fig_png))
                system(sprintf("cp %s %s/%s.%s.top30.genes.MatchNormalOrganoid.cumulative.samples.oncoplots.png", fig_png, ExportDir, OrderID, CODE ) ) # export
                #image size
                hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                if( hm_size[1] > 1600 )
                {
                    resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                    system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                }
                # fig list 
                FIGUER_LIST <- rbind(
                    FIGUER_LIST,
                    data.frame( 
                        report_category = "oncoplot_2",
                        variant_call    = "Matchced Normal Organoid",
                        cancer_code     = CODE,
                        preset_name     = "top30",
                        sample_group    = "",
                        png             = sprintf("data/cum_oncoplot_top30_MatchORG_%s.png", CODE)
                    )
                )
                system(sprintf("cp %s %s/cum_oncoplot_top30_MatchORG_%s.png", fig_png, report_data_dir, CODE) )
                message(sprintf("|----->> Matched-Normal ORGANOID Variants ... oncoplot(1) created and PNG copied for report : %s", CODE) )
                system(sprintf("rm %s", tmp_png))

                # preset gene oncoplot -----------------------------------------
                CancerCodePresetNameList <- PresetList[which(PresetList$cancer_code == CODE), "preset_name"]
                for( PSN in CancerCodePresetNameList )
                {
                    oncoplot_preset <- WES_analysis.Draw.OncoPlots(
                        targetGeneCountTable = WES_analysis.Get.Oncoplot.Target.Genes( 
                                                    selectGroupName  = "all", 
                                                    selectCancerType = CODE , 
                                                    VariantData      = data_SomaticVariants %>% filter( variant_call_mode == "ORG_N" ), 
                                                    topGenes         = NULL, 
                                                    geneSets         = ClientPresetGenes[[ PSN ]], 
                                                    useBoth          = FALSE, 
                                                    SampleInfo       = CumOncoplotSampleInfo
                                                ),
                        SampleInfo           = CumOncoplotSampleInfo,
                        SampleOrderByIdGroup = TRUE,
                        includeAllSamples    = TRUE
                    )

                    tmp_png <- sprintf("%s/cum_oncoplot_preset_%s_MatchORG_%s_tmp.png", ResultDir, PSN, CODE ) 
                    fig_png <- sprintf("%s/cum_oncoplot_preset_%s_MatchORG_%s.png",     ResultDir, PSN, CODE ) 
                    png( tmp_png, width=20, height=20, units="in", res=200)
                    draw(oncoplot_preset)
                    dev.off()
                    system(sprintf("convert %s -trim -bordercolor white -border 20x20 %s", tmp_png, fig_png))
                    system(sprintf("cp %s %s/%s.%s.preset.%s.MatchNormalOrganoid.cumulative.samples.oncoplots.png", fig_png, ExportDir, OrderID, CODE, PSN) ) # export

                    # image size
                    hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
                    hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
                    if( hm_size[1] > 1600 )
                    {
                        resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
                        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
                    }
                    # fig list 
                    FIGUER_LIST <- rbind(
                        FIGUER_LIST,
                        data.frame( 
                            report_category = "oncoplot_2",
                            variant_call    = "Matched Normal Organoid",
                            cancer_code     = CODE,
                            preset_name     = PSN,
                            sample_group    = "",
                            png             = sprintf("data/cum_oncoplot_preset_%s_MatchORG_%s.png", PSN, CODE)
                        )
                    )
                    system(sprintf("cp %s %s/cum_oncoplot_preset_%s_MatchORG_%s.png", fig_png, report_data_dir, PSN, CODE) )
                    message(sprintf("|----->> Matched-Normal ORGANOID Variants ... oncoplot(1) created and PNG copied for report : %s", CODE) )
                    system(sprintf("rm %s", tmp_png))
                }
            }
        }
    }else{
        message("|----->> Matched-Normal TISSUE Variants ... no data to analysis. skip.")
    }
    #---------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------

#---/ ANALYSIS : INDIVIDUAL MATCHING /--------------------------------------------------------------
    message( sprintf('\n|--->> Run "INDIVIDUAL MATCHING" analysis <<%s| \n', paste(rep("-", 100), collapse="") ))
    
    indiv_match_heatmap <- WES_analysis.Indiv.Match.Heatmap( IndivMatchResult= data_IndivMatch %>% filter( match_res == "matched"), SampleInfo=sinfo )

    tmp_png <- sprintf("%s/indiv_match_heatmap_tmp.png", ResultDir ) 
    fig_png <- sprintf("%s/indiv_match_heatmap.png",     ResultDir ) 
    png( tmp_png, width=25, height=25, units="in", res=200)
    draw(indiv_match_heatmap)
    dev.off()
    system(sprintf("convert %s -trim -bordercolor white -border 100x20 %s", tmp_png, fig_png))
    # image size
    hm_size <- system(sprintf("identify %s | awk '{print $3}'", fig_png), intern=T)
    hm_size <- as.numeric(unlist(strsplit(hm_size, "x")))
    if( hm_size[1] > 1600 )
    {
        resizeRatio <- paste0(round(1600/hm_size[1] * 100 , 0), "%")
        system(sprintf("convert %s -resize %s %s", fig_png, resizeRatio, fig_png))
    }
    # fig list 
    FIGUER_LIST <- rbind(
        FIGUER_LIST,
        data.frame( 
            report_category = "individual_match",
            variant_call    = "",
            cancer_code     = "",
            preset_name     = "",
            sample_group    = "",
            png             = sprintf("data/indiv_match_heatmap.png")
        )
    )
    system(sprintf("cp %s %s/indiv_match_heatmap.png", fig_png, report_data_dir) )
    message(sprintf("|----->> Individual Matching ... Heatmap created and PNG copied for report : %s", SeqFolderID) )
    system(sprintf("rm %s", tmp_png))
#---------------------------------------------------------------------------------------------------

#---/ FIGURE LIST TABLE SAVE /----------------------------------------------------------------------
    saveRDS( FIGUER_LIST, file=sprintf("%s/data_figures_list.rds", ResultDir))
    system(sprintf("cp %s/data_figures_list.rds %s/data_figures_list.rds", ResultDir, report_data_dir))  
#---------------------------------------------------------------------------------------------------

#---/ REPORT RMD /----------------------------------------------------------------------------------
    # standard------------------------------------------------------------------
    system(sprintf("cp %s/DEFAULT_wes_standard_analysis_report_v2.Rmd %s/%s.Standard_Analysis_Report.Rmd", resource_dir, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_ORDER_ID/%s/g' %s/%s.Standard_Analysis_Report.Rmd", OrderID, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_SEQ_FOLDER_ID/%s/g' %s/%s.Standard_Analysis_Report.Rmd", SeqFolderID, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_LOW_QUAL_SEQ_ID/NULL/g' %s/%s.Standard_Analysis_Report.Rmd", WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_EXCLUDE_SEQ_ID/%s/g' %s/%s.Standard_Analysis_Report.Rmd", ifelse(is.null(ExcludeSeqID), "NULL", paste0('"',ExcludeSeqID,'"')), WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_INCLUDE_SEQ_ID/%s/g' %s/%s.Standard_Analysis_Report.Rmd", ifelse(is.null(IncludeSeqID), "NULL", paste0('"',IncludeSeqID,'"')), WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_CLIENT_ID/%s/g' %s/%s.Standard_Analysis_Report.Rmd", ClientID, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_VARIANT_CALL_MODE/%s/g' %s/%s.Standard_Analysis_Report.Rmd", TonlyCallMode, WES_REPORT_DIR, OrderID))

    # advanced -----------------------------------------------------------------
    system(sprintf("cp %s/DEFAULT_wes_advanced_analysis_report_v2.Rmd %s/%s.Advanced_Analysis_Report.Rmd" , resource_dir, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_ORDER_ID/%s/g' %s/%s.Advanced_Analysis_Report.Rmd", OrderID, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_SEQ_FOLDER_ID/%s/g' %s/%s.Advanced_Analysis_Report.Rmd", SeqFolderID, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_LOW_QUAL_SEQ_ID/NULL/g' %s/%s.Advanced_Analysis_Report.Rmd", WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_EXCLUDE_SEQ_ID/%s/g' %s/%s.Advanced_Analysis_Report.Rmd", ifelse(is.null(ExcludeSeqID), "NULL", paste0('"',ExcludeSeqID,'"')), WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_INCLUDE_SEQ_ID/%s/g' %s/%s.Advanced_Analysis_Report.Rmd", ifelse(is.null(IncludeSeqID), "NULL", paste0('"',IncludeSeqID,'"')), WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_CLIENT_ID/%s/g' %s/%s.Advanced_Analysis_Report.Rmd", ClientID, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_VARIANT_CALL_MODE/%s/g' %s/%s.Advanced_Analysis_Report.Rmd", TonlyCallMode, WES_REPORT_DIR, OrderID))
    
    # export data list pdf -----------------------------------------------------
    system(sprintf("cp %s/DEFAULT.Data_List.Rmd %s/%s.Data_List.Rmd", resource_dir, WES_REPORT_DIR, OrderID))
    system(sprintf("sed -i 's/DEFAULT_ORDER_ID/%s/g' %s/%s.Data_List.Rmd", OrderID, WES_REPORT_DIR, OrderID))

    # report config 
    if(       is.null(NgsLib)        ){ NGS_LIB <- "Twist Exome 2.0"  
    }else if( NgsLib == "twist"      ){ NGS_LIB <- "Twist Exome 2.0" 
    }else if( NgsLib == "sureselect" ){ NGS_LIB <- "SureSelect V6" 
    }else{ NGS_LIB <- NgsLib }
    system(sprintf("cp %s/wes_report_config.yaml %s/wes_report_config.yaml", resource_dir, WES_REPORT_DIR))
    system(sprintf("sed -i 's/NGS_LIBRARY_HERE/%s/g' %s/wes_report_config.yaml", NGS_LIB, WES_REPORT_DIR))


    # export file list ---------------------------------------------------------
    #> seq_folder seqid sample_name tonly match_ts match_org tonly_gf, germline
    SampleID_Tonly    <- sinfo$seq_id
    SampleID_MatchTS  <- sinfo %>% filter( sample_group %in% analysis_run_info$MatchTS$SampleGroup  ) %>% .$seq_id
    SampleID_MatchORG <- sinfo %>% filter( sample_group %in% analysis_run_info$MatchORG$SampleGroup ) %>% .$seq_id

    ExportSampleList <- data.frame(
        SeqFolderID = SeqFolderID,
        SeqID       = sinfo$seq_id,
        SampleName  = sinfo$sample_name
    )
    ExportSampleList$Tonly    <- SampleID_Tonly[ match(ExportSampleList$SeqID, SampleID_Tonly) ]
    ExportSampleList$MatchTS  <- SampleID_MatchTS[ match(ExportSampleList$SeqID, SampleID_Tonly) ]
    ExportSampleList$MatchORG <- SampleID_MatchORG[ match(ExportSampleList$SeqID, SampleID_Tonly) ]

    ExportSampleList$Tonly    <- sapply(ExportSampleList$Tonly,    function(w) ifelse( is.na(w), "NO", "YES") )
    ExportSampleList$MatchTS  <- sapply(ExportSampleList$MatchTS,  function(w) ifelse( is.na(w), "NO", "YES") )
    ExportSampleList$MatchORG <- sapply(ExportSampleList$MatchORG, function(w) ifelse( is.na(w), "NO", "YES") )
    
    ExportSampleList$Germline <- "NO"
    ExportSampleList$TonlyGF  <- "NO"

    MetaDir <- sprintf("%s/%s/meta", BaseDir, SeqFolderID)
    write.table( ExportSampleList, sprintf("%s/%s_export_data_list.tsv", MetaDir, SeqFolderID), quote=FALSE, col.names=TRUE,  row.names=FALSE, sep="\t")
    write.table( ExportSampleList, sprintf("%s/%s_export_data.list",     MetaDir, SeqFolderID), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" " )

    # copy data for client export
    # system(sprintf("/storage/home/kangsm/runScripts/NGS_service.ExportData_v2.sh --baseDir %s --seqFolder %s --checksum md5", BaseDir, SeqFolderID))

#---------------------------------------------------------------------------------------------------

## FINAL REPORT

    # /mnt/d/01.NGS_RUO_SERVICE/ReportResources/NGS_service.PrepareFinalReport.sh \
    #     --InputPdf "advanced_report.pdf" \
    #     --LocalRun "true" \
    #     --ReportType "advanced" \
    #     --NgsApplication "wes" \
    #     --OrderId "REPORT-v2-TEST" \
    #     --RemoveTmp "false"

    
    # /storage/home/kangsm/runScripts/NGS_service.PrepareFinalReport.sh 
    #     --InputPdf "advanced_report.pdf" \
    #     --LocalRun "false" \
    #     --ReportType "advanced" \
    #     --NgsApplication "wes" \
    #     --OrderId "REPORT-v2-TEST" \
    #     --RemoveTmp "false"
#---------------------------------------------------------------------------------------------------

  












    