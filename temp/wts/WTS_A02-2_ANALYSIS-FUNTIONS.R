
#---| PACKAGES |-------------------------------------------------------------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("Hmisc"))
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#===| EXPRESSION-VALUES-MATRIX |=====================================================================================================================#
#' @description create expression-value matrix from saved expression values.
#' @param SeqIdList seq-id list.
#' @param ValueType expression-value type. 
#' @param DataSource data source. (saved data location)
#' @param dbHost database host ( if DataSource = 'db' )
#' @param dbUser database user id ( if DataSource = 'db' )
#' @param dbPasswd database user password ( if DataSource = 'db' )
#' @param dbName database name ( if DataSource = 'db' )
#' @param CountAsInteger convert count value into integer (only for ValueType = 'count')
#' @param Convert2Matrix convert expression-value list into matrix format
#' @export
    Create.Expression.Values.Matrix <- function( SeqIdList, ValueType=c("count","tpm","fpkm"), DataSource=c("db","rds"), 
        dbHost, dbUser, dbPasswd, dbName, CountAsInteger=TRUE , Convert2Matrix=TRUE 
    )
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("RMySQL")) 
        source("/data/wts/params/ruo_wts_db.R")
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
        on.exit(dbDisconnect(dbCon))
        #-----------------------------------------------------------------------
        SEQ_IDS <- paste(paste0("'", SeqIdList, "'"), collapse=",")
        #-----------------------------------------------------------------------
        if( DataSource == "db" )
        {
            message("|---> selected data-source = 'database'")
            message(sprintf("|---> expression-value type = '%s'", ValueType))
            if( ValueType == "tpm" ){ 
                EV <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,tpm FROM expr_values WHERE seq_id IN (%s)", SEQ_IDS)) %>% dplyr::rename( expr=3 )
            }else if( ValueType == "fpkm" ){
                EV <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,fpkm FROM expr_values WHERE seq_id IN (%s)", SEQ_IDS)) %>% dplyr::rename( expr=3 )
            }else{
                EV <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,read_count FROM expr_values WHERE seq_id IN (%s)", SEQ_IDS)) %>% dplyr::rename( expr=3 )
            }                
            #-------------------------------------------------------------------
            EV_List <- setNames( EV$expr, nm=EV$seq_id )
            EV_List <- lapply( EV_List, function(ev) unlist(strsplit(ev, ";")) )
            EV_List <- lapply( EV_List, function(ev)
            {
                names(ev) <- sapply(ev, function(y) unlist(strsplit(y, ","))[1])
                if( all( ValueType == "count" & CountAsInteger) )
                {    
                    nev <- sapply(ev, function(y) round(as.numeric(unlist(strsplit(y, ","))[2]),0))
                }else{
                    nev <- sapply(ev, function(y) as.numeric(unlist(strsplit(y, ","))[2]))
                }
                return(nev)
            })
        }else{
            message("|---> selected data-source = 'local RDS file'")
            message(sprintf("|---> expression-value type = '%s'", ValueType))
            if( ValueType == "tpm" ){ 
                EV <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,tpm_rds_path FROM expr_values WHERE seq_id IN (%s)", SEQ_IDS)) %>% dplyr::rename( path=3 )
            }else if( ValueType == "fpkm" ){
                EV <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,fpkm_rds_path FROM expr_values WHERE seq_id IN (%s)", SEQ_IDS)) %>% dplyr::rename( path=3 )
            }else{
                EV <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,read_count_rds_path FROM expr_values WHERE seq_id IN (%s)", SEQ_IDS)) %>% dplyr::rename( path=3 )
            } 
            #-------------------------------------------------------------------
            EV_List <- list()
            for( k in 1:nrow(EV) )
            {
                EV_List <- c(EV_List, readRDS( EV[k, 'path'] ))
            }
            if( all( ValueType == "count" & CountAsInteger) )
            {
                EV_List <- lapply( EV_List, function(ev) round(ev, 0) )
            }
        } 
        #-----------------------------------------------------------------------
        if( Convert2Matrix )
        {
            message("|---> expression-value output format = 'matrix'")
            EXPR_VALUES <- as.data.frame( EV_List )
        }else{
            message("|---> expression-value output format = 'list'")
            EXPR_VALUES <- EV_List
        }
        #-----------------------------------------------------------------------
        return(EXPR_VALUES)
    }
#====================================================================================================================================================#

#===| SAMPLE CLUSTERING ANALYSIS MODULES |===========================================================================================================#
#---/ UMAP /---------------------------------------------------------------------------------------#
#' @description 
#' @param ExprValueMatrix input expression-value matrix
#' @param SeqIdList            sample-list for umap-clustering analysis
#' @param SeqIdGroup           sample group. if not provided, name of SeqIdList will be used. if no names in SeqIdList, analysis will be run as No Groups. 
#' @param LowExpressionFilter  filtering low expression level features or not. cut-off is pre-defined on each value-type.
#' @param ValueType            expression-value type. one of 'count', 'tpm' and 'fpkm'. default = 'tpm'
#' @param GeneFilter           select genes (features) for analysis. NULL = use all genes. "deg" = use only DEGs. 0~100 = quantile % of variance sorted order.
#' @param DegList              DEG list vector (ens-geneid, ponly used when GeneFilter = "deg" )
#' @param ComponentN           UMAP-analysis components, 2 or 3. default = 3
#' @param UmapNeighbors        UMAP-analysis neighbors #N. 
#' @param GroupColors          sample group colors. if NULL, pre-defined colors will be used.
#' @export 
    run.UMAP.Clustering <- function( ExprValueMatrix, SeqIdList, SeqIdGroup=NULL, LowExpressionFilter=TRUE, ValueType="tpm", 
        GeneFilter=NULL, DegList=NULL, ComponentN=3, UmapNeighbors=NULL, GroupColors=NULL
    )
    {
        suppressPackageStartupMessages(library("umap"))
        # sample group colors --------------------------------------------------
        SampleGroupColors <- c(
            "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#80B1D3","#A65628","#F781BF","#999999","#66C2A5","#8DA0CB",
            "#A6D854","#FFD92F","#E5C494","#FFFF33","#FDB462","#B3DE69","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#FCCDE5"
        )
        #-----------------------------------------------------------------------
        SeqIdList <- unlist(strsplit(SeqIdList, ","))
        #-----------------------------------------------------------------------
        evm <- ExprValueMatrix[, SeqIdList]
        #-----------------------------------------------------------------------
        if( LowExpressionFilter )
        {
            message("|---> performing low-expression-filtering...")
            if( ValueType == "tpm" ){ cutoff = 0.5 }else if( ValueType == "fpkm" ){ cutoff = 1 }else{ cutoff = 10 }
            evm <- evm[apply( evm, 1, function(z) !all( z <= cutoff )), ]
        }
        #-----------------------------------------------------------------------
        if( is.null(GeneFilter) )
        {
            message("|---> no gene-filter. use all genes.")
            umapInput <- evm
        }else if( GeneFilter == "deg" ){
            message("|---> use only DEGs.")
            umapInput <- evm[ DegList, ]
        }else{
            message(sprintf("|---> use only variance top %s percentage genes", GeneFilter*100))
            GeneVariance <- rev(sort(apply( log2(evm+1), 1, var)))
            umapInput    <- evm[names(GeneVariance[ GeneVariance >= quantile(GeneVariance, 1-GeneFilter) ]), ]
        }
        # UMAP config ----------------------------------------------------------
        umapConfig                      <- umap.defaults
        umapConfig$n_components         <- ComponentN
        umapConfig$metric               <- "cosine"
        umapConfig$n_epochs             <- 1000
        umapConfig$negative_sample_rate <- ifelse( length(SeqIdList) > 10, 5, length(SeqIdList)/2 )
        umapConfig$random_state         <- 1234
        umapConfig$transform_state      <- 123456  
        if( is.null(UmapNeighbors) )
        { umapConfig$n_neighbors <- as.integer(sqrt(length(SeqIdList))) }else{ umapConfig$n_neighbors <- UmapNeighbors } 
        
        # Run UMAP Analysis -----------------------------------------------------------------------#
            message(sprintf("|---> run UMAP analysis : componenets = %s", ComponentN))
            umap_res <- suppressWarnings(umap(t(as.matrix(umapInput)), config=umapConfig))
        #------------------------------------------------------------------------------------------#

        umap_analysis_result <- as.data.frame(umap_res$layout) 
        if( ComponentN == 2 )
        { 
            umap_analysis_result   <- umap_analysis_result %>% dplyr::rename(X=1, Y=2) 
            umap_analysis_result$Z <- NA
        }else{ 
            umap_analysis_result   <- umap_analysis_result %>% dplyr::select(c(1:3)) %>% dplyr::rename(X=1, Y=2, Z=3) 
        }
        # sample groups -------------------------------------------------------#
        if( is.null(SeqIdGroup) )
        {
            if( !is.null(names(SeqIdList)) )
            {
                umap_analysis_result$group <- (names(SeqIdList)[ match( rownames(umap_analysis_result), SeqIdList) ] )
            }else{
                umap_analysis_result$group <- SeqIdList
            }
        }else{
            umap_analysis_result$group = SeqIdGroup
        }
        # add colors ----------------------------------------------------------#
        if( !is.null(GroupColors) )
        {
            if( length(which(umap_analysis_result$group %nin% names(GroupColors))) > 0 )
            {
                message("|---!> input group colors are not sufficient to represent all sample groups. pre-defiend colors will be used.")
                names(SampleGroupColors)    <- unique(umap_analysis_result$group)
                umap_analysis_result$colors <- SampleGroupColors[ umap_analysis_result$group ]
            }else{
                umap_analysis_result$colors <- GroupColors[ umap_analysis_result$group ]
            }
        }else{
            names(SampleGroupColors)    <- unique(umap_analysis_result$group)
            umap_analysis_result$colors <- SampleGroupColors[ umap_analysis_result$group ]
        }
        #-----------------------------------------------------------------------
        umap_analysis_result$seq_id              <- rownames(umap_analysis_result)
        umap_analysis_result$analysis_method     <- "UMAP"
        umap_analysis_result$analysis_components <- umapConfig$n_components
        umap_analysis_result$distance_metric     <- umapConfig$metric

        #-----------------------------------------------------------------------
        return(umap_analysis_result)
    }
#--------------------------------------------------------------------------------------------------#
#---/ PCA /----------------------------------------------------------------------------------------#
#' @description run PCA analysis usin expression value matrix
#' @param ExprValueMatrix      input expression-value matrix
#' @param SeqIdList            sample-list for umap-clustering analysis
#' @param ValueType            expression-value type. one of 'count', 'tpm' and 'fpkm'. default = 'tpm'
#' @param SeqIdGroup           sample group. if not provided, name of SeqIdList will be used. if no names in SeqIdList, analysis will be run as No Groups. 
#' @param LowExpressionFilter  filtering low expression level features or not. cut-off is pre-defined on each value-type.
#' @param Log2                 perform log2 transformation of expression value matrix
#' @param PcaComponenetsN      PCA analysis component #N. default = 5
#' @param GroupColors          sample group colors. SHOULD be named vector. if NULL, pre-defined colors will be used.
#' @export 
    run.PCA.Clustering <- function( 
        ExprValueMatrix, SeqIdList, SeqIdGroup, ValueType="tpm", LowExpressionFilter=TRUE, Log2=TRUE, PcaComponenetsN=5, GroupColors=NULL
    )
    {
        suppressPackageStartupMessages(library("FactoMineR"))
        # sample group colors --------------------------------------------------
        SampleGroupColors <- c(
            "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#80B1D3","#A65628","#F781BF","#999999","#66C2A5","#8DA0CB",
            "#A6D854","#FFD92F","#E5C494","#FFFF33","#FDB462","#B3DE69","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#FCCDE5"
        )
        # sample id list -------------------------------------------------------
        SeqIdList <- unlist(strsplit(SeqIdList, ","))
        # exprssion value matrix -----------------------------------------------
        evm <- ExprValueMatrix[, SeqIdList]
        # low expression filtering ---------------------------------------------
        if( LowExpressionFilter )
        {
            message("|---> performing low-expression-filtering...")
            if( ValueType == "tpm" ){ cutoff = 0.5 }else if( ValueType == "fpkm" ){ cutoff = 1 }else{ cutoff = 10 }
            evm <- evm[apply( evm, 1, function(z) !all( z <= cutoff )), ]
        }
        # log-transformation ---------------------------------------------------
        if( Log2 ){ evm <- log2(evm + 1) }
        
        # RUN PCA ANALYSIS ------------------------------------------------------------------------#
            pca_res <- PCA( t(evm), graph=FALSE, scale.unit=TRUE, ncp=PcaComponenetsN )
            # sample stats ----------------------------------------------------#
            pca_sample_stats <- pca_res$ind$coord %>% as.data.frame() %>% 
                dplyr::select(c(1:3)) %>%
                dplyr::rename(X=1, Y=2, Z=3) %>% 
                mutate(
                    seq_id               = rownames(.), 
                    analysis_method      = "PCA",
                    analysis_components = 5,
                    distance_metric      = ""
                )
            # feature stats ---------------------------------------------------#
            pca_feature_stats <- data.frame( 
                comp        = 1:PcaComponenetsN,
                var_contrib = rownames(pca_res$var$contrib)[1:PcaComponenetsN],
                var_pct     = round(pca_res$eig[1:PcaComponenetsN, "percentage of variance"],2),
                var_cum_pct = round(pca_res$eig[1:PcaComponenetsN, "cumulative percentage of variance"],2),
                X           = pca_res$var$coord[1:PcaComponenetsN, "Dim.1"],
                Y           = pca_res$var$coord[1:PcaComponenetsN, "Dim.2"],
                Z           = pca_res$var$coord[1:PcaComponenetsN, "Dim.3"]
            )
        #------------------------------------------------------------------------------------------#

        # add sample group ----------------------------------------------------#
        if( is.null(SeqIdGroup) )
        {
            if( !is.null(names(SeqIdList)) )
            {
                pca_sample_stats$group <- (names(SeqIdList)[ match( pca_sample_stats$seq_id, SeqIdList) ] )
            }else{
                pca_sample_stats$group <- SeqIdList
            }
        }else{
            pca_sample_stats$group = SeqIdGroup
        }
        # add colors ----------------------------------------------------------#
        if( !is.null(GroupColors) )
        {
            if( length(which(pca_sample_stats$group %nin% names(GroupColors))) > 0 )
            {
                message("|---!> input group colors are not sufficient to represent all sample groups. pre-defiend colors will be used.")
                names(SampleGroupColors) <- unique(pca_sample_stats$group)
                pca_sample_stats$colors  <- SampleGroupColors[ pca_sample_stats$group ]
            }else{
                pca_sample_stats$colors <- GroupColors[ pca_sample_stats$group ]
            }
        }else{
            names(SampleGroupColors) <- unique(pca_sample_stats$group)
            pca_sample_stats$colors  <- SampleGroupColors[ pca_sample_stats$group ]
        }
        # final result list ---------------------------------------------------#
        PCA_RESULT <- list(
            sample_stats  = pca_sample_stats,
            feature_stats = pca_feature_stats,
            analysis_object = pca_res
        )
        #----------------------------------------------------------------------#
        return(PCA_RESULT)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Sample Relative Distance /-------------------------------------------------------------------#
#' @description calculate distances between expression values profiles
#' @param SeqFolderID     REQUIRED. analysis batch id (SEQ_FOLDER)
#' @param ExprMatrix      REQUIRED. Input expression values matrix.
#' @param DistanceMetric  REQUIRED. distance metric. supprted metric = 'cosine', 'pearson', 'spearman', 'kendall', 'euclidean', 'manhattan',' supremum'
#' @param DbImport        import result into database. default = TRUE
#' @export 
    Expression.Profile.Distance <- function( SeqFolderID=NULL, ExprMatrix=NULL, DistanceMetric=NULL, DbImport = TRUE )
    {
        if( is.null(SeqFolderID)    ){ stop("|---!!! No analysis batch id (SeqFolder). REQUIRED. Check again. STOPPED.")   }
        if( is.null(ExprMatrix)     ){ stop("|---!!! No input epression values profiles. REQUIRED. Check again. STOPPED.") }
        if( is.null(DistanceMetric) ){ stop("|---!!! No distance metric. REQUIRED. Check again. STOPPED.")                 }
        #-----------------------------------------------------------------------
        suppressPackageStartupMessages(library("lsa"))
        # scaling function -----------------------------------------------------
        ScaleDown <- function( x, y )
        {
            cv <- (x-min(x))*((max(y)-min(y))/(max(x)-min(x))) + min(y)
            cv <- round(cv, 6)
            return(cv)
        }
        # low expression remove ------------------------------------------------
        mat <- ExprMatrix[rowSums(ExprMatrix) != 0 , ]
        mat <- mat[apply( mat, 1, function(y) length(which( y <= 0.5 )) != ncol(mat) ), ]
        # sample pairwise list -------------------------------------------------
        cal_dist_list <- as.data.frame(t(combn(colnames(mat), 2)))
        # calculate distance ---------------------------------------------------
        if( DistanceMetric == "cosine" ){  
            cal_dist_list$distance <- apply( cal_dist_list[,1:2], 1, function(pf) round(1-cosine(mat[ ,pf][,1],  mat[ ,pf][,2]), 6) ) 
        }else if( DistanceMetric == "pearson" ){
            cal_dist_list$distance <- apply( cal_dist_list[,1:2], 1, function(pf) round(1-cor( x=mat[ ,pf][,1], y=mat[ ,pf][,2], method="pearson"), 6) )
        }else if( DistanceMetric == "spearman" ){
            cal_dist_list$distance <- apply( cal_dist_list[,1:2], 1, function(pf) round(1-cor( x=mat[ ,pf][,1], y=mat[ ,pf][,2], method="spearman"), 6) )
        }else if( DistanceMetric == "kendall" ){
            cal_dist_list$distance <- apply( cal_dist_list[,1:2], 1, function(pf) round(1-cor( x=mat[ ,pf][,1], y=mat[ ,pf][,2], method="kendall"),6) )
        }else if( DistanceMetric == "euclidean" ){
            cal_dist_list$distance <- apply( cal_dist_list[,1:2], 1, function(pf) sqrt(sum((mat[ ,pf][,1] - mat[ ,pf][,2])^2)) )
            cal_dist_list$distance <- ScaleDown(c(0,cal_dist_list$distance), c(0,1))[-1]
        }else if( DistanceMetric == "manhattan" ){
            cal_dist_list$distance <- apply( cal_dist_list[,1:2], 1, function(pf) sum(abs(mat[ , pf][,1] - mat[ , pf][,2])) )
            cal_dist_list$distance <- ScaleDown(c(0,cal_dist_list$distance), c(0,1))[-1]
        }else if( DistanceMetric == "supremum" ){
            cal_dist_list$distance <- apply( cal_dist_list[,1:2], 1, function(pf) max(abs(mat[ ,pf][,1] - mat[ ,pf][,2])) )
            cal_dist_list$distance <- ScaleDown(c(0,cal_dist_list$distance), c(0,1))[-1]
        }else{
            stop("|---!!! DistanceMetric is not supported. Check again. STOPPED.")
        }
        #----------------------------------------------------------------------#
        colnames(cal_dist_list) <- c("seq_id_1","seq_id_2", "distance")
        #----------------------------------------------------------------------#
        cal_dist_list$seq_folder <- SeqFolderID
        cal_dist_list$metric     <- DistanceMetric
        #----------------------------------------------------------------------#  
        if( DbImport )
        {
            source("/data/wts/params/ruo_wts_db.R")
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM sample_dist_calculation WHERE seq_folder = '%s' AND metric = '%s'", SEQ_FOLDER, DistanceMetric ) )
            UpdateData      <- dbWriteTable(dbCon, name="sample_dist_calculation", value=cal_dist_list, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
        }
        #----------------------------------------------------------------------#
        return( cal_dist_list )
    }
#--------------------------------------------------------------------------------------------------#
#---/ Convert Sample Distance Data Into Matrix Format /--------------------------------------------# 
#' @description convert sample distance calculation result table into matrix 
#' @param DistanceResultTable  REQUIRED. input result table
#' @param SampleOrderList      sample order vector. if NULL, sample id will sorted by alphabetical order
#' @export 
    Convert.Distance.Table.To.Matrix <- function( DistanceResultTable=NULL, SampleOrderList=NULL )
    {
        if( is.null(DistanceResultTable) ){ stop("|---!!! No distance calculation result table found. REQUIRED. Check again. STOPPED.") }
        #-----------------------------------------------------------------------
        sample_list        <- union( DistanceResultTable$seq_id_1, DistanceResultTable$seq_id_2 )
        names(sample_list) <- sample_list
        #-----------------------------------------------------------------------
        if( !is.null(SampleOrderList) )
        {
            sample_list <- sample_list[ SampleOrderList ]
        }else{
            sample_list <- sort(sample_list)
        }
        #-----------------------------------------------------------------------
        res_mat <- matrix( rep(0, length(sample_list)*length(sample_list)), ncol=length(sample_list) )
        rownames(res_mat) = colnames(res_mat) = sample_list
        for( i in 1:nrow(DistanceResultTable) )
        {
            res_mat[ DistanceResultTable[i,"seq_id_1"], DistanceResultTable[i,"seq_id_2"] ] <- DistanceResultTable[i, "distance"]
            res_mat[ DistanceResultTable[i,"seq_id_2"], DistanceResultTable[i,"seq_id_1"] ] <- DistanceResultTable[i, "distance"]
        }
        labels        <- names(SampleOrderList)
        names(labels) <- SampleOrderList
        rownames(res_mat) <- labels[rownames(res_mat)]
        colnames(res_mat) <- labels[colnames(res_mat)]
        #-----------------------------------------------------------------------
        return(res_mat)
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| DEG ANALYSIS MODULES |=========================================================================================================================#
#---/ DEG Analysis Info Table /--------------------------------------------------------------------#
#' @description prepare DEG-analysis-info table for DEG analysis
#' @param SeqFolderId           sample batch id ( = SEQ_FOLDER )
#' @param SampleInfoTable       sample info table ( = SINFO )
#' @param SeqIdList             sample id ( = seq id )
#' @param SeqIdGroup            sample group id
#' @param ControlTreatPairList  REQUIRED. ctrl and treat 'SeqIdGroup' pair in list format. eg. list(list(ctrl=..., treat=...),...)
#' @param GenomeAssembly        reference genome assembly version. either 'hg19' or 'hg38'. default = "hg19"
#' @param GetInfoFromDb         load exist DEG-analysis-info from database. default = TRUE. if FALSE, DEG-analysis-info will be made newly.
#' @param CreateInfo            create DEG-analysis-info newly based on sample info. default = FALSE
#' @param DbImport              import DEG-analysis-info into database. default = FALSE
#' @export 
    Prepapre.WTS.Analysis.Info.Table <- function( SeqFolderId=NULL, SampleInfoTable=NULL, SeqIdList=NULL, SeqIdGroup=NULL, ControlTreatPairList=list(), 
        GenomeAssembly="hg19", GetInfoFromDb=TRUE, CreateInfo=TRUE, DbImport=TRUE 
    )
    {
        suppressPackageStartupMessages(library("RMySQL"))
        suppressPackageStartupMessages(library("Hmisc"))
        source("/data/wts/params/ruo_wts_db.R")
        # check control-treat pair info ----------------------------------------
        if( length(ControlTreatPairList) == 0 )
        { stop("|---!!! control : treat paired group list is REQURIED but NOT FOUND. please check again. STOPPED.") }
        # select method --------------------------------------------------------
        if( !GetInfoFromDb )
        {
            message("|--->> create DEG analysis info using user input sample info.")
            if( is.null(SampleInfoTable) )
            {
                if( is.null(SeqIdList) ) { stop("|---!!! no seq ids. please check or change method to get data from database. STOPPED.") 
                }else{
                    SampleList <- SeqIdList
                    if( is.null(SeqIdGroup) ) {
                        if( !is.null(names(SeqIdList)) ){ SampleGroup <- names(SeqIdList) }else{ 
                            SampleGroup       <- SeqIdList 
                            names(SampleList) <- SeqIdList
                        }
                    }else{ names(SampleList) <- SeqIdGroup }                
                }
            }else{
                SampleList <- SampleInfoTable$seq_id
                names(SampleList) <- SampleInfoTable$sample_group
            }
            # 
            ANALYSIS_INFO_TABLE <- data.frame()
            for( k in 1:length(ControlTreatPairList) )
            {
                SeqIDs_CTRL     <- SampleList[ names(SampleList) %in% unlist(strsplit(ControlTreatPairList[[k]]$ctrl, ",")) ]
                SeqIDs_TREAT    <- SampleList[ names(SampleList) %in% unlist(strsplit(ControlTreatPairList[[k]]$treat, ",")) ]
                SeqIDs_CTRL     <- SeqIDs_CTRL[SeqIDs_CTRL %nin% ControlTreatPairList[[k]]$rm_ctrl ]
                SeqIDs_TREAT    <- SeqIDs_TREAT[SeqIDs_TREAT %nin% ControlTreatPairList[[k]]$rm_treat ]
                GroupName_CTRL  <- ifelse(is.null(ControlTreatPairList[[k]]$ctrl_group_name), ControlTreatPairList[[k]]$ctrl, ControlTreatPairList[[k]]$ctrl_group_name)
                GroupName_TREAT <- ifelse(is.null(ControlTreatPairList[[k]]$treat_group_name), ControlTreatPairList[[k]]$treat, ControlTreatPairList[[k]]$treat_group_name)

                ANALYSIS_INFO_TABLE <- rbind(ANALYSIS_INFO_TABLE, 
                    rbind(
                        data.frame(
                            seq_folder      = SeqFolderId,
                            analysis_id     = k,
                            analysis_group  = "CTRL",
                            sample_id       = SeqIDs_CTRL,
                            sample_group    = GroupName_CTRL,
                            genome_assembly = GenomeAssembly
                        ),
                        data.frame(
                            seq_folder      = SeqFolderId,
                            analysis_id     = k,
                            analysis_group  = "TREAT",
                            sample_id       = SeqIDs_TREAT,
                            sample_group    = GroupName_TREAT,
                            genome_assembly = GenomeAssembly
                        ) 
                    )
                )
            }
            ANALYSIS_INFO_TABLE$analysis_id <- paste0("AID", formatC(ANALYSIS_INFO_TABLE$analysis_id, digits=0, flag="0", width=3))
            #
            if( DbImport )
            {
                dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM deg_analysis_info WHERE seq_folder = '%s'", SeqFolderId))
                UpdateData      <- dbWriteTable(dbCon, name="deg_analysis_info", value=ANALYSIS_INFO_TABLE, row.names=FALSE, append=TRUE)
                dbDisconnect(dbCon)
            }
        }else{
            if( CreateInfo )
            {
                message("|--->> create DEG analysis info using sample info in database.")
                # get sample info from database --------------------------------
                dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_INFO )
                SampleInfo <- dbGetQuery( dbCon, sprintf("SELECT seq_folder,seq_id,sample_group FROM seqid_info WHERE seq_folder = '%s'", SeqFolderId))
                dbDisconnect(dbCon)
                # check sample info --------------------------------------------
                if( nrow(SampleInfo) == 0 )
                {
                    stop(sprintf("|---!!! SEQ_FOLDER_ID : %s ... >>> no samples info in database. please check again. STOPPED.", SeqFolderId))
                }else{
                    SampleList  <- SampleInfo$seq_id
                    SampleGroup <- SampleInfo$sample_group
                    names(SampleList) <- SampleGroup
                }
                # create DEG analysis info -------------------------------------
                ANALYSIS_INFO_TABLE <- data.frame()
                for( k in 1:length(ControlTreatPairList) )
                {
                    SeqIDs_CTRL     <- SampleList[ names(SampleList) %in% unlist(strsplit(ControlTreatPairList[[k]]$ctrl, ",")) ]
                    SeqIDs_TREAT    <- SampleList[ names(SampleList) %in% unlist(strsplit(ControlTreatPairList[[k]]$treat, ",")) ]
                    SeqIDs_CTRL     <- SeqIDs_CTRL[SeqIDs_CTRL %nin% ControlTreatPairList[[k]]$rm_ctrl ]
                    SeqIDs_TREAT    <- SeqIDs_TREAT[SeqIDs_TREAT %nin% ControlTreatPairList[[k]]$rm_treat ]
                    GroupName_CTRL  <- ifelse(is.null(ControlTreatPairList[[k]]$ctrl_group_name), ControlTreatPairList[[k]]$ctrl, ControlTreatPairList[[k]]$ctrl_group_name)
                    GroupName_TREAT <- ifelse(is.null(ControlTreatPairList[[k]]$treat_group_name), ControlTreatPairList[[k]]$treat, ControlTreatPairList[[k]]$treat_group_name)

                    ANALYSIS_INFO_TABLE <- rbind(ANALYSIS_INFO_TABLE, 
                        rbind(
                            data.frame(
                                seq_folder      = SeqFolderId,
                                analysis_id     = k,
                                analysis_group  = "CTRL",
                                sample_id       = SeqIDs_CTRL,
                                sample_group    = GroupName_CTRL,
                                genome_assembly = GenomeAssembly
                            ),
                            data.frame(
                                seq_folder      = SeqFolderId,
                                analysis_id     = k,
                                analysis_group  = "TREAT",
                                sample_id       = SeqIDs_TREAT,
                                sample_group    = GroupName_TREAT,
                                genome_assembly = GenomeAssembly
                            ) 
                        )
                    )
                }
                ANALYSIS_INFO_TABLE$analysis_id <- paste0("AID", formatC(ANALYSIS_INFO_TABLE$analysis_id, digits=0, flag="0", width=3))
                # import into database -----------------------------------------
                if( DbImport )
                {
                    dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                    DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM deg_analysis_info WHERE seq_folder = '%s'", SeqFolderId))
                    UpdateData      <- dbWriteTable(dbCon, name="deg_analysis_info", value=ANALYSIS_INFO_TABLE, row.names=FALSE, append=TRUE)
                    dbDisconnect(dbCon)
                }
            }else{
                message("|--->> load DEG analysis info from database.")
                dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                ANALYSIS_INFO_TABLE <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_info WHERE seq_folder = '%s'", SeqFolderId))
                dbDisconnect(dbCon)
            }
        }
        #-----------------------------------------------------------------------
        return(ANALYSIS_INFO_TABLE)
    }
#--------------------------------------------------------------------------------------------------#
#---/ DEG Analysis Meta Table /--------------------------------------------------------------------#
#' @description create DEG analysis set meta table. (final format = list) 
#' @param SeqFolderId       seq-folder. ( = analysis batch, project batch )
#' @param LoadFromDb        load analysis set info table from database. default = TRUE. ( prioroty = 1 )
#' @param CreateFromObject  create analysis set list from analysis set info table input. default = FALSE ( prioroty = 2 )
#' @param AnalysisInfoTable analysis set info table input. required when CreateFromObject = TRUE
#' @param ControlGroupID    group-id of control group
#' @param TreatGroupIDs     group-ids of treat group
#' @param ExcludeSeqIDs     seq-ids that should be removed from analysis
#' @param ExcludeGroupIDs   group-ids that should be removed from analysis
#' @param GenomeAssembly    genome assembly version (need for DEG-analysis gene mapping)
#' @export 
    Create.DEG.Analysis.Set.Table.List <- function( SeqFolderId=NULL, LoadFromDb=TRUE, CreateFromObject=FALSE, AnalysisInfoTable=NULL,  
        ControlGroupID=NULL, TreatGroupIDs=NULL, ExcludeSeqIDs=NULL, ExcludeGroupIDs=NULL, GenomeAssembly="hg19"
    )
    {
        suppressPackageStartupMessages(library("RMySQL"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        source("/data/wts/params/ruo_wts_db.R")
        #-----------------------------------------------------------------------
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        #-----------------------------------------------------------------------
        if( LoadFromDb )
        {
            if( is.null(SeqFolderId) ){ 
                stop("|---!!! no SeqFolderId. REQUIRED when LoadFromDb = TRUE. please check again. STOPPED.")
            }else{
                dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
                AnalysisInfoTable <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_info WHERE seq_folder = '%s'", SeqFolderId))
                dbDisconnect(dbCon)
                ##
                ANALYSIS_SET_LIST <- lapply( indexing(AnalysisInfoTable$analysis_id), function(aid) AnalysisInfoTable %>% filter( analysis_id == aid ) )
                ##
            }
        }else if( CreateFromObject ){
            if( !is.data.frame(AnalysisInfoTable) )
            {
                stop("|---!!! AnalysisInfoTable is not 'data.frame'. SHOULD BE 'data.frame' when CreateFromObject = TRUE. please check again. STOPPED.")
            }else{
                ##
                ANALYSIS_SET_LIST <- lapply( indexing(AnalysisInfoTable$analysis_id), function(aid) AnalysisInfoTable %>% filter( analysis_id == aid ) )
                ##
            }
        }else{
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_INFO )
            sinfo <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SeqFolderId))
            dbDisconnect(dbCon)
            # exclude groups ---------------------------------------------------
            if( !is.null(ExcludeGroupIDs) )
            {
                ExcludeGroupIDs <- unlist(strsplit(ExcludeGroupIDs, ","))
                sinfo <- sinfo %>% filter( sample_group %nin% ExcludeGroupIDs) 
            }
            # exclude seq-ids --------------------------------------------------
            if( !is.null(ExcludeSeqIDs) )
            { 
                ExcludeSeqIDs <- unlist(strsplit(ExcludeSeqIDs, ","))
                sinfo <- sinfo %>% filter( seq_id %nin% ExcludeSeqIDs) 
            }
            # Group Check ------------------------------------------------------
            if( length(unique(sinfo$sample_group)) < 2 ){ stop("No more than 2 Groups") }        
            # CTRL-GROUP -------------------------------------------------------
            if( is.null(ControlGroupID) ){ ControlGroupID <- unique(sinfo$sample_group)[1] 
            }else{
                if( ControlGroupID %nin% sinfo$sample_group )
                { stop("CTRL-GROUP is not in sample-group. please check again. stopped.") }else{ ControlGroupID <- ControlGroupID }
            }
            # TREAT-GROUPS -----------------------------------------------------
            if( is.null(TreatGroupIDs) )
            { 
                TreatGroupIDs <- sinfo %>% filter( sample_group != ControlGroupID ) %>% .$sample_group %>% unique() 
            }else{
                TreatGroupIDs <- unlist(strsplit( TreatGroupIDs, "," ))
                if( all(TreatGroupIDs %in% sinfo$sample_group)  )
                {
                    message("|--->> all treat-group-ids are exist in sample info. use all.")
                }else{
                    message("|--->> some treat-group-ids are NOT EXIST in sample info. use only avaiable ones.")
                    TreatGroupIDs <- TreatGroupIDs[ TreatGroupIDs %in% sinfo$sample_group ]
                }
            }
            # create analysis set list -----------------------------------------
            TreatGroupIndex   <- indexing(TreatGroupIDs)
            ANALYSIS_SET_LIST <- lapply( TreatGroupIndex, function(idx) 
            {
                analysisSet <- rbind(
                    data.frame(
                        seq_folder      = SeqFolderId,
                        analysis_id     = "",
                        analysis_group  = "CTRL",
                        sampel_id       = sinfo[which(sinfo$sample_group == ControlGroupID), "seq_id"],
                        sampel_group    = ControlGroupID,
                        genome_assembly = GenomeAssembly
                    ),
                    data.frame(
                        seq_folder      = SeqFolderId,
                        analysis_id     = "",
                        analysis_group  = "TREAT",
                        sampel_id       = sinfo[which(sinfo$sample_group == idx), "seq_id"],
                        sampel_group    = idx,
                        genome_assembly = GenomeAssembly
                    )
                )
                return(analysisSet)
            })
            for( k in 1:length(ANALYSIS_SET_LIST) )
            { 
                ASL <- ANALYSIS_SET_LIST[[k]]
                ASL$analysis_id <- paste0("AID", formatC(k, digits=0, flag="0", width=3)) 
                ANALYSIS_SET_LIST[[k]] <- ASL
                rm(ASL)
            }
        }
        # 
        ANALYSIS_SET_LIST <- lapply( ANALYSIS_SET_LIST, function(ASL){
            ASL$analysis_group <- factor(ASL$analysis_group, levels=c("CTRL","TREAT"))
            return(ASL)
        })
        return(ANALYSIS_SET_LIST)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Calculate DESeq2 Score /---------------------------------------------------------------------#
#' @description calculate DESeq2 DEG scores 
#' @param DESeqAnalysisResult DESeq2 analysis result table 
#' @param ResultSort sort result by score ( descending absolute score value )
#' @export 
    Calculate.DESeq.Scores <- function( DESeqAnalysisResult, ResultSort=TRUE )
    {
        suppressPackageStartupMessages(library("dplyr"))
        #-----------------------------------------------------------------------
        max_pvalue_score <- max(-log10(DESeqAnalysisResult[which(DESeqAnalysisResult$pvalue != 0 & !is.na(DESeqAnalysisResult$pvalue)), 'pvalue'])+2)
        #-----------------------------------------------------------------------
        DESeqAnalysisResult$score <- apply( DESeqAnalysisResult[,c("log2FoldChange","stat","pvalue")], 1, function(w)
        {
            if( is.na(w[3])     ){ pvalueScore <- 0 
            }else if( w[3] == 0 ){ pvalueScore <- max_pvalue_score
            }else{                 pvalueScore <- -log10(w[3]) 
            }
            score <- sqrt( abs(w[1]) * abs(w[2]) * pvalueScore )
            if( all(!is.na(w[1]) & w[1] < 0) ){ score = score * (-1) }
            score <- round(score, 2)
            return(score)
        })
        #-----------------------------------------------------------------------
        if( ResultSort )
        {
            DESeqAnalysisResult <- DESeqAnalysisResult %>% arrange( dplyr::desc(abs(score)) )
        }
        #-----------------------------------------------------------------------
        return(DESeqAnalysisResult)
    }
#--------------------------------------------------------------------------------------------------#
#---/ DESeq2 Run /---------------------------------------------------------------------------------#
#' @description DEG analysis using DESeq2
#' @param ExprCountMatrix rnaseq integer count matrix. major input. 
#' @param ExprTpmMatrix tpm value matrix for filtering. if NULL, will be created from data-source.
#' @param AnalysisDesign analysis design table.                   
#' @param FilterByCount filtering low expression based on count. default min=10, min.total=15
#' @param FilterByTPM filtering low expression based on TPM. default 0.5
#' @param DegAssignMethod DEG assign method. either 'ratio' or 'cutoff'. default = 'cutoff'. if the number of DEGs are more than 10% of gene space size, 10% ratio will be forcely used
#' @param PvalueCutOff p-value cut-off. default = 0.05
#' @param FDRCutOff FDR cut-off. default = 0.05
#' @param logFoldChangeCutOff fold-change cut-off. default = log2(1.3)
#' @param DegRatio DEG assign ratio when DegAssignMethod = 'ratio'. max ratio = 10% (0.1)
#' @param MatrixCreationDataSource tpm matrixc creation data-source. 'db' or 'rds'. default = 'db' 
#' @param WriteParamsFile    write deg analysis param table as TSV file. default = TRUE
#' @param ImportParamsIntoDb import deg analysis param table into database. default = TRUE 
#' @export 
    run.DESeq2.DEG.Analysis <- function( ExprCountMatrix, ExprTpmMatrix=NULL, AnalysisDesign, FilterByCount=TRUE, FilterByTPM=TRUE, DegAssignMethod="cutoff",
        PvalueCutOff=0.05, FDRCutOff=0.05, logFoldChangeCutOff=log2(1.3), DegRatio=0.1, MatrixCreationDataSource="db", 
        WriteParamsFile=TRUE, ParamsFileSaveDir=NULL, ImportParamsIntoDb = TRUE, GeneInfo=GENE_INFO
    )
    {
        suppressPackageStartupMessages(library("DESeq2"))
        suppressPackageStartupMessages(library("edgeR"))
        suppressPackageStartupMessages(library("dplyr"))
        # count-based filtering ------------------------------------------------
        if( FilterByCount )
        {
            message("|---> filtering based on read-counts...")
            FilterDesign <- model.matrix( ~ AnalysisDesign$analysis_group )
            keepGenes    <- filterByExpr( ExprCountMatrix, design=FilterDesign, group=AnalysisDesign$analysis_group, min.count=10, min.total.count=15 )
        }
        # tpm-based filtering --------------------------------------------------
        if( FilterByTPM )
        {
            message("|---> filtering based on TPM...")
            if( is.null(ExprTpmMatrix) )
            { 
                message(sprintf("|---> no input TPM-matrix. create from data-source. selected data-source = '%s'", MatrixCreationDataSource))
                if( MatrixCreationDataSource == "db" )
                {    
                    suppressPackageStartupMessages(library("RMySQL"))
                    ExprTpmMatrix <- Create.Expression.Values.Matrix( SeqIdList = AnalysisDesign$sample_id, ValueType="tpm", DataSource="db", Convert2Matrix=TRUE )
                }else{
                    ExprTpmMatrix <- Create.Expression.Values.Matrix( SeqIdList = AnalysisDesign$sample_id, ValueType="tpm", DataSource="rds", Convert2Matrix=TRUE )
                }
            }
            rmGenes <- apply(ExprTpmMatrix, 1, function(w) length(which(w < 0.5)) == length(w) )
        }else{
            rmGenes = c()
        }
        # DESeq2 Run -----------------------------------------------------------
        message("|---> run DESeq2 analysis ...")
        CountMat  <- ExprCountMatrix[, AnalysisDesign$sample_id]
        dds       <- suppressMessages(DESeqDataSetFromMatrix(countData=CountMat, colData=AnalysisDesign, design= ~ analysis_group))
        dds       <- suppressMessages(estimateSizeFactors(dds))
        dds       <- suppressMessages(estimateDispersions(dds))
        dds_res   <- suppressMessages(DESeq(dds, betaPrior=FALSE))
        deseq_res <- suppressMessages(results(dds_res, pAdjustMethod = "BH", independentFiltering=TRUE))
        #-----------------------------------------------------------------------
        message("|---> add gene-filtering results ...")
        DESEQ_RES <- data.frame( ens_geneid = rownames(deseq_res), deseq_res )
        if( FilterByCount ){ DESEQ_RES$keep <- keepGenes[ DESEQ_RES$ens_geneid ] }else{ DESEQ_RES$keep = TRUE }
        if( FilterByTPM   ){ DESEQ_RES[ which(rmGenes), "keep" ] <- FALSE }
        DESEQ_RES$ens_geneid = sapply(DESEQ_RES$ens_geneid, function(y) unlist(strsplit(y, "\\."))[1] )
        DESEQ_RES <- DESEQ_RES %>% left_join(., GeneInfo[,c("ens_geneid","symbol","entrez","locus_type","hgnc_id")], by="ens_geneid")
        # deseq2 score ---------------------------------------------------------
        message("|---> calculate DEG-scores ...")
        DESEQ_RES <- Calculate.DESeq.Scores( DESeqAnalysisResult=DESEQ_RES, ResultSort=TRUE )
        DESEQ_RES_GENES <- DESEQ_RES %>% filter( locus_type == "gene with protein product" ) %>% filter( baseMean != 0 )
        # DEG tagging ----------------------------------------------------------
        GeneSpaceSize <- DESEQ_RES_GENES %>% filter( keep ) %>% nrow(.)
        message("|---> make DEG-tags ...")
        if( DegAssignMethod == "cutoff" )
        {
            message("|---> DEG assign method : cut-off values")
            DESEQ_RES_GENES$deg <- ""
            DESEQ_RES_GENES[which( 
                DESEQ_RES_GENES$log2FoldChange >= logFoldChangeCutOff &
                DESEQ_RES_GENES$pvalue         <= PvalueCutOff        
            ), "deg" ] <- "UP"
            DESEQ_RES_GENES[which( 
                DESEQ_RES_GENES$log2FoldChange <= -logFoldChangeCutOff &
                DESEQ_RES_GENES$pvalue         <= PvalueCutOff        
            ), "deg" ] <- "DOWN"
        }else{
            message("|---> DEG assign method : gene space ratio")
            deg_numbers      <- as.integer(GeneSpaceSize * DegRatio)
            deg_score_cutoff <- rev(sort(abs(DESEQ_RES_GENES$score)))[deg_numbers]
            DESEQ_RES_GENES$deg    <- ""
            DESEQ_RES_GENES[which( abs(DESEQ_RES_GENES$score) >= deg_score_cutoff & DESEQ_RES_GENES$log2FoldChange > 0 ), "deg"] <- "UP"
            DESEQ_RES_GENES[which( abs(DESEQ_RES_GENES$score) >= deg_score_cutoff & DESEQ_RES_GENES$log2FoldChange < 0 ), "deg"] <- "DOWN"
        }
        # deg cut-off params ---------------------------------------------------
            initialDegRatio <- round(length(which(DESEQ_RES_GENES$deg != ""))/GeneSpaceSize , 3)
            deg_cutoff_param_table <- data.frame(
                seq_folder     = unique(AnalysisDesign$seq_folder  ),
                analysis_id     = unique(AnalysisDesign$analysis_id ),
                cutoff_method  = DegAssignMethod,
                cutoff_logfc   = ifelse( DegAssignMethod == "cutoff", logFoldChangeCutOff, NA),
                cutoff_pval    = ifelse( DegAssignMethod == "cutoff", PvalueCutOff, NA),
                cutoff_fdr     = ifelse( DegAssignMethod == "cutoff", FDRCutOff, NA),
                cutoff_ratio   = ifelse( DegAssignMethod == "cutoff", NA, DegRatio),
                gene_space     = GeneSpaceSize,
                deg_ratio_1st  = initialDegRatio,
                up_genes_1st   = length(which(DESEQ_RES_GENES$deg == "UP")),
                down_genes_1st = length(which(DESEQ_RES_GENES$deg == "DOWN"))
            )
        if( length(which(DESEQ_RES_GENES$deg != "")) > GeneSpaceSize * 0.11 )
        {
            # adjust DEGs : too many DEGs --------------------------------------
            message(sprintf("|---> too many DEGs. DEG ratio = %s", round(length(which(DESEQ_RES_GENES$deg != ""))/GeneSpaceSize, 2) ))
            message("|---> reducing DEGs to 10% of gene-space-size...")
            deg_numbers      <- as.integer(GeneSpaceSize * 0.1)
            deg_score_cutoff <- rev(sort(abs(DESEQ_RES_GENES$score)))[deg_numbers]
            DESEQ_RES_GENES$deg    <- ""
            DESEQ_RES_GENES[which( abs(DESEQ_RES_GENES$score) >= deg_score_cutoff & DESEQ_RES_GENES$log2FoldChange > 0 ), "deg"] <- "UP"
            DESEQ_RES_GENES[which( abs(DESEQ_RES_GENES$score) >= deg_score_cutoff & DESEQ_RES_GENES$log2FoldChange < 0 ), "deg"] <- "DOWN"
            # add 2nd adjust deg cut-off params
            deg_cutoff_param_table$cutoff_adjust  <- 1
            deg_cutoff_param_table$adjust_method  <- "ratio"
            deg_cutoff_param_table$adjust_cutoff  <- 0.1
            deg_cutoff_param_table$up_genes_2nd   <- length(which(DESEQ_RES_GENES$deg == "UP"))
            deg_cutoff_param_table$down_genes_2nd <- length(which(DESEQ_RES_GENES$deg == "DOWN"))
        }else if( length(which(DESEQ_RES_GENES$deg != "")) == 0 )
        {
            message(sprintf("|---> no DEGs. adjust deg cut-off to ratio = 0.05"))
            message("|---> default foldchange and p-value cut-off will be also applied after adjustment with ratio.")
            deg_numbers          <- as.integer(GeneSpaceSize * 0.05)
            deg_score_cutoff     <- rev(sort(abs(DESEQ_RES_GENES$score)))[deg_numbers]
            default_logfc_cutoff <- log2(1.3)
            default_pval_cutoff  <- 0.05
            DESEQ_RES_GENES$deg    <- ""
            DESEQ_RES_GENES[which( abs(DESEQ_RES_GENES$score) >= deg_score_cutoff & 
                DESEQ_RES_GENES$log2FoldChange >= default_logfc_cutoff & DESEQ_RES_GENES$pvalue <= default_pval_cutoff ), "deg"] <- "UP"
            DESEQ_RES_GENES[which( abs(DESEQ_RES_GENES$score) >= deg_score_cutoff & 
                DESEQ_RES_GENES$log2FoldChange <= -default_logfc_cutoff &  DESEQ_RES_GENES$pvalue <= default_pval_cutoff ), "deg"] <- "DOWN"
            # add 2nd adjust deg cut-off params
            deg_cutoff_param_table$cutoff_adjust  <- 1
            deg_cutoff_param_table$adjust_method  <- "ratio"
            deg_cutoff_param_table$adjust_cutoff  <- 0.05
            deg_cutoff_param_table$up_genes_2nd   <- length(which(DESEQ_RES_GENES$deg == "UP"))
            deg_cutoff_param_table$down_genes_2nd <- length(which(DESEQ_RES_GENES$deg == "DOWN"))       
        }else{
            message(sprintf("|---> no more adjustment."))
            deg_cutoff_param_table$cutoff_adjust  <- 0
            deg_cutoff_param_table$adjust_method  <- NA
            deg_cutoff_param_table$adjust_cutoff  <- NA
            deg_cutoff_param_table$up_genes_2nd   <- NA
            deg_cutoff_param_table$down_genes_2nd <- NA
        }
        # write params as file -------------------------------------------------
        if( WriteParamsFile )
        {
            if( is.null(ParamsFileSaveDir) )
            {
                PARAM_FILE_NAME <- sprintf("./%s.%s.DEG.anaysis.params.tsv", unique(AnalysisDesign$seq_folder), unique(AnalysisDesign$analysis_id ) )
            }else{
                PARAM_FILE_NAME <- sprintf("%s/%s.%s.DEG.anaysis.params.tsv", ParamsFileSaveDir, unique(AnalysisDesign$seq_folder), unique(AnalysisDesign$analysis_id ) )
            }
            write.table( deg_cutoff_param_table, PARAM_FILE_NAME, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
        }
        # import params into database ------------------------------------------
        if( ImportParamsIntoDb )
        {
            source("/data/wts/params/ruo_wts_db.R")
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData <- dbGetQuery(dbCon, sprintf("DELETE FROM deg_analysis_params WHERE seq_folder = '%s' AND analysis_id = '%s'", 
                unique(AnalysisDesign$seq_folder), unique(AnalysisDesign$analysis_id) )
            )
            UpdateData <- dbWriteTable(dbCon, name="deg_analysis_params", value=deg_cutoff_param_table, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
        }
        #-----------------------------------------------------------------------
        DESEQ_RES_GENES$genome_assembly <- unique(AnalysisDesign$genome_assembly)
        #-----------------------------------------------------------------------
        DESEQ_RES_GENES <- DESEQ_RES_GENES %>% dplyr::rename( base_mean=baseMean, log2_foldchange=log2FoldChange, lfc_se=lfcSE )
        #-----------------------------------------------------------------------
        message("|---> DESeq2 DEG analysis DONE.")
        return(DESEQ_RES_GENES)
    } 
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| PATHWAY ANALYSIS MODULES |=====================================================================================================================#
#---/ Load DEG Analysis Result From Database or RDS /----------------------------------------------#
#' @description load DEG analysis result based on analysis batch and analysis-id
#' @param SeqFolder    REQUIRED. analysis batch id ( = seq_folder )
#' @param AnalysisId   REQUIRED. analysis id. 
#' @param DataSource   data source. either 'db' or 'rds'. default = 'db'
#' @param GeneInfo     REQUIRED. gene information table. 
#' @param GlobalRdsDir RDS base dir. default = /data/wts/RDS
#' @export 
    load.DEG.Analysis.Profile <- function( SeqFolder=NULL, AnalysisId=NULL, DataSource='db', GeneInfo=NULL, GlobalRdsDir="/data/wts/RDS" )
    {
        # seq-folder check -----------------------------------------------------
        if( is.null(SeqFolder) ){ stop("|---!!! no seq-folder-id. REQUIRED. please check again. STOPPED.") }
        # analysis id check ----------------------------------------------------
        if( is.null(AnalysisId) ){ stop("|---!!! no analysis-id. REQUIRED. please check again. STOPPED.") }
        # gene info ------------------------------------------------------------
        if( is.null(GeneInfo) ){ stop("|---!!! no gene info table. REQUIRED. please check again. STOPPED.") }
        #-------------------------------------------------------------------------------------------
        message(sprintf("|---> %s - %s ", SeqFolder, AnalysisId))
        #-------------------------------------------------------------------------------------------
        if( DataSource == "db" )
        {
            # load data from database ------------------------------------------
            message("|---> load DEG analysis result data from database.")
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            deg_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_result WHERE seq_folder = '%s' AND analysis_id = '%s'", SeqFolder, AnalysisId))
            dbDisconnect(dbCon)

            deg_score         <- sapply(unlist(strsplit(deg_res$deg_score, ";")),  function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
            names(deg_score)  <- sapply(unlist(strsplit(deg_res$deg_score, ";")),  function(z) unlist(strsplit(z, ","))[1] ) 
            deg_logfc         <- sapply(unlist(strsplit(deg_res$deg_logfc, ";")),  function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
            names(deg_logfc)  <- sapply(unlist(strsplit(deg_res$deg_logfc, ";")),  function(z) unlist(strsplit(z, ","))[1] ) 
            deg_pvalue        <- sapply(unlist(strsplit(deg_res$deg_pvalue, ";")), function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
            names(deg_pvalue) <- sapply(unlist(strsplit(deg_res$deg_pvalue, ";")), function(z) unlist(strsplit(z, ","))[1] ) 
            deg_tag           <- sapply(unlist(strsplit(deg_res$deg_tag, ";")),    function(z) unlist(strsplit(z, ","))[2] ) 
            names(deg_tag)    <- sapply(unlist(strsplit(deg_res$deg_tag, ";")),    function(z) unlist(strsplit(z, ","))[1] ) 
            
            DEG_RES <- data.frame(
                seq_folder      = SeqFolder,
                analysis_id     = AnalysisId,
                ens_geneid      = names(deg_score),
                log2_foldchange = deg_logfc,
                pvalue          = deg_pvalue,
                score           = deg_score,
                deg             = deg_tag,
                keep            = TRUE,
                GeneInfo[match(names(deg_score), GeneInfo$ens_geneid), c("symbol","entrez","locus_type")] 
            )
            DEG_RES[which(DEG_RES$deg == "noDEG"), "deg"] = ""
        }else{
            # load data from RDS -----------------------------------------------
            message("|---> load DEG analysis result data from RDS.")
            DEG_RES <- readRDS(sprintf("%s/RDS_DegAnalysis/%s.%s.DEG.analysis.result.table.rds", GlobalRdsDir, SeqFolder, AnalysisId))
        }
        #-----------------------------------------------------------------------
        return(DEG_RES)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Load PathwaySet Data from Database and Convert Into Gene List /------------------------------#
#' @description prepare pathway set genes list 
#' @param PathwaySetID   pathway-set-id. PathwaySetID has higher priority than PathwaySetName.
#' @param PathwaySetName pathway-set-name
#' @param Identifier     gene identifier. either 'entrez' or 'symbol'. default = 'entrez'
#' @param DataSource     pathwayset data source. either "db" or "rds". default = 'db'
#' @param PathSubset     pathway subset by data-source. (optional), if multiple subsets, insert ',' seprarated string as input.
#' @export 
    load.Pathway.Gene.List <- function( PathwaySetID=NULL, PathwaySetName=NULL, Identifier="entrez", DataSource="db", PathSubset=NULL )
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("RMySQL"))
        source("/data/wts/params/ruo_wts_db.R")
        dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db="public_database" )
        on.exit(dbDisconnect(dbCon))
        #-----------------------------------------------------------------------
        if( all(is.null(PathwaySetID) & is.null(PathwaySetName)) )
        { 
            stop("|---!!! STOPPED. no both pathway set id and pathway set name. at least one of the pathway set identifier shoud be exist.") 
        }else{
            if( is.null(PathwaySetID) )
            {
                PathwaySetID <- dbGetQuery(dbCon, sprintf("SELECT set_id FROM pathset_info WHERE set_name = '%s'", PathwaySetName))$set_id
                if( length( PathwaySetID) == 0 ){ stop("|---!!! STOPPED. no pathway set has such name" ) }
            }
        }
        #-----------------------------------------------------------------------
        if( DataSource == "rds" )
        {
            RDS_FILE <- dbGetQuery(dbCon, sprintf("SELECT rds_path FROM pathset_info WHERE set_id = '%s'", PathwaySetID))$rds_path
            if( file.exists(RDS_FILE) )
            {
                PSL <- readRDS(RDS_FILE)
            }else{
                stop("|---!!! no pathway sets. please check set-id or set-name again. STOPPED.")
            }  
        }else{
            PSL <- dbGetQuery(dbCon, sprintf("SELECT path_id,pathway,data_source,ext_id,entrez,hgnc_symbol FROM pathwaysets WHERE set_id = '%s'", PathwaySetID))
            #---------------------------------------------------------------
            if( nrow(PSL) == 0 ) { stop("|---!!! no pathway sets. please check set-id or set-name again. STOPPED.") }
        }
        #-----------------------------------------------------------------------
        if( !is.null(PathSubset) )
        {
            PathSubset <- unlist(strsplit(PathSubset, ","))
            #-------------------------------------------------------------------
            if( length(intersect(tolower(PathSubset), tolower(PSL$data_source))) > 0 )
            {
                message(sprintf("|---> subset identifier found. only subset will be used. found subset = %s", PathSubset))
                PSL <- PSL %>% filter( data_source %in% PathSubset )
            }else{
                message("|---> subset identifier NOT found. ALL pathway set will be used.")
            }
        }
        #-----------------------------------------------------------------------
        if( Identifier == "symbol" ){ pathwayList <- as.list(PSL$hgnc_symbol) }else{ pathwayList <- as.list(PSL$entrez) }
        names(pathwayList) <- PSL$path_id 
        pathwayList <- lapply( pathwayList, function(pw) unlist(strsplit(pw, ";")))
        pathwayInfo <- PSL[,c("path_id","pathway","data_source","ext_id")]
        #-----------------------------------------------------------------------
        PATHWAY_SET_LIST <- list(
            pathwayGenes = pathwayList,
            pathwayInfo  = pathwayInfo
        )
        #-----------------------------------------------------------------------
        return(PATHWAY_SET_LIST)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Calulcate Hypergeometric-Test and Enrichment Factor /----------------------------------------#
#' @description 
#' @param PathwayGenesList pathway gene list
#' @param InputGeneSet input gene set (eg. DEGs)
#' @param GeneSpace total gene list (background all genes)
#' @param Threads cpu cores for parallel analysis
#' @param MinPathSize minimum pathway set genes size 
#' @param MaxPathSize maximum pathway set genes size
#' @export  
    Calculate.Pathway.Gene.Enrich <- function( PathwayGenesList, InputGeneSet, GeneSpace, Threads, MinPathSize=10, MaxPathSize=500 )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("parallel"))
        #-----------------------------------------------------------------------
        GeneSpaceSize <- length(GeneSpace)
        InputSetSize  <- length(InputGeneSet)
        # pathway set filtering ------------------------------------------------
        PathListFilterd <- mclapply( PathwayGenesList, function(pw) pw[ pw %in% GeneSpace ], mc.cores=Threads )
        PathSizeCheck   <- mclapply( PathListFilterd, function(pw) ifelse( length(pw) < MinPathSize | length(pw) > MaxPathSize , FALSE, TRUE ), mc.cores=Threads )
        PF_Pathway      <- PathListFilterd[ which( PathSizeCheck == TRUE ) ]
        #-----------------------------------------------------------------------
        PathEnrichResultList <- mclapply( PF_Pathway, function(PW) 
        {
            PW_HITS     <- length( intersect( PW, InputGeneSet ) )
            PW_BG_RATIO <- length( intersect( PW, GeneSpace )    )
            # hypergeometric p-value -------------------------------------------
            hyperGeometricPval <- phyper(PW_HITS, PW_BG_RATIO, (GeneSpaceSize-PW_BG_RATIO), InputSetSize, lower.tail=FALSE)
            # enrichment factor ------------------------------------------------
            enrichFactor <- (PW_HITS/InputSetSize) / (PW_BG_RATIO/GeneSpaceSize)
            #-------------------------------------------------------------------
            enrich_res <- data.frame(
                pathway_size  = length(PW),
                pvalue        = hyperGeometricPval,
                enrich_factor = enrichFactor,
                enrich_genes  = PW_HITS,
                gene_list     = paste( intersect( PW, InputGeneSet ), collapse=";" )
            )
            #-------------------------------------------------------------------
            return(enrich_res)
        }, mc.cores = Threads)
        #-----------------------------------------------------------------------
        return(PathEnrichResultList)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Run Pathway Analysis /-----------------------------------------------------------------------#
#' @description
#' @param DegAnalysisResult deg analysis result table
#' @param PathwaySetID      pathway-set-id. default = 'pw01' ( = consensuspathdb )
#' @param PathwaySetID      pathway-set-name. default = 'cpdb' ( = consensuspathdb )
#' @param PathSubset        pathway subset by data-source. (optional)
#' @param AddMergeDegRun    use 'UP' and 'DOWN' merged gene list as input geneset. default = TRUE
#' @param GeneIdentifier    gene identifier. either 'entrez' or 'symbol'. default = 'entrez'
#' @param Threads           cpu cores for parallel processing. default = 15. 
#' @export 
    run.DEG.Pathway.Analysis <- function( DegAnalysisResult, PathwaySetID="pw01", PathwaySetName="cpdb", PathSubset=NULL, 
        AddMergeDegRun=TRUE, GeneIdentifier="entrez", Threads=15 
    )
    { 
        # get pathway set ------------------------------------------------------
        message("|--->>> Loading pathwayset data ...")
        PATHWAY_SET <- load.Pathway.Gene.List( PathwaySetID=PathwaySetID, PathwaySetName=PathwaySetName, Identifier=GeneIdentifier, DataSource="db", PathSubset=PathSubset )
        PathwayList <- PATHWAY_SET$pathwayGenes
        PathwayInfo <- PATHWAY_SET$pathwayInfo
        # analysis input gene set index ----------------------------------------
        PathwayAnalysisIndex <- list( "UP"="UP", "DOWN"="DOWN" )
        deg_stats  <- as.data.frame(table(DegAnalysisResult$deg)) %>% dplyr::rename( degtype=1, genes=2 ) %>% 
            mutate( degtype=as.character(degtype) ) %>% filter( degtype %in% c("UP","DOWN") ) 
        analyzable <- deg_stats %>% filter( genes > 10 ) %>% .$degtype
        PathwayAnalysisIndex <- PathwayAnalysisIndex[ analyzable ]
        if( AddMergeDegRun )
        {
            if( sum(deg_stats$genes) > 10 ){
                PathwayAnalysisIndex <- c(PathwayAnalysisIndex, list("ALL_DEGs"=c("UP","DOWN")) )
                message("|--->>> Analysis using merged DEG (Up+Down) is INCLUDED.")
            }else{
                message("|--->>> Analysis using merged DEG (Up+Down) is EXCLUDED due to low DEGs.")
            }
        }else{
            message("|--->>> Analysis using merged DEG (Up+Down) is EXCLUDED due to analysis-config.")
        }
        if( length(PathwayAnalysisIndex) >= 1 )
        {       
            # gene identifier apply -------------------------------------------- 
            if( GeneIdentifier == "symbol" )
            {
                PathwayAnalysisIndex <- lapply(PathwayAnalysisIndex, function(gt)
                { DegAnalysisResult %>% filter( keep ) %>% filter( !is.na(symbol) ) %>% filter( deg %in% gt ) %>% .$symbol %>% unique() })
                GENE_SPACE <- DegAnalysisResult %>% filter( keep ) %>% filter( !is.na(symbol) ) %>% .$symbol %>% unique()
                message("|--->>> 'Gene Symbol' is used for Gene-Indentifier.")
                message(sprintf("|--->>> Total %s gene as gene-space.", length(GENE_SPACE)))
            }else{
                PathwayAnalysisIndex <- lapply(PathwayAnalysisIndex, function(gt)
                { DegAnalysisResult %>% filter( keep ) %>% filter( !is.na(entrez) ) %>% filter( deg %in% gt ) %>% .$entrez %>% unique() })
                GENE_SPACE <- DegAnalysisResult %>% filter( keep ) %>% filter( !is.na(entrez) ) %>% .$entrez %>% unique()
                message("|--->>> 'EntrezID' is used for Gene-Indentifier.")
                message(sprintf("|--->>> Total %s gene as gene-space.", length(GENE_SPACE)))
            }
            # run Over-Representation-Analysis (enrichment analysis) -----------
            message("|--->>> perform hypergeometric-test and enrichment-factor calculation ....")
            PathwayAnalysisResultListSet <- mclapply( PathwayAnalysisIndex, function(gn) 
            {
                PathwayAnlysisResList <- Calculate.Pathway.Gene.Enrich( 
                    PathwayGenesList = PathwayList, 
                    InputGeneSet     = gn, 
                    GeneSpace        = GENE_SPACE, 
                    Threads          = Threads, 
                    MinPathSize      = 10, 
                    MaxPathSize      = 500 
                )
                PathwayAnlysisResult <- ldply(PathwayAnlysisResList, .id="path_id") %>% mutate(path_id=as.character(path_id)) 
                # scores -------------------------------------------------------            
                max_pvalue_score <- max(-log10(PathwayAnlysisResult[which(PathwayAnlysisResult$pvalue != 0 & !is.na(PathwayAnlysisResult$pvalue)), 'pvalue'])+2)
                PathwayAnlysisResult$score <- apply( PathwayAnlysisResult[,c("pvalue","enrich_factor")], 1, function(w)
                {
                    if( is.na(w[1])     ){ pvalueScore <- 0 
                    }else if( w[1] == 0 ){ pvalueScore <- max_pvalue_score
                    }else{                 pvalueScore <- -log10(w[1]) 
                    }
                    score <- sqrt( pvalueScore * w[2] )
                    score <- round(score, 2)
                    return(score)
                })
                # add pathway info ---------------------------------------------
                PathwayAnlysisResult <- data.frame(
                    PathwayAnlysisResult,
                    PathwayInfo[match(PathwayAnlysisResult$path_id, PathwayInfo$path_id), c("pathway","data_source","ext_id") ]
                ) %>% arrange( dplyr::desc(score) )
                # filter no gene-enriched pathawys -----------------------------
                PathwayAnlysisResult <- PathwayAnlysisResult %>% filter( enrich_genes >= 3 )
                #---------------------------------------------------------------
                return(PathwayAnlysisResult)
            }, mc.cores = length(PathwayAnalysisIndex) )
            # post processing --------------------------------------------------
            if( is.null( PathSubset ) ){ SubsetName <- "" }else{ SubsetName <- PathSubset }
            PATHWAY_ANALYSIS_RESULT <- data.frame()
            for( k in 1:length(PathwayAnalysisResultListSet) )
            {
                PATHWAY_ANALYSIS_RESULT <- rbind(
                    PATHWAY_ANALYSIS_RESULT,
                    data.frame(
                        seq_folder           = SEQ_FOLDER,
                        deg_analysis_id      = unique(DegAnalysisResult$analysis_id),
                        input_geneset        = names(PathwayAnalysisResultListSet)[k],
                        analysis_pathset_id  = PathwaySetID,
                        analysis_subset_name = SubsetName,
                        PathwayAnalysisResultListSet[[k]]                    
                    )
                )
            }
            message("|--->>> pathway analysis result table is created.")
        }else{
            PATHWAY_ANALYSIS_RESULT <- data.frame(seq_folder=unique(DegAnalysisResult$seq_folder), deg_analysis_id=unique(DegAnalysisResult$analysis_id), 
                input_geneset=NA, analysis_pathset_id=NA, analysis_subset_name=NA, path_id=NA, pathway_size=NA, pvalue=NA, enrich_factor=NA, enrich_genes=NA, gene_list=NA, score=NA, pathway=NA, data_source=NA, ext_id=NA
            )
        }
        #-----------------------------------------------------------------------
        return(PATHWAY_ANALYSIS_RESULT)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Load Pathway Analysis Result From Database for Plots /---------------------------------------#
#' @description load pathway analysis result based on analysis batch and analysis-id
#' @param SeqFolder    REQUIRED. analysis batch id ( = seq_folder )
#' @param AnalysisId   REQUIRED. analysis id. 
#' @param DataSource   data source. either 'db' or 'rds'. default = 'db'
#' @param RdsDir       RDS dir. default = /data/wts/RDS/RDS_PathwayAnalysis
#' @export 
    load.Pathway.Analysis.Profile <- function( SeqFolder=NULL, AnalysisId=NULL, DataSource='db', RdsDir="/data/wts/RDS/RDS_PathwayAnalysis" )
    {
        # seq-folder check -----------------------------------------------------
        if( is.null(SeqFolder) ){ stop("|---!!! no seq-folder-id. REQUIRED. please check again. STOPPED.") }
        # analysis id check ----------------------------------------------------
        if( is.null(AnalysisId) ){ stop("|---!!! no analysis-id. REQUIRED. please check again. STOPPED.") }
        # gene info ------------------------------------------------------------
        #if( is.null(GeneInfo) ){ stop("|---!!! no gene info table. REQUIRED. please check again. STOPPED.") }
        #-------------------------------------------------------------------------------------------
        message(sprintf("|---> %s - %s ", SeqFolder, AnalysisId))
        #-------------------------------------------------------------------------------------------
        if( DataSource == "db" )
        {
            # load data from database ------------------------------------------
            message("|---> load pathway analysis result data from database.")
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            ora_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM pathway_analysis_top200_result WHERE seq_folder = '%s' AND deg_analysis_id = '%s'", SeqFolder, AnalysisId))
            dbDisconnect(dbCon)
        }else{
            # load data from RDS -----------------------------------------------
            message("|---> load DEG analysis result data from RDS.")
            ora_res <- readRDS(sprintf("%s/%s.%s.pathway.analysis.top200.result.table.rds", RdsDir, SeqFolder, AnalysisId))
        }
        #-----------------------------------------------------------------------
        return(ora_res)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Make Network Node and Edge Data List of Top-Ranked Results /---------------------------------#
#' @description create nodes/edges network data of top-ranked pathway analysis results
#' @param PathwayAnalysisResult      REQUIRED. pathway analysis reuslt. 
#' @param MultiPathInvolveThreshold  pathway numbers cut-off of  multi-pathway involved genes. default = 3 
#' @param PathwaySetId               pathway set id
#' @param TopRank_N                  top-rank #N of pathway analysis reuslt to show on network.
#' @param UseMergedAnalysisResult    use 'ALL_DEGs' resutl instead of UP/DOWN genes result combine. default = FALSE
#' @param DegAnalysisResult          deg analysis result for up/down gene identification. REQUIRD WHEN UseMergedAnalysisResult = TRUE
#' @param GeneInfo                   REQUIRED. gene information. 
#' @eport 
    create.TopRank.Pathways.and.Genes.NetData <- function( PathwayAnalysisResult, MultiPathInvolveThreshold=3, PathwaySetId=NULL,
        TopRank_N=NULL, UseMergedAnalysisResult=FALSE, DegAnalysisResult=NULL, GeneInfo=NULL
    )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        # gene information -----------------------------------------------------
        if( is.null(GeneInfo) ){ stop("|---!!! No gene information table found. REQUIRED. please check again. STOPPED.") }
        # pathway set id check -------------------------------------------------
        if( is.null(PathwaySetId) ){ stop("|---!!! No pathway set id selected. REQUIRED. please check again. STOPPED.") }
        # top-rank analysis result ---------------------------------------------
        if( is.null(TopRank_N) )
        {
            stop("|---!!! TopRank #N is not specified. REQUIRED. please check again. STOPPED.")
        }else{
            if( UseMergedAnalysisResult )
            {
                TopRankRes <- PathwayAnalysisResult %>% dplyr::filter( analysis_pathset_id == PathwaySetId, input_geneset == "ALL_DEGs"  ) %>% 
                    arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N)
            }else{
                TopRankRes_up   <- PathwayAnalysisResult %>% dplyr::filter( analysis_pathset_id == PathwaySetId, input_geneset == "UP"   ) %>% 
                    arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N)
                TopRankRes_down <- PathwayAnalysisResult %>% dplyr::filter( analysis_pathset_id == PathwaySetId, input_geneset == "DOWN" ) %>% 
                    arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N)    
            }
        }
        # 
        if( UseMergedAnalysisResult )
        {
            if( is.null(DegAnalysisResult) )
            {
                stop("|---!!! No DEG analysis result found. Required if DegAnalysisResult = TRUE. please check again. STOPPED.")
            }else{
                genes_path <- lapply(with(TopRankRes, structure( get('gene_list'), names=get('pathway'))), function(gn) unlist(strsplit(gn, ";")))
                genes_path <- lapply( genes_path, function(gn) 
                {
                    gs <- GeneInfo[match(gn, GeneInfo$entrez), "symbol"] 
                    gs <- unique(gs[ !is.na(gs) ])
                    return(gs)
                })
                genes_counts <- table(unlist(genes_path))
                NET_GENES    <- names(genes_counts[genes_counts >= MultiPathInvolveThreshold ])
                deg_net      <- lapply( genes_path,   function(y) y[y %in% NET_GENES]   )
                deg_net      <- deg_net[ which(lapply( deg_net , length ) != 0 ) ]

                DEG_NET <- data.frame()
                for( k in 1:length(deg_net) )
                { DEG_NET <- rbind( DEG_NET, data.frame( pathway = names(deg_net)[k], gene = deg_net[[k]], direct = "" ) ) }

                DEG_NET$direct <- sapply( DEG_NET$gene, function(gn) 
                {
                    if(       DegAnalysisResult[which(DegAnalysisResult$symbol == gn), 'deg'] == "UP"   ){ drct <- "up_gene" 
                    }else if( DegAnalysisResult[which(DegAnalysisResult$symbol == gn), 'deg'] == "DOWN" ){ drct <- "down_gene" 
                    }else{ drct = "" }
                    return(drct)
                }) 
            }
        }else{
            # UP genes top-rank ----------------------------------------------------
            up_genes_path <- lapply(with(TopRankRes_up, structure( get('gene_list'), names=get('pathway'))), function(gn) unlist(strsplit(gn, ";")))
            up_genes_path <- lapply( up_genes_path, function(gn) 
            {
                gs <- GeneInfo[match(gn, GeneInfo$entrez), "symbol"] 
                gs <- unique(gs[ !is.na(gs) ])
                return(gs)
            })
            # DOWN genes top-rank --------------------------------------------------
            down_genes_path <- lapply(with(TopRankRes_down, structure( get('gene_list'), names=get('pathway'))), function(gn) unlist(strsplit(gn, ";")))
            down_genes_path <- lapply( down_genes_path, function(gn) 
            {
                gs <- GeneInfo[match(gn, GeneInfo$entrez), "symbol"] 
                gs <- unique(gs[ !is.na(gs) ])
                return(gs)
            })
            # UP/DOWN gene counts and threshold filtering --------------------------
            up_genes_counts   <- table(unlist(up_genes_path))
            NET_UP_GENES      <- names(up_genes_counts[up_genes_counts >= MultiPathInvolveThreshold ])
            down_genes_counts <- table(unlist(down_genes_path))
            NET_DOWN_GENES    <- names(down_genes_counts[down_genes_counts >= MultiPathInvolveThreshold ])
            # pathway filtering ----------------------------------------------------
            up_net   <- lapply( up_genes_path,   function(y) y[y %in% NET_UP_GENES]   )
            up_net   <- up_net[ which(lapply( up_net , length ) != 0 ) ]
            down_net <- lapply( down_genes_path, function(y) y[y %in% NET_DOWN_GENES] )
            down_net <- down_net[ which(lapply( down_net , length ) != 0 ) ]
            # up genes network -----------------------------------------------------
            UP_NET <- data.frame()
            if( length(up_net) > 0 ){
                for( k in 1:length(up_net) )
                { UP_NET <- rbind( UP_NET, data.frame( pathway = names(up_net)[k], gene = up_net[[k]], direct = "up_gene" ) ) }
            }
            # down genes network ---------------------------------------------------
            DOWN_NET <- data.frame()
            if( length(down_net) > 0 ){
                for( k in 1:length(down_net) )
                { DOWN_NET <- rbind( DOWN_NET, data.frame( pathway = names(down_net)[k], gene = down_net[[k]], direct = "down_gene" ) ) }
            }
            # network --------------------------------------------------------------
            DEG_NET <- rbind(UP_NET, DOWN_NET)
        }
        # nodes/edges ----------------------------------------------------------
        if( nrow(DEG_NET) > 0 ){
            NET_nodes <- rbind(
                data.frame(id=unique(DEG_NET$pathway), type='pathway'),
                unique( DEG_NET[,c('gene','direct')] %>% dplyr::rename( id=1, type=2) )
            )
            NET_edges <- unique(DEG_NET[,c('pathway','gene')] %>% dplyr::rename( from=1, to=2))
            # network data results -------------------------------------------------
            NET_RESULT <- list( node = NET_nodes, edge = NET_edges )
        }else{
            NET_RESULT <- data.frame()
        }
        return(NET_RESULT)
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| GSEA  MODULES |================================================================================================================================#
#---/ Run GSEA /-----------------------------------------------------------------------------------#
#' @description run GSEA (Gene Set Enrich Analysis) analysis
#' @param DegAnalysisResult deg analysis result table
#' @param AnalysisID        deg analysis id 
#' @param PathwaySetID      pathway-set-id. default = 'pw01' ( = consensuspathdb )
#' @param PathwaySetName    pathway-set-name. default = 'cpdb' ( = consensuspathdb )
#' @param PathSubset        pathway subset by data-source. (optional)
#' @param GeneIdentifier    gene identifier. either 'entrez' or 'symbol'. default = 'entrez'
#' @param RankMetric        input profile rank-metric. either 'score' or 'stat'. default = 'score'
#' @param Threads           cpu cores for parallel processing. default = 15.
#' @param MinPathSize       minimum pathway set genes size
#' @param MaxPathSize       maximum pathway set genes size
#' @export 
    run.GSEA.Analysis <- function( DegAnalysisResult, AnalysisID=NULL, PathwaySetID="pw01", PathwaySetName="cpdb", PathSubset=NULL, 
        GeneIdentifier="entrez", RankMetric="score", Threads=15, MinPathSize=10, MaxPathSize=500
    )
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("fgsea"))
        #-----------------------------------------------------------------------
        if( is.null(AnalysisID) ){ stop("|---!!! no analysis-id. STOPPED.")}
        # get pathway set ------------------------------------------------------
        PATHWAY_SET <- load.Pathway.Gene.List( PathwaySetID=PathwaySetID, PathwaySetName=PathwaySetName, Identifier=GeneIdentifier, DataSource="db", PathSubset=PathSubset )
        PathwayList <- PATHWAY_SET$pathwayGenes
        PathwayInfo <- PATHWAY_SET$pathwayInfo
        # gene space -----------------------------------------------------------
        if( GeneIdentifier == "symbol" )
        {
            PF_DegRes  <- DegAnalysisResult %>% filter( keep ) %>% filter( !is.na(symbol) ) %>% unique() 
            GENE_SPACE <- PF_DegRes %>% .$symbol %>% unique()
        }else{
            PF_DegRes  <- DegAnalysisResult %>% filter( keep ) %>% filter( !is.na(entrez) ) %>% unique()
            GENE_SPACE <- PF_DegRes %>% .$entrez %>% unique()
        }
        # input profile rank ---------------------------------------------------
        if( RankMetric == "stat" )
        {
            PF_DegRes      <- PF_DegRes %>% filter( stat != 0 ) %>% arrange(dplyr::desc(stat)) 
            GENE_EXPR_RANK <- with(PF_DegRes, structure(get('stat'), names=get(GeneIdentifier))) 
        }else{
            PF_DegRes      <- PF_DegRes %>% filter( score != 0 ) %>% arrange(dplyr::desc(score)) 
            GENE_EXPR_RANK <- with(PF_DegRes, structure(get('score'), names=get(GeneIdentifier)))
        }
        # pathway filtering ----------------------------------------------------
        PathListFilterd <- mclapply( PathwayList,     function(pw) pw[ pw %in% GENE_SPACE ], mc.cores=Threads )
        PathSizeCheck   <- mclapply( PathListFilterd, function(pw) ifelse( length(pw) < MinPathSize | length(pw) > MaxPathSize , FALSE, TRUE ), mc.cores=Threads )
        PF_PATHWAY_LIST <- PathListFilterd[ which( PathSizeCheck == TRUE ) ]
        # run GSEA -------------------------------------------------------------        
        GSEA_results = suppressWarnings(fgseaMultilevel(
            pathways    = PF_PATHWAY_LIST,
            stats       = GENE_EXPR_RANK/max(abs(GENE_EXPR_RANK)),
            minSize     = MinPathSize, 
            maxSize     = MaxPathSize,
            eps         = 1e-60,
            nPermSimple = 10000
        ))
        # leading edge genes ---------------------------------------------------
        GSEA_results$genes <- unlist(lapply( GSEA_results$leadingEdge, function(gn) paste(gn, collapse=";")) )
        GSEA_results$hits  <- unlist(lapply( GSEA_results$leadingEdge, length ) )
        # scores ---------------------------------------------------------------
        GSEA_results <- GSEA_results %>% filter( !is.na(pval) )
        if( length(which(GSEA_results$pval == 0)) > 0 )
        { max_pval_score <- max(-log10(GSEA_results[which(GSEA_results$pval != 0 ), "pval"])) }
        GSEA_results$score <- apply(GSEA_results[,c("pval","NES")], 1, function(w)
        {
            if( w[1] == 0 ){ pval_score <- max_pval_score }else{ pval_score <- -log10(w[1]) }
            score <- sqrt( pval_score * abs(w[2]) )
            if( w[2] < 0 ){ score <- score * (-1) }
            return(score)
        })
        GSEA_results <- GSEA_results %>% dplyr::rename( path_id=pathway )
        # final result table ---------------------------------------------------
        GSEA_RESULTS <- data.frame(
            seq_folder  = SEQ_FOLDER,
            analysis_id = AnalysisID,
            PathwayInfo[match(GSEA_results$path_id, PathwayInfo$path_id), c("pathway","data_source","ext_id")],
            GSEA_results[,c("path_id","pval","padj","ES","NES","score","log2err","hits","genes")]
        ) %>% as.data.frame() %>% arrange(dplyr::desc(abs(score))) 
        #-----------------------------------------------------------------------
        return(GSEA_RESULTS)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Load GSEA Result From Database for Plots /---------------------------------------------------#
#' @description load GSEA result based on analysis batch and analysis-id
#' @param SeqFolder    REQUIRED. analysis batch id ( = seq_folder )
#' @param AnalysisId   REQUIRED. analysis id. 
#' @param DataSource   data source. either 'db' or 'rds'. default = 'db'
#' @param RdsDir       RDS dir. default = /data/wts/RDS/RDS_PathwayAnalysis
#' @export 
    load.GSEA.Profile <- function( SeqFolder=NULL, AnalysisId=NULL, DataSource='db', RdsDir="/data/wts/RDS/RDS_GSEA" )
    {
        # seq-folder check -----------------------------------------------------
        if( is.null(SeqFolder) ){ stop("|---!!! no seq-folder-id. REQUIRED. please check again. STOPPED.") }
        # analysis id check ----------------------------------------------------
        if( is.null(AnalysisId) ){ stop("|---!!! no analysis-id. REQUIRED. please check again. STOPPED.") }
        # gene info ------------------------------------------------------------
        #if( is.null(GeneInfo) ){ stop("|---!!! no gene info table. REQUIRED. please check again. STOPPED.") }
        #-------------------------------------------------------------------------------------------
        message(sprintf("|---> %s - %s ", SeqFolder, AnalysisId))
        #-------------------------------------------------------------------------------------------
        if( DataSource == "db" )
        {
            # load data from database ------------------------------------------
            message("|---> load pathway analysis result data from database.")
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            ora_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM gsea_top200_result WHERE seq_folder = '%s' AND analysis_id = '%s'", SeqFolder, AnalysisId))
            dbDisconnect(dbCon)
        }else{
            # load data from RDS -----------------------------------------------
            message("|---> load DEG analysis result data from RDS.")
            ora_res <- readRDS(sprintf("%s/%s.%s.GSEA.top200.result.table.rds", RdsDir, SeqFolder, AnalysisId))
        }
        #-----------------------------------------------------------------------
        return(ora_res)
    }
#--------------------------------------------------------------------------------------------------#
#====================================================================================================================================================#

#===| GO ENRICH ANALYSIS MODULES |===================================================================================================================#
#---/ Get GO-Ancesotors using EBI API /------------------------------------------------------------#
#' @description get GO ancestor chart as graph using QuickGO API
#' @param InputGO  input GO term ID
#' @param root     end pooint GO term. either 'gobp', 'gocc', or 'gomf'. default = 'gobp'
#' @param GoTerms  REQUIRED. Gencurix pathwayset list containing root
#' @export
    Get.GO.Ancestor.Graph <- function( InputGO, root="gobp", GoTerms=NULL )
    {
        # input check ----------------------------------------------------------
        if( is.null(GoTerms) )
        { stop("|---!!! no pathwayset list containing root. REQUIRED. please check again. STOPPED.") }
        # packages -------------------------------------------------------------
        suppressPackageStartupMessages(library(httr))
        suppressPackageStartupMessages(library(jsonlite))
        suppressPackageStartupMessages(library(xml2))
        # root GO term ---------------------------------------------------------
        GOBP="GO%3A0008150"
        GOCC="GO%3A0005575"
        GOMF="GO%3A0003674"
        # QuickGO URL ----------------------------------------------------------
        API_Base_Graph1 = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/graph?startIds=GO%3A"
        API_Base_Graph2 = "&stopIds="
        # assign root GO term --------------------------------------------------
        if( root =="gocc" ){ ROOT <- GOCC }else if( root == "gomf" ){ ROOT <- GOMF }else{ ROOT <- GOBP }
        # GO term ID -----------------------------------------------------------
        INPUT_GO <- unlist(strsplit(InputGO, ":"))[2]
        # get data using API ---------------------------------------------------
        REQUEST_URL   <- paste0( API_Base_Graph1, INPUT_GO, API_Base_Graph2, ROOT )
        GO_GRAPH_JSON <- GET(REQUEST_URL, accept("application/json"))
        stop_for_status(GO_GRAPH_JSON)
        # convert data into node/edge list -------------------------------------
        GRAPH_LIST <- content(GO_GRAPH_JSON)
        # check data exist -----------------------------------------------------
        if( GRAPH_LIST$numberOfHits == 0 )
        {
            data_edges = data_nodes = data.frame()
        }else{
            data_edges <- as.data.frame(do.call(rbind, GRAPH_LIST$results[[1]]$edges))
            data_nodes <- as.data.frame(do.call(rbind, GRAPH_LIST$results[[1]]$vertices))
            for( i in 1:ncol(data_edges) ){ data_edges[,i] <- unlist(data_edges[,i]) }
            for( i in 1:ncol(data_nodes) ){ data_nodes[,i] <- unlist(data_nodes[,i]) }
            data_nodes <- data_nodes %>% filter( id %in% GoTerms )
            data_edges <- data_edges %>% filter( subject %in% GoTerms, object %in% GoTerms )
        }
        # final result ---------------------------------------------------------
        RESULT <- list( edge = data_edges, node = data_nodes )
        return(RESULT)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Create Network of Top-Enriched GO /----------------------------------------------------------#
#' @description make network data (node/edge) of top-ranked GO enrichment analysis result
#' @param GoPathwayAnalysisResult  GO enrich analysis result 
#' @param GoTerms                  GO term table 
#' @param TopRank_N                Top rank #N to create network
#' @param AnalysisGeneType         Input DEG type of GO enrich analysis  
#' @param GoRoot                   GO term root ( = ancestor end point )
#' @param IncludeAllTerms          shown on all terms on network (even not in top-ranked result list)
#' @parma Threads                  parallel cpu cores #N
#' @eport 
    make.TopRank.GO.To.Network <- function( GoPathwayAnalysisResult, GoTerms=NULL,
        TopRank_N=NULL, AnalysisGeneType="ALL_DEGs", GoRoot="gobp", IncludeAllTerms=FALSE, Threads=NULL
    )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("parallel"))
        # top-rank analysis result ---------------------------------------------
        if( is.null(TopRank_N) )
        {
            stop("|---!!! TopRank #N is not specified. REQUIRED. please check again. STOPPED.")
        }else{
            TopRankRes <- GoPathwayAnalysisResult %>% dplyr::filter( input_geneset == AnalysisGeneType ) %>% 
                arrange(dplyr::desc(score)) %>% dplyr::slice(1:TopRank_N)
        }
        if( nrow(TopRankRes) > 0 ){
            # go term id check -----------------------------------------------------
            if( all(is.na(TopRankRes$ext_id) | TopRankRes$ext_id == "") )
            {
                stop( "|---!!! all GO terms has no ID. please check analysis result. STOPPED." )
            }else{
                if( length(grep("^GO", TopRankRes$ext_id)) != nrow(TopRankRes) )
                {
                    if( IncludeAllTerms )
                    { 
                        stop( "|---!!! some terms has no ID and 'IncludeAllTerms'=TRUE. STOPPED." ) 
                    }else{
                        TopRankRes <- TopRankRes[grep("^GO", TopRankRes$ext_id), ]
                    }
                }
            }
            # parallel processing param --------------------------------------------
            if( is.null(Threads) ){ THREADS <- nrow(TopRankRes)/3 }else{ THREADS <- Threads }
            # get GO term ancestor tree data ---------------------------------------
            GoTreeNetworkData <- mclapply( TopRankRes$ext_id, function(GTID) Get.GO.Ancestor.Graph(InputGO=GTID, root=GoRoot, GoTerms=GoTerms), mc.cores=THREADS )
            # prepare node / edge tables -------------------------------------------
            EDGE_DATA <- ldply(lapply( GoTreeNetworkData, function(gn) gn$edge )) %>% unique() %>%
                mutate( 
                    analysis_id       = unique(GoPathwayAnalysisResult$deg_analysis_id),
                    seq_folder        = unique(GoPathwayAnalysisResult$seq_folder),
                    analysis_deg_type = AnalysisGeneType
                )
            NODE_DATA <- ldply(lapply( GoTreeNetworkData, function(gn) gn$node )) %>% unique() %>% 
                mutate( 
                    analysis_id       = unique(GoPathwayAnalysisResult$deg_analysis_id),
                    seq_folder        = unique(GoPathwayAnalysisResult$seq_folder),
                    analysis_deg_type = AnalysisGeneType
                )
            # final result ---------------------------------------------------------
            NETWORK_DATA <- list( node = NODE_DATA, edge = EDGE_DATA )
            #-----------------------------------------------------------------------
        }else{
            NETWORK_DATA <- data.frame()
        }
        return(NETWORK_DATA)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Make Network Node and Edge List of GO /------------------------------------------------------#
#' @description 
#' @param GoPathwayAnalysisResult  GO enrich analysis result 
#' @param GoTerms                  GO term table 
#' @param TopRank_N                Top rank #N to create network
#' @param UseAllDegResult          If TRUE, use 'ALL_DEGs' as DEG type. Otherwise, make UP and DOWN network component separately and merge
#' @param GoRoot                   GO term root ( = ancestor end point )
#' @param IncludeAllTerms          shown on all terms on network (even not in top-ranked result list)
#' @parma Threads                  parallel cpu cores #N
#' @export 
    create.TopRank.GO.Network.Data <- function( GoPathwayAnalysisResult, GoTerms=NULL, TopRank_N=10, UseAllDegResult=FALSE,
        GoRoot="gobp", IncludeAllTerms=FALSE, Threads=NULL
    )
    {
        suppressPackageStartupMessages(library("dplyr"))
        #-----------------------------------------------------------------------
        if( UseAllDegResult )
        {
            message("|---> use ALL_DEGs analysis result. This network cannot add profiles Up/Down regulation of genes")
            NetworkDataList <- make.TopRank.GO.To.Network(
                GoPathwayAnalysisResult = GoPathwayAnalysisResult,
                GoTerms                 = GoTerms,
                TopRank_N               = TopRank_N,
                AnalysisGeneType        = "ALL_DEGs",
                GoRoot                  = GoRoot,
                IncludeAllTerms         = IncludeAllTerms,
                Threads                 = Threads
            )
        }else{
            message("|---> Up/Down network will be created separately, then will be merged as one.")
            NetDataUp <- make.TopRank.GO.To.Network(
                GoPathwayAnalysisResult = GoPathwayAnalysisResult,
                GoTerms                 = GoTerms,
                TopRank_N               = TopRank_N,
                AnalysisGeneType        = "UP",
                GoRoot                  = GoRoot,
                IncludeAllTerms         = IncludeAllTerms,
                Threads                 = Threads
            )
            NetDataDown <- make.TopRank.GO.To.Network(
                GoPathwayAnalysisResult = GoPathwayAnalysisResult,
                GoTerms                 = GoTerms,
                TopRank_N               = TopRank_N,
                AnalysisGeneType        = "DOWN",
                GoRoot                  = GoRoot,
                IncludeAllTerms         = IncludeAllTerms,
                Threads                 = Threads
            )
            if( class(NetDataUp) != "data.frame" ){
                if( class(NetDataDown) != "data.frame" ){
                    NetworkDataList <- list(
                        node = rbind(
                            NetDataUp$node   %>% mutate( analysis_deg_type = "UP_DOWN" ),
                            NetDataDown$node %>% mutate( analysis_deg_type = "UP_DOWN" )
                        ) %>% unique(),
                        edge = rbind(
                            NetDataUp$edge   %>% mutate( analysis_deg_type = "UP_DOWN" ),
                            NetDataDown$edge %>% mutate( analysis_deg_type = "UP_DOWN" )
                        ) %>% unique()
                    )
                }else{
                    NetworkDataList <- NetDataUp
                }
            }else{
                if( class(NetDataDown) != "data.frame" ){
                    NetworkDataList <- NetDataDown
                }else{
                    NetworkDataList <- data.frame()
                }
            }
        }
        #-----------------------------------------------------------------------
        return(NetworkDataList)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Load GO Enrich Result from Database for Plots /----------------------------------------------#
#' @description load GO enrich analysis result based on analysis batch and analysis-id
#' @param SeqFolder    REQUIRED. analysis batch id ( = seq_folder )
#' @param AnalysisId   REQUIRED. analysis id. 
#' @param DataSource   data source. either 'db' or 'rds'. default = 'db'
#' @param RdsDir       RDS dir. default = /data/wts/RDS/RDS_PathwayAnalysis
#' @export 
    load.GO.Enrich.Analysis.Profile <- function( SeqFolder=NULL, AnalysisId=NULL, DataSource='db', RdsDir="/data/wts/RDS/RDS_GoEnrichAnalysis" )
    {
        # seq-folder check -----------------------------------------------------
        if( is.null(SeqFolder) ){ stop("|---!!! no seq-folder-id. REQUIRED. please check again. STOPPED.") }
        # analysis id check ----------------------------------------------------
        if( is.null(AnalysisId) ){ stop("|---!!! no analysis-id. REQUIRED. please check again. STOPPED.") }
        # gene info ------------------------------------------------------------
        #if( is.null(GeneInfo) ){ stop("|---!!! no gene info table. REQUIRED. please check again. STOPPED.") }
        #-------------------------------------------------------------------------------------------
        message(sprintf("|---> %s - %s ", SeqFolder, AnalysisId))
        #-------------------------------------------------------------------------------------------
        if( DataSource == "db" )
        {
            # load data from database ------------------------------------------
            message("|---> load pathway analysis result data from database.")
            dbCon <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            ora_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM go_enrich_analysis_top200_result WHERE seq_folder = '%s' AND deg_analysis_id = '%s'", SeqFolder, AnalysisId))
            dbDisconnect(dbCon)
        }else{
            # load data from RDS -----------------------------------------------
            message("|---> load DEG analysis result data from RDS.")
            ora_res <- readRDS(sprintf("%s/%s.%s.Go.enrichment.analysis.top200.result.table.rds", RdsDir, SeqFolder, AnalysisId))
        }
        #-----------------------------------------------------------------------
        return(ora_res)
    }
#--------------------------------------------------------------------------------------------------#
#---/ Create GOplot Object /-----------------------------------------------------------------------#
# #' @description create plotting object for GOplot functions. 
# #' @param GoEnrichProfile  REQUIRED. GO enrich analysis result 
# #' @param DegProfile       REQUIRED. DEG analysis result 
# #' @param GoTermTable      GO-Term information table. if NULL, pathway names of GoEnrichProfile will be used.
# #' @param GeneInfo         REQUIRED. gene information table to convert entrez into gene symbol.
# #' @param GoCategory       GO category. default = 'BP' 
# #' @export
#     create.GOplot.Object <- function( GoEnrichProfile=NULL, DegProfile=NULL, 
#         GoTermTable=NULL, GeneInfo=NULL, GoCategory="BP"
#     )
#     {
#         suppressPackageStartupMessages(library("GOplot"))
#         suppressPackageStartupMessages(library("Hmisc"))
#         #-----------------------------------------------------------------------
#         if( is.null(GoEnrichProfile) ){ stop("|---!!! No GO enrich analysis result found. REQUIRED. please check again. STOPPED.")}
#         if( is.null(DegProfile)      ){ stop("|---!!! No GO enrich analysis result found. REQUIRED. please check again. STOPPED.")}
#         if( is.null(GeneInfo)        ){ stop("|---!!! No GO enrich analysis result found. REQUIRED. please check again. STOPPED.")}
#         # terms data -----------------------------------------------------------
#         GoEnrichProfile$Genes <- sapply( GoEnrichProfile$gene_list, function(gl) 
#         {
#             symbol_list <- GeneInfo[match(unlist(strsplit(gl, ";")), GeneInfo$entrez), "symbol"]
#             symbol_list <- symbol_list[ !is.na(symbol_list) ]
#             symbols <- paste(symbol_list, collapse=", ")
#             return(symbols)
#         })
#         if( is.null(GoTermTable) )
#         {
#             GoEnrichProfile$Term <- tolower(GoEnrichProfile$pathway)
#             GoEnrichProfile$Term <- gsub("^gobp_", "", GoEnrichProfile$Term)
#             GoEnrichProfile$Term <- gsub("^gocc_", "", GoEnrichProfile$Term)
#             GoEnrichProfile$Term <- gsub("^gomf_", "", GoEnrichProfile$Term)
#             GoEnrichProfile$Term <- gsub("_", "", tolower(GoEnrichProfile$Term))
#         }else{
#             GoEnrichProfile$Term <- GoTermTable[match(GoEnrichProfile$ext_id, GoTermTable$GO_TermID), "GO_Name"]
#         }
#         # genes data -----------------------------------------------------------
#         DegProfile <- DegProfile %>% filter( deg != "" )
#         # maek object ----------------------------------------------------------
#         goplot_obj <- suppressWarnings(circle_dat(
#             terms  = data.frame(
#                 Category = GoCategory,
#                 ID       = GoEnrichProfile$ext_id,
#                 Term     = GoEnrichProfile$Term,
#                 Genes    = GoEnrichProfile$Genes,
#                 adj_pval = GoEnrichProfile$pvalue
#             ),
#             genes = data.frame(
#                 ID = DegProfile$symbol,
#                 logFC = DegProfile$log2_foldchange,
#                 adj.P.val = DegProfile$pvalue
#             )
#         ))
#         #-----------------------------------------------------------------------
#         return(goplot_obj)
#     }
#--------------------------------------------------------------------------------------------------#


