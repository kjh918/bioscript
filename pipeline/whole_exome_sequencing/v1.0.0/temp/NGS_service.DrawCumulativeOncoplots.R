

#---| PACKAGES |---------------------------------------------------------------#
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("reshape2"))
    suppressPackageStartupMessages(library("scales"))
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages(library("cowplot"))
    suppressPackageStartupMessages(library("magick"))
    suppressPackageStartupMessages(library("igraph"))
    suppressPackageStartupMessages(library("ComplexHeatmap"))
    suppressPackageStartupMessages(library("circlize"))
    suppressPackageStartupMessages(library("lubridate"))
    suppressPackageStartupMessages(library("Cairo"))
    options(bitmapType='cairo')
    grDevices::X11.options(type='cairo')
    options(device='x11')
#------------------------------------------------------------------------------#

#---| ARGUMENTS |--------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"), action="store", default=NA, type="character", help="BASE_DIR"),
        make_option(c("--SEQ_FOLDER"), action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--EXCLUDE_SEQ_ID"), action="store", default=NA, type="character", help="Exclude SeqIDs"),
        make_option(c("--CLIENT_ID"), action="store", default=NA, type="character", help="Client ID")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR     <- ARGS$BASE_DIR
    SEQ_FOLDER   <- ARGS$SEQ_FOLDER
    ExcludeSeqID <- ARGS$EXCLUDE_SEQ_ID
    CLIENT_ID    <- ARGS$CLIENT_ID
#------------------------------------------------------------------------------#

#---| DRAW ONCOPLOT FUNCTION |-------------------------------------------------#
#' @param targetGeneCountTable 'oncoPlotTargetGenes' function result table
#' @param SampleInfo           sample info table (from DB) 
#' @param SampleOrderByIdGroup sort by sample group id
#' @param includeAllSamples    include all samples even have no variants
#' @export
    drawOncoPlotsExtra <- function( 
        targetGeneCountTable, 
        SampleInfo, 
        SampleOrderByIdGroup = TRUE, 
        includeAllSamples    = FALSE 
    )
    {
        varColors <- c( 
            "Missense_Mutation"      = "#4DAF4A", 
            "Nonsense_Mutation"      = "#E41A1C", 
            "Frame_Shift_Ins"        = "#FF7F00", 
            "Frame_Shift_Del"        = "#8DD3C7",
            "In_Frame_Ins"           = "#984EA3", 
            "In_Frame_Del"           = "#377EB8", 
            "Nonstop_Mutation"       = "#A65628",
            "Translation_Start_Site" = "#FFFF33", 
            "Splice_Site"            = "#F781BF",
            "Multi_Variants"         = "#4c4c4c" 
        )
        
        if( length(unique(SampleInfo$sample_info)) >= 3 )
        {
            tissueColors        <- RColorBrewer::brewer.pal( length(unique(SampleInfo$sample_info)), "Dark2")
            names(tissueColors) <- unique(SampleInfo$sample_info)
        }else{
            tissueColors        <- c("#1B9E77","#D95F02")[1:length(unique(SampleInfo$sample_info))] 
            names(tissueColors) <- unique(SampleInfo$sample_info)
        }
    
        if( SampleOrderByIdGroup )
        {
            uniqSampleGroup <- SampleInfo %>% filter( sample_group != "" ) %>% .$sample_group %>% unique()
    
            if( length(uniqSampleGroup) > 8 )
            {
                groupColors <- rainbow(length(uniqSampleGroup))
            }else if( length(uniqSampleGroup) >= 3 & length(uniqSampleGroup) <= 8 ){
                groupColors <- RColorBrewer::brewer.pal(length(uniqSampleGroup), "Set2")
            }else if( length(uniqSampleGroup) == 0 ){
                groupColors <- "#BC80BD"
            }else{
                groupColors <- c("#BC80BD","#FFED6F")[1:length(uniqSampleGroup)]
            }

            if( length(uniqSampleGroup) == 0 )
            {
                names(groupColors) <- ""
                SampleOrderByIdGroup = FALSE
            }else{
                names(groupColors) <- uniqSampleGroup
            }
            
        }
        
        annotR_data <- targetGeneCountTable %>% 
            group_by(HGNC_SYMBOL, mutations) %>% 
            summarise( varCount = length(seq_id) ) %>% 
            data.frame() %>% 
            reshape(., idvar="HGNC_SYMBOL", timevar="mutations", direction="wide")
    
        annotR_data[is.na(annotR_data)] <- 0 
        colnames(annotR_data) <- gsub("^varCount.", "", colnames(annotR_data))
        rownames(annotR_data) <- annotR_data[,1]
        annotR_data <- annotR_data[, c("HGNC_SYMBOL", names(varColors)[which(names(varColors) %in% colnames(annotR_data))]) ]
    
        annotC_data <- targetGeneCountTable %>% 
            group_by(seq_id, mutations) %>% 
            summarise(length(HGNC_SYMBOL)) %>% 
            data.frame() %>% 
            reshape(., idvar="seq_id", timevar="mutations", direction="wide")
    
        colnames(annotC_data) <- gsub("^length.HGNC_SYMBOL..", "", colnames(annotC_data))
        rownames(annotC_data) <- annotC_data[,1]
        annotC_data[is.na(annotC_data)] <- 0
        annotC_data <- annotC_data[, c("seq_id", names(varColors)[which(names(varColors) %in% colnames(annotC_data))]) ]
    
        oncoPlotMatrix <- targetGeneCountTable %>% data.frame() %>% reshape(., idvar="HGNC_SYMBOL", timevar="seq_id", direction="wide")
    
        rownames(oncoPlotMatrix) <- oncoPlotMatrix[,1]
        colnames(oncoPlotMatrix) <- gsub("^mutations.", "", colnames(oncoPlotMatrix))
    
        if( SampleOrderByIdGroup )
        {
            geneOrder <- apply(oncoPlotMatrix[,-1], 1, function(y) length(y[!is.na(y)])) %>% sort() %>% rev()
    
            topOneGeneTissueOrders <- targetGeneCountTable %>% 
                filter( HGNC_SYMBOL == (apply(oncoPlotMatrix[,-1], 1, function(y) length(y[!is.na(y)])) %>% 
                sort() %>% rev() %>% head(1) %>% names())) %>%
                left_join(., SampleInfo, by='seq_id') %>% .$sample_info %>% table() %>% sort() %>% rev() 
    
            colOrders <- data.frame()
    
            for( n in 1:length(topOneGeneTissueOrders) )
            {
                colOrders <- rbind( 
                    colOrders,
                    rbind(
                        SampleInfo %>% filter( sample_info == names(topOneGeneTissueOrders)[n], sample_group %nin% uniqSampleGroup ),
                        SampleInfo %>% filter( sample_info == names(topOneGeneTissueOrders)[n], sample_group %in% uniqSampleGroup ) %>% arrange(sample_group, desc(sample_origin))
                    )
                )
            }
    
            if( !includeAllSamples ){ colOrders <- colOrders %>% filter(seq_id %in% targetGeneCountTable$seq_id ) }
    
            colOrders <- colOrders %>% dplyr::rename( id=seq_id )
        }else{
    
            geneOrder <- apply(oncoPlotMatrix[,-1], 1, function(y) length(y[!is.na(y)])) %>% sort() %>% rev()
    
            colOrderMat <- oncoPlotMatrix
            for( i in 2:ncol(colOrderMat))  
            {
                colOrderMat[!is.na(colOrderMat[,i]),i] <- as.numeric(1)
                colOrderMat[ is.na(colOrderMat[,i]),i] <- as.numeric(0)
            }
    
            colOrders <- t(colOrderMat[names(geneOrder),-1])
            colOrders <- data.frame(id=rownames(colOrders), colOrders, check.names=F)
            colOrders <- colOrders %>% arrange(across(names(geneOrder), desc))
        }
    
        if( includeAllSamples )
        {
            noVarSeqIDs <- colOrders$id[ colOrders$id %nin% colnames(oncoPlotMatrix) ]
    
            if( length(noVarSeqIDs) > 0 )
            {
                noVarMat           <- matrix( rep(NA, length(noVarSeqIDs)*nrow(oncoPlotMatrix)), ncol=length(noVarSeqIDs) )
                colnames(noVarMat) <- noVarSeqIDs
                rownames(noVarMat) <- rownames(oncoPlotMatrix)
                oncoPlotMatrix     <- cbind(oncoPlotMatrix, noVarMat)
    
                add_annotC <- data.frame(noVarSeqIDs, matrix(rep(0,length(noVarSeqIDs)*(ncol(annotC_data)-1)), nrow=length(noVarSeqIDs)))
                colnames(add_annotC) <- colnames(annotC_data)
    
                annotC_data           <- rbind( annotC_data, add_annotC )
                rownames(annotC_data) <- annotC_data$seq_id
            }
        }
    
        sampleFreq <- targetGeneCountTable %>% 
            group_by(HGNC_SYMBOL) %>% 
            summarise( sampleFreq = round(length(seq_id)/(ncol(oncoPlotMatrix)-1) *100, 0) ) %>%
            mutate(sampleFreq = gsub("$", "%", sampleFreq))
    
        annotR <- rowAnnotation(
            varRatio   = anno_text( sampleFreq[match(names(geneOrder), sampleFreq$HGNC_SYMBOL), ] %>% .$sampleFreq , gp=gpar(fontsize=9) ),
            SAMPLES    = anno_barplot( 
                annotR_data[names(geneOrder),-1], 
                border     = FALSE, 
                extend     = 0.5, 
                width      = unit(3, "cm"), 
                axis_param = list( labels_rot = 90 ),
                gp         = gpar( fill = varColors[colnames(annotR_data)[2:ncol(annotR_data)]], col =  varColors[colnames(annotR_data)[2:ncol(annotR_data)]] )
            ),
            annotation_name_gp = gpar(fontsize = 7), gap = unit(2, "mm")
        )        
    
        annotC <- HeatmapAnnotation(
            VARIANTS = anno_barplot( 
                annotC_data[colOrders$id,-1], 
                border = FALSE , 
                height = unit(1.5, "cm"), 
                extend = 0.5, 
                axis_param = list( labels_rot = -90 ),
                gp = gpar( fill = varColors[colnames(annotC_data)[2:ncol(annotC_data)]], col  = varColors[colnames(annotC_data)[2:ncol(annotC_data)]] )
            ),
            TISSUE = SampleInfo[match(colOrders$id, SampleInfo$seq_id), "sample_info"  ], 
            GROUP  = SampleInfo[match(colOrders$id, SampleInfo$seq_id), "sample_group" ],
            ORIGIN = SampleInfo[match(colOrders$id, SampleInfo$seq_id), "sample_origin"],
            col    = list( TISSUE = tissueColors, GROUP = groupColors, ORIGIN = c("TS" = "#1F78B4", "ORG" = "#A6CEE3") ),
            gap    = unit(1.2, "mm"),
            annotation_name_gp      = gpar(fontsize = 7), 
            annotation_name_side    = 'left', 
            annotation_legend_param = list( title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=7) )        
        )
    
        inputMat <- as.matrix(oncoPlotMatrix[names(geneOrder), colOrders$id])
    
        if( ncol(inputMat) <= 40 ){ cellWidth  = 4.5 ; colFontsize = 8 }else{ cellWidth = 3.2 ; colFontsize = 7 }
        if( nrow(inputMat) >= 35 ){ cellHeight = 3 ; rowFontsize = 6 }else{ cellHeight = 4 ; rowFontsize = 9 } 
    
        oncoPlotHeatMap <- Heatmap( inputMat, 
            cluster_rows         = FALSE, 
            cluster_columns      = FALSE, 
            col                  = varColors, 
            na_col               = "#f9f9f9",
            border_gp            = gpar(col = "#727272", lwd =0.5),
            rect_gp              = gpar(col = "#727272", lwd = 0.3),
            row_names_gp         = gpar(fontsize = rowFontsize), 
            column_names_gp      = gpar(fontsize = colFontsize), 
            row_names_side       = "left",
            width                = ncol(inputMat)*unit(cellWidth, 'mm'), 
            height               = nrow(inputMat)*unit(cellHeight, 'mm'),
            column_labels        = SampleInfo[match(colOrders$id, SampleInfo$seq_id), "sample_name"], 
            column_names_rot     = 60,
            heatmap_legend_param = list( 
                title        = 'VARIANT CLASS', 
                title_gp     = gpar(fontsize=10),
                labels_gp    = gpar(fontsize=8), 
                legend_width = unit(70,"mm") 
            ),
            right_annotation     = annotR, 
            top_annotation       = annotC
        )

        return(oncoPlotHeatMap)
    }
#------------------------------------------------------------------------------#


#---| DATA |-------------------------------------------------------------------#
    if( !is.na(ExcludeSeqID) ){ ExcludeSeqID = unlist(strsplit(ExcludeSeqID, ",")) }
    ###
    dbCon     = dbConnect(dbDriver("MySQL"), host='192.168.0.34', port=3306, user='gcx', password='gencurix!!', db = 'gcx_ngs_service' )
    sinfo     = dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    order_id  = dbGetQuery(dbCon, sprintf("SELECT ngs_order_id FROM seqfolder_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    var_pf    = dbGetQuery(dbCon, sprintf("SELECT * FROM gcx_wes.variants_somatic WHERE seq_folder = '%s'", SEQ_FOLDER))
    dbDisconnect(dbCon)
    ###
    sinfo[is.na(sinfo)] <- ""
    ###
    if( all(!is.na(ExcludeSeqID)) )
    { 
        varSinfo <- sinfo   %>% filter( seq_id %nin% ExcludeSeqID )
        var_pf   <- var_pf  %>% filter( seq_id %nin% ExcludeSeqID )
    }else{
        varSinfo <- sinfo
    }
    ###
    load("/storage/home/kangsm/shinyWeb/resources/PresetGenes.Rdata")
    source("/storage/home/kangsm/shinyWeb/resources/Fun_QC_Report.R")
    source("/storage/home/kangsm/shinyWeb/resources/Fun_VARIANTS.R")
    ###
    analysisVisGroupTag <- c(
        "TS_N"  = "variant call : matched normal TISSUE",
        "ORG_N" = "variant call : matched normal ORGANOID",
        "Tonly" = "variant call : Tumor-only",
        "TonlyGermlineFiltered" = "variant call : Tumor-only Germline Filtered"
    )
    ###
    analysis_profiles_raw <- var_pf %>% select(c("variant_call_mode","seq_id")) %>% unique() %>%
        left_join(., varSinfo[,c("seq_id","sample_info","sample_tissue")], by="seq_id") 
    analysis_profiles <- analysis_profiles_raw %>% 
        group_by(variant_call_mode, sample_info, sample_tissue) %>% reframe( analysis=length(seq_id) )
    ###
    cancerCodeTissue <- lapply(indexing(varSinfo$sample_info), function(cc) analysis_profiles %>% filter( sample_info == cc ) %>% .$sample_tissue %>% unique() )
    cancerTypePresetList <- lapply(cancerCodeTissue, function(st) 
    {
        genelist = ifelse( st %in% names(c11_preset_genes), st, "pan_cancers" )
        genelist = c(genelist, "mutation_panel")
        return(genelist)
    })
#------------------------------------------------------------------------------#

#---| DRAW ONCOPLOTS |---------------------------------------------------------#
    QUERY = sprintf("SELECT variants_somatic.variant_call_mode,COUNT(DISTINCT(variants_somatic.seq_id)),seqid_info.sample_info FROM variants_somatic INNER JOIN gcx_ngs_service.seqid_info ON variants_somatic.seq_id = seqid_info.seq_id WHERE seqid_info.sample_info in (%s) AND seqid_info.panel = 'WES' AND seqid_info.client_facility_id = '%s' GROUP BY variant_call_mode, sample_info;", paste(paste0("'", names(cancerCodeTissue) , "'"), collapse=","), CLIENT_ID )
    ###
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', port=3306, user='gcx', password='gencurix!!', db = 'gcx_wes' )
    cumSinfo_raw = dbGetQuery(dbCon, QUERY)
    cumSinfo_filtered = cumSinfo_raw %>% dplyr::rename(counts=2) %>% filter( counts >= 10 ) 
    dbDisconnect(dbCon)
    ###
    if( nrow(cumSinfo_filtered) > 0 )
    {
        for( i in 1:nrow(cumSinfo_filtered) )
        {
            cumVarData = CUM_ONCOPLOT_DATASETS2( 
                TISSUE_INFO       = cumSinfo_filtered$sample_info[i], 
                VARIANT_CALL_MODE = cumSinfo_filtered$variant_call_mode[i],
                CLIENT_ID         = CLIENT_ID
            )
            CUM_PRESETS = cancerTypePresetList[[ cumSinfo_filtered$sample_info[i] ]]
            CUM_MODE    = analysisVisGroupTag[ cumSinfo_filtered$variant_call_mode[i] ]
            #------------------------------------------------------------------#
            matchedTissue  = cancerCodeTissue[[ cumSinfo_filtered$sample_info[i] ]]
            matchedSamples = cumVarData$cumVarPF$seq_id %>% unique()
            OCP_TOP30 <- drawOncoPlotsExtra( 
                targetGeneCountTable = oncoPlotTargetGenes( 
                    selectCancerType = matchedTissue, 
                    VariantsMAF      = cumVarData$cumVarPF, 
                    topGenes         = 30, 
                    SampleInfo       = cumVarData$cumSinfo %>% filter( seq_id %in% matchedSamples )
                ),  
                SampleInfo           = cumVarData$cumSinfo %>% filter( seq_id %in% matchedSamples ), 
                SampleOrderByIdGroup = TRUE, 
                includeAllSamples    = TRUE
            )
            ##
            OCP_PRESET1 <- drawOncoPlotsExtra( 
                targetGeneCountTable = oncoPlotTargetGenes( 
                    selectCancerType = matchedTissue, 
                    VariantsMAF      = cumVarData$cumVarPF, 
                    geneSets         = c11_preset_genes[[ CUM_PRESETS[1] ]], 
                    SampleInfo       = cumVarData$cumSinfo %>% filter( seq_id %in% matchedSamples )
                ),  
                SampleInfo           = cumVarData$cumSinfo %>% filter( seq_id %in% matchedSamples ), 
                SampleOrderByIdGroup = TRUE, 
                includeAllSamples    = TRUE
            )
            ##
            OCP_PRESET2 <- drawOncoPlotsExtra( 
                targetGeneCountTable = oncoPlotTargetGenes( 
                    selectCancerType = matchedTissue, 
                    VariantsMAF      = cumVarData$cumVarPF, 
                    geneSets         = c11_preset_genes[[ CUM_PRESETS[2] ]], 
                    SampleInfo       = cumVarData$cumSinfo %>% filter( seq_id %in% matchedSamples )
                ),  
                SampleInfo           = cumVarData$cumSinfo %>% filter( seq_id %in% matchedSamples ), 
                SampleOrderByIdGroup = TRUE, 
                includeAllSamples    = TRUE
            )
            #------------------------------------------------------------------#
            
            cumulative_sample_oncoplot_dir = sprintf("%s/%s/meta/cumulative_sample_oncoplots", BASE_DIR, SEQ_FOLDER)
            if( !dir.exists(cumulative_sample_oncoplot_dir) ){ system(sprintf("mkdir -p %s", cumulative_sample_oncoplot_dir)) }

            plotName1 = sprintf("%s/%s.%s.top30.genes.%s.cumulative.samples.oncoplots", cumulative_sample_oncoplot_dir, order_id$ngs_order_id, cumSinfo_filtered$sample_info[i], cumSinfo_filtered$variant_call_mode[i])
            
            plotName2 = sprintf("%s/%s.%s.preset.%s.%s.cumulative.samples.oncoplots", cumulative_sample_oncoplot_dir, order_id$ngs_order_id, cumSinfo_filtered$sample_info[i], CUM_PRESETS[1], cumSinfo_filtered$variant_call_mode[i])
            
            plotName3 = sprintf("%s/%s.%s.preset.%s.%s.cumulative.samples.oncoplots", cumulative_sample_oncoplot_dir, order_id$ngs_order_id, cumSinfo_filtered$sample_info[i], CUM_PRESETS[2], cumSinfo_filtered$variant_call_mode[i])


            png(sprintf("%s_raw.png",plotName1), width=36, height=27, units="in", res=200 )
            draw(OCP_TOP30)
            dev.off()

            png(sprintf("%s_raw.png",plotName2), width=36, height=27, units="in", res=200 )
            draw(OCP_PRESET1)
            dev.off()

            png(sprintf("%s_raw.png",plotName3), width=36, height=27, units="in", res=200 )
            draw(OCP_PRESET2)
            dev.off()

            system(sprintf("convert %s_raw.png -trim -bordercolor white -border 20x20 %s.png", plotName1, plotName1))
            system(sprintf("convert %s_raw.png -trim -bordercolor white -border 20x20 %s.png", plotName2, plotName2))
            system(sprintf("convert %s_raw.png -trim -bordercolor white -border 20x20 %s.png", plotName3, plotName3))

        }
    }
#------------------------------------------------------------------------------#