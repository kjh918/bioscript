
#---| Pre-Defined Values |-----------------------------------------------------#
    statColumnNames <- factor( 
        c(  'seq_id'            = "seq_id"         , 'sample_name'            = "SampleName"           ,
            'Missense_Mutation' = "Missense"       , 'Nonsense_Mutation'      = "Nonsense"             , 
            'Frame_Shift_Ins'   = "FrameShift INS" , 'Frame_Shift_Del'        = "FrameShift DEL"       ,  
            'In_Frame_Ins'      = "InFrame INS"    , 'In_Frame_Del'           = "InFrame DEL"          ,  
            'Splice_Site'       = "SpliceSite"     , 'Translation_Start_Site' = "Translation StartSite",
            'Nonstop_Mutation'  = "NonStop"        , "3'Flank"                = "3-Flank"              ,  
            "3'UTR"             = "3-UTR"          , "5'Flank"                = "5-Flank"              , 
            "5'UTR"             = "5-UTR"          , "IGR"                    = "IGR"                  ,   
            "Intron"            = "Intron"         , "RNA"                    = "RNA"                  ,  
            "Silent"            = "Silent"         , "Splice Region"          = "Splice Region"        ,  
            "Targeted Region"   = "Targeted Region"
        ),
        levels = c( 
            "SampleName",     "Missense",    "Nonsense",    "NonStop",       "FrameShift INS", 
            "FrameShift DEL", "InFrame INS", "InFrame DEL", "SpliceSite",    "Translation StartSite",
            "Silent",         "Intron",      "IGR",         "3-UTR",         "5-UTR",
            "3-Flank",        "5-Flank",     "RNA",         "Splice Region", "Targeted Region",
            "seq_id" 
        )
    )    
    statColumnNames <- statColumnNames[order(statColumnNames)]
#------------------------------------------------------------------------------#

#---| VARIANTS STATS TABLE RE-FORMATTING |-------------------------------------#
#' @param VAR_STATS     variant class summary table ( from DB )
#' @param VARIANT_GROUP variant group-type for selection 
#' @param SAMPLE_INFO   sample information ( from DB )
#' @export
    VARIANT_STATS_REFORMAT <- function( VAR_STATS, VARIANT_GROUP, SAMPLE_INFO )
    {
        require(reshape2)
        require(dplyr)
        varStat_rf = VAR_STATS %>%
            filter( variant_group == VARIANT_GROUP ) %>% 
            mutate( SampleName = SAMPLE_INFO[match(seq_id, SAMPLE_INFO$seq_id), "sample_name"] ) %>%
            group_by( SampleName, Variant_Classification ) %>% 
            summarise( Variants = sum(stats) ) %>% data.frame() %>%
            reshape(., idvar = "SampleName", timevar = "Variant_Classification", direction = "wide" )

        if( VARIANT_GROUP %in% c("non_synonymous","filtered_somatic_only") )
        {
            ColumnList <- c( "Variants.Missense_Mutation","Variants.Nonsense_Mutation","Variants.Nonstop_Mutation", 
                             "Variants.Frame_Shift_Ins",  "Variants.Frame_Shift_Del",  "Variants.In_Frame_Ins", 
                             "Variants.In_Frame_Del",     "Variants.Splice_Site",      "Variants.Translation_Start_Site" )
            missingCols <- ColumnList[ ColumnList %nin% colnames(varStat_rf) ]
            
            if( length(missingCols) > 0 ) 
            {
                extraStats <- as.data.frame(matrix( rep( NA, length(missingCols)*nrow(varStat_rf) ), ncol=length(missingCols) ))
                colnames(extraStats) <- missingCols
                varStat_rf <- cbind(varStat_rf, extraStats)
            }   

            varStat_rf = varStat_rf %>%
                dplyr::rename( 
                    'Missense'              = Variants.Missense_Mutation, 
                    'Nonsense'              = Variants.Nonsense_Mutation, 
                    'NonStop'               = Variants.Nonstop_Mutation, 
                    'FrameShift INS'        = Variants.Frame_Shift_Ins, 
                    'FrameShift DEL'        = Variants.Frame_Shift_Del, 
                    'InFrame INS'           = Variants.In_Frame_Ins, 
                    'InFrame DEL'           = Variants.In_Frame_Del, 
                    'SpliceSite'            = Variants.Splice_Site, 
                    'Translation StartSite' = Variants.Translation_Start_Site
                ) %>%
                dplyr::select( c(
                    "SampleName","Missense","Nonsense","NonStop",
                    "FrameShift INS","FrameShift DEL","InFrame INS","InFrame DEL",
                    "SpliceSite","Translation StartSite"
                ) )
        }else{
            colnames(varStat_rf) = gsub("Variants.", ""  , colnames(varStat_rf))
            colnames(varStat_rf) = gsub("_"        , " " , colnames(varStat_rf))
        }
        varStat_rf[is.na(varStat_rf)] = 0 
        return(varStat_rf)
    }
#------------------------------------------------------------------------------#

#---| VARIANT-STAT-TABLE : NON-SYNONYMOUS |----------------------------------------#
#' @param VAR_STAT_RF         re-formed variant class summary table 
#' @param SAMPLE_ORDER        sample orders
#' @param SAMPLE_INFO         sample information ( from DB )
#' @param SAMPLE_ORDER        sample orders
#' @param VARIANT_CALL_METHOD variant calling method. TS_N, ORG_N, or Tonly.
#' @export
    VAR_SUM_TABLE_NonSynonymous <- function( VAR_STAT_RF, SAMPLE_ORDER, SAMPLE_INFO, VARIANT_CALL_METHOD )
    {
        VarSumNonSyn <- VAR_STAT_RF[match(SAMPLE_ORDER, VAR_STAT_RF$SampleName), ] %>%
            filter( !is.na(SampleName) )
        
        GrIdx <- SAMPLE_INFO %>% filter( sample_name %in% VarSumNonSyn$SampleName ) %>% .$sample_group %>% indexing()
        SampleGroupRows <- ldply( GrIdx, function(s) 
        {
            vsinfo     <- SAMPLE_INFO[match(VarSumNonSyn$SampleName, SAMPLE_INFO$sample_name), ]       
            sinfo.rows <- which(vsinfo$sample_group == s)
            return( data.frame(group = s, start = head(sinfo.rows, 1), end = tail(sinfo.rows, 1)) )
        })
        
        if( VARIANT_CALL_METHOD == "TS_N"  ){ VCM = "Matched Normal (Tissue)"   }
        if( VARIANT_CALL_METHOD == "ORG_N" ){ VCM = "Matched Normal (Organoid)" }
        if( VARIANT_CALL_METHOD == "Tonly" ){ VCM = "Tumor-only" }
        if( VARIANT_CALL_METHOD == "TonlyGermlineFiltered" ){ VCM = "Tumor-only and Germline Filtered" }

        VarSumNonSynTable <- VarSumNonSyn %>%
            kable( align=c("l",rep("c", 9)), row.names=F, escape=FALSE, 
                caption = sprintf("Non-Synonymous Variants : %s", VCM ) 
            ) %>%
            kable_classic( "striped", full_width=T, font_size=12, position='center', html_font='Roboto Condensed') %>%
            row_spec(0, align='c', bold=T, background="#D0E4A4") %>%
            column_spec(1, width="6cm")

        if( nrow(SampleGroupRows) > 0 )
        {
            for( n in 1:nrow(SampleGroupRows) )
            { 
                VarSumNonSynTable <- VarSumNonSynTable %>% 
                    group_rows(group_label = SampleGroupRows$group[n], SampleGroupRows$start[n], SampleGroupRows$end[n] ) 
            }
        }
        
        return(VarSumNonSynTable)
    }
#------------------------------------------------------------------------------#

#---| VARIANT-STAT-TABLE : SYNONYMOUS |----------------------------------------#
#' @param VAR_STAT_RF         re-formed variant class summary table 
#' @param SAMPLE_INFO         sample information ( from DB )
#' @param STAT_COLNAMES       pre-defined ordered column names ( statColumnNames )
#' @param SAMPLE_ORDER        sample orders
#' @param VARIANT_CALL_METHOD variant calling method. TS_N, ORG_N, or Tonly.
#' @export
    VAR_SUM_TABLE_Synonymous <- function( VAR_STAT_RF, STAT_COLNAMES, SAMPLE_ORDER, SAMPLE_INFO, VARIANT_CALL_METHOD )
    {
        VarSumSyn <- VAR_STAT_RF %>% dplyr::rename( sample_name = 1 ) 

        colnames(VarSumSyn) <- STAT_COLNAMES[colnames(VarSumSyn)]
        vss_colNameOrders    <- factor(colnames(VarSumSyn), levels=STAT_COLNAMES)

        VarSumSyn <- VarSumSyn[match(SAMPLE_ORDER, VarSumSyn$SampleName), ] %>%
            filter( !is.na(SampleName) )
        
        GrIdx <- SAMPLE_INFO %>% filter( sample_name %in% VarSumSyn$SampleName ) %>% .$sample_group %>% indexing()
        SampleGroupRows <- ldply( GrIdx, function(s) 
        {
            vsinfo     <- SAMPLE_INFO[match(VarSumSyn$SampleName, SAMPLE_INFO$sample_name), ]       
            sinfo.rows <- which(vsinfo$sample_group == s)
            return( data.frame(group = s, start = head(sinfo.rows, 1), end = tail(sinfo.rows, 1)) )
        })

        if( VARIANT_CALL_METHOD == "TS_N"  ){ VCM = "Matched Normal (Tissue)"   }
        if( VARIANT_CALL_METHOD == "ORG_N" ){ VCM = "Matched Normal (Organoid)" }
        if( VARIANT_CALL_METHOD == "Tonly" ){ VCM = "Tumor-only" }
        if( VARIANT_CALL_METHOD == "TonlyGermlineFiltered" ){ VCM = "Tumor-only and Germline Filtered" }

        VarSumSynTable <- VarSumSyn %>%
            relocate( order(vss_colNameOrders) ) %>%
            kable( align=c("l", rep("c", 10)), row.names=F, escape=FALSE, 
                caption = sprintf("Synonymous Variants : %s", VCM )
            ) %>%
            kable_classic("striped", full_width=T, font_size=12, position='center', html_font='Roboto Condensed') %>%
            row_spec(0, align='c', bold=T, background="#D0E4A4") %>%
            column_spec(1, width="6cm")

        if( !is.null(SampleGroupRows) )
        {
            for( n in 1:nrow(SampleGroupRows) )
            { 
                VarSumSynTable <- VarSumSynTable %>% 
                    group_rows(group_label = SampleGroupRows$group[n], SampleGroupRows$start[n], SampleGroupRows$end[n] ) 
            }
        } 

        return(VarSumSynTable)
    }
#------------------------------------------------------------------------------#

#---| VARIANT-STAT-TABLE : FILTER-PASS NON-SYNONYMOUS |------------------------#
#' @param VAR_STAT_RF         re-formed variant class summary table 
#' @param SAMPLE_ORDER        sample orders
#' @param SAMPLE_INFO         sample information ( from DB )
#' @param VARIANT_CALL_METHOD variant calling method. TS_N, ORG_N, or Tonly.
#' @export
    VAR_SUM_TABLE_NonSynonymousPF <- function( VAR_STAT_RF, SAMPLE_ORDER, SAMPLE_INFO, VARIANT_CALL_METHOD )
    {
        VarSumNonSyn_PF <- VAR_STAT_RF[match(SAMPLE_ORDER, VAR_STAT_RF$SampleName), ] %>%
            filter( !is.na(SampleName) )

        GrIdx <- SAMPLE_INFO %>% filter( sample_name %in% VarSumNonSyn_PF$SampleName ) %>% .$sample_group %>% indexing()
        SampleGroupRows <- ldply( GrIdx, function(s) 
        {
            vsinfo     <- SAMPLE_INFO[match(VarSumNonSyn_PF$SampleName, SAMPLE_INFO$sample_name), ]       
            sinfo.rows <- which(vsinfo$sample_group == s)
            return( data.frame(group = s, start = head(sinfo.rows, 1), end = tail(sinfo.rows, 1)) )
        })

        if( VARIANT_CALL_METHOD == "TS_N"  ){ VCM = "Matched Normal (Tissue)"   }
        if( VARIANT_CALL_METHOD == "ORG_N" ){ VCM = "Matched Normal (Organoid)" }
        if( VARIANT_CALL_METHOD == "Tonly" ){ VCM = "Tumor-only" }
        if( VARIANT_CALL_METHOD == "TonlyGermlineFiltered" ){ VCM = "Tumor-only and Germline Filtered" }

        VarSumNonSyn_PF_TB <- VarSumNonSyn_PF %>%
            kable( align=c("l",rep("c", 9)), row.names=F, escape=FALSE, 
                caption = sprintf("Filtered Non-Synonymous Variants : %s", VCM )
            ) %>%
            kable_classic( "striped", full_width=T, font_size=12, position='center', html_font='Roboto Condensed') %>%
            row_spec(0, align='c', bold=T, background="#D0E4A4") %>%
            column_spec(1, width="6cm")

        if( !is.null(SampleGroupRows) )
        {
            for( n in 1:nrow(SampleGroupRows) )
            { 
                VarSumNonSyn_PF_TB <- VarSumNonSyn_PF_TB %>% 
                    group_rows(group_label = SampleGroupRows$group[n], SampleGroupRows$start[n], SampleGroupRows$end[n] ) 
            }
        } 

        return(VarSumNonSyn_PF_TB)
    }
#------------------------------------------------------------------------------#

#---| PRESET GENE VARIANTS TABLE |---------------------------------------------#
#' @param PF_VARS             somatic variants table ( from DB )
#' @param CANCER_TISSUE       cancer tissue 
#' @param GENE_LIST           preset gene list ( hgnc symbols )
#' @param SAMPLE_INFO         sample information ( from DB )
#' @param PRESET_NAME         gene preset name
#' @param VARIANT_CALL_METHOD variant calling method. TS_N, ORG_N, or Tonly.
#' @export
    PRESET_GENE_VARIANTS_TABLE <- function( PF_VARS, CANCER_TISSUE, GENE_LIST, SAMPLE_INFO, PRESET_NAME, VARIANT_CALL_METHOD )
    {
        seqIdList <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% .$seq_id

        cancerType <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% .$sample_info %>% unique()

        if( VARIANT_CALL_METHOD == "TS_N"  ){ VCM = "Matched Normal (Tissue)"   }
        if( VARIANT_CALL_METHOD == "ORG_N" ){ VCM = "Matched Normal (Organoid)" }
        if( VARIANT_CALL_METHOD == "Tonly" ){ VCM = "Tumor-only" }
        if( VARIANT_CALL_METHOD == "TonlyGermlineFiltered" ){ VCM = "Tumor-only and Germline Filtered" }

        presetGeneVarsTable <- PF_VARS %>% 
            filter( variant_call_mode == VARIANT_CALL_METHOD, seq_id %in% seqIdList, GENE_SYMBOL %in% GENE_LIST ) %>%
            mutate( SampleName = SAMPLE_INFO[match(seq_id, SAMPLE_INFO$seq_id), c("sample_name")], varDNA=paste(Refseq_TRANSCRIPT,HGVSc, sep=":")) %>%
            arrange( seq_id, Chromosome, Start_Position) %>% 
            dplyr::select(c("SampleName","GENE_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) %>%
            dplyr::rename(
                'Gene'            = GENE_SYMBOL,
                'DNA Change'      = varDNA, 
                'PROTEIN Change'  = HGVSp_Short,
                'Depth'           = DEPTH, 
                'Tumor Alt Count' = TUMOR_ALT_COUNT
            ) %>%
            kable( align=c("r", rep("c", 6)), row.names = FALSE, escape = FALSE, 
                caption = sprintf("%s ( Preset Genes : %s ) - Variant Call : %s", cancerType, PRESET_NAME, VCM ) 
            ) %>%
            kable_paper(full_width = T, font_size=11, html_font = 'Roboto Condensed', position = 'center' ) %>%
            column_spec(1, width = '3cm') %>%
            row_spec(0, align='center', font_size = 10, bold=T, background = "#D0E4A4") %>% 
            collapse_rows(1)

        return(presetGeneVarsTable)
    }
#------------------------------------------------------------------------------#

#---| PRESET GENE VARIANTS HEATMAP |-------------------------------------------#
#' @param PF_VARS             somatic variants table ( from DB )
#' @param CANCER_TISSUE       cancer tissue 
#' @param GENE_LIST           preset gene list ( hgnc symbols )
#' @param SAMPLE_INFO         sample information ( from DB )
#' @param VARIANT_CALL_METHOD variant calling method. TS_N, ORG_N, or Tonly.
#' @export
    PRESET_GENE_VARIANTS_HEATMAP <- function( PF_VARS, CANCER_TISSUE, GENE_LIST, SAMPLE_INFO, VARIANT_CALL_METHOD )
    {
        seqIdList <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% mutate( sample_name = as.character(sample_name) )
        
        cancerType <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% .$sample_info %>% unique()

        psg_HM <- PF_VARS %>% 
            filter( variant_call_mode == VARIANT_CALL_METHOD, seq_id %in% seqIdList$seq_id, GENE_SYMBOL %in% GENE_LIST ) %>%
            mutate( SampleName = SAMPLE_INFO[match(seq_id, SAMPLE_INFO$seq_id), c("sample_name")]) %>%
            group_by( GENE_SYMBOL, SampleName ) %>%
            reframe( Variant_Class = paste(unique(Variant_Classification), collapse=";") ) %>%
            as.data.frame()
        
        multiVars <- grep(";", psg_HM$Variant_Class)
        if( length(multiVars) > 0 ){ psg_HM[ multiVars, "Variant_Class" ] <- "Multi_Variants" }

        PSG_HM <- reshape(psg_HM , idvar = "GENE_SYMBOL", timevar = "SampleName", direction = "wide")
        colnames(PSG_HM) <- gsub("^Variant_Class.", "", colnames(PSG_HM))

        MissIDs <- seqIdList$sample_name[seqIdList$sample_name %nin% colnames(PSG_HM)]
        if( length(MissIDs) > 0 )
        {
            addMat <- matrix( rep(NA, length(MissIDs)*nrow(PSG_HM)), ncol = length(MissIDs) )
            colnames(addMat) <- MissIDs
            PSG_HM <- cbind(PSG_HM, addMat)
        }

        PSG_HM <- PSG_HM[, c("GENE_SYMBOL", seqIdList$sample_name)]
        rownames(PSG_HM) <- PSG_HM$GENE_SYMBOL

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

        rowFontsize = 10
        colFontsize = 8
        cellWidth   = 6
        cellHeight  = 6

        HM <- Heatmap( as.matrix(data.frame(PSG_HM[,-1])), 
            cluster_rows = TRUE, cluster_columns = FALSE, 
            col = varColors, na_col="#ffffff", rect_gp = gpar(col = "#5b5b5b"),
            row_names_gp = gpar(fontsize = rowFontsize), column_names_gp = gpar(fontsize = colFontsize), 
            row_names_side = "left",
            width = ncol(PSG_HM[,-1])*unit(cellWidth, 'mm'), height = nrow(PSG_HM[,-1])*unit(cellHeight, 'mm'),
            #column_labels = SampleInfo[match(colOrders$id, SampleInfo$seq_id), "sample_name"], 
            column_names_rot = 90, column_names_side = "top",
            heatmap_legend_param = list( title='VARIANT CLASS', title_gp=gpar(fontsize=10),labels_gp=gpar(fontsize=10), legend_width=unit(70,"mm") )
        )

        return(HM)
    }
#------------------------------------------------------------------------------#

#---| GROUP SELECT (ONCOPLOT) |------------------------------------------------#
#' @param Group      group ID ( = group name , "sample_group" in SampleInfo )
#' @param SampleInfo sample info table (from DB)
#' @export
    SelectGroup  <- function( Group = "all", SampleInfo )
    { 
        if( Group == "all" ) { return(unique(SampleInfo$sample_group)) }else{ return(Group) } 
    }
#------------------------------------------------------------------------------#

#---| CANCER SELECT (ONCOPLOT) |-----------------------------------------------#
#' @param cn         cancer type ( "sample_tissue" in SampleInfo )
#' @param SampleInfo sample info table (from DB)
#' @export
    SelectCancer <- function( cn = "all", SampleInfo )
    { 
        if( cn == "all" ){ return(unique(SampleInfo$sample_tissue)) }else{ return(cn) } 
    }
#------------------------------------------------------------------------------#

#---| ONCOPLOT TARGET GENE TABLE |---------------------------------------------#
#' @param selectGroupName  sample group name
#' @param selectCancerType cancer type 
#' @param VariantsMAF      filtered variants list table (from DB)
#' @param topGenes         #N of genes to draw plot
#' @param geneSets         gene symbol list to fraw plot (HGNC symbol)
#' @param useBoth          apply topGenes and geneSets arguments both 
#' @param SampleInfo       sample info table (from DB)
#' @export
    oncoPlotTargetGenes <- function( 
        selectGroupName  = "all", 
        selectCancerType = "all", 
        VariantsMAF, 
        topGenes         = 30, 
        geneSets         = NULL, 
        useBoth          = FALSE, 
        SampleInfo 
    )
    {
        SelectSeqID <- SampleInfo %>% 
            filter( 
                sample_group %in% SelectGroup(Group=selectGroupName,SampleInfo), 
                sample_tissue %in% SelectCancer(cn=selectCancerType,SampleInfo) 
            ) %>% .$seq_id
        
        opTargetGenes <- VariantsMAF %>% 
            filter( seq_id %in% SelectSeqID ) %>% 
            group_by( GENE_SYMBOL ) %>% 
            summarise( length(GENE_SYMBOL) ) %>%
            dplyr::rename( varCount=2 ) %>% 
            arrange( desc(varCount) )
        
        if( is.null(geneSets) )
        {
            if( nrow(opTargetGenes %>% filter( varCount >= opTargetGenes$varCount[topGenes])) > 30 )
            { 
                otg <- opTargetGenes %>% 
                    filter( varCount >  opTargetGenes$varCount[topGenes]) %>% 
                    .$GENE_SYMBOL
            }else{
                otg <- opTargetGenes %>% 
                    filter( varCount >= opTargetGenes$varCount[topGenes]) %>% 
                    .$GENE_SYMBOL
            }
                
            oncoPlotGenes <- VariantsMAF %>% 
                filter( seq_id %in% SelectSeqID, GENE_SYMBOL %in% otg ) %>% 
                .$GENE_SYMBOL %>% 
                table() %>% 
                sort() %>% 
                rev()     
        }else{
            geneSetMAF <- VariantsMAF %>% 
                filter( seq_id %in% SelectSeqID, GENE_SYMBOL %in% geneSets )
                
            if( useBoth )
            {
                gsMafGenes <- geneSetMAF %>% 
                    group_by(GENE_SYMBOL) %>% 
                    summarise(length(GENE_SYMBOL)) %>% 
                    dplyr::rename(varCount=2) %>% 
                    arrange(desc(varCount)) 

                oncoPlotGenes <- geneSetMAF %>% 
                    filter( GENE_SYMBOL %in% (gsMafGenes %>% filter( varCount >= gsMafGenes$varCount[topGenes] ) %>% .$GENE_SYMBOL) ) %>%
                    .$GENE_SYMBOL %>% 
                    table() %>% 
                    sort() %>% 
                    rev() 
            }else{
                oncoPlotGenes <- geneSetMAF %>% 
                    .$GENE_SYMBOL %>% 
                    table() %>% 
                    sort() %>% 
                    rev() 
            }    
        }
        
        oncoPlotGenesCount <- VariantsMAF %>% 
            filter( seq_id %in% SelectSeqID, GENE_SYMBOL %in% names(oncoPlotGenes) ) %>% 
            group_by(seq_id, GENE_SYMBOL) %>% 
            summarise( mutations = paste(unique(Variant_Classification), collapse=";")) 
            
        oncoPlotGenesCount$mutations <- sapply(oncoPlotGenesCount$mutations, function(mt) ifelse( length(grep(";",mt)) == 1, "Multi_Variants", mt ))
            
        return(oncoPlotGenesCount)
    }
#------------------------------------------------------------------------------#

#---| DRAW ONCOPLOT |----------------------------------------------------------#
#' @param targetGeneCountTable 'oncoPlotTargetGenes' function result table
#' @param SampleInfo           sample info table (from DB) 
#' @param SampleOrderByIdGroup sort by sample group id
#' @param includeAllSamples    include all samples even have no variants
#' @export
    drawOncoPlots <- function( 
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
            group_by(GENE_SYMBOL, mutations) %>% 
            summarise( varCount = length(seq_id) ) %>% 
            data.frame() %>% 
            reshape(., idvar="GENE_SYMBOL", timevar="mutations", direction="wide")
    
        annotR_data[is.na(annotR_data)] <- 0 
        colnames(annotR_data) <- gsub("^varCount.", "", colnames(annotR_data))
        rownames(annotR_data) <- annotR_data[,1]
        annotR_data <- annotR_data[, c("GENE_SYMBOL", names(varColors)[which(names(varColors) %in% colnames(annotR_data))]) ]
    
        annotC_data <- targetGeneCountTable %>% 
            group_by(seq_id, mutations) %>% 
            summarise(length(GENE_SYMBOL)) %>% 
            data.frame() %>% 
            reshape(., idvar="seq_id", timevar="mutations", direction="wide")
    
        colnames(annotC_data) <- gsub("^length.GENE_SYMBOL..", "", colnames(annotC_data))
        rownames(annotC_data) <- annotC_data[,1]
        annotC_data[is.na(annotC_data)] <- 0
        annotC_data <- annotC_data[, c("seq_id", names(varColors)[which(names(varColors) %in% colnames(annotC_data))]) ]
    
        oncoPlotMatrix <- targetGeneCountTable %>% data.frame() %>% reshape(., idvar="GENE_SYMBOL", timevar="seq_id", direction="wide")
    
        rownames(oncoPlotMatrix) <- oncoPlotMatrix[,1]
        colnames(oncoPlotMatrix) <- gsub("^mutations.", "", colnames(oncoPlotMatrix))
    
        if( SampleOrderByIdGroup )
        {
            geneOrder <- apply(oncoPlotMatrix[,-1], 1, function(y) length(y[!is.na(y)])) %>% sort() %>% rev()
    
            topOneGeneTissueOrders <- targetGeneCountTable %>% 
                filter( GENE_SYMBOL == (apply(oncoPlotMatrix[,-1], 1, function(y) length(y[!is.na(y)])) %>% 
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
            group_by(GENE_SYMBOL) %>% 
            summarise( sampleFreq = round(length(seq_id)/(ncol(oncoPlotMatrix)-1) *100, 0) ) %>%
            mutate(sampleFreq = gsub("$", "%", sampleFreq))
    
        annotR <- rowAnnotation(
            varRatio   = anno_text( sampleFreq[match(names(geneOrder), sampleFreq$GENE_SYMBOL), ] %>% .$sampleFreq , gp=gpar(fontsize=9) ),
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
    
        if( ncol(inputMat) <= 40 ){ cellWidth  = 4.5 ; colFontsize = 8 }else{ cellWidth = 2.25 ; colFontsize = 6 }
        if( nrow(inputMat) >= 35 ){ cellHeight = 3 ; rowFontsize = 6 }else{ cellHeight = 4 ; rowFontsize = 9 } 
    
        oncoPlotHeatMap <- Heatmap( inputMat, 
            cluster_rows         = FALSE, 
            cluster_columns      = FALSE, 
            col                  = varColors, 
            na_col               = "#f9f9f9",
            border_gp            = gpar(col = "#727272", lwd = 0.5),
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

#---| DRAW VAF-PLOT |----------------------------------------------------------#
#' @param SeqID_1     sample 1 ( default = Tissue sample , X-axis )
#' @param SeqID_2     sampel 2 
#' @param rawMAFList      unfiltered non-synonymous variants MAF format table ( from Rdata ) 
#' @param filteredMAF somatic variants only MAF foramt table ( from DB )
#' @param sampleInfo  sample information ( from DB )
#' @param groupName   group ID
    DRAW_VAF_PLOT <- function( SeqID_1=sampleTS, SeqID_2, rawMAFList, filteredMAF, sampleInfo = sinfo, groupName )
    {
        vf = merge(
            x = rawMAFList[[SeqID_1]] %>% mutate(VAF = round(TUMOR_ALT_COUNT/DEPTH, 3) ) %>% dplyr::select(c("vkey","seq_id","VAF")),   
            y = rawMAFList[[SeqID_2]] %>% mutate(VAF = round(TUMOR_ALT_COUNT/DEPTH, 3) ) %>% dplyr::select(c("vkey","seq_id","VAF")),
            by = 'vkey', all.x = TRUE, all.y = TRUE
        ) %>% dplyr::rename(s1=2, s1_vaf=3, s2=4, s2_vaf=5) %>% mutate(s1_fv="", s2_fv="") %>% mutate(s1 = SeqID_1, s2 = SeqID_2)
        #--------------------------------------------------------------------------------------------------------------------------------------------#
        vf$s1_fv = vf$s2_fv = 1
        vf[which(vf$vkey %in% (filteredMAF %>% filter(seq_id == SeqID_1) %>% .$vkey) ), "s1_fv"] = 4
        vf[which(vf$vkey %in% (filteredMAF %>% filter(seq_id == SeqID_2) %>% .$vkey) ), "s2_fv"] = 2
        vf[is.na(vf$s1_vaf), "s1_vaf"] = 0
        vf[is.na(vf$s2_vaf), "s2_vaf"] = 0
        
        vf$tag = apply( vf[,c("s1_fv","s2_fv")], 1, function(z)
        {
            if( sum(z) == 6 ){ tg = "cm" }
            if( sum(z) == 5 ){ tg = "s1" }
            if( sum(z) == 3 ){ tg = "s2" }
            if( sum(z) == 2 ){ tg = "no" }
            return(tg)
        })
        vf[which(vf$tag == "s1" & vf$s2_vaf != 0 ), "tag" ] = "no"
        vf[which(vf$tag == "s2" & vf$s1_vaf != 0 ), "tag" ] = "no"
        #--------------------------------------------------------------------------------------------------------------------------------------------#
        plotTitle = paste0(
            gsub(paste0(gsub("\\-","_", groupName), "_"), "", sampleInfo[which(sampleInfo$seq_id == SeqID_2), "sample_name"]), " / ", 
            gsub(paste0(gsub("\\-","_", groupName), "_"), "", sampleInfo[which(sampleInfo$seq_id == SeqID_1), "sample_name"])
        )
        #--------------------------------------------------------------------------------------------------------------------------------------------#
        gp_vaf = ggplot(vf, aes(x = s1_vaf, y = s2_vaf)) + geom_point( color = "#eaeaea", alpha = 0.8 ) +
                geom_point(data=(vf %>% filter(tag == "s1")), aes(x = s1_vaf, y = s2_vaf), color = "#e6a722", size = 1.2, alpha = 0.8) + 
                geom_point(data=(vf %>% filter(tag == "s2")), aes(x = s1_vaf, y = s2_vaf), color = "#69aadb", size = 1.2, alpha = 0.8) + 
                geom_point(data=(vf %>% filter(tag == "cm")), aes(x = s1_vaf, y = s2_vaf), color = "#e64522", size = 1.8) +
                theme(
                    axis.line.x  = element_line(colour = "#6f6d6d"), axis.line.y  = element_line(colour = "#6f6d6d"),
                    axis.text.x  = element_text(face = "bold", size = 9),
                    axis.text.y  = element_text(face = "bold", size = 9),
                    axis.title.x = element_text(face = "bold", size = 8, margin = margin(t = 5)),
                    axis.title.y = element_text(face = "bold", size = 8, margin = margin(r = 5)),
                    panel.background = element_rect(fill='white'), 
                    panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'), legend.position  = "none",
                    plot.title=element_text(face="bold", size=8),
                    plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
                ) +
                labs( x = gsub(paste0(groupName, "_"), "", sampleInfo[which(sampleInfo$seq_id == SeqID_1), "sample_name"]), 
                    y = gsub(paste0(groupName, "_"), "", sampleInfo[which(sampleInfo$seq_id == SeqID_2), "sample_name"]) ) +
                coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
                ggtitle(plotTitle) 
        return(gp_vaf)
    }
#------------------------------------------------------------------------------#

#---| INDIVIDUAL MATCHING HEATMAP |--------------------------------------------#
#' @param IND_MATCH_RES individual matching results ( from DB )
#' @param SAMPLE_INFO   sample information ( from DB )
#' @export 
    INDIVIDUAL_MATCH_HEATMAP <- function( IND_MATCH_RES, SAMPLE_INFO )
    {      
        IndivMatchGroups <- unique(SAMPLE_INFO$sample_group)
    
        IMR <- matrix( rep(0, length(SAMPLE_INFO$seq_id)*length(SAMPLE_INFO$seq_id)), ncol=length(SAMPLE_INFO$seq_id))
        rownames(IMR) = colnames(IMR) = SAMPLE_INFO$seq_id
    
        for( i in 1:nrow(IND_MATCH_RES) )
        {
        IMR[ IND_MATCH_RES[i,3], IND_MATCH_RES[i,4] ] = IND_MATCH_RES[i,6]
        IMR[ IND_MATCH_RES[i,4], IND_MATCH_RES[i,3] ] = IND_MATCH_RES[i,6]
        }
        for( k in 1:ncol(IMR)){ IMR[k,k] = 1 }
        
        colnames(IMR) = rownames(IMR) = SAMPLE_INFO[match(colnames(IMR), SAMPLE_INFO$seq_id), 'sample_name']

        groupColors <- RColorBrewer::brewer.pal(length(IndivMatchGroups), "Paired")
        groupColors <- groupColors[ 1:length(IndivMatchGroups) ]
        names(groupColors) <- IndivMatchGroups

        if( ncol(IMR) > 20 ){ unitSize  <- 3  }else{ unitSize  <- 5   }
        if( ncol(IMR) > 20 ){ fontRatio <- 1  }else{ fontRatio <- 1.1 }

        annotC <- HeatmapAnnotation(
            GROUP                   = SAMPLE_INFO[match(colnames(IMR), SAMPLE_INFO$sample_name), 'sample_group'],
            col                     = list( GROUP = groupColors ),
            gap                     = unit(1.4, "mm"),
            annotation_name_gp      = gpar(fontsize = 8*fontRatio), 
            annotation_name_side    = 'left', 
            annotation_legend_param = list( title_gp=gpar(fontsize=8*fontRatio), labels_gp=gpar(fontsize=8*fontRatio) )  
        )

        col_fun = colorRamp2( c(0,0.5,1), c( "#f4f4f4", "#d0e4a4","#375523"), transparency=0 )

        IM_HM = Heatmap( as.matrix(IMR), 
            col                  = col_fun, 
            cluster_columns      = FALSE, 
            cluster_rows         = FALSE, 
            row_names_gp         = gpar(fontsize = 8*fontRatio), 
            column_names_gp      = gpar(fontsize = 8*fontRatio), 
            show_heatmap_legend  = TRUE,
            border_gp            = gpar(col = "#d6d6d6", lwd=0.7), 
            rect_gp              = gpar(col = "#d6d6d6", lwd=0.5),
            width                = ncol(IMR)*unit(unitSize, 'mm'), 
            height               = nrow(IMR)*unit(unitSize, 'mm'),
            heatmap_legend_param = list(
                title            = 'Correlation', 
                title_gp         = gpar(fontsize = 10*fontRatio), 
                labels_gp        = gpar(fontsize = 8*fontRatio),
                legend_direction = 'horizontal', 
                legend_width     = unit(20*fontRatio,"mm"), 
                legend_height    = unit( 3*fontRatio,"mm")
            ),
            top_annotation       = annotC
        )

        return(IM_HM)
    }
#------------------------------------------------------------------------------#

#---| INDIVIDUAL MATCHING TABLE |----------------------------------------------#
#' @param IND_MATCH_RES individual matching results ( from DB )
#' @param SAMPLE_INFO   sample information ( from DB )
#' @export
   INDIVIDUAL_MATCH_TABLE <- function( IND_MATCH_RES, SAMPLE_INFO )
    {
        IND_MATCH_RES$group <- apply(IND_MATCH_RES[,c("seqid_1","seqid_2")], 1, function(y)
        { paste(unique(SAMPLE_INFO[which(SAMPLE_INFO$seq_id %in% y), "sample_group"]), collapse=";") })
        IMRT <- IND_MATCH_RES %>% 
            mutate( 'Sample A' = SAMPLE_INFO[match(seqid_1, SAMPLE_INFO$seq_id), "sample_name"] ) %>%
            mutate( 'Sample B' = SAMPLE_INFO[match(seqid_2, SAMPLE_INFO$seq_id), "sample_name"] ) %>%
            dplyr::select( c("group","Sample A","Sample B", "cor", "match_res") ) %>%
            dplyr::rename(
                'Group'       = group,
                'Correlation' = cor,
                'Matching'    = match_res
            ) %>%
            kbl( escape = FALSE, align=c('c', 'c','c','c'), row.names=F) %>% 
            kable_classic("striped", full_width = T, font_size=11, html_font = 'Roboto Condensed', position = 'center' ) %>%
            row_spec(0, align='center', font_size = 10, bold=T, background = "#D0E4A4") %>%
            column_spec(1, width='8cm') %>% column_spec(2, width='8cm') %>% 
            column_spec(4, width='2cm', background=ifelse( IND_MATCH_RES$match_res == "matched", "#326a1a", "#a72e2e"), color='white') %>%
            collapse_rows(1)
        return(IMRT)
    }
#------------------------------------------------------------------------------#

#---| CONCORDANCE TERM |-------------------------------------------------------#
#' @export 
    CONCORDANCE_TERM <- function()
    {
        concordanceTerm <- data.frame(rbind(
            c(" ", "Concordance" , " Normalized Tanimoto coefficient(Tm) with possible maximum Tm" ),
            c(" ", "MatchRate1"  , " Ratio of common variants to ID1 sample filtered variants"     ),
            c(" ", "MatchRate2"  , " Ratio of common variants to ID2 sample filtered variants"     )
        )) %>%
            kable( align = c("c","c","l"), row.names = FALSE, escape = FALSE, col.names = NULL) %>%
            kable_classic( full_width = FALSE, font_size = 12, position = "left", html_font = "Roboto Condensed") %>%
            column_spec(2, bold = TRUE, width = "4cm", color = "#000000", background = c("#fad9d2","#faedd2","#e1eef7") )
        return(concordanceTerm)
    }
#------------------------------------------------------------------------------#

#---| GET TANIMOTO |-----------------------------------------------------------#
#' @param SET1 input set 1
#' @param SET2 input set 2
#' @export 
    GET_TM = function( SET1,SET2)
    {
        if( length(SET1) != 0 & length(SET2) != 0 )
        {
            a   = length(SET1)
            b   = length(SET2)
            mxT = ifelse( a >= b, b/a, a/b )
            Tm  = round(length(intersect(SET1,SET2))/length(union(SET1,SET2))/mxT, 3)
        }else{
            Tm =0
        }
        return(Tm)
    }
#------------------------------------------------------------------------------#

#---| GET COMMON RATIO |-------------------------------------------------------#
#' @param SET1 input set 1
#' @param SET2 input set 2TRUECONDORDANCE_RATIO_TABLE
#' @param p    match rate direction. 1 = set1, 2 = set2
#' @export 
    GET_COMMON_RATE = function( SET1, SET2, p = 1 )
    {
        if( length(SET1) != 0 & length(SET2) != 0 )
        {
            commonV = intersect( SET1, SET2 )
            if( p == 1 )
            {
                CR = round(length(commonV)/length(SET1) * 100 , 2)
            }else{
                CR = round(length(commonV)/length(SET2) * 100 , 2)
            }
        }else{
            CR = 0
        }
        return(CR)
    }
#------------------------------------------------------------------------------#

#---| GROUP SAMPLES VARIANT CONCORDANCE TABLE |--------------------------------#
#' @param VAR_PF      somatic variants table ( from DB )
#' @param GROUP_ID    group ID
#' @param SAMPLE_INFO sample information ( from DB )
#' @export 
    CONDORDANCE_RATIO_TABLE <- function( VAR_PF, GROUP_ID, SAMPLE_INFO, VARIANT_CALL_METHOD )
    {
        VAR_PF_VCM <- VAR_PF %>% filter( variant_call_mode == VARIANT_CALL_METHOD)
        groupSinfo  <- SAMPLE_INFO %>% filter( sample_group == GROUP_ID )
        pairSeqIDs  <- VAR_PF_VCM %>% 
            filter( variant_call_mode == VARIANT_CALL_METHOD, seq_id %in% groupSinfo$seq_id ) %>% 
            .$seq_id %>% unique()
        if( length(pairSeqIDs) < 2)
        {
            return( data.frame() %>% kable() )
            breaks
        }
        scorePairs <- t(combn(pairSeqIDs, 2))
        # groupSinfo  <- SAMPLE_INFO %>% filter( sample_group == GROUP_ID )
        # scorePairs  <- t(combn(groupSinfo$seq_id, 2))
        conScoreRes <- data.frame()

        for( K in 1:nrow(scorePairs) )
        {
            v1 <- VAR_PF_VCM %>% filter( seq_id == scorePairs[K, 1] ) %>% .$vkey
            v2 <- VAR_PF_VCM %>% filter( seq_id == scorePairs[K, 2] ) %>% .$vkey
            tm <- GET_TM( SET1 = v1, SET2 = v2 )
            rate1 <- GET_COMMON_RATE( SET1 = v1, SET2 = v2, p = 1 )
            rate2 <- GET_COMMON_RATE( SET1 = v1, SET2 = v2, p = 2 )
            conScoreRes <- rbind(
                conScoreRes,
                data.frame(
                    'ID1' = SAMPLE_INFO[which(SAMPLE_INFO$seq_id == scorePairs[K,1]), "sample_name"], 
                    'ID2' = SAMPLE_INFO[which(SAMPLE_INFO$seq_id == scorePairs[K,2]), "sample_name"],
                    'ID1 Variants' = length(v1),
                    'ID2 Variants' = length(v2),
                    'Common Variants'  = length(intersect(v1, v2)),
                    'MatchRate1 (%)'   = rate1,
                    'MatchRate2 (%)'   = rate2,
                    'Concordance'        = tm,
                    check.names = FALSE
                )
            )
        }

        if( VARIANT_CALL_METHOD == "TS_N"  ){ VCM = "matched Normal TISSUE"   }
        if( VARIANT_CALL_METHOD == "ORG_N" ){ VCM = "matched Normal ORGANOIDS" }
        if( VARIANT_CALL_METHOD == "Tonly" ){ VCM = "Tumor-only" }
        if( VARIANT_CALL_METHOD == "TonlyGermlineFiltered" ){ VCM = "Tumor-only and Germline Filtered" }
        
        # sampleBgColors <- apply(groupSinfo[,c("seq_id","sample_origin","matched_normal") ], 1, function(sm)
        # {
        #     if( sm[2] == "TS"  ){ CL <- ifelse( sm[3] == "3", "#DAE8FC", "#F8CECC") }
        #     if( sm[2] == "ORG" ){ CL <- ifelse( sm[3] == "2", "#D5E8D4", "#FFDCB8") }
        #     return(CL)
        # })
        # names(sampleBgColors) <- groupSinfo$sample_name

        conScoreResTB <- conScoreRes %>% 
            kable( align=c("c"), row.names = FALSE, escape = FALSE, 
                caption=sprintf("\t %s - variant calling : %s ", GROUP_ID, VCM)
            ) %>%
            kable_classic("striped", full_width = T, font_size=11, html_font = 'Roboto Condensed', position = 'center' ) %>%
            row_spec(0, align='center', font_size = 10, bold=T, background = "#D0E4A4") %>%
            column_spec(1, width='5cm') %>% 
            column_spec(2, width='5cm') %>%
            # column_spec(1, width='5cm', background=sampleBgColors[conScoreRes$ID1]) %>% 
            # column_spec(2, width='5cm', background=sampleBgColors[conScoreRes$ID2]) %>% 
            column_spec(8, width='2cm', background=ifelse( conScoreRes$Concordance > 0.5, "#326a1a", "#e6a722"), color='white') %>%
            column_spec(6, width='2cm', background=ifelse(  conScoreRes$MatchRate1 > 50,  "#326a1a", "#e6a722"), color='white') %>%
            column_spec(7, width='2cm', background=ifelse(  conScoreRes$MatchRate2 > 50,  "#326a1a", "#e6a722"), color='white')

        return(conScoreResTB)
    }
#------------------------------------------------------------------------------#

#---| VARIANT CONCORDANCE VENN-DIAGRAM |---------------------------------------#
#' @export    
    DRAW_VENN <- function( VENN_LIST, SAMPLE_INFO, NAME, RM_GROUPID=TRUE )
    {
        require(ggVennDiagram)

        vennFaceColoss <- c(
            "1"="#80B1D3",   "2"="#F8CECC", "1/2"="#82B366",
            "3"="#CCEBC5", "1/3"="#B85450", "2/3"="#D79B00", "1/2/3"="#BC80BD",
            "4"="#E1D5E7", "1/4"="#66C2A5", "2/4"="#FCCDE5",   "3/4"="#A6D854",
            "1/2/4" = "#8DA0CB", "1/3/4"="#377EB8","2/3/4"="#FF7F00","1/2/3/4"="#FFED6F"
        )
        SetLabelColors <- c( "TS1"="#B85450", "TS3"="#006EAF", "ORG1"="#2D7600", "ORG2"="#cf7e3b", "ORG0"="#2d2d2d" )

        if( length(VENN_LIST) < 2 )
        {
            VennPlot = ggplot() + theme_void()
        }else{

            VN = ggVennDiagram::Venn(VENN_LIST)
            Vdata <- process_data(VN)

            valueLabel <- venn_regionlabel(Vdata)
            valueLabel <- valueLabel %>% mutate(pct = round(count/sum(count)*100, 1)) %>%
                mutate(lbl = paste0(count,"\n(", pct, "%)"))

            setLabel <- venn_setlabel(Vdata)
            if( nrow(setLabel) == 4 ){ setLabel$X = c(0.15,0.3, 0.6, 0.82) }
            if( nrow(setLabel) == 3 ){ setLabel$X = c(  -2,  5,   2 ) }
            if( nrow(setLabel) == 2 ){ setLabel$X = c(3,3) ; setLabel$Y = c(-1.5,5.5) }

            if( NAME == "id" )
            {
                gid           = unique(SAMPLE_INFO[match(setLabel$name, SAMPLE_INFO$seq_id), "sample_group"])
                setLabel$name = SAMPLE_INFO[match(setLabel$name, SAMPLE_INFO$seq_id), "sample_name"]
            }else{
                gid           = unique(SAMPLE_INFO[match(setLabel$name, SAMPLE_INFO$seq_id), "sample_group"])
            }
            SL <- SAMPLE_INFO[match(setLabel$name, SAMPLE_INFO$sample_name), ] %>% mutate(set_lbl=paste0(sample_origin, matched_normal))
            sl <- SL$set_lbl ; names(sl) = SL$sample_name
            setLabel$clr <- sl[ setLabel$name ]

            if( RM_GROUPID )
            { setLabel$name = gsub( paste0(gsub("\\-","_", gid),"_"), "", setLabel$name) }

            VennPlot <- ggplot() +
                # Circle Region
                geom_polygon( data = venn_regionedge(Vdata), aes(X, Y, fill = id, group = id), show.legend = FALSE, alpha=0.3 ) +
                scale_fill_manual(values=vennFaceColoss) +
                #scale_fill_manual(values=c("1"="#B85450", "2"="#6C8EBF", "1/2"="#82B366")) +
                # Edge
                geom_path(data = venn_setedge(Vdata), aes(X, Y, group = id), color="#FFFFFF", linewidth = 1, show.legend = FALSE) +
                # SetLabelbel in bold
                geom_text(data = setLabel, aes(X, Y, label = name, col = clr), #x=c(3,3), y=c(-1.5,5.5),
                        fontface = "bold", size = 3, hjust=0.5, show.legend = FALSE
                ) +
                scale_color_manual(values=SetLabelColors) +
                # Values
                geom_label(data = valueLabel, aes(X, Y, label = lbl), fontface = "bold", size=3.5, alpha = 0.5 ) +
                theme_void() 
            
            if( length(VENN_LIST) == 2 ){ 
                vp <- VennPlot + coord_flip() 
                EmptyPlot <- ggplot() + theme_void()
                VennPlot  <- cowplot::plot_grid(EmptyPlot, vp, EmptyPlot, ncol=1, rel_heights =c(0.2, 0.6, 0.2) )
            }
        }
        return(VennPlot)
    }
#------------------------------------------------------------------------------#

#---| VARIANT CONCORDANCE VENN-DIAGRAM |---------------------------------------#
#' @param VAR_PF                somatic variants table ( from DB )
#' @param SEQ_ID_1              input seqID 1
#' @param SEQ_ID_2              input seqID 2
#' @param SAMPLE_INFO           sample information ( from DB )
#' @param GROUP_NAME            group name 
#' @param removeGroupNamePrefix remove group name prefix from setname ( default = T )
#' @export
    CONCORDANCE_VENN <- function( VAR_PF, SEQ_ID_1, SEQ_ID_2, SAMPLE_INFO, GROUP_NAME, removeGroupNamePrefix = FALSE, VARIANT_CALL_METHOD )
    {
        require(ggvenn)
        color.Venn <- c( NORMAL="#1BA1E2", TS="#FA6800", ORG="#60A917" )

        vl <- list(
            s1 = VAR_PF %>% filter(seq_id == SEQ_ID_1) %>% .$vkey,
            s2 = VAR_PF %>% filter(seq_id == SEQ_ID_2) %>% .$vkey
        )

        if( removeGroupNamePrefix )
        {
            GN = gsub("\\-", "_", GROUP_NAME)
            vennListName <- gsub( paste0(GN,"_"), "", SAMPLE_INFO[match(c(SEQ_ID_1, SEQ_ID_2), SAMPLE_INFO$seq_id), "sample_name"] )
            vennListName <- gsub("ORG_", "ORG\n", vennListName)
        }else{
            vennListName <- SAMPLE_INFO[match(c(SEQ_ID_1, SEQ_ID_2), SAMPLE_INFO$seq_id), "sample_name"] 
        }
        
        names(vl)  <- vennListName

        conVarVenn <- DRAW_VENN( VENN_LIST = vl, SAMPLE_INFO = SAMPLE_INFO, NAME="name" )
        # sampleType <- SAMPLE_INFO[match(c(SEQ_ID_1, SEQ_ID_2), SAMPLE_INFO$seq_id), "sample_origin"]
        # VennColor  <- color.Venn[ sampleType ]
        # names(VennColor) = NULL

        # conVarVenn <- ggvenn(vl,
        #     stroke_color  = "white", 
        #     set_name_size = 2.5, 
        #     text_size     = 3, 
        #     fill_color    = VennColor, 
        #     auto_scale    = T
        # ) + 
        # theme( plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) #+
        # #coord_cartesian(xlim = c(-1.8, 1.5), ylim = c(-2.5, 2.5))
        # conVarVenn$layers[[3]]$data$x <- c(-0.5, 0.5)

        return(conVarVenn)
    }
#------------------------------------------------------------------------------#

#---| VARIANT CLASS PORTION CHANGE BAR PLOT |----------------------------------#
#' @param VAR_PF                somatic variants table ( from DB )
#' @param SEQ_ID_1              input seqID 1
#' @param SEQ_ID_2              input seqID 2
#' @param VARIANT_GROUP
#' @param SAMPLE_INFO           sample information ( from DB )
#' @param showLegend            diplay legend or not. default = TRUE
#' @param GROUP_NAME            group name 
#' @param removeGroupNamePrefix remove group name prefix from setname ( default = T )
#' @export
    VARIANT_CLASS_PORTION_BARPLOT <- function( VAR_STATS, SEQ_ID_1, SEQ_ID_2, VARIANT_GROUP, SAMPLE_INFO,
        showLegend = T, GROUP_NAME, removeGroupNamePrefix )
    {   
        require(ggalluvial)
        color.Variants <- c( 
            "Missense_Mutation" = "#4DAF4A", "Nonsense_Mutation"      = "#E41A1C", "Frame_Shift_Ins" = "#FF7F00", 
            "Frame_Shift_Del"   = "#8DD3C7", "In_Frame_Ins"           = "#984EA3", "In_Frame_Del"    = "#377EB8", 
            "Nonstop_Mutation"  = "#A65628", "Translation_Start_Site" = "#FFFF33", "Splice_Site"     = "#F781BF", 
            "Multi_Variants"    = "#4c4c4c" 
        )
        sid <- SAMPLE_INFO[match(c(SEQ_ID_1,SEQ_ID_2), SAMPLE_INFO$seq_id ), "sample_name"]    
        
        varClassPortionBarPlot <- rbind(
            VAR_STATS %>% filter( seq_id == SEQ_ID_1, variant_group == VARIANT_GROUP ) %>% 
                group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>%
                mutate( varN_PCT = varN/sum(varN)*100 , ID=SEQ_ID_1 ) ,
            VAR_STATS %>% filter( seq_id == SEQ_ID_2, variant_group == VARIANT_GROUP ) %>% 
                group_by(Variant_Classification) %>% summarise(varN = sum(stats)) %>%
                mutate( varN_PCT = varN/sum(varN)*100 , ID=SEQ_ID_2 )
        ) %>% 
            mutate( SID = SAMPLE_INFO[match(ID, SAMPLE_INFO$seq_id ), "sample_name"] ) 
        
        if( removeGroupNamePrefix )
        {  
            sid = gsub( paste0(gsub("\\-", "_", GROUP_NAME), "_"), "", sid )
            varClassPortionBarPlot <- varClassPortionBarPlot %>% 
                mutate( SID = gsub( paste0(gsub("\\-", "_", GROUP_NAME), "_"), "", SID ) )
        }
        
        varClassPortionBarPlot <- varClassPortionBarPlot %>%
            mutate( SID = factor(SID , levels = sid) ) %>%
            ggplot(aes(x=SID, y=varN_PCT)) +
                geom_flow( aes(alluvium = Variant_Classification), 
                    alpha=0.9, linetype="dotted", fill="#FFFFFF", color="#999999", curve_type="linear", width=0.5
                ) +
                geom_col( aes(fill = Variant_Classification), width=0.5, color="#999999" ) +
                scale_fill_manual( values = color.Variants ) +
                geom_hline( yintercept = 0, col = "#3d3d3d" ) +
                labs( x = "", y = "% of Variant Class") +
                theme(
                    axis.line.y  = element_line(colour = "#3d3d3d"),
                    axis.text.x  = element_text(size = 9, angle=60, vjust=1, hjust=1),
                    axis.text.y  = element_text(size = 10),
                    axis.title.y = element_text(size = 10, margin = margin(r = 3)),
                    axis.ticks.x = element_blank(),
                    legend.text  = element_text(size=9),
                    legend.title = element_text(size=10),
                    legend.key.size = unit(4, 'mm'),
                    panel.background = element_rect(fill="#FFFFFF"), 
                    panel.grid.major = element_line(color = '#c1c1c1', linetype = 'dotted'), 
                    plot.margin  = unit(c(1, 0.5, 1, 0.2), "cm")
                )
            
        if( !showLegend ){ varClassPortionBarPlot = varClassPortionBarPlot + theme(legend.position="none") }
        
        return(varClassPortionBarPlot)
    }
#------------------------------------------------------------------------------#

#---| GENOMIC CHANGE PLOT |----------------------------------------------------#
#' @param VAR_PF                somatic variants table ( from DB )
#' @param SEQ_ID_1              input seqID 1
#' @param SEQ_ID_2              input seqID 2
#' @param SAMPLE_INFO           sample information ( from DB )
#' @param showLegend            diplay legend or not. default = TRUE
#' @param GROUP_NAME            group name 
#' @param removeGroupNamePrefix remove group name prefix from setname ( default = T )
#' @export
    GENOMIC_CHANGE_PLOT <- function( VAR_PF, SEQ_ID_1, SEQ_ID_2, SAMPLE_INFO, showLegend = T, GROUP_NAME, removeGroupNamePrefix )
    {
        color.AlleleChange <- c( 
            "A > C" = "#4DAF4A", "A > G" = "#E41A1C", "A > T" = "#FF7F00", 
            "C > A" = "#8DD3C7", "C > G" = "#984EA3", "C > T" = "#377EB8", 
            "G > A" = "#A65628", "G > C" = "#FFFF33", "G > T" = "#F781BF", 
            "T > A" = "#4c4c4c", "T > C" = "#F8CECC", "T > G" = "#E1D5E7"
        )
        
        sid <- SAMPLE_INFO[match(c(SEQ_ID_1,SEQ_ID_2), SAMPLE_INFO$seq_id ), "sample_name"]

        GenomicChange <- rbind(
            VAR_PF %>% filter( seq_id == SEQ_ID_1, Variant_Type == "SNP" ) %>% 
                mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2)), 
            VAR_PF %>% filter( seq_id == SEQ_ID_2, Variant_Type == "SNP" ) %>% 
                mutate( VAC = paste0( Reference_Allele, " > ", Tumor_Seq_Allele2))
        ) %>% 
            group_by( seq_id, VAC ) %>% summarise( GDC = length(VAC) ) %>%
            group_by( seq_id ) %>% mutate( GDC_PCT = GDC/sum(GDC)*100 ) %>% 
            mutate( SID = SAMPLE_INFO[match(seq_id, SAMPLE_INFO$seq_id ), "sample_name"] )
        
        if( removeGroupNamePrefix )
        {  
            sid = gsub( paste0(gsub("\\-", "_", GROUP_NAME), "_"), "", sid )
            GenomicChange <- GenomicChange %>% 
                mutate( SID = gsub( paste0(gsub("\\-", "_", GROUP_NAME), "_"), "", SID ) )
        }
        
        GenomicChangePlot <- GenomicChange %>% 
            mutate( SID = factor(SID , levels = sid) ) %>%
            ggplot( aes(x=SID, y=GDC_PCT) ) +
                geom_flow( aes(alluvium = VAC), 
                    alpha=0.9, linetype="dotted", fill="#FFFFFF", color="#999999", curve_type="linear", width=0.5
                ) +
                geom_col( aes(fill = VAC), width=0.5, color="#999999") +
                scale_fill_manual(values=color.AlleleChange) +
                geom_hline(yintercept=0, col="#3d3d3d") +
                labs( x = "", y = "% of Genomic Change", fill = "Genomic\nChange") +
                theme(
                    axis.line.y  = element_line(colour = "#3d3d3d"),
                    axis.text.x  = element_text(size = 9, angle=60, vjust=1, hjust=1),
                    axis.text.y  = element_text(size = 10),
                    axis.title.y = element_text(size = 10, margin = margin(r = 3)),
                    axis.ticks.x = element_blank(),
                    legend.text  = element_text(size=9),
                    legend.title = element_text(size=10),
                    legend.key.size = unit(4, 'mm'),
                    panel.background = element_rect(fill="#FFFFFF"), 
                    panel.grid.major = element_line(color = '#c1c1c1', linetype = 'dotted'), 
                    plot.margin  = unit(c(1, 0.2, 1, 0.2), "cm")
                )  
        return(GenomicChangePlot)
    }
#------------------------------------------------------------------------------#

#---| VAF PLOT LEGEND |--------------------------------------------------------#
#' @export 
    VAF_PLOT_LEGEND <- function()
    {
        VafPlotLegend <- data.frame(rbind(
            c(" ", "X,Y-axis"  , " Variant Allele Frequency"              ),
            c(" ", "RED"       , " Common Filtered Variants"              ),
            c(" ", "ORANGE"    , " Sample 1 Only Filtered Variants (TS)"  ),
            c(" ", "BLUE"      , " Sample 2 Only Filtered Variants (ORG)" ),
            c(" ", "GRAY"      , " Unfilitered Variants (Sample1+Sample2)")
        )) %>% 
            kable( align = c("c","c","l"), row.names = FALSE, escape = FALSE, col.names = NULL, caption = "VAF correlation plot legend") %>%
            kable_classic(full_width = FALSE, font_size = 12, position = "left", html_font = "Roboto Condensed") %>%
            column_spec(2, color = c("#000000","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF"), width = "1.5cm", background = c("#FFFFFF","#e64522","#e6a722","#69aadb","#a3a3a3")) %>%
            column_spec(3, width = "7cm")
        return(VafPlotLegend)
    }
#------------------------------------------------------------------------------#

#---| CUMULATIVE SOMATIC VARIANTS |--------------------------------------------#
#' @param SEQ_FOLDERS oncoplot target dataset ids
#' @export
    CUM_ONCOPLOT_DATASETS <- function( SEQ_FOLDERS = NULL )
    {
        dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', port=3306, user='gcx', password='gencurix!!', db = 'gcx_wes' ) 
        on.exit(dbDisconnect(dbCon))
        cumVarPF_seqFolders = dbGetQuery(dbCon, "SELECT seq_folder FROM variants_somatic")
        if(length(SEQ_FOLDERS) == 0 )
        {
            allSeqFolders = paste(paste0( "'", unique(cumVarPF_seqFolders$seq_folder), "'" ), collapse=",")
        }else{
            allSeqFolders = paste(paste0( "'", unique(SEQ_FOLDERS), "'" ), collapse=",")
        }
        cum_sinfo = dbGetQuery(dbCon, sprintf("SELECT * FROM gcx_ngs_service.seqid_info WHERE seq_folder IN (%s)", allSeqFolders))
        cum_sinfo = cum_sinfo %>% filter( matched_normal != 2 )
        cum_var_pf = dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,GENE_SYMBOL,Variant_Classification FROM variants_somatic WHERE seq_folder IN (%s)", allSeqFolders))
        res = list(
            cumSinfo = cum_sinfo,
            cumVarPF = cum_var_pf
        )
        return(res)
    }
#------------------------------------------------------------------------------#

#---| HLA-TYPING RESULT LIST |-------------------------------------------------#
#' @param HLA_RES      HLA type result table ( from DB )
#' @param SAMPLE_INFO  sample information ( from DB )
#' @param SAMPLE_ORDER sample orders
#' @export 
    HLA_TYPE_RESULT <- function( HLA_RES, SAMPLE_INFO, SAMPLE_ORDER )
    {
        tissueList <- unique(SAMPLE_INFO$sample_tissue)
        tidx       <- sapply(tissueList, function(y) list(y))  

        hla_tissue_res <- lapply(tidx, function(tx) {
            mhc.tissue = HLA_RES %>% dplyr::select(c('seq_id','locus','allele1_hla_la','allele2_hla_la')) %>% 
                left_join( . , SAMPLE_INFO[,c('seq_id','sample_name','sample_tissue')], by='seq_id' ) %>%
                dplyr:: filter(sample_tissue %in% tx)
        
            mhc.tissue$HLA = apply(mhc.tissue[,c('allele1_hla_la','allele2_hla_la')], 1, function(y) paste(y, collapse="#") )
        
            MHC.tissue = mhc.tissue %>% 
                dplyr::select(c('sample_name','locus','HLA')) %>% 
                mutate( HLA = gsub("#","\\<br\\>", HLA) ) %>% 
                mutate( HLA = gsub("\\*","\\* ", HLA) ) %>% 
                filter( locus %nin% c("E","F","G") ) %>%
                reshape( idvar='sample_name', timevar='locus', direction='wide' ) %>% 
                dplyr::rename( 'Sample_ID' = 1 )

            if( nrow(MHC.tissue) > 0 )
            {
                colnames(MHC.tissue) = gsub("^HLA.D","D", colnames(MHC.tissue))
                colnames(MHC.tissue) = gsub("\\.","-", colnames(MHC.tissue))
                MHC.tissue = MHC.tissue[match(SAMPLE_ORDER, MHC.tissue[,1]), ] 
                MHC.tissue = MHC.tissue[!is.na(MHC.tissue$Sample_ID), ] 
            }else{
                MHC.tissue = data.frame()
            }
            
            return(MHC.tissue)
        })
        return(hla_tissue_res)
    }
#------------------------------------------------------------------------------#

#---| HLA-TYPING RESULT TABLE |------------------------------------------------#
#' @param HLA_RES_LIST HLA result list ( derived from 'HLA_TYPE_RESULT' )
#' @param TISSUE       cancer tissue name
#' @param SAMPLE_INFO  ample information ( from DB )
#' @param NO_MATCH     optitype adn hla-la mis-matched alleles check
#' @export
    HLA_RESULT_TABLE <- function( HLA_RES_LIST, TISSUE, SAMPLE_INFO, NO_MATCH )
    {
        tissueName = SAMPLE_INFO %>% filter( sample_tissue == TISSUE ) %>% .$sample_info %>% unique()
        hla.sinfo  = SAMPLE_INFO %>% filter( sample_name %in% HLA_RES_LIST[[ TISSUE ]]$Sample_ID )

        if( length(which(!is.na(hla.sinfo$sample_group))) > 0 )
        {
            hla.groups = unique(hla.sinfo$sample_group)
            hla.groups = hla.groups[!is.na(hla.groups)]
            is.hla.groups = TRUE
        }else{
            hla.groups = NA
            is.hla.groups = FALSE
        }
        htr = HLA_RES_LIST[[ TISSUE ]]

        if( nrow(htr) == 0 )
        {
            next
        }else{
            if( is.hla.groups )
            {
                no.group.samples = htr$Sample_ID[ which(htr$Sample_ID %nin% (hla.sinfo %>% filter(sample_group %in% hla.groups ) %>% .$sample_name) ) ]
                
                htr2 = rbind(
                    htr %>% filter( Sample_ID %in%  no.group.samples),
                    htr %>% filter( Sample_ID %nin% no.group.samples)
                )
                hla.res = htr2 %>%
                    kbl( escape = FALSE, align=c('r', rep('c',10)), row.names=F, caption = tissueName ) %>% 
                    kable_classic("striped", full_width = T, font_size=8, html_font = 'Roboto Condensed' ) %>%
                    add_header_above(c(" "=1,"MHC-I"=3, "MHC-II"=7), font_size=15, bold=T) %>% 
                    row_spec(0, align='center', font_size = 12, bold=T, background = "#D0E4A4") %>% 
                    column_spec(1, width='3cm')
                
                for( k in 1:length(hla.groups) )
                {
                    hla.group.rows = which(htr2$Sample_ID %in% (hla.sinfo %>% filter(sample_group == hla.groups[k] ) %>% .$sample_name) )
                    hla.res = hla.res %>% group_rows(group_label= hla.groups[k], min(hla.group.rows), max(hla.group.rows) )
                }
                
                is.noMatch = ifelse( length(which(htr2$Sample_ID %in% NO_MATCH$sample_name)) > 0, TRUE, FALSE )
                if( is.noMatch )
                {
                    hla.res = hla.res %>% 
                    column_spec(1, color=sapply(htr2$Sample_ID, function(y) ifelse( y %in% NO_MATCH$sample_name, "#e64522", "#000000"))) 
                }       
                
            }else{
                hla.res = htr %>%
                    kbl( escape = FALSE, align=c('r', rep('c',10)), row.names=F, caption = tissueName ) %>% 
                    kable_classic("striped", full_width = T, font_size=8, html_font = 'Roboto Condensed' ) %>%
                    add_header_above(c(" "=1,"MHC-I"=3, "MHC-II"=7), font_size=15, bold=T) %>% 
                    row_spec(0, align='center', font_size = 12, bold=T, background = "#D0E4A4") %>% 
                    column_spec(1, width='3cm')
                
                is.noMatch = ifelse( length(which(htr$Sample_ID %in% NO_MATCH$sample_name)) > 0, TRUE, FALSE )
                if( is.noMatch )
                {
                    hla.res = hla.res %>% 
                    column_spec(1, color=sapply(htr$Sample_ID, function(y) ifelse( y %in% noMatch$sample_name, "#e64522", "#000000"))) 
                }    
            }
        }
        return(hla.res)
    }
#------------------------------------------------------------------------------#

#---| HLA-NO MATCHED SAMPLE RESULT TABLE |-------------------------------------#
#' @param NO_MATCH_RES no matched result table
#' @export 
    HLA_NO_MATCH_RES_TABLE <- function( NO_MATCH_RES )
    {
        NoMatchResTb <- NO_MATCH_RES %>% mutate( 
            sw1 = paste(allele1_hla_la, allele2_hla_la, sep="#"), 
            sw2 = paste(allele1_optitype, allele2_optitype, sep="#")
        ) %>%
            dplyr::select(c('sample_name','locus','sw1','sw2')) %>% 
            mutate(locus = gsub("^","HLA-",locus)) %>%
            mutate( sw1 = gsub("#","\\<br\\>", sw1) ) %>% 
            mutate( sw2 = gsub("#","\\<br\\>", sw2) ) %>%
            mutate( sw1 = gsub("\\*","\\* ", sw1) ) %>% 
            mutate( sw2 = gsub("\\*","\\* ", sw2) ) %>%
            dplyr::rename('Sample_ID'=1, 'Locus'=2, 'HLA-LA'=3, OptiType=4) %>%
                kbl( escape = FALSE, align=c('c', rep('c',3)), row.names=F) %>% 
                kable_classic("striped", full_width = F, font_size=10, html_font = 'Roboto Condensed', position = 'left' ) %>%
                row_spec(0, align='center', font_size = 12, bold=T, background = "#D0E4A4") %>%
                column_spec(1, width='6cm') %>% 
                column_spec(2, width='3cm') %>% column_spec(3:4, width='4cm')
        return(NoMatchResTb)
    }
#------------------------------------------------------------------------------#


#---| TMB TABLE |--------------------------------------------------------------#
#' @param TMB_RES      tmb result table ( from DB )
#' @param SAMPLE_INFO  sample information ( from DB )
#' @export 
    TMB_TABLE <- function( TMB_RES, SAMPLE_INFO, TISSUE )
    {
        TMB_RES$variant_call = sapply( TMB_RES$variant_call_mode, function(y) 
        {
            if( y == "Tumor-only" ){ vc = "Tumor-Only" 
            }else if( y == "TonlyGermlineFiltered" ){ vc = "Tumor-Only" 
            }else{ vc = "Matched-Normal" }
            return(vc)
        })
        TMB_RES$status       = sapply( TMB_RES$TMB_value, function(y) ifelse( as.numeric(y) > 10, "TMB High", ""))
        tmbRes <- TMB_RES %>%
            left_join(., SAMPLE_INFO[,c("seq_id","sample_name","sample_info","sample_group")], by='seq_id') %>%
            arrange( sample_info, seq_id ) %>%
            dplyr::select(c("sample_info","sample_group","sample_name","TMB_value","status","variants","variant_call")) %>%
            dplyr::rename(
                "Tissue"        = sample_info,
                "Group"         = sample_group,
                "Sample ID"     = sample_name,
                "TMB (Mut/Mb)"  = TMB_value,
                "TMB Status"    = status,
                "Variants"      = variants,
                "Variant Call"  = variant_call
            ) 
        if( TISSUE != "all" ){ 
            tmbRes = tmbRes %>% filter( Tissue == TISSUE)
            tmb_caption = TISSUE
        }else{
            tmb_caption = ""
        }
        tmbRes = tmbRes %>% arrange(Group)
        tmbResTable <- tmbRes %>% arrange(Group) %>%
            kbl( escape = FALSE, align=c('c','c','l','c','c','c'), row.names=F, caption = tmb_caption ) %>% 
            kable_classic( full_width = T, font_size=11, html_font = 'Roboto Condensed' ) %>%
            row_spec(0, align='center', font_size = 10, bold=T, background = "#D0E4A4") %>% 
            column_spec(4, color = sapply(tmbRes[,4], function(y) ifelse( as.numeric(y) > 10, "#e64522", "#000000"))) %>%
            column_spec(5, background = sapply(tmbRes[,5], function(y) ifelse( y == "TMB High", "#e64522", "#FFFFFF" )), color = "#FFFFFF") %>%
            collapse_rows(c(1,2))
        return(tmbResTable)
    }     
#------------------------------------------------------------------------------#

#---| TMB GROUP-POINT-PLOT |---------------------------------------------------#
#' @param TMB_RES     tmb result table ( from DB )
#' @param SAMPLE_INFO sample information ( from DB )
#' @param GROUP_NAME  group name 
#' @export
    TMB_PLOT <- function( TMB_RES, SAMPLE_INFO, GROUP_NAME )
    {
        tmbRes <- TMB_RES %>%
            left_join(., SAMPLE_INFO, by='seq_id') %>%
            filter( sample_group == GROUP_NAME , matched_normal != "3" ) %>%
            mutate( sample_origin = factor(sample_origin, levels = c("TS","ORG")))
        tmbRes$NT = sapply( tmbRes$matched_normal, function(y) ifelse( y %in% c(0,1), "CANCER", "NORMAL") )

        if( max( tmbRes$TMB_value ) < 10 ){ maxY = ifelse( max( tmbRes$TMB_value ) > 5, 11, 5 ) 
        }else{ maxY = max( tmbRes$TMB_value ) + 2 }

        typeColor <- c( "TS"="#FA6800", "ORG"="#60A917" )

        tmbPlot <- ggplot() +
            geom_point( data=tmbRes %>% filter(matched_normal==1), aes( x = sample_origin, y = TMB_value, col = sample_origin, shape = NT ), size=5, alpha=0.7 ) +
            geom_point( data=tmbRes %>% filter(matched_normal==2), aes( x = sample_origin, y = TMB_value, col = sample_origin, shape = NT ), size=5, alpha=0.7 ) +
            coord_cartesian( ylim = c(0, maxY)) +
            scale_colour_manual(values=typeColor) + 
            scale_shape_manual(values=c(19, 15)) +
            geom_hline( yintercept = 0, col = "#3d3d3d" ) + 
            labs( x = "", y = "TMB (Mut/Mb)", color = "TYPE", shape = "GROUP" ) +
            theme(
                axis.line.y  = element_line(colour = "#3d3d3d"),
                axis.text.x  = element_text(size = 12),
                axis.text.y  = element_text(size = 12),
                axis.title.y = element_text(size = 12, margin = margin(r = 3)),
                axis.ticks.x = element_blank(),
                legend.text  = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(4, 'mm'),
                panel.background = element_rect(fill="#FFFFFF"), 
                panel.grid.major = element_line(color = '#c1c1c1', linetype = 'dotted'), 
                plot.margin  = unit(c(1, 0.5, 1, 0.2), "cm")
            ) + 
            ggtitle( GROUP_NAME )

        return(tmbPlot)
    }   
#------------------------------------------------------------------------------#

#---| TMB CANCER-TYPE-PLOT |---------------------------------------------------#
#' @param TMB_ALL_RES
#' @param SAMPLE_INFO
#' @param CANCER_TYPE
#' @export
    TMB_BOXPLOT <- function( TMB_ALL_RES, SAMPLE_INFO, CANCER_TYPE, MASK_MAX = false )
    {
        if( CANCER_TYPE == "all" ){ CANCER_TYPE = unique(SAMPLE_INFO$sample_info) }
        tmbRes <- TMB_ALL_RES %>%
            left_join(., SAMPLE_INFO, by='seq_id') %>%
            filter( sample_info %in% CANCER_TYPE)
        tmbRes$RunMode = sapply( tmbRes$variant_call_mode, function(y) ifelse( y == 'tonly', "Tumor-only", "Matched-Normal") )

        if( max( tmbRes$TMB_value ) < 10 ){ maxY = ifelse( max( tmbRes$TMB_value ) > 5, 11, 5 ) 
        }else{ maxY = max( tmbRes$TMB_value ) + 2 }

        if( MASK_MAX )
        {
            if( maxY > 40 ){ maxY = 45 }
            tmbRes[which(tmbRes$TMB_value > 45), "TMB_value"] = 45
        }
        typeColor <- c( "TS"="#FA6800", "ORG"="#60A917" )

        tmbBoxplot <- ggplot(tmbRes, aes( x = sample_info, y = TMB_value) ) +
            geom_hline( yintercept = 10, col = "#14426c", linetype = 'dotted', alpha=0.5 ) + 
            geom_violin( aes( fill = sample_info ), alpha=0.3, col = "#FFFFFF", trim = FALSE, bw = 1.5 ) +
            geom_hline( yintercept = 0, col = "#3d3d3d" ) + 
            geom_jitter( aes( shape = RunMode, col = sample_origin ), width = 0.2, size = 2, alpha=0.9  ) + 
            labs( x = "", y = "TMB (Mut/Mb)", color = "GROUP", fill = "CANCER", shape = "Variant Call" ) +
            coord_cartesian( ylim = c(0, maxY)) +
            theme(
                axis.line.y  = element_line(colour = "#3d3d3d"),
                axis.text.x  = element_text(size = 12),
                axis.text.y  = element_text(size = 12),
                axis.title.y = element_text(size = 12, margin = margin(r = 3)),
                axis.ticks.x = element_blank(),
                legend.text  = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(4, 'mm'),
                panel.background = element_rect(fill="#FFFFFF"), 
                #panel.grid.major = element_line(color = '#c1c1c1', linetype = 'dotted'), 
                plot.margin  = unit(c(1, 0.5, 1, 0.2), "cm")
            ) +
            scale_colour_manual(values=typeColor) +
            scale_shape_manual(values=c(16,13))

        return(tmbBoxplot)
    }
#------------------------------------------------------------------------------#

#---| MSI TABLE |--------------------------------------------------------------#
#' @param MSI_RES     msi-status result ( from DB )
#' @param SAMPLE_INFO sample information ( from DB )
#' @param TISSUE      filter result by specified tissue. default = "all" (no filter)
#' @export 
    MSI_TABLE <- function( MSI_RES, SAMPLE_INFO, TISSUE="all" )
    {
        msiRes <- MSI_RES %>%
            left_join(., SAMPLE_INFO[,c("seq_id","sample_name","sample_info","sample_group","matched_normal")], by='seq_id') %>%
            arrange( sample_info, sample_group, desc(matched_normal), seq_id ) 
        msiRes$status = apply( msiRes[,c("msi_run_mode","mantis_msi_status","msisensor2_msi_status")], 1, function(y)
        {
            if( y[1] == "Tumor-Only" )
            {
                sts = ifelse( y[3] == "MSI-H", "MSI-H", "MSS" )
            }else{
                sts = ifelse( y[2] == "MSI-H" & y[3] == "MSI-H", "MSI-H" , "MSS" )
            }
            return(sts)
        })

        msiRes <- msiRes %>%
                dplyr::select(c("sample_info","sample_group","sample_name","mantis_score","msisensor2_score","status","msi_run_mode")) %>%
                dplyr::rename(
                    "Tissue"     = sample_info,
                    "Group"      = sample_group,
                    "Sample ID"  = sample_name,
                    "MANTIS"     = mantis_score,
                    "MSIsensor2" = msisensor2_score,
                    "Status"     = status,
                    "Run Mode"   = msi_run_mode
                )
        if( TISSUE != "all" )
        { 
            msiRes <- msiRes %>% filter( Tissue == TISSUE ) 
            msi_caption <- TISSUE
        }else{
            msi_caption <- ""
        }

        mantis_background = sapply( msiRes$MANTIS, function(y) { 
            if( y == "-" ){ bgclr = "#FFFFFF" }else{ bgclr = ifelse( as.numeric(y) > 0.4 , "#e64522", "#FFFFFF") }
            return(bgclr) })
        mantis_color = sapply( msiRes$MANTIS, function(y) { 
            if( y == "-" ){ mclr = "#000000" }else{ mclr = ifelse( as.numeric(y) > 0.4 , "#FFFFFF", "#000000") } 
            return(mclr) })
        sensor_background = sapply( msiRes$MSIsensor2, function(y) { 
            bgclr = ifelse( as.numeric(y) > 20 , "#e64522", "#FFFFFF")  
            return(bgclr) })
        sensor_color = sapply( msiRes$MSIsensor2, function(y) { 
            bgclr = ifelse( as.numeric(y) > 20 , "#FFFFFF", "#000000")  
            return(bgclr) })
        status_background = sapply( msiRes$Status, function(y) ifelse( y == "MSI-H", "#e64522","#FFFFFF") )
        status_color = sapply( msiRes$Status, function(y) ifelse( y == "MSI-H", "#FFFFFF","#000000") )

        msiResTable <- msiRes %>%
            kbl( escape = FALSE, align=c('c','c','l','c','c','c','c'), row.names=F, caption = msi_caption ) %>% 
            kable_classic( full_width = T, font_size=12, html_font = 'Roboto Condensed' ) %>%
            row_spec(0, align='center', font_size = 11, bold=T, background = "#D0E4A4") %>% 
            column_spec(4, background = mantis_background, color = mantis_color ) %>%
            column_spec(5, background = sensor_background, color = sensor_color ) %>%
            column_spec(6, background = status_background, color = status_color ) %>%
            collapse_rows(c(1,2)) 
        return(msiResTable)
    }     
#------------------------------------------------------------------------------#

#---| BLANK PLOTS |------------------------------------------------------------#
    empty_ts_n_label   = "No two or more variant results\nwith matched-normal TISSUE"
    empty_org_n_label  = "No two or more variant results\nwith matched-normal ORGANOID"
    empty_tonly_label  = "No two or more variant results\nfrom Tumor-only variant call"
    empty_no_res_label = "No Results fill condition"
    EMPTY_PLOT <- ggplot() + theme_void()
    EMPTY_PLOT_TS_N  <- ggplot() + annotate("text", x=1, y=1, size=4, label=empty_ts_n_label) + theme_void()
    EMPTY_PLOT_ORG_N <- ggplot() + annotate("text", x=1, y=1, size=4, label=empty_org_n_label) + theme_void()
    EMPTY_PLOT_TONLY <- ggplot() + annotate("text", x=1, y=1, size=4, label=empty_tonly_label) + theme_void()
    EMPTY_PLOT_NO_RES <- ggplot() + annotate("text", x=1, y=1, size=4, label=empty_no_res_label) + theme_void()

#------------------------------------------------------------------------------#

#---| CUMULATIVE SOMATIC VARIANTS2 |--------------------------------------------#
#' @param SEQ_FOLDERS oncoplot target dataset ids
#' @export
    CUM_ONCOPLOT_DATASETS2 <- function( TISSUE_INFO = NULL, VARIANT_CALL_MODE = "Tonly", CLIENT_ID=NULL, PANEL="WES")
    {
        dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', port=3306, user='gcx', password='gencurix!!', db = 'gcx_wes' ) 
        on.exit(dbDisconnect(dbCon))
        if( is.null(TISSUE_INFO) )
        {
            cum_sinfo <- dbGetQuery(dbCon, "SELECT * FROM gcx_ngs_service.seqid_info")
        }else{
            cum_sinfo <- dbGetQuery(dbCon, sprintf("SELECT * FROM gcx_ngs_service.seqid_info WHERE sample_info = '%s'", TISSUE_INFO))
        }
        allseqID = paste(paste0("'", cum_sinfo$seq_id, "'"), collapse=",")

        cum_var_pf = dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,GENE_SYMBOL,Variant_Classification FROM variants_somatic WHERE variant_call_mode = '%s' AND seq_id IN (%s)", VARIANT_CALL_MODE, allseqID))
        cum_sinfo = cum_sinfo %>% filter( panel %in% PANEL )
        if( is.null(CLIENT_ID) )
        {
            res = list(
                cumSinfo = cum_sinfo,
                cumVarPF = cum_var_pf
            )
        }else{
            cum_sinfo = cum_sinfo %>% filter(client_facility_id == CLIENT_ID)
            cum_var_pf = cum_var_pf %>% filter( seq_folder %in% cum_sinfo$seq_folder )
            res = list(
                cumSinfo = cum_sinfo,
                cumVarPF = cum_var_pf
            )
        }
        return(res)
    }
#------------------------------------------------------------------------------#
