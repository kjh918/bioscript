
### WES ANALYSIS REPORT MODULES : COMMON MODULES ###

#---| PRE-DEFINED DATA |------------------------------------------------------------------------------------------------------------------------------
#---/ variant stat table columns /------------------------------------------------------------------
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
#
#---/ sample info table /---------------------------------------------------------------------------
    WES_report.Sample.Info.Table <- function( SampleInfo=NULL, FontSize=10 )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        #-----------------------------------------------------------------------
        if( is.null(SampleInfo) ){ stop("No sample information table. STOPPED.") }
        # sample orders --------------------------------------------------------
        SampleOrderTable <- rbind(
            SampleInfo[which(SampleInfo$sample_group != ""), c("seq_id","sample_group","sample_info","sample_name","sample_label","sample_tissue")] %>% 
                arrange( sample_info, sample_group, seq_id ),
            SampleInfo[which(SampleInfo$sample_group == ""), c("seq_id","sample_group","sample_info","sample_name","sample_label","sample_tissue")] %>% 
                arrange( sample_info, seq_id) %>% mutate(sample_group = "No Grouped")
        )
        # sample group index ---------------------------------------------------
        SampleGroupIndex <- ldply( lapply( indexing( SampleOrderTable$sample_group), function(sg)
            data.frame( group=sg, start=head(which( SampleOrderTable$sample_group == sg ), 1), end=tail(which( SampleOrderTable$sample_group == sg ), 1) )
        ))
        # output table ---------------------------------------------------------
        SinfoTableOutput <- SampleOrderTable %>% mutate( INDEX=" " ) %>% dplyr::select(c(
            "INDEX", "sample_name","sample_label","sample_tissue","sample_info" 
        )) %>% dplyr::rename(' '=1, 'SAMPLE ID'=2, 'SAMPLE NAME'=3, 'SAMPLE TISSUE'=4, 'CANCER CODE'=5 ) %>%
            kbl( escape=FALSE, align=c('c','l','c','c','c')) %>%
            kable_styling(bootstrap_options="condensed", font_size=FontSize, full_width = TRUE ) %>%
            row_spec(0, align ='c', bold=TRUE, background="#82B366", color="#FFFFFF") %>%
            column_spec(c(1), width_min="1em",  width="2em" ) %>% 
            column_spec(c(2), width_min="20em", width="23em") %>% 
            column_spec(c(3), width_min="15em", width="20em")  
        for( k in 1:nrow(SampleGroupIndex) )
        {
            SinfoTableOutput <- SinfoTableOutput %>% pack_rows( 
                group_label = SampleGroupIndex[k, "group"],  
                start_row   = SampleGroupIndex[k, "start"], 
                end_row     = SampleGroupIndex[k, "end"  ], 
                indent      = FALSE,
                background  = '#D5E8D4',
                label_row_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #82B366;"
            )
        }
        return(SinfoTableOutput)
    }
#
#---/ sample info table (Organoid Science) /--------------------------------------------------------
    WES_report.Sample.Info.Table.OrganoidScience <- function( SampleInfo=NULL, FontSize=10 )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        #-----------------------------------------------------------------------
        if( is.null(SampleInfo) ){ stop("No sample information table. STOPPED.") }
        # sample orders --------------------------------------------------------
        SampleOrderTable <- rbind(
            SampleInfo[which(SampleInfo$sample_group != ""), c("seq_id","sample_group","sample_info","sample_name","sample_label","sample_tissue","sample_origin")] %>% 
                arrange( sample_info, sample_group, seq_id ),
            SampleInfo[which(SampleInfo$sample_group == ""), c("seq_id","sample_group","sample_info","sample_name","sample_label","sample_tissue","sample_origin")] %>% 
                arrange( sample_info, seq_id) %>% mutate(sample_group = "No Grouped")
        )
        SampleOrderTable$sample_origin <- gsub("TS", "Tissue",    SampleOrderTable$sample_origin)
        SampleOrderTable$sample_origin <- gsub("ORG", "Organoid", SampleOrderTable$sample_origin)
        SampleOrderTable$sample_tissue <- capitalize(SampleOrderTable$sample_tissue)
        # sample group index ---------------------------------------------------
        SampleGroupIndex <- ldply( lapply( indexing( SampleOrderTable$sample_group), function(sg)
            data.frame( group=sg, start=head(which( SampleOrderTable$sample_group == sg ), 1), end=tail(which( SampleOrderTable$sample_group == sg ), 1) )
        ))
        # output table ---------------------------------------------------------
        SinfoTableOutput <- SampleOrderTable %>% mutate( INDEX=" " ) %>% dplyr::select(c(
            "INDEX", "sample_name","sample_label","sample_origin","sample_tissue","sample_info" 
        )) %>% dplyr::rename(' '=1, 'SAMPLE ID'=2, 'SAMPLE NAME'=3, 'SAMPLE TYPE'=4, 'SAMPLE TISSUE'=5, 'CANCER CODE'=6 ) %>%
            kbl( escape=FALSE, align=c('c','l','c','c','c')) %>%
            kable_styling(bootstrap_options="condensed", font_size=FontSize, full_width = TRUE ) %>%
            row_spec(0, align ='c', bold=TRUE, background="#82B366", color="#FFFFFF") %>%
            column_spec(c(1), width_min="1em",  width="2em" ) %>% 
            column_spec(c(2), width_min="20em", width="23em") %>% 
            column_spec(c(3), width_min="15em", width="20em")  
        for( k in 1:nrow(SampleGroupIndex) )
        {
            SinfoTableOutput <- SinfoTableOutput %>% pack_rows( 
                group_label = SampleGroupIndex[k, "group"],  
                start_row   = SampleGroupIndex[k, "start"], 
                end_row     = SampleGroupIndex[k, "end"  ], 
                indent      = FALSE,
                background  = '#D5E8D4',
                label_row_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #82B366;"
            )
        }
        return(SinfoTableOutput)
    }
#
#---/ ngs qc stat summary table /-------------------------------------------------------------------
    WES_report.NGS.QC.Stat.Summary.Table <- function( DataQC=NULL , SampleInfo=NULL, LowQualSeqID=NULL )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        #-----------------------------------------------------------------------
        ReportFeatures <- c('seq_id','rds_pf','pct_q30_pf','pct_gc_content','pct_rds_dup','pct_base_on_target','depth_target_mean','pct_base_target_30x','pct_base_target_50x','pct_base_target_100x')
        ReportColumns <- c("sample_name","sample_origin","rds_pf","pct_q30_pf","pct_gc_content","pct_rds_dup","pct_base_on_target","depth_target_mean","pct_base_target_30x","pct_base_target_50x","pct_base_target_100x")
        #-----------------------------------------------------------------------
        if( is.null(DataQC)     ){ stop("No QC data. STOPPED.") }
        if( is.null(SampleInfo) ){ stop("No sample information table. STOPPED.") }
        # sample orders --------------------------------------------------------
        SampleOrderTable <- rbind(
            SampleInfo[which(SampleInfo$sample_group != ""), c("seq_id","sample_group","sample_info","sample_name","sample_label","sample_tissue")] %>% 
                arrange( sample_info, sample_group, seq_id ),
            SampleInfo[which(SampleInfo$sample_group == ""), c("seq_id","sample_group","sample_info","sample_name","sample_label","sample_tissue")] %>% 
                arrange( sample_info, seq_id) %>% mutate(sample_group = "No Grouped")
        )
        # sample group index ---------------------------------------------------
        SampleGroupIndex <- ldply( lapply( indexing( SampleOrderTable$sample_group), function(sg)
            data.frame( group=sg, start=head(which( SampleOrderTable$sample_group == sg ), 1), end=tail(which( SampleOrderTable$sample_group == sg ), 1) )
        ))
        # qc data table --------------------------------------------------------
        InputData <- DataQC[match(SampleOrderTable$seq_id, DataQC$seq_id), ] %>% 
            left_join(. , sinfo %>% dplyr::select(c('seq_id','sample_name','sample_group','sample_origin','sample_tissue')), by='seq_id')
            
        QcTableOutput <- InputData %>%
            dplyr::select(all_of(ReportColumns)) %>%
            mutate( rds_pf = formatC(round(rds_pf), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
            mutate_at(c(4:11), round,1) %>%
            dplyr::rename(
                'SAMPLE ID'            = sample_name, 
                'SAMPLE TYPE'          = sample_origin, 
                'TOTAL READS'          = rds_pf, 
                'Q30 RATE'             = pct_q30_pf, 
                'GC RATIO'             = pct_gc_content,
                'DUPLICATES RATE'      = pct_rds_dup, 
                'ON-TARGET RATE'       = pct_base_on_target, 
                'TARGET DEPTH'         = depth_target_mean, 
                'TARGET COVERAGE 30X'  = pct_base_target_30x, 
                'TARGET COVERAGE 50X'  = pct_base_target_50x, 
                'TARGET COVERAGE 100X' = pct_base_target_100x
            ) %>% 
            kbl( escape=FALSE, row.names=FALSE, align=c('l','c','r',rep('c', 8)) ) %>%
            kable_styling(bootstrap_options="condensed", full_width=TRUE, font_size=10, position="center") %>%
            row_spec(0, align='c', bold=T, background = "#82B366", color="#FFFFFF") %>% 
            column_spec(c(1),    width_min='15em', width='25em') %>% 
            column_spec(c(2),    width_min='2em',  width='2em' ) %>% 
            column_spec(c(3),    width_min='5em',  width='7em' ) %>% 
            column_spec(c(4:6),  width_min='2em',  width='3em' ) %>%
            column_spec(c(7),    width_min='6em',  width='7em' ) %>%
            column_spec(c(8:11), width_min='3em',  width='3em' ) 

        for( k in 1:nrow(SampleGroupIndex) )
        {
            QcTableOutput <- QcTableOutput %>% pack_rows( 
                group_label = SampleGroupIndex[k, "group"],  
                start_row   = SampleGroupIndex[k, "start"], 
                end_row     = SampleGroupIndex[k, "end"  ], 
                indent      = TRUE,
                background  = '#D5E8D4',
                label_row_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #82B366;"
            )
        }

        if( !is.null(LowQualSeqID) )
        {
            QcTableOutput <- QcTableOutput %>%
                column_spec(1, color=sapply(InputData$seq_id, function(y) ifelse( y %in% LowQualSeqID, "#ff0000", "#000000"))) 
        }

        return(QcTableOutput)
    }
#
#---/ wes alignment stats barplot /-----------------------------------------------------------------
    WES_analysis.Alignment.Stats.Barplot <- function( batchID_List=SeqFolderID, SampleInfo=sinfo )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("scales"))
        suppressPackageStartupMessages(library("RMySQL"))
        # seq folder id check --------------------------------------------------
        if( is.null(batchID_List) ){ stop("|---!!! No seq folder id found. REQUIRD. please check again. STOPPED.")}
        #-----------------------------------------------------------------------
        StatData <- NGS_report.Load.Stats( batchID_List =batchID_List, Src="report" )
        #-----------------------------------------------------------------------
        if( !is.null(SampleInfo) )
        { 
            StatData <- StatData %>% filter( seq_id %in% SampleInfo$seq_id )
            StatData$seq_id <- SampleInfo[match(StatData$seq_id, SampleInfo$seq_id), "sample_name"] 
            StatData$seq_id <- factor( StatData$seq_id, levels=SampleInfo$sample_name[length(SampleInfo$sample_name):1] )
        }else{
            StatData$seq_id <- factor( StatData$seq_id, levels=unique(StatData$seq_id)[length(unique(StatData$seq_id)):1] )
        }
        #-----------------------------------------------------------------------
        align_qc <- StatData %>% 
            dplyr::select(c("seq_id","rds_pf","rds_pf_align","rds_pf_align_hq")) %>% 
            mutate(rds_unalign=rds_pf-rds_pf_align, rds_align_nohq=rds_pf_align-rds_pf_align_hq) %>%
            mutate(pct_align_hq = round(rds_pf_align_hq/rds_pf*100, 1), pct_unalign = round( rds_unalign/rds_pf*100, 1 ))

        align_qc2 <- rbind(
            align_qc %>% 
                dplyr::select(c(1,2,4,7)) %>% 
                mutate(Read_Group="Aligned PF & HQ", rds_pf=rds_pf_align_hq/2) %>% 
                dplyr::rename(txt_x=rds_pf, reads=rds_pf_align_hq, pct=pct_align_hq) %>% 
                mutate(pct=paste0(pct,"%")),
            align_qc %>% 
                dplyr::select(c(1,2,6)) %>% mutate(pct="", Read_Group="Aligned PF", rds_pf=0) %>% 
                dplyr::rename(txt_x=rds_pf, reads=rds_align_nohq),
            align_qc %>% 
                dplyr::select(c(1,2,5)) %>% 
                mutate(pct="", Read_Group="Unaligned", rds_pf=0) %>% 
                dplyr::rename(txt_x=rds_pf, reads=rds_unalign)
        ) %>% 
            mutate(Read_Group=factor(Read_Group, levels=c("Unaligned","Aligned PF","Aligned PF & HQ"))) 
        
        align.plot <- ggplot(align_qc2, aes(x = reads, y = seq_id, fill=Read_Group)) + 
            geom_bar(position="stack", stat="identity", width =0.7) +
            labs( x="Reads", y="") + 
            theme(
                axis.text.x  = element_text(colour="black", size=9, hjust=0.5),
                axis.text.y  = element_text(family='Nunito',colour="black", size=8, hjust=1),
                axis.title.x = element_text(family='Nunito',face='bold', colour="black",size=10 , vjust=1),
                axis.title.y = element_text(family='Nunito',face='bold', colour="black",size=10 , vjust=1),
                axis.line.x  = element_line(colour = "#6f6d6d"), 
                axis.line.y  = element_line(colour = "#6f6d6d"),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'),
                legend.text  = element_text(family='Nunito'),
                legend.title = element_text(family='Nunito'),
                legend.position = "bottom"
            ) +
            geom_text(aes(x=txt_x, y=seq_id, label = pct),size=3, color=c('black'), hjust = 0.25) + 
            scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
            scale_fill_manual(values=c("Unaligned"='#B85450', "Aligned PF"='#D79B00',"Aligned PF & HQ"='#82B366' )) +
            coord_cartesian(xlim = c(0, max(align_qc2$reads)*1.2)) 

        return(align.plot)
    }
#
#---/ wes target coverage plot /--------------------------------------------------------------------
    WES_analysis.Target.Coverage.Plot <- function( batchID_List=SeqFolderID, SampleInfo=sinfo, X_LIMIT=NULL )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("scales"))
        suppressPackageStartupMessages(library("RMySQL"))
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        # seq folder id check --------------------------------------------------
        if( is.null(batchID_List) ){ stop("|---!!! No seq folder id found. REQUIRD. please check again. STOPPED.")}
        #-----------------------------------------------------------------------
        StatData <- NGS_report.Load.Stats( batchID_List =batchID_List, Src="report" )
        #-----------------------------------------------------------------------
        if( !is.null(SampleInfo) )
        { 
            StatData <- StatData %>% filter( seq_id %in% SampleInfo$seq_id )
            StatData$seq_id <- SampleInfo[match(StatData$seq_id, SampleInfo$seq_id), "sample_name"] 
            StatData$seq_id <- factor( StatData$seq_id, levels=SampleInfo$sample_name[length(SampleInfo$sample_name):1] )
        }else{
            StatData$seq_id <- factor( StatData$seq_id, levels=unique(StatData$seq_id)[length(unique(StatData$seq_id)):1] )
        }
        #-----------------------------------------------------------------------
        StatData$seq_id <- as.character(StatData$seq_id)
        INDEX           <- indexing( StatData$seq_id )
        
        mosd_cov <- ldply(lapply(INDEX, function(y) {
            df = StatData %>% filter( seq_id == y )
            data.frame( 
                depth = as.numeric(unlist(strsplit(df$mosdepth_cov_x_depth, ";" ))),
                cov   = as.numeric(unlist(strsplit(df$mosdepth_cov_y_cov,   ";" )))
            ) %>% 
            mutate( seq_id = y )
        })) %>% 
            filter( depth != 0 ) 

        cov.plot <- ggplot(mosd_cov) + 
            geom_point( aes(x=depth, y=cov ), size=0, color="white" ) +
            stat_smooth( aes(x=depth, y=cov, group=seq_id, color=seq_id), method = "lm", formula = y ~ poly(x, 21), se = FALSE, alpha=0.8 ) + 
            labs( x = "Coverage Depth", y = "Target Percentage") + 
            theme(
                axis.text.x  = element_text(family='Nunito', colour="black", size=12, hjust=0.5),
                axis.text.y  = element_text(family='Nunito', colour="black", size=12, hjust=0.5),
                axis.title.x = element_text(family='Nunito', face='bold', colour="black",size=14 , vjust=1),
                axis.title.y = element_text(family='Nunito', face='bold', colour="black",size=14 , vjust=1),
                axis.line.x  = element_line(colour = "#6f6d6d"), 
                axis.line.y  = element_line(colour = "#6f6d6d"),
                plot.margin  = unit(c(1, 1, 1, 1), "cm"),
                legend.position  = "bottom", 
                legend.text      = element_text(family='Roboto Condensed',size=10),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted')
            )+
            scale_color_discrete(name="") 
        
        if( !is.null(X_LIMIT) ){ cov.plot <- cov.plot +  coord_cartesian(xlim=c(0, X_LIMIT)) }

        return(cov.plot)
    }
#
#---/ variant stat data re-format /-----------------------------------------------------------------
    WES_report.Reformat.Variant.Stats <- function( VariantStatData, VariantGroup, SampleInfo=NULL, VariantCallMethod=NULL )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("reshape2"))
        #-----------------------------------------------------------------------
        if( is.null(SampleInfo) ){      
               VariantStatData$SampleName <- VariantStatData$seq_id
        }else{ VariantStatData$SampleName <- SampleInfo[match(VariantStatData$seq_id, SampleInfo$seq_id), "sample_name"] }
        #-----------------------------------------------------------------------
        if(       VariantCallMethod == "tumor_only"                  ){ VarCallMethod <- "Tonly" 
        }else if( VariantCallMethod == "match_normal_tissue"         ){ VarCallMethod <- "TS_N" 
        }else if( VariantCallMethod == "match_normal_organoid"       ){ VarCallMethod <- "ORG_N"
        }else if( VariantCallMethod == "tumor_only_gerline_filtered" ){ VarCallMethod <- "TonlyGermlineFiltered" }
        #-----------------------------------------------------------------------
        VS = VariantStatData %>%
            filter( variant_group == VariantGroup ) %>% 
            filter( variant_call_mode == VarCallMethod ) %>% 
            group_by( SampleName, Variant_Classification ) %>% 
            reframe( Variants = sum(stats) ) %>% data.frame() %>%
            reshape(., idvar = "SampleName", timevar = "Variant_Classification", direction = "wide" )
        #-----------------------------------------------------------------------
        if( VariantGroup %in% c("non_synonymous","filtered_somatic_only") )
        {
            ColumnList <- c( "Variants.Missense_Mutation","Variants.Nonsense_Mutation","Variants.Nonstop_Mutation", 
                             "Variants.Frame_Shift_Ins",  "Variants.Frame_Shift_Del",  "Variants.In_Frame_Ins", 
                             "Variants.In_Frame_Del",     "Variants.Splice_Site",      "Variants.Translation_Start_Site" )
            missingCols <- ColumnList[ ColumnList %nin% colnames(VS) ]
            #-------------------------------------------------------------------
            if( length(missingCols) > 0 ) 
            {
                extraStats <- as.data.frame(matrix( rep( NA, length(missingCols)*nrow(VS) ), ncol=length(missingCols) ))
                colnames(extraStats) <- missingCols
                VS <- cbind(VS, extraStats)
            }   
            #-------------------------------------------------------------------
            colnames(VS) = gsub("Variants.", ""  , colnames(VS))
            # VS <- VS %>%
            #     dplyr::rename( 
            #         'Missense'              = Variants.Missense_Mutation, 
            #         'Nonsense'              = Variants.Nonsense_Mutation, 
            #         'Nonstop'               = Variants.Nonstop_Mutation, 
            #         'FrameShift INS'        = Variants.Frame_Shift_Ins, 
            #         'FrameShift DEL'        = Variants.Frame_Shift_Del, 
            #         'InFrame INS'           = Variants.In_Frame_Ins, 
            #         'InFrame DEL'           = Variants.In_Frame_Del, 
            #         'SpliceSite'            = Variants.Splice_Site, 
            #         'Translation StartSite' = Variants.Translation_Start_Site
            #     ) %>%
            #     dplyr::select( c(
            #         "SampleName","Missense","Nonsense","NonStop",
            #         "FrameShift INS","FrameShift DEL","InFrame INS","InFrame DEL",
            #         "SpliceSite","Translation StartSite"
            #     ) )
        }else{
            colnames(VS) = gsub("Variants.", ""  , colnames(VS))
            colnames(VS) = gsub("_"        , " " , colnames(VS))
        }
        VS[is.na(VS)] = 0 
        return(VS)
    }
#
#---/ variant stat summary table /------------------------------------------------------------------
    WES_report.Variant.Stat.Summary.Table <- function( VariantStatTable, SampleOrder=NULL, SampleInfo=NULL, 
        StatColumns=statColumnNames
    )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        # sample order ---------------------------------------------------------
        if( is.null(SampleOrder) )
        {
            SampleOrderTable <- rbind(
                SampleInfo[which(SampleInfo$sample_group != ""), c("seq_id","sample_group","sample_info","sample_name")] %>% 
                    arrange( sample_info, sample_group, seq_id ),
                SampleInfo[which(SampleInfo$sample_group == ""),  c("seq_id","sample_group","sample_info","sample_name")] %>% 
                    arrange(sample_info, seq_id) %>% mutate(sample_group = "No Grouped")
            )
        }else{
            if( length(grep("^WES_", SampleOrder)) == nrow(SampleInfo) )
            {
                SampleOrderTable <- SampleInfo[match(SampleOrder, SampleInfo$seq_id), c("seq_id","sample_group","sample_info","sample_name")]
            }else{
                SampleOrderTable <- SampleInfo[match(SampleOrder, SampleInfo$sample_name), c("seq_id","sample_group","sample_info","sample_name")]
            }
        }
        # sample group index ---------------------------------------------------
        SampleGroupIndex <- ldply( lapply( indexing( SampleOrderTable$sample_group), function(sg)
            data.frame( group=sg, start=head(which( SampleOrderTable$sample_group == sg ), 1), end=tail(which( SampleOrderTable$sample_group == sg ), 1) )
        ))
        # row order ------------------------------------------------------------
        if( length(grep("^WES_", VariantStatTable$SampleName)) == nrow(VariantStatTable) )
        {      VariantStatTable <- VariantStatTable[match(SampleOrderTable$seq_id,      VariantStatTable$SampleName), ]
        }else{ VariantStatTable <- VariantStatTable[match(SampleOrderTable$sample_name, VariantStatTable$SampleName), ] }
        # column order ---------------------------------------------------------
        colnames(VariantStatTable)[1] <- "sample_name"
        colnames(VariantStatTable)    <- StatColumns[colnames(VariantStatTable)]
        VariantStatTable <- VariantStatTable %>% relocate( order(factor(colnames(VariantStatTable), levels=StatColumns)) )
        #-----------------------------------------------------------------------
        VarSumTableOutput <- VariantStatTable %>% 
            kbl( escape=FALSE, row.names=FALSE, align=c("l", rep("c", (ncol(VariantStatTable)-1) )) ) %>%
            kable_styling(bootstrap_options = "condensed", full_width=TRUE, font_size=10, position="center") %>%
            row_spec(0, align = 'c', bold = TRUE, background = "#82B366", color="#FFFFFF" ) %>%
            column_spec(1, width_min="10em", width="10em")
        #-----------------------------------------------------------------------
        for( k in 1:nrow(SampleGroupIndex) )
        {
            VarSumTableOutput <- VarSumTableOutput %>%
                pack_rows( 
                    group_label=SampleGroupIndex[k, "group"],  start_row=SampleGroupIndex[k, "start"], end_row=SampleGroupIndex[k, "end"], 
                    label_row_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #82B366;",
                    background = "#D5E8D4"
                )
        }
        return(VarSumTableOutput)
    }
#
#---/ preset gene title /---------------------------------------------------------------------------
    WES_report.Preset.Gene.Variants.Title.Table <- function( CancerCode=NULL, PresetName=NULL, VariantCallMethod=NULL )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        PresetVarTitle <- data.frame(matrix(c("CANCER CODE", CancerCode, "PRESET", PresetName, "VARIANT CALL METHOD", " ", VariantCallMethod), byrow=T, nrow=1)) %>%
            kbl( escape=FALSE, row.names=FALSE, align=c( rep('c',6), 'l'), col.names=NULL ) %>%
            kable_styling(bootstrap_options = "condensed", full_width=TRUE, font_size=10) %>%
            row_spec(1, bold=TRUE) %>%
            column_spec(c(1), width_min="4em", width="4em", background="#D5E8D4", border_left=TRUE, 
                extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #D5E8D4"  ) %>%
            column_spec(c(3),   width_min="3em",  width="3em",  background="#D5E8D4", extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #D5E8D4;" ) %>%
            column_spec(c(5),   width_min="8em",  width="9em",  background="#D5E8D4", extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #D5E8D4;" ) %>%
            column_spec(c(6),   width_min="1em",  width="1em",  extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #D5E8D4;") %>%
            column_spec(c(2,4), width_min="5em",  width="5em",  extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #D5E8D4;") %>%
            column_spec(c(7)  , width_min="10em", width="15em", border_right=TRUE, extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #D5E8D4;" )
        return(PresetVarTitle)
    }
#
#---/ preset gene variant table /-------------------------------------------------------------------
    WES_analysis.PresetGeneVariantsTable <- function( VariantData, CancerCode, SampleInfo, VariantCallMethod, PresetGeneList )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        cancer_sinfo <- SampleInfo %>% filter( sample_info == CancerCode )
        #-----------------------------------------------------------------------
        if(       VariantCallMethod == "tumor_only"                  ){ VarCallMethod <- "Tonly"
        }else if( VariantCallMethod == "match_normal_tissue"         ){ VarCallMethod <- "TS_N" 
        }else if( VariantCallMethod == "match_normal_organoid"       ){ VarCallMethod <- "ORG_N"
        }else if( VariantCallMethod == "tumor_only_gerline_filtered" ){ VarCallMethod <- "TonlyGermlineFiltered" }
        #-----------------------------------------------------------------------
        cancer_sinfo[which(cancer_sinfo$sample_group == ""), "sample_group"] <- "No Grouped"
        cancer_sinfo_group <- sort(unique(cancer_sinfo$sample_group))
        if( cancer_sinfo_group[[1]] == "No Grouped" ){ cancer_sinfo_group = cancer_sinfo_group[ c(2:length(cancer_sinfo_group), 1) ]}
        #-----------------------------------------------------------------------
        var_data           <- data.frame()
        for( k in 1:length(cancer_sinfo_group) )
        {
            group_sample_list <- cancer_sinfo %>% filter( sample_group == cancer_sinfo_group[k] ) %>% arrange( seq_id )
            var_data <- rbind(
                var_data, 
                VariantData %>% filter( 
                    variant_call_mode ==  VarCallMethod, 
                    seq_id           %in% group_sample_list$seq_id,
                    HGNC_SYMBOL      %in% PresetGeneList
                ) %>% arrange( seq_id, Chromosome, Start_Position)
            )
        }
        #-----------------------------------------------------------------------
        if( nrow(var_data) > 0 )
        {
            preset_var_res <- var_data %>% mutate(
                sample_name  = cancer_sinfo[match(seq_id, cancer_sinfo$seq_id), "sample_name"], 
                sample_group = cancer_sinfo[match(seq_id, cancer_sinfo$seq_id), "sample_group"],
                varDNA       = paste(Refseq_TRANSCRIPT,HGVSc, sep=":")
            ) 
        }else{
            preset_var_res <- data.frame()
        }
        #-----------------------------------------------------------------------
        return(preset_var_res)
    }
#
#---/ preset gene variant table (report) /----------------------------------------------------------
    WES_report.PresetGeneVariantsTable <- function( PresetGeneVarTable )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        sample_group_names <- unique(PresetGeneVarTable$sample_group)

        TableOutput <- PresetGeneVarTable %>% dplyr::select(
                c("sample_name","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")
            ) %>% dplyr::rename(
                'SAMPLE ID'=sample_name, 'GENE'= HGNC_SYMBOL, 'DNA CHANGE'=varDNA, 'PROTEIN CHANGE'=HGVSp_Short, 'TUMOR ALT COUNT'=TUMOR_ALT_COUNT
            ) %>% kbl( align='c', row.names = FALSE, escape = FALSE ) %>%
            kable_styling(bootstrap_options='condensed', full_width=TRUE, font_size=10 ) %>%
            row_spec(0, background="#82B366", color="#FFFFFF", align='c') %>%
            row_spec(1:nrow(PresetGeneVarTable), extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #cccccc;" ) %>% 
            column_spec(1,        width_min="15em", width="17em", bold=TRUE ) %>%
            column_spec(c(2,5,7), width_min="3em",  width="4em"  ) %>%
            column_spec(c(3,4),   width_min="6em",  width="6em"  ) %>%
            column_spec(c(6),     width_min="4em",  width="5em"  ) %>%
            collapse_rows(1)

        for( SGN in sample_group_names )
        {
            TableOutput <- TableOutput %>% pack_rows( 
                group_label = SGN, 
                start_row   = head(which(PresetGeneVarTable$sample_group == SGN), 1), 
                end_row     = tail(which(PresetGeneVarTable$sample_group == SGN), 1), 
                indent      = FALSE,
                label_row_css = "border-bottom: 1px solid; border-color: #82B366;",
                background = "#D5E8D4"
            ) 
        }
        
        return(TableOutput)
    }
#
#---/ preset gene variant heatmap /-----------------------------------------------------------------    
    WES_analysis.Preset.Gene.Variants.Heatmap <- function( PresetVariantsTable, CancerSampleInfo )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("ComplexHeatmap"))
        suppressPackageStartupMessages(library("circlize"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------
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
        #-----------------------------------------------------------------------
        var_mat <- PresetVariantsTable %>% group_by( HGNC_SYMBOL, sample_name ) %>%
            reframe( Variant_Class = paste(unique(Variant_Classification), collapse=";") ) %>% as.data.frame()
        multiVars <- grep(";", var_mat$Variant_Class)
        if( length(multiVars) > 0 ){ var_mat[ multiVars, "Variant_Class" ] <- "Multi_Variants" }
        #-----------------------------------------------------------------------
        VarMat <- reshape(var_mat , idvar = "HGNC_SYMBOL", timevar = "sample_name", direction = "wide")
        colnames(VarMat) <- gsub("^Variant_Class.", "", colnames(VarMat))
        #-----------------------------------------------------------------------
        MissIDs <- CancerSampleInfo[which(CancerSampleInfo$sample_name %nin% colnames(VarMat)), "sample_name" ]
        if( length(MissIDs) > 0 )
        {
            addMat <- matrix( rep(NA, length(MissIDs)*nrow(VarMat)), ncol = length(MissIDs) )
            colnames(addMat) <- MissIDs
            VarMat <- cbind(VarMat, addMat)
            VarMat <- VarMat[, c("HGNC_SYMBOL", CancerSampleInfo$sample_name)]
        }
        #-----------------------------------------------------------------------
        rownames(VarMat) <- VarMat$HGNC_SYMBOL
        #-----------------------------------------------------------------------
        HM <- Heatmap( as.matrix(data.frame(VarMat[,-1])), 
            cluster_rows = TRUE, cluster_columns = FALSE, 
            col = varColors, na_col="#ffffff", rect_gp = gpar(col = "#5b5b5b"),
            row_names_gp = gpar(fontsize = 11), column_names_gp = gpar(fontsize = 10), 
            row_names_side = "left",
            width = ncol(VarMat[,-1])*unit(6, 'mm'), height = nrow(VarMat[,-1])*unit(6, 'mm'),
            column_names_rot = 90, column_names_side = "top",
            heatmap_legend_param = list( title='VARIANT CLASS', title_gp=gpar(fontsize=9),labels_gp=gpar(fontsize=9), legend_width=unit(70,"mm") )
        )
        #-----------------------------------------------------------------------
        return(HM)
    }
#
#---/ calculate tanimoto /--------------------------------------------------------------------------
    WES_analysis.Get.Tanimoto <- function( SET1, SET2 )
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
#
#---/ calculate common ratio /----------------------------------------------------------------------
    WES_analysis.Get.Common.Rate <- function( SET1, SET2, p = 1 )
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
#
#---/ variant concordance venn-diagram /------------------------------------------------------------
    WES_analysis.Variant.Concordance.VennDiagram <- function( VariantList )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("ggVennDiagram"))
        #-----------------------------------------------------------------------
        vennFaceColoss <- c(
            "1"="#80B1D3",   "2"="#F8CECC", "1/2"="#82B366",
            "3"="#CCEBC5", "1/3"="#B85450", "2/3"="#D79B00", "1/2/3"="#BC80BD",
            "4"="#E1D5E7", "1/4"="#66C2A5", "2/4"="#FCCDE5",   "3/4"="#A6D854",
            "1/2/4" = "#8DA0CB", "1/3/4"="#377EB8","2/3/4"="#FF7F00","1/2/3/4"="#FFED6F"
        ) 
        SetLabelColors <- c( "TS1"="#B85450", "TS3"="#006EAF", "ORG1"="#2D7600", "ORG2"="#cf7e3b", "ORG0"="#2d2d2d" )
        #-----------------------------------------------------------------------
        VN    <- ggVennDiagram::Venn(VariantList)
        Vdata <- process_data(VN)

        valueLabel <- venn_regionlabel(Vdata)
        valueLabel <- valueLabel %>% mutate(pct = round(count/sum(count)*100, 1)) %>% mutate(lbl = paste0(count,"\n(", pct, "%)"))

        setLabel <- venn_setlabel(Vdata)
        if( nrow(setLabel) == 4 ){ setLabel$X = c(0.15,0.3, 0.6, 0.82) }
        if( nrow(setLabel) == 3 ){ setLabel$X = c(  -2,  5,   2 ) }
        if( nrow(setLabel) == 2 ){ setLabel$X = c(3,3) ; setLabel$Y = c(-1.5,5.5) }

        VP <- ggplot() +
            # Circle Region
            geom_polygon( data = venn_regionedge(Vdata), aes(X, Y, fill = id, group = id), show.legend = FALSE, alpha=0.3 ) +
            scale_fill_manual(values=vennFaceColoss) +
            # Edge
            geom_path(data = venn_setedge(Vdata), aes(X, Y, group = id), color="#FFFFFF", linewidth = 1, show.legend = FALSE) +
            # SetLabelbel in bold
            geom_text(data = setLabel, aes(X, Y, label = name), fontface = "bold", size = 4, hjust=0.5, show.legend = FALSE) +
            # Values
            geom_label(data = valueLabel, aes(X, Y, label = lbl), fontface = "bold", size=3.5, alpha = 0.5, fill="white" ) +
            theme_void() 
        #-----------------------------------------------------------------------
            if( length(VariantList) == 2 ){ 
                VP        <- VP + coord_flip() 
                EmptyPlot <- ggplot() + theme_void()
                VennPlot  <- cowplot::plot_grid(EmptyPlot, VP, EmptyPlot, ncol=1, rel_heights =c(0.2, 0.6, 0.2) )
            }else{
                VennPlot  <- VP
            }
        #-----------------------------------------------------------------------
        return(VennPlot)
    }
#
#---/ variant calss compare barplot /---------------------------------------------------------------
    WES_analysis.Variant.Class.Compare.Barplot <- function( VariantClassCompareData )
    {   
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("ggalluvial"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------
        color.Variants <- c( 
            "Missense_Mutation" = "#4DAF4A", "Nonsense_Mutation"      = "#E41A1C", "Frame_Shift_Ins" = "#FF7F00", 
            "Frame_Shift_Del"   = "#8DD3C7", "In_Frame_Ins"           = "#984EA3", "In_Frame_Del"    = "#377EB8", 
            "Nonstop_Mutation"  = "#A65628", "Translation_Start_Site" = "#FFFF33", "Splice_Site"     = "#F781BF", 
            "Multi_Variants"    = "#4c4c4c" 
        )
        #-----------------------------------------------------------------------
        VariantClassCompareData$sample_name <- gsub("_", "\n", VariantClassCompareData$sample_name)
        sampleNameFactor <- unique(VariantClassCompareData$sample_name)
        if( "TS" %in% sampleNameFactor )
        {
            sampleNameFactorLevel <- c("TS", sampleNameFactor[sampleNameFactor != "TS"] )
        }else{
            sampleNameFactorLevel <- sort(sampleNameFactor)
        }
        VariantClassCompareData$sample_name <- factor( VariantClassCompareData$sample_name, levels = sampleNameFactorLevel )
        #-----------------------------------------------------------------------
        VarClassBarplot <- VariantClassCompareData %>% ggplot(aes(x=sample_name, y=varN_PCT)) +
            geom_flow( aes(alluvium = Variant_Classification), alpha=0.8, linetype="dotted", fill="#FFFFFF", color="#c1c1c1", curve_type="linear", width=0.5 ) +
            geom_col( aes(fill = Variant_Classification), width=0.5, color="#FFFFFF" ) +
            scale_fill_manual( values = color.Variants ) +
            scale_y_continuous( expand=c(0,0)) +
            geom_hline( yintercept = 0, col = "#3d3d3d" ) +
            labs( x = "", y = "Variant Class (%)") +
            theme(
                axis.line.y  = element_line(colour = "#3d3d3d"),
                axis.text.x  = element_text(size = 10, vjust=0.5, hjust=0.5, face="bold"),
                axis.text.y  = element_text(size = 10, face="bold"),
                axis.title.y = element_text(size = 12, margin = margin(r = 3)),
                axis.ticks.x = element_blank(),
                legend.text  = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(4, 'mm'),
                panel.background = element_rect(fill="#FFFFFF"), 
                panel.grid.major = element_blank(), 
                plot.margin  = unit(c(1, 0.5, 1, 0.2), "cm")
            )
        return(VarClassBarplot)
    }
#
#---/ genomic change barplot /----------------------------------------------------------------------
    WES_analysis.Genomic.Change.Barplot <- function( GenomicChangeData )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------
        color.AlleleChange <- c( 
            "A > C" = "#4DAF4A", "A > G" = "#E41A1C", "A > T" = "#FF7F00", 
            "C > A" = "#8DD3C7", "C > G" = "#984EA3", "C > T" = "#377EB8", 
            "G > A" = "#A65628", "G > C" = "#FFFF33", "G > T" = "#F781BF", 
            "T > A" = "#4c4c4c", "T > C" = "#F8CECC", "T > G" = "#E1D5E7"
        )
        #-----------------------------------------------------------------------
        GenomicChangeData$sample_name <- gsub("_", "\n", GenomicChangeData$sample_name)
        sampleNameFactor <- unique(GenomicChangeData$sample_name)
        if( "TS" %in% sampleNameFactor )
        {
            sampleNameFactorLevel <- c("TS", sampleNameFactor[sampleNameFactor != "TS"] )
        }else{
            sampleNameFactorLevel <- sort(sampleNameFactor)
        }
        GenomicChangeData$sample_name <- factor( GenomicChangeData$sample_name, levels = sampleNameFactorLevel )
        #-----------------------------------------------------------------------
        GenomicChangePlot <- GenomicChangeData %>% ggplot( aes(x=sample_name, y=GDC_PCT) ) +
            geom_flow( aes(alluvium = VAC), alpha=0.9, linetype="dotted", fill="#FFFFFF", color="#c1c1c1", curve_type="linear", width=0.5 ) +
            geom_col( aes(fill = VAC), width=0.5, color="#FFFFFF") +
            scale_fill_manual(values=color.AlleleChange) +
            scale_y_continuous( expand=c(0,0)) +
            geom_hline(yintercept=0, col="#3d3d3d") +
            labs( x = "", y = "Genomic Change (%)", fill = "Genomic Change") +
            theme(
                axis.line.y  = element_line(colour = "#3d3d3d"),
                axis.text.x  = element_text(size = 10, vjust=0.5, hjust=0.5, face="bold"),
                axis.text.y  = element_text(size = 10, face="bold"),
                axis.title.y = element_text(size = 12, margin = margin(r = 3)),
                axis.ticks.x = element_blank(),
                legend.text  = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(4, 'mm'),
                panel.background = element_rect(fill="#FFFFFF"), 
                panel.grid.major = element_blank(), 
                plot.margin  = unit(c(1, 0.2, 1, 0.2), "cm")
            )  
        return(GenomicChangePlot)
    }
#
#---/ vaf correlation scatter plot /----------------------------------------------------------------
    WES_analysis.Vaf.Correlation.Plot <- function( VafCorData )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------
        VafCorData[is.na(VafCorData$X.vaf), "X.vaf"] = 0
        VafCorData[is.na(VafCorData$Y.vaf), "Y.vaf"] = 0

        X.only  <- VafCorData[which(!is.na(VafCorData$X.mut) &  is.na(VafCorData$Y.mut) & VafCorData$Y.vaf == 0 ), ]
        Y.only  <- VafCorData[which( is.na(VafCorData$X.mut) & !is.na(VafCorData$Y.mut) & VafCorData$X.vaf == 0 ), ]
        XY.both <- VafCorData[which(!is.na(VafCorData$X.mut) & !is.na(VafCorData$Y.mut)), ]

        X.label <- unique(VafCorData$X[!is.na(VafCorData$X)])
        Y.label <- unique(VafCorData$Y[!is.na(VafCorData$Y)])

        vafcorplot <- ggplot() + 
            theme(
                axis.line.x  = element_line(colour = "#6f6d6d"), 
                axis.line.y  = element_line(colour = "#6f6d6d"),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted')
            ) +
            geom_point(data=VafCorData, aes(x=X.vaf, y=Y.vaf), alpha=0.6, pch=19, size=2,   color="#eaeaea") +
            geom_point(data=X.only,     aes(x=X.vaf, y=Y.vaf), alpha=0.8, pch=21, size=3,   fill="#2D7600", color="#FFFFFF") +
            geom_point(data=Y.only,     aes(x=X.vaf, y=Y.vaf), alpha=0.8, pch=21, size=3,   fill="#006EAF", color="#FFFFFF") +
            geom_point(data=XY.both,    aes(x=X.vaf, y=Y.vaf), alpha=0.9, pch=21, size=4,     fill="#B85450", color="#FFFFFF") +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            labs( x=X.label, y=Y.label ) +
            theme(
                axis.text.x  = element_text(size = 9),
                axis.text.y  = element_text(size = 9),
                axis.title.x = element_text(face = "bold", size = 13),
                axis.title.y = element_text(face = "bold", size = 13),
                legend.position  = "none",
                plot.margin  = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
            )
        return(vafcorplot)
    }
#
#---/ oncoplot target gene /------------------------------------------------------------------------
    WES_analysis.Get.Oncoplot.Target.Genes <- function( 
        selectGroupName  = "all", 
        selectCancerType = "all", 
        VariantData, 
        topGenes         = 30, 
        geneSets         = NULL, 
        useBoth          = FALSE, 
        SampleInfo 
    )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        # sample group filter --------------------------------------------------
        if( selectGroupName == "all" ){ 
            hm_sinfo <- SampleInfo 
        }else{ 
            hm_SampleGroup <- unlist(strsplit(selectGroupName, ","))
            hm_sinfo <- SampleInfo %>% filter( sample_group %in% hm_SampleGroup) 
        }
        # cancer code filter ---------------------------------------------------
        if( selectCancerType == "all" ){ 
            hm_sinfo <- hm_sinfo 
        }else{ 
            hm_CancerCode <- unlist(strsplit(selectCancerType, ","))
            hm_sinfo <- hm_sinfo %>% filter( sample_info %in% hm_CancerCode ) 
        }
        # oncoplot target sample id --------------------------------------------
        SelectSeqID <- hm_sinfo$seq_id
        # gene-variant sort table ----------------------------------------------
        opTargetGenes <- VariantData %>% filter( seq_id %in% SelectSeqID ) %>% 
            group_by( HGNC_SYMBOL ) %>% reframe( varCount=length(HGNC_SYMBOL) ) %>% arrange( desc(varCount) )
        # gene set filter ------------------------------------------------------            
        if( is.null(geneSets) )
        {
            if( is.null(topGenes) )
            {
                oncoPlotGenes <- opTargetGenes %>% .$HGNC_SYMBOL %>% table() %>% sort() %>% rev()   
            }else{
                topGeneList   <- opTargetGenes %>% filter( varCount >= opTargetGenes$varCount[topGenes]) %>% .$HGNC_SYMBOL
                if( length(topGeneList) > 30 )
                { 
                    topGeneList   <- opTargetGenes %>% filter( varCount >= opTargetGenes$varCount[topGenes] + 1 ) %>% .$HGNC_SYMBOL
                } 
                oncoPlotGenes <- opTargetGenes %>% filter( HGNC_SYMBOL %in% topGeneList ) %>% .$HGNC_SYMBOL %>% table() %>% sort() %>% rev()
            }
        }else{
            geneSetGenes <- opTargetGenes %>% filter( HGNC_SYMBOL %in% geneSets ) 
            if( useBoth )
            {
                if( nrow(geneSetGenes) >= topGenes )
                {
                    topGeneList   <- geneSetGenes %>% filter( varCount >= geneSetGenes$varCount[topGenes]) %>% .$HGNC_SYMBOL
                    oncoPlotGenes <- geneSetGenes %>% filter( HGNC_SYMBOL %in% topGeneList ) %>% .$HGNC_SYMBOL %>% table() %>% sort() %>% rev()
                }else{
                    oncoPlotGenes <- geneSetGenes %>% .$HGNC_SYMBOL %>% table() %>% sort() %>% rev()
                }
            }else{
                oncoPlotGenes <- geneSetGenes %>% .$HGNC_SYMBOL %>% table() %>% sort() %>% rev()
            }
        }
       
        oncoPlotGenesCount <- VariantData %>% 
            filter( seq_id %in% SelectSeqID, HGNC_SYMBOL %in% names(oncoPlotGenes) ) %>% 
            group_by(seq_id, HGNC_SYMBOL) %>% 
            reframe( mutations = paste(unique(Variant_Classification), collapse=";")) 
            
        oncoPlotGenesCount$mutations <- sapply(oncoPlotGenesCount$mutations, function(mut) ifelse( length(unlist(strsplit(mut, ";"))) > 1, "Multi_Variants", mut ))
            
        return(oncoPlotGenesCount)
    }
#
#---/ draw oncoplot /-------------------------------------------------------------------------------
    WES_analysis.Draw.OncoPlots <- function( 
        targetGeneCountTable, 
        SampleInfo, 
        SampleOrderByIdGroup = TRUE, 
        includeAllSamples    = FALSE 
    )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("ComplexHeatmap"))
        suppressPackageStartupMessages(library("circlize"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------

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

        if( includeAllSamples )
        {
            noVarSeqIDs <- SampleInfo[which(SampleInfo$seq_id %nin% targetGeneCountTable$seq_id), "seq_id"]
            #noVarSeqIDs <- colOrders$id[ colOrders$id %nin% colnames(oncoPlotMatrix) ]
    
            if( length(noVarSeqIDs) > 0 )
            {
                targetGeneCountTableExtra <- data.frame( seq_id = noVarSeqIDs, HGNC_SYMBOL = NA, mutations = NA)
                targetGeneCountTable <- rbind(as.data.frame(targetGeneCountTable), targetGeneCountTableExtra)

                # noVarMat           <- matrix( rep(NA, length(noVarSeqIDs)*nrow(oncoPlotMatrix)), ncol=length(noVarSeqIDs) )
                # colnames(noVarMat) <- noVarSeqIDs
                # rownames(noVarMat) <- rownames(oncoPlotMatrix)
                # oncoPlotMatrix     <- cbind(oncoPlotMatrix, noVarMat)
    
                # add_annotC <- data.frame(noVarSeqIDs, matrix(rep(0,length(noVarSeqIDs)*(ncol(annotC_data)-1)), nrow=length(noVarSeqIDs)))
                # colnames(add_annotC) <- colnames(annotC_data)
    
                # annotC_data           <- rbind( annotC_data, add_annotC )
                # rownames(annotC_data) <- annotC_data$seq_id
            }
        }
        
        annotR_data <- targetGeneCountTable %>% 
            group_by(HGNC_SYMBOL, mutations) %>% 
            reframe( varCount = length(seq_id) ) %>% 
            data.frame() %>% 
            reshape(., idvar="HGNC_SYMBOL", timevar="mutations", direction="wide") %>%
            filter( !is.na(HGNC_SYMBOL) ) 

        annotR_data[is.na(annotR_data)] <- 0 
        colnames(annotR_data) <- gsub("^varCount.", "", colnames(annotR_data))
        rownames(annotR_data) <- annotR_data[,1]
        annotR_data <- annotR_data[, c("HGNC_SYMBOL", names(varColors)[which(names(varColors) %in% colnames(annotR_data))]) ]
    
        annotC_data <- targetGeneCountTable %>% 
            group_by(seq_id, mutations) %>% 
            reframe(length(HGNC_SYMBOL)) %>% 
            data.frame() %>% 
            reshape(., idvar="seq_id", timevar="mutations", direction="wide")
    
        colnames(annotC_data) <- gsub("^length.HGNC_SYMBOL..", "", colnames(annotC_data))
        rownames(annotC_data) <- annotC_data[,1]
        annotC_data[is.na(annotC_data)] <- 0
        annotC_data <- annotC_data[, c("seq_id", names(varColors)[which(names(varColors) %in% colnames(annotC_data))]) ]
    
        oncoPlotMatrix <- targetGeneCountTable %>% data.frame() %>% reshape(., idvar="HGNC_SYMBOL", timevar="seq_id", direction="wide")
        oncoPlotMatrix <- oncoPlotMatrix[which(!is.na(oncoPlotMatrix$HGNC_SYMBOL)), ]

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
    

        sampleFreq <- targetGeneCountTable %>% 
            group_by(HGNC_SYMBOL) %>% 
            reframe( sampleFreq = round(length(seq_id)/(ncol(oncoPlotMatrix)-1) *100, 0) ) %>%
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
            col    = list( TISSUE = tissueColors, GROUP = groupColors, ORIGIN = c("TS" = "#1F78B4", "ORG" = "#A6CEE3", "BLOOD"="#F8CECC") ),
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
#
#---/ cumulative oncoplot data load from database /-------------------------------------------------
    WES_analysis.Get.Cumulative.Oncoplot.Data <- function( CancerCode, ClientID, VariantCallMethod, SampleLimit=10 )    
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("RMySQL"))
        #-----------------------------------------------------------------------
        source("/storage/home/kangsm/myScripts/Default_Scripts/ruo_wes_db.R")
        dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_info )
        cancercode_seqid_all <- dbGetQuery(dbCon, sprintf("SELECT seq_id from seqid_info WHERE sample_info = '%s' AND client_facility_id = '%s'", CancerCode, ClientID)) 
        dbDisconnect(dbCon)
        
        dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
        cancercode_vardata_seqid_all <- dbGetQuery(dbCon, sprintf("SELECT DISTINCT(seq_id) FROM variants_somatic WHERE variant_call_mode = '%s'", VariantCallMethod ))
        dbDisconnect(dbCon)

        oncoplot2_target_samples <- cancercode_vardata_seqid_all$seq_id[ cancercode_vardata_seqid_all$seq_id %in% cancercode_seqid_all$seq_id ]

        if( length(oncoplot2_target_samples) < SampleLimit )
        {
            res <- NULL
        }else{
            querySeqID <- paste(paste0("'", oncoplot2_target_samples, "'"), collapse=",")

            dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
            var_data <- dbGetQuery(dbCon, sprintf("SELECT seq_id,HGNC_SYMBOL,Variant_Classification FROM variants_somatic WHERE variant_call_mode = '%s' AND seq_id IN (%s)", VariantCallMethod, querySeqID))
            dbDisconnect(dbCon)

            dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_info )
            var_sinfo <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_id IN (%s)", querySeqID))
            dbDisconnect(dbCon)

            res <- list(var_sinfo = var_sinfo, var_data = var_data )
        }
        return(res)
    }
#
#---/ individual matching heatmap /-----------------------------------------------------------------
    WES_analysis.Indiv.Match.Heatmap <- function( IndivMatchResult, SampleInfo )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("ComplexHeatmap"))
        suppressPackageStartupMessages(library("RColorBrewer"))
        suppressPackageStartupMessages(library("circlize"))
        #-----------------------------------------------------------------------
        ColorPresets <- c(
            "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
            "#A65628", "#F781BF", "#999999", "#66C2A5", "#8DA0CB", "#FCCDE5",
            "#A6D854", "#FFD92F", "#E5C494", "#FFFFB3", "#FB8072", "#80B1D3",
            "#FDB462", "#B3DE69", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
        )
        #-----------------------------------------------------------------------
        grouped_sinfo   <- SampleInfo %>% filter( sample_group != "" ) 
        ungrouped_sinfo <- SampleInfo %>% filter( sample_group == "" ) 

        if( nrow(grouped_sinfo) > 0   ){ GROUPED_SINFO   <- grouped_sinfo %>% arrange( sample_info, sample_group, seq_id ) }else{ GROUPED_SINFO = data.frame() }
        if( nrow(ungrouped_sinfo) > 0 ){ UNGROUPED_SINFO <- ungrouped_sinfo %>% mutate( sample_group = "No Grouped" ) %>% arrange( sample_info, seq_id ) }else{ UNGROUPED_SINFO = data.frame() }
        
        im_sinfo <- as.data.frame(rbind( GROUPED_SINFO, UNGROUPED_SINFO ))
            
        #-----------------------------------------------------------------------
        IMR <- matrix( rep(0, length(im_sinfo$seq_id)*length(im_sinfo$seq_id)), ncol=length(im_sinfo$seq_id))
        rownames(IMR) = colnames(IMR) = im_sinfo$seq_id
        #-----------------------------------------------------------------------
        if( nrow(IndivMatchResult) < 1 )
        {
            for( i in 1:nrow(im_sinfo) ){ IMR[i,i] <- 1 }
        }else{
            for( i in 1:nrow(IndivMatchResult) )
            {
                IMR[ IndivMatchResult[i,3], IndivMatchResult[i,4] ] = IndivMatchResult[i,6]
                IMR[ IndivMatchResult[i,4], IndivMatchResult[i,3] ] = IndivMatchResult[i,6]
            }
            for( i in 1:nrow(im_sinfo) )
            {
                IMR[i, i] <- 1
            }
        }
        #-----------------------------------------------------------------------
        rownames(IMR) <- im_sinfo[match(rownames(IMR), im_sinfo$seq_id), "sample_name"]
        colnames(IMR) <- im_sinfo[match(colnames(IMR), im_sinfo$seq_id), "sample_name"]
        #-----------------------------------------------------------------------
        IndivMatchSampleGroup    <- unique(im_sinfo$sample_group)
        SampleGroupColors        <- ColorPresets[1:length(IndivMatchSampleGroup)]
        names(SampleGroupColors) <- IndivMatchSampleGroup
        #-----------------------------------------------------------------------
        if( ncol(IMR) > 20 ){ unitSize  <- 3  }else{ unitSize  <- 5   }
        if( ncol(IMR) > 20 ){ fontRatio <- 1  }else{ fontRatio <- 1.1 }
        #-----------------------------------------------------------------------
        annotC <- HeatmapAnnotation(
            GROUP                   = im_sinfo[match(colnames(IMR), im_sinfo$sample_name), 'sample_group'],
            col                     = list( 'GROUP' = SampleGroupColors ),
            gap                     = unit(1.4, "mm"),
            annotation_name_gp      = gpar(fontsize = 8*fontRatio), 
            annotation_name_side    = 'left', 
            annotation_legend_param = list( title_gp=gpar(fontsize=8*fontRatio), labels_gp=gpar(fontsize=8*fontRatio) )  
        )
        #-----------------------------------------------------------------------
        col_fun = colorRamp2( c(0,0.5,1), c( "#FFFFFF", "#cde0c1","#344728"), transparency=0 )
        #-----------------------------------------------------------------------
        
        IndivMatchHeatmap <- Heatmap( as.matrix(IMR), 
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
        return(IndivMatchHeatmap)
    }
#
#---/ individual matching table /-------------------------------------------------------------------
    WES_report.Indiv.Match.Table <- function( IndivMatchResult, SampleInfo )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        indiv_res <- IndivMatchResult %>% mutate(
            seqid_1 = SampleInfo[match(seqid_1, SampleInfo$seq_id), "sample_name" ],
            seqid_2 = SampleInfo[match(seqid_2, SampleInfo$seq_id), "sample_name" ]
        ) %>% mutate(
            group_1 = SampleInfo[match(seqid_1, SampleInfo$sample_name), "sample_group" ],
            group_2 = SampleInfo[match(seqid_2, SampleInfo$sample_name), "sample_group" ]
        ) %>% filter( group_1 == group_2 ) 

        sample_orders = lapply(unique(SampleInfo$sample_info), function(y) {
            SampleInfo %>% filter( sample_info == y ) %>% arrange(sample_group) %>% dplyr::select(c("sample_info","sample_group")) %>% unique()
        }) %>% ldply()
        
        INDIV_RES = data.frame()
        for( j in 1:nrow(sample_orders) ){
            INDIV_RES = rbind(INDIV_RES,
                indiv_res %>% filter( group_1 == sample_orders[j, "sample_group"], group_2 == sample_orders[j, "sample_group"]) 
            )
        }

        INDIV_RES$group_match  <- apply( INDIV_RES[c("group_1","group_2")], 1, function(y) 
        {
            if( all(y[1] == "" & y[2] == "") ) 
            {
                gm <- "unmatched"
            }else{
                gm <- ifelse( y[1] == y[2], "matched", "unmatched") 
            }
            return(gm)
        })

        indiv_res_table <- INDIV_RES %>% dplyr::select(
            c("seqid_1","seqid_2","group_1","group_2","cor","match_res","group_match")
        ) %>% dplyr::rename(
            'Sample-A'=seqid_1, 'Sample-B'=seqid_2, 'Sample-A Group'=group_1, 'Sample-B Group'=group_2, 'Correlation'=cor, 'Sample Matching'=match_res, 'Group Matching'=group_match
        ) %>% 
            kbl( escape=FALSE, align='c' ) %>% 
            kable_styling(bootstrap_options='condensed', full_width=TRUE, font_size=10 ) %>%
            row_spec(0, background="#9673A6", color="#FFFFFF") %>%
            column_spec(c(1,2), width_min="12em", width="15em" ) %>%
            column_spec(c(3,4), width_min="16em", width="8em" ) %>%
            column_spec(c(5),   width_min="3em", width="4em",
                background=ifelse( INDIV_RES$match_res == "matched", "#82B366", "#D79B00"), color="#FFFFFF"
            ) %>%
            column_spec(c(6),   width_min="4em", width="5m", 
                background=ifelse( INDIV_RES$match_res == "matched", "#D5E8D4", "#F8CECC"),
                color=ifelse( INDIV_RES$match_res == "matched", "#326a1a", "#a72e2e")
            ) %>%
            column_spec(c(7),   width_min="4em", width="5m", 
                background=ifelse( INDIV_RES$group_match == "matched", "#82B366", "#B85450"),
                color=ifelse( INDIV_RES$group_match == "matched", "#D5E8D4", "#F8CECC")
            ) 
        return(indiv_res_table)
    }
#
#---/ variant concordance result table /------------------------------------------------------------
    WES_report.Variant.Concordance.Table <- function( ConcordanceResultTable )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        ColNames <- c("SAMPLE-A","SAMPLE-B","SAMPLE-A","SAMPLE-B","COMMON","SAMPLE-A","SAMPLE-B","CONCORDANCE") 
        #-----------------------------------------------------------------------
        res <- ConcordanceResultTable[,-1] 
        OutputTable <- res %>% #dplyr::select(!c("Sample.Group")) %>%
            kbl( escape=FALSE, col.names=ColNames, align='c', row.names=FALSE ) %>%
            kable_styling(bootstrap_options='condensed', full_width=TRUE, font_size=10 ) %>%
            row_spec(0, align='c', bold=T, extra_css="border-top: 1px solid; border-bottom: 1px solid; border-color: #6C8EBF;" ) %>% 
            row_spec(c(1:nrow(res)), extra_css="border-top: 0px solid; border-bottom: 1px solid; border-color: #DAE8FC;" ) %>%
            column_spec(c(1,2),    width_min='10em', width='15em') %>% 
            column_spec(c(3,4,5), width_min='4em',  width='5em' ) %>%
            column_spec(c(6),      width_min='4em',  width='5em', background=ifelse(res$MatchRate.ID1 > 50, "#82B366", "#D79B00"), color="#FFFFFF" ) %>%
            column_spec(c(7),      width_min='4em',  width='5em', background=ifelse(res$MatchRate.ID2 > 50, "#82B366", "#D79B00"), color="#FFFFFF" ) %>%
            column_spec(c(8),      width_min='4em',  width='5em', background=ifelse(res$Concordance.Score > 0.5, "#82B366", "#D79B00"), color="#FFFFFF" ) %>%
            add_header_above(header = c(" "=2, "VARIANTS"=3, "MATCH RATE"=2,  " "=1), font_size=10, bold = TRUE, line_sep=15 ) 

        return(OutputTable)          
    }
#
#---/ vaf-plot title /-----------------------------------------------------------------------------
    WES_report.Vaf.Corr.Plot.Title.Table <- function( CancerCode=NULL, SampleGroup=NULL, VariantCallMethod=NULL )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        PresetVarTitle <- data.frame(matrix(c("CANCER CODE", CancerCode, "SAMPLE GROUP", SampleGroup, "VARIANT CALL METHOD", " ", VariantCallMethod), byrow=T, nrow=1)) %>%
            kbl( escape=FALSE, row.names=FALSE, align=c( rep('c',6), 'l'), col.names=NULL ) %>%
            kable_styling(bootstrap_options = "condensed", full_width=TRUE, font_size=10) %>%
            row_spec(1, bold=TRUE) %>%
            column_spec(c(1), width_min="4em", width="4em", background="#FFE6CC", border_left=TRUE, 
                extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #FFE6CC;"  ) %>%
            column_spec(c(3),   width_min="4em",  width="5em",  background="#FFE6CC", extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #FFE6CC;" ) %>%
            column_spec(c(5),   width_min="8em",  width="9em",  background="#FFE6CC", extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #FFE6CC;" ) %>%
            column_spec(c(6),   width_min="1em",  width="1em",  extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #FFE6CC;") %>%
            column_spec(c(2),   width_min="3em",  width="3em",  extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #FFE6CC;") %>%
            column_spec(c(4),   width_min="7em",  width="7em",  extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #FFE6CC;") %>%
            column_spec(c(7)  , width_min="10em", width="15em", border_right=TRUE, extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #FFE6CC;" )
        return(PresetVarTitle)
    }
#
#---/ hla-type table prepare /----------------------------------------------------------------------
    WES_report.HLA.Type.Table.Prepare <- function( HlaResult, SampleInfo, CancerCode )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("reshape2"))
        # hla sample info ------------------------------------------------------
        hla_sinfo <- rbind(
            SampleInfo %>% filter( sample_info == CancerCode, sample_group != "" ) %>% arrange(sample_group, seq_id),
            SampleInfo %>% filter( sample_info == CancerCode, sample_group == "" ) %>% mutate( sample_group = "No Grouped" ) %>% arrange(seq_id)
        )
        # hla result -----------------------------------------------------------
        hla_res       <- HlaResult %>% filter( seq_id %in% hla_sinfo$seq_id ) %>% dplyr::select(c('seq_id','locus','allele1_hla_la','allele2_hla_la'))
        hla_res$HLA   <- apply(hla_res[,c('allele1_hla_la','allele2_hla_la')], 1, function(y) paste(y, collapse="#") )
        hla_res_table <- hla_res %>% dplyr::select(c('seq_id','locus','HLA')) %>% 
            mutate( HLA = gsub("#","\\<br\\>", HLA) ) %>% 
            mutate( HLA = gsub("\\*","\\* ", HLA) ) %>% 
            filter( locus %nin% c("E","F","G") ) %>%
            reshape( idvar='seq_id', timevar='locus', direction='wide' ) 
        colnames(hla_res_table) <-  gsub("^HLA.D","D", colnames(hla_res_table))
        colnames(hla_res_table) <-  gsub("\\.","-", colnames(hla_res_table))

        hla_res_table <- hla_res_table[match(hla_sinfo$seq_id, hla_res_table$seq_id), ]
        
        # not matched sample ---------------------------------------------------
        no_matched_samples        <- HlaResult %>% filter( seq_id %in% hla_sinfo$seq_id, concordance == "No" ) 
        if( nrow(no_matched_samples) > 0 )
        {
            no_matched_samples$sample_name <- hla_sinfo[match(no_matched_samples$seq_id, hla_sinfo$seq_id), "sample_name"]
        }else{
            no_matched_samples <- data.frame()
        }

        # sample groups --------------------------------------------------------
        SampleGroupList <- unique(hla_sinfo$sample_group)

        # sample group index ---------------------------------------------------
        SampleGroupIndex <- ldply(lapply(indexing(SampleGroupList), function(sg)
        {
            data.frame(
                sample_group = sg, 
                start_row    = head(which( hla_res_table$seq_id %in% hla_sinfo[which(hla_sinfo$sample_group == sg), "seq_id"] ), 1),
                end_row      = tail(which( hla_res_table$seq_id %in% hla_sinfo[which(hla_sinfo$sample_group == sg), "seq_id"] ), 1)
            )
        }))
        
        # sample name modification ---------------------------------------------
        hla_res_table$sample_name <- hla_sinfo[match(hla_res_table$seq_id, hla_sinfo$seq_id), "sample_name"]
        for( sg in SampleGroupList[ SampleGroupList != "No Grouped"] )
        {
            sample_group_tag               <- paste0(gsub("\\-", "_", sg), "_")
            hla_res_table$sample_name      <- gsub(sample_group_tag, "", hla_res_table$sample_name)
            no_matched_samples$sample_name <- gsub(sample_group_tag, "", no_matched_samples$sample_name)
        }
        
        # result list 
        HLA_RES <- list(
            result_table       = hla_res_table,
            sample_group_index = SampleGroupIndex,
            no_match           = no_matched_samples
        )
        
        return(HLA_RES)
    }
#
#---/ hla-type result table /-----------------------------------------------------------------------
    WES_report.HLA_Type.Table <- function( HlaResultTable, HlaResultGroupIndex, HlaResultNoMatch )
    {    
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        OutputTable <- HlaResultTable %>% dplyr::select(c(12,2:11)) %>%
            dplyr::rename('SAMPLE NAME'=1) %>%
            kbl( escape = FALSE, align=c('c'), row.names=F ) %>% 
            kable_styling(bootstrap_options='condensed', full_width = TRUE, font_size=9 ) %>%
            row_spec(0, font_size = 11, bold=TRUE ) %>% 
            column_spec(c(1),    width_min="10em", width="12em" ) %>%
            column_spec(c(2:11), width_min="4em",  width="8em"  ) %>%
            add_header_above(c(" "=1,"MHC-I"=3, "MHC-II"=7), font_size=12, bold=TRUE, line_sep=15)
        #-----------------------------------------------------------------------
        for( k in 1:nrow(HlaResultGroupIndex) )
        {
            OutputTable <- OutputTable %>% pack_rows(
                group_label = HlaResultGroupIndex[k, "sample_group"],
                start_row   = HlaResultGroupIndex[k, "start_row"],
                end_row     = HlaResultGroupIndex[k, "end_row"],
                indent      = FALSE,
                background  = "#E1D5E7"
            )
        }
        #-----------------------------------------------------------------------
        if( nrow(HlaResultNoMatch) > 0 )
        {
            NoMatchLabelColors <- sapply( HlaResultTable$seq_id, function(w) ifelse( w %in% HlaResultNoMatch$seq_id , "#C73500", "#000000") )   
            OutputTable        <- OutputTable %>% column_spec(1, color=NoMatchLabelColors) 
        }
        #-----------------------------------------------------------------------
        return(OutputTable)
    }
#
#---/ no hla-match sample table /-------------------------------------------------------------------
    WES_report.HLA_Type.No.Match.Samples.Table <- function( HlaResultNoMatch, SampleInfo )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        no_match_table <- HlaResultNoMatch %>% mutate(
            seq_id = SampleInfo[match(seq_id, SampleInfo$seq_id), "sample_name"],
            sw1    = paste(allele1_hla_la,   allele2_hla_la,   sep="#"), 
            sw2    = paste(allele1_optitype, allele2_optitype, sep="#")
        ) %>% dplyr::select(c('seq_id','locus','sw1','sw2')) %>% 
            mutate( locus = gsub("^","HLA-",locus)    ) %>%
            mutate( sw1   = gsub("#","\\<br\\>", sw1) ) %>% 
            mutate( sw2   = gsub("#","\\<br\\>", sw2) ) %>%
            mutate( sw1   = gsub("\\*","\\* ", sw1)   ) %>% 
            mutate( sw2   = gsub("\\*","\\* ", sw2)   ) %>%
            dplyr::rename('SAMPLE NAME'=1, 'LOCUS'=2, 'HLA-LA'=3, 'OPTITYPE'=4) %>%
            kbl( escape = FALSE, align='c', row.names=F ) %>% 
            kable_styling(bootstrap_options="condensed", full_width = FALSE, font_size=10, position = 'left' ) %>%
            row_spec(0, align='center', font_size = 10, bold=TRUE, background = "#9673A6", color="#FFFFFF" ) %>%
            column_spec(1,   width_min="10em", width="12em" ) %>%
            column_spec(2,   width_min="3em",  width="4em"  ) %>%
            column_spec(3:4, width_min="12em", width="15em" ) 

        return(no_match_table)
    }
#
#---/ msi result table /----------------------------------------------------------------------------
    WES_report.MSI.Result.Table <- function( MsiResult, SampleInfo, CancerCode )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("reshape2"))
        # hla sample info ------------------------------------------------------
        msi_sinfo <- rbind(
            SampleInfo %>% filter( sample_info == CancerCode, sample_group != "" ) %>% arrange(sample_group, seq_id),
            SampleInfo %>% filter( sample_info == CancerCode, sample_group == "" ) %>% mutate( sample_group = "No Grouped" ) %>% arrange(seq_id)
        )
        # msi data -------------------------------------------------------------
        msi_data <- MsiResult %>% filter( msi_run_mode == "Matched-Normal", seq_id %in% msi_sinfo$seq_id )
        if( nrow(msi_data) > 0 )
        {
            msi_data <- rbind(
                msi_data,
                MsiResult %>% filter( msi_run_mode == "Tumor-Only", seq_id %in% msi_sinfo$seq_id ) %>% filter( seq_id %nin% msi_data$seq_id )
            )
        }else{
            msi_data <- MsiResult %>% filter( msi_run_mode == "Tumor-Only", seq_id %in% msi_sinfo$seq_id )
        }
        msi_data <- msi_data[match(msi_sinfo$seq_id, msi_data$seq_id), ]
        msi_data[which(msi_data$mantis_score == "-"), "mantis_score"] = ""
        # sample groups 
        SampleGroupList <- unique(msi_sinfo$sample_group)

        # sample group index ---------------------------------------------------
        SampleGroupIndex <- ldply(lapply(indexing(SampleGroupList), function(sg)
        {
            data.frame(
                sample_group = sg, 
                start_row    = head(which( msi_data$seq_id %in% msi_sinfo[which(msi_sinfo$sample_group == sg), "seq_id"] ), 1),
                end_row      = tail(which( msi_data$seq_id %in% msi_sinfo[which(msi_sinfo$sample_group == sg), "seq_id"] ), 1)
            )
        }))
        # sample name modification ---------------------------------------------
        msi_data$sample_name <- msi_sinfo[match(msi_data$seq_id, msi_sinfo$seq_id), "sample_name"]
        for( sg in SampleGroupList[ SampleGroupList != "No Grouped" ] )
        {
            sample_group_tag      <- paste0(gsub("\\-", "_", sg), "_")
            msi_data$sample_name  <- gsub(sample_group_tag, "", msi_data$sample_name)
        }
        # msi status column ----------------------------------------------------
        msi_data$status <- apply( msi_data[,c("msi_run_mode","msisensor2_msi_status","mantis_msi_status")], 1, function(z) ifelse( z[1] == "Matched-Normal", z[3], z[2] ) )

        # output table ---------------------------------------------------------
        StatusBackgroundColors <- sapply( msi_data$status, function(w) ifelse( w == "MSS", "#FFFFFF", "#B85450") )
        StatusLabelColors      <- sapply( msi_data$status, function(w) ifelse( w == "MSS", "#000000", "#FFFFFF") )
        sw1_label_colors       <- sapply( msi_data$msisensor2_score, function(w) ifelse( as.numeric(w) > 20,  "#C73500", "#000000") )
        sw2_label_colors       <- suppressWarnings(sapply( msi_data$mantis_score, function(w) ifelse( w == "", "#000000", ifelse( as.numeric(w) > 0.4, "#C73500", "#000000") )))

        OutputTable <- msi_data %>% dplyr::select(c("sample_name","msisensor2_score","mantis_score","status","msi_run_mode")) %>%
            dplyr::rename(
                'SAMPLE NAME'=sample_name, 'MSIsensor2'=msisensor2_score, 'MANTIS'=mantis_score, 'MSI STATUS'=status, 'VARIANT CALL'=msi_run_mode
            ) %>%
            kbl( escape = FALSE, align=c('c'), row.names=FALSE ) %>% 
            kable_styling(bootstrap_options='condensed', full_width = TRUE, font_size=10 ) %>%
            row_spec(0, font_size = 12, bold=TRUE, extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #D79B00;" ) %>% 
            column_spec(c(1), width_min="12em", width="15em" ) %>%
            column_spec(c(2), width_min="6em",  width="8em", color=sw1_label_colors  ) %>%
            column_spec(c(3), width_min="6em",  width="8em", color=sw2_label_colors  ) %>%
            column_spec(c(4), width_min="5em",  width="8em", background=StatusBackgroundColors, color=StatusLabelColors  ) %>%
            column_spec(c(5), width_min="5em",  width="8em"  ) %>%
            add_header_above(c(" "=1,"SCORE"=2, " "=2), font_size=12, bold=TRUE, line_sep=15)

        for( k in 1:nrow(SampleGroupIndex) )
        {
            OutputTable <- OutputTable %>% pack_rows(
                group_label = SampleGroupIndex[k, "sample_group"],
                start_row   = SampleGroupIndex[k, "start_row"],
                end_row     = SampleGroupIndex[k, "end_row"],
                indent      = FALSE,
                background  = "#FFE6CC",
                label_row_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #D79B00;"
            )
        }
        
        return(OutputTable)
    }
#
#---/ tmb result table /----------------------------------------------------------------------------
    WES_report.TMB.Result.Table <- function( TmbResult, SampleInfo, CancerCode, TonlyMode="Tonly" )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("reshape2"))
        # tmb sample info ------------------------------------------------------
        tmb_sinfo <- rbind(
            SampleInfo %>% filter( sample_info == CancerCode, sample_group != "" ) %>% arrange(sample_group, seq_id),
            SampleInfo %>% filter( sample_info == CancerCode, sample_group == "" ) %>% mutate( sample_group = "No Grouped" ) %>% arrange(seq_id)
        )
        # tmb data -------------------------------------------------------------
        if( TonlyMode == "TonlyGermlineFiltered" ){ TONLY_MODE = "TonlyGermlineFiltered" }else{ TONLY_MODE = "Tumor-only" }
        tmb_data <- TmbResult %>% filter(  variant_call_mode == "NT-paired", seq_id %in% tmb_sinfo$seq_id )
        
        if( nrow(tmb_data) > 0 )
        {
            tmb_data <- rbind(
                tmb_data,
                TmbResult %>% filter(  variant_call_mode %in% TONLY_MODE, seq_id %in% tmb_sinfo$seq_id ) %>% filter( seq_id %nin% tmb_data$seq_id )
            )
        }else{
            tmb_data <- TmbResult %>% filter(  variant_call_mode %in% TONLY_MODE, seq_id %in% tmb_sinfo$seq_id )
        }
        tmb_data <- tmb_data[match(tmb_sinfo$seq_id, tmb_data$seq_id), ]

        # sample groups 
        SampleGroupList <- unique(tmb_sinfo$sample_group)

        # sample group index ---------------------------------------------------
        SampleGroupIndex <- ldply(lapply(indexing(SampleGroupList), function(sg)
        {
            data.frame(
                sample_group = sg, 
                start_row    = head(which( tmb_data$seq_id %in% tmb_sinfo[which(tmb_sinfo$sample_group == sg), "seq_id"] ), 1),
                end_row      = tail(which( tmb_data$seq_id %in% tmb_sinfo[which(tmb_sinfo$sample_group == sg), "seq_id"] ), 1)
            )
        }))
        # sample name modification ---------------------------------------------
        tmb_data$sample_name <- tmb_sinfo[match(tmb_data$seq_id, tmb_sinfo$seq_id), "sample_name"]
        for( sg in SampleGroupList[ SampleGroupList != "No Grouped" ] )
        {
            sample_group_tag      <- paste0(gsub("\\-", "_", sg), "_")
            tmb_data$sample_name  <- gsub(sample_group_tag, "", tmb_data$sample_name)
        }
        # msi status column ----------------------------------------------------
        tmb_data$status <- sapply( tmb_data$TMB_value, function(z) ifelse( as.numeric(z) > 10 , "TMB High", "" ) )

        # output table ---------------------------------------------------------
        StatusBackgroundColors <- sapply( tmb_data$status, function(w) ifelse( w == "TMB High", "#B85450", "#FFFFFF") )
        StatusLabelColors      <- sapply( tmb_data$status, function(w) ifelse( w == "TMB High", "#FFFFFF", "#000000") )
        ValueLabelColors       <- sapply( tmb_data$status, function(w) ifelse( w == "TMB High", "#C73500", "#000000") )   

        OutputTable <- tmb_data %>% dplyr::select(c("sample_name","TMB_value","variants","status","variant_call_mode")) %>%
            dplyr::rename(
                'SAMPLE NAME'=sample_name, 'TMB (Mut/Mb)'=TMB_value, 'VARIANTS'=variants, 'TMB STATUS'=status, 'VARIANT CALL'=variant_call_mode
            ) %>%
            kbl( escape = FALSE, align=c('c'), row.names=F ) %>% 
            kable_styling(bootstrap_options='condensed', full_width = TRUE, font_size=10 ) %>%
            row_spec(0, font_size = 12, bold=TRUE, extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;" ) %>% 
            column_spec(c(1), width_min="12em", width="15em" ) %>%
            column_spec(c(2), width_min="6em",  width="8em", color=ValueLabelColors  ) %>%
            column_spec(c(3), width_min="6em",  width="8em"  ) %>%
            column_spec(c(4), width_min="5em",  width="8em", background=StatusBackgroundColors, color=StatusLabelColors  ) %>%
            column_spec(c(5), width_min="5em",  width="8em"  ) 

        for( k in 1:nrow(SampleGroupIndex) )
        {
            OutputTable <- OutputTable %>% pack_rows(
                group_label = SampleGroupIndex[k, "sample_group"],
                start_row   = SampleGroupIndex[k, "start_row"],
                end_row     = SampleGroupIndex[k, "end_row"],
                indent      = FALSE,
                background  = "#DAE8FC",
                label_row_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #6C8EBF;"
            )
        }
        
        return(OutputTable)
    }
#
#### APPENDIX ######################################################################################
#
#---/ APPENDIX : wes variant filters - GENE /-------------------------------------------------------
    WES_report_appendix.Variant.Filter.Genes <- function()
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        FilterGenes <- data.frame(rbind(
            c("GENE LOCUS GROUP", " ", "HGNC Protein-Coding Genes Only"            ),
            c( rep(" ", 3)),
            c("WHITELIST GENES" , " ", "14 Eseential Genes of Cancer NGS Panel Test")
        )) %>% 
            kbl( align=c('c','c','l'), row.names=FALSE, escape=FALSE, col.names=NULL) %>%
            kable_styling(bootstrap_options="condensed", full_width=FALSE, font_size = 10, position="left") %>%
            column_spec(1, width_min = "15em", width = "15em", background = "#6C8EBF" ) %>%
            column_spec(2, width_min = "1em", width = "1em" ) %>%
            column_spec(3, width_min = "48em", width = "48em" ) %>%
            row_spec(c(1,3), bold = TRUE, font_size=10, extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #DAE8FC" ) %>%
            row_spec(c(2), background = "#FFFFFF", font_size=3, extra_css = "border-bottom: 0px solid" ) %>%
            pack_rows(group_label="GENE FILTERS", start_row=1, end_row=3, indent=FALSE, 
                label_row_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #DAE8FC;" )
        return(FilterGenes)
    }
#
#---/ APPENDIX : wes variant essential genes /------------------------------------------------------
    WES_report_appendix.Essential.Gene.List <- function()
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        EssentialGenes <- c("ALK","BRAF","BRCA1","BRCA2","EGFR","ERBB2","IDH1","IDH2","KIT","KRAS","MYC","MYCN","NRAS","PDGFRA", "","" )
        EGTable <- data.frame(matrix(EssentialGenes, ncol=8, byrow=T)) %>% 
            kbl( align = 'c', row.names = FALSE, escape = FALSE, col.names=NULL ) %>% 
            kable_styling(bootstrap_options="condensed", full_width=FALSE, font_size=10, position="left") %>%
            column_spec(c(1:8), width_min = "8em",  width = "8em", border_right=TRUE, border_left=TRUE,
                extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #cccccc;") %>%
            pack_rows(group_label="SOLID TUMOR NGS PANEL ESSENTIAL GENES", start_row=1, end_row=2, indent=FALSE,
                label_row_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #cccccc;")
        return(EGTable)
    }
#
#---/ APPENDIX : wes variant filters - FREQUENCY /--------------------------------------------------
    WES_report_appendix.Variant.Filter.Frequency <- function()
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        FilterFreq <- data.frame(rbind(
            c("MINIMUM READ DEPTH (DP)"  , " ", "DP >= 10"                         ),
            c(rep(" ",3)),
            c("ALTERED ALLELE DEPTH (AD)", " ", "AD >= 2"                          ),
            c(rep(" ",3)),
            c("ALLELE FREQUENCY (AF)"    , " ", "AF >= 2% (SNV), AF >= 5% (INDEL)" ),
            c(rep(" ",3)),
            c("1000 Genomes FREQUENCY"   , " ", "Phase3 AF & EastAsian AF <= 1%"   ),
            c(rep(" ",3)),
            c("gnomAD FREQUENCY"         , " ", "AF & EastAsian AF <= 1%"          )
        )) %>%
            kbl( align = c('c','c','l'), row.names = FALSE, escape = FALSE, col.names=NULL ) %>%
            kable_styling(bootstrap_options="condensed", full_width=FALSE, font_size=10, position="left") %>%
            column_spec(c(1), width_min = "15em",  width = "15em", background = "#6C8EBF", color="#FFFFFF" ) %>%
            column_spec(2, width_min = "1em", width = "1em" ) %>%
            column_spec(3, width_min = "25em", width = "50em" ) %>%
            row_spec(c(1,3,5,7,9), extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #DAE8FC") %>% 
            row_spec(c(2,4,6,8),  background = "#FFFFFF", extra_css = "border-bottom: 0px solid;" ) %>% 
            pack_rows(group_label="READ DEPTH FILTERS", start_row=1, end_row=3, indent=FALSE,
                label_row_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #DAE8FC;" ) %>%
            pack_rows(group_label="ALLELE FREQUENCY FILTERS", start_row=5, end_row=5, indent=FALSE,
                label_row_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #DAE8FC;" ) %>%
            pack_rows(group_label="POPULATION FREQUENCY FILTERS", start_row=7, end_row=9, indent=FALSE,
                label_row_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #DAE8FC;" ) 
        return(FilterFreq)
    }
#
#---/ APPENDIX : wes analysis terminology : ENGLISH /-----------------------------------------------
    WES_report.appendix.Terms.Eng <- function()
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        TermTable <- data.frame(rbind(
            c("TOTAL READ",      " ",  "The total number of short sequences produced by next-generation sequencing (NGS) process."                      ),
            c("Q30 RATE",        " ",  "The quality score that indicates the probability of an incorrect base call. Q30 score means 99.9% accuracy."    ),
            c("GC RATIO",        " ",  "The percentage of Guanine(G) or Cytosine(C) bases in sequences. Extreme GC ratio can make sequencing difficult."),
            c("DUPLICATES RATE", " ",  "The percentage of mapped reads that are duplicates."                                                            ),
            c("ON-TARGET RATE",  " ",  "The percentage of sequencing reads that map to the target region."                                              ),
            c("TARGET DEPTH",    " ",  "Also known as Sequencing Depth. The number of times a specific genomic region is sequenced."                    ),
            c("TARGET COVERAGE", " ",  "The average number of reads that cover a target region."                                                        )
        )) %>%
            kbl( escape=FALSE, col.names=c("TERMS"," ","DEFINITION"), align=c('c','c','l') ) %>%
            kable_styling(bootstrap_options="condensed", full_width=TRUE, font_size=10 ) %>%
            row_spec(0, background = "#6C8EBF", color = "#FFFFFF", bold=TRUE, font_size=12 ) %>%
            column_spec(1, background = "#DAE8FC", width_min="15em", width="15em" ) %>%
            column_spec(2, width_min="1em", width="1em") %>%  
            column_spec(3, width_min="30em", width="50em")
        #-----------------------------------------------------------------------
        return(TermTable)
    }
#
#---/ APPENDIX : wes reference, database, swtools /-------------------------------------------------
    WES_report.appendix.ref.swtool.standard <- function()
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        resource_table <- data.frame(rbind(
            c("REFERENCE GENOME VERSION",   " ", "hg19 (GRCh37)"                                                                             ),
            c("REFERENCE GENOME FASTA"  ,   " ", "human_g1k_v37_decoy,fasta"                                                                 ),
            rep(" ", 3),
            c("KNOWN SNP DATABASE",         " ", "dbSNP v138, 1000G phase1"                                                                  ),
            c("KNOWN INDEL DATABASE",       " ", "hg19 known indels, 1000G phase1, Mills and 1000G goldstandard"                             ),
            c("GERMLINE VARIANTS DATABASE", " ", "gnomAD v2.1, KOVA hg19"                                                                    ),
            c("GENE ANNOTATION",            " ", "RefSeq (2020-10-26), Ensembl release 110.1, HGNC (2023-10)"                                ),
            c("FUNCTIONAL ANNOTATION",      " ", "dbSNP v154, ClinVar (2023-12-09), SIFT v2.2.2, POLYPHEN v5.2.2, COSMIC v99, Alpha-Missense"),
            rep(" ", 3),
            c("FASTQ QC",                   " ", "fastqc v0.12.1, fastp v0.23.4, Fastq-Screen v0.15.3, MultiQC v1.18"                        ),
            c("ALIGNMENT",                  " ", "bwa mem v0.7.17, picard v3.1.0, GATK3 v3.8, GATK4 v4.4.0.0"                                ),
            c("NGS QC",                     " ", "GATK4 v4.4.0.0, mosdepth v0.3.6, alfred v0.2.6, MultiQC v1.18"                             ),
            c("VARIANT CALL",               " ", "Mutect2 v4.4.0.0, GATK4 v4.4.0.0"                                                          ),
            c("ANNOTATION",                 " ", "ensembl-vep release 110.1, vcf2maf v1.6.21"                                                )
        )) %>% 
            kbl( escape=FALSE, align=c('l'), col.names=NULL ) %>%
            kable_styling(bootstrap_options="condensed", full_width=TRUE, font_size=10 ) %>%
            column_spec(1, background = "#DAE8FC", width_min="15em", width="15em" ) %>% 
            column_spec(2, width_min="1em", width="1em" ) %>%  
            column_spec(3, width_min="30em", width="50em" ) %>%            
            row_spec(c(3,9), background="#FFFFFF", extra_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #6C8EBF;" ) %>%
            row_spec(c(1:2,4:8,10:14), extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;" ) %>%
            pack_rows(group_label="REFERENCES", start_row=1,  end_row=2, 
                label_row_css="border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;"  ) %>%
            pack_rows(group_label="DATABASES",  start_row=4,  end_row=8,
                label_row_css="border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;"  ) %>%
            pack_rows(group_label="ANALYSIS TOOLS", start_row=10, end_row=14,
                label_row_css="border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;"  )
        return(resource_table)
    }
    WES_report.appendix.ref.swtools.advanced <- function()
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        resource_table <- data.frame(rbind(
            c("REFERENCE GENOME VERSION",   " ", "hg19 (GRCh37)"                                                                             ),
            c("REFERENCE GENOME FASTA"  ,   " ", "human_g1k_v37_decoy,fasta"                                                                 ),
            rep(" ", 3),
            c("KNOWN SNP DATABASE",         " ", "dbSNP v138, 1000G phase1"                                                                  ),
            c("KNOWN INDEL DATABASE",       " ", "hg19 known indels, 1000G phase1, Mills and 1000G goldstandard"                             ),
            c("GERMLINE VARIANTS DATABASE", " ", "gnomAD v2.1, KOVA hg19"                                                                    ),
            c("GENE ANNOTATION",            " ", "RefSeq (2020-10-26), Ensembl release 110.1, HGNC (2023-10)"                                ),
            c("FUNCTIONAL ANNOTATION",      " ", "dbSNP v154, ClinVar (2023-12-09), SIFT v2.2.2, POLYPHEN v5.2.2, COSMIC v99, Alpha-Missense"),
            rep(" ", 3),
            c("INDIVIDUAL MATCHING",        " ", "NGSCheckMate v1.0.1"                                                                       ),
            c("HLA TYPING",                 " ", "OptiType v1.3.3, HLA-LA v1.0.3"                                                            ),
            c("MICROSATELLITE INSTABILITY", " ", "MANTIS v1.0.5, MSIsensor2 v0.1"                                                            ),
            c("ANALYSIS AND PLOTS",         " ", "R v4.3.1"                                                                                  )
        )) %>% 
            kbl( escape=FALSE, align=c('l'), col.names=NULL ) %>%
            kable_styling(bootstrap_options="condensed", full_width=TRUE, font_size=10 ) %>%
            column_spec(1, background = "#DAE8FC", width_min="10em", width="10em") %>% 
            column_spec(2, width_min="1em", width="1em" ) %>%  
            column_spec(3, width_min="30em", width="50em", bold=TRUE  ) %>%            
            row_spec(c(3,9), background="#FFFFFF", extra_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #6C8EBF;" ) %>%
            row_spec(c(1:2,4:8,10:13), extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;" ) %>%
            pack_rows(group_label="REFERENCES", start_row=1,  end_row=2, 
                label_row_css="border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;"  ) %>%
            pack_rows(group_label="DATABASES",  start_row=4,  end_row=8,
                label_row_css="border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;"  ) %>%
            pack_rows(group_label="ANALYSIS TOOLS", start_row=10, end_row=13,
                label_row_css="border-top: 0px solid; border-bottom: 1px solid; border-color: #6C8EBF;"  )
        return(resource_table)
    }
#
#---/ APPENDIX : organoid science preset genes list /-----------------------------------------------
    WES_report.appendix.client.preset.genes <- function( ClientID=NULL, cnf=cnf )
    {
        suppressPackageStartupMessages(library("plyr"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        suppressPackageStartupMessages(library("RMySQL"))
        #-----------------------------------------------------------------------
        if( is.null(ClientID) ){ stop("no client id. required.")}
        #-----------------------------------------------------------------------
        dbCon  <- dbConnect(dbDriver("MySQL"), host=cnf$db$host, user=cnf$db$user, port=cnf$db$port, password=cnf$db$pw, db=cnf$db$db_data )
        PresetFromDB <- dbGetQuery(dbCon, sprintf("SELECT * FROM preset_genes WHERE client_id = '%s'", ClientID))
        dbDisconnect(dbCon)
        #-----------------------------------------------------------------------
        pstb <- rbind(
            PresetFromDB %>% filter( cancer_code != "" ) %>% arrange( cancer_code ),
            PresetFromDB %>% filter( cancer_code == "" ) 
        ) 
        pstb$genes <- gsub(";", ", ", pstb$genes )
        
        preset_table <- pstb %>% dplyr::select(c("preset_name","cancer_code","genes")) %>% 
            kbl( escape=FALSE, align=c('c','c','l'), col.names=c("PRESET NAME","CANCER","GENES") ) %>% 
            kable_styling(bootstrap_options="condensed", full_width=TRUE, font_size=10 ) %>%
            row_spec(0, background="#6C8EBF", color="#FFFFFF", font_size=11) %>%
            row_spec(c(1:nrow(pstb)), extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #6C8EBF;" ) %>%
            column_spec(1, background="#DAE8FC", width_min="10em", width="12em") %>%
            column_spec(2, width_min="2em", width="2em", border_right=TRUE,
                extra_css = "border-right: 1px solid; border-right-color: #DAE8FC;" ) %>%
            column_spec(3, width_min="35em")             
        return(preset_table)
    }
#
#---/ APPENDIX : concordance terms /----------------------------------------------------------------
    WES_report.Concordance.Terms <- function()
    {
        ConcordTermTable <- as.data.frame(rbind(
            c("CONCORDANCE", " ", "Normalized Tanimoto-coefficient (Tm)"),
            c(" ",           " ", ": general tanimoto-coefficient divided with possible maximum tanimoto-coefficient"),
            c(" ",           " ", "  ( between 0 ~ 1 )"),
            rep(" ",3),
            rep(" ",3),
            rep(" ",3),
            rep(" ",3),
            c("MATCH RATE" , " ", "Ratio of common variants with variants of sample-A or sample-B. (percentage) ")
        )) %>% 
            kbl( align = c('c','c','l'), row.names = FALSE, escape = FALSE, col.names=NULL ) %>%
            kable_styling(bootstrap_options="condensed", full_width=FALSE, font_size=11, position="left") %>%
            column_spec(1, width_min = "12em", width = "12em", background = "#6C8EBF", color="#FFFFFF" ) %>%
            column_spec(2, width_min = "1em",  width = "1em"  ) %>%
            column_spec(3, width_min = "50em", width = "50em" ) %>%
            row_spec(1,   extra_css = "border-top: 1px solid; border-bottom: 0px solid; border-color: #DAE8FC;") %>% 
            row_spec(2,   background="#FFFFFF", font_size=10, extra_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #DAE8FC;") %>% 
            row_spec(3,   background="#FFFFFF", font_size=10, extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #DAE8FC;") %>%
            row_spec(4,   background="#FFFFFF", font_size=15, extra_css = "border-top: 1px solid; border-bottom: 0px solid; border-color: #DAE8FC;") %>%
            row_spec(5:7, background="#FFFFFF", font_size=15, extra_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #DAE8FC;") %>%
            row_spec(8,   extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #DAE8FC;") 
        return(ConcordTermTable)
    }
#
#---/ vaf correlation plot legnd /------------------------------------------------------------------
    WES_report.Vaf.Corr.Plot.Legnd <- function()
    {
        ColumnBackgroundColor <- c("#FFFFFF","#FFFFFF","#B85450","#FFFFFF","#2D7600","#FFFFFF","#006EAF","#FFFFFF","#adadad")
        ColumnLabelColor      <- c("#000000","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")
        VafPlotLegend <- as.data.frame(rbind(
            c("X,Y-axis",  " ", "Variant Allele Frequency (VAF)"                  ),
            rep(" ",3),
            c("RED dot",   " ", "Common Somatic Variants"                         ),
            rep(" ",3),
            c("GREEN dot", " ", "Sample-A Somatic Variants (generally TS sample)" ),
            rep(" ",3),
            c("BLUE dot ", " ", "Sample-B SOmatic Variants (generally ORG sample)"),
            rep(" ",3),
            c("GRAY dot ", " ", "Unfiltered Raw Variants"                         )
        )) %>% 
            kbl( align = c('c','c','l'), row.names = FALSE, escape = FALSE, col.names=NULL ) %>%
            kable_styling(bootstrap_options="condensed", full_width=FALSE, font_size=10, position="left") %>%
            column_spec(1, width_min = "10em", width = "10em", background = ColumnBackgroundColor, color=ColumnLabelColor ) %>%
            column_spec(2, width_min = "1em",  width = "1em"  ) %>%
            column_spec(3, width_min = "50em", width = "50em" ) %>%
            row_spec(c(2,4,6,8),   extra_css = "border-top: 0px solid; border-bottom: 0px solid; border-color: #e0e0e0;" ) %>%
            row_spec(c(1,3,5,7,9), extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-color: #e0e0e0;")
            
        return(VafPlotLegend)
    }
#
#---/ msi-h cut off /-------------------------------------------------------------------------------
    WES_report.MSI.Cutoff <- function()
    {
        msi_cutoff <- as.data.frame(rbind(
            c("Tumor Only Variant Call",     " ", "MSIsensor2 score > 20"),
            c("Matched Normal Variant Call", " ", "Both MSIsensor2 score > 20 and MANTIS score > 0.4")
        )) %>% 
            kbl( align = c('c','c','l'), row.names = FALSE, escape = FALSE, col.names=NULL ) %>%
            kable_styling(bootstrap_options="condensed", full_width=FALSE, font_size=10, position="left") %>%
            column_spec(1, width_min = "15em",  width = "20em", background="#FFE6CC" ) %>%
            column_spec(2, width_min = "1em",  width = "1em"  ) %>%
            column_spec(3, width_min = "20em", width = "25em" ) %>%
            row_spec(c(1:2), extra_css = "border-top: 1px solid; border-bottom: 1px solid; border-color: #adadad;" ) 
        return(msi_cutoff)
    }

#
#### DATA EXPORT ###################################################################################
#
#---/ FASTQ-FILES-INFO /----------------------------------------------------------------------------
    WES_report.Export.File.List.FASTQ <- function()
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        fastqFiles <- as.data.frame(rbind(
            c(" ", "[ Sample ID ]_R1.fastq.gz", " ", "Raw FASTQ Read 1 File"      ),
            c(" ", "[ Sample ID ]_R2.fastq.gz", " ", "Raw FASTQ Read 2 File"      ),
            c(" ", "FASTQ.md5sum"             , " ", "MD5 checksum of FASTQ Fles" )
        )) %>% 
            kbl( align=c("l"), row.names=F, col.names=c(" ","FILENAME"," ","DESCRIPTION"), escape=FALSE ) %>%
            kable_classic_2('striped', full_width=F, font_size = 11, position='left' ) %>%
            row_spec(0, background = "#D5E8D4", bold=TRUE) %>%
            column_spec(1, width_min="1em",  width="1em") %>%
            column_spec(2, width_min="36em", width="36em") %>%
            column_spec(3, width_min="1em",  width="1em") %>%
            column_spec(4, width_min="25em", width="25em") 
        return(fastqFiles)
    }
#
#---/ BAM-FILES-INFO |------------------------------------------------------------------------------
    WES_report.Export.File.List.BAM <- function()
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        bamFiles <- as.data.frame(rbind(
            c(" ", "[ Sample ID ].bam"     , " ", "Aligned BAM File"         ),
            c(" ", "[ Sample ID ].bam.bai" , " ", "BAM Index File"           ),
            c(" ", "BAM.md5sum"            , " ", "MD5 checksum of BAM Files")
        )) %>% 
            kbl( align=c("l"), row.names=F, col.names=c(" ","FILENAME"," ","DESCRIPTION"), escape=FALSE ) %>%
            kable_classic_2('striped', full_width=F, font_size = 11, position='left' ) %>%
            row_spec(0, background = "#D5E8D4", bold=TRUE) %>%
            column_spec(1, width_min="1em",  width="1em") %>%
            column_spec(2, width_min="36em", width="36em") %>%
            column_spec(3, width_min="1em",  width="1em") %>%
            column_spec(4, width_min="25em", width="25em")
        return(bamFiles)
    }
#
#---/ VCF-FILES-INFO /------------------------------------------------------------------------------
    WES_report.Export.File.List.VCF <- function( TONLY=TRUE, TS_N=FALSE, ORG_N=FALSE, TONLY_GF=FALSE, GERMLINE=FALSE  )
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        vcf_types = c()
        if( TONLY    ){ vcf_types = c(vcf_types, 1) }
        if( TS_N     ){ vcf_types = c(vcf_types, 2) }
        if( ORG_N    ){ vcf_types = c(vcf_types, 3) }
        if( TONLY_GF ){ vcf_types = c(vcf_types, 4) }
        if( GERMLINE ){ vcf_types = c(vcf_types, 5) }
        #----------------------------------------------------------------------#
        vcfFiles <- as.data.frame(rbind(
            c(" ", "[ Sample ID ].tumor.only.vcf"                   , " ", "Annotated VCF file. Variant Calling : Tumor-only"                       ),
            c(" ", "[ Sample ID ].matched.normal.TISSUE.vcf"        , " ", "Annotated VCF file. Variant Calling : Matched Normal - TISSUE"          ),
            c(" ", "[ Sample ID ].matched.normal.ORGANOID.vcf"      , " ", "Annotated VCF file. Variant Calling : Matched Normal - ORGANOID"        ),
            c(" ", "[ Sample ID ].tumor.only.germline.filtered.vcf" , " ", "Annotated VCF file. Variant Calling : Tumor-only and Germline Filtered" ),
            c(" ", "[ Sample ID ].germline.vcf"                     , " ", "VCF file. Variant Calling : DeepVariant Germline Call"                  ),
            c(" ", "VCF.md5sum"                                     , " ", "MD5 checksum of VCF Fles"                                               )
        )) %>% slice(c(vcf_types,6)) %>%
            kbl( align=c("l"), row.names=F, col.names=c(" ","FILENAME"," ","DESCRIPTION"), escape=FALSE ) %>%
            kable_classic_2('striped', full_width=F, font_size = 11, position='left' ) %>%
            row_spec(0, background = "#D5E8D4", bold=TRUE) %>%
            column_spec(1, width_min="1em",  width="1em") %>%
            column_spec(2, width_min="36em", width="36em") %>%
            column_spec(3, width_min="1em",  width="1em") %>%
            column_spec(4, width_min="25em", width="25em") 
        return(vcfFiles)
    }
#
#---/ MAF-FILES-INFO /------------------------------------------------------------------------------
    WES_report.Export.File.List.MAF <- function( TONLY=TRUE, TS_N=FALSE, ORG_N=FALSE, TONLY_GF=FALSE, GERMLINE=FALSE )
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        maf_types = c()
        if( TONLY    ){ maf_types = c(maf_types, 1, 2) }
        if( TS_N     ){ maf_types = c(maf_types, 3, 4) }
        if( ORG_N    ){ maf_types = c(maf_types, 5, 6) }
        if( TONLY_GF ){ maf_types = c(maf_types, 7) }
        if( GERMLINE ){ maf_types = c(maf_types, 8) }
        #----------------------------------------------------------------------#
        mafFiles <- as.data.frame(rbind(
            c(" ", "[ Sample ID ].tumor.only.maf"                                   , " ", "MAF Converted from Tumor-only VCF File"                      ),
            c(" ", "[ Sample ID ].tumor.only.somatic.variants.maf"                  , " ", "Somatic Variants of Tumor-only MAF"                          ),
            c(" ", "[ Sample ID ].matched.normal.TISSUE.maf"                        , " ", "MAF Converted from Matched Normal TISSUE VCF"                ),
            c(" ", "[ Sample ID ].matched.normal.TISSUE.somatic.variants.maf"       , " ", "Somatic Variants of Matched Normal TISSUE MAF"               ),
            c(" ", "[ Sample ID ].matched.normal.ORGANOID.maf"                      , " ", "MAF Converted From Matched Normal ORGANOID VCF"              ),
            c(" ", "[ Sample ID ].matched.normal.ORGANOID.somatic.variants.maf"     , " ", "Somatic Variants of Matched Normal ORGANOID MAF"             ),
            c(" ", "[ Sample ID ].tumor.only.germline.filtered.somatic.variants.maf", " ", "Tumor-only called and Germline Filtered Somatic Variants MAF"),
            c(" ", "[ Sample ID ].germline.maf"                                     , " ", "MAF Converted From DeepVariant Germline Call VCF"            ),
            c(" ", "MAF.md5sum"                                                     , " ", "MD5 checksum of MAF Files"                                   )
        )) %>% slice(c(maf_types,9)) %>%
            kbl( align=c("l"), row.names=F, col.names=c(" ","FILENAME"," ","DESCRIPTION"), escape=FALSE ) %>%
            kable_classic_2('striped', full_width=F, font_size = 11, position='left' ) %>%
            row_spec(0, background = "#D5E8D4", bold=TRUE) %>%
            column_spec(1, width_min="1em",  width="1em") %>%
            column_spec(2, width_min="36em", width="36em") %>%
            column_spec(3, width_min="1em",  width="1em") %>%
            column_spec(4, width_min="25em", width="25em")  
        return(mafFiles)
    }
#
#---/ REPORT-FILES-INFO /---------------------------------------------------------------------------
    WES_report.Export.File.List.REPORT <- function( STANDARD=TRUE, ADVANCED=TRUE, ONCOPLOT_PNG=TRUE, MAF_XLSX=TRUE, NONHUMAN=FALSE )
    {
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("Hmisc"))
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        #-----------------------------------------------------------------------
        report_contents = c()
        if( STANDARD     ){ report_contents = c(report_contents, 1, 2) }
        if( ADVANCED     ){ report_contents = c(report_contents, 3, 4) }
        if( ONCOPLOT_PNG ){ report_contents = c(report_contents, 5   ) }
        if( ONCOPLOT_PNG ){ report_contents = c(report_contents, 6   ) }
        if( NONHUMAN     ){ report_contents = c(report_contents, 7   ) }
        report_contents = c(report_contents, 8)
        #----------------------------------------------------------------------#
        reportFiles <- as.data.frame(rbind(
            c( " ", "[ Order ID ].Standard_Analysis_Report.pdf"                                      , " ", "NGS QC and Statistics Report. PDF format" ),
            c( " ", "[ Order ID ].Standard_Analysis_Report.html"                                     , " ", "NGS QC and Statistics Report. HTML format"),
            c( " ", "[ Order ID ].Advanced_Analysis_Report.pdf"                                      , " ", "WES Analysis Resport. PDF format"         ),
            c( " ", "[ Order ID ].Advanced_Analysis_Report.html"                                     , " ", "WES Analysis Resport. HTML format"        ),
            c( " ", "[ Order ID ].[ Tissue ].[ Geneset ].[ Variant Call ].cumulative.samples.oncoplots.png" , " ", "Cumlative Samples Oncoplots"              ),
            c( " ", "[ Order ID ].Preset.Genes.Variants.[ Variant Call ].xlsx"                       , " ", "Somatic Variants List. MS Excel format"   ),
            c( " ", "[ Order ID ].NonHumanSamples.pdf"                                               , " ", "Non-Human Sample Check Result. PDF format"),
            c( " ", "[ Order ID ].Data_List.pdf"                                                     , " ", "Data List Description. PDF format"        )
        )) %>% slice(report_contents) %>%
            kbl( align=c("l"), row.names=F, col.names=c(" ","FILENAME"," ","DESCRIPTION"), escape=FALSE ) %>%
            kable_classic_2('striped', full_width=F, font_size = 11, position='left' ) %>%
            row_spec(0, background = "#D5E8D4", bold=TRUE) %>%
            column_spec(1, width_min="1em",  width="1em") %>%
            column_spec(2, width_min="36em", width="36em") %>%
            column_spec(3, width_min="1em",  width="1em") %>%
            column_spec(4, width_min="25em", width="25em")  
        return(reportFiles)
    }
#



