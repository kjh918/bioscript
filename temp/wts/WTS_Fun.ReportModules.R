
# "cnf" REQUIRED. 
# cnf = read_yaml( config_file )

#---| INFO TABLES |----------------------------------------------------------------------------------------------------------------------------------#
#' @description WTS sample information table
#' @param SampleInfoTable  sample information table, if NULL, will be loaded from database.
#' @param SeqFolderId      REQUIRED. analysis batch id (seq_folder)
#' @param RmSeqId          sample id to eliminate from report
#' @param FontSize         font size. default = 10
#' @export 
    WTS_report.Sample.Info.Table <- function( SampleInfoTable=NULL, SeqFolderId=NULL, RmSeqId=NULL, FontSize=10,
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_info
    )
    {
        require(dplyr)
        require(Hmisc)
        require(knitr)
        require(kableExtra)
        #-----------------------------------------------------------------------
        if( is.null(SampleInfoTable) )
        {
            if( is.null(SeqFolderId) ){ stop("No SeqFolderId for loading sample info from database. STOPPED.") }else{
                require(RMySQL)
                dbCon        <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
                sample_info  <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SeqFolderId))
                dbDisconnect(dbCon)
            }
        }else{
            sample_info <- SampleInfoTable
        }
        # seq-id filter --------------------------------------------------------
        if( !is.null(RmSeqId) )
        {
            rm_seq_ids <- unlist(strsplit(RmSeqId, ","))
            sample_info <- sample_info %>% filter( seq_id %nin% rm_seq_ids )
        }
        #-----------------------------------------------------------------------
        sample_groups <- unique(sample_info$sample_group)
        #======================================================================#
        sinfo_table <- sample_info[,c("sample_group", "sample_name","sample_label","sample_tissue","sample_info")] %>% dplyr::rename(
            'Sample Group ID'= sample_group, 
            'Sample ID'      = sample_name,
            'Sample Name'    = sample_label,
            'Sample Tissue'  = sample_tissue,
            'Sample Info'    = sample_info
        ) %>% kbl( escape = FALSE, align = c('l','l','l','c','c') ) %>% 
            kable_styling( bootstrap_options = "condensed", full_width = TRUE, font_size = FontSize ) 
        #-----------------------------------------------------------------------
        for( SG in sample_groups )
        { sinfo_table <- sinfo_table %>% row_spec( max(which(sample_info$sample_group == SG)), extra_css = "border-bottom: 1px dotted; border-bottom-color: #4c4c4c" ) }
        #-----------------------------------------------------------------------
        sinfo_table <- sinfo_table %>% 
            row_spec( 0, font_size = FontSize+1, bold = TRUE, background = "#82B366", color = "#FFFFFF", extra_css = "border-bottom: 1px solid; border-bottom-color: #4c4c4c" ) %>%
            column_spec( 1, extra_css = "border-bottom: 1px dotted; border-bottom-color: #4c4c4c" ) %>%
            collapse_rows(1, valign = "middle")
        #======================================================================#
        return(sinfo_table)
    }
#--------------------------------------------------------------------------------------------------#
#' @description create RNAseq QC data summary table
#' @param SeqFolderId  REQUIRED. analysis batch id (seq folder).
#' @param SampleInfo   sample info table. if NULL, will be loaded from database.  
#' @param FontSize     font size. default = 10
#' @export
    WTS_report.NGS.Stat.Table <- function( SeqFolderId=NULL, SampleInfo=NULL, FontSize=10,
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_info
    )
    {
        require(dplyr)
        require(knitr)
        require(kableExtra)
        # SeqFolderId check ----------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! No seq_folder id. REQUIRED. Please check again. STOPPED.") }
        # header color ---------------------------------------------------------
        #if( is.null(HeaderColor) ){ header_color <- "#cacfab" }else{ header_color <- HeaderColor }
        # qc data --------------------------------------------------------------
        QCData <- NGS_report.Load.Stats( SeqFolderID = SeqFolderId, Src="report" )
        # sample info check ----------------------------------------------------
        if( is.null(SampleInfo) )
        {
            require(RMySQL)
            dbCon        <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
            sample_info  <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SeqFolderId))
            dbDisconnect(dbCon)
        }else{
            sample_info <- SampleInfo
        }
        #-----------------------------------------------------------------------
        qc_data <- QCData %>% mutate( 
            seqc_total_rds        = formatC(round(seqc_total_rds), format = "f", big.mark = ",", drop0trailing = TRUE),
            star_total_rds        = formatC(round(star_total_rds), format = "f", big.mark = ",", drop0trailing = TRUE)
        #) %>% mutate(
        #    map_rate = paste0( star_pct_uniq_mapped, " / ", star_pct_multi_mapped )
        ) %>% dplyr::select(c(
            "seq_id","seqc_total_rds","star_total_rds","star_pct_uniq_mapped","star_pct_multi_mapped","pct_pf_q30_rate","pct_gc_content",
            "median_5to3bias","pct_coding_base","pct_utr_base","pct_intronic_base","pct_intergenic_base"
        )) 
        qc_data$seq_id <- sample_info[match(qc_data$seq_id, sample_info$seq_id), "sample_name" ]
        #-----------------------------------------------------------------------
        qc_stats <- qc_data %>% dplyr::rename(
            'SAMPLE ID'       = seq_id,
            'TOTAL'           = seqc_total_rds,
            'MAPPED'         = star_total_rds,
            'UNIQUE'          = star_pct_uniq_mapped,
            'MULTI'           = star_pct_multi_mapped,
            'Q30 RATIO'       = pct_pf_q30_rate,
            'GC RATIO'        = pct_gc_content,
            '5-3 BIAS'        = median_5to3bias,
            'CODING'          = pct_coding_base,
            'UTR'             = pct_utr_base,
            'INTRON'          = pct_intronic_base,
            'INTERGENIC'      = pct_intergenic_base
        )
        #======================================================================#
        qc_stat_table <- qc_stats %>% kbl( escape=FALSE, align='c', valign='middle' ) %>% 
            kable_styling( bootstrap_options = "condensed", full_width = TRUE, font_size = FontSize ) %>%
            column_spec(c(1),     width_min="5em",  width="8em") %>%
            column_spec(c(2:3),   width_min="3em",  width="3em") %>%
            column_spec(c(4),     width_min="4em",  width="6em") %>%
            column_spec(c(5:8),   width_min="4em",  width="5em") %>%
            column_spec(c(9:12),  width_min="4em",  width="4em") %>%
            row_spec( 0, font_size = FontSize, bold = TRUE, extra_css = "border-bottom: 1px solid; border-bottom-color: #4c4c4c" ) %>% 
            add_header_above(header = c(" "=1, "READS"=2, "ALIGN RATIO"=2,  " "=3, "BASE RATIO ON GENE STRUCTURE"=4), 
                font_size = FontSize-1, escape = FALSE,  line_sep = 8, bold = TRUE
            ) 
        #-----------------------------------------------------------------------
        SampleGroups <- unique(sample_info$sample_group)
        for( SG in SampleGroups ){
            groupSampleRows <- which(qc_stats[,1] %in% sample_info[which(sample_info$sample_group == SG), "sample_name" ] )
            qc_stat_table   <- qc_stat_table %>% pack_rows(
                group_label=SG, 
                start_row = min(groupSampleRows), 
                end_row =max(groupSampleRows),
                escape = FALSE,
                indent = FALSE,
                label_row_css = "border-bottom: 1px solid; border-bottom-color: #D5E8D4; background: #eaeaea"
            ) 
        }
        #======================================================================#
        return(qc_stat_table)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| MAPPING STATS BARPLOT |------------------------------------------------------------------------------------------------------------------------#
#' @description draw STAR alignment statistics barplot
#' @param SeqFolderId  REQUIRED. analysis batch id.
#' @param SampleInfo   sample info table. if NULL, seq-id will be used as sample name.
#' @note  data downloaded from database using 'NGS_report.Load.Stat' function in "NGS_Fun.ReportModules_v2.R"
#' @export 
    WTS_plot.Align.Stats.Barplot <- function( SeqFolderId=NULL, SampleInfo=NULL )
    {
        require(dplyr)
        require(ggplot2)
        require(scales)
        require(Cairo)
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        # SeqFolderId check ----------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! No seq_folder id. REQUIRED. Please check again. STOPPED.") }
        # colors ---------------------------------------------------------------
        AlignGroupColors <- c("Unique_Mapped"="#6C8EBF", "Multi_Mapped"="#82B366", "Unmapped"="#D79B00")
        #-----------------------------------------------------------------------
        WTS_Stats <- NGS_report.Load.Stats( SeqFolderID = SeqFolderId, Src="bam" )
        #-----------------------------------------------------------------------
        if( !is.null(SampleInfo) )
        { WTS_Stats$seq_id <- SampleInfo[match(WTS_Stats$seq_id, SampleInfo$seq_id), "sample_name"] }
        #-----------------------------------------------------------------------
        star_align_stats <- rbind(
            data.frame( Group="Unique_Mapped", WTS_Stats[, c("seq_id","star_uniq_mapped_rds") ] %>% dplyr::rename( rds=2 ) ),
            data.frame( Group="Multi_Mapped",  WTS_Stats[, c("seq_id","star_multi_mapped_rds")] %>% dplyr::rename( rds=2 ) ),
            data.frame( Group="Unmapped",      WTS_Stats[, c("seq_id","star_unmapped")        ] %>% dplyr::rename( rds=2 ) )
        )
        #-----------------------------------------------------------------------
        star_align_stats$Group <- factor(star_align_stats$Group, levels=c("Unmapped","Multi_Mapped","Unique_Mapped"))
        #-----------------------------------------------------------------------
        star_align_stats$seq_id <- factor(unique(star_align_stats$seq_id), levels=unique(star_align_stats$seq_id)[length(unique(star_align_stats$seq_id)):1])
        # plot range -----------------------------------ooo---------------------
        MaxLimitRDS <- star_align_stats %>% group_by( seq_id ) %>% reframe( total_rds = sum(rds) ) %>% .$total_rds %>% max() * 1.2
        #======================================================================#
        align_stat_barplot <- ggplot(star_align_stats, aes( x=rds, y=seq_id, fill=Group) ) + 
            geom_bar(stat="identity", position="stack", width=0.8) +
            coord_cartesian( xlim = c(0, MaxLimitRDS) ) +
            scale_fill_manual( values = AlignGroupColors ) +
            scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6), expand = c(0,0)) +
            labs( x = "Reads", y = "", fill = "") +
            theme(
                panel.background = element_rect(fill = '#FFFFFF'),
                panel.grid.major = element_line(color = '#e0e0e0', linetype = 'dotted'),
                axis.line.x  = element_line(colour = "#6f6d6d"),
                axis.ticks.y = element_blank(),
                axis.text.x  = element_text(family="Helvetica", face="bold", size=10, hjust=0.5, vjust=0.5, margin = margin(t = 0)),
                axis.text.y  = element_text(family="Helvetica", face="bold", size=10, hjust=1,   vjust=0.5, margin = margin(r = 0)),
                axis.title   = element_text(family="Helvetica", face="bold", size=12, hjust=0.5, vjust=0.5, margin = margin(t = 10, r = 10)),
                legend.text  = element_text(family="Helvetica", size=8),
                plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
            )
        #======================================================================#
        return(align_stat_barplot)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| MAPPING STATS BARPLOT or PIE-CHART |-----------------------------------------------------------------------------------------------------------#
#' @description draw STAR alignment statistics barplot
#' @param SeqFolderId   REQUIRED. analysis batch id.
#' @param SampleInfo    sample info table. if NULL, seq-id will be used as sample name.
#' @param PlotType      plot type. barplot = 'bar', piechart = 'pie'. default = 'bar'
#' @param SeqID         draw plot for specific sample only when PlotType = 'pie'. if NULL, draw all samples.
#' @param PlotGridCols  column numbers of plot_grid when PlotType = 'pie'. 
#' @export 
    WTS_plot.Align.Region.Ratio <- function( SeqFolderId=NULL, SampleInfo=NULL, PlotType="bar", 
        SeqID=NULL, PlotGridCols=2 
    )
    {
        require(ggplot2)
        require(scales)
        require(cowplot)
        require(Cairo)
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        # data check -----------------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! No seq folder id found. REQUIRD. please check again. STOPPED.")}
        # colors ---------------------------------------------------------------
        AlignGroupColors <- c("Coding"="#82B366", "UTR"="#6C8EBF", "Intron"="#9673A6", "Intergenic"="#D79B00")
        #-----------------------------------------------------------------------
        WTS_Stats <- NGS_report.Load.Stats( SeqFolderID = SeqFolderId, Src="bam" )
        #-----------------------------------------------------------------------
        if( !is.null(SampleInfo) )
        {
            WTS_Stats$seq_id <- SampleInfo[match(WTS_Stats$seq_id, SampleInfo$seq_id), "sample_name"]
        }
        #-----------------------------------------------------------------------
        region_align_stats <- rbind(
            data.frame( Group="Coding",     WTS_Stats[, c("seq_id","pct_coding_base")     ] %>% dplyr::rename( pct=2 ) ),
            data.frame( Group="UTR",        WTS_Stats[, c("seq_id","pct_utr_base")        ] %>% dplyr::rename( pct=2 ) ),
            data.frame( Group="Intron",     WTS_Stats[, c("seq_id","pct_intronic_base")   ] %>% dplyr::rename( pct=2 ) ),
            data.frame( Group="Intergenic", WTS_Stats[, c("seq_id","pct_intergenic_base") ] %>% dplyr::rename( pct=2 ) )
        )
        #-----------------------------------------------------------------------
        region_align_stats$Group <- factor(region_align_stats$Group, levels=c("Intergenic","Intron","UTR","Coding"))
        #-----------------------------------------------------------------------
        region_align_stats$seq_id <- factor(unique(region_align_stats$seq_id), levels=unique(region_align_stats$seq_id)[length(unique(region_align_stats$seq_id)):1])
        #-----------------------------------------------------------------------
        if( PlotType == "pie" )
        {
            # PIE CHART =======================================================#
            if( !is.null(SeqID) )
            {
                Input.plot <- region_align_stats %>% filter( seq_id == SeqID )
                plot_index <- indexing(as.character(Input.plot$seq_id))
                empty_plot <- ggplot()+theme_void() 
            }else{
                plot_index <- indexing(as.character(region_align_stats$seq_id))
            }
            #-------------------------------------------------------------------
            region_plot_list <- lapply( plot_index, function(idx)
            {
                stat_R <- region_align_stats %>% filter( seq_id == idx ) %>% mutate(label = paste0(Group, "\n", round(pct,1), "%"))

                gp <- ggplot(stat_R, aes( x=pct, y=seq_id, fill = Group) ) +  
                    geom_bar( stat="identity") +
                    coord_radial("x", expand= FALSE) +
                    scale_fill_manual( values = AlignGroupColors ) +
                    geom_text( aes( x = pct/2 + c(0, cumsum(pct)[-length(pct)]), label = label), fontface="bold", size=4, color="#FFFFFF") +
                    ggtitle( unique(stat_R$seq_id ) ) +
                    theme(
                        panel.background = element_rect(fill = '#FFFFFF'),
                        axis.title       = element_blank(),
                        axis.ticks.y     = element_blank(),
                        axis.text.y      = element_blank(),
                        #axis.ticks.length.x = unit( 2, "mm"),
                        axis.text.x      = element_text(family="Helvetica", face="bold", size=12, hjust=0.5, vjust=0.5 ),
                        plot.title       = element_text(face="bold", size=14, hjust=0.5),
                        legend.position  = "none",
                        plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
                    )
                return(gp)
            })
            if( length(region_plot_list) == 1 ){ region_plot_list = c(region_plot_list, list( empty_plot )) }
            #-------------------------------------------------------------------
            RegionRatioPlot <- plot_grid( plotlist=region_plot_list, ncol=PlotGridCols )    
            #==================================================================#
        }else{
            # BARPLOT =========================================================#
            gene_region_only_stats <- region_align_stats %>% filter( Group != "Intergenic" ) %>% group_by( seq_id ) %>% reframe( gene_pct = sum(pct) ) %>% as.data.frame()
            #-------------------------------------------------------------------
            RegionRatioPlot <- ggplot(region_align_stats, aes( x=pct, y=seq_id) ) + 
                geom_bar( aes(fill = Group), stat="identity", position="stack", width=0.8) +
                coord_cartesian( xlim = c(0, 100) ) +
                geom_bar( data = gene_region_only_stats, aes( x=gene_pct, y=seq_id ),  fill=NA, color = "#B85450", stat="identity", position="stack", width=0.8, linewidth=1) +
                scale_fill_manual( values = AlignGroupColors ) +
                scale_x_continuous(expand = c(0,0)) +
                labs( x = "Percentage (%)", y = "", fill = "") +
                theme(
                    panel.background = element_rect(fill = '#FFFFFF'),
                    panel.grid.major = element_line(color = '#e0e0e0', linetype = 'dotted'),
                    axis.line.x  = element_line(colour = "#6f6d6d"),
                    axis.ticks.y = element_blank(),
                    axis.text.x  = element_text(family="Helvetica", face="bold", size=10, hjust=0.5, vjust=0.5, margin = margin(t = 0)),
                    axis.text.y  = element_text(family="Helvetica", face="bold", size=10, hjust=1,   vjust=0.5, margin = margin(r = 0)),
                    axis.title   = element_text(family="Helvetica", face="bold", size=12, hjust=0.5, vjust=0.5, margin = margin(t = 10, r = 10)),
                    legend.text  = element_text(family="Helvetica", size=8),
                    plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
                )
            #==================================================================#
        }        
        return(RegionRatioPlot)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| COVERAGE BIAS HISTOGRAM |----------------------------------------------------------------------------------------------------------------------#
#' @description create RNAseq coverage bias plot
#' @param SeqFolderId  REQUIRED. analysis batch id ( seq_folder )
#' @param SampleInfo   REQUIRED. sample info table
#' @param SaveDir      folder for png save. default = data_figures 
#' @note  data will be downloaded from database automatically, using database connection params in global report params.     
#' @note  result object has no plot data. this function only generate plots and save them.
#' @export 
    WTS_plot.Coverage.Bias.Histogram <- function( SeqFolderId=NULL, SampleInfo=NULL, SaveDir="data_figures",
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_data
    )
    {
        require(ggplot2)
        require(RMySQL)
        require(cowplot)
        require(Cairo)
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #-----------------------------------------------------------------------
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        #-----------------------------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! No analysis batch ID ( 'seq_folder' ). REQUIRED. please check again. STOPPED.") }
        if( is.null(SampleInfo)  ){ stop("|---!!! No sample info table REQUIRED. please check again. STOPPED.")             }
        # read data ------------------------------------------------------------
        dbCon         <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        cov_hist_data <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,cov_bias_hist_value FROM qc_report WHERE seq_folder = '%s'", SeqFolderId))
        dbDisconnect(dbCon)
        #-----------------------------------------------------------------------
        cov_hist_data$sample_name  <- SampleInfo[match(cov_hist_data$seq_id, SampleInfo$seq_id), "sample_name" ]
        cov_hist_data$sample_group <- SampleInfo[match(cov_hist_data$seq_id, SampleInfo$seq_id), "sample_group"]
        #-----------------------------------------------------------------------
        sample_group_index <- indexing( cov_hist_data$sample_group )
        #-----------------------------------------------------------------------
        for( k in 1:length(sample_group_index) )
        {
            rep_index <- indexing( SampleInfo[which(SampleInfo$sample_group == sample_group_index[k]), "sample_name" ] )
            rep_bias_plot_list <- lapply(rep_index, function(sid) 
            {
                hist_data <- cov_hist_data %>% filter( sample_name == sid )
                plot_data <- data.frame(
                    index = 0:100,
                    value = as.numeric(unlist(strsplit(hist_data$cov_bias_hist_value, ";")))
                )
                #==============================================================#
                bias_plot <- ggplot(plot_data, aes( x = index, y = value ) ) +
                    geom_line( linewidth = 0.6, color = "#b4d1a3", alpha=1 ) +
                    geom_point( size = 1.1, pch=19, color = "#82b366", alpha=0.7 ) +
                    coord_cartesian(xlim=c(0,100)) +
                    labs( x = "Normalized Position", y = "Normalized Coverage" ) +
                    ggtitle( sid ) +
                    geom_hline( yintercept = 1, linewidth = 0.5, linetype="dashed", color="#D79B00") +
                theme(
                    panel.background = element_rect(fill = '#FFFFFF', color = "#2d2d2d"),
                    panel.grid.major = element_line(color = '#e0e0e0', linetype = 'dotted'),
                    axis.line        = element_line(colour = "#2d2d2d"),
                    axis.text        = element_text(family="Helvetica", face="bold", size=8),
                    axis.title       = element_text(family="Helvetica", face="bold", size=9),
                    plot.title       = element_text(family="Helvetica", face="bold", size=10, hjust=0.5),
                    legend.position  = "none",
                    plot.margin      = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
                )
                #==============================================================#    
                return(bias_plot)  
            })
            rep_lv_bias_plot <- plot_grid( plotlist = rep_bias_plot_list, nrow = 1 )
            #-------------------------------------------------------------------
            rep_lv_plot_width  <- 2.5 * length(rep_index)
            rep_lv_plot_height <- 2
            #==================================================================#
            ggsave( rep_lv_bias_plot, file=sprintf("%s/coverage_bias_plot_%s.png", SaveDir, k), 
                width = rep_lv_plot_width, height = rep_lv_plot_height, unit = 'in', dpi = 150
            )
            #==================================================================#
        }
        return(NULL)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| PCA FEATURE CONTRIBUTION TABLE |---------------------------------------------------------------------------------------------------------------#
#' @description PCA feature contribution stats summary table
#' @param SeqFolderID  REQUIRED. analysis batch id ( seq_folder )
#' @param FontSize     font size. default = 10.
#' @export 
    WTS_report.PCA.Feature.Contrib.Table <- function( SeqFolderId=NULL, FontSize=10,
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_data
    )
    {
        require(RMySQL)
        require(dplyr)
        require(knitr)
        require(kableExtra)
        #-----------------------------------------------------------------------
        dbCon      <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        pca_stats  <- dbGetQuery(dbCon, sprintf("SELECT * FROM pca_analysis_feature_stats WHERE seq_folder = '%s'", SeqFolderId ))
        dbDisconnect(dbCon)
        #======================================================================#
        pca_feature_table <- pca_stats %>% dplyr::select(
            c("comp","var_contrib","var_pct","var_cum_pct","gene_symbol","gene_type")
        ) %>% dplyr::rename(
            'PCA Component'               = comp,
            'Feature ID'                  = var_contrib,
            'Feature Contribution'        = var_pct,
            'Cumulative Contribution'     = var_cum_pct,
            'Gene Symbol'                 = gene_symbol,
            'Gene Type'                   = gene_type
        ) %>% kbl( escape = FALSE, align = 'c' ) %>% 
            kable_styling( bootstrap_options = "condensed", full_width = TRUE, font_size = FontSize ) %>%
            column_spec(c(1),     width_min="3em",  width="4em") %>%
            column_spec(c(2),     width_min="4em",  width="4em") %>%
            column_spec(c(3:4),   width_min="8em",  width="10em") %>%
            column_spec(c(5),     width_min="4em",  width="8em") %>%  
            column_spec(c(6),     width_min="4em",  width="12em") %>%
            row_spec( 0, font_size = FontSize, bold = TRUE, background = "#82B366", color = "#FFFFFF", 
                extra_css = "border-bottom: 1px solid; border-bottom-color: #4c4c4c" 
            )
        #======================================================================#
        return(pca_feature_table)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| ANALYSIS INFO SUMMARY TABLE |------------------------------------------------------------------------------------------------------------------#
#' @description  analysis sample groups and sample names (or sample ids)
#' @param AnalysisInfo  REQUIRED. information table of analysis pair sample groups
#' @param AnalysisId    REQUIRED. analysis id
#' @param SampleInfo    sample info table. if NULL, seq-id will be used as sample name.
#' @export
    WTS_report.Analysis.Info.Summary <- function( AnalysisInfo=NULL, AnalysisId=NULL, SampleInfo=NULL )
    {
        require(dplyr)
        require(knitr)
        require(kableExtra)
        #-----------------------------------------------------------------------
        if(is.null(AnalysisInfo) ){ stop("No analysis info table. REQUIRED. check again. STOPPED.") }
        if(is.null(AnalysisId)   ){ stop("No analysis id. REQUIRED. check again. STOPPED.") }
        #-----------------------------------------------------------------------
        ainfo_set <- AnalysisInfo %>% filter( analysis_id == AnalysisId )
        #-----------------------------------------------------------------------
        if( !is.null(SampleInfo) )
        { ainfo_set$sample_name <- SampleInfo[match(ainfo_set$sample_id, SampleInfo$seq_id), "sample_name"] }
        #-----------------------------------------------------------------------
        ainfo_table <- rbind(
            data.frame( 
                index  = "TARAGET", 
                group  = unique(ainfo_set[which(ainfo_set$analysis_group == "TREAT"), "sample_group"]),
                sample = paste(ainfo_set[which(ainfo_set$analysis_group == "TREAT"), "sample_name"], collapse=", ")
            ),
            data.frame( 
                index  = "CONTROL", 
                group  = unique(ainfo_set[which(ainfo_set$analysis_group == "CTRL"), "sample_group"]),
                sample = paste(ainfo_set[which(ainfo_set$analysis_group == "CTRL"), "sample_name"], collapse=", ")
            )
        ) %>% as.data.frame()
        colnames(ainfo_table) <- NULL
        #======================================================================#
        ainfo_table_out <- ainfo_table %>% kbl( escape = FALSE, align = c('c','c','l'), col.names = c(" ","Sample Group", "Sample Name") ) %>% 
            kable_styling( bootstrap_options = "condensed", full_width = TRUE, font_size = 9 ) %>% 
            column_spec(1, background = c("#D79B00","#6C8EBF"), color = c("#FFFFFF"), width_min = "2em", width = "3em" ) %>%
            column_spec(2, bold = TRUE, width_min = "8em", width = "8em", background = c("#FFE6CC","#DAE8FC") ) %>% 
            column_spec(3, width_min = "12em", width = "25em") %>% 
            row_spec(0, align = 'c')
        #======================================================================#
        return(ainfo_table_out)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| DEG CUT-OFF STATS TABLE |----------------------------------------------------------------------------------------------------------------------#
#' @description DEG cut-off and adjustment stats table
#' @param SeqFolderId REQURIED. analysis batch id ( seq_folder )
#' @param AnalysisId  REQUIRED. analysis id
#' @export 
    WTS_report.DEG.Stats.Table <- function( SeqFolderId=NULL, AnalysisId=NULL,
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_data
    )
    {
        require(RMySQL)
        require(dplyr)
        require(knitr)
        require(kableExtra)
        # seq-folder check -----------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! no seq-folder-id. REQUIRED. please check again. STOPPED.") }
        # analysis id check ----------------------------------------------------
        if( is.null(AnalysisId) ){ stop("|---!!! no analysis-id. REQUIRED. please check again. STOPPED.") }
        # load data from database ----------------------------------------------
        dbCon     <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        deg_stats <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_params WHERE seq_folder = '%s' AND analysis_id = '%s'", SeqFolderId, AnalysisId))
        dbDisconnect(dbCon)
        #-----------------------------------------------------------------------
        if( deg_stats$cutoff_adjust == 1 )
        {
            param_1st <- ifelse(
                deg_stats$cutoff_method == "cutoff",
                paste0( "|logFC| > ", deg_stats$cutoff_logfc, ", P-value < ", deg_stats$cutoff_pval ),
                paste0( "DEG Ratio < ", deg_stats$cutoff_ratio )
            )
            #==================================================================#
            deg_param_table <- data.frame( check.names=FALSE,
                UP   = deg_stats$up_genes_2nd,
                DOWN = deg_stats$down_genes_2nd,
                'DEG Params' = sprintf("DEG Ratio < %s", deg_stats$adjust_cutoff ),
                'Gene Space' = deg_stats$gene_space,
                Adjustment = "YES",
                'cutoff params' = param_1st,
                up   = deg_stats$up_genes_1st,
                down = deg_stats$down_genes_1st
            ) %>% kbl( escape = FALSE, align = c('c') ) %>% 
                kable_styling( bootstrap_options = "condensed", full_width = TRUE, font_size = 10 ) %>% 
                column_spec(c(1,2,7,8), width_min="2em", width="3em") %>%
                column_spec(c(4,5), width_min="6em", width="7em") %>%
                column_spec(c(3,6), width_min="8em", width="12em") %>%
                add_header_above(header = c("FINAL DEGs"=3, " "=2, "Initial DEGs "=3 ), font_size = 11, escape = FALSE,  line_sep = 8, bold = TRUE )
            #==================================================================#
        }else{
            param_1st <- ifelse(
                deg_stats$cutoff_method == "cutoff",
                paste0( "|logFC| > ", deg_stats$cutoff_logfc, ", P-value < ", deg_stats$cutoff_pval ),
                paste0( "DEG Ratio < ", deg_stats$cutoff_ratio )
            )
            #==================================================================#
            deg_param_table <- data.frame(
                UP = deg_stats$up_genes_1st,
                DOWN = deg_stats$down_genes_1st,
                'DEG Params' = param_1st,
                'Gene Space' = deg_stats$gene_space,
                Adjustment = "NO",
                'cutoff params' = "-",
                up   = "-",
                down = "-"
            ) %>% kbl( escape = FALSE, align = c('c') ) %>% 
                kable_styling( bootstrap_options = "condensed", full_width = TRUE, font_size = 10 ) %>% 
                column_spec(c(1,2,7,8), width_min="2em", width="3em") %>%
                column_spec(c(4,5), width_min="6em", width="7em") %>%
                column_spec(c(3,6), width_min="8em", width="12em") %>%
                add_header_above(header = c("FINAL DEGs"=3, " "=2, "Initial DEGs "=3 ), font_size = 11, escape = FALSE,  line_sep = 8, bold = TRUE )
            #==================================================================#
        }
        return(deg_param_table)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| TOP10 DEG LIST TABLE |-------------------------------------------------------------------------------------------------------------------------#
#' @description top 10 deg gene list ( up and down each )
#' @param SeqFolderId  REQUIRED. analysis batch id ( seq_folder )
#' @param AnalysisId   REQUIRED. analysis id
#' @param GeneInfo     REQUIRED. gene information table.
#' @export
    WTS_report.Top10.DEG.List.Table <- function( SeqFolderId=NULL, AnalysisId=NULL, GeneInfo=NULL,        
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_data
    )
    {
        require(RMySQL)
        require(dplyr)
        require(knitr)
        require(kableExtra)
        # seq-folder check -----------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! no seq-folder-id. REQUIRED. please check again. STOPPED.") }
        # analysis id check ----------------------------------------------------
        if( is.null(AnalysisId) ){ stop("|---!!! no analysis-id. REQUIRED. please check again. STOPPED.") }
        # gene info ------------------------------------------------------------
        if( is.null(GeneInfo) ){ stop("|---!!! no gene info table. REQUIRED. please check again. STOPPED.") }
        # load data from database ----------------------------------------------
        dbCon   <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        deg_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_result WHERE seq_folder = '%s' AND analysis_id = '%s'", SeqFolderId, AnalysisId))
        dbDisconnect(dbCon)
        #-----------------------------------------------------------------------
        deg_score         <- sapply(unlist(strsplit(deg_res$deg_score, ";")),  function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
        names(deg_score)  <- sapply(unlist(strsplit(deg_res$deg_score, ";")),  function(z) unlist(strsplit(z, ","))[1] ) 
        deg_logfc         <- sapply(unlist(strsplit(deg_res$deg_logfc, ";")),  function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
        names(deg_logfc)  <- sapply(unlist(strsplit(deg_res$deg_logfc, ";")),  function(z) unlist(strsplit(z, ","))[1] ) 
        deg_pvalue        <- sapply(unlist(strsplit(deg_res$deg_pvalue, ";")), function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
        names(deg_pvalue) <- sapply(unlist(strsplit(deg_res$deg_pvalue, ";")), function(z) unlist(strsplit(z, ","))[1] ) 
        deg_tag           <- sapply(unlist(strsplit(deg_res$deg_tag, ";")),    function(z) unlist(strsplit(z, ","))[2] ) 
        names(deg_tag)    <- sapply(unlist(strsplit(deg_res$deg_tag, ";")),    function(z) unlist(strsplit(z, ","))[1] ) 
        #-----------------------------------------------------------------------
        DEG_RES <- data.frame(
            ens_geneid      = names(deg_score),
            log2_foldchange = deg_logfc,
            pvalue          = deg_pvalue,
            score           = deg_score,
            deg             = deg_tag,
            GeneInfo[match(names(deg_score), GeneInfo$ens_geneid), c("symbol","entrez","locus_type","hgnc_id")] 
        )
        #-----------------------------------------------------------------------
        DEG_RES$locus_type <- gsub("gene with protein product", "protein", DEG_RES$locus_type )
        DEG_RES$locus_type <- gsub("RNA, long non-coding",      "lncRNA",  DEG_RES$locus_type )
        DEG_RES[is.na(DEG_RES)] <- ""
        #-----------------------------------------------------------------------
        list_up   <- DEG_RES %>% filter( deg == "UP"   ) %>% arrange( dplyr::desc(score) ) %>% slice(1:10) %>% dplyr::select(c("symbol","ens_geneid","entrez","locus_type"))     
        list_down <- DEG_RES %>% filter( deg == "DOWN" ) %>% arrange( dplyr::desc(score) ) %>% slice(1:10) %>% dplyr::select(c("symbol","ens_geneid","entrez","locus_type")) 
        rownames(list_up) = rownames(list_down) = NULL
        if( nrow(list_up) < 10 ){
            need_rows <- 10-nrow(list_up)
            for( j in 1:need_rows ){ list_up <- rbind( list_up, c("&nbsp;", "&nbsp;", "&nbsp;", "&nbsp;") ) }
        }
        if( nrow(list_down) < 10 ){
            need_rows2 <- 10-nrow(list_down)
            for( j in 1:need_rows2 ){ list_down <- rbind( list_down, c("&nbsp;", "&nbsp;", "&nbsp;", "&nbsp;") ) }
        }
        #-----------------------------------------------------------------------
        top10_deg_list <- data.frame( check.names = FALSE, Rank = 1:10, list_up, blk = "&nbsp;" , list_down ) 
        colnames(top10_deg_list) <- NULL
        top10_deg_list <- top10_deg_list[apply(top10_deg_list[,-1], 1, function(x) ifelse(all(x=="&nbsp;"), FALSE, TRUE) ), ]
        #======================================================================#
        if( nrow(top10_deg_list) > 0 ){
            deg_list_column_names <- c("Rank","Gene","Ensembl","Entrez","Gene Type"," ","Gene","Ensembl","Entrez","Gene Type")
            top10_deg_table <- top10_deg_list %>% kbl( escape = FALSE, align = 'c', col.names = deg_list_column_names ) %>%
                kable_styling( bootstrap_options = "condensed", full_width = TRUE, font_size = 9 ) %>% 
                column_spec(c(1),    width_min="2em",  width="2em") %>%
                column_spec(c(2),    width_min="4em",  width="4em", background = "#F8CECC", color = "#C73500", bold = TRUE ) %>%
                column_spec(c(7),    width_min="4em",  width="4em", background = "#D5E8D4", color = "#2D7600", bold = TRUE ) %>%
                column_spec(c(3,8),  width_min="5em",  width="5em") %>%
                column_spec(c(4,9),  width_min="3em",  width="3em") %>%  
                column_spec(c(5,10), width_min="5em",  width="7em") %>%
                column_spec(c(6),    width_min="1em",  width="1em", extra_css = "border-bottom: 0px solid; border-top: 0px solid" ) %>%
                add_header_above(header = c(" "=1, "UP-Regulated"=4, " "=1, "DOWN-Regulated"=4), font_size = 11, escape = FALSE,  line_sep = 8, bold = TRUE )       
        }else{
            top10_deg_table <- '<br><div style="color: #000000; font-size: 12px; font-weight: bold;">No DEGs.</div><br>'
        }
        #======================================================================#
        return(top10_deg_table)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| PATHWAY ANALYSIS TOP10 PATHWAYS |--------------------------------------------------------------------------------------------------------------#
#' @description 
#' @param SeqFolderId     REQUIRED. analysis batch id ( seq_folder )
#' @param AnalysisId      REQUIRED. analysis id
#' @param AnalysisMethod  REQUIRED. analysis method. c("ORA","GO","GSEA"). 
#' @param PathwaySetID    REQUIRED. pathway set id.
#' @param DegType         pathway analysis input deg type. c("UP","DOWN","ALL_DEGs"). if AnalysisMethod == "GSEA", this param will be ignored.
#' @param FontSize        table font size. default = 10
#' @param 
#' @export
    WTS_report.Top10.Pathway.List.Table <- function( SeqFolderId=NULL, AnalysisId=NULL, AnalysisMethod=NULL, PathwaySetID=NULL, DegType=NULL, FontSize=10,
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_data    
    )
    {
        require(RMySQL)
        require(dplyr)
        require(knitr)
        require(kableExtra)
        # seq-folder check -----------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! No SeqFolderID. REQUIRED.") }
        # analysis id check ----------------------------------------------------
        if( is.null(AnalysisId) ){ stop("|---!!! No AnalysisID. REQUIRED.") }
        # pathway analysis input deg type check --------------------------------
        if( is.null(AnalysisMethod) ){
            stop("|---!!! No AnalysisMethod. REQUIRED.")
        }else{
            if( AnalysisMethod != "GSEA" )
            { if( is.null(DegType) ){ stop("|---!!! No DEG-Type. REQURIED when AnalysisMethod is not 'GSEA'.") } }
        }
        # pathway set id check -------------------------------------------------
        if ( is.null(PathwaySetID) ){ stop("|---!!! No PathwaySetID. REQUIRED.")}
        # load data from database ----------------------------------------------
        dbCon   <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        if( AnalysisMethod == "ORA" ){  
            top10_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM pathway_analysis_top200_result WHERE seq_folder = '%s' AND deg_analysis_id = '%s' AND input_geneset = '%s'", SeqFolderId, AnalysisId, DegType))
        }else if( AnalysisMethod == "GO" ){
            top10_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM go_enrich_analysis_top200_result WHERE seq_folder = '%s' AND deg_analysis_id = '%s' AND input_geneset = '%s'", SeqFolderId, AnalysisId, DegType))
            top10_res$pathway <- gsub("^GOBP_", "", top10_res$pathway)
            top10_res$pathway <- gsub("_", " ", top10_res$pathway)
            top10_res$ext_id  <- gsub("\\:", "", top10_res$ext_id)
        }else if( AnalysisMethod == "GSEA" ){
            top_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM gsea_top200_result WHERE seq_folder = '%s' AND analysis_id = '%s'", SeqFolderId, AnalysisId))
            top10_res <- rbind(
                top_res %>% filter( NES > 0, analysis_pathset_id == PathwaySetID ) %>% arrange(dplyr::desc(score)) %>% dplyr::slice(1:5),
                top_res %>% filter( NES < 0, analysis_pathset_id == PathwaySetID ) %>% arrange(dplyr::desc(abs(score))) %>% dplyr::slice(1:5)
            )
        }
        dbDisconnect(dbCon)
        # color ----------------------------------------------------------------
        if(       is.null( DegType)     ){ BgColor <- "#E1D5E7"
        }else if( DegType == "UP"       ){ BgColor <- "#F8CECC" 
        }else if( DegType == "DOWN"     ){ BgColor <- "#D5E8D4" 
        }else if( DegType == "ALL_DEGs" ){ BgColor <- "#FFDCB8" }
        # filter by pathwaysetid -----------------------------------------------        
        top10_res <- top10_res %>% filter( analysis_pathset_id == PathwaySetID )
        #=======================================================================
        if( nrow(top10_res) == 0 )
        {
            top10_table <- data.frame(v1="No Results") %>% kbl( escape = FALSE, align = c('l'), col.names=NULL ) %>% 
                kable_styling( bootstrap_options = "condensed", full_width  =TRUE, font_size = FontSize ) 
        }else{
            if( AnalysisMethod != "GSEA" )
            {
                top10_table <- top10_res %>% arrange(dplyr::desc(score)) %>% slice(1:10) %>% mutate( 
                        Rank = 1:10, 
                        Genes = paste0( enrich_genes, "/", pathway_size ),
                        pvalue = as.character(signif(pvalue, 2)),
                        enrich_factor = as.character(round(enrich_factor, 2))
                    ) %>% dplyr::select(
                        c("Rank", "pathway", "Genes", "score", "pvalue","enrich_factor","data_source","ext_id" )
                    ) %>% dplyr::rename( 
                        Pathway = pathway, Score = score, 'P-value'=pvalue, 'Enrichment Factor'=enrich_factor, 'Pathway Source'=data_source, 'Source PathID'=ext_id
                    ) %>% kbl( escape = FALSE, align = c('c','l', rep('c', 6)) ) %>% 
                    kable_styling( bootstrap_options = "condensed", full_width  =TRUE, font_size = FontSize ) %>%
                    row_spec(0, background = BgColor ) %>%
                    column_spec(1, width_min="1em", width="1em" ) %>%
                    column_spec(2, width_min="15em", width="17em" ) %>%
                    column_spec(c(3,4,5,6), width_min="2em", width="3em" ) %>%
                    column_spec(c(7), width_min="4em", width="4em" ) %>%
                    column_spec(c(8), width_min="5em", width="5em" ) 
            }else{
                top10_table <- top10_res %>% mutate( 
                        Rank  = c(1:5, 1:5), 
                        pval  = as.character(signif(pval, 2)),
                        padj  = as.character(signif(padj, 2)),
                        score = as.character(round(score, 2)),
                        ES    = as.character(round(ES, 3)),
                        NES   = as.character(round(NES, 3))
                    ) %>% dplyr::select(
                        c("Rank", "pathway", "hits", "score", "pval", "padj","ES","NES","data_source","ext_id" )
                    ) %>% dplyr::rename( 
                        Pathway = pathway, Score = score, Genes=hits, 'P-value'=pval, FDR=padj, 'Pathway Source'=data_source, 'Source PathID'=ext_id
                    ) %>% kbl( escape = FALSE, align = c('c','l', rep('c', 6)) ) %>% 
                    kable_styling( bootstrap_options = "condensed", full_width  =TRUE, font_size = FontSize ) %>%
                    row_spec(0, background = "#E1D5E7" ) %>%
                    column_spec(1, width_min="1em", width="1em", background=c(rep("#F8CECC",5),rep("#D5E8D4",5)) ) %>%
                    column_spec(2, width_min="15em", width="18em" ) %>%
                    column_spec(c(3,4,5,6,7,8), width_min="3em", width="4em" ) %>%
                    column_spec(c(9,10), width_min="4em", width="5em" ) 
            }
        }
        #=======================================================================
        return(top10_table)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---| FUGION GENES and SPLICE VARIANTS |-------------------------------------------------------------------------------------------------------------#
#' @description create consensus fusion gene result
#' @param SeqFolderId       REQUIRED. analysis batch id ( seq_folder )
#' @param SampleInfo        sample information table
#' @param FusionGeneResult  fusion gene result raw-data
#' @export 
    WTS_report.FusionGenes.Table <- function( SeqFolderId=NULL, SampleInfo=NULL, FusionGeneResult=NULL )
    {
        require(dplyr)
        require(RMySQL)
        # function -------------------------------------------------------------
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        # input check ----------------------------------------------------------
        if (is.null(SeqFolderId)      ){ stop("|---!!! No SeqFolderID found. REQURIED.") }
        if (is.null(SampleInfo)       ){ stop("|---!!! No sample info table found. REQURIED.") }
        if (is.null(FusionGeneResult) ){ stop("|---!!! No Fusion-Gene analysis result data found. REQURIED.") }
        # filtering result -----------------------------------------------------
        fusion_gene_data <- data.frame( FusionGeneResult, sample_group = SampleInfo[match(FusionGeneResult$seq_id, SampleInfo$seq_id), "sample_group" ] )

        sample_group_index <- unique(fusion_gene_data$sample_group)
        FUSION_GENE_RESULT <- data.frame()

        for( sgid in sample_group_index )
        {
            SelCols    <- c("geneA","geneA_strand","geneA_breakpoint","geneB","geneB_strand","geneB_breakpoint","type_arriba","analysis_method")
            sg_fg_data <- fusion_gene_data %>% filter( sample_group == sgid ) %>% mutate( fgidx = paste( geneA, geneB, sep="_" ) )

            if( "starfusion,arriba" %in% sg_fg_data$analysis_method )
            {
                sidx       <- indexing( SampleInfo[which(SampleInfo$sample_group == sgid), "seq_id"] )
                as_fg_list <- Reduce(intersect, lapply(sidx, function(sx) sg_fg_data[which(sg_fg_data$seq_id == sx & sg_fg_data$analysis_method == "starfusion,arriba"), "fgidx"] ))
                as_fg_list <- as_fg_list[ grep("IGL@", as_fg_list, invert=T) ]
                if( length(as_fg_list) > 0 ){
                        as_fg_res <- unique(sg_fg_data[which(sg_fg_data$fgidx %in% as_fg_list & sg_fg_data$analysis_method == "starfusion,arriba"), SelCols ]) 
                        as_fg_res <- as_fg_res %>% mutate( index = 1:nrow(as_fg_res) )
                }else{ as_fg_res <- data.frame() }
            
                fg_res <- as_fg_res %>% mutate( sample_group = sgid )
            }else{
                sidx       <- indexing( SampleInfo[which(SampleInfo$sample_group == sgid), "seq_id"] )
                ab_fg_list <- Reduce(intersect, lapply(sidx, function(sx) sg_fg_data[which(sg_fg_data$seq_id == sx & sg_fg_data$analysis_method == "arriba"), "fgidx"] ))
                ab_fg_list <- ab_fg_list[ grep("IGL@", ab_fg_list, invert=T) ]

                if( length(ab_fg_list) > 0 ){
                        ar_fg_res <- unique(sg_fg_data[which(sg_fg_data$fgidx %in% ab_fg_list & sg_fg_data$analysis_method == "arriba"), SelCols ])
                        ar_fg_res <- ar_fg_res %>% mutate( index = 1:nrow(ar_fg_res) )
                }else{ ar_fg_res <- data.frame() }

                fg_res <- ar_fg_res %>% mutate( sample_group = sgid )
            }

            FUSION_GENE_RESULT <- rbind(FUSION_GENE_RESULT, fg_res)
        }
        return(FUSION_GENE_RESULT)
    }  
#--------------------------------------------------------------------------------------------------#
#' @description create printable table of fusion-gene result
#' @param FusionGeneReultTable  fusion-gene result table formatted to report (load from database)
#' @param SampleGroupId         sample group id
#' @param FontSize              font size. default = 9 
#' @export 
    WTS_report.FusionGene.Result.Table <- function( FusionGeneReultTable=NULL, SampleGroupId=NULL, FontSize=9 )
    {
        require(dplyr)
        require(knitr)
        require(kableExtra)
        #-----------------------------------------------------------------------
        fg_res <- FusionGeneReultTable %>% filter( sample_group == SampleGroupId )
        if( nrow(fg_res) > 0 )
        {
            fg_res$analysis_method <- gsub("arriba", "AR", fg_res$analysis_method)
            fg_res$analysis_method <- gsub("starfusion", "SF", fg_res$analysis_method)
            fg_res$geneA <- gsub(",", " ", fg_res$geneA)
            fg_res$geneB <- gsub(",", " ", fg_res$geneB)

            fg_res_table <- fg_res %>% dplyr::select(
                c("geneA","geneA_strand","geneA_breakpoint","geneB","geneB_strand","geneB_breakpoint","type_arriba","analysis_method","cancer_related",
                "kinase","oncogene","tumorSuppressor","transcrptionFactor","non_cancer_disease")
            )
            colnames(fg_res_table) <- NULL
            ReportVisColnames <- c("Gene","Strand","BreakPoint","Gene","Strand","BreakPoint","Fusion Type","Swtool","CANCER","K","OC","TSG","TF","NCD") 
            ResN <- nrow(fg_res_table)

            fg_res_table_output <- fg_res_table %>% kbl( escape = FALSE, col.names = ReportVisColnames, align = c(rep("c",13), "l"), vlaign='middle') %>%
                kable_styling( bootstrap_options = "condensed", full_width  =TRUE, font_size = FontSize ) %>% 
                column_spec(c(1,4), width_min="3em", width="4em" ) %>%
                column_spec(c(2,5), width_min="1em", width="1em" ) %>%
                column_spec(c(3,6), width_min="3em", width="3em" ) %>%
                column_spec(c(7,8), width_min="2em", width="3em" ) %>%
                column_spec(c(9), width_min="2em", width="2em" ) %>%
                column_spec(c(10,11,12,13), width_min="1em", width="1em" ) %>%
                add_header_above(header = c("Gene A"=3, "Gene B"=3, " "=8), font_size = FontSize, escape = FALSE,  line_sep = 10, bold = TRUE ) %>%
                footnote( general = "SF: STAR-Fusion, AR: arriba, CANCER: Cancer-related, K: Kinase, OC: Oncogene, TSG: TumorSuppressorGene, TF: TranscriptionFactor, NCD: Non-Cancer Disease")       
        }else{
            NoRes           <- data.frame("  ", "No Fusion-Genes Found.", " " )
            colnames(NoRes) <- NULL
            fg_res_table_output <- NoRes %>% kbl( escape = FALSE, align = 'l') %>% 
                kable_styling( bootstrap_options = "condensed", full_width  =TRUE, font_size = 11  ) %>% 
                column_spec(c(1), width_min="1em", width="1em" ) %>%
                column_spec(c(2), width_min="10em", width="20em" ) %>%
                row_spec(1, bold=TRUE)
        }
        return(fg_res_table_output)
    }
#--------------------------------------------------------------------------------------------------#
#' @description create table of splice variants result for report
#' @param SpliceVariantData  splice variant analysis result table ( load from database )
#' @param SampleInfo         sample information table
#' @param FontSize           font size. default = 9
#' @export
    WTS_report.Splice.Variants.Table <- function( SpliceVariantData=NULL, SampleInfo=NULL, FontSize=9 )
    {
        require(dplyr)
        require(knitr)
        require(kableExtra)
        #-----------------------------------------------------------------------
        SpliceVariantData$sample_name <- SampleInfo[match(SpliceVariantData$seq_id, SampleInfo$seq_id), "sample_name"]
        SpliceVariantData$GENE        <- sapply(SpliceVariantData$genes, function(gn) unlist(strsplit(gn, "\\^"))[1] )
        SpliceVariantData$ens_geneid  <- sapply(SpliceVariantData$genes, function(gn) unlist(strsplit(gn, "\\^"))[2] )
        #-----------------------------------------------------------------------
        sv_res <- SpliceVariantData %>% dplyr::select(c("sample_name", "variant_name","GENE", "intron","uniq_mapped","multi_mapped","ens_geneid")) %>%
            dplyr::rename('SAMPLE ID'=sample_name, 'VARIANT'=variant_name, 'LOCATION'=intron, 'UNIQUE'=uniq_mapped, 'MULTI'=multi_mapped, 'ENSEMBL ID'=ens_geneid ) %>%
            kbl( escape = FALSE, align='c' ) %>%
            kable_styling( bootstrap_options = "condensed", full_width  =TRUE, font_size = FontSize ) %>% 
            column_spec(c(1),   width_min="8em", width="12em" ) %>%
            column_spec(c(2,7), width_min="6em", width="8em" ) %>%
            column_spec(c(3),   width_min="4em", width="4em" ) %>%
            column_spec(c(4),   width_min="10em", width="14em" ) %>%
            column_spec(c(5,6), width_min="3em", width="4em" ) %>%
            add_header_above(header = c(" "=4, "Mapped Reads"=2, " "=1), font_size = 10, escape = FALSE,  line_sep = 15, bold = TRUE ) 
        #-----------------------------------------------------------------------
        return(sv_res)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#' @description WTS SWTOOLS AND VERSIONS
#' @export 
    Appendix_swtools.versions.WTS <- function()
    {
        suppressPackageStartupMessages(library("knitr"))
        suppressPackageStartupMessages(library("kableExtra"))
        suppressPackageStartupMessages(library("dplyr"))
        #-----------------------------------------------------------------------
        swtool_ver <- data.frame(
            c(
                "FASTQ QC",             "&nbsp;",
                "FASTQ to BAM Process", "&nbsp;",
                "BAM QC",               "&nbsp;",
                "QUANTIFICATION",       "     ",
                "FUSION GENES",         "&nbsp;",
                "SPLICE VARIANTS",      "&nbsp;",
                "STANDARD ANALYSIS",    "&nbsp;",
                "ADVANCED ANALYSIS"
            ),
            c(" "),
            c(
                "fastqc (v0.12.1), fastp (v0.23.4)",      "&nbsp;",
                "STAR (v2.7.11)",                         "&nbsp;",
                "GATK4 (v4.4.0.0), RNA-SeQC (v2.4.2)",    "&nbsp;",
                "RSEM (v1.3.3), RNA-SeQC (v2.4.2)",       "      ",
                "STAR-Fusion (v0.13.0), Arriba (v2.4.0)", "&nbsp;",
                "CTAT Splicing (v0.0.2)",                 "&nbsp;",
                "R (v4.3.1), DESeq2 (v1.40.2), FactoMineR (v2.11), umap (v0.2.10)",           "&nbsp;",
                "R (v4.3.1), fgsea (v1.26.0)"                                   
            )
        )  
        colnames(swtool_ver) <- NULL

        swtool_ver_table <- swtool_ver %>% kbl( escape = FALSE, align = c('c','c','l') ) %>% 
            kable_styling( bootstrap_options = "condensed", full_width = TRUE, position = "left" ) %>%
            column_spec( 1, width_min = "10em", width = "15em", background = "#6C8EBF", color = "#FFFFFF" ) %>%
            column_spec( 2, width_min = "1em", width = "1em" ) %>%
            column_spec( 3, width_min = "25em", width = "50em" ) %>%
            row_spec( c(2,4,6,8,10,12,14), font_size = 4, background = "#FFFFFF", extra_css = "border-bottom: 0px solid; border-top: 0px solid" ) %>%
            row_spec( c(1,3,5,7,9,11,13,15), bold = TRUE, font_size = 12, extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-bottom-color: #DAE8FC" )
            
        return(swtool_ver_table)
    } 
#----------------------------------------------------------------------------------------------------------------------------------------------------#













# #---| DEG RESULT PLOTS |
# #' @description 
# #' @param 
# #' @param 
# #' @export  
#     WTS_report.DEG.Plots <- function( SeqFolderId=NULL, AnalysisId=NULL, AnalysisInfo=NULL,
#         db_host     = cnf$db$host,
#         db_user     = cnf$db$user,
#         db_password = cnf$db$pw,
#         db_port     = cnf$db$port,
#         db_name     = cnf$db$db_data
#     )
#     {
#         require(RMySQL)
#         require(dplyr)
#         require(ggplot2)
#         require(GGally)
#         # seq-folder check -----------------------------------------------------
#         if( is.null(SeqFolderId) ){ stop("|---!!! no seq-folder-id. REQUIRED. please check again. STOPPED.") }
#         # analysis id check ----------------------------------------------------
#         if( is.null(AnalysisId) ){ stop("|---!!! no analysis-id. REQUIRED. please check again. STOPPED.") }
#         #-----------------------------------------------------------------------
#         # load data from database ----------------------------------------------
#         dbCon   <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
#         deg_res <- dbGetQuery(dbCon, sprintf("SELECT * FROM deg_analysis_result WHERE seq_folder = '%s' AND analysis_id = '%s'", SeqFolderId, AnalysisId))
#         dbDisconnect(dbCon)
#         ##
#         deg_score         <- sapply(unlist(strsplit(deg_res$deg_score, ";")),  function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
#         names(deg_score)  <- sapply(unlist(strsplit(deg_res$deg_score, ";")),  function(z) unlist(strsplit(z, ","))[1] ) 
#         deg_logfc         <- sapply(unlist(strsplit(deg_res$deg_logfc, ";")),  function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
#         names(deg_logfc)  <- sapply(unlist(strsplit(deg_res$deg_logfc, ";")),  function(z) unlist(strsplit(z, ","))[1] ) 
#         deg_pvalue        <- sapply(unlist(strsplit(deg_res$deg_pvalue, ";")), function(z) as.numeric(unlist(strsplit(z, ","))[2]) ) 
#         names(deg_pvalue) <- sapply(unlist(strsplit(deg_res$deg_pvalue, ";")), function(z) unlist(strsplit(z, ","))[1] ) 
#         deg_tag           <- sapply(unlist(strsplit(deg_res$deg_tag, ";")),    function(z) unlist(strsplit(z, ","))[2] ) 
#         names(deg_tag)    <- sapply(unlist(strsplit(deg_res$deg_tag, ";")),    function(z) unlist(strsplit(z, ","))[1] ) 
        
#         DEG_RES <- data.frame(
#             ens_geneid      = names(deg_score),
#             log2_foldchange = deg_logfc,
#             pvalue          = deg_pvalue,
#             score           = deg_score,
#             deg             = deg_tag
#         )
#         # volcano plot ---------------------------------------------------------
#         VolcanoInput         <- DEG_RES[which(!is.na(DEG_RES$pvalue)), ]
#         VolcanoInput$y_value <- abs(VolcanoInput$score)
#         max_y                <- max(VolcanoInput$y_value) * 1.3
#         max_x                <- as.integer(max(abs(VolcanoInput$log2_foldchange))) * 1.2 
#         FoldChangeGuide      <- VolcanoInput %>% filter( deg %in% c("UP","DOWN") ) %>% .$log2_foldchange %>% abs() %>% min()
#         ScoreGuide           <- VolcanoInput %>% filter( deg %in% c("UP","DOWN") ) %>% .$y_value %>% min()
#         VolcanoInput$deg_type <- sapply( VolcanoInput$deg, function(g)
#         { 
#             if( g == "UP" ){ gt <- "UP-Regulation" }else if( g == "DOWN" ){ gt <- "DOWN-Regulation" }else{ gt <- "Non-DEGs" }
#             return(gt)
#         })
#         VolcanoInput$deg_type <- factor(VolcanoInput$deg_type, levels = c("UP-Regulation","DOWN-Regulation","Non-DEGs"))
#         GENE_COLORS <- c( "UP-Regulation"="#9f3122", "DOWN-Regulation"="#347019", "Non-DEGs"="#d6d6d6" )

#         VPLOT <- ggplot( VolcanoInput, aes(x=log2_foldchange, y=y_value, color=deg_type) ) + 
#             geom_point( size = 2, pch=19, alpha=0.6 ) +
#             labs( x="FoldChange", y="Score Metric", color="Genes" ) +
#             scale_color_manual( values = GENE_COLORS ) +
#             coord_cartesian( xlim = c(-max_x, max_x), ylim = c(0, max_y) ) +
#             geom_hline(yintercept=ScoreGuide, linetype=2, color='#1F78B4') +
#             geom_vline(xintercept=c( -FoldChangeGuide, FoldChangeGuide ), linetype=2, color='#1F78B4') +
#             theme(
#                 panel.background = element_rect(fill='white'), panel.grid.major = element_line(color = "#d6d6d6", linetype = 'dotted'),
#                 axis.text        = element_text(family="Helvetica", colour="#2d2d2d", size=10),
#                 axis.title       = element_text(family="Helvetica", colour="#2d2d2d", size=11),
#                 axis.line        = element_line(colour="#2d2d2d"),
#                 legend.title     = element_text(family="Helvetica", colour="#2d2d2d", size=10),
#                 legend.text      = element_text(family="Helvetica", colour="#2d2d2d", size=9),
#                 plot.margin      = unit(c(0.2,0.2,0.2,0.2), "cm")  
#             )
#         #-----------------------------------------------------------------------
        
#         # correlation scatter
#         # SeqIdList <- c(
#         #     sort(AnalysisInfo[which(AnalysisInfo$analysis_id == AnalysisId & AnalysisInfo$analysis_group == "CTRL"), "sample_id"]),
#         #     sort(AnalysisInfo[which(AnalysisInfo$analysis_id == AnalysisId & AnalysisInfo$analysis_group == "TREAT"), "sample_id"])
#         # )
#         # SEQ_IDS <- paste(paste0("'", SeqIdList, "'"), collapse=",")
#         # #-----------------------------------------------------------------------
#         # dbCon   <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
#         # EV <- dbGetQuery(dbCon, sprintf("SELECT seq_folder,seq_id,read_count FROM expr_values WHERE seq_id IN (%s)", SEQ_IDS)) %>% dplyr::rename( expr=3 )
#         # dbDisconnect(dbCon)
#         # EV_List <- setNames( EV$expr, nm=EV$seq_id )
#         # EV_List <- lapply( EV_List, function(ev) unlist(strsplit(ev, ";")) )
#         # EV_List <- lapply( EV_List, function(ev)
#         # {
#         #     names(ev) <- sapply(ev, function(y) unlist(strsplit(y, ","))[1])
#         #     nev <- sapply(ev, function(y) round(as.numeric(unlist(strsplit(y, ","))[2]),0))
#         #     return(nev)
#         # })
#         # EXPR_VALUES <- as.data.frame( EV_List )
#         # CPM <- edgeR::cpm(EXPR_VALUES[, SeqIdList], log=TRUE, normalized.lib.sizes=TRUE, prior.count=1) %>% as.data.frame()
#         # CPM$clr = "cpm"
#         # ggpairs( CPM, columns = 1:4, aes(color = clr), diag = list( continuous = "blankDiag" ), upper = list( continuous = "points") )

#         return(VPLOT)
#     }
            





