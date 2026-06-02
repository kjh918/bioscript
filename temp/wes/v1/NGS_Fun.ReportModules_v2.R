
# "cnf" REQUIRED. 
# cnf = read_yaml( config_file )

#---/ LOAD NGS STATS FROM DATABASE /-----------------------------------------------------------------------------------------------------------------#
#' @description load stats data from database
#' @param SeqFolderID
#' @export 
    NGS_report.Load.Stats <- function( SeqFolderID=NULL, Src="report",
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = as.numeric(cnf$db$port),
        db_name     = cnf$db$db_data
    )
    {
        require(RMySQL)
        #-----------------------------------------------------------------------
        if( Src == "bam" ){ dbTable <- "qc_bam" }else{ dbTable <- "qc_report" }
        #-----------------------------------------------------------------------
        dbCon  <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        qcStat <- dbGetQuery(dbCon, sprintf("SELECT * FROM %s WHERE seq_folder = '%s'", dbTable, SeqFolderID))
        dbDisconnect(dbCon)
        return(qcStat)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ Order Info Table /-----------------------------------------------------------------------------------------------------------------------------#
#' @description create "ORDER INFORMATION TABLE"
#' @param SeqFolderId  analysis batch id ( = seq_folder )
#' @param ReportDate   report generation date. if NULL, today() will be used.
#' @param HeaderColor  category strip color. if NULL, default color will be used.
#' @export 
    NGS_report.Order.Info.Table <- function( SeqFolderId=NULL, ReportDate=NULL, HeaderColor=NULL, KOR=FALSE,
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_info
    )
    {
        require(lubridate)
        require(knitr)
        require(kableExtra)
        require(RMySQL)
        #-----------------------------------------------------------------------
        seq_folder_ids <- paste(paste0("'", unlist(strsplit(SeqFolderId, ",")), "'"), collapse=",")
        dbCon        <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        batch_info   <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqfolder_info WHERE seq_folder IN (%s)", seq_folder_ids))
        client_info  <- dbGetQuery(dbCon, sprintf("SELECT * FROM facility_info WHERE facility_id = '%s'", unique(batch_info$facility_id) ))
        sample_info  <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder IN (%s)", seq_folder_ids))
        dbDisconnect(dbCon)
        # analysis batch id check ----------------------------------------------
        if( is.null(SeqFolderId) ){ stop("NO SEQ_FODLER ID. REQUIRED.") }
        # header color ---------------------------------------------------------
        if( is.null(HeaderColor) ){ header_color <- "#82B366" }else{ header_color <- HeaderColor }
        # report date ----------------------------------------------------------
        if( is.null(ReportDate) ){ report_date <- today() }else{ report_date <- as_date(ReportDate) }
        # NA value replace -----------------------------------------------------
        if( is.na(unique(sample_info$date_sample_in)) ){ sample_in_date <- " " }else{ sample_in_date <- as_date(unique(sample_info$date_sample_in)) }
        if( is.na(unique(sample_info$date_sample_qc)) ){ sample_qc_date <- " " }else{ sample_qc_date <- as_date(unique(sample_info$date_sample_qc)) }

        # order info table -----------------------------------------------------
        if( KOR )
        {
            inst_name <- gsub("_", " ", client_info$facility_name_kor)
        }else{
            inst_name <- gsub("_", " ", client_info$facility_name)
        }

        order_info_table <- data.frame(
            c(  "CLIENT INSTITUTE"      , " " ,
                "TOTAL SAMPLES"         , " " ,      
                "DATE of SAMPLE RECEIPT", " " ,
                "DATE of SAMPLE QC"     , " " , 
                "DATE of NGS REPORT"    , " " ,
                "ORDER ID"              , " " ,
                "ANALYSIS INSTITUTE"
            ),
            c(" "),
            c(  paste0(" "," ", inst_name)                                , " ",
                paste0(" "," ", nrow(sample_info))                        , " ",      
                paste0(" "," ", as.character(sample_in_date))             , " ",
                paste0(" "," ", as.character(sample_qc_date))             , " ", 
                paste0(" "," ", as.character(report_date))                , " ",
                paste0(" "," ", unique(batch_info$ngs_order_id))          , " ",
                paste0(" "," ", "GENCURIX Inc.")                                                             
            ),
            c(" ")
        )
        colnames(order_info_table) <- NULL
        # create kable for report rmd ------------------------------------------ 
        ORDER_INFO_TABLE <- order_info_table %>% kbl( escape = FALSE, align = c('c','c','l','c') ) %>% 
            kable_styling( bootstrap_options = "basic", full_width = TRUE, position = "left" ) %>%
            column_spec( 1, width_min = "10em", width = "15em", background = header_color, color = "#FFFFFF" ) %>%
            column_spec( 2, width_min = "1em", width = "1em" ) %>%
            column_spec( 3, width_min = "25em", width = "50em" ) %>%
            row_spec( c(1,3,5,7,9,11,13), bold = TRUE, font_size = 12, extra_css = "border-top: 0px solid; border-bottom: 1px dotted; border-color: #D5E8D4" ) %>%
            row_spec( c(2,4,6,8,10,12), background = "#FFFFFF", extra_css = "border-bottom: 0px solid" )
        #-----------------------------------------------------------------------
        return(ORDER_INFO_TABLE)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
#---/ Sequencing Info Table /------------------------------------------------------------------------------------------------------------------------#
#' @description create "SEQUENCING INFORMATION TABLE"
#' @param SeqFolderId  analysis batch id ( = seq_folder )
#' @param Application  ngs application. c("WES","WTS","WGS","TSO"). default = "WES". manual input values are possible.
#' @param Material     ngs material. c("DNA","RNA","miRNA","DNA,RNA"). default = "DNA". manual input values are possible.
#' @param NgsLibrary   ngs construction library. c("twist","sureselect","rnaseq"). default = "twist". manual input values are possible.
#' @param ReadType     ngs read type. c("PE","SE"). default = "PE", manual input values are possible.
#' @param ReadLength   ngs read length. default = 151. manual input values are possible.
#' @param HeaderColor  category strip color. if NULL, default color will be used.
#' @export 
    NGS_report.NGS.Info.Table <- function( SeqFolderId=NULL, Application=NULL, Material=NULL, NgsLibrary=NULL,
        ReadType=NULL, ReadLength=NULL, Platform=NULL, HeaderColor=NULL,
        db_host     = cnf$db$host,
        db_user     = cnf$db$user,
        db_password = cnf$db$pw,
        db_port     = cnf$db$port,
        db_name     = cnf$db$db_info
    )
    {
        require(lubridate)
        require(knitr)
        require(kableExtra)
        require(RMySQL)
        dbCon        <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_password, db=db_name )
        batch_info   <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqfolder_info WHERE seq_folder = '%s'", SeqFolderId))
        client_info  <- dbGetQuery(dbCon, sprintf("SELECT * FROM facility_info WHERE facility_id = '%s'", batch_info$facility_id ))
        sample_info  <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SeqFolderId))
        dbDisconnect(dbCon)
        # analysis batch id check ----------------------------------------------
        if( is.null(SeqFolderId) ){ stop("NO SEQ_FODLER ID. REQUIRED.") }
        # header color ---------------------------------------------------------
        if( is.null(HeaderColor) ){ header_color <- "#82B366" }else{ header_color <- HeaderColor }
        # ngs library ----------------------------------------------------------
        if( is.null(NgsLibrary) ){ ngs_library <- "Twist Exome 2.0" 
        }else{ 
            if(       NgsLibrary == "twist"      ){ ngs_library <- "Twist Exome 2.0"
            }else if( NgsLibrary == "sureselect" ){ ngs_library <- "Agilent SureSelect V6"
            }else if( NgsLibrary == "rnaseq"     ){ ngs_library <- "Illumina RNA Prep with Enrichment"
            }else{                                  ngs_library <- NgsLibrary }
        }
        # ngs application ------------------------------------------------------
        if( is.null(Application) ){ ngs_application <- "Whole Exome Sequencing (WES)"
        }else{
            if(       Application == "WES"  ){ ngs_application <- "Whole Exome Sequencing (WES)" 
            }else if( Application == "WTS"  ){ ngs_application <- "Whole Transcriptome Sequencing (WTS, RNAseq)"
            }else if( Application == "WGS"  ){ ngs_application <- "Whole Genome Sequencing (WGS)" 
            }else if( Application == "TSO"  ){ ngs_application <- "TruSight Oncology 500 (TSO-500)" 
            }else{                             ngs_application <- Application }
        }
        # ngs material ---------------------------------------------------------
        if( is.null(Material) ){ ngs_material <- "DNA"
        }else{
            if(       Material == "DNA"    ){ ngs_material <- "DNA" 
            }else if( Material == "RNA"    ){ ngs_material <- "RNA" 
            }else if( Material == "miRNA"  ){ ngs_material <- "micro RNA" 
            }else if( Material == "DNA,RNA"){ ngs_material <- "DNA, RNA"
            }else{                            ngs_material <- Material }
        }
        # ngs seq read type ----------------------------------------------------
        if(       ReadType == "SE" ){ ngs_type <- "Single-end (SE)" 
        }else if( ReadType == "PE" ){ ngs_type <- "Paired-end (PE)"
        }else{                        ngs_type <- ReadType }
        # ngs read length ------------------------------------------------------
        if( is.null(ReadLength) ){ ngs_read_length <- "151" 
        }else{                     ngs_read_length <- ReadLength }
        # ngs platform ---------------------------------------------------------
        if( is.null(Platform) ){ ngs_platform <- "ILLUMINA Novaseq" 
        }else{                   ngs_platform <- Platform }
        # ngs info table -------------------------------------------------------
        ngs_seq_info_table <- data.frame(
            c(  "NGS APPLICATION" , "&nbsp;",
                "NGS MATERIAL"    , "&nbsp;",      
                "NGS LIBRARY"     , "&nbsp;",
                "NGS READ TYPE"   , "&nbsp;", 
                "READ LENGTH"     , "&nbsp;",
                "NGS PLATFORM" 
            ),
            c("&nbsp;"),
            c(  paste0("&nbsp;","&nbsp;", ngs_application), "&nbsp;",
                paste0("&nbsp;","&nbsp;", ngs_material)   , "&nbsp;",      
                paste0("&nbsp;","&nbsp;", ngs_library)    , "&nbsp;",
                paste0("&nbsp;","&nbsp;", ngs_type)       , "&nbsp;", 
                paste0("&nbsp;","&nbsp;", ngs_read_length), "&nbsp;",
                paste0("&nbsp;","&nbsp;", ngs_platform)  
            ),
            c("&nbsp;")
        )
        colnames(ngs_seq_info_table) <- NULL
        # create kable for report rmd ------------------------------------------
        NGS_SEQ_INFO_TABLE <- ngs_seq_info_table %>% kbl( escape = FALSE, align = c('c','c','l','c') ) %>% 
            kable_styling( bootstrap_options = "basic", full_width = TRUE, position = "left" ) %>%
            column_spec( 1, width_min = "10em", width = "15em", background = header_color, color = "#FFFFFF" ) %>%
            column_spec( 2, width_min = "1em", width = "1em" ) %>%
            column_spec( 3, width_min = "25em", width = "50em" ) %>%
            row_spec( c(1,3,5,7,9,11), bold = TRUE, font_size = 12, extra_css = "border-top: 0px solid; border-bottom: 1px dotted; border-color: #D5E8D4" ) %>%
            row_spec( c(2,4,6,8,10), background = "#FFFFFF", extra_css = "border-bottom: 0px solid" )
        #-----------------------------------------------------------------------
        return(NGS_SEQ_INFO_TABLE)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ GC-Ratio Plot Module /-------------------------------------------------------------------------------------------------------------------------#
#' @description create GC-ratio barplot 
#' @param SeqFolderId  REQURIED. analysis batch id (=seq_folder id)
#' @param SampleInfo   sample info table. if NULL, seq-id will be used as sample name.
#' @export
    NGS_report.GC.Contents.Barplot <- function( SeqFolderId=NULL, SampleInfo=NULL )
    {
        require(ggplot2)
        # seq folder id check --------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! No seq folder id found. REQUIRD. please check again. STOPPED.")}
        #-----------------------------------------------------------------------
        StatData <- NGS_report.Load.Stats( SeqFolderID = SeqFolderId, Src="report" )
        StatData$pct_gc_content <- round(StatData$pct_gc_content, 1)
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
        StatData <- StatData %>% dplyr::rename( 'Sample_ID' = seq_id, 'GC_Ratio' = pct_gc_content )
        #-----------------------------------------------------------------------
        gc_contents_plot <- ggplot( StatData, aes( x = GC_Ratio, y = Sample_ID ) ) +
            geom_bar( position="stack", stat="identity", width=0.8, fill="#D79B00" ) +
            coord_cartesian( xlim = c(0,100) ) +
            labs( x = "GC Ratio (%)", y = "" ) +
            scale_x_continuous( breaks = c(0,25,50,75,100), label=c("0","25","50","75","100"), expand=c(0,0) ) +
            geom_text(aes(label = paste0(GC_Ratio,"%")), position = position_stack(vjust = .5), size=4, color=c('white')) +
            theme(
                panel.background = element_rect(fill = '#FFFFFF'),
                panel.grid.major = element_line(color = '#e0e0e0', linetype = 'dotted'),
                axis.line.x  = element_line(colour = "#6f6d6d"),
                axis.ticks.y = element_blank(),
                axis.text.x  = element_text(family="Nunito", face="bold", size=10, hjust=0.5, vjust=0.5, margin = margin(t = 0)),
                axis.text.y  = element_text(family="Nunito", face="bold", size=8, hjust=1,   vjust=0.5, margin = margin(r = 0)),
                axis.title   = element_text(family="Nunito", face="bold", size=10, hjust=0.5, vjust=0.5, margin = margin(t = 10, r = 10)),
                legend.text  = element_text(family='Nunito'),
                legend.title = element_text(family='Nunito'),
                plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
            )
        #-----------------------------------------------------------------------
        return(gc_contents_plot)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
#---/ Duplicates-Ratio Plot Module /-----------------------------------------------------------------------------------------------------------------#
#' @description create duplicates-ratio barplot
#' @param SeqFolderId  REQURIED. analysis batch id (=seq_folder id)
#' @param SampleInfo   sample info table. if NULL, seq-id will be used as sample name.
#' @param Application  NGS application type. 'WES' for DNAseq and 'WTS' for RNAseq. default = "WES"
#' @export 
    NGS_report.Duplicates.Ratio.Barplot <- function( SeqFolderId=NULL, SampleInfo=NULL, Application=NULL )
    {
        require(RMySQL)
        require(ggplot2)
        require(dplyr)
        # seq folder id check --------------------------------------------------
        if( is.null(SeqFolderId) ){ stop("|---!!! No seq folder id found. REQUIRD. please check again. STOPPED.")}
        #-----------------------------------------------------------------------
        bar_colors <- c("Unique Reads"="#DAE8FC", "Duplicates"="#B85450")
        #-----------------------------------------------------------------------
        StatData <- NGS_report.Load.Stats( SeqFolderID = SeqFolderId, Src="report" )
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
        if( is.null(Application) ){ 
            StatData <- StatData[,c("seq_id","pct_duplicates")] %>% dplyr::rename( dup = 2 ) 
        }else{
            if(       Application == "WTS" ){ StatData <- StatData[,c("seq_id","pct_duplicates")] %>% dplyr::rename( dup = 2 ) 
            }else if( Application == "WES" ){ StatData <- StatData[,c("seq_id","pct_rds_dup")]    %>% dplyr::rename( dup = 2 ) }
        }
        #-----------------------------------------------------------------------
        StatData$dup = round(StatData$dup, 1)
        PlotData <- rbind(
            data.frame(sample_id = StatData$seq_id, group = "Unique Reads",  pct = 100-StatData$dup),
            data.frame(sample_id = StatData$seq_id, group = "Duplicates",    pct = StatData$dup    )
        )
        #-----------------------------------------------------------------------
        dup_plot <- ggplot( PlotData, aes( x = pct, y = sample_id, fill = factor(group, levels = c("Unique Reads","Duplicates")), group = sample_id ) ) + 
            geom_bar( stat = "identity" ) + 
            scale_x_continuous( breaks = c(0,25,50,75,100), label=c("0","25","50","75","100"), expand=c(0,0) ) +
            labs( x = "Duplicates Ratio (%)", y = "", fill = "" ) +
            scale_fill_manual( values = bar_colors ) +
            theme(
                panel.background = element_rect(fill = '#FFFFFF'),
                panel.grid.major = element_line(color = '#e0e0e0', linetype = 'dotted'),
                axis.line.x  = element_line(colour = "#6f6d6d"),
                axis.ticks.y = element_blank(),
                axis.text.x  = element_text(family="Nunito", face="bold", size=10, hjust=0.5, vjust=0.5, margin = margin(t = 0)),
                axis.text.y  = element_text(family="Nunito", face="bold", size=8, hjust=1,   vjust=0.5, margin = margin(r = 0)),
                axis.title   = element_text(family="Nunito", face="bold", size=10, hjust=0.5, vjust=0.5, margin = margin(t = 10, r = 10)),
                plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                legend.position = "bottom",
                legend.text  = element_text(family='Nunito'),
                legend.title = element_text(family='Nunito')
            ) +
            geom_text( 
                data = data.frame(sample_id = StatData$seq_id, group = "Unique Reads", pct = 100-StatData$dup, lbl= paste0(StatData$dup, "%") ), 
                aes(label = lbl), position = position_stack(vjust = .8), size=4, color="#B85450"
            ) 
        #-----------------------------------------------------------------------
        return(dup_plot)
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
#---/ Empty Plot Module /----------------------------------------------------------------------------------------------------------------------------#
#' @description draw blank plot for spacing
#' @note no any argument
#' @export 
    NGS_report.Blank.Plot <- function(){ ggplot() + theme_void() }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ Page Breaks V1 /-------------------------------------------------------------------------------------------------------------------------------#
#' @description page breaks with gencurix logo image (version 1)
#' @export 
    NGS_report.Page.Breaks <- function()
    {
        cat('<div style="page-break-after: always;"></div>')
        cat('<div><img style="float: right;" src="resources/gcx_report_logo.png" width="80" height="20"></div>')
        cat('<br>')
        cat('<hr>')
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
#---/ Page Breaks V2 /-------------------------------------------------------------------------------------------------------------------------------#
#' @description page breaks without gencurix logo image (version 2)
#' @export 
    NGS_report.Page.Breaks2 <- function()
    {
        cat('<div style="page-break-after: always;"></div>')
        cat('<br>')
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#

#---/ Appendix : swtools and version /---------------------------------------------------------------------------------------------------------------#
#' @description WES SWTOOLS AND VERSIONS
#' @export 
    Appendix_swtools.versions.WES <- function()
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
                "VARIANT CALL",         "     ",
                "ANNOTATION",           "&nbsp;",
                "HLA TYPING",           "&nbsp;",
                "MSI",                  "&nbsp;",
                "INDIVIDUAL MATCHING",  "&nbsp;",
                "CONTAMINATION QC"
            ),
            c(" "),
            c(
                "fastqc (v0.12.1), fastp (v0.23.4)",                                     "&nbsp;",
                "bwa (v0.7.17), Picard (v3.1.0), GATK3 (v3.8), GATK4 (v4.4.0.0)",        "&nbsp;",
                "GATK4 (v4.4.0.0), mosdepth (v0.3.6), alfred (v0.2.6), multiqc (v1.18)", "&nbsp;",
                "Mutect2 in GATK4 (v4.4.0.0)",                                           "      ",
                "ensembl-vep (release 110.1), vcf2maf (v1.6.21)",                        "&nbsp;",
                "Optitype (v1.3.3), HLA-LA (v1.0.3)",                                    "&nbsp;",
                "MANTIS (v1.0.5), MSIsensor2 (v0.1)",                                    "&nbsp;",
                "NGSCheckMate (v1.0.1)",                                                 "&nbsp;",
                "Fastq-Screen (v0.15.3)"
            )
        )  
        colnames(swtool_ver) <- NULL

        swtool_ver_table <- swtool_ver %>% kbl( escape = FALSE, align = c('c','c','l') ) %>% 
            kable_styling( bootstrap_options = "condensed", full_width = TRUE, position = "left" ) %>%
            column_spec( 1, width_min = "10em", width = "15em", background = "#6C8EBF", color = "#FFFFFF" ) %>%
            column_spec( 2, width_min = "1em", width = "1em" ) %>%
            column_spec( 3, width_min = "25em", width = "50em" ) %>%
            row_spec( c(2,4,6,8,10,12,14,16), font_size = 4, background = "#FFFFFF", extra_css = "border-bottom: 0px solid; border-top: 0px solid" ) %>%
            row_spec( c(1,3,5,7,9,11,13,15,17), bold = TRUE, font_size = 12, extra_css = "border-top: 0px solid; border-bottom: 1px solid; border-bottom-color: #DAE8FC" )
            
        return(swtool_ver_table)
    } 
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

#---/ Appendix : References and Additional Databases /---------------------------------------------#
#' @description WES Reference and Additional Databases
#' @export 


