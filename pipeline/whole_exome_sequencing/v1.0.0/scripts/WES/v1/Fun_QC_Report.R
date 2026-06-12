
### GCX WES STANDARD REPORT : NGS-QC FUNCTIONS ###

#---| INDEXING |---------------------------------------------------------------#
#' @param X vector to do indexing
#' @export 
    indexing <- function(X) { sapply(unique(X), function(y) list(y)) }
#------------------------------------------------------------------------------#

#---| ORDER INFORMATION TABLE GENERATION |-------------------------------------#
#' @param CLIENT_NAME client name (in English)
#' @param ORDER_ID    order ID 
#' @param SAMPLE_INFO sample information ( from DB )
#' @export 
    make.ORDER_INFO_TABLE <- function( CLIENT_NAME, ORDER_ID, SAMPLE_INFO )
    {
        orderInfoTable <- data.frame(rbind(
            c("Client Insititute"     , sprintf("   %s",CLIENT_NAME)                                      ),
            c("Total Samples"         , paste0("   ",nrow(SAMPLE_INFO), " samples")                       ),
            c("Date of Sample Receipt", paste0("   ",as.character(as_date(unique(SAMPLE_INFO$date_sample_in)))) ),
            c("Date of Sample QC"     , paste0("   ",as.character(as_date(unique(SAMPLE_INFO$date_sample_qc)))) ),
            c("Date of NGS Report"    , paste0("   ",as.character(today()))                               ),
            c("Order ID"              , paste0("   ",ORDER_ID)                                            ),
            c("Analysis Institute"    , "   GENCURIX"                                                     )
        )) %>% 
        knitr::kable( align=c("c","l"), row.names=F, col.names=NULL, escape=FALSE ) %>%
        kable_styling(full_width=T, bootstrap_options = "basic", font_size = 14, position='left', html_font = 'Noto Sans KR') %>%
        column_spec(1, width="5cm", background = "#D0E4A4") 
        return(orderInfoTable)
    }
#------------------------------------------------------------------------------#

#---| SAMPLE INFORMATION TABLE GENERATION |------------------------------------#
#' @param SAMPLE_INFO sample information ( from DB )
#' @export 
    make.SAMPLE_INFO_TABLE <- function( SAMPLE_INFO )
    {
        sampleInfoTable <- SAMPLE_INFO %>% 
            dplyr::select(c("sample_label","sample_name","sample_info","sample_origin")) %>%
            dplyr::rename('Sample Label'=1, 'Sample ID'=2, 'Tissue Info'=3, 'Sample Type'=4) %>%
            knitr::kable( align=c("l","l","c","c"), row.names=F, escape=FALSE ) %>%
            kable_classic("striped", full_width=T,font_size = 12, position='center', html_font = 'Roboto Condensed') %>%
            row_spec(0, align='center', background = "#D0E4A4")
        return(sampleInfoTable)
    }
#------------------------------------------------------------------------------#

#---| NGS QC RESULT TABLE GENERATION |-----------------------------------------#
#' @param QC_RES           'qc_report' table data ( from DB )
#' @param QC_GROUP         grouping sample's group IDs 
#' @param TAG_SAMPLE_NAMES samples to highlighten in red color. default = NULL
#' @export 
    make.QC_SUMMARY_TABLE <- function( QC_RES, QC_GROUP, TAG_SAMPLE_NAMES=NULL )
    {
        qcResTable = QC_RES %>% 
            dplyr::select(c(11,13,2:10)) %>% 
            mutate( rds_pf = formatC(round(rds_pf), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
            mutate_at(c(4:11), round,1) %>% 
            dplyr::rename(
                'Sample ID'=sample_name, 
                'Sample Origin'=sample_origin, 
                'Total Reads'=rds_pf, 
                'Q30 (%)'=pct_q30_pf, 
                'GC (%)'=pct_gc_content,
                'Duplicates (%)'=pct_rds_dup, 
                'On Target  Bases (%)'=pct_base_on_target, 
                'Target Depth (X)'=depth_target_mean, 
                'Target Coverage  30X'=pct_base_target_30x, 
                'Target Coverage  50X'=pct_base_target_50x, 
                'Target Coverage 100X'=pct_base_target_100x
            ) %>%
            knitr::kable( align=c("l","c","r",rep("c", 8)), row.names=F, escape=FALSE ) %>%
            kable_classic("striped", full_width=T, font_size = 11, position='center', html_font = 'Roboto Condensed') %>%
            row_spec(0, align='c', bold=T, background = "#D0E4A4") 

        qcGroupID <- unique(QC_GROUP$sample_group)
        for( G in 1:length(qcGroupID) ){
            groupSampleRows = which(QC_RES$sample_name %in% (QC_GROUP %>% filter( sample_group == qcGroupID[G] ) %>% .$sample_name) )
            qcResTable      = qcResTable %>% group_rows(group_label=qcGroupID[G], min(groupSampleRows), max(groupSampleRows)) 
        }

        qcResTable <- qcResTable %>%
            column_spec(c(1),        width='3.5cm') %>% 
            column_spec(c(6,7,9:11), width='1.8cm') %>%
            column_spec(c(4:5),      width='1.3cm') %>%
            column_spec(c(8),        width='2cm'  )
        
        if( !is.null(TAG_SAMPLE_NAMES) )
        {
            qcResTable <- qcResTable %>%
                column_spec(1, color=sapply(QC_RES$sample_name, function(y) ifelse( y %in% TAG_SAMPLE_NAMES, "#e64522", "#000000"))) 
        }

        return(qcResTable)
    }
#------------------------------------------------------------------------------#

#---| DRAW GC-CONTENTS BAR PLOT |----------------------------------------------#
#' @param QC_RES       'qc_report' table data ( from DB )
#' @param SAMPLE_INFO  sample information ( from DB )
#' @param SAMPLE_ORDER sample orders
#' @export 
    draw.GC_CONTENTS_PLOT <- function( QC_RES, SAMPLE_INFO, SAMPLE_ORDER )
    {
        gc.plot <- QC_RES %>% dplyr::select(c('seq_id','pct_gc_content')) %>% 
            left_join(. , SAMPLE_INFO %>% dplyr::select(c('seq_id','sample_name')), by='seq_id') %>%
            mutate(sample_name = factor(sample_name, levels=SAMPLE_ORDER[length(SAMPLE_ORDER):1])) %>%
            dplyr::rename(Sample_ID=sample_name,'pct.GC'=pct_gc_content) %>% 
            ggplot(aes( x = pct.GC, y=Sample_ID)) + 
                geom_bar(position="stack", stat="identity", width =0.7, fill='#b3c935') +
                labs( y ="", x = "GC (%)") + 
                coord_cartesian(xlim = c(0,100)) + 
                theme(
                    axis.text.x  = element_text(colour="black", size=9, hjust=0.5),
                    axis.text.y  = element_text(family='Roboto Condensed', colour="black", size=8, hjust=1),
                    axis.title.x = element_text(face='bold', colour="black",size=10 , vjust=1),
                    axis.title.y = element_text(face='bold', colour="black",size=10 , vjust=1),
                    axis.line.x  = element_line(colour = "#6f6d6d"), 
                    axis.line.y  = element_line(colour = "#6f6d6d"),
                    #plot.margin  = unit(c(1, 1, 1, 1), "cm"),
                    #plot.title   = element_text(family="Roboto Condensed", colour="black",size=FONT_SIZE, hjust=0.5),
                    #legend.position  = "none", 
                    panel.background = element_rect(fill='white'), 
                    panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted')
                )+
                geom_text(aes(label = paste0(round(pct.GC,1),"%")), position = position_stack(vjust = .5),size=3, color=c('white')) +
                ggtitle("GC Contents")
        return(gc.plot)
    }
#------------------------------------------------------------------------------#

#---| DRAW DUPLICATES BAR PLOT |-----------------------------------------------#
#' @param QC_RES       'qc_report' table data ( from DB )
#' @param SAMPLE_INFO  sample information ( from DB )
#' @param SAMPLE_ORDER sample orders
#' @export 
    draw.DUPLICATES_PLOT <- function( QC_RES, SAMPLE_INFO, SAMPLE_ORDER )
    {
        dup_qc <- QC_RES %>% dplyr::select(c('seq_id','pct_rds_dup','rds_pf_align')) %>%
            mutate( Uniq.Reads=rds_pf_align*(1-(pct_rds_dup/100)), Dup.Reads=rds_pf_align*(pct_rds_dup/100) , pct_align="" ) %>%
            dplyr::rename( "pct.Dup.Reads"=pct_rds_dup)

        dup_qc2 <- rbind(
            dup_qc %>% 
                dplyr::select(c("seq_id","Uniq.Reads","pct_align")) %>% 
                dplyr::rename( reads=2, pct=3) %>% mutate(total=0, Read_Group='Uniq.Reads'),
            dup_qc %>% 
                dplyr::select(c("seq_id","Dup.Reads","pct.Dup.Reads","rds_pf_align")) %>% 
                dplyr::rename( reads=2, pct=3, total=4) %>% 
                mutate(Read_Group='Dup.Reads', pct=paste0(round(pct,1), "%")) 
        ) %>% 
            left_join(. , SAMPLE_INFO %>% dplyr::select(c('seq_id','sample_name')), by='seq_id') %>% 
            mutate(sample_name = factor(sample_name, levels=SAMPLE_ORDER[length(SAMPLE_ORDER):1])) %>%
            dplyr::rename( 'Sample_ID'=sample_name )

        dup.plot <- ggplot(dup_qc2, aes(x = reads, y = Sample_ID, fill=Read_Group)) + 
            geom_bar(position="stack", stat="identity", width =0.7) +
            labs( x="Reads", y="") +
            theme(
                axis.text.x  = element_text(colour="black", size=9, hjust=0.5),
                axis.text.y  = element_text(family='Roboto Condensed',colour="black", size=8, hjust=1),
                axis.title.x = element_text(face='bold', colour="black",size=10 , vjust=1),
                axis.title.y = element_text(face='bold', colour="black",size=10 , vjust=1),
                axis.line.x  = element_line(colour = "#6f6d6d"), 
                axis.line.y  = element_line(colour = "#6f6d6d"),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted')
            ) +
            geom_text(aes(x=total, y=Sample_ID, label = pct),size=3, color=c('black'), hjust = -0.1) + 
            scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
            scale_fill_manual(values=c('Dup.Reads'='#c63664','Uniq.Reads'='#b3c935')) +
            ggtitle("Duplicates") +
            coord_cartesian(xlim = c(0, max(dup_qc2$total)*1.1))
        return(dup.plot)
    }
#------------------------------------------------------------------------------#

#---| DRAW ALIEGNMENT STATS BAR PLOT |-----------------------------------------#
#' @param QC_RES       'qc_report' table data ( from DB )
#' @param SAMPLE_INFO  sample information ( from DB )
#' @param SAMPLE_ORDER sample orders
#' @export 
    draw.ALIGNMENT_PLOT <- function( QC_RES, SAMPLE_INFO, SAMPLE_ORDER )
    {
        align_qc <- QC_RES %>% 
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
            mutate(Read_Group=factor(Read_Group, levels=c("Unaligned","Aligned PF","Aligned PF & HQ"))) %>%
            left_join(. , SAMPLE_INFO %>% dplyr::select(c('seq_id','sample_name')), by='seq_id') %>% 
            mutate(sample_name = factor(sample_name, levels=SAMPLE_ORDER[length(SAMPLE_ORDER):1])) %>%
            dplyr::rename( 'Sample_ID'=sample_name )
        
        align.plot <- ggplot(align_qc2, aes(x = reads, y = Sample_ID, fill=Read_Group)) + 
            geom_bar(position="stack", stat="identity", width =0.7) +
            labs( x="Reads", y="") + 
            theme(
                axis.text.x  = element_text(colour="black", size=9, hjust=0.5),
                axis.text.y  = element_text(family='Roboto Condensed',colour="black", size=8, hjust=1),
                axis.title.x = element_text(face='bold', colour="black",size=10 , vjust=1),
                axis.title.y = element_text(face='bold', colour="black",size=10 , vjust=1),
                axis.line.x  = element_line(colour = "#6f6d6d"), 
                axis.line.y  = element_line(colour = "#6f6d6d"),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted')
            ) +
            geom_text(aes(x=txt_x, y=Sample_ID, label = pct),size=3, color=c('black'), hjust = 0.25) + 
            scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
            scale_fill_manual(values=c("Unaligned"='#c63664', "Aligned PF"='#375523',"Aligned PF & HQ"='#b3c935' )) +
            coord_cartesian(xlim = c(0, max(align_qc2$reads)*1.2)) +
            ggtitle("Alignments")
        return(align.plot)
    }
#------------------------------------------------------------------------------#

#---| DRAW COVERAGE BAR PLOT |-------------------------------------------------#
#' @param QC_RES       'qc_report' table data ( from DB )
#' @param SAMPLE_INFO  sample information ( from DB )
#' @param SAMPLE_ORDER sample orders
#' @param INDEX        index list format of sample 
#' @param X_LIMIT      X-axis maximum limit ( default = NULL )
#' @export 
    draw.COVERAGE_PLOT <- function( QC_RES, SAMPLE_INFO, SAMPLE_ORDER, INDEX, X_LIMIT=NULL )
    {
        mosd_cov <- ldply(lapply(INDEX, function(y) {
            df = QC_RES %>% filter( seq_id == y)
            data.frame( 
                depth = as.numeric(unlist(strsplit(df$mosdepth_cov_x_depth, ";" ))),
                cov   = as.numeric(unlist(strsplit(df$mosdepth_cov_y_cov, ";" )))
            ) %>% 
            mutate( seq_id = y )
        })) %>% 
            filter( depth != 0 ) %>% 
            left_join(. , SAMPLE_INFO %>% dplyr::select(c('seq_id','sample_name')), by='seq_id') %>% 
            mutate(sample_name = factor(sample_name, levels=SAMPLE_ORDER[length(SAMPLE_ORDER):1])) %>%
            dplyr::rename( 'Sample_ID'=sample_name )

        cov.plot <- ggplot(mosd_cov, aes(x=depth, y=cov, group=Sample_ID, color=Sample_ID)) + 
            geom_line( linewidth = 0.5 ) +
            labs( x = "Coverage Depth", y = "Target (%)") + 
            theme(
                axis.text.x  = element_text(colour="black", size=12, hjust=0.5),
                axis.text.y  = element_text(colour="black", size=12, hjust=0.5),
                axis.title.x = element_text(face='bold', colour="black",size=14 , vjust=1),
                axis.title.y = element_text(face='bold', colour="black",size=14 , vjust=1),
                axis.line.x  = element_line(colour = "#6f6d6d"), 
                axis.line.y  = element_line(colour = "#6f6d6d"),
                plot.margin  = unit(c(1, 1, 1, 1), "cm"),
                legend.position  = "bottom", 
                legend.text      = element_text(family='Roboto Condensed',size=9),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted')
            )+
            scale_color_discrete(name="") +
            ggtitle("Target Coverage Depth")
        
        if( !is.null(X_LIMIT) ){ cov.plot <- cov.plot +  coord_cartesian(xlim=c(0, X_LIMIT)) }

        return(cov.plot)
    }
#------------------------------------------------------------------------------#

#---| SEQ-INFO COMMON TABLE |--------------------------------------------------#
#' @export 
    SEQUENCING_INFO_TABLE <- function(PANEL="WES", READ_LENGTH=151)
    {
        require(lubridate)
        if( PANEL == "WES" )
        {
            seqInfo <- data.frame(rbind(
                c("Application"       , "   Whole Exome Sequencing (WES)"   ),
                c("Material"          , "   DNA (gDNA)"                     ),
                c("Library"           , "   Twist Exome 2.0"                ),
                c("Type of Read"      , "   Paired-end"                     ),
                c("Read Length"       , paste0("   ",READ_LENGTH)                 ),
                c("Platform"          , "   Illumina Novaseq"               ),
                c("Date of NGS Report", paste0("   ",as.character(today())) ),
                c("Analysis Institute", "   GENCURIX"                       )
            )) 
        }
        if( PANEL == "WES2" )
        {
            seqInfo <- data.frame(rbind(
                c("Application"       , "   Whole Exome Sequencing (WES)"   ),
                c("Material"          , "   DNA (gDNA)"                     ),
                c("Library"           , "   SuerSelect V6"                ),
                c("Type of Read"      , "   Paired-end"                     ),
                c("Read Length"       , paste0("   ",READ_LENGTH)                 ),
                c("Platform"          , "   Illumina Novaseq"               ),
                c("Date of NGS Report", paste0("   ",as.character(today())) ),
                c("Analysis Institute", "   GENCURIX"                       )
            )) 
        }
        if( PANEL == "TSO" )
        {
            seqInfo <- data.frame(rbind(
                c("Application"       , "   TruSight Oncology 500 (TSO500)" ),
                c("Material"          , "   DNA / RNA"                      ),
                c("Type of Read"      , "   Paired-end"                     ),
                c("Platform"          , "   Illumina NextSeq 550Dx"         ),
                c("Date of NGS Report", paste0("   ",as.character(today())) ),
                c("Analysis Institute", "   GENCURIX"                       )
            ))    
        }

        seqInfoTable <- seqInfo %>% 
            knitr::kable( align=c("c","l"), row.names=F, col.names=NULL, escape=FALSE ) %>%
            kable_styling(full_width=T, bootstrap_options = "basic", font_size = 14, position='left', html_font = 'Noto Sans KR') %>%
            column_spec(1, width="5cm", background = "#D0E4A4") 

        return(seqInfoTable)
    }
#------------------------------------------------------------------------------#

#---| FASTQ-FILES-INFO |-------------------------------------------------------#
#' @export
    FASTQ_FILES <- function()
    {
        fastqFiles <- as.data.frame(rbind(
            c("   [ Sample ID ]_R1.fastq.gz", "   Raw-FASTQ Read 1 File"      ),
            c("   [ Sample ID ]_R2.fastq.gz", "   Raw-FASTQ Read 2 File"      ),
            c("   FASTQ.md5sum"             , "   MD5 checksum of FASTQ Fles")
        )) %>% 
            knitr::kable( align=c("l","l"), row.names=F, col.names=c("FileName","Description"), escape=FALSE ) %>%
            kable_classic("striped", full_width=T,font_size = 12, position='left', html_font = 'Roboto Condensed') %>%
            row_spec(0, background = "#D0E4A4") %>%
            column_spec(1, width = '12cm') 
        return(fastqFiles)
    }
#------------------------------------------------------------------------------#

#---| BAM-FILES-INFO |---------------------------------------------------------#
#' @export
    BAM_FILES <- function()
    {
        bamFiles <- as.data.frame(rbind(
            c("   [ Sample ID ].bam"     , "   Aligned BAM File"       ),
            c("   [ Sample ID ].bam.bai" , "   BAM Index File"           ),
            c("   BAM.md5sum"            , "   MD5 checksum of BAM Files")
        )) %>% 
            knitr::kable( align=c("l","l"), row.names=F, col.names=c("FileName","Description"), escape=FALSE ) %>%
            kable_classic("striped", full_width=T,font_size = 12, position='left', html_font = 'Roboto Condensed') %>%
            row_spec(0, background = "#D0E4A4") %>%
            column_spec(1, width = '12cm') 
        return(bamFiles)
    }
#------------------------------------------------------------------------------#

#---| VCF-FILES-INFO |---------------------------------------------------------#
#' @export
    VCF_FILES <- function(TONLY=TRUE, TS_N=TRUE, ORG_N=TRUE, GERMLINE=FALSE)
    {
        vcf_types = c()
        if( TONLY ){ vcf_types = c(vcf_types, 1) }
        if( TS_N  ){ vcf_types = c(vcf_types, 2) }
        if( ORG_N ){ vcf_types = c(vcf_types, 3) }
        if( GERMLINE ){ vcf_types = c(vcf_types, 4) }
        #----------------------------------------------------------------------#
        vcfFiles <- as.data.frame(rbind(
            c("   [ Sample ID ].tumor.only.vcf",              "   Annotated VCF file. Variant Calling : Tumor-only"),
            c("   [ Sample ID ].matched.normal.TISSUE.vcf",   "   Annotated VCF file. Variant Calling : Matched Normal - TISSUE"),
            c("   [ Sample ID ].matched.normal.ORGANOID.vcf", "   Annotated VCF file. Variant Calling : Matched Normal - ORGANOID"),
            c("   [ Sample ID ].germline.vcf",                "   VCF file. Variant Calling : DeepVariant Germline Call"),
            c("   VCF.md5sum" ,                               "   MD5 checksum of VCF Fles"             )
        )) %>% slice(c(vcf_types,5)) %>%
            knitr::kable( align=c("l","l"), row.names=F, col.names=c("FileName","Description"), escape=FALSE ) %>%
            kable_classic("striped", full_width=T,font_size = 12, position='left', html_font = 'Roboto Condensed') %>%
            row_spec(0, background = "#D0E4A4") %>%
            column_spec(1, width = '12cm') 
        return(vcfFiles)
    }
#------------------------------------------------------------------------------#

#---| MAF-FILES-INFO |---------------------------------------------------------#
#' @export
    MAF_FILES <- function( TONLY=TRUE, TS_N=TRUE, ORG_N=TRUE, GERMLINE=FALSE )
    {
        maf_types = c()
        if( TONLY ){ maf_types = c(maf_types, 1, 2) }
        if( TS_N  ){ maf_types = c(maf_types, 3, 4) }
        if( ORG_N ){ maf_types = c(maf_types, 5, 6) }
        if( GERMLINE ){ maf_types = c(maf_types, 7) }
        #----------------------------------------------------------------------#

        mafFiles <- as.data.frame(rbind(
            c("   [ Sample ID ].tumor.only.maf",                                    "   MAF Converted from Tumor-only VCF File"               ),
            c("   [ Sample ID ].tumor.somatic.variants.only.maf",                   "   Filtered Variants Only of Tumor-only MAF"             ),
            c("   [ Sample ID ].matched.normal.TISSUE.maf",                         "   MAF Converted from Matched Normal TISSUE VCF"         ),
            c("   [ Sample ID ].matched.normal.TISSUE.somatic.variants.only.maf",   "   Filtered Variants Only of Matched Normal TISSUE MAF"  ),
            c("   [ Sample ID ].matched.normal.ORGANOID.maf",                       "   MAF Converted From Matched Normal ORGANOID VCF"       ),
            c("   [ Sample ID ].matched.normal.ORGANOID.somatic.variants.only.maf", "   Filtered Variants Only of Matched Normal ORGANOID MAF"),
            c("   [ Sample ID ].germline.maf",                                      "   MAF Converted From DeepVariant Germline Call VCF"     ),
            c("   MAF.md5sum",                                                      "   MD5 checksum of MAF Files"                            )
        )) %>% slice(c(maf_types,8)) %>%
            knitr::kable( align=c("l","l"), row.names=F, col.names=c("FileName","Description"), escape=FALSE ) %>%
            kable_classic("striped", full_width=T,font_size = 12, position='left', html_font = 'Roboto Condensed') %>%
            row_spec(0, background = "#D0E4A4") %>%
            column_spec(1, width = '12cm') 
        return(mafFiles)
    }
#------------------------------------------------------------------------------#

#---| REPORT-FILES-INFO |------------------------------------------------------#
#' @export
    REPORT_FILES <- function( STANDARD=TRUE, ADVANCED=TRUE, ONCOPLOT_PNG=TRUE, MAF_XLSX=TRUE, NONHUMAN=FALSE )
    {
        report_contents = c()
        if( STANDARD     ){ report_contents = c(report_contents, 1, 2) }
        if( ADVANCED     ){ report_contents = c(report_contents, 3, 4) }
        if( ONCOPLOT_PNG ){ report_contents = c(report_contents, 5   ) }
        if( ONCOPLOT_PNG ){ report_contents = c(report_contents, 6   ) }
        if( NONHUMAN     ){ report_contents = c(report_contents, 7   ) }
        report_contents = c(report_contents, 8)
        #----------------------------------------------------------------------#
        reportFiles <- as.data.frame(rbind(
            c( "   [ Order ID ].Standard_Analysis_Report.pdf"  , "   NGS QC and Statistics Report. PDF format" ),
            c( "   [ Order ID ].Standard_Analysis_Report.html" , "   NGS QC and Statistics Report. HTML format"),
            c( "   [ Order ID ].Advanced_Analysis_Report.pdf"  , "   WES Analysis Resport. PDF format" ),
            c( "   [ Order ID ].Advanced_Analysis_Report.html" , "   WES Analysis Resport. HTML format"),
            c( "   [ Order ID ].[Tissue].[Geneset].[VariantType].cumulative.samples.oncoplots.png" , "   Cumlative Samples Oncoplots"),
            c( "   [ Order ID ].Preset.Genes.Variants.[ Variant Call ].xlsx" , "   Somatic Variants List. MS Excel format"),
            c( "   [ Order ID ].NonHumanSamples.pdf"           , "   Non-Human Sample Check Result. PDF format"),
            c( "   [ Order ID ].Data_List.pdf"                 , "   Data List Description. PDF format")
        )) %>% slice(report_contents) %>%
            knitr::kable( align=c("l","l"), row.names=F, col.names=c("FileName","Description"), escape=FALSE ) %>%
            kable_classic("striped", full_width=T,font_size = 12, position='left', html_font = 'Roboto Condensed') %>%
            row_spec(0, background = "#D0E4A4") %>%
            column_spec(1, width = '12cm') 
        return(reportFiles)
    }
#------------------------------------------------------------------------------#

#---| APPENDIX TERM KOREAN |---------------------------------------------------#
#' @export 
    APPENDIX_TERM_KOR <- function()
    {
        termList <- c('Total Reads','Q30 (%)','GC (%)','Duplication (%)',
                      'On Target (%)', 'Target Depth (X)','Target Coverage n-X' )
        defList <- c(
        '연기서열분석(시퀀싱.sequencing)라이브러리에 포함된 DNA단편에서생성한 분석량에 대한 염기쌍 정보, 염기서열 분석으로 나온 전체 시퀀스',
        '각 염기(base)별 신뢰성을 표시한 품질 지표로, 해당 염기서열 결과가 0.1%의 오류일 확률',
        'Target영역의 전체 염기중, 구아닌(G)와 사이토신(C)의 비율. GC %가 높을수록 DNA 밀도가 높으며 변성이 잘 이루어지지 않음',
        '중복된 Read의 비율',
        'Target영역에 mapping된 base의 비율',
        '시퀀싱하였을 때 출력된 데이터에서 Target영역에 대한 각 염기 위치별 평균 Read수를 의미',
        'Target영역에 n개 이상의 Read가 있는 비율'
        )

        appendix_term_kor <- data.frame( TERMS=termList, DEFINITION=defList ) %>%
            kable( align=c('c','l'), row.names=F, escape=FALSE ) %>%
            kable_styling( full_width=T, bootstrap_options='condensed', font_size = 12 ) %>%
            column_spec(c(1), width='3cm') %>%
            row_spec(c(0), background = "#D0E4A4")
        return(appendix_term_kor)
    }
#------------------------------------------------------------------------------#

#---| APPENDIX TERM ENGLISH |---------------------------------------------------#
#' @export 
    APPENDIX_TERM_ENG <- function()
    {
        termList <- c('Total Reads','Q30 (%)','GC (%)','Duplication (%)',
                      'On Target (%)', 'Target Depth (X)','Target Coverage n-X' )
        defList <- c(
        "Total Bases of DNA fragments (reads) generated from sequencing library. The entire sequences produced from sequencing.",
        "Sequenced base quality indicator. The probability that the sequence result is an error of 0.1%.",
        "Ratio of Guanine(G) and Cytosine(C) bases in the target DNA region. The higher the GC%, the higher the DNA density and the less likely it is to be denatured.",
        "Rate of duplicated reads",
        "Proportion of base mapped to target region",
        "The average number of reads for each nucleotide position in the target region",
        "Ratio of n or more reads in target region"
        )

        appendix_term_eng <- data.frame( TERMS=termList, DEFINITION=defList ) %>%
            kable( align=c('c','l'), row.names=F, escape=FALSE ) %>%
            kable_styling( full_width=T, bootstrap_options='condensed', font_size = 12 ) %>%
            column_spec(c(1), width='3cm') %>%
            row_spec(c(0), background = "#D0E4A4")
        return(appendix_term_eng)
    }
#------------------------------------------------------------------------------#

#---| APPENDIX REFERENCES |----------------------------------------------------#
#' @export
    APPENDIX_REFERENCE <- function()
    {
        appendix_ref <- as.data.frame(rbind(
            c("Reference Genome Version "   , "  hg19 (GRCh37)"                                                ),
            c("Reference Genome File "      , "  human_g1k_v37_decoy"                                          ),
            c(""                            , "  "                                                             ),
            c("Known SNP Databases "        , "  dbSNP(v138), 1000G(phase1)"                                   ),
            c("Known Indel Databases  "     , "  hg19 known indels, 1000G(phase1), Mills and 1000G goldstandar"),
            c("Germline Variants Databases ", "  gnomAD(v2.1)"                                                 ),
            c("Gene Annotation "            , "  RefSeq(2020-10-26), Ensembl(v110), HGNC(2023-10)"             ),
            c("Functional Annotation"       , "  dbSNP(v154), ClinVar(2023-12-09), SIFT(v2.2.2), POLYPHEN(v5.2.2), COSMIC(v99), Alpha-Missense")
        )) %>% kable( align=c('c','l'), row.names=F, col.names=c('Cetegory','Database and Version'), escape=FALSE ) %>%
            kable_classic( full_width=F, font_size = 12, html_font = 'Roboto Condensed', position = 'left' ) %>% 
            column_spec(c(1), width='4.5cm') %>% column_spec(c(2), width='10cm') %>% row_spec(0, background = "#D0E4A4") %>%
            group_rows(group_label="Reference Genome", 1, 2) %>%
            group_rows(group_label="Database References", 4, 8)
        return(appendix_ref)
    }
#------------------------------------------------------------------------------#

#---| APPENDIX REFERENCES |----------------------------------------------------#
#' @export
    APPENDIX_SWTOOLS <- function()
    {
        appendix_swtools <- as.data.frame(rbind(
            c("fastqc"      ,"  v0.12.1"      ),
            c("fastp"       ,"  v0.23.4"      ),
            c("BWA mem"     ,"  v0.7.17"      ),
            c("Picard"      ,"  v3.1.0"       ),
            c("GATK3"       ,"  v3.8"         ),
            c("GATK4"       ,"  v4.4.0.0"     ),
            c("Ensembl-VEP" ,"  release 110.1"),
            c("vcf2maf"     ,"  v1.6.21"      ),
            c("mosdepth"    ,"  v0.3.6"       ),
            c("alfred"      ,"  v0.2.6"       ),
            c("multiqc"     ,"  v1.18"        ),
            c("OptiType"    ,"  v1.3.3"       ),
            c("HLA-LA"      ,"  v1.0.3"       ),
            c("NGSCheckMate","  v1.0.1"       ),
            c("Fastq-Screen","  v0.15.3"      ),
            c("MANTIS"      ,"  v1.0.5"       ),
            c("MSIsensor2"  ,"  v0.1"         ),
            c("R"           ,"  v4.3.1"       )
        )) %>% kable( align=c('c','l'), row.names=F, col.names=c('Swtool Name','Version'), escape=FALSE ) %>%
            kable_classic( full_width=F, font_size = 12, html_font = 'Roboto Condensed', position = 'left' ) %>% 
            column_spec(c(1), width='4.5cm') %>% column_spec(c(2), width='10cm') %>% row_spec(0, background = "#D0E4A4")
        return(appendix_swtools)
    }
#------------------------------------------------------------------------------#

#---| PAGE BREAKS |------------------------------------------------------------#
#' @export 
    PAGE_BREAKS = function()
    {
        cat('<div style="page-break-after: always;"></div>')
        cat('<div><img style="float: right;" src="resources/gcx_report_logo.png" width="80" height="20"></div>')
        cat('<br>')
        cat('<hr>')
    }
#------------------------------------------------------------------------------#

#---| VARIANT FILTERING  CUT-OFF APPENDIX |------------------------------------#
#' @export 
    APPENDIX_VAR_FILTER.GENE_LIST <- function()
    {
        cutoff_GENE <- data.frame(rbind(
            c("Gene Locus Group", "   HGNC Protein-Coding Genes Only"            ),
            c("Whitelist Genes" , "   14 Eseential Genes of Cancer NGS Panel Test")
        )) %>% 
            kable( align = c("c","l"), row.names = FALSE, escape = FALSE, col.names = NULL, caption = "Gene Filtering" ) %>%
            kable_classic(full_width = FALSE, font_size = 12, position = "left", html_font = "Roboto Condensed") %>%
            column_spec(1, bold = TRUE, width = "5cm", background = "#D0E4A4") %>%
            column_spec(2, width = "12.5cm")
        return(cutoff_GENE)
    }
#' @export
    APPENDIX_VAR_FILTER.AF <- function()
    {
        cutoff_AF <- data.frame(rbind(
            c("Overall Depth"         , "  10 >= Depth "                     ),
            c("Altered Allele Depth"  , "   2 >= Altered Allele Depth"       ),
            c("Allele Frequency (AF)" , "  AF >= 2% (SNV), AF >= 5% (INDEL)" )
        )) %>% 
            kable( align = c("c","l"), row.names = FALSE, escape = FALSE, col.names = NULL, caption = "Allele Frequency Filtering") %>%
            kable_classic(full_width = FALSE, font_size = 12, position = "left", html_font = "Roboto Condensed") %>%
            column_spec(1, bold = TRUE, width = "5cm", background = "#D0E4A4") %>%
            column_spec(2, width = "12.5cm")
        return(cutoff_AF)
    }
#' @export
    APPENDIX_VAR_FILTER.POP_AF <- function()
    {
        cutoff_POPAF <- data.frame(rbind(
            c("1000 Genomes Frequency", "   1000 Genomes Phase3 Allele Frequency & 1000G EastAsian Alelle Frequency <= 1%" ),
            c("gnomAD Frequency"      , "   gnomAD Exome Alelle Frequency & gnomAD Exome EastAsian Alelle Frequency <= 1%" )
        )) %>% 
            kable( align = c("c","l"), row.names = FALSE, escape = FALSE, col.names = NULL, caption = "Population Frequency Filtering") %>%
            kable_classic(full_width = FALSE, font_size = 12, position = "left", html_font = "Roboto Condensed") %>%
            column_spec(1, bold = TRUE, width = "5cm", background = "#D0E4A4") %>%
            column_spec(2, width = "12.5cm")
        return(cutoff_POPAF)
    }
#' @export
    ESSENTIAL_GENE_LIST <- function(gcx_preset_genes)
    {
        essentialGeneList <- data.frame(matrix(c(gcx_preset_genes$ngs_essential_genes, ""), ncol=5)) %>% 
            mutate( id = c("Solid Tumor","NGS Panel", "Essential Genes")) %>%
            relocate(c(6,1:5)) %>% 
            kable( align = c(rep("c",6)), row.names = FALSE, escape = FALSE, col.names=NULL, caption="Solid Tumor NGS Panel Test Essential Genes" ) %>% 
            kable_classic(full_width = FALSE, font_size = 12, position = "left", html_font = "Roboto Condensed") %>%
            column_spec(1, bold = TRUE, width = "5cm", background = "#D0E4A4") %>%
            column_spec(2:6, width = "2.5cm")
        return(essentialGeneList)
    }
#--------------------------------------------------------------------------------------------------#   

#---| PRESET GENESET LIST TABLE |------------------------------------------------------------------#
#' @export 
    preset_C11_names <- c(
        "hCRC"="colorectal", "hLC"="lung", "hPC"="prostate", "HCC"="liver", "hGC"="gastric",
        "Other Cancers"="pan_cancers", "All Cancers"="mutation_panel"
    )
    ORGANOID_SCIENCE_PRESET_GENES_LIST <- function(c11_preset_genes)
    {    
        presetGnTb <- data.frame("Preset Name"=preset_C11_names, check.names=F) %>%
            left_join(., ldply(lapply(c11_preset_genes, function(pg) paste(pg, collapse=" ")), .id="Preset Name"), by="Preset Name") %>% 
            mutate( "Cacner Types"=names(preset_C11_names) ) %>%
            dplyr::select(c(1,3,2)) %>% dplyr::rename(Genes=3) %>%
            kable( align=c("r", "c", "l"), row.names = FALSE, escape = FALSE, 
                caption = "ORGANOID SCIENCE Preset Genes" 
            ) %>%
            kable_classic("striped", full_width = T, font_size=12, html_font = 'Roboto Condensed', position = 'center' ) %>%
            row_spec(0, align='center', font_size = 11, bold=T, background = "#D0E4A4") %>%
            column_spec(1, width = "2.5cm") %>% 
            column_spec(2, width = "3.5cm")
        return(presetGnTb)
    }
#--------------------------------------------------------------------------------------------------#


