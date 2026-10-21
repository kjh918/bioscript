

options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("parallel"))
library(Cairo)  
grDevices::X11.options(type='cairo')
options(device='x11')
options(stringsAsFactors=FALSE) 

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR        <- "/data/wes"
    GENOME_ASSEMBLY <- "hg19"
    WES_LIB_KIT     <- "twist.exome.2.0"
    SEQ_ID          <- ""
    THREADS         <- 15
    CNVKIT_RUN_MODE <- "Tonly"
    VCF_TAG         <- "mutect2.NT"
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"),    action="store", default=NA, type="character", help="BASE_DIR"),
        make_option(c("--SEQ_FOLDER"),  action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--SEQ_ID"),      action="store", default=NA, type="character", help="SEQ_ID"),
        make_option(c("--ASSEMBLY"),    action="store", default="hg19", type="character", help="GENOME_ASSEMBLY"),
        make_option(c("--THREADS"),     action="store", default=15, type="double", help="THREADS"),
        make_option(c("--CNVKIT_RUN_MODE"),  action="store", default="Tonly", type="character", help="CNVKIT RUN MODE"),
        make_option(c("--VCF_TAG"),  action="store", default="mutect2.NT", type="character", help="VCF TAG")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR        <- ARGS$BASE_DIR
    SEQ_FOLDER      <- ARGS$SEQ_FOLDER
    SEQ_ID          <- ARGS$SEQ_ID
    GENOME_ASSEMBLY <- ARGS$ASSEMBLY
    CNVKIT_RUN_MODE <- ARGS$CNVKIT_RUN_MODE
    VCF_TAG         <- ARGS$VCF_TAG
    THREADS         <- as.numeric(ARGS$THREADS)
#--------------------------------------------------------------------------------------------------#

#---| FUNCTIONS |----------------------------------------------------------------------------------#
    #' @description convert cnvkit call result to gene-level copy number alteration result format
    #' @param CNV_REF_GENES   gene list table from refFlat file which used in CNVkit annotation
    #' @param CNVKIT_CALL_RES cnvkit cnv call result .cns file
    #' @param PLOIDY          integer formatted ploidy for adjust copy number ( from ASCAT or PureCN ) 
    #' @param THREADS         multi-core n-cpu
    #' @param ADJUST_CN       if adjust CN with input ploidy or not
    #' @export
        GENE_LEVEL_CNA <- function( CNV_REF_GENES, CNVKIT_CALL_RES, PLOIDY, THREADS, ADJUST_CN=TRUE )
        {
            if( ADJUST_CN )
            {
                CNVKIT_CALL_RES$adjust_cn <- round( 2^CNVKIT_CALL_RES$log2 * PLOIDY )
            }else{
                CNVKIT_CALL_RES$adjust_cn <- CNVKIT_CALL_RES$cn
            }

            CallRes <- ldply(mclapply( indexing( 1:nrow(CNVKIT_CALL_RES) ), function(k) 
            {
                data.frame(
                    gene = unique(unlist(strsplit(CNVKIT_CALL_RES[ k, "gene" ], ","))),
                    chr  = CNVKIT_CALL_RES[ k, "chromosome" ],
                    cn   = CNVKIT_CALL_RES[ k, "adjust_cn" ]
                )
            }, mc.cores = THREADS )) %>% filter( gene != "-" )

            GENE_CNA <- merge(
                x     = CallRes, 
                y     = CNV_REF_GENES[,c("gene","symbol","name","location","entrez","ens_geneid","status")], 
                by    = "gene", 
                all.x = T 
            ) %>% filter( !is.na(symbol) ) %>% unique()

            GENE_CNA$absolute_cna <- sapply( GENE_CNA$cn, function(n) 
            {
                if( n == 0         ){ cna1 <- "DELETION"      }
                if( n == 1         ){ cna1 <- "LOSS"          }
                if( n == 2         ){ cna1 <- "NEUTRAL"       }
                if( n  > 2 & n < 8 ){ cna1 <- "GAIN"          }
                if( n >= 8         ){ cna1 <- "AMPLIFICATION" }
                return(cna1)
            })
            
            GENE_CNA$ploidy_cna <- sapply( GENE_CNA$cn, function(n) 
            {
                if( n == 0                     ){ cna2 <- "DELETION"      }
                if( n >= 1 & n < PLOIDY        ){ cna2 <- "LOSS"          }
                if( n == PLOIDY                ){ cna2 <- "NEUTRAL"       }
                if( n  > PLOIDY & n < PLOIDY*4 ){ cna2 <- "GAIN"          }
                if( n >= PLOIDY*4              ){ cna2 <- "AMPLIFICATION" }
                return(cna2)
            })

            GIDX <- indexing( GENE_CNA$gene )

            GENE_CNA_rev <- mclapply( GIDX, function(gn) 
            {
                gcn <- GENE_CNA %>% filter( gene == gn ) 

                if( nrow(gcn) > 1 )
                {   
                    colOrders = colnames(gcn)
                    copies = gcn$cn

                    maxAbsCopy = ifelse( length(copies[ copies >= 2 ]) > 0, max(copies[ copies >= 2 ]), NA )
                    minAbsCopy = ifelse( length(copies[ copies  < 2 ]) > 0, min(copies[ copies  < 2 ]), NA )

                    if( is.na(maxAbsCopy) )
                    { 
                        maxAbsCNA = NA 
                    }else{
                        if( maxAbsCopy == 2 ){ maxAbsCNA = "NEUTRAL" }else{ maxAbsCNA = ifelse( maxAbsCopy >= 8, "AMPLIFICATION", "GAIN" ) }
                    }

                    if( is.na(minAbsCopy) )
                    {
                        minAbsCNA = NA
                    }else{
                        minAbsCNA = ifelse( minAbsCopy == 0 , "DELETION", "LOSS" )
                    }

                    maxPldCopy = ifelse( length(copies[ copies >= PLOIDY ]) > 0, max(copies[ copies >= PLOIDY ]), NA )
                    minPldCopy = ifelse( length(copies[ copies  < PLOIDY ]) > 0, min(copies[ copies  < PLOIDY ]), NA )

                    if( is.na(maxPldCopy) )
                    { 
                        maxPldCNA = NA
                    }else{
                        if( maxPldCopy == PLOIDY)
                        { maxPldCNA = "NEUTRAL" }else{ maxPldCNA = ifelse( maxPldCopy >= PLOIDY*4, "AMPLIFICATION", "GAIN") }
                    }

                    if( is.na(minPldCopy) )
                    {
                        minPldCNA = NA
                    }else{
                        minPldCNA = ifelse( minPldCopy == 0 , "DELETION", "LOSS" ) 
                    }

                    abs_cna_res = c(maxAbsCNA, minAbsCNA)
                    if( length(abs_cna_res[!is.na(abs_cna_res)]) < 2 ){ abs_cna_res = abs_cna_res[!is.na(abs_cna_res)] }

                    ploidy_cna_res = c(maxPldCNA, minPldCNA)
                    if( length(ploidy_cna_res[!is.na(ploidy_cna_res)]) < 2 ){ ploidy_cna_res = ploidy_cna_res[!is.na(ploidy_cna_res)] }

                    gcn2 = data.frame(
                        unique(gcn[,c("gene","chr","symbol","name","location","entrez","ens_geneid","status")]),
                        cn           = c(maxAbsCopy, minAbsCopy),
                        absolute_cna = abs_cna_res,
                        ploidy_cna   = ploidy_cna_res
                    )

                    gcn = gcn2[ which(!is.na(gcn2$cn)), colOrders ]
                }

                return(gcn)

            }, mc.cores = THREADS)

            GENE_CNA_RESULTS = ldply(GENE_CNA_rev)[,-1] 

            return(GENE_CNA_RESULTS)
        }
    #--------------------------------------------------------------------------#
    #' @description make index list
    #' @export
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
    #--------------------------------------------------------------------------#
    #' @description draw CNA class barplots 
    #' @param CNA_STATS CNA class stats table
    #' @export
        DRAW_CNA_CLASS_BARPLOT <- function( CNA_STATS, X_LABEL_ROTATE = TRUE )
        {
            CNA_STATS %>% arrange( desc(RATIO) ) %>%  
                mutate( CNA_CLASS = gsub(",",", ", CNA_CLASS) ) %>%
                mutate( CNA_CLASS = factor(CNA_CLASS, levels = CNA_CLASS)) %>%
                ggplot( aes(x = CNA_CLASS, y = RATIO, fill = CNA_CLASS) ) +
                    geom_bar(width = 0.75, stat = "identity", position = "dodge") +
                    geom_hline( yintercept = 0, color = "#6f6d6d" ) + 
                    coord_cartesian( ylim = c(0,100) ) + 
                    labs( x = "", y = "Genes (%)", fill = "CNA Class" ) +
                    theme(
                        panel.background = element_rect(fill='white'), 
                        panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth = 0.5),
                        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                        plot.title   = element_text(face = "bold", size = 16),
                        axis.line.y  = element_line(colour = "#6f6d6d"),
                        axis.text.y  = element_text(face = "bold", size = 14),
                        axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 5)),
                        axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
                        axis.ticks.x = element_blank(), 
                        legend.key.size = unit(4, 'mm')
                    ) + 
                    ggtitle( unique(CNA_STATS$CNA_GROUP) )
            if( X_LABEL_ROTATE )
            {
                CNA.ClassPlot <- CNA.ClassPlot + theme( axis.text.x  = element_text(face = "bold", size = 10, angle = 50, vjust= 1, hjust = 1) )
            }else{
                CNA.ClassPlot <- CNA.ClassPlot + theme( axis.text.x  = element_text(face = "bold", size = 12) )
            }
            return(CNA.ClassPlot)
        }
    #--------------------------------------------------------------------------#
    #' @description add cumulative chromosome coordinates for plots 
    #' @param CNR              log-ratio of bin data ( data of .cnr file of cnvkit )
    #' @param CHR_LENGTH_TABLE chromosome length data ( from kangsm RDS )
    #' @export   
        addChrCoords <- function( CNR, CHR_LENGTH_TABLE ) 
        {
            CNR$start2 <- apply( CNR[,c("chromosome","start")], 1, function(Z) 
            { as.numeric(Z[2]) + as.numeric(CHR_LENGTH_TABLE[which(CHR_LENGTH_TABLE$chr == Z[1]), "add_X"]) })

            CNR$end2 <- apply( CNR[,c("chromosome","end")],   1, function(Z) 
            { as.numeric(Z[2]) + as.numeric(CHR_LENGTH_TABLE[which(CHR_LENGTH_TABLE$chr == Z[1]), "add_X"]) })

            CNR$mid2 <- apply( CNR[,c("start2","end2")], 1, function(Z) median(Z) )

            return(CNR)
        }
    #--------------------------------------------------------------------------#
    #' @description CNV result concordance
    #' @param CNR cnvkit .cnr      file data
    #' @param CNS cnvkit .call.cns file data
    #' @export   
        MERGE.METRIC.TO.CNR <- function( CNR, METRIC )
        {
            CNR <- CNR %>% filter( gene %nin% c("Antitarget", "-")  )

            CNR$BAF = CNR$LOG2_CALL = CNR$SEGMENT = CNR$LOG2_GAP = CNR$MEAN = CNR$MEDIAN = CNR$MODE = CNR$MSE = CNR$CN = NA
            
            CNR_COLS    <- c("CN","MEAN","MEDIAN","MODE","MSE","SEGMENT","BAF","LOG2_CALL","LOG2_GAP")
            METRIC_COLS <- c("cn","mean","median","mode","mse","segment","baf_call","log2_call","log2_gap") 

            METRIC

            for( k in 1:nrow(METRIC) )
            {   
                CNR[ which( 
                    CNR$start      >= METRIC[k, "start"] & 
                    CNR$end        <= METRIC[k, "end"]   & 
                    CNR$chromosome == METRIC[k, "chromosome"]
                ), 
                CNR_COLS ] <- METRIC[ k, METRIC_COLS] 
            }
            return(CNR)
        }
    #--------------------------------------------------------------------------#
    #' @description draw cnv plot contains logRatio segments and B-allele frequencies (BAF) 
    #' @param SEGMENT          modified cnvkit .cnr data
    #' @param CHR_LENGTH_TABLE chromosome length data ( from kangsm RDS )
    #' @param PARAMS           parameters of input (used in 'GENE_LEVEL_CNA' function) 
    #' @param CHR              selected chromsome ( default = all chromosomes )
    #' @param GENOME_ASSEMBLY  genome assembly version for coordinates
    #' @param CNA_TYPE         coloring method. 'cn' or 'ploidy' ( cn = AbsoluteCN , ploidy = PloidyCN )
    #' @export    
        DRAW_CNV_SEGMENT_PLOT <- function( SEGMENT , CHR_LENGTH_TABLE, PARAMS, CHR = "all", GENOME_ASSEMBLY, CNA_TYPE )
        {
            # CNV COLORS
            CNA.COLORS <- c("DELETION"="#1a5387", "LOSS"="#8ca9c3", "NEUTRAL"="#4c4c4c", "GAIN"="#cf9890", "AMPLIFICATION"="#9f3122")
            # SEGMENT DATA FIX
            SEGMENT <- SEGMENT[which(!is.na(SEGMENT$CN)), ]
            SEGMENT$CNA_ABS <- sapply( SEGMENT$CN, function(n) 
            {
                if( n == 0         ){ cna1 <- "DELETION"      }
                if( n == 1         ){ cna1 <- "LOSS"          }
                if( n == 2         ){ cna1 <- "NEUTRAL"       }
                if( n  > 2 & n < 8 ){ cna1 <- "GAIN"          }
                if( n >= 8         ){ cna1 <- "AMPLIFICATION" }
                return(cna1)
            })
            SEGMENT$CNA_PLOIDY <- sapply( SEGMENT$CN, function(n) 
            {
                if( n == 0                                           ){ cna2 <- "DELETION"      }
                if( n >= 1 & n < PARAMS$PLOIDY_INT                   ){ cna2 <- "LOSS"          }
                if( n == PARAMS$PLOIDY_INT                           ){ cna2 <- "NEUTRAL"       }
                if( n  > PARAMS$PLOIDY_INT & n < PARAMS$PLOIDY_INT*4 ){ cna2 <- "GAIN"          }
                if( n >= PARAMS$PLOIDY_INT*4                         ){ cna2 <- "AMPLIFICATION" }
                return(cna2)
            })
            log2Outlier = 3
            SEGMENT$logRatio = SEGMENT$log2 + SEGMENT$LOG2_GAP       
            SEGMENT[which(SEGMENT$logRatio >  log2Outlier), "logRatio"] <-  log2Outlier
            SEGMENT[which(SEGMENT$logRatio < -log2Outlier), "logRatio"] <- -log2Outlier
            SEGMENT[which(SEGMENT$BAF > 1), "BAF"] <- 1
            SEGMENT$baf1 = SEGMENT$BAF + 1.5
            SEGMENT$baf2 = 1 - SEGMENT$BAF + 1.5
            # PLOT PARAMETERS
            if( CHR == "all"     ){ CHROMOSOME = unique(SEGMENT$chromosome) }else{ CHROMOSOME = CHR }
            if( CHR == "all"     ){ PS    = 0.1 }else{ PS    = 0.5 }
            if( CHR == "all"     ){ chrFS = 5   }else{ chrFS = 12  }
            if( CHR == "all"     ){ ChromosomeTitle = "All Chromosomes" }else{ ChromosomeTitle = CHR }
            if( CNA_TYPE == "cn" ){ PlotType = "AbsoluteCN" }else{ PlotType = "PloidyCN" }

            breakPoint.Y    <- c(-1,-0,1,1.5,2,2.5)
            limit.Y         <- c(-1.5, 3)
            label.Y         <- c("", sprintf("CN = %s",PARAMS$PLOIDY_INT), "",  "0.0","0.5","1.0")
            breakPoint.X    <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$add_X
            label.X         <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$chr
            chrEndPoint.X   <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$cum_chr_length
            chrStartPoint.X <- ifelse( CHR == "all", 0, chrLengthTable %>% filter( chr %in% CHR ) %>% .$add_X %>% head(1))
            
            if( GENOME_ASSEMBLY == "hg19" & CHR != "all" )
            { label.X <- paste0("chr", label.X ) ; ChromosomeTitle <- paste0("chr", ChromosomeTitle ) }

            # CNV PLOT
            cnvPlot <- SEGMENT %>% filter( chromosome %in% CHROMOSOME ) %>% ggplot() +
                coord_cartesian( ylim = limit.Y ) +
                theme(
                    panel.background = element_rect(fill='white'),
                    axis.line.y      = element_line(colour = "#7a7a7a"),
                    axis.ticks.x     = element_blank(),
                    axis.text.x      = element_text(hjust = 0, angle = 30, vjust=1, size=chrFS ),
                    legend.position  = "none"
                ) +
                labs( x = "Chromosome", y = "", color = "CNA" ) +
                geom_hline( yintercept = 0, col = "#7a7a7a", alpha = 0.5 ) 
            
            if( CNA_TYPE == "cn" )
            {
                cnvPlot <- cnvPlot + geom_point( aes( x = mid2, y = logRatio/3 , color = CNA_ABS ), size = PS, alpha = 0.15 )
            }else{
                cnvPlot <- cnvPlot + geom_point( aes( x = mid2, y = logRatio/3 , color = CNA_PLOIDY ), size = PS, alpha = 0.15 )
            }
            
            cnvPlot <- cnvPlot  +
                geom_point( aes( x = mid2, y = logRatio/3 , color = CNA_PLOIDY ), size = PS, alpha = 0.15 ) +
                scale_color_manual(values = CNA.COLORS) +
                geom_hline( yintercept = c(1.5,2,2.5), color = "#b7b7b7", alpha = 0.3, linetype = 'dashed' ) +
                geom_point( aes( x = mid2, y = baf1 ), color = "#a3a3a3", shape = 23, size = PS*1.5, alpha = 0.4 ) +
                geom_point( aes( x = mid2, y = baf2 ), color = "#a3a3a3", shape = 23, size = PS*1.5, alpha = 0.4 ) +
                scale_y_continuous( breaks = breakPoint.Y, labels = label.Y ) +
                scale_x_continuous( breaks = breakPoint.X, labels = label.X ) +
                geom_vline( xintercept = chrStartPoint.X, col = "#7a7a7a", alpha = 0.5, linewidth = 0.2, linetype = 'dashed') +
                geom_vline( xintercept = chrEndPoint.X,   col = "#7a7a7a", alpha = 0.5, linewidth = 0.2, linetype = 'dashed' ) +
                ggtitle( paste0( PlotType, " : ", ChromosomeTitle) )



        }
#--------------------------------------------------------------------------------------------------#

#---| PRE-DEFINED DATA |---------------------------------------------------------------------------#
    # CHR LENGTH
    load("/storage/home/kangsm/myDB/rds/chr_length.Rdata")                     # CHR_LIST,hg19_chr_length,hg38_chr_length"
    #--------------------------------------------------------------------------#
    # CNV REF GENES
    load("/storage/home/kangsm/myDB/rds/cnv_refGenes_list_hg19_hg38.Rdata")    # CNV_REF_GENES
    # CANCER GENES
    load("/storage/home/kangsm/myDB/rds/PresetGenes.Rdata")                    # c11_preset_genes, gcx_preset_genes
    #--------------------------------------------------------------------------#
    if( GENOME_ASSEMBLY == "hg38" ){ chrLengthTable <- hg38_chr_length    }else{ chrLengthTable <- hg19_chr_length }
    if( GENOME_ASSEMBLY == "hg38" ){ cnv_ref_genes  <- CNV_REF_GENES$hg38 }else{ cnv_ref_genes  <- CNV_REF_GENES$hg19 }
#--------------------------------------------------------------------------------------------------#

#---| FILE TAGS |----------------------------------------------------------------------------------#
    if( CNVKIT_RUN_MODE %in% c("nt","NT") )
    { 
        RUN_MODE_TAG <- ifelse( VCF_TAG == "mutect2.ORG.NT", "NT.ORG_N", "NT.TS_N" )
    }else{
        RUN_MODE_TAG <- "Tonly"
    }
#--------------------------------------------------------------------------------------------------#

#---| RE-CALL CNVKIT RESULT |----------------------------------------------------------------------#
    CNVKIT_RES_DIR     <- sprintf("%s/%s/%s/cnv/cnvkit", BASE_DIR, SEQ_FOLDER, SEQ_ID )
    CNVKIT_SUMMARY_DIR <- sprintf("%s/%s.%s.summary",    CNVKIT_RES_DIR, SEQ_ID, RUN_MODE_TAG )
    if( ! dir.exists(CNVKIT_SUMMARY_DIR) ){ system(sprintf("mkdir -p %s",CNVKIT_SUMMARY_DIR )) }
    #--------------------------------------------------------------------------#
    CNVKIT_RES_FILE   <- sprintf("%s/%s.%s.call.cns",    CNVKIT_RES_DIR, SEQ_ID, RUN_MODE_TAG)
    CNVKIT_PARAM_FILE <- sprintf("%s/%s.%s.call.params", CNVKIT_RES_DIR, SEQ_ID, RUN_MODE_TAG)
    #--------------------------------------------------------------------------#
    CNVKIT_RES    <- read.delim( CNVKIT_RES_FILE )
    CNVKIT_PARAMS <- read.delim( CNVKIT_PARAM_FILE )
    #--------------------------------------------------------------------------#
    GeneCNV <- GENE_LEVEL_CNA(
        CNV_REF_GENES   = cnv_ref_genes,
        CNVKIT_CALL_RES = CNVKIT_RES,
        PLOIDY          = CNVKIT_PARAMS$PLOIDY_INT,
        THREADS         = THREADS,
        ADJUST_CN       = TRUE
    )
    #--------------------------------------------------------------------------#
    write.table( GeneCNV, 
        sprintf("%s/%s.%s.CNVkit.GeneLevel.CopyNumberAlterations.txt", CNVKIT_RES_DIR, SEQ_ID, RUN_MODE_TAG),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
#--------------------------------------------------------------------------------------------------#

#---| GENE COPY-NUMBER DISTRIBUTION |--------------------------------------------------------------#     
    Gene_CN_stats <- as.data.frame( table(GeneCNV$cn) ) %>% dplyr::rename( CN = 1, GENES = 2) %>% 
        mutate( GENES_pct = round(GENES/sum(GENES)*100, 1) , ID = SEQ_ID ) 
    #--------------------------------------------------------------------------#
    CN.DistributePlot <- ggplot(Gene_CN_stats) +
        geom_line( aes( x = CN, y = GENES_pct, group = ID), linewidth = 1, color = "#82B366", alpha = 0.5 ) +
        geom_point( aes( x = CN, y = GENES_pct), size = 5, color = "#B85450" ) +
        coord_cartesian( ylim = c(0,100) ) +
        labs( x = "Estimated CN", y = "Genes (%)" ) +
        theme(
            panel.background = element_rect(fill='white'), 
            panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth = 0.5),
            plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            plot.title   = element_text(face = "bold", size = 16),
            axis.line.x  = element_line(colour = "#6f6d6d"),
            axis.line.y  = element_line(colour = "#6f6d6d"),
            axis.text.x  = element_text(face = "bold", size = 12),
            axis.text.y  = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 5)),
            axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
            axis.ticks.x = element_blank()
        ) + 
        ggtitle( "Distribution of Estimated CopyNumbers" )
    #--------------------------------------------------------------------------#
    ggsave(CN.DistributePlot, width=5, height=5, unit='in', dpi=150, type='cairo',
        file = sprintf("%s/%s.%s.Gene.Copynumber.Distribution.png", CNVKIT_SUMMARY_DIR, SEQ_ID, RUN_MODE_TAG )  )
#--------------------------------------------------------------------------------------------------#

#---| COPY-NUMBER-ALTERATIONS STATS PLOT |---------------------------------------------------------#
    cnaStats <- GeneCNV %>% group_by(gene) %>% 
        reframe( CN      = paste(cn, collapse=","), 
                 ABS_CNA = paste(absolute_cna, collapse=","),
                 PLD_CNA = paste(ploidy_cna, collapse=",")
        ) 
    #--------------------------------------------------------------------------#
    cnaStat_CN <- cnaStats %>% group_by( CN ) %>% reframe( CNA = length(gene) ) %>% arrange(CN) %>% 
        mutate( ID = SEQ_ID, RATIO = round(CNA/sum(CNA)*100, 1), CNA_GROUP = "CopyNumber"  ) %>% dplyr::rename( CNA_CLASS = 1) 
    cnaStat_AbsCNA <- cnaStats %>% group_by( ABS_CNA ) %>% reframe( CNA = length(gene) ) %>% 
        mutate( ID = SEQ_ID, RATIO = round(CNA/sum(CNA)*100, 1), CNA_GROUP = "AbsoluteCN CNA" ) %>% dplyr::rename( CNA_CLASS = 1) 
    cnaStat_PldCN  <- cnaStats %>% group_by( PLD_CNA ) %>% reframe( CNA = length(gene) ) %>% 
        mutate( ID = SEQ_ID, RATIO = round(CNA/sum(CNA)*100, 1), CNA_GROUP = "PloidyCN CNA"  ) %>% dplyr::rename( CNA_CLASS = 1) 
    #--------------------------------------------------------------------------#
    CNA.ClassPlotList <- list(
        CN         = DRAW_CNA_CLASS_BARPLOT( CNA_STATS = cnaStat_CN,     X_LABEL_ROTATE = TRUE ),
        AbsoluteCN = DRAW_CNA_CLASS_BARPLOT( CNA_STATS = cnaStat_AbsCNA, X_LABEL_ROTATE = TRUE  ),
        PloidyCN   = DRAW_CNA_CLASS_BARPLOT( CNA_STATS = cnaStat_PldCN,  X_LABEL_ROTATE = TRUE  )
    )
    CNA.ClassPlots <- cowplot::plot_grid( plotlist = CNA.ClassPlotList, ncol = 1 )
    #--------------------------------------------------------------------------#
    ggsave(CNA.ClassPlots, width=5, height=12, unit='in', dpi=150, type='cairo',
        file = sprintf("%s/%s.%s.CNA.Class.Distribution.png", CNVKIT_SUMMARY_DIR, SEQ_ID, RUN_MODE_TAG )  )
#--------------------------------------------------------------------------------------------------#

#---| TSO500 GENES CNA |---------------------------------------------------------------------------#
    GeneCNV_TSO <- GeneCNV %>% filter( gene %in% gcx_preset_genes$tso500 )
    write.table( GeneCNV_TSO, 
        sprintf("%s/%s.%s.CNVkit.GeneLevel.TSO500.CopyNumberAlterations.txt", CNVKIT_RES_DIR, SEQ_ID, RUN_MODE_TAG),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
    #--------------------------------------------------------------------------#
    cnaStats_TSO <- GeneCNV_TSO %>% 
        group_by(gene) %>% 
        reframe( CN      = paste(cn, collapse=","), 
                 ABS_CNA = paste(absolute_cna, collapse=","),
                 PLD_CNA = paste(ploidy_cna, collapse=",")
        ) 
    #--------------------------------------------------------------------------#
    cnaStat_CN_tso <- cnaStats %>% group_by( CN ) %>% reframe( CNA = length(gene) ) %>% arrange(CN) %>% 
        mutate( ID = SEQ_ID, RATIO = round(CNA/sum(CNA)*100, 1), CNA_GROUP = "CopyNumber"  ) %>% dplyr::rename( CNA_CLASS = 1) 
    cnaStat_AbsCNA_tso <- cnaStats %>% group_by( ABS_CNA ) %>% reframe( CNA = length(gene) ) %>% 
        mutate( ID = SEQ_ID, RATIO = round(CNA/sum(CNA)*100, 1), CNA_GROUP = "AbsoluteCN CNA" ) %>% dplyr::rename( CNA_CLASS = 1) 
    cnaStat_PldCN_tso  <- cnaStats %>% group_by( PLD_CNA ) %>% reframe( CNA = length(gene) ) %>% 
        mutate( ID = SEQ_ID, RATIO = round(CNA/sum(CNA)*100, 1), CNA_GROUP = "PloidyCN CNA"  ) %>% dplyr::rename( CNA_CLASS = 1) 
    #--------------------------------------------------------------------------#
    CNA.TSO.ClassPlotList <- list(
        CN         = DRAW_CNA_CLASS_BARPLOT( CNA_STATS = cnaStat_CN_tso,     X_LABEL_ROTATE = TRUE ),
        AbsoluteCN = DRAW_CNA_CLASS_BARPLOT( CNA_STATS = cnaStat_AbsCNA_tso, X_LABEL_ROTATE = TRUE  ),
        PloidyCN   = DRAW_CNA_CLASS_BARPLOT( CNA_STATS = cnaStat_PldCN_tso,  X_LABEL_ROTATE = TRUE  )
    )
    CNA.TSO.ClassPlots <- cowplot::plot_grid( plotlist = CNA.TSO.ClassPlotList, ncol = 1 )
    #--------------------------------------------------------------------------#
    ggsave(CNA.TSO.ClassPlots, width=5, height=12, unit='in', dpi=150, type='cairo',
        file = sprintf("%s/%s.%s.CNA.Class.TSO500.Genes.Distribution.png", CNVKIT_SUMMARY_DIR, SEQ_ID, RUN_MODE_TAG )  )
#--------------------------------------------------------------------------------------------------#

#---| CNA CLASS RESULT FILE SAVE |-----------------------------------------------------------------#
    CNA_CLASS_RES     <- rbind( cnaStat_CN, cnaStat_AbsCNA, cnaStat_PldCN )
    CNA_CLASS_TSO_RES <- rbind( cnaStat_CN_tso, cnaStat_AbsCNA_tso, cnaStat_PldCN_tso )
    #--------------------------------------------------------------------------#
    write.table( as.data.frame(CNA_CLASS_RES), 
        sprintf("%s/%s.%s.CNVkit.Gene.CNA.Class.Distribution.txt", CNVKIT_SUMMARY_DIR, SEQ_ID, RUN_MODE_TAG),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
    write.table( as.data.frame(CNA_CLASS_TSO_RES), 
        sprintf("%s/%s.%s.CNVkit.TSO500.CNA.Class.Distribution.txt", CNVKIT_SUMMARY_DIR, SEQ_ID, RUN_MODE_TAG),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
#--------------------------------------------------------------------------------------------------#

#---| SEGMENT PLOTS |------------------------------------------------------------------------------#
    CNVKIT_CNR_FILE <- sprintf("%s/%s.%s.cnr", CNVKIT_RES_DIR, SEQ_ID, RUN_MODE_TAG)  
    CNR <- read.delim( CNVKIT_CNR_FILE ) %>% filter( gene %nin% c("Antitarget", "-") )
    # SegmentMetrics ----------------------------------------------------------#
    SegMetricFile <- sprintf("%s/%s.%s.segmetrics.txt", CNVKIT_RES_DIR, SEQ_ID, RUN_MODE_TAG)  
    SegMetric     <- read.delim( SegMetricFile ) 
    # Modified CNR ------------------------------------------------------------#
    SM <- SegMetric[,-4] %>% mutate( key = paste(chromosome, start, sep="_") ) %>% 
        left_join(. ,
            CNVKIT_RES %>% mutate( key = paste(chromosome, start, sep="_") ) %>% 
                select(c("key","cn","log2","baf")) %>% dplyr::rename( log2_call = log2, baf_call = baf ),
                by = "key"
        ) %>%
        mutate( segment = paste0("seg", 1:nrow(SM))) %>% 
        select(c(1:4,6:14,16,17,18,19)) %>%
        mutate( log2_gap = log2_call - log2 ) %>%
        mutate( cn = round(2^(log2_call) * CNVKIT_PARAMS$PLOIDY_INT) )
    #--------------------------------------------------------------------------#
    SegRes <- MERGE.METRIC.TO.CNR( CNR = CNR, METRIC = SM )
    SegRes <- addChrCoords( CNR = SegRes, CHR_LENGTH_TABLE = chrLengthTable ) 
    #--------------------------------------------------------------------------#
    CHR_LIST <- c(1:22, "X", "Y")
    if( GENOME_ASSEMBLY == "hg38" ){ CHR_LIST <- paste0("chr", CHR_LIST) }
    #--------------------------------------------------------------------------#
    for( CHR in c("all", CHR_LIST) )
    {
        if( CHR == "all" ){ FW = 12 ; FH = 3, }
        
        AbsCNA.CNVplot <- DRAW_CNV_SEGMENT_PLOT(
            SEGMENT = SegRes,
            CHR_LENGTH_TABLE = chrLengthTable,
            PARAMS = CNVKIT_PARAMS,
            CHR = CHR,
            GENOME_ASSEMBLY = GENOME_ASSEMBLY,
            CNA_TYPE = "cn"
        )
        ggsave( AbsCNA.CNVplot,  width=FW, height=FH, unit='in', dpi=150, type='cairo',
            file = sprintf("%s/%s.%s.AbsoluteCN.%s.CNVplot.png", CNVKIT_SUMMARY_DIR, SEQ_ID, RUN_MODE_TAG, CHR )  
        )
        #----------------------------------------------------------------------#
        PloidyCNA.CNVplot <- DRAW_CNV_SEGMENT_PLOT(
            SEGMENT = SegRes,
            CHR_LENGTH_TABLE = chrLengthTable,
            PARAMS = CNVKIT_PARAMS,
            CHR = CHR,
            GENOME_ASSEMBLY = GENOME_ASSEMBLY,
            CNA_TYPE = "ploidy"
        )
        ggsave( PloidyCNA.CNVplot,  width=FW, height=FH, unit='in', dpi=150, type='cairo',
            file = sprintf("%s/%s.%s.PloidyCN.%s.CNVplot.png", CNVKIT_SUMMARY_DIR, SEQ_ID, RUN_MODE_TAG, CHR )  
        )
    }
#--------------------------------------------------------------------------------------------------#

### HERE AFTER, SOURCE CODES FOR SEQC2 VALIDATION

# #--------------------------------------------------------------------------------------------------#
#     #' @description compare two CNA result and get concordance
#     #' @param GeneLevelCNV1 result of 'GENE_LEVEL_CNA' (input 1)
#     #' @param GeneLevelCNV2 result of 'GENE_LEVEL_CNA' (input 2)
#     #' @param PARAM1        parameters of input 1 (used in 'GENE_LEVEL_CNA' function) 
#     #' @param PARAM2        parameters of input 2 (used in 'GENE_LEVEL_CNA' function) 
#     #' @export
#         CNV_CONCORDANCE <- function( GeneLevelCNV1, GeneLevelCNV2, PARAM1, PARAM2 )
#         {
#             res1 <- GeneLevelCNV1 %>% group_by(gene) %>% 
#                 reframe(
#                     cn         = paste(sort(cn),           collapse=";"), 
#                     abs_cna    = paste(sort(absolute_cna), collapse=";"),
#                     ploidy_cna = paste(sort(ploidy_cna),   collapse=";")
#                 ) %>% as.data.frame()
            
#             res2 <- GeneLevelCNV2 %>% group_by(gene) %>% 
#                 reframe(
#                     cn         = paste(sort(cn),           collapse=";"), 
#                     abs_cna    = paste(sort(absolute_cna), collapse=";"),
#                     ploidy_cna = paste(sort(ploidy_cna),   collapse=";")
#                 ) %>% as.data.frame()

#             GeneSpace <- unique(intersect( GeneLevelCNV1$gene, GeneLevelCNV2$gene ))

#             RES <- data.frame(
#                 gene        = GeneSpace,
#                 cn1         = res1[match(GeneSpace, res1$gene), "cn"],
#                 cn2         = res2[match(GeneSpace, res2$gene), "cn"],
#                 abs_cna1    = res1[match(GeneSpace, res1$gene), "abs_cna"],
#                 abs_cna2    = res2[match(GeneSpace, res2$gene), "abs_cna"],
#                 ploidy_cna1 = res1[match(GeneSpace, res1$gene), "ploidy_cna"],
#                 ploidy_cna2 = res2[match(GeneSpace, res2$gene), "ploidy_cna"]        
#             )

#             CON_RES <- data.frame(
#                 genespace = length(GeneSpace),
#                 concordance_CN = round(sum(apply(RES[,c("cn1","cn2")], 1, function(y) y[1] == y[2]))/length(GeneSpace)*100, 1),
#                 concordance_AbsoluteCNA = round(sum(apply(RES[,c("abs_cna1","abs_cna2")], 1, function(y) y[1] == y[2]))/length(GeneSpace)*100, 1),
#                 concordance_PloidyCNA   = round(sum(apply(RES[,c("ploidy_cna1","ploidy_cna2")], 1, function(y) y[1] == y[2]))/length(GeneSpace)*100, 1),
#                 input1_genes = length(unique(GeneLevelCNV1$gene)),
#                 input2_genes = length(unique(GeneLevelCNV2$gene)),
#                 input1_cn_range = paste(range(GeneLevelCNV1$cn), collapse="-"),
#                 input2_cn_range = paste(range(GeneLevelCNV2$cn), collapse="-"),
#                 input1_ploidy=PARAM1$PLOIDY,
#                 input2_ploidy=PARAM2$PLOIDY,
#                 input1_ploidy_int=PARAM1$PLOIDY_INT,
#                 input2_ploidy_int=PARAM2$PLOIDY_INT,
#                 input1_purity=PARAM1$PURITY,
#                 input2_purity=PARAM2$PURITY
#             )

#             return(CON_RES)
#         }
# #--------------------------------------------------------------------------------------------------#

# #---| VALIDATION |---------------------------------------------------------------------------------#
#     ascat_valid = pcn_NT_valid = pcn_Tonly_valid = data.frame()

#     for( ID in c("FD1_T","FD2_T","FD3_T","NV1_T","NV2_T","NV3_T") )
#     {
#         ansid = unlist(strsplit(ID, ""))
#         ID2 = paste0(ansid[1], ansid[2], "_", ansid[3])
#         message(sprintf("ID2 = %s", ID2))

#         ASCAT_CALL_FILE <- sprintf("/storage/home/kangsm/analysis/cnv_seqc2/hg38/%s/cnv/cnvkit/%s.NT.TS_N.ASCAT.call.cns", ID, ID)
#         ASCAT_PARAMS    <- sprintf("/storage/home/kangsm/analysis/cnv_seqc2/hg38/%s/cnv/cnvkit/%s.NT.TS_N.call.ASCAT.params", ID, ID)

#         PURECN_NT_CALL_FILE <- sprintf("/storage/home/kangsm/analysis/cnv_seqc2/hg38/%s/cnv/cnvkit/%s.NT.TS_N.PureCN.NT.call.cns", ID, ID)
#         PURECN_NT_PARAMS    <- sprintf("/storage/home/kangsm/analysis/cnv_seqc2/hg38/%s/cnv/cnvkit/%s.NT.TS_N.call.PureCN.params", ID, ID)

#         PURECN_TONLY_CALL_FILE <- sprintf("/storage/home/kangsm/analysis/cnv_seqc2/hg38/%s/cnv/cnvkit/%s.NT.TS_N.PureCN.Tonly.call.cns", ID, ID)
#         PURECN_TONLY_PARAMS    <- sprintf("/storage/home/kangsm/analysis/cnv_seqc2/hg38/%s/cnv/cnvkit/%s.NT.TS_N.call.PureCN.Tonly.params", ID, ID)

#         ANSWERS_FILE <- sprintf("/storage/home/kangsm/myDB/benchmarks/seqc2/cnv/cnvkit/WES_%s.call.cns", ID2)
#         ANSWER_PARAMS <- data.frame( PLOIDY=3, PLOIDY_INT=3, PURITY=1 )

#         seqc2.cnv <- GENE_LEVEL_CNA(
#             CNV_REF_GENES   = cnv_ref_genes,
#             CNVKIT_CALL_RES = read.delim( ANSWERS_FILE ),
#             PLOIDY          = 3,
#             THREADS         = 15,
#             ADJUST_CN       = FALSE
#         )
#         # ASCAT paramter
#         ascat.param.cnv <- GENE_LEVEL_CNA(
#             CNV_REF_GENES   = cnv_ref_genes,
#             CNVKIT_CALL_RES = read.delim( ASCAT_CALL_FILE ),
#             PLOIDY          =  read.delim( ASCAT_PARAMS )$PLOIDY_INT,
#             THREADS         = 15,
#             ADJUST_CN       = TRUE
#         )
#         ascat_res = CNV_CONCORDANCE( GeneLevelCNV1 = seqc2.cnv, GeneLevelCNV2 = ascat.param.cnv, PARAM1 = ANSWER_PARAMS, PARAM2 = read.delim( ASCAT_PARAMS ) )
#         ascat_res = data.frame( id = ID, ascat_res )
#         ascat_valid = rbind(ascat_valid, ascat_res)
#         message(sprintf("%s ASCAT done.", ID))
#         #--------------------------------------------------------------------------#
#         # PureCN NT-mode parameter
#         pcn.NT.param.cnv <- GENE_LEVEL_CNA(
#             CNV_REF_GENES   = cnv_ref_genes,
#             CNVKIT_CALL_RES = read.delim( PURECN_NT_CALL_FILE ),
#             PLOIDY          = read.delim( PURECN_NT_PARAMS )$PLOIDY_INT,
#             THREADS         = 15,
#             ADJUST_CN       = TRUE
#         )
#         pcnNT_res = CNV_CONCORDANCE( GeneLevelCNV1 = seqc2.cnv, GeneLevelCNV2 = pcn.NT.param.cnv, PARAM1 = ANSWER_PARAMS, PARAM2 = read.delim( PURECN_NT_PARAMS ) )
#         pcnNT_res = data.frame( id = ID, pcnNT_res )
#         pcn_NT_valid = rbind(pcn_NT_valid, pcnNT_res)
#         message(sprintf("%s PureCN NT done.", ID))
#         #--------------------------------------------------------------------------#
#         # PureCN Tonly-mode parameter
#         pcn.Tonly.param.cnv <- GENE_LEVEL_CNA(
#             CNV_REF_GENES   = cnv_ref_genes,
#             CNVKIT_CALL_RES = read.delim( PURECN_TONLY_CALL_FILE ),
#             PLOIDY          = read.delim( PURECN_TONLY_PARAMS )$PLOIDY_INT,
#             THREADS         = 15,
#             ADJUST_CN       = TRUE
#         )
#         pcnTonly_res = CNV_CONCORDANCE( GeneLevelCNV1 = seqc2.cnv, GeneLevelCNV2 = pcn.Tonly.param.cnv, PARAM1 = ANSWER_PARAMS, PARAM2 = read.delim( PURECN_TONLY_PARAMS ) )
#         pcnTonly_res = data.frame( id = ID, pcnTonly_res )
#         pcn_Tonly_valid = rbind(pcn_Tonly_valid, pcnTonly_res)
#          message(sprintf("%s PureCN Tonly done.", ID))
#         #--------------------------------------------------------------------------#
#         message(sprintf("%s FINISHED.", ID))
#     }
# #--------------------------------------------------------------------------------------------------#
    



