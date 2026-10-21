
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
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR        <- "/data/wes"
    GENOME_ASSEMBLY <- "hg19"
    WES_LIB_KIT     <- "twist.exome.2.0"
    SEQ_ID          <- ""
    THREADS         <- 15
    RUN_MODE        <- "Tonly"
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"),    action="store", default=NA, type="character", help="BASE_DIR"),
        make_option(c("--SEQ_FOLDER"),  action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--SEQ_ID"),      action="store", default=NA, type="character", help="SEQ_ID"),
        make_option(c("--ASSEMBLY"),    action="store", default="hg19", type="character", help="GENOME_ASSEMBLY"),
        make_option(c("--THREADS"),     action="store", default=15, type="double", help="THREADS"),
        make_option(c("--RUN_MODE"),    action="store", default="Tonly", type="character", help="CNVKIT RUN MODE")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR        <- ARGS$BASE_DIR
    SEQ_FOLDER      <- ARGS$SEQ_FOLDER
    SEQ_ID          <- ARGS$SEQ_ID
    GENOME_ASSEMBLY <- ARGS$ASSEMBLY
    RUN_MODE        <- ARGS$RUN_MODE
    THREADS         <- as.numeric(ARGS$THREADS)
#--------------------------------------------------------------------------------------------------#

#---| PRE-DEFINED DATA |---------------------------------------------------------------------------#
    # CHR LENGTH
    load("/storage/home/kangsm/myDB/rds/chr_length.Rdata")                     # CHR_LIST,hg19_chr_length,hg38_chr_length"
    #--------------------------------------------------------------------------#
    # CNV REF GENES
    load("/storage/home/kangsm/myDB/rds/cnv_refGenes_list_hg19_hg38.Rdata")    # CNV_REF_GENES
    # CANCER GENES
    load("/storage/home/kangsm/myDB/rds/PresetGenes.Rdata")                    # c11_preset_genes, gcx_preset_genes
#--------------------------------------------------------------------------------------------------#

#---| FUNCTIONS |----------------------------------------------------------------------------------#
    #' @description add cumulative chromosome coordinates for plots 
    #' @param CNR              log-ratio of bin data ( data of .cnr file of cnvkit )
    #' @param CHR_LENGTH_TABLE chromosome length data ( from kangsm RDS )
    #' @export   
        addChrCoords <- function(CNR, CHR_LENGTH_TABLE) 
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
        CNS2CNR <- function( CNR, CNS, CHR_LENGTH_TABLE )
        {
            CNR %>% filter( gene != "Antitarget"  )

            CNR$cn        <- NA
            CNR$log2.call <- NA
            
            for( k in 1:nrow(CNS) )
            {   
                CNR[ which( 
                    CNR$start      >= CNS[k, "start"] & 
                    CNR$end        <= CNS[k, "end"]   & 
                    CNR$chromosome == CNS[k, "chromosome"]
                ), 
                c("cn","log2.call") ] <- CNS[ k, c("cn","log2") ] 
            }

            CNR_Coords <- addChrCoords( CNR, CHR_LENGTH_TABLE )
            CNR_Coords <- CNR_Coords %>% filter( !is.na(cn) )

            return(CNR_Coords)
        }
    #--------------------------------------------------------------------------#
    #' @description CNVkit result cnv-plot
    #' @param CNVKIT_RESULT cnvkit result ( cns2cnr results )
    #' @param PARAMS        cbvkit call parameter
    #' @param CHR_LENGTH_TABLE chromosome length data ( from kangsm RDS )
    #' @param CHR              plotting target chromosome ( default = 'all' )
    #' @param GENOME_ASSEMBLY  reference genome assembly version hg19 or hg38.
    #' @export
        DRAW_CNVPLOT <- function( CNVKIT_RESULT, CHR_LENGTH_TABLE, PARAMS, CHR="all", GENOME_ASSEMBLY )
        {
            logR.Colors <- c("Del-"="#14426c", "Del"="#8ca9c3", "Neu"="#4c4c4c", "Amp"="#cf9890", "Amp+"="#7f271b")
            
            if( CHR == "all" ){ CHROMOSOME = unique(CNVKIT_RESULT$chromosome) }else{ CHROMOSOME = CHR }
            if( CHR == "all" ){ PS = 0.1 }else{ PS = 0.5 }
            if( CHR == "all" ){ CHR_TITLE = "All Chromosomes" }else{ CHR_TITLE = CHR }

            Y.bp    <- c(-4,-2,0,2,4)
            Y.end   <- c(-5,5)
            X.bp    <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$add_X
            X.lbl   <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$chr
            X.chr   <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$cum_chr_length
            X.start <- ifelse( CHR == "all", 0, chrLengthTable %>% filter( chr %in% CHR ) %>% .$add_X %>% head(1))
            if( GENOME_ASSEMBLY == "hg19" & CHR != "all" )
            {
                X.lbl     <- paste0("chr", X.lbl    )
                CHR_TITLE <- paste0("chr", CHR_TITLE)
            }

            cnvPlot <- CNVKIT_RESULT %>% filter( chromosome %in% CHROMOSOME ) %>% 
                ggplot() +
                    coord_cartesian(ylim = Y.end ) +
                    scale_y_continuous( breaks = Y.bp+(PARAMS$logR_SHIFT), labels = as.character(Y.bp) ) +
                    scale_x_continuous( breaks = X.bp, labels = X.lbl) +
                    theme(
                        panel.background = element_rect(fill='white'),
                        axis.line.y      = element_line(colour = "#7a7a7a"),
                        axis.ticks.x     = element_blank(),
                        axis.text.x      = element_text(hjust = 0, angle = 30, vjust=1, size=5 ),
                        legend.position  = "none"
                    ) + 
                    labs( x = "Chromosome", y = "logRatio", color = "CNA" ) +
                    geom_hline( yintercept = 0+(PARAMS$logR_SHIFT), col = "#7a7a7a", alpha = 0.5 ) +
                    geom_point( aes( x = mid2, y = log2, color = cna ), size = PS, alpha = 0.3 ) +
                    scale_color_manual(values = logR.Colors) +
                    geom_vline( xintercept = X.start, col = "#7a7a7a", alpha = 0.5, linewidth = 0.2, linetype = 'dashed') +
                    geom_vline( xintercept = X.chr, col = "#7a7a7a", alpha = 0.5, linewidth = 0.2, linetype = 'dashed' ) +
                    ggtitle( CHR_TITLE )
                    
            return(cnvPlot)
        }
    #--------------------------------------------------------------------------#
    indexing = function(x){ sapply(unique(x), function(z) list(z)) }
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS & FILES |----------------------------------------------------------------------------#
    if( RUN_MODE %in% c("nt","NT") ){ ANALYSIS_MODE = "NT" }else{ ANALYSIS_MODE = "Tonly" }
    #--------------------------------------------------------------------------#
    CNVKIT_RES_DIR     <- sprintf("%s/%s/%s/cnv/cnvkit", BASE_DIR, SEQ_FOLDER, SEQ_ID)   
    CNVKIT_SUMMARY_DIR <- sprintf("%s/%s.%s.summary", CNVKIT_RES_DIR, SEQ_ID, ANALYSIS_MODE)
    if( ! dir.exists(CNVKIT_SUMMARY_DIR) ){ system(sprintf("mkdir -p %s", CNVKIT_SUMMARY_DIR)) }
    #--------------------------------------------------------------------------#
    CNR_FILE   <- sprintf("%s/%s.%s.cnr",      CNVKIT_RES_DIR, SEQ_ID, ANALYSIS_MODE)
    CNS_FILE   <- sprintf("%s/%s.%s.call.cns", CNVKIT_RES_DIR, SEQ_ID, ANALYSIS_MODE)
    RUN_PARAMS <- sprintf("%s/%s.call.params", CNVKIT_RES_DIR, SEQ_ID)
#--------------------------------------------------------------------------------------------------#

#---| READ CNVKIT RESULTS |------------------------------------------------------------------------#
    resCNR <- read.delim( CNR_FILE   )
    resCNS <- read.delim( CNS_FILE   )
    PARAMS <- read.delim( RUN_PARAMS )
    #--------------------------------------------------------------------------#
    if( GENOME_ASSEMBLY == "hg38" ){ chrLengthTable <- hg38_chr_length }else{ chrLengthTable <- hg19_chr_length }
    #--------------------------------------------------------------------------#
    resCNVKIT <- CNS2CNR( CNR = resCNR, CNS = resCNS, CHR_LENGTH_TABLE = chrLengthTable )
    #--------------------------------------------------------------------------#
    resCNVKIT$cna <- sapply( resCNVKIT$cn, function(n) 
    {
        if( n == 0         ){ cna <- "Del-" }
        if( n == 1         ){ cna <- "Del"  }
        if( n == 2         ){ cna <- "Neu"  }
        if( n  > 2 & n < 8 ){ cna <- "Amp"  }
        if( n >= 8         ){ cna <- "Amp+" }
        return(cna)
    })
    #--------------------------------------------------------------------------#
    resCNVKIT <- resCNVKIT %>% filter( abs(log2) < 5 )
#--------------------------------------------------------------------------------------------------#

#---| CNV-PLOTS |----------------------------------------------------------------------------------#
    chrs4plot <- c(1:22, "X","Y")
    if( GENOME_ASSEMBLY == "hg38" )
    {
        chrs4plot <- c("all", paste0("chr", chrs4plot))
    }else{
        chrs4plot <- c("all", chrs4plot)
    }
    #--------------------------------------------------------------------------#
    for( chr in chrs4plot )
    {
        for( cnaType in c(1,2) )
        {
            outputName <- sprintf("%s/%s.%s.chromosome.%s.cnvplot.png", CNVKIT_SUMMARY_DIR, SEQ_ID, ANALYSIS_MODE, chr)
            outWidth   <- ifelse( chr == "all", 16, 8 )

            CNVPLOT <- DRAW_CNVPLOT(
                CNVKIT_RESULT    = resCNVKIT,
                CHR_LENGTH_TABLE = chrLengthTable,
                PARAMS           = PARAMS,
                GENOME_ASSEMBLY  = GENOME_ASSEMBLY,
                CHR              = chr
            )

            ggsave( CNVPLOT, file = outputName, width = outWidth, height = 3, unit='in', dpi=150, type='cairo')
        }
    }
#--------------------------------------------------------------------------------------------------#

#---| GENE LEVEL CNV-RESULTS |---------------------------------------------------------------------#
    if( GENOME_ASSEMBLY == "hg38" )
    { cnv_genes <- CNV_REF_GENES$hg38 }else{ cnv_genes <- CNV_REF_GENES$hg19 }
    #--------------------------------------------------------------------------#
    geneCN <- ldply(mclapply( indexing( 1:nrow(resCNS) ), function(k) 
    {
        data.frame(
            gene = unique(unlist(strsplit(resCNS[ k, "gene" ], ","))),
            chr  = resCNS[ k, "chromosome" ],
            cn   = resCNS[ k, "cn" ]
        )
    }, mc.cores = THREADS )) %>% filter( gene != "-" )
    #--------------------------------------------------------------------------#
    geneCNA <- merge(geneCN, cnv_genes[,c("gene","symbol","name","location","entrez","ens_geneid","status")], by="gene", all.x=T ) %>%
        filter( !is.na(symbol)) 
    #--------------------------------------------------------------------------#
    geneCNA$cna <- sapply( geneCNA$cn, function(n) 
    {
        if( n == 0         ){ cna <- "Del-" }
        if( n == 1         ){ cna <- "Del"  }
        if( n == 2         ){ cna <- "Neu"  }
        if( n  > 2 & n < 8 ){ cna <- "Amp"  }
        if( n >= 8         ){ cna <- "Amp+" }
        return(cna)
    })
    #--------------------------------------------------------------------------#
    GENE_CNA <- ldply(mclapply( indexing(geneCNA$gene), function(gn)
    {
        gnCN <- geneCNA %>% filter( gene == gn )
        gcns <- unique(gnCN$cn)

        if( length(gcns) > 1 )
        {
            maxCN = ifelse( length(gcns[gcns >= 2]) > 0, max(gcns[gcns >= 2]), NA )
            minCN = ifelse( length(gcns[gcns  < 2]) > 0, min(gcns[gcns  < 2]), NA )
        
            if( is.na(maxCN) )
            {
                maxCNA = NA
            }else{
                if( maxCN == 2 )
                {
                    maxCNA = "Neu"
                }else{
                    maxCNA = ifelse( maxCN >= 8, "Amp+", "Amp" )
                }
            }
            if( is.na(minCN) )
            {
                minCNA = NA
            }else{
                minCNA = ifelse( minCN == 0, "Del-", "Del" )
            }
            call.cn  = c(minCN, maxCN)
            call.cn  = call.cn[call.cn != 2]
            call.cna = c(minCNA, maxCNA)
            call.cna = call.cna[call.cna != "Neu"]

            gnCN$cn  = paste(call.cn[!is.na(call.cn)],   collapse=",")
            gnCN$cna = paste(call.cna[!is.na(call.cna)], collapse=",")
        }
        return(unique(gnCN))
    }, mc.cores = THREADS ))[,-1] %>%
        select(c("gene","chr","location","cn","cna","symbol","entrez","ens_geneid","status")) %>%
        dplyr::rename( 
            GENE=gene,
            CHR=chr,
            CN=cn,
            HGNC_SYMBOL=symbol,
            LOCATION=location,
            ENTREZ=entrez,
            ENS_GENEID=ens_geneid,
            CNA=cna,
            HGNC_STATUS=status
        )
    #--------------------------------------------------------------------------#
    GENE_dual <- GENE_CNA[grep(",", GENE_CNA$CNA), ]
    GENE_amp  <- GENE_CNA %>% filter( CNA %in% c("Amp","Amp+") ) %>% mutate(CN = as.numeric(CN)) %>% arrange(desc(CN))
    GENE_del  <- GENE_CNA %>% filter( CNA %in% c("Del","Del-") ) %>% mutate(CN = as.numeric(CN)) %>% arrange(desc(CN))
    #--------------------------------------------------------------------------#
    GENE_AMP <- rbind( GENE_dual, GENE_amp )
    GENE_DEL <- rbind( GENE_dual, GENE_del )
    #--------------------------------------------------------------------------#
    GENE_PANEL_CNV <- GENE_CNA %>% 
        filter( HGNC_SYMBOL %in% Reduce(union, gcx_preset_genes[c("tso500","ngs_essential_genes")]) ) %>%
        filter( CNA != "Neu" )
    #--------------------------------------------------------------------------#
    GENE_dual$amp   <- sapply(GENE_dual$CNA, function(q) unlist(strsplit(q, ","))[2] )
    GENE_dual$del   <- sapply(GENE_dual$CNA, function(q) unlist(strsplit(q, ","))[1] )    
    GENE_PANEL_CNV2 <- rbind(
        GENE_dual %>% 
            filter( HGNC_SYMBOL %in%  Reduce(union, gcx_preset_genes[c("tso500","ngs_essential_genes")]) ) %>%
            filter( HGNC_SYMBOL %nin% gcx_preset_genes$hla_blacklist ) %>%
            filter( amp == "Amp+" | del == "Del-" ) %>%
            select(c(1:ncol(GENE_PANEL_CNV))),
        GENE_PANEL_CNV %>% 
            filter( HGNC_SYMBOL %nin% gcx_preset_genes$hla_blacklist ) %>% 
            filter( CNA %in% c("Amp+","Del-") )
    )
    #--------------------------------------------------------------------------#
    write.table( GENE_CNA, 
        sprintf("%s/%s.%s.CNV.Genes.Raw.txt",CNVKIT_SUMMARY_DIR, SEQ_ID, ANALYSIS_MODE), 
        quote=F, col.names=T, row.names=F, sep="\t" 
    )
    write.table( GENE_AMP, 
        sprintf("%s/%s.%s.CNV.Genes.Amplification.txt",CNVKIT_SUMMARY_DIR, SEQ_ID, ANALYSIS_MODE), 
        quote=F, col.names=T, row.names=F, sep="\t" 
    )
    write.table( GENE_DEL, 
        sprintf("%s/%s.%s.CNV.Genes.Deletion.txt",CNVKIT_SUMMARY_DIR, SEQ_ID, ANALYSIS_MODE), 
        quote=F, col.names=T, row.names=F, sep="\t" 
    )
    write.table( GENE_PANEL_CNV, 
        sprintf("%s/%s.%s.CNV.CancePanelGenes.txt",CNVKIT_SUMMARY_DIR, SEQ_ID, ANALYSIS_MODE), 
        quote=F, col.names=T, row.names=F, sep="\t" 
    )
    write.table( GENE_PANEL_CNV2, 
        sprintf("%s/%s.%s.CNV.CancePanelGenes.Siginificant.txt",CNVKIT_SUMMARY_DIR, SEQ_ID, ANALYSIS_MODE), 
        quote=F, col.names=T, row.names=F, sep="\t" 
    )
#--------------------------------------------------------------------------------------------------#


    # CN0_LINE  = 4
    # CN14_LINE = 11 

    # SEG %>% filter( chromosome %in% CHR) %>% ggplot() + 
    #     coord_cartesian(ylim=c(-4,12)) +
    #     theme(
    #         panel.background = element_rect(fill='white'),
    #         axis.line.y      = element_line(colour = "#7a7a7a"),
    #         axis.ticks.x     = element_blank(),
    #         axis.text.x      = element_text(hjust = 0 )
    #     ) + 
    #     labs( x = "Chromosome", y = "CN" ) +
    #     geom_point( aes( x = mid2, y = log2 ), size = 0.6, col = "#d6d6d6", alpha = 0.4 ) +   
    #     geom_hline( yintercept = 0, col = "#7a7a7a", alpha = 0.3 ) +
    #     geom_vline( xintercept = ChrLength %>% filter( chr %in% CHR ) %>% .$cum_chr_length, linewidth = 0.8, col = "#7a7a7a", alpha = 0.2 ) +
    #     geom_hline( yintercept = c( CN0_LINE:CN14_LINE ), linetype = "dashed", linewidth = 0.25, col = "#60A917", alpha = 0.4 ) +
    #     geom_segment( data = ANS_CNVKIT %>% filter( chromosome %in% CHR), 
    #         aes( x = start2, y = cn*0.5 + CN0_LINE, xend = end2, yend = cn*0.5 + CN0_LINE ), col = "#C73500", size = 3, alpha = 0.6 ) +
    #     geom_segment( data = CALL %>% filter( chromosome %in% CHR), 
    #         aes(x = start2, y = cn*0.5 + CN0_LINE, xend = end2, yend = cn*0.5 + CN0_LINE ), col = "#006EAF", size = 3, alpha = 0.6) +
    #     scale_y_continuous( breaks = c(0, CN0_LINE:CN14_LINE ), 
    #         labels=c(sprintf("Ploidy %s", PLOIDY),"0","2","4","6","8","10","12","14") ) +
    #     scale_x_continuous( breaks = ChrLength %>% filter( chr %in% CHR ) %>% .$add_X, 
    #         labels = ChrLength %>% filter( chr %in% CHR ) %>% .$chr) 


