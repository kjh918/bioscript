

#
#

#---| PACKAGES |-----------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("RMySQL"))
#--------------------------------------------------------------------------------------------------#

#---| SET DEFAULT VALUES |-------------------------------------------------------------------------#
    BASE_DIR                = "/data/wts"
    DRAW_STAR_FUSION_CIRCOS = TRUE
    GENOME_ASSEMBLY         = "hg19"
    DB_IMPORT               = TRUE
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"),                action="store", default="/data/wts", type="character", help="BASE_DIR"),   
        make_option(c("--SEQ_FOLDER"),              action="store", default=NA,          type="character", help="SEQ_FOLDER"),
        make_option(c("--SEQ_ID"),                  action="store", default=NA,          type="character", help="$SEQ_ID"),
        make_option(c("--DRAW_STAR_FUSION_CIRCOS"), action="store", default=TRUE,        type="logical",   help="DRAW CIRCOS-PLOT from STAR-FUSION"),
        make_option(c("--GENOME_ASSEMBLY"),         action="store", default=NA,          type="character", help="GENOME_ASSEMBLY"),
        make_option(c("--DB_IMPORT"),               action="store", default=TRUE,        type="logical",   help="IMPORT RESULT INTO DATABASE")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR                <- ARGS$BASE_DIR
    SEQ_FOLDER              <- ARGS$SEQ_FOLDER
    SEQ_ID                  <- ARGS$SEQ_ID
    DRAW_STAR_FUSION_CIRCOS <- ARGS$DRAW_STAR_FUSION_CIRCOS
    GENOME_ASSEMBLY         <- ARGS$GENOME_ASSEMBLY
    DB_IMPORT               <- ARGS$DB_IMPORT
#--------------------------------------------------------------------------------------------------#

#---| FUNCTIONS |----------------------------------------------------------------------------------#
#' @description filtering, mergeing, summarization of arriba & star-fusion results 
#' @import Hmisc
#' @import plyr
#' @import dplyr
#' @param ARRIBA_RESULT      arriba analysis result
#' @param STAR_FUSION_RESULT star-fusion analysis result
#' @export
    FusionGeneResultSummarize <- function( ARRIBA_RESULT, STAR_FUSION_RESULT )
    {
        require(Hmisc)
        require(plyr)
        require(dplyr)
        
        res_ar <- ARRIBA_RESULT %>% mutate( fusionKey = paste(breakpoint1, breakpoint2, sep="_") ) %>% 
            dplyr::rename(gene1=1) %>% filter( confidence %in% c("high","medium"))
        res_sf <- STAR_FUSION_RESULT %>% mutate(
            geneA = sapply(LeftGene, function(y) unlist(strsplit(y, "\\^"))[1]),
            geneB = sapply(RightGene, function(z) unlist(strsplit(z, "\\^"))[1]), 
            geneA_strand = sapply(LeftBreakpoint, function(q) tail(unlist(strsplit(q, "\\:")), 1) ),
            geneB_strand = sapply(RightBreakpoint, function(w) tail(unlist(strsplit(w, "\\:")), 1) ),
            LeftBreakpoint2  = gsub("^chr", "", gsub("\\:\\-$", "", gsub("\\:\\+$", "", LeftBreakpoint ))) ,
            RightBreakpoint2 = gsub("^chr", "", gsub("\\:\\-$", "", gsub("\\:\\+$", "", RightBreakpoint))) 
        ) %>% mutate( fusionKey = paste(LeftBreakpoint2, RightBreakpoint2, sep="_") ) %>% 
            dplyr::rename(FusionName=1)

        CommonFusion <- intersect( res_ar$fusionKey, res_sf$fusionKey )
        SfOnlyFusion <- res_sf$fusionKey[ res_sf$fusionKey %nin% CommonFusion ]
        ArOnlyFusion <- res_ar$fusionKey[ res_ar$fusionKey %nin% CommonFusion ]

        CommonFusionReport <- res_sf %>% 
            filter( fusionKey %in% CommonFusion ) %>% 
            select(c("FusionName","geneA","geneB","geneA_strand","geneB_strand","FFPM","JunctionReadCount","SpanningFragCount","fusionKey")) %>%
            dplyr::rename( junctionReads=JunctionReadCount, spanningReads=SpanningFragCount ) %>%
            left_join(
                ., 
                ( res_ar %>% filter( fusionKey %in% CommonFusion ) %>% 
                    select(c("gene1","gene2","strand1.gene.fusion.","strand2.gene.fusion.","breakpoint1","breakpoint2","type","confidence","split_reads1","split_reads2","discordant_mates","site1","site2","fusionKey")) %>%
                    dplyr::rename( gene1_strand='strand1.gene.fusion.', gene2_strand='strand2.gene.fusion.', 
                        breakPoint1=breakpoint1, breakPoint2=breakpoint2, splitReads1=split_reads1, splitReads2=split_reads2 )
                ), 
                by = "fusionKey"
            ) %>% select(!c("fusionKey"))
        if( nrow(CommonFusionReport) == 0 ){ CommonFusionReport= "No_Common_FusionGenes"}

        SfOnlyReport <- res_sf %>% filter( fusionKey %in% SfOnlyFusion ) %>% 
            select(c("FusionName","geneA","geneB","geneA_strand","geneB_strand","FFPM","JunctionReadCount","SpanningFragCount"))
        if( nrow(SfOnlyReport) == 0 ){ SfOnlyReport= "No_StarFusion-Only_FusionGenes"}

        ArOnlyReport <- res_ar %>% filter( fusionKey %in% ArOnlyFusion ) %>% 
            select(c("gene1","gene2","strand1.gene.fusion.","strand2.gene.fusion.","breakpoint1","breakpoint2","type","confidence","split_reads1","split_reads2","discordant_mates","site1","site2"))
        if( nrow(ArOnlyReport) == 0 ){ ArOnlyReport= "No_Arriba-Only_FusionGenes"}

        CommonFusionArribaFormat <- res_ar %>% filter( fusionKey %in% CommonFusion ) %>% select(!c("fusionKey"))
        if( nrow(CommonFusionArribaFormat) == 0 ){ CommonFusionArribaFormat= "No_Common_FusionGenes"}

        CommonFusionStarFusionFormat <- res_sf %>% filter( fusionKey %in% CommonFusion ) %>% select(!c("fusionKey"))
        if( nrow(CommonFusionStarFusionFormat) == 0 ){ CommonFusionStarFusionFormat= "No_Common_FusionGenes"}
        

        FUSION_RESULTS <- list(
            Arriba_StarFusion_Common = CommonFusionReport,
            Arriba_Only              = ArOnlyReport,
            StarFusion_Only          = SfOnlyReport,
            Arriba_CommonFusion      = CommonFusionArribaFormat,
            StarFusion_CommonFusion  = CommonFusionStarFusionFormat
        )

        return(FUSION_RESULTS)
    }
#--------------------------------------------------------------------------------------------------#

#---| READ FUSION-GENES and SPLICE-VARIANTS RESULT |-----------------------------------------------#
    message(">>> Read Fusion-Genes and Splice-Variants Results.")
    #--------------------------------------------------------------------------#
    ArribaResDir     <- sprintf( "%s/%s/%s/fusionGenes_spliceVariants/arriba", BASE_DIR, SEQ_FOLDER, SEQ_ID )
    StarFusionResDir <- sprintf( "%s/%s/%s/fusionGenes_spliceVariants/star_fusion", BASE_DIR, SEQ_FOLDER, SEQ_ID )
    SpliceVarResDir  <- sprintf( "%s/%s/%s/fusionGenes_spliceVariants/splice_variants", BASE_DIR, SEQ_FOLDER, SEQ_ID )
    #--------------------------------------------------------------------------#
    ResArriba     <- read.delim(sprintf("%s/%s.arriba.fusion.genes.tsv", ArribaResDir, SEQ_ID))
    ResStarFusion <- read.delim(sprintf("%s/star-fusion.fusion_predictions.abridged.tsv", StarFusionResDir))
    ResSpliceVars <- read.delim(sprintf("%s/%s.cancer.introns", SpliceVarResDir,SEQ_ID))
#--------------------------------------------------------------------------------------------------#
    
#---| FILTERING and SUMMARY |----------------------------------------------------------------------#
    message(">>> Summary Fusion-Genes Result.")
    #--------------------------------------------------------------------------#
    # FusionGenes
    FusionGeneResultList <- FusionGeneResultSummarize( ARRIBA_RESULT=ResArriba, STAR_FUSION_RESULT=ResStarFusion )
    #--------------------------------------------------------------------------#
    message(">>> Summary Splice-Variants Result.")
    #--------------------------------------------------------------------------#
    # SpliceVariants
    cancerIntrons_SV    <- ResSpliceVars %>% filter( uniq_mapped >= 5, !is.na(variant_name) )
    cancerIntrons_NonSV <- ResSpliceVars %>% filter( uniq_mapped >= 5,  is.na(variant_name) )
    if( nrow(cancerIntrons_SV)    == 0 ){ cancerIntrons_SV    = "No_SplicieVariants" }
    if( nrow(cancerIntrons_NonSV) == 0 ){ cancerIntrons_NonSV = "No_CancerIntrons" }
    SpliceVariantsList <- list(
        Cancer_SpliceVariants = cancerIntrons_SV,
        Cancer_Introns = cancerIntrons_NonSV
    )
#--------------------------------------------------------------------------------------------------#

#---| SAVE RESULTS |-------------------------------------------------------------------------------#
    FUSION_SPLICEVARS <- c(FusionGeneResultList[1:3], SpliceVariantsList)
    #--------------------------------------------------------------------------#
    ResultDir <- sprintf( "%s/%s/%s/fusionGenes_spliceVariants", BASE_DIR, SEQ_FOLDER, SEQ_ID )
    #--------------------------------------------------------------------------#
    # all results summary
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    message(">>> Save Results.")
    cat(">>> ALL RESULTS INTEGRATED XLSX-FILE.....")
    #--------------------------------------------------------------------------#
    openxlsx::write.xlsx( FUSION_SPLICEVARS, sprintf("%s/%s.FusionGenes.SpliceVariants.CancerIntrons.Summary.xlsx", ResultDir, SEQ_ID), rowNames=FALSE, overwrite=TRUE )
    #--------------------------------------------------------------------------#
    message("DONE.")
    cat(">>> COMMON FUSION-GENES.....")
    #--------------------------------------------------------------------------#
    # common-fusion-genes
    if( class(FusionGeneResultList$Arriba_StarFusion_Common) != "data.frame" )
    {
        write(
            "No Common FusionGenes Found.",
            sprintf("%s/%s.FusionGenes.Arriba.StarFusion.Common_NoCommonFound.tsv", ResultDir, SEQ_ID)
        )
        RESULT_COMMON_FUSION <- data.frame()
    }else{
        write.table( FusionGeneResultList$Arriba_StarFusion_Common, 
            sprintf("%s/%s.FusionGenes.Arriba.StarFusion.Common.tsv", ResultDir, SEQ_ID), 
            quote=F, row.names=F, col.names=T, sep="\t"
        )
        RESULT_COMMON_FUSION <- FusionGeneResultList$Arriba_StarFusion_Common %>% 
            dplyr::select(c("geneA","geneB","FFPM","confidence","gene1_strand","gene2_strand","breakPoint1","breakPoint2","type",
                "junctionReads","spanningReads","splitReads1","splitReads2","discordant_mates","site1","site2")) %>%
            dplyr::rename(
                FFPM_starfusion=FFPM, confidence_arriba=confidence, geneA_strand=gene1_strand, geneB_strand=gene2_strand, 
                geneA_breakpoint=breakPoint1, geneB_breakpoint=breakPoint2, type_arriba=type, 
                junction_rds_starfusion=junctionReads, span_rds_starfusion=spanningReads, junction_rds_arriba=splitReads1, span_rds_arriba=splitReads2,
                discordant_rds=discordant_mates, geneA_site_arriba=site1, geneB_site_arriba=site2
            ) %>% mutate( analysis_method="starfusion,arriba")
    }
    #--------------------------------------------------------------------------#
    message("DONE.")
    cat(">>> ARRIBA-ONLY FUSION-GENES.....")
    #--------------------------------------------------------------------------#
    # arriba-only fusion genes 
    if( class(FusionGeneResultList$Arriba_Only) != "data.frame" )
    {
        write( 
            "No Arriba FusionGenes Found.",
            sprintf("%s/%s.FusionGenes.Arriba.Only_NotFound.tsv", ResultDir, SEQ_ID)
        )
        RESULT_ARRIBA_ONLY <- data.frame()
    }else{
        write.table( FusionGeneResultList$Arriba_Only, 
            sprintf("%s/%s.FusionGenes.Arriba.Only.tsv", ResultDir, SEQ_ID), 
            quote=F, row.names=F, col.names=T, sep="\t"
        )
        RESULT_ARRIBA_ONLY <- FusionGeneResultList$Arriba_Only %>% mutate(
            FFPM_starfusion=NA, junction_rds_starfusion=NA, span_rds_starfusion=NA
        ) %>% dplyr::select(
            c("gene1","gene2","FFPM_starfusion","confidence","strand1.gene.fusion.","strand2.gene.fusion.","breakpoint1","breakpoint2","type",
            "junction_rds_starfusion","span_rds_starfusion","split_reads1","split_reads2","discordant_mates","site1","site2")
        ) %>% dplyr::rename(
            geneA=gene1, geneB=gene2, confidence_arriba=confidence, geneA_strand=strand1.gene.fusion., geneB_strand=strand2.gene.fusion.,
            geneA_breakpoint=breakpoint1, geneB_breakpoint=breakpoint2, type_arriba=type, junction_rds_arriba=split_reads1, span_rds_arriba=split_reads2,
            discordant_rds=discordant_mates, geneA_site_arriba=site1, geneB_site_arriba=site2
        ) %>% mutate( analysis_method="arriba")
    }
    #--------------------------------------------------------------------------#
    message("DONE.")
    cat(">>> STAR-FUSION-ONLY FUSION-GENES.....")
    #--------------------------------------------------------------------------#
    # star-fusion-only fusion genes 
    if( class(FusionGeneResultList$StarFusion_Only) != "data.frame" )
    {
        write(
            "No STAR-Fusion FusionGenes Found.",
            sprintf("%s/%s.FusionGenes.StarFusion.Only_NotFound.tsv", ResultDir, SEQ_ID)
        )
        RESULT_STARFUSION_ONLY <- data.frame()
    }else{
        write.table( FusionGeneResultList$StarFusion_Only, 
            sprintf("%s/%s.FusionGenes.StarFusion.Only.tsv", ResultDir, SEQ_ID), 
            quote=F, row.names=F, col.names=T, sep="\t"
        )
        RESULT_STARFUSION_ONLY <- FusionGeneResultList$StarFusion_Only %>% mutate(
            confidence_arriba=NA, geneA_breakpoint=NA, geneB_breakpoint=NA, type_arriba=NA, junction_rds_arriba=NA, span_rds_arriba=NA,discordant_rds=NA, 
            geneA_site_arriba=NA, geneB_site_arriba=NA
        ) %>% dplyr::select(
            c("geneA","geneB","FFPM","confidence_arriba","geneA_strand","geneB_strand","geneA_breakpoint","geneB_breakpoint","type_arriba",
            "JunctionReadCount","SpanningFragCount","junction_rds_arriba","span_rds_arriba","discordant_rds","geneA_site_arriba","geneB_site_arriba")
        ) %>% dplyr::rename(
            FFPM_starfusion=FFPM, junction_rds_starfusion=JunctionReadCount, span_rds_starfusion=SpanningFragCount
        ) %>% mutate( analysis_method="starfusion")
    }
    #--------------------------------------------------------------------------#
    message("DONE.")
    cat(">>> SPLICE-VARIANTS.....")
    #--------------------------------------------------------------------------#
    # splicing-variants  
    if( class(SpliceVariantsList$Cancer_SpliceVariants) != "data.frame" )
    {
        write(
            "No Splice-Variants Found.",
            sprintf("%s/%s.Cancer.Splicing.Variants_NotFound.tsv", ResultDir, SEQ_ID),
        )
        RESULT_SPLICE_VARIANTS <- data.frame()
    }else{
        write.table( SpliceVariantsList$Cancer_SpliceVariants, 
            sprintf("%s/%s.Cancer.Splicing.Variants.tsv", ResultDir, SEQ_ID), 
            quote=F, row.names=F, col.names=T, sep="\t"
        )
        RESULT_SPLICE_VARIANTS <- SpliceVariantsList$Cancer_SpliceVariants %>%
            mutate( group = "cancer_splice_variant")
    }
    #--------------------------------------------------------------------------#
    message("DONE.")
    cat(">>> CANCER-INTRONS.....")
    #--------------------------------------------------------------------------#
    # cancer introns  
    if(  class(SpliceVariantsList$Cancer_Introns) != "data.frame" )
    {
        write(
            "No Cancer-Intron Found.",
            sprintf("%s/%s.Cancer.Introns_NotFound.tsv", ResultDir, SEQ_ID), 
        )
        RESULT_CANCER_INTRON <- data.frame()
    }else{
        write.table( SpliceVariantsList$Cancer_Introns, 
            sprintf("%s/%s.Cancer.Introns.tsv", ResultDir, SEQ_ID), 
            quote=F, row.names=F, col.names=T, sep="\t"
        )
        RESULT_CANCER_INTRON <- SpliceVariantsList$Cancer_Introns %>% 
        mutate( group = "cancer_intron" )
    }
    #--------------------------------------------------------------------------#
    message("DONE.")
    #--------------------------------------------------------------------------#  
    # common-fusiongenes arriba/star-fusion format
    if( class(FusionGeneResultList$Arriba_CommonFusion) != "data.frame" )
    {
        write(
            "No Common FusionGenes Found.",
            sprintf("%s/%s.FusionGenes.Common.Arriba.Format_NotFound.tsv", ResultDir, SEQ_ID)
        )
        write(
            "No Common FusionGenes Found.",
            sprintf("%s/%s.FusionGenes.Common.StarFusion.Format_NotFound.tsv", ResultDir, SEQ_ID)
        )
    }else{    
        # arriba format
        cat(">>> COMMON FUSION-GENES IN ARRIBA-FORMAT FOR PLOT.....")
        write.table( FusionGeneResultList$Arriba_CommonFusion %>% dplyr::rename('#gene1'=1, 'strand1(gene/fusion)'=3, 'strand2(gene/fusion)'=4), 
            sprintf("%s/%s.FusionGenes.Common.Arriba.Format.tsv", ResultDir, SEQ_ID), 
            quote=F, row.names=F, col.names=T, sep="\t"
        )
        message("DONE.")
        # star-fusion format
        cat(">>> COMMON FUSION-GENES IN STAR_FUSION-FORMAT FOR PLOT.....")
        write.table( FusionGeneResultList$StarFusion_CommonFusion %>% dplyr::rename('#FusionName'=1), 
            sprintf("%s/%s.FusionGenes.Common.StarFusion.Format.tsv", ResultDir, SEQ_ID), 
            quote=F, row.names=F, col.names=T, sep="\t"
        )
        message("DONE.")
    }
    #--------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#---| DATABASE IMPORT |----------------------------------------------------------------------------#
    if( DB_IMPORT )
    {
        FUSION_GENE_RESULTS <- rbind( RESULT_COMMON_FUSION, RESULT_ARRIBA_ONLY, RESULT_STARFUSION_ONLY ) %>% 
            mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID )

        SPLICE_VARIANT_RESULTS <- rbind( RESULT_SPLICE_VARIANTS, RESULT_CANCER_INTRON ) %>% 
            mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID )
   


        if( nrow(FUSION_GENE_RESULTS) > 0 )
        {
            source("/data/wts/params/ruo_wts_db.R")
            dbCon            <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData1 <- dbGetQuery(dbCon, sprintf("DELETE FROM fusion_genes WHERE seq_folder = '%s' AND seq_id = '%s'", SEQ_FOLDER, SEQ_ID))
            UpdateData1      <- dbWriteTable(dbCon, name="fusion_genes", value=FUSION_GENE_RESULTS, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
        }
        if( nrow(SPLICE_VARIANT_RESULTS) > 0 )
        {
            dbCon            <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=3306, password=DB_PW, db=DB_NAME_DATA )
            DeleteExtstData2 <- dbGetQuery(dbCon, sprintf("DELETE FROM splice_variants WHERE seq_folder = '%s' AND seq_id = '%s'", SEQ_FOLDER, SEQ_ID))
            UpdateData2      <- dbWriteTable(dbCon, name="splice_variants", value=SPLICE_VARIANT_RESULTS, row.names=FALSE, append=TRUE)
            dbDisconnect(dbCon)
        }
    }
#--------------------------------------------------------------------------------------------------#

#---| STAR-FUSION CIRCOS-PLOT |--------------------------------------------------------------------#
    if( DRAW_STAR_FUSION_CIRCOS )
    {
        message("\n>>> CREATE Circos-Plots....")
        #----------------------------------------------------------------------#
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        suppressPackageStartupMessages(library("chimeraviz"))
        source("/storage/home/kangsm/myScripts/R_FUNCTIONAL_MODULES/FunMod_StarFusionResults.R")
        #----------------------------------------------------------------------#
        InputResultFile <- sprintf("%s/%s.FusionGenes.Common.StarFusion.Format.tsv", ResultDir, SEQ_ID)
        if( !file.exists(InputResultFile) )
        {
            "No STAR-Fusion format Common FusionGene result file found. SKIP circos-plot drawing. "
        }else{
            fusionGenesList = import_starfusion(
                filename       = sprintf("%s/%s.FusionGenes.Common.StarFusion.Format.tsv", ResultDir, SEQ_ID), 
                genome_version = GENOME_ASSEMBLY
            )
            #------------------------------------------------------------------#
            png( sprintf("%s/%s.FusionGenes.Common.Circos.Plot.png", ResultDir, SEQ_ID), width=5, height=5, units="in", res=200, type="cairo" )
            plot_circle(fusionGenesList)
            dev.off()
            message(">>> CREATE Circos-Plots....DONE.")
        }
        #----------------------------------------------------------------------#        
    }else{
        message("\n>>> CREATE Circos-Plots SKIPPED.")
    }
#--------------------------------------------------------------------------------------------------#

    













