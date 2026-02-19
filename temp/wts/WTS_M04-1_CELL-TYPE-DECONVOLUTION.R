
#---| PACKAGES |------------------------------------------------------------------------------------
    suppressPackageStartupMessages(library("SCdeconR"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("openxlsx"))
    suppressPackageStartupMessages(library("Cairo"))
    suppressPackageStartupMessages(library("ggplot2"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
#---------------------------------------------------------------------------------------------------

#---| FUNCTIONS |-----------------------------------------------------------------------------------
    SCdeconR.Predicted.Proportion.Barplot <- function( prop, sampleOrder=NULL )
    {
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("Cairo"))
        grDevices::X11.options(type='cairo')
        options(device='x11')
        options(bitmapType='cairo')
        #
        prop <- as.data.frame(t(prop))
        prop <- prop %>% mutate(sampleid = rownames(prop)) %>% reshape2::melt(id.var = "sampleid",
            value.name = "ct_prop", variable.name = "ct") %>% as.data.frame()
        #
        if( !is.null(sampleOrder) )
        {
            prop$sampleid = factor( prop$sampleid, levels=sampleOrder )
        }
        #        
        CellTypeColor = c(
            "Epithelial cells"  = "#B85450",
            "Monocytes"         = "#6C8EBF",
            "Fibroblast"        = "#60A917",
            "T cell"            = "#1BA1E2",
            "NK cells"          = "#006EAF",
            "Endothelial cells" = "#FF7F0E",
            "Pericytes"         = "#FFBF00",
            "Adipocyte"         = "#9673A6"
        )
        gp <- ggplot() + geom_col(data=prop, aes(x = sampleid, y = ct_prop, fill = ct)) + 
            labs(x="", y = "Predicted proportion", fill = "Cell types") + 
            scale_y_continuous(expand=c(0,0)) +
            scale_fill_manual(values=CellTypeColor) +
            theme(
                panel.background = element_rect(fill='white'),
                axis.line.x      = element_line(color="#2d2d2d"), 
                axis.line.y      = element_line(color="#2d2d2d"),
                panel.grid.major = element_line(color='#dcdbdb', linetype='dotted'),
                axis.title       = element_text(size = 20, face = "bold", family="Nunito"), 
                axis.text.x      = element_text(size = 15, family="Nunito", angle=90, hjust=1 ), 
                axis.text.y      = element_text(size = 15, family="Nunito" ), 
                axis.ticks.x     = element_blank(),
                legend.title     = element_text(size = 15, face = "bold", family="Nunito"),
                legend.text      = element_text(size = 12),
                plot.margin      = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
            )
        return(gp)
    }
#---------------------------------------------------------------------------------------------------

#---| LOAD DATA |-----------------------------------------------------------------------------------
    refdata         <- readRDS("/storage/home/kangsm/myDB/rds/SCdeconR_refdata.RDS")
    hg19_GTF        <- readRDS("/storage/home/kangsm/myDB/rds/hg19_gtf.RDS")
    FeatureGenesGTF <- readRDS("/storage/home/kangsm/myDB/rds/SCdeconR_FeatureGenes_hg19.RDS")
    phenodata       <- readRDS("/storage/home/kangsm/myDB/rds/SCdeconR_ref_phenodata.RDS")
    source("/storage/home/kangsm/runScripts/WTS_Fun.AnalysisModules.R")
#---------------------------------------------------------------------------------------------------

#---| ANALYSIS INFO |-------------------------------------------------------------------------------
    BatchID <- "WTS_25_02"
#---------------------------------------------------------------------------------------------------


#---| INPUT DATA |
    dbCon <- dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='gcx_ngs_service' )
    SampleInfo <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", BatchID))
    dbDisconnect(dbCon)
    CountMatrix <- Create.Expression.Values.Matrix( SeqIdList=SampleInfo$seq_id, ValueType="count", DataSource="rds", CountAsInteger=TRUE, Convert2Matrix=TRUE )
    colnames(CountMatrix) <- gsub(" ", "_", SampleInfo[match(colnames(CountMatrix), SampleInfo$seq_id), "sample_name"])
#---------------------------------------------------------------------------------------------------


#---| RUN ANALYSIS |-------------------------------------------------------------------------------#

    CountMatrixRev <- CountMatrix[FeatureGenesGTF$gene_id, ] %>% as.data.frame() %>% mutate( gene = FeatureGenesGTF$gene )
    InputRNaseq    <- CountMatrixRev %>% group_by( gene ) %>% summarise_all(sum) %>% as.data.frame()
    rownames(InputRNaseq) <- InputRNaseq$gene

    deconvolution_result <- scdecon(
        bulk              = InputRNaseq[,-1], 
        ref               = GetAssayData(refdata, layer = "data", assay = "SCT"),
        phenodata         = phenodata, 
        filter_ref        = FALSE, 
        decon_method      = "OLS", 
        norm_method_sc    = "LogNormalize",
        norm_method_bulk  = "TMM", 
        trans_method_sc   = "none", 
        trans_method_bulk = "log2", 
        marker_strategy   = "all",
        min_pct_ct        = 0.03,
        lfc_markers       = log2(1.3)
    )
#--------------------------------------------------------------------------------------------------#

#---| SAVE RESULTS
    OutDir <- sprintf("/data/wts/%s/analysis/deconvolution", BatchID)
    if( !dir.exists(OutDir) ){ system(sprintf("mkdir -p %s", OutDir)) }

    PredictTable  <- deconvolution_result$prediction 
    PredictResult <- data.frame(CELL_TYPE=rownames(PredictTable), PredictTable) 
    OutputPredictTable <- sprintf("%s/%s_cell.type.deconvolution.result.xlsx", OutDir, BatchID)
    write.xlsx( PredictResult, file=OutputPredictTable, rowName=FALSE, overwrite=TRUE )

    deconv_barplot <- SCdeconR.Predicted.Proportion.Barplot(prop=deconvolution_result[[1]], sampleOrder=NULL)
    OutputBarplot  <- sprintf("%s/%s_cell.type.deconvolution.png", OutDir, BatchID)
    ggsave(deconv_barplot, file=OutputBarplot, width=14, height=8, unit='in', dpi=150 )
#---------------------------------------------------------------------------------------------------

