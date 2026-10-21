
#### packages
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("RMySQL")
        library("data.table")
        library("parallel")
        library("optparse")
        library("ggplot2")
        library("Cairo")
        library("SingleCellExperiment")
        library("Seurat")
        library("patchwork")
        library("SingleR")
        library("celldex")
    })
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
    Rounding <- function(x, digits = 0) {   
        posneg <- sign(x)
        z <- abs(x) * 10^digits
        z <- floor(z + 0.5)
        v <- z / 10^digits * posneg
        return(v)
    }
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
####

#### PARAMS 
BaseDir <- "/data/scRNAseq/02.Public-Datasets"
BatchID <- "002.Tabula-Muris-Senis"
WorkDir <- sprintf("%s/%s/analysis/scRNAseq", BaseDir, BatchID)
sample_list <- c("10X_P4_2","10X_P4_7","10X_P8_13")
# sample_list <- list.files(sprintf("%s/%s/scRNAseq", BaseDir, BatchID))
DataTax <- "mouse"
####

#### CELL-RANGER QC RESULTS
    cellranger_qc_res <- lapply( indexing(sample_list), function(sid) 
        data.table::fread(sprintf("%s/%s/scRNAseq/%s/outs/metrics_summary.csv", BaseDir, BatchID, sid)) 
    ) %>% ldply(.id="Sample_ID") %>% mutate(Sample_ID=as.character(Sample_ID))  
    colnames(cellranger_qc_res) <- gsub(" ", "_", colnames(cellranger_qc_res))
    write.table( cellranger_qc_res, sprintf("%s/scRNAseq_cellranger.qc.output.tsv", WorkDir), quote=F, col.names=T, row.names=F, sep="\t" )
####

#### READ DATA 
    # data dir list
    data_dirs_list <- lapply(indexing(sample_list), function(sid) sprintf("%s/%s/scRNAseq/%s/outs/filtered_feature_bc_matrix", BaseDir, BatchID, sid))
    # read data as seurat-object
    sce_list <- lapply(indexing(sample_list), function(sid) {
        counts <- Read10X(data.dir = data_dirs_list[[sid]])
        sce    <- CreateSeuratObject( counts=counts, project=sid, min.cells=3, min.features=200 )
        sce$sample <- sid
        return(sce)
    })
####

#### BASIC QC
    # calculate mitochondrial ratio
    sce_list <- lapply(sce_list, function(sce) { 
        sce[["pct.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-|^MT-") 
        return(sce)
    })
    # cell numbers check
    sce_cell_numbers <- sapply(sce_list, ncol)
    # qc-summary (before filter)
    sce_qc_summary_before <- lapply(names(sce_list), function(sid) {
        md <- sce_list[[sid]]@meta.data
        data.frame(
            sample = sid,
            cells  = nrow(md),
            median_nFeature = median(md$nFeature_RNA),
            median_nCount   = median(md$nCount_RNA),
            median_percent_mt = median(md$pct.mt)
        )
    }) %>% bind_rows()
    # qc-plots (before filter)
    p_qc.violin.before <- cowplot::plot_grid(
        plotlist = lapply( names(sce_list), function(sid) {
            p.viloin <- VlnPlot( sce_list[[sid]], layer="counts", features=c("nFeature_RNA","nCount_RNA","pct.mt"), ncol=3, pt.size=0 ) + 
                plot_annotation(title = sid) 
            return(p.viloin)
        }),
        nrow=1
    )
    p_feat.scatter.before <- cowplot::plot_grid(
        plotlist = lapply( names(sce_list), function(sid) {
            p.feat.scatter <- FeatureScatter( sce_list[[sid]], feature1="nCount_RNA", feature2="nFeature_RNA" ) + 
                ggtitle(sid) + theme(plot.title=element_text(face="bold", size=14, family="Roboto"))
            return(p.feat.scatter)
        }),
        nrow=1
    )
    p_qc.summary.before <- cowplot::plot_grid( p_qc.violin.before, p_feat.scatter.before, nrow=2 )
    ggsave( p_qc.summary.before, 
        file=sprintf("%s/scRNAseq_qc.summary.before.filtering.png", WorkDir),
        width=8*length(sce_list), height=12, unit='in', dpi=150
    )

    # qc-simulation
    message(sprintf("[ %s ] run QC threshold simulation...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_qc_simulation <- lapply(names(sce_list), function(sid) {
        md <- sce_list[[sid]]@meta.data
        data.frame(
            sample = sid,
            total  = nrow(md),
            keep_mt10_feat300 = sum( md$nFeature_RNA > 300 & md$pct.mt < 10 ),
            keep_mt15_feat300 = sum( md$nFeature_RNA > 300 & md$pct.mt < 15 ),
            keep_mt20_feat300 = sum( md$nFeature_RNA > 300 & md$pct.mt < 20 ),
            keep_mt10_feat400 = sum( md$nFeature_RNA > 400 & md$pct.mt < 10 ),
            keep_mt15_feat400 = sum( md$nFeature_RNA > 400 & md$pct.mt < 15 ),
            keep_mt20_feat400 = sum( md$nFeature_RNA > 400 & md$pct.mt < 20 ),
            keep_mt10_feat500 = sum( md$nFeature_RNA > 500 & md$pct.mt < 10 ),
            keep_mt15_feat500 = sum( md$nFeature_RNA > 500 & md$pct.mt < 15 ),
            keep_mt20_feat500 = sum( md$nFeature_RNA > 500 & md$pct.mt < 20 )
        )
    }) %>% bind_rows()
    sce_qc_simulation_pct <- sce_qc_simulation %>% mutate(
        across( starts_with("keep_"), ~ round(.x / total * 100 , 2))
    )
    openxlsx::write.xlsx(
        list(
            retain_cell = sce_qc_simulation,
            retain_pct  = sce_qc_simulation_pct
        ),
        sprintf("%s/scRNAseq_qc.simulation.results.xlsx", WorkDir)
    )

    ### !!!!
    ### sce_qc_simulation_pct 결과보고 적정 mt 및 feature threshold 설정
    ### 
    PCT_MT_CUTOFF <- 15     # <- 수동설정
    FEATURE_CUTFF <- 300    # <- 수동설정
    ###

    ### qc-pass subset
    sce_list_qc <- lapply(sce_list, function(sce) { subset( sce, subset = nFeature_RNA > FEATURE_CUTFF & pct.mt < PCT_MT_CUTOFF ) })
    ###

    # qc-summary (after filter)
    sce_qc_summary_after <- lapply(names(sce_list_qc), function(sid) {
        md <- sce_list_qc[[sid]]@meta.data
        data.frame(
            sample = sid,
            cells  = nrow(md),
            median_nFeature   = median(md$nFeature_RNA),
            median_nCount     = median(md$nCount_RNA),
            median_percent_mt = median(md$pct.mt)
        )
    }) %>% bind_rows()

    # qc-summary before/after table output
    openxlsx::write.xlsx(
        list(
            before = sce_qc_summary_before,
            after  = sce_qc_summary_after
        ),
        sprintf("%s/scRNAseq_qc.result.stats.xlsx", WorkDir)
    )
    # qc-plots (after filter) 
    p_qc.violin.after <- cowplot::plot_grid(
        plotlist = lapply( names(sce_list_qc), function(sid) {
            p.viloin <- VlnPlot( sce_list_qc[[sid]], layer="counts", features=c("nFeature_RNA","nCount_RNA","pct.mt"), ncol=3, pt.size=0 ) + 
                plot_annotation(title = sid) 
            return(p.viloin)
        }),
        nrow=1
    )
    p_feat.scatter.after <- cowplot::plot_grid(
        plotlist = lapply( names(sce_list_qc), function(sid) {
            p.feat.scatter <- FeatureScatter( sce_list_qc[[sid]], feature1="nCount_RNA", feature2="nFeature_RNA" ) + 
                ggtitle(sid) + theme(plot.title=element_text(face="bold", size=14, family="Roboto"))
            return(p.feat.scatter)
        }),
        nrow=1
    )
    p_qc.summary.after <- cowplot::plot_grid( p_qc.violin.after, p_feat.scatter.after, nrow=2 )
    ggsave( p_qc.summary.after, 
        file=sprintf("%s/scRNAseq_qc.summary.after.filtering.png", WorkDir),
        width=8*length(sce_list_qc), height=12, unit='in', dpi=150
    )
####

#### NORMALIZATION

    # 1. LogNormalize : 기본 정규화
    # 2. SCTransform  : depth/variance 보정 정규화 

    # run normalization 
    sce_list_norm <- lapply(sce_list_qc, function(sce) {
        # LogNormalize
        sce <- NormalizeData(sce)
        sce <- FindVariableFeatures(sce)
        sce <- ScaleData(sce)
        # SCTransform
        sce <- SCTransform( sce, vars.to.regress="pct.mt", verbose=TRUE )
        return(sce)
    })
    # save intermediate seurat object
    saveRDS( sce_list_norm, file=sprintf("%s/scRANseq_filter.pass.normalized.SCE.onject.list.RDS", WorkDir))
####

#### INTEGRATION

    # integration : feature selection
    # | nfeatures | 장점                   | 단점                            | 추천 상황                            |
    # |      2000 | 보수적, 빠름, noise 적음 | rare cell type signal 부족 가능 | 단순한 샘플, cell type 다양성 낮음     |
    # |      3000 | 균형 좋음               | 약간 더 많은 noise 포함 가능      | 일반적인 SCT integration 기본 후보    |
    # |      5000 | 다양한 signal 포착 가능  | batch/noise도 더 들어갈 수 있음  | cell type이 매우 다양하거나 atlas 수준 |

    # sample-name fix before integration
    sce_list_norm <- lapply(names(sce_list_norm), function(sid) {
        sce <- sce_list_norm[[sid]]
        sce <- RenameCells(sce, add.cell.id=sid)
        return(sce)
    })

    ### !!!!
    ###
    INTEGRATION_nFEATURES <- 3000 # <- 수동설정
    ###

    # run integration
    integrate_features <- SelectIntegrationFeatures(
        object.list = sce_list_norm, 
        nfeatures   = INTEGRATION_nFEATURES
    )
    sce_list_norm <- PrepSCTIntegration(
        object.list     = sce_list_norm, 
        anchor.features = integrate_features
    )
    integrate_anchors <- FindIntegrationAnchors(
        object.list          = sce_list_norm, 
        normalization.method = "SCT", 
        anchor.features      = integrate_features
    )
    SCE_DATA <- suppressWarnings(IntegrateData(
        anchorset            = integrate_anchors,
        normalization.method = "SCT"
    ))
    # save intermediate seurat object
    saveRDS( SCE_DATA, file=sprintf("%s/scRNAseq_integrated.SCE.object.RDS", WorkDir))
####

#### PCA & UMAP
    
    ## PCA ##
    SCE_DATA    <- RunPCA( SCE_DATA, assay="integrated" )
    # elbow plot
    p_pca.elbow <- ElbowPlot(SCE_DATA, ndims = 50)
    ggsave( p_pca.elbow, 
        file=sprintf("%s/scRNAseq_PCA.elbow.plot.png", WorkDir), 
        width=8, height=8, unit='in', dpi=150
    )

    ### !!! Elbow plot을 보고 umap에 사용할 dims 결정 
    ###
    UMAP_DIMS        <- 1:20                    # <- 수동설정
    RESOLUTION_RANGE <- c(0.2, 0.4, 0.6, 0.8)   # <- 수동설정
    ###

    ## UMAP ##
    # DefaultAssay(SCE_DATA) <- "integrated"
    SCE_DATA <- RunUMAP( SCE_DATA, assay="integrated", dims=UMAP_DIMS )
    SCE_DATA <- FindNeighbors( SCE_DATA, dims=UMAP_DIMS )
    SCE_DATA <- FindClusters( SCE_DATA, resolution=RESOLUTION_RANGE )

    # Clusters by Resolution ranges
    ResolutionList         <- paste0("integrated_snn_res.", RESOLUTION_RANGE)
    clusters_by_resolution <- sapply( ResolutionList, function(x) length(unique(SCE_DATA[[x]][, 1])) )
    
    # umap-plot by resolution
    p_umap.cluster.by.resolution <- cowplot::plot_grid(
        plotlist = lapply( ResolutionList, function(rs) { 
            p.umap <- DimPlot( SCE_DATA, reduction = "umap", group.by=rs, label=TRUE )
            return(p.umap)
        }),
        nrow=1
    )
    ggsave( p_umap.cluster.by.resolution, 
        file=sprintf("%s/scRNAseq_UMAP.clusters.by.resolution.png", WorkDir), 
        width=8*length(ResolutionList)*1.1, height=8, unit='in', dpi=150
    )

    ### !!!! Resolution 별로 생성한 UMAP plot 을 기반으로 Clustering Resolution 설정
    ####
    SELECTED_RESOLUTION <- "integrated_snn_res.0.4"  # <- 수동설정
    ####

    # apply selected cluster-resolution
    Idents(SCE_DATA) <- SELECTED_RESOLUTION
    SCE_DATA$seurat_clusters <- Idents(SCE_DATA)
    
    # remove test resolution-data (columns)
    for( RS in ResolutionList ){ SCE_DATA[[ RS ]] <- NULL }
    
    # umap-plot of selected cluster-resolution
    p_umap.by.cluster <- DimPlot( SCE_DATA, reduction="umap", group.by="seurat_clusters", label=TRUE )
    p_umap.by.sample  <- DimPlot( SCE_DATA, reduction="umap", group.by="sample" )
    p_umap <- cowplot::plot_grid( p_umap.by.cluster, p_umap.by.sample, nrow=1 )   
    ggsave( p_umap, file=sprintf("%s/scRNAseq_UMAP.selected.clusters.png", WorkDir), 
        width=16*1.1, height=8, unit='in', dpi=150
    )

    #
    sample_cluster_table <- table( SCE_DATA$sample, SCE_DATA$seurat_clusters )
    umap_cluster_summary <- list(
        cluster_cells = sample_cluster_table,
        cluster_ratio = round( prop.table(sample_cluster_table, margin = 1) * 100, 1 ),
        sample_ratio  = round( prop.table(sample_cluster_table, margin = 2) * 100, 1 )
    )

    openxlsx::write.xlsx( umap_cluster_summary, sprintf("%s/scRNAseq_UMAP.clusters.metric.xlsx", WorkDir) )

####

#### SingleR Annotation

    # SingleR References : Main/Immune-Focused 2종류 동시적용 기본으로함
    if( DataTax == "human"){
        ref_main   <- HumanPrimaryCellAtlasData()
        ref_immune <- MonacoImmuneData()
    }else if( DataTax == "mouse" ){
        ref_main   <- MouseRNAseqData()
        ref_immune <- ImmGenData()
    }else{
        stop("No Taxonomy info.", call.=FALSE)
    } 
    singler_ref_list       <- list(SingleR_main=ref_main, SingleR_immune=ref_immune)
    singler_ref_list_label <- lapply(singler_ref_list, function(x) x$label.main)

    # sce copy for singler annotation
    sce_query <- suppressMessages(as.SingleCellExperiment( SCE_DATA, assay="SCT" ))

    # cell-level
    message(sprintf("[ %s ] annotate cells using SingleR...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_singler_cell_annot <- SingleR(test=sce_query, ref=singler_ref_list, labels=singler_ref_list_label )
    # cluster-level 
    message(sprintf("[ %s ] annotate clusters using SingleR...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_singler_cluster_annot <- SingleR(test=sce_query, ref=singler_ref_list, labels=singler_ref_list_label, clusters=sce_query$seurat_clusters )

    sce_singler_cluster_annot_res <- data.frame(
        cluster              = as.character(rownames(sce_singler_cluster_annot)), 
        SingleR_cluster      = sce_singler_cluster_annot$labels,
        SingleR_cluster_conf = sce_singler_cluster_annot$pruned.labels,
        SingleR_cluster_ref  = sce_singler_cluster_annot$reference
    )

    # add SingleR annotation to main SCE-object
    SCE_DATA$SingleR_cell         <- sce_singler_cell_annot$labels
    SCE_DATA$SingleR_cell_ref     <- sce_singler_cell_annot$reference
    SCE_DATA$SingleR_cluster      <- sce_singler_cluster_annot_res[match(SCE_DATA$seurat_clusters, sce_singler_cluster_annot_res$cluster), "SingleR_cluster"     ]
    SCE_DATA$SingleR_cluster_conf <- sce_singler_cluster_annot_res[match(SCE_DATA$seurat_clusters, sce_singler_cluster_annot_res$cluster), "SingleR_cluster_conf"]
    SCE_DATA$SingleR_cluster_ref  <- sce_singler_cluster_annot_res[match(SCE_DATA$seurat_clusters, sce_singler_cluster_annot_res$cluster), "SingleR_cluster_ref" ]

    # SingleR annotation result summary table list
    sce_singler_annot_res_cells <- table( cluster = SCE_DATA$seurat_clusters, SingleR = SCE_DATA$SingleR_cell, useNA = "ifany" )
    sce_singler_annot_res_prop  <- round( prop.table(sce_singler_annot_res_cells, margin = 1) * 100, 1 )

    sce_singler_annot_res_cluster_cells <- sce_singler_cluster_annot_res %>% left_join(
        data.frame( cluster=rownames(sce_singler_annot_res_cells), as.data.frame.matrix(singler_annot_res_cells), check.names = FALSE ),
        by="cluster"
    )
    sce_singler_annot_res_cluster_prop  <- sce_singler_cluster_annot_res %>% left_join(
        data.frame( cluster = rownames(sce_singler_annot_res_prop), as.data.frame.matrix(singler_annot_res_prop), check.names = FALSE ),
        by="cluster"
    )
    sce_singler_annot_summary <- list(
        Cell_annot    = sce_singler_annot_res_cells, Cell_props = sce_singler_annot_res_cluster_cells,
        Cluster_annot = sce_singler_annot_res_cluster_cells, Cluster_props = sce_singler_annot_res_cluster_prop
    )
    openxlsx::write.xlsx( sce_singler_annot_summary, sprintf("%s/scRNAseq_SingleR.annotation.summary.xlsx", WorkDir) )

####

#### CELL-MARKER-2.0 ANNOTATION 
    
    #### 
    SELECTED_CELLMARKER_TISSUE <- c("Lung","Liver","Spleen")
    ####
    
    message(sprintf("[ %s ] cell-marker selected tissue : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(SELECTED_CELLMARKER_TISSUE, collapse=", ")))
    message(sprintf("[ %s ] filter CellMarker by selected tissues....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cellmarker_data <- readRDS(sprintf("/storage/home/kangsm/myDB/rds/CellMarker2.0_curated.%s.single.cell.markers.RDS", DataTax))
    # if( DataTax == "mouse" ){
    #     cellmarker_data <- readRDS("/storage/home/kangsm/myDB/rds/CellMarker2.0_curated.mouse.single.cell.markers.RDS")
    # }else if( DataTax == "human" ){
    #     cellmarker_data <- readRDS("/storage/home/kangsm/myDB/rds/CellMarker2.0_curated.human.single.cell.markers.RDS")
    # }
    message(sprintf("[ %s ] annotate clusters using CellMarker...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))   
    cellmarker_set <- FilterCellMarkerSet(cellmarker_data=cellmarker_data, genespace=rownames(SCE_DATA), selected_tissues=SELECTED_CELLMARKER_TISSUE)
    # annotation
    SCE_DATA <- AnnotateCellMarkerSet(sce=SCE_DATA, cellmarker_set=cellmarker_set, is_seurat=TRUE )
####

#### UMAP PLOT BY ANNOTATYED CELL TYPES
    message(sprintf("[ %s ] create UMAP plots with annotations....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    umap_annot_plot <- suppressWarnings(cowplot::plot_grid( 
        plotlist=list(
            DimPlot( SCE_DATA, reduction="umap", group.by="SingleR_cluster_conf", label=TRUE, repel=TRUE, label.size=6 ) + ggtitle("SingleR") + 
                theme(
                    legend.text  = element_text(size = 14, family="Arial"),
                    legend.title = element_text(size = 15, family="Arial"),
                    plot.title  = element_text(face="bold", size=18, family="Roboto")
                ),
            DimPlot( SCE_DATA, reduction="umap", group.by="CellMarkerAnnot", label=TRUE, repel=TRUE, label.size=6 ) + ggtitle("CellMarker 2.0") +
                theme(
                    legend.text  = element_text(size = 14, family="Arial"),
                    legend.title = element_text(size = 15, family="Arial"),
                    plot.title  = element_text(face="bold", size=18, family="Roboto")
                )
        ), 
        nrow=1
    ))
    ggsave(umap_annot_plot, file=sprintf("%s/scRNAseq_anootated.umap.plot.png", WorkDir), width=20, height=10, unit='in', dpi=150  )
####

#### ANNOTATION SUMMARY TABLE
    message(sprintf("[ %s ] integrate and summarize annotation results....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    integrated_annot_summary <- SCE_DATA@meta.data %>% dplyr::select(
        c("cluster"="seurat_clusters","singler_cell"="SingleR_cluster_conf", "cellmarker_tissue"="CellMarker_tissue", "cell_marker_cell"="CellMarker_cell")
    ) %>% unique() %>% arrange(as.numeric(cluster))
    rownames(integrated_annot_summary) <- NULL
    write.table(integrated_annot_summary, sprintf("%s/scRNAseq_cluster.annotation.summary.tsv", WorkDir), quote=F, col.names=T, row.names=F, sep="\t")
####

#### CLUSTER MARKERS
    message(sprintf("[ %s ] extract cluster markers....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    Idents(SCE_DATA) <- "seurat_clusters"
    SCE_DATA <- PrepSCTFindMarkers(SCE_DATA, assay="SCT")
    message(sprintf("[ %s ] positive-only markers....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_cluster_markers_pos <- FindAllMarkers(
        object          = SCE_DATA,
        assay           = "SCT",
        slot            = "data",
        only.pos        = TRUE,
        min.pct         = 0.2,
        logfc.threshold = 0.25,
        test.use        = "wilcox"
    )
    message(sprintf("[ %s ] all positive & negative markers....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_cluster_markers_all <- FindAllMarkers(
        object          = SCE_DATA,
        assay           = "SCT",
        slot            = "data",
        only.pos        = FALSE,
        min.pct         = 0.2,
        logfc.threshold = 0.25,
        test.use        = "wilcox"
    )
    write.table( sce_cluster_markers_pos, sprintf("%s/scRNAseq_cluster.markers.positive.only.tsv", WorkDir), quote=F, col.names=T, row.names=F, sep="\t" )
    write.table( sce_cluster_markers_all, sprintf("%s/scRNAseq_cluster.markers.all.tsv",           WorkDir), quote=F, col.names=T, row.names=F, sep="\t" )
####

#### SAVE ANALYSIS-FINISHED SCE-OBJECTS
    saveRDS(SCE_DATA, file=sprintf("%s/scRNAseq_analyzed.SCE.objects.RDS", WorkDir))
####




