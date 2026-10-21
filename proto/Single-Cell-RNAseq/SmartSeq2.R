
#### packages #
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
        library("scater")
        library("scran")
        library("scuttle")
        library("igraph")
        library("SingleR")
        library("celldex")
        library("rtracklayer")
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
####
####
    setwd("/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/analysis/SmartSeq2")
    BaseDir <- "/data/scRNAseq/02.Public-Datasets"
    BatchID <- "004.GSE151334.PRJNA635600.cell-line-smartseq2"
    WorkDir <- "/data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/analysis/SmartSeq2"
    DataTax        <- "human"
    GenomeAssembly <- "hg38_ercc.plus"
    #
    sinfo <- fread(sprintf("%s/%s/meta/metadata_SmartSeq2.tsv", BaseDir, BatchID))
    
####
#### 
    gtf_file   <- sprintf("/storage/references_and_index/RNAseq/%s/genes/genes.gtf", GenomeAssembly)
    gene_annot <- AnnotDataFromGTF(gtf_file=gtf_file, tax=DataTax)
####
#### featureCount matrix
    count_file_check <- mclapply( indexing(sinfo$SeqID), function(SeqID) {
        data_file <- sprintf("%s/%s/SmartSeq2/%s/quant/featureCounts/featureCounts.tsv", BaseDir, BatchID, SeqID)
        check_res <- ifelse( file.exists(data_file), TRUE, FALSE)
        return(check_res)
    }, mc.cores=16)
    # >>>   table(unlist(count_file_check))

    count_list <- mclapply(indexing(sinfo$SeqID), function(SeqID) {
        data_dir <- sprintf("%s/%s/SmartSeq2/%s/quant/featureCounts", BaseDir, BatchID, SeqID)        
        cnt <- fread(sprintf("%s/featureCounts.tsv", data_dir), skip=1)
        out <- cnt[, c(1,6,7)]
        colnames(out)[3] <- SeqID
        return(out)
    }, mc.cores=16)

    count_mat <- SmartSeq.CountListToMatrix(smartseq_count_list=count_list)
####

#### Initial Basic QC
    qc_smartseq2 <- SmartSeq.BasicQC( count_mat=count_mat, gene_annot=gene_annot, include_ercc=TRUE, draw_plots=TRUE, out_dir=WorkDir ) 
    
    # 초기 QC 상태를 보고 필터링 컷오프 설정
    # ERCC의 경우 필요하면 QC에 포함 (대게 10~20% 이상이면 제거)

    ###
    THRESHOLD_GENE  <- 1000     # 300, 500
    THRESHOLD_COUNT <- 30000    # 10000 
    THRESHOLD_MT    <- 20       # 15, 20, 25 등 <<-- MAJOR FILTERING FACTOR
    THRESHOLD_RIBO  <- 25       # 15, 20, 25 등 
    FILTER_BY_ERCC  <- TRUE     #
    THRESHOLD_ERCC  <- 20       # 10, 20 등 (대게 20 이상 제거)
    ###

    # qc simulation
    if( !FILTER_BY_ERCC ){
        qc_simulation_res <- data.frame(
            gene = THRESHOLD_GENE, count = THRESHOLD_COUNT, mt = THRESHOLD_MT, ercc = "not_used",
            total_cell  = nrow(qc_smartseq2),
            retain_cell = qc_smartseq2 %>% filter( 
                detected_genes >= THRESHOLD_GENE, total_counts >= THRESHOLD_COUNT, pct.mt <= THRESHOLD_MT,
                pct.ribo <= THRESHOLD_RIBO
            ) %>% nrow(.)
        )
    }else{
        qc_simulation_res <- data.frame(
            gene = THRESHOLD_GENE, count = THRESHOLD_COUNT, mt = THRESHOLD_MT, ercc = THRESHOLD_ERCC,
            total_cell  = nrow(qc_smartseq2),
            retain_cell = qc_smartseq2 %>% filter( 
                detected_genes >= THRESHOLD_GENE, total_counts >= THRESHOLD_COUNT, pct.mt <= THRESHOLD_MT,
                pct.ribo <= THRESHOLD_RIBO, pct.ercc <= THRESHOLD_ERCC 
            ) %>% nrow(.)            
        )
    }
    qc_simulation_res <- qc_simulation_res %>% mutate( pct.retain = round(retain_cell/total_cell * 100, 3) )
    print(qc_simulation_res)
    
    ## 여러 필터링 조건으로 시뮬레이션 후 결정
    #write.table( qc_simulation_res, sprintf("%s/smartseq2_qc.filter.simulation.tsv", WorkDir), quote=F, row.names=F, col.names=T, sep="\t")
    
    # qc-pass-cells status tags
    qc_smartseq2$qc_pass <- "pass"
    qc_smartseq2$qc_pass[qc_smartseq2$detected_genes < THRESHOLD_GENE ] <- "fail"
    qc_smartseq2$qc_pass[qc_smartseq2$total_counts   < THRESHOLD_COUNT] <- "fail"
    qc_smartseq2$qc_pass[qc_smartseq2$pct.mt         > THRESHOLD_MT   ] <- "fail"
    if( FILTER_BY_ERCC ) { qc_smartseq2$qc_pass[qc_smartseq2$pct.ercc > THRESHOLD_ERCC ] <- "fail" }

    # perform low-qc cells filtering
    keep_cells   <- qc_smartseq2$cell_id[qc_smartseq2$qc_pass == "pass"]
    count_mat_qc <- count_mat[, keep_cells, drop = FALSE]
    qc_smartseq2_pass <- qc_smartseq2[qc_smartseq2$qc_pass == "pass", ]

    # qc-pass cells plots
    png(sprintf("%s/smartseq2_qc.pass.png", WorkDir), width=21, height=14, units="in", res=150 )
    par(mfrow = c(2,3))
    hist( qc_smartseq2_pass$log10_total_counts,   breaks=50, main="Library size",    xlab="log10(total counts + 1)"   )
    hist( qc_smartseq2_pass$log10_detected_genes, breaks=50, main="Detected genes",  xlab="log10(detected genes + 1)" )
    plot( x=qc_smartseq2_pass$log10_total_counts, y=qc_smartseq2_pass$log10_detected_genes,
        pch=19, cex=1, col = ifelse(qc_smartseq2$pct.mt > 20, "#E41A1C", "#7a7a7a"),
        xlab = "log10(total counts + 1)", ylab = "log10(detected genes + 1)", main = "SmartSeq2 Basic QC"
    )
    legend( "bottomright", legend=c("pct.mt <= 20", "pct.mt > 20"), col=c("#7a7a7a", "#E41A1C"), pch=19, bty="n" )
    hist( qc_smartseq2_pass$pct.mt,   breaks=50, main="Mitochondrial %", xlab="pct.mt"   )
    hist( qc_smartseq2_pass$pct.ribo, breaks=50, main="Ribosomal %",     xlab="pct.ribo" )
    if( FILTER_BY_ERCC ){ hist( qc_smartseq2_pass$pct.ercc, breaks=50, main="ERCC-spikein %", xlab="pct.ercc" ) }
    dev.off()

    write.table( qc_smartseq2_pass, sprintf("%s/smartseq2_qc.pass.sample.info.tsv", WorkDir), quote=F, row.names=F, col.names=T, sep="\t")

####

#### create SCE-object + scran-normalization

    # sce object
    sce <- SingleCellExperiment( 
        assays  = list(counts = count_mat_qc), 
        colData = DataFrame( qc_smartseq2_pass[colnames(count_mat_qc), ] )
    )
    rowData(sce)$gene_id   <- rownames(sce)
    rowData(sce)$gene      <- gene_annot[match(rownames(sce), gene_annot$gene_id), "gene"]
    rowData(sce)$gene_type <- gene_annot[match(rownames(sce), gene_annot$gene_id), "gene_type"]

    # ERCC scaling
    if( FILTER_BY_ERCC )
    {
        is_ercc <- grepl("^ERCC", rowData(sce)$gene_id)
        sce     <- splitAltExps( sce, ifelse(is_ercc, "ERCC", "endogenous") )
        sce     <- computeSpikeFactors( sce, assay.type = "counts", spikes = "ERCC" )
    }

    # Low-expressed genes filter
    # 필터링 시뮬레이션 후 적정선 설정 
    ### 
    MIN_CELLS <- c(1, 2, 3, 5, 10, 15, 20)
    ###

    gene_filter_simulation <- data.frame(
        min_cells      = MIN_CELLS,
        retained_genes = sapply(MIN_CELLS, function(x) { sum(rowSums(counts(sce) > 0) >= x) })
    ) %>% mutate( 
        pct.retained = round(retained_genes / nrow(counts(sce)) * 100, 3)
    )   
    print(gene_filter_simulation)

    ###
    LOW_EXPRESS_GENE_CELL_CUTOFF <- 3
    ###

    keep_genes <- rowSums(counts(sce) > 0) >= LOW_EXPRESS_GENE_CELL_CUTOFF
    sce        <- sce[keep_genes, ]
    count_mat_filtered <- counts(sce)


    # endogenous scaling & normalization
    set.seed(1)
    sce_quick_clust <- quickCluster(sce)
    sce <- computeSumFactors(sce, clusters=sce_quick_clust)
    sce <- logNormCounts(sce)
    
    qc_smartseq2_pass$size_factor <- round(sizeFactors(sce), 4)
    write.table( qc_smartseq2_pass, sprintf("%s/smartseq2_qc.pass.sample.info.tsv", WorkDir), quote=F, row.names=F, col.names=T, sep="\t")

    png(sprintf("%s/smartseq2_scran.normlization.size.factor.png", WorkDir), width=14, height=7, units="in", res=150 )
    par(mfrow = c(1,2))
    plot( x=qc_smartseq2_pass$total_counts, y=qc_smartseq2_pass$size_factor,
        log="xy", pch=19, cex=1, 
        xlab="Total counts", ylab="scran size factor", main="Size factor vs total counts"
    )
    plot( x=qc_smartseq2_pass$detected_genes, y=qc_smartseq2_pass$size_factor,
        log="y", pch=19, cex=1, 
        xlab="Detected genes", ylab="scran size factor", main="Size factor vs detected genes"
    )   
    dev.off()

    ### SAVE PROCESSED DATA 
    saveRDS( qc_smartseq2_pass,  file=sprintf("%s/smartseq2_qc.pass.sample.info.RDS",        WorkDir))
    saveRDS( count_mat_filtered, file=sprintf("%s/smartseq2_qc.pass.raw.count.matrix.RDS",   WorkDir))
    saveRDS( gene_annot,         file=sprintf("%s/smartseq2_gene.annotation.RDS",            WorkDir))
    saveRDS( sce,                file=sprintf("%s/smartseq2_qc.pass.normalized.SCE.obj.RDS", WorkDir))

####

#### PCA
    set.seed(1)
    sce_vm  <- modelGeneVar(sce)
    sce_hvg <- getTopHVGs( sce_vm, n=3000 )
    saveRDS( sce_vm,  file=sprintf("%s/smartseq2_model.gene.var.RDS", WorkDir))
    saveRDS( sce_hvg, file=sprintf("%s/smartseq2_hvg.3000.RDS",        WorkDir))

    set.seed(1)
    # run PCA
    sce <- runPCA( sce, subset_row=sce_hvg, exprs_values="logcounts", ncomponents=50, name="pca" )
    
    # elbow plot 확인 및 적정 PCA comp cutoff 설정
    percent_var <- attr(reducedDim(sce, "pca"), "percentVar")

    png(sprintf("%s/smartseq2_pca.plot.png", WorkDir), width=14, height=7, units="in", res=150 )
    par(mfrow = c(1,2))
    plot( percent_var[1:50], type="b", pch=19, 
        xlab="PC", ylab="Variance explained (%)", main="SmartSeq2 PCA elbow plot"
    )
    plot( x=reducedDim(sce, "pca")[, 1], y=reducedDim(sce, "pca")[, 2],
        pch=19, cex=1, col = ifelse(sce$pct.mt > 10, "#E41A1C", "#7a7a7a"),
        xlab="PC1", ylab="PC2", main="SmartSeq2 PCA"
    )
    dev.off()

    ### elbow plot 보고 매뉴얼 결정
    PC_N   <- 8   # elbow plot 에서 정한 PC
    SNN_K  <- 50  # 
    UMAP_K <- 50  # neighbor 수
    ###
####

#### UMAP & Clustering
    # clustering 용 객체 생성
    pcadim <- reducedDim(sce, "pca")[, 1:PC_N]
    reducedDim(sce, "pcadim") <- pcadim
    for( snn_k in c(20,50,70,100) ){
        for( res in seq(0.1,1,0.1) ){
            test_snn    <- bluster::makeSNNGraph( pcadim, k=snn_k )
            test_leiden <- igraph::cluster_leiden(test_snn, objective_function="modularity", resolution=res, n_iterations=10 )
            memb <- factor(igraph::membership(test_leiden))
            cat( "snn_k =", snn_k, "resolution =", res, "clusters =", nlevels(memb), "sizes =", paste(sort(table(memb)), collapse = ", "), "\n" )
        }
    }

    ###
    SNN_K              <- 70
    CLUSTER_RESOLUTION <- 0.2
    UMAP_K             <- 348 
    ###

    graph_snn      <- bluster::makeSNNGraph( pcadim, k=SNN_K )
    cluster_leiden <- igraph::cluster_leiden(graph_snn, objective_function="modularity", resolution=CLUSTER_RESOLUTION, n_iterations=10 )
    sce$cluster    <- factor(igraph::membership(cluster_leiden))
    colLabels(sce) <- sce$cluster

    #cluster_walktrap <- igraph::cluster_walktrap(graph_snn)
    #sce$cluster      <- factor( igraph::membership(cluster_walktrap) )
    #colLabels(sce)   <- sce$cluster
    
    # run umap
    sce <- scater::runUMAP( sce, dimred="pcadim", n_neighbors=UMAP_K, name="umap" )

    png(sprintf("%s/smartseq2_umap.plot.png", WorkDir), width=7, height=7, units="in", res=150 )
    par(mfrow = c(1,1))
    plot( x=reducedDim(sce, "umap")[, 1], y=reducedDim(sce, "umap")[, 2],
        pch=19, cex=0.8, col=as.integer(sce$cluster),
        xlab="UMAP1", ylab="UMAP2", main="SmartSeq2 UMAP Clusters"
    )
    #legend( "topright", legend=levels(sce$cluster), col=seq_along(levels(sce$cluster)), pch=19, bty = "n" )
    dev.off()

    cluster_qc <- aggregate( 
        cbind(total_counts, detected_genes, pct.mt, pct.ribo) ~ cluster, data = as.data.frame(colData(sce)), FUN = median 
    ) %>% left_join(
        data.frame( cluster = names(table(sce$cluster)), n_cells = as.integer(table(sce$cluster)) ) %>% mutate(
            pct_cluster = round(n_cells/ncol(sce) * 100, 3)
        ),
        by="cluster"
    )
    write.table( cluster_qc, sprintf("%s/smartseq2_umap.cluster.metric.tsv", WorkDir), quote=F, row.names=F, col.names=T, sep="\t")
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
    singler_ref_list <- list(SingleR_main=ref_main, SingleR_immune=ref_immune)
    singler_ref_list_label <- lapply(singler_ref_list, function(x) x$label.main)

    # SingleR 용 sce 따로 생성 (gene-level aggregation 때문)
    snglr_logcnt <- as.matrix(logcounts(sce))
    snglr_cnt    <- as.matrix(counts(sce))
    snglr_gene   <- rowData(sce)$gene 
    
    message(sprintf("[ %s ] aggregate log-counts by gene-level....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    logcnt_df     <- data.frame(gene=snglr_gene, snglr_logcnt, check.names=FALSE ) %>% data.table()
    logcnt_df_rev <- logcnt_df[, lapply(.SD, mean, na.rm=TRUE), by=gene] %>% as.data.frame()
    rownames(logcnt_df_rev) <- logcnt_df_rev$gene
    #logcnt_df_rev <- aggregate(. ~ gene, data=logcnt_df, FUN=mean )

    message(sprintf("[ %s ] aggregate counts by gene-level....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cnt_df     <- data.frame(gene=snglr_gene, snglr_cnt, check.names=FALSE ) %>% data.table()
    cnt_df_rev <- cnt_df[, lapply(.SD, mean, na.rm=TRUE), by=gene] %>% as.data.frame()
    rownames(cnt_df_rev) <- cnt_df_rev$gene   
    #cnt_df_rev <- aggregate(. ~ gene, data=cnt_df, FUN=sum )

    message(sprintf("[ %s ] create shadow sce-object for SingleR annotation...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_singler <- SingleCellExperiment(
        assays  = list( counts=as.matrix(cnt_df_rev[,-1]), logcounts=as.matrix(logcnt_df_rev[,-1]) ),
        colData = colData(sce)
    )

    # SingleR annotation : cell-level, cluster-level 모두 수행
    # cell-level
    message(sprintf("[ %s ] annotate cells using SingleR...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_singler_cell_annot <- SingleR(test=sce_singler, ref=singler_ref_list, labels=singler_ref_list_label,
        assay.type.test = "logcounts", assay.type.ref  = "logcounts"
    )
    # cluster-level 
    message(sprintf("[ %s ] annotate clusters using SingleR...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    sce_singler_cluster_annot <- SingleR(test=sce_singler, ref=singler_ref_list, labels=singler_ref_list_label, clusters=sce$cluster,
        assay.type.test = "logcounts", assay.type.ref  = "logcounts"
    )

    sce_singler_cluster_annot_res <- data.frame(
        cluster              = as.character(rownames(sce_singler_cluster_annot)), 
        SingleR_cluster      = sce_singler_cluster_annot$labels,
        SingleR_cluster_conf = sce_singler_cluster_annot$pruned.labels,
        SingleR_cluster_ref  = sce_singler_cluster_annot$reference
    )
    # add SingleR annotation to main SCE-object
    sce$SingleR_cell     <- sce_singler_cell_annot$labels
    sce$SingleR_cell_ref <- sce_singler_cell_annot$reference
    sce$SingleR_cluster      <- sce_singler_cluster_annot_res[match(sce$cluster, sce_singler_cluster_annot_res$cluster), "SingleR_cluster"     ]
    sce$SingleR_cluster_conf <- sce_singler_cluster_annot_res[match(sce$cluster, sce_singler_cluster_annot_res$cluster), "SingleR_cluster_conf"]
    sce$SingleR_cluster_ref  <- sce_singler_cluster_annot_res[match(sce$cluster, sce_singler_cluster_annot_res$cluster), "SingleR_cluster_ref" ]

    # SingleR annotation result summary table list
    singler_annot_res_cells <- table( cluster = sce$cluster, SingleR = sce$SingleR_cell, useNA = "ifany" )
    singler_annot_res_prop  <- round( prop.table(singler_annot_res_cells, margin = 1) * 100, 1 )

    singler_annot_res_cluster_cells <- sce_singler_cluster_annot_res %>% left_join(
        data.frame( cluster = rownames(singler_annot_res_cells), as.data.frame.matrix(singler_annot_res_cells), check.names = FALSE ),
        by="cluster"
    )
    singler_annot_res_cluster_prop  <- sce_singler_cluster_annot_res %>% left_join(
        data.frame( cluster = rownames(singler_annot_res_prop), as.data.frame.matrix(singler_annot_res_prop), check.names = FALSE ),
        by="cluster"
    )
    singler_annot_summary <- list(
        Cell_annot = singler_annot_res_cells, Cell_props = singler_annot_res_prop,
        Cluster_annot = singler_annot_res_cluster_cells, Cluster_props = singler_annot_res_cluster_prop
    )
    openxlsx::write.xlsx( singler_annot_summary, sprintf("%s/smartseq2_SingleR.annotation.summary.xlsx", WorkDir) )

####

#### CELL-MARKER-2.0 ANNOTATION 
    
    # CELL-MARKER 2.0 : cluster annotation 수행하며 tissue 를 지정해서 실행 (샘플소스 tissue로 지정)
    #### 
    SELECTED_CELLMARKER_TISSUE <- c("Kidney","Breast","Skin")
    ####
    # cellmarker_data <- readRDS(sprintf("/storage/home/kangsm/myDB/rds/CellMarker2.0_curated.%s.single.cell.markers.RDS", DataTax))
    
    
    if( DataTax == "mouse" ){
        cellmarker_data <- readRDS("/storage/home/kangsm/myDB/rds/CellMarker2.0_curated.mouse.single.cell.markers.RDS")
    }else if( DataTax == "human" ){
        cellmarker_data <- readRDS("/storage/home/kangsm/myDB/rds/CellMarker2.0_curated.human.single.cell.markers.RDS")
        cellmarker_data <- readRDS("/storage/home/kangsm/myDB/rds/CellMarker2.0_curated.human.markers.RDS")
    }
       
    cellmarker_set <- FilterCellMarkerSet(
        cellmarker_data = cellmarker_data, 
        genespace       = rowData(sce)$gene, 
        min_marker_n    = 5,
        selected_tissues=SELECTED_CELLMARKER_TISSUE        
    )
    # annotation
    sce <- AnnotateCellMarkerSet(sce=sce, cellmarker_set=cellmarker_set, is_seurat=FALSE)
####

#### UMAP PLOT BY ANNOTATYED CELL TYPES

    umap_annot_plot <- cowplot::plot_grid( 
        plotlist=list(
            scater::plotReducedDim(sce, dimred = "umap", colour_by = "SingleR_cluster_conf", point_size=1 ) + ggtitle("SingleR") + 
                guides(colour=guide_legend(title = "SingleR Annot")) +
                theme(
                    legend.text  = element_text(size = 14, family="Arial"),
                    legend.title = element_text(size = 15, family="Arial"),
                    plot.title  = element_text(face="bold", size=18, family="Roboto")
                ),
            scater::plotReducedDim(sce, dimred = "umap", colour_by = "CellMarkerAnnot", point_size=1 ) + ggtitle("CellMarker 2.0") +
                guides(colour=guide_legend(title = "CellMarker Annot")) +
                theme(
                    legend.text  = element_text(size = 14, family="Arial"),
                    legend.title = element_text(size = 15, family="Arial"),
                    plot.title  = element_text(face="bold", size=18, family="Roboto")
                )
        ), 
        nrow=1, rel_widths=c(1,2)
    )
    ggsave(umap_annot_plot, file=sprintf("%s/smartseq2_anootated.umap.plot.png", WorkDir), width=20, height=10, unit='in', dpi=150  )

####

#### CLUSTER MARKERS

    sce_cluster_markers <- findMarkers( sce, groups=sce$cluster, assay.type="logcounts", direction="up")
    marker_gene_annot   <- gene_annot %>% dplyr::select(c("gene_id","gene","gene_type"))

    sce_cluster_top_marker_list <- lapply(names(sce_cluster_markers), function(cl) {
        tab <- as.data.frame(sce_cluster_markers[[cl]])
        tab$gene_id <- rownames(tab)
        tab$cluster <- cl
        tab <- tab[order(tab$Top), ]
        selected_top <- tab$Top[50]

        topMarker <- tab %>% filter( Top <= selected_top ) %>% left_join( marker_gene_annot, by="gene_id") %>% dplyr::select(
            c("cluster","gene","gene_id","Top","p.value", "FDR", "summary.logFC","gene_type")
        ) %>% unique()
        return(topMarker)
    }) 
    saveRDS(sce_cluster_top_marker_list, file=sprintf("%s/smartseq2_cluster.markers.RDS", WorkDir))
   
####

#### 
    sinfo2 = as.data.frame(sinfo)
    colData(sce)$cell <- as.character(sinfo2[match(colData(sce)$cell_id, sinfo2$SeqID), "Cell"])
    saveRDS(sce, file=sprintf("%s/smartseq2_analyzed.SCE.objects.RDS", WorkDir))

    pseudobulk_sce <- scuttle::aggregateAcrossCells( sce, ids=colData(sce)$cell, use.assay.type="counts" )
    pseudobulk_counts <- assay(pseudobulk_sce, "counts")

    pseudobulk_gene_symbol <- rowData(pseudobulk_sce)$gene
    names(pseudobulk_gene_symbol) <- rownames(pseudobulk_sce)

    pseudobulk_counts_df   <- data.frame(gene = pseudobulk_gene_symbol, as.matrix(pseudobulk_counts), check.names = FALSE )
    pseudobulk_counts_df   <- pseudobulk_counts_df[!is.na(pseudobulk_counts_df$gene) & pseudobulk_counts_df$gene != "", ] %>% data.table()
    pseudobulk_counts_gene <- pseudobulk_counts_df[, lapply(.SD, mean, na.rm=TRUE), by=gene] %>% as.data.frame() 
    rownames(pseudobulk_counts_gene) <- pseudobulk_counts_gene$gene
    pseudobulk_counts_mat <- as.matrix(pseudobulk_counts_gene[, -1])

    saveRDS(pseudobulk_sce, file=sprintf("%s/smartseq2_pseudo.bulk.counts.SCE.objects.RDS", WorkDir))
    saveRDS(pseudobulk_counts_gene, file=sprintf("%s/smartseq2_pseudo.bulk.counts.df.RDS", WorkDir))
    saveRDS(pseudobulk_counts_mat,  file=sprintf("%s/smartseq2_pseudo.bulk.counts.matrix.RDS", WorkDir))

####

#### Variant Calling 전처리용 samplelist (optional)
    # 들어가야 하는 필드
    # 1. individual id    <- RGSM
    # 2. seq-id           <- RGID & RGLB
    # 3. GenomeAssembly  
    # 4. RGPL
    # 5. RGPU
    # 6. RGCN              
    # 7.... tissue/cell type <- 보조 데이터 추가 얼마든지 (6번까지가 필수)
    var_call_sinfo <- data.frame(
        IndivID = sinfo$Cell,
        SeqID   = sinfo$SeqID,
        GenomeAssembly = "hg38_ercc.plus",
        RGPL = "ILLUMINA",
        RGPU = "SMARTSEQ2",
        RGCN = "GENCURIX"
    )

    write.table(var_call_sinfo,
        sprintf("%s/%s/meta/variant.call.SmartSeq2.sample.info", BaseDir, BatchID),
        quote=F, col.names=T, row.names=F, sep=" "
    )

    for( INDIV in unique(var_call_sinfo$IndivID) ){
        var_call_sinfo_indiv <- var_call_sinfo %>% filter( IndivID == INDIV )
        write.table(var_call_sinfo_indiv[,c("IndivID","SeqID")],
            sprintf("%s/%s/meta/variant.call.SmartSeq2.sample.%s.list", BaseDir, BatchID, INDIV),
            quote=F, col.names=T, row.names=F, sep=" "
        )
    }
    write.table(
        var_call_sinfo[,c("IndivID","GenomeAssembly")] %>% unique(),
        sprintf("%s/%s/meta/variant.call.SmartSeq2.sample.list", BaseDir, BatchID),
            quote=F, col.names=T, row.names=F, sep=" "
    )


