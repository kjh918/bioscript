

#### Gene Annotation Table From GTF
    AnnotDataFromGTF <- function( gtf_file, tax )
    {
        message(sprintf("[ %s ] read GTF....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
        gtf <- rtracklayer::import(gtf_file) %>% as.data.frame() %>% data.table()
        gtf_gene <- gtf[type == "gene"]
        message(sprintf("[ %s ] convert into annotation table....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
        if( tax == "mouse" ){
            gene_annot <- gtf_gene %>% dplyr::select(c(
                "gene_id","gene_name","gene_type","mgi_id","chrom"="seqnames"
            )) %>% mutate(chrom=as.character(chrom))
        }else if( tax == "human" ){
            gene_annot <- gtf_gene %>% dplyr::select(c(
                "gene_id","gene_name","gene_type","hgnc_id","chrom"="seqnames"
            )) %>% mutate(chrom=as.character(chrom))
        }else{
            gene_annot <- gtf_gene %>% dplyr::select(c(
                "gene_id","gene_name","gene_type","chrom"="seqnames"
            )) %>% mutate(chrom=as.character(chrom))
        }
        message(sprintf("[ %s ] fix ERCC annotations....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
        gene_annot <- gene_annot %>% mutate(
            gene_name = case_when( 
                grepl("ERCC-", gene_id) ~ gene_id,
                TRUE ~ gene_name
            ),
            gene_type = case_when(
                grepl("ERCC-", gene_id) ~ "ERCC-spikein",
                TRUE ~ gene_type
            )       
        )
        gene_annot <- as.data.frame(gene_annot) %>% dplyr::rename(gene=gene_name)
        return(gene_annot)
    }
####

#### Cell-Marker 2.0 Annotation
    # cell-marker filtering
    FilterCellMarkerSet <- function( cellmarker_data, genespace, selected_tissues, min_marker_n ){
        if( !is.null(selected_tissues) ){
            cellmarker_use <- cellmarker_data %>% filter( tissue %in% selected_tissues ) %>% filter( gene %in% genespace )
        }else{
            cellmarker_use <- cellmarker_data %>% filter( gene %in% genespace )
        }
        cellmarker_use_stat <- cellmarker_use %>% group_by(index) %>% reframe(n_markers=length(unique(gene))) %>% as.data.frame() 
        cellmarker_index_filtered <- cellmarker_use_stat %>% filter( n_markers >= min_marker_n ) %>% .$index
        cellmarker_set <- cellmarker_use %>% filter( index %in% cellmarker_index_filtered ) %>% 
            #distinct(index, gene) %>%
            group_by(index, tissue, tissue_type, cell_type, cancer_type, cell_name) %>%
            reframe( cellmarker_genes = list(unique(gene)), n_cellmarker = n_distinct(gene)  ) %>% 
            as.data.frame()
        return(cellmarker_set)
    }
    # annotation 
    AnnotateCellMarkerSet <- function( sce, cellmarker_set, is_seurat=TRUE ) {
        # exprmat
        if( is_seurat )
        {
            message(sprintf("[ %s ] cell-marker annotation with droplet scRNAseq data....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
            sce2 <- sce
            sce2[["RNA"]] <- JoinLayers(sce2[["RNA"]])
            expr_gene_mat <- GetAssayData(sce2, assay="RNA", layer="data")
            message(sprintf("[ %s ] calculate scors....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
            module_score_mat <- sapply(cellmarker_set$cellmarker_genes, function(marker_genes) {
                genes_use <- intersect(marker_genes, rownames(expr_gene_mat))
                if (length(genes_use) < 3) { return(rep(NA, ncol(expr_gene_mat))) }
                colMeans(expr_gene_mat[genes_use, , drop = FALSE])
            })
            colnames(module_score_mat) <- cellmarker_set$index
            module_score_df            <- as.data.frame(module_score_mat)
            module_score_df$cluster    <- sce$seurat_clusters
            cluster_module_score <- module_score_df %>% group_by(cluster) %>% summarise(
                across( where(is.numeric), ~ { 
                    v <- .x
                    if (all(is.na(v))) { NA } else { mean(v, na.rm = TRUE) }
                })
            )
            cluster_module_score_rev   <- reshape2::melt( cluster_module_score, 
                id.vars = "cluster", variable.name = "CellMarker_index", value.name = "module_score"
            )
            annot_result <- cluster_module_score_rev %>% group_by(cluster) %>% arrange(dplyr::desc(module_score)) %>% 
                reframe( Tissue_CellType=CellMarker_index[1], Score=round(module_score[1],3) ) %>% as.data.frame()

            sce$CellMarkerAnnot   <- as.character(annot_result[match(as.character(sce$seurat_clusters), annot_result$cluster), "Tissue_CellType"])
            sce$CellMarker_tissue <- sapply(sce$CellMarkerAnnot, function(y) unlist(strsplit(y, ":"))[1])
            sce$CellMarker_cell   <- sapply(sce$CellMarkerAnnot, function(y) gsub("_", " ", unlist(strsplit(y, ":"))[2]))
            sce$CellMarkerScore   <- as.numeric(annot_result[match(as.character(sce$seurat_clusters), annot_result$cluster), "Score"])
            sce$CellMarkerAnnot   <- gsub("_", " ", sce$CellMarkerAnnot)
            sce$CellMarkerAnnot   <- gsub(":", " : ", sce$CellMarkerAnnot)
        }else{
            message(sprintf("[ %s ] cell-marker annotation with SmartSeq2 data....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
            cm_expr     <- logcounts(sce)   
            gene_symbol <- rowData(sce)$gene
            names(gene_symbol) <- rownames(sce)
            expr_df   <- data.frame(gene = gene_symbol, as.matrix(cm_expr), check.names = FALSE )
            expr_df   <- expr_df[!is.na(expr_df$gene) & expr_df$gene != "", ] %>% data.table()
            message(sprintf("[ %s ] aggregate expression by gene-level....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
            expr_gene <- expr_df[, lapply(.SD, mean, na.rm=TRUE), by=gene] %>% as.data.frame() 
            rownames(expr_gene) <- expr_gene$gene
            expr_gene_mat <- as.matrix(expr_gene[, -1])
            message(sprintf("[ %s ] calculate scors....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
            module_score_mat <- sapply(cellmarker_set$cellmarker_genes, function(marker_genes) {
                genes_use <- intersect(marker_genes, rownames(expr_gene_mat))
                if (length(genes_use) < 3) { return(rep(NA, ncol(expr_gene_mat))) }
                colMeans(expr_gene_mat[genes_use, , drop = FALSE])
            })
            colnames(module_score_mat) <- cellmarker_set$index
            rownames(module_score_mat) <- colnames(sce)
            module_score_df            <- as.data.frame(module_score_mat)
            module_score_df$cluster    <- sce$cluster
            #cluster_module_score       <- aggregate( . ~ cluster, data = module_score_df, FUN = mean )
            cluster_module_score <- module_score_df %>% group_by(cluster) %>% summarise(
                across( where(is.numeric), ~ { 
                    v <- .x
                    if (all(is.na(v))) { NA } else { mean(v, na.rm = TRUE) }
                })
            )
            cluster_module_score_rev   <- reshape2::melt( cluster_module_score, 
                id.vars = "cluster", variable.name = "CellMarker_index", value.name = "module_score"
            )   
            annot_result <- cluster_module_score_rev %>% group_by(cluster) %>% arrange(dplyr::desc(module_score)) %>% 
            reframe( Tissue_CellType=CellMarker_index[1], Score=round(module_score[1],3) ) %>% as.data.frame()

            colData(sce)$CellMarkerAnnot   <- as.character(annot_result[match(sce$cluster, annot_result$cluster), "Tissue_CellType"])
            colData(sce)$CellMarker_tissue <- sapply(colData(sce)$CellMarkerAnnot, function(y) unlist(strsplit(y, ":"))[1])
            colData(sce)$CellMarker_cell   <- sapply(colData(sce)$CellMarkerAnnot, function(y) gsub("_", " ", unlist(strsplit(y, ":"))[2]))
            colData(sce)$CellMarkerScore   <- as.numeric(annot_result[match(sce$cluster, annot_result$cluster), "Score"])
            colData(sce)$CellMarkerAnnot   <- gsub("_", " ", colData(sce)$CellMarkerAnnot)
            colData(sce)$CellMarkerAnnot   <- gsub(":", " : ", colData(sce)$CellMarkerAnnot)
        }
        return(sce)
    }
####

#### SmartSeq2 Functions
    #
    SmartSeq.CountListToMatrix <- function( smartseq_count_list )
    {
        message(sprintf("[ %s ] convert individual count list into merged matrix format....", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))       
        count_list_mats <- Reduce(function(x, y) { merge(x, y, by = c("Geneid", "Length"), all = TRUE) }, smartseq_count_list)
        count_mat <- as.data.frame(count_list_mats[, !c("Geneid", "Length")])
        rownames(count_mat) <- count_list_mats$Geneid
        count_mat <- as.matrix(count_mat)
        storage.mode(count_mat) <- "integer"
        return(count_mat)
    }
    #
    SmartSeq.BasicQC <- function( count_mat, gene_annot, include_ercc=TRUE, draw_plots=TRUE, out_dir=NULL )
    {
        qc_smartseq2 <- data.frame(
            cell_id        = colnames(count_mat),
            total_counts   = colSums(count_mat),
            detected_genes = colSums(count_mat > 0)
        )
        count_mat_genes <- gene_annot[match(rownames(count_mat), gene_annot$gene_id), ]
        # mitochondrial
        mito <- rbind(
            count_mat_genes %>% filter(grepl("chrM|MT",   chrom    )),
            count_mat_genes %>% filter(grepl("^MT-|^mt-", gene_name))
        )
        # ribosomal
        ribo <- count_mat_genes %>% filter(grepl("^Rpl|^Rps|^RPL|^RPS", gene_name))
        # ercc-spikein
        if(include_ercc) { ercc <- count_mat_genes %>% filter(grepl("^ERCC", gene_name)) }

        # qc values
        qc_smartseq2$log10_total_counts   <- Rounding(log10(qc_smartseq2$total_counts   + 1), 3)
        qc_smartseq2$log10_detected_genes <- Rounding(log10(qc_smartseq2$detected_genes + 1), 3)
        qc_smartseq2$pct.mt   <- Rounding(colSums(count_mat[mito$gene_id, , drop = FALSE]) / colSums(count_mat) * 100, 3)
        qc_smartseq2$pct.ribo <- Rounding(colSums(count_mat[ribo$gene_id, , drop = FALSE]) / colSums(count_mat) * 100, 3)
        if(include_ercc){
            qc_smartseq2$pct.ercc <- Rounding(colSums(count_mat[ercc$gene_id, , drop = FALSE]) / colSums(count_mat) * 100, 3)
        }else{
            qc_smartseq2$pct.ercc <- NA
        }

        # plots
        if(draw_plots)
        {
            if(is.null(out_dir)){ out_dir = "./" }
            # density_libsize <- density(qc_smartseq2$log10_total_counts,   na.rm=TRUE)
            # density_gene    <- density(qc_smartseq2$log10_detected_genes, na.rm=TRUE)
            # density_mt      <- density(qc_smartseq2$pct.mt,               na.rm=TRUE)
            # density_ribo    <- density(qc_smartseq2$pct.ribo,             na.rm=TRUE)
            png(sprintf("%s/smartseq2_qc.initial.status.summary.png", out_dir), width=21, height=14, units="in", res=150 )
            par(mfrow = c(2,3))
            # plot( density_libsize, main="Library size",   xlab="log10(total counts + 1)",   yalb="Density" )
            # plot( density_gene,    main="Detected genes", xlab="log10(detected genes + 1)", yalb="Density" )
            hist( qc_smartseq2$log10_total_counts,   breaks=50, main="Library size",   xlab="log10(total counts + 1)"   )
            hist( qc_smartseq2$log10_detected_genes, breaks=50, main="Detected genes", xlab="log10(detected genes + 1)" )
            plot( x=qc_smartseq2$log10_total_counts, y=qc_smartseq2$log10_detected_genes,
                pch=19, cex=1, col = ifelse(qc_smartseq2$pct.mt > 20, "#E41A1C", "#7a7a7a"),
                xlab = "log10(total counts + 1)", ylab = "log10(detected genes + 1)", 
                main = "SmartSeq2 Basic QC"
            )
            legend( "bottomright", legend=c("pct.mt <= 20", "pct.mt > 20"), col=c("#7a7a7a", "#E41A1C"), pch=19, bty="n" )
            # plot( density_mt,   breaks=50, main="Mitochondrial %", xlab="pct.mt"   ,ylab = "Density")
            # plot( density_ribo, breaks=50, main="Ribosomal %",     xlab="pct.ribo" ,ylab = "Density")
            hist( qc_smartseq2$pct.mt,   breaks=50, main="Mitochondrial %", xlab="pct.mt"   )
            hist( qc_smartseq2$pct.ribo, breaks=50, main="Ribosomal %",     xlab="pct.ribo" )
            
            if(include_ercc){ hist( qc_smartseq2$pct.ercc, breaks=50, main="ERCC-spikein %",  xlab="pct.ercc" ) }
            dev.off()
        }
        write.table( qc_smartseq2, sprintf("%s/smartseq2_qc.initial.status.summary.tsv", out_dir), quote=F, row.names=F, col.names=T, sep="\t")
        return(qc_smartseq2)

    }
    #
####


