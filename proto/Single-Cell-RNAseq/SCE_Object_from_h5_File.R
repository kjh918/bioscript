options(stringsAsFactors=FALSE) 
suppressPackageStartupMessages({
    library(rhdf5)
    library(Seurat)
    library(Matrix)
    library(plyr)
    library(dplyr)
})
SysTime = function(){format(Sys.time(), "%Y-%m-%d %H:%M:%S")}

DataSetList <- list.files("/storage/home/kangsm/myDB/single-cell-RNAseq-objects/TISCH2", pattern="Seurat.Object.RDS")
DataSetList <- gsub("_Seurat.Object.RDS", "", DataSetList)


TISCH_DATASET_NAME <- "PRAD_GSE172301"

TISCH2_DIR   <- "/storage/home/kangsm/myDB/single-cell-RNAseq-objects/TISCH2"
SCdeconR_DIR <- "/storage/references_and_index/RNAseq/deconvolution_reference"

for( TISCH_DATASET_NAME in DataSetList ) {

    message(sprintf("[ %s ] %s processing start.", SysTime(), TISCH_DATASET_NAME))
   
    h5_file         <- sprintf("%s/%s_expression.h5", TISCH2_DIR, TISCH_DATASET_NAME)
    meta_file       <- sprintf("%s/%s_CellMetainfo_table.tsv", TISCH2_DIR, TISCH_DATASET_NAME)
    sce_output_file <- sprintf("%s/%s_Seurat.Object.RDS", TISCH2_DIR, TISCH_DATASET_NAME)
    
    ## data structure 확인
    # info <- h5ls(h5_file, recursive = TRUE)
    # print(info)

    expr        <- Read10X_h5(filename=h5_file, use.names=TRUE, unique.features=TRUE )
    zero_counts <- sparseMatrix( i=integer(0), j=integer(0), dims=dim(expr), dimnames=dimnames(expr) )
    sce         <- suppressWarnings(CreateSeuratObject(counts=zero_counts, project=TISCH_DATASET_NAME, min.cells=0, min.features=0))
    LayerData(object=sce, assay="RNA", layer="data") <- expr

    meta   <- read.delim(meta_file, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
    common <- intersect(colnames(sce), rownames(meta))
    sce    <- subset(sce, cells=common)
    meta   <- meta[colnames(sce), , drop = FALSE]
    sce    <- AddMetaData( sce, metadata = meta )
    saveRDS(sce, file=sce_output_file)

    phenodata <- meta %>% tibble::rownames_to_column("cellid") %>% dplyr::select(c(
        "cellid","celltypes"="Celltype (major-lineage)","subjectid"="Patient"
    ))
    rownames(phenodata) <- phenodata$cellid
    SCdeconR_ref_output <- sprintf("%s/%s.SCdeconR.reference.RDS", SCdeconR_DIR, TISCH_DATASET_NAME)
    scdecon_ref <- list(
        seurat_data = sce,
        pheno_data  = phenodata,
        ussage = "assay='RNA', layer='data'"
    )
    saveRDS(scdecon_ref, file=SCdeconR_ref_output)

    message(sprintf("[ %s ] %s processing done. %s", SysTime(), TISCH_DATASET_NAME, paste0(rep("=", 50), collapse="")))
}


## 
    ref_data = readRDS("/data/wts/RUO_Service/WTS-26-03/analysis/deconvolution/HTAPP_deconv.ref_brca.liver.RDS")
    SCdeconR_ref_output <- sprintf("%s/%s.SCdeconR.reference.RDS", SCdeconR_DIR, "BRCA_metastasis_human_liver")

    scdecon_ref <- list(
        seurat_data = ref_data$seurat_data,
        pheno_data  = ref_data$pheno_data,
        ussage = "assay='SCT', layer='data'"
    )
    saveRDS(scdecon_ref, file=SCdeconR_ref_output)
