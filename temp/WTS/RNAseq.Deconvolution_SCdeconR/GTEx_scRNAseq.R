
library(plyr)
library(dplyr)
library(Matrix)
library(Seurat)
library(setwidth)
library(data.table)
library(openxlsx)
library(ribiosUtils)

work_dir = "/storage/home/kangsm/analysis/gcx.Organoid_Science/RNAseq_deconvolution/GTEx_snRNAseq"
setwd(work_dir)

exp_design = fread("ExpDesign-E-ANND-2.tsv")
exp_design_rev = exp_design[,c(1,4,12,14,16,26,28)]
colnames(exp_design_rev) = c("assay","gtex_id","tissue","tissue_site","protocol","cell_type","author_cell_type")
exp_design_rev = exp_design_rev %>% arrange( cell_type, author_cell_type )
cell_type_summary1 = unique(exp_design_rev[,c("tissue","tissue_site","cell_type","author_cell_type")])
cell_type_summary2 = cell_type_summary1 %>% arrange( tissue, cell_type, author_cell_type)
write.xlsx(list(cell_type_summary1,cell_type_summary2) , "GTEx_Cell.Type.Summary.v1.0.xlsx")

ExpDesign <- read.xlsx("GTEx_Cell.Type.Summary.v1.1.xlsx", sheet=2) ###

#---------------------------------------------------------------------------------------------------

mtx_file      = "E-ANND-2.aggregated_filtered_normalised_counts.mtx"
cells_file    = "E-ANND-2.aggregated_filtered_normalised_counts.mtx_cols"
features_file = "E-ANND-2.aggregated_filtered_normalised_counts.mtx_rows"
GTEX_MTX      = ReadMtx(mtx=mtx_file, cells=cells_file, features=features_file)
GeneMap = fread("E-ANND-2.aggregated_filtered_normalised_counts.mtx_rows", header=FALSE) %>% as.data.frame()
colnames(GeneMap) = c("ens_geneid", "gene")


gtex_cell_type_markers <- fread("E-ANND-2.marker_genes_inferred_cell_type_-_ontology_labels.tsv")
TopCellTypeMarkers     <- gtex_cell_type_markers %>% filter( rank %in% seq(0,249,1) )

gcx_cell_type <- unique(ExpDesign$gcx_deconv_type4)
gcx_cell_type <- gcx_cell_type[!is.na(gcx_cell_type)]


    EXPR_DF  = data.frame()
    FEATURES = data.frame()
    BARCODES = c()
    META     = data.frame()
    
    for( i in 1:length(gcx_cell_type) )
    {
        selected_cells <- unique(exp_design_rev[which(exp_design_rev$cell_type %in% ExpDesign[which(ExpDesign$gcx_deconv_type4 == gcx_cell_type[i]), "cell_type"]), "assay"]$assay)

        selected_genes <- unique(TopCellTypeMarkers[which(TopCellTypeMarkers$cluster %in% ExpDesign[which(ExpDesign$gcx_deconv_type4 == gcx_cell_type[i]), "cell_type"]), ]$genes)
        selected_genes_rev <- GeneMap[match(selected_genes, GeneMap$ens_geneid),"gene"]

        selected_expr <- GTEX_MTX[selected_genes_rev, selected_cells]

        expr_df <- matrix2longdf( mat=selected_expr, longdf.colnames = c("gene", "sampleid", "counts") )

        EXPR_DF  <- rbind(EXPR_DF, expr_df)
        BARCODES <- c(BARCODES, selected_cells)
        FEATURES <- rbind(FEATURES, data.frame(g1=selected_genes_rev, g2=selected_genes_rev, g3="Gene Expression"))

        META <- rbind(META, data.frame(cellid=selected_cells, celltype=gcx_cell_type[i], subjectid="gcx1", cohort="gcx_gtex1"))  
        print(gcx_cell_type[i])
    }


library(Matrix)
library(dplyr)
library(tidyr)

mat <- df %>% pivot_wider(id_cols = gene,values_from = count,names_from = cell,
values_fill = list(count=0)) %>% tibble::column_to_rownames("gene")

mat = Matrix(as.matrix(mat),sparse=TRUE)


