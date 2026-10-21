
# Cell-Type-Specific DE analysis -------------------------------------------------------------------
    library(ENIGMA)
    load("exampleDatasets.Rdata")

    egm <- create_ENIGMA(bulk=Nulk, ref=Ref, ref_type="aggre")
    egm <- get_cell_proportion(egm, method="RLR")
    egm <- ENIGMA_trace_norm(egm, solver="proximalpoint", pos=FALSE, preprocess="none", alphar=0.1)

    p <- 100
    y <- as.numeric(gl(2, p/2)) - 1
    DEG <- FindCSE_DEG(egm, y)

    DEG$CellType1 ..
#---------------------------------------------------------------------------------------------------

# github deconvolution
    res = GeneSigTest(object,filtering=TRUE)
    head(res$call)
    head(res$pval)
    egm = res$egm # the filtered ENIGMA object
#---------------------------------------------------------------------------------------------------


# Estimate CSE on ASD patients dataset -------------------------------------------------------------
    # sample data : https://doi.org/10.5281/zenodo.7462559

    # load datasets and build ENIGMA class
    library(ENIGMA)
    load("brain_data.Rdata")
    load("/storage/home/kangsm/analysis/gcx.Organoid_Science/RNAseq_deconvolution/ENIGMA_setup/ENIGMA_analysis_data/brain_analysis/brain_data.Rdata")
    egm = create_ENIGMA(bulk = Bulk, ref = profile,ref_type = "aggre")

    # calculate cell type fraction and run ENIGMA trace model
    # Here I run two models with different preprocessing methods: sqrt or log
    egm = get_cell_proportion(egm, method = "RLR")
    egm = ENIGMA_trace_norm(egm, alpha = 0.1,model_tracker = TRUE,model_name = "sqrt",preprocess = "sqrt")
    egm = ENIGMA_trace_norm(egm, alpha = 0.1,model_tracker = TRUE,model_name = "log", preprocess = "log")

    # compare with ground truth
    cse = sce2array(egm,norm_output = FALSE,model_name = "log")
    ## ENIGMA would construct CSE into raw gene expression space so need to perform log transformation for visualization
    umap_tr = umap::umap(apply(log2(cse+1), 1, as.vector))
    ggpubr::ggscatter(data.frame(umap = umap_tr$layout, cell = rep(colnames(profile), each = ncol(Bulk))), x = "umap.1", y = "umap.2", color = "cell", size = 1) + ggtitle('Ground truth')

    # Then compare with ground truth, compare on the sample-level
    celltype <- dimnames(CSE_array)[[1]]
    gene <- intersect(dimnames(cse)[[1]],dimnames(CSE_array)[[2]])

    cor_sample_egm <- cor_sample_bulk <- cor_sample_ref <- NULL
    celltype_list <- NULL
    for(k in celltype){
    for(i in 1:dim(CSE_array)[3]){
        cor_sample_egm <- c(cor_sample_egm, cor(cse[gene,i,k],CSE_array[k,gene,i],method = "sp"))
        cor_sample_bulk <- c(cor_sample_bulk, cor(Bulk[gene,i],CSE_array[k,gene,i],method = "sp"))
        cor_sample_ref <- c(cor_sample_ref, cor(profile[gene,k],CSE_array[k,gene,i],method = "sp"))
        celltype_list <- c(celltype_list,k)
    }
    }

    dat <- data.frame(cor = c(cor_sample_egm,cor_sample_bulk,cor_sample_ref), celltype = rep(celltype_list,3), method = c(rep("ENIGMA",length(cor_sample_egm)),rep("Bulk",length(cor_sample_egm)),rep("Reference",length(cor_sample_egm))))
    dat$method <- as.factor(dat$method)
    dat$method <- factor(dat$method,level = c("Reference","Bulk","ENIGMA"))
    ggplot2::ggplot(dat,aes(x = celltype, y = cor, fill = method)) + ggplot2::geom_boxplot()
#---------------------------------------------------------------------------------------------------

# Identify cell type-specific cell states ----------------------------------------------------------


