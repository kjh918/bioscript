
### TUTMORIALS
# https://liuy12.github.io/SCdeconR/articles/Cell_type_deconvolution_and_differential_analysis.html

library(SCdeconR)
library(plotly)
# register backend for parallel processing
library(doFuture)
registerDoFuture()
plan("multisession", workers = 4)
options(doFuture.rng.onMisuse="ignore")

# ref1 + ref2
    ref_list_all <- c(
        paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample1"), 
        paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample2")
    )
    phenodata_list_all <- c(
        paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample1.txt"), 
        paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample2.txt")
    )

# ref1
    ref_list_1 <- c( paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample1"), paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample1") ) 
    phenodata_list_1 <- c( paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample1.txt"), paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample1.txt") )
        
# ref2
    ref_list_2 <- c( paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample2"), paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample2")  )
    phenodata_list_2 <- c( paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample2.txt"), paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample2.txt") )


ref_data_all <- construct_ref(
    ref_list        = ref_list_all, 
    phenodata_list  = phenodata_list_all, 
    data_type       = "cellranger",
    method          = "harmony", 
    group_var       = "subjectid", 
    vars_to_regress = "percent_mt",
    verbose         = TRUE
)
saveRDS(ref_data_all, file="/storage/home/kangsm/myDB/rds/SCdeconR_refdata.RDS")

ref_data_1 <- construct_ref(
    ref_list        = ref_list_1, 
    phenodata_list  = phenodata_list_1, 
    data_type       = "cellranger",
    method          = "harmony", 
    group_var       = "subjectid", 
    vars_to_regress = "percent_mt",
    verbose         = TRUE
)

ref_data_2 <- construct_ref(
    ref_list        = ref_list_2, 
    phenodata_list  = phenodata_list_2, 
    data_type       = "cellranger",
    method          = "harmony", 
    group_var       = "subjectid", 
    vars_to_regress = "percent_mt",
    verbose         = TRUE
)



# construct integrated reference using harmony algorithm
refdata  <- construct_ref(ref_list = ref_list, phenodata_list = phenodata_list, data_type = "cellranger",
    method = "harmony", group_var = "subjectid", nfeature_rna = 50, vars_to_regress = "percent_mt",
    verbose = FALSE)
saveRDS(refdata,  file="/storage/home/kangsm/myDB/rds/SCdeconR_refdata.RDS" )

# create phenodata
phenodata <- data.frame(cellid = colnames(refdata), celltypes = refdata$celltype, subjectid = refdata$subjectid)
saveRDS(phenodata, file="/storage/home/kangsm/myDB/rds/SCdeconR_ref_phenodata.RDS")

# Here is the number of cells per cell-type in refdata
table(refdata$celltype)

### generate two sets of bulk data that have different cell-type proportions.
# equal proportions across cell-types
prop1 <- data.frame(celltypes = unique(refdata$celltype), proportion = rep(0.125, 8))
# generate 20 artificial bulk samples with the above cell-type proportions
bulk_sim1 <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"), phenodata = phenodata,
    num_mixtures = 20, prop = prop1, replace = TRUE)

# high proportion of epithelial cells
prop2 <- data.frame(celltypes = unique(refdata$celltype), proportion = c(0.8, 0.1, 0.1, rep(0, 5)))
# generate 20 artificial bulk samples with the above cell-type proportions
bulk_sim2 <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"), phenodata = phenodata,
    num_mixtures = 2, prop = prop2, replace = TRUE)

# combine the two datasets
bulk_sim <- list(cbind(bulk_sim1[[1]], bulk_sim2[[1]]), cbind(bulk_sim1[[2]], bulk_sim2[[2]]))

# deconvolution
decon_res <- scdecon(bulk = bulk_sim2[[1]], ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
    phenodata = phenodata, filter_ref = TRUE, decon_method = "OLS", norm_method_sc = "LogNormalize",
    norm_method_bulk = "TMM", trans_method_sc = "none", trans_method_bulk = "log2", marker_strategy = "all")
#
prop_barplot(prop = decon_res[[1]], interactive = FALSE)


# Feature Gene SPace
    FeatureGenes <- Features(refdata)
    # INTS11 -> CPSF3L = INTS11

    FeatureGenes_rev <- c(FeatureGenes, "CPSF3L")
    # gtf <- as.data.frame(rtracklayer::import('/storage/references_and_index/hg19/gtf/hg19.gtf'))
    # saveRDS(gtf, file="/storage/home/kangsm/myDB/rds/hg19_gtf.RDS")
    # hg38_gtf=as.data.frame(rtracklayer::import('/storage/references_and_index/hg38/gtf/hg38.gtf'))
    # saveRDS(hg38_gtf, file="/storage/home/kangsm/myDB/rds/hg38_gtf.RDS")
    gtf <- readRDS("/storage/home/kangsm/myDB/rds/hg19_gtf.RDS")

    FeatureGenesGTF  <- data.frame(
        gene = FeatureGenes,
        hg19_gtf[match(FeatureGenes, hg19_gtf$gene_name),c("gene_id","gene_type","gene_status","gene_name")]
    )
    
    hg38_gtf=readRDS("/storage/home/kangsm/myDB/rds/hg38_gtf.RDS")
    hg19_gtf=readRDS("/storage/home/kangsm/myDB/rds/hg19_gtf.RDS")
    hg19_gtf_gene <- unique(hg19_gtf[,c("gene_id","gene_type","gene_status","gene_name")])
    hg19_gtf_gene$ens_geneid = sapply(hg19_gtf_gene$gene_id, function(y) unlist(strsplit(y, "\\."))[1])

    for( missRow in which(is.na(FeatureGenesGTF$gene_id)) )
    {
        miss_gene <- FeatureGenesGTF$gene[ missRow ]
        hg38_geneid <- hg38_gtf %>% filter( gene_name == miss_gene ) %>% .$gene_id %>% unique()
        if( length(hg38_geneid) == 1 )
        {
            hg19_info <- hg19_gtf_gene %>% filter( ens_geneid == unlist(strsplit( hg38_geneid, "\\." ))[1] )
            if( nrow(hg19_info) == 1) 
            {
                FeatureGenesGTF[missRow, "gene_id"    ] <- hg19_info$gene_id
                FeatureGenesGTF[missRow, "gene_type"  ] <- hg19_info$gene_type
                FeatureGenesGTF[missRow, "gene_status"] <- hg19_info$gene_status
                FeatureGenesGTF[missRow, "gene_name"  ] <- hg19_info$gene_name
                #print( FeatureGenesGTF$gene[ missRow ] )
            }else{
                message(sprintf("hg38 geneid exist, but no in hg19 : %s", FeatureGenesGTF$gene[ missRow ]))
                next
            }
        }else{
            message(sprintf("not in hg38 : %s", FeatureGenesGTF$gene[ missRow ]))
            next
        }
    }

    FeatureMap <- openxlsx::read.xlsx()
    FeatureMap <- FeatureMap[which(!is.na(FeatureMap$gene_id)), ]
    saveRDS(FeatureMap, file="/storage/home/kangsm/myDB/rds/SCdeconR_FeatureGenes_hg19.RDS")
    
#--------------------------------------------------------------------------------------------------#

### ANALYSIS
    suppressPackageStartupMessages(library("SCdeconR"))
    #suppressPackageStartupMessages(library("doFuture"))
    suppressPackageStartupMessages(library("Cairo"))
    suppressPackageStartupMessages(library("ggplot2"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
    #options(doFuture.rng.onMisuse="ignore")

    #registerDoFuture()
    #plan("multisession", workers = 4)

    # Load Reference Data
    refdata <- readRDS("/storage/home/kangsm/myDB/rds/SCdeconR_refdata.RDS")
    
    # FeauterGeneSpace Matrix
    hg19_GTF        <- readRDS("/storage/home/kangsm/myDB/rds/hg19_gtf.RDS")
    FeatureGenesGTF <- readRDS("/storage/home/kangsm/myDB/rds/SCdeconR_FeatureGenes_hg19.RDS")

    cell_list <- c(
        "HDFn_mRNA","HDFn_total",
        "HCC827_rep1","HCC827_rep2","HCC827_rep3",
        "MDAMB231_rep1","MDAMB231_rep2","MDAMB231_rep3",
        "MDAMB453_rep1","MDAMB453_rep2","MDAMB453_rep3",
        "NCIH596_rep1","NCIH596_rep2","NCIH596_rep3",
        "SKBR3_rep1","SKBR3_rep2","SKBR3_rep3",
        "T47D_rep1","T47D_rep2","T47D_rep3"
    )

    cnt_list <- lapply( cell_list, function(y) readRDS(sprintf("/data/wts/RDS/RDS_ReadCount/%s.ReadCount.rds", y)) )

    CNT = data.frame(cnt_list[[1]])
    for( i in 2:length(cnt_list) ){ CNT = data.frame( CNT, cnt_list[[i]])}

    CNT_rev <- CNT[FeatureGenesGTF$gene_id, ] %>% as.data.frame() %>% mutate( gene = FeatureGenesGTF$gene )
    InputRNaseq <- CNT_rev %>% group_by( gene ) %>% summarise_all(sum) %>% as.data.frame()
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
    suppressPackageStartupMessages(library("Cairo"))
    grDevices::X11.options(type='cairo')
    options(device='x11')
    options(bitmapType='cairo')
    prop_barplot(prop = deconvolution_result[[1]], interactive = FALSE)

    prop = deconvolution_result[[1]]

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
    
    
    
    
    prop$sampleid = gsub("HCC827_rep","HCC-827_R", prop$sampleid)
    prop$sampleid = gsub("MDAMB231_rep","MDA-MB-231_R", prop$sampleid)
    prop$sampleid = gsub("MDAMB453_rep","MDA_MB-453_R", prop$sampleid)
    prop$sampleid = gsub("NCIH596_rep","NCI-H596_R", prop$sampleid)
    prop$sampleid = gsub("SKBR3_rep","SK-BR-3_R", prop$sampleid)
    prop$sampleid = gsub("T47D_rep","T-47D_R", prop$sampleid)
    prop$sampleid = gsub("HDFn_mRNA","HDFn_R1", prop$sampleid)
    prop$sampleid = gsub("HDFn_total","HDFn_R2", prop$sampleid)

    prop$sampleid = factor( prop$sampleid, levels=c(
        "HDFn_R1","HDFn_R2","MDA-MB-231_R1", "MDA-MB-231_R2", "MDA-MB-231_R3","NCI-H596_R1","NCI-H596_R2",   "NCI-H596_R3",
        "HCC-827_R1",    "HCC-827_R2" ,   "HCC-827_R3","SK-BR-3_R1",    "SK-BR-3_R2",    "SK-BR-3_R3",
        "MDA_MB-453_R1", "MDA_MB-453_R2" ,"MDA_MB-453_R3","T-47D_R1"  ,    "T-47D_R2"   ,   "T-47D_R3"
    ))

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
    ggsave(gp, file="SCdeconR_setup_CCLE_test.png", width=16, height=8, unit='in', dpi=150 )


    # /storage/home/kangsm/analysis/gcx.Organoid_Science/RNAseq_deconvolution/test_organoid_rnaseq/GSE277147
    valid_data_dir = "/storage/home/kangsm/analysis/gcx.Organoid_Science/RNAseq_deconvolution/test_organoid_rnaseq/GSE277147"
    valid_res_files = list.files(valid_data_dir, pattern="^GSM")


    i=1
    VALID_DATA = fread(sprintf("%s/%s", valid_data_dir, valid_res_files[k]), header=F)
    colnames(VALID_DATA) = c("gene", unlist(strsplit(valid_res_files[k], "\\."))[1] )

    for( k in 2:length(valid_res_files) )
    {
        dat = fread(sprintf("%s/%s", valid_data_dir, valid_res_files[k]), header=F)
        colnames(dat) = c("gene", unlist(strsplit(valid_res_files[k], "\\."))[1] )
        VALID_DATA = merge(VALID_DATA, dat, by="gene")
    }
    VALID_DATA = as.data.frame(VALID_DATA)
    rownames(VALID_DATA) = VALID_DATA$gene
    InputRNaseq2 <- VALID_DATA[ which(VALID_DATA$gene %in% FeatureGenesGTF$gene) , ]

    deconvolution_result2 <- scdecon(
        bulk              = InputRNaseq2[,-1], 
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

    deconv_barplot = SCdeconR.Predicted.Proportion.Barplot(prop=deconvolution_result2[[1]], sampleOrder=NULL)

    ggsave(deconv_barplot, file="SCdeconR_setup_GSE277147_test.png", width=10, height=8, unit='in', dpi=150 )



    
    cnt1 <- readRDS("/data/wts/RDS/RDS_ReadCount/HDFn_mRNA.ReadCount.rds")
    cnt2 <- readRDS("/data/wts/RDS/RDS_ReadCount/HDFn_total.ReadCount.rds")
    cnt3 <- readRDS("/data/wts/RDS/RDS_ReadCount/HCC827_rep1.ReadCount.rds")
    
    CNT$gene <- hg19_GTF[match(rownames(CNT), hg19_GTF$gene_id), "gene_name"]
    CNT_FeatureSpace <- CNT %>% filter( gene %in% FeatureGenesGTF$gene_name )
    CNT_FeatureSpace_Rev <- CNT_FeatureSpace %>% group_by( gene ) %>% reframe(
        HDFn_mRNA=sum(HDFn_mRNA), HDFn_total=sum(HDFn_total), HCC827=sum(HCC827_rep1)
    ) %>% as.data.frame()
    rownames(CNT_FeatureSpace_Rev) <- CNT_FeatureSpace_Rev$gene

    CNT_Rev <- CNT %>% group_by( gene ) %>% reframe(
        HDFn_mRNA=sum(HDFn_mRNA), HDFn_total=sum(HDFn_total), HCC827=sum(HCC827_rep1)
    ) %>% as.data.frame()
    rownames(CNT_Rev) <- CNT_Rev$gene


    rownames(CNT_FeatureSpace_Rev) <- CNT_FeatureSpace_Rev$gene


    CNT <- data.frame(s1=cnt1, s2=cnt2, s3=cnt3)
    colnames(CNT) <- c("h1","h2")
    CNT$gene <- rownames(CNT)
    CNT$gene <- hg19_GTF[match(CNT$gene, hg19_GTF$gene_id), "gene_name"]
    CNT_FeatureSpace <- CNT %>% filter( gene %in% FeatureGenesGTF$gene_name )

    rownames(CNT) <- CNT$gene
    rownames(CNT_FeatureSpace) <- CNT_FeatureSpace$gene

    CNT_FeatureSpace_Rev <- CNT_FeatureSpace %>% group_by(gene) %>% reframe( h1=sum(h1), h2=sum(h2) ) %>% as.data.frame()
    rownames(CNT_FeatureSpace_Rev) <- CNT_FeatureSpace_Rev$gene

    


    # Load Bulk-RNAseq Data
    bulk_sim1 <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"), phenodata = phenodata,
        num_mixtures = 20, prop = prop1, replace = TRUE)

    # 
    # Perform Deconvolution
    deconvolution_result <- scdecon(
        bulk = bulk_sim2[[1]], 
        ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
        phenodata = phenodata, 
        filter_ref = TRUE, 
        decon_method = "OLS", 
        norm_method_sc = "LogNormalize",
        norm_method_bulk = "TMM", 
        trans_method_sc = "none", 
        trans_method_bulk = "log2", 
        marker_strategy = "all"
    )

    deconvolution_result <- scdecon(
        bulk = CNT_Rev[,-1], 
        ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
        phenodata = phenodata, 
        filter_ref = TRUE, 
        decon_method = "OLS", 
        norm_method_sc = "LogNormalize",
        norm_method_bulk = "TMM", 
        trans_method_sc = "none", 
        trans_method_bulk = "log2", 
        marker_strategy = "all",
        min_pct_ct = 0.03,
        lfc_markers = log2(1.3)
    )
    prop_barplot(prop = deconvolution_result[[1]], interactive = FALSE)






MTX_FILE='/storage/apps/R/4.3.1/lib/R/library/SCdeconR/extdata/refdata/sample2/matrix.mtx.gz'
CELL_FILE='/storage/apps/R/4.3.1/lib/R/library/SCdeconR/extdata/refdata/sample2/barcodes.tsv.gz'
FEATURE_FILE='/storage/apps/R/4.3.1/lib/R/library/SCdeconR/extdata/refdata/sample2/features.tsv.gz'

mm = Seurat::ReadMtx(mtx=MTX_FILE, cells=CELL_FILE, features=FEATURE_FILE)

gene cell count
a    c1    1
a    c2    1
a    c3    4
b    c1    3
b    c2    1
b    c3    1
f    c1    3
d    c8    9
e    c11   1

library(Matrix)
library(dplyr)
library(tidyr)

mat <- df %>% pivot_wider(id_cols = gene,values_from = count,names_from = cell,
values_fill = list(count=0)) %>% tibble::column_to_rownames("gene")

mat = Matrix(as.matrix(mat),sparse=TRUE)


