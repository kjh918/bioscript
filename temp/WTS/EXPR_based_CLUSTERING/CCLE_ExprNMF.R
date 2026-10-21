
#---| PACKAGES |-----------------------------------------------------------------------------------3
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("umap"))
    suppressPackageStartupMessages(library("openxlsx"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("scales"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("NMF"))
    suppressPackageStartupMessages(library("ConsensusClusterPlus"))
    suppressPackageStartupMessages(library("parallel"))
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    ccle_data_dir = "/storage/home/kangsm/myDB/ccle_rnaseq"
    res_dir       = "/storage/home/kangsm/analysis/gcx.rnaseq.based.clustering"

#---| LOAD CELL LINE INFO |------------------------------------------------------------------------#

    #cinfo <- read.xlsx(sprintf("%s/Cell_Lines_Info_Master_v1.0_20220321.xlsx", ccle_data_dir))
    load(sprintf("%s/CCLE_RNAseq_v2_TPM_Matrix.Rdata", ccle_data_dir))          # TpmMat , REF 
    cinfo = read.xlsx(sprintf("%s/CCLE_Cell_Info_rev.xlsx", ccle_data_dir))

    cinfo_stat = cinfo %>% group_by(tissue, disease) %>% reframe( cells = length(cellid) ) %>% 
        arrange( tissue, desc(cells) )
    CINFO = cinfo %>% filter( tissue %nin% c("Blood","Lymph"))
    SOLID_TUMOR_TPM = TpmMat[,CINFO$cellid]
    
    SOLID_TUMOR_TPM_FILTER  = SOLID_TUMOR_TPM[which(apply(SOLID_TUMOR_TPM, 1, function(y) mad(y) ) != 0), ]
    SOLID_TUMOR_TPM_FILTER2 = SOLID_TUMOR_TPM_FILTER[apply(SOLID_TUMOR_TPM_FILTER, 1, function(y) quantile(y, 0.95) > 1 ), ]

    SOLID_TUMOR_TPM_Z = t(apply(SOLID_TUMOR_TPM_FILTER2, 1, function(z) (z-median(z))/mad(z)))

#---| GENE FILTERING |-----------------------------------------------------------------------------#

    indexing = function(x){ sapply(unique(x), function(z) list(z)) }
    gidx = indexing(rownames(SOLID_TUMOR_TPM_FILTER2))    
    DIDX = indexing(CINFO$disease)
    TIDX = indexing(CINFO$tissue)

    MarkerList_disease = lapply( DIDX, function(ts) 
    {
        message( sprintf("diease = %s", ts))
        cid1 = CINFO %>% filter( disease == ts ) %>% .$cellid
        cid2 = CINFO %>% filter( disease != ts ) %>% .$cellid
        
        ts_markers = mclapply( gidx,  function(y) 
        {
            gv1 = as.matrix(SOLID_TUMOR_TPM_FILTER2)[y, cid1]
            gv2 = as.matrix(SOLID_TUMOR_TPM_FILTER2)[y, cid2]
            res = data.frame( 
                T_test = t.test(gv1,gv2)$p.value,
                W_test = wilcox.test(gv1,gv2)$p.value,
                Mean1 = mean(gv1),
                Mean2 = mean(gv2),
                Median1 = median(gv1),
                Median2 = median(gv2),
                Var1    = var(gv1)
            )
            return(res)
        }, mc.cores= 60)
        ts_markers = ldply(ts_markers, .id="gene") %>% mutate( gene = as.character(gene) ) 
        ts_markers = ts_markers %>%
            mutate( mean_gap = Mean2 - Mean1, mean_ratio = Mean1/Mean2, median_gap = Median2 - Median1, median_ratio = Median1/Median2 )
        
        ts_markers = ts_markers %>% filter( T_test < 0.05 , W_test < 0.05 )
        return(ts_markers)
    })    

    MarkerList_tissue = lapply( TIDX, function(ts) 
    {
        message( sprintf("diease = %s", ts))
        cid1 = CINFO %>% filter( tissue == ts ) %>% .$cellid
        cid2 = CINFO %>% filter( tissue != ts ) %>% .$cellid
        
        ts_markers = mclapply( gidx,  function(y) 
        {
            gv1 = as.matrix(SOLID_TUMOR_TPM_FILTER2)[y, cid1]
            gv2 = as.matrix(SOLID_TUMOR_TPM_FILTER2)[y, cid2]
            res = data.frame( 
                T_test = t.test(gv1,gv2)$p.value,
                W_test = wilcox.test(gv1,gv2)$p.value,
                Mean1 = mean(gv1),
                Mean2 = mean(gv2),
                Median1 = median(gv1),
                Median2 = median(gv2),
                Var1    = var(gv1)
            )
            return(res)
        }, mc.cores= 60)
        ts_markers = ldply(ts_markers, .id="gene") %>% mutate( gene = as.character(gene) ) 
        ts_markers = ts_markers %>%
            mutate( mean_gap = Mean2 - Mean1, mean_ratio = Mean1/Mean2, median_gap = Median2 - Median1, median_ratio = Median1/Median2 )
        
        ts_markers = ts_markers %>% filter( T_test < 0.05 , W_test < 0.05 )
        return(ts_markers)
    }) 

    MarkerList_disease2 = lapply( DIDX, function(tsn) 
    {
        message( sprintf("diease = %s", tsn))
        cid1 = CINFO %>% filter( disease == tsn ) %>% .$cellid
        cid2 = CINFO %>% filter( disease != tsn ) %>% .$cellid
        
        ts_markers = mclapply( gidx,  function(y) 
        {
            gv1 = as.matrix(SOLID_TUMOR_TPM_FILTER2)[y, cid1]
            gv2 = as.matrix(SOLID_TUMOR_TPM_FILTER2)[y, cid2]
            res = data.frame( 
                T_test = t.test(gv1,gv2)$p.value,
                W_test = wilcox.test(gv1,gv2)$p.value,
                Mean1 = mean(gv1),
                Mean2 = mean(gv2),
                Median1 = median(gv1),
                Median2 = median(gv2)
            )
            return(res)
        }, mc.cores= 60)
        ts_markers = ldply(ts_markers, .id="gene") %>% mutate( gene = as.character(gene) ) 
        ts_markers = ts_markers %>%
            mutate( mean_gap = Mean2 - Mean1, mean_ratio = Mean1/Mean2, median_gap = Median2 - Median1, median_ratio = Median1/Median2 )
        
        ts_markers = ts_markers %>% filter( Mean2 == 0 | Median2 == 0 )
        return(ts_markers)
    })    

    MD1 = lapply( MarkerList_disease, function(mkd) 
    {
        mkd %>% filter( T_test < 0.005 , W_test < 0.005, mean_ratio > 1.5, median_ratio > 1.5 ) %>% dplyr::select(!c("Var1"))
    })

    MD2 = lapply( MarkerList_tissue, function(mkd) 
    {
        mkd %>% filter( T_test < 0.005 , W_test < 0.005, mean_ratio > 1.5, median_ratio > 1.5 ) %>% dplyr::select(!c("Var1"))
    })

    MD1_stat = mclapply( gidx, function(gn)
    {
        chk = lapply( MD1, function(y) gn %in% y$gene )
        res = data.frame(matrix(unlist(chk), nrow=1))
        colnames(res) = names(MD1)
        return(res)
    }, mc.cores=40)
    MD1_stat = ldply(MD1_stat, .id="gene") %>% mutate(gene = as.character(gene))

    MD2_stat = mclapply( gidx, function(gn)
    {
        chk = lapply( MD2, function(y) gn %in% y$gene )
        res = data.frame(matrix(unlist(chk), nrow=1))
        colnames(res) = names(MD2)
        return(res)
    }, mc.cores=40)
    MD2_stat = ldply(MD2_stat, .id="gene") %>% mutate(gene = as.character(gene))

    MD1_filter = lapply( DIDX, function(ds) 
    {
        ds_name    = ds
        no_ds_name = CINFO[which(CINFO$disease != ds_name), "disease"] %>% unique()
        ts_name    = CINFO[which(CINFO$disease == ds_name), "tissue"] %>% unique()
        no_ts_name = CINFO[which(CINFO$tissue != ts_name), "tissue"] %>% unique()
        ds_group_name    = CINFO[which(CINFO$tissue == ts_name), "disease"] %>% unique()
        no_ds_group_name = CINFO[which(CINFO$tissue != ts_name), "disease"] %>% unique()

        mk = MD1[[ds_name]]
        if( nrow(mk) > 0 )
        {
            mk_status = data.frame()
            for( k in 1:nrow(mk) )
            {
                mks = data.frame(
                    gene = mk$gene[k],
                    Matched.Disease = sum(MD1_stat[which(MD1_stat$gene == mk$gene[k]), ds_name]),
                    NonMatched.Disease = sum(MD1_stat[which(MD1_stat$gene == mk$gene[k]), no_ds_name]),
                    TissueMatched.DiseaseGroup = sum(MD1_stat[which(MD1_stat$gene == mk$gene[k]), ds_group_name]),
                    TissueNonMatched.DiseaseGroup = sum(MD1_stat[which(MD1_stat$gene == mk$gene[k]), no_ds_group_name]),
                    Matched.Tissue = sum(MD2_stat[which(MD2_stat$gene == mk$gene[k]), ts_name]),
                    NonMatched.Tissue = sum(MD2_stat[which(MD2_stat$gene == mk$gene[k]), no_ts_name])
                )
                mk_status = rbind(mk_status, mks)
            }
        }else{
            mk_status = data.frame()
        }
        message(ds)
        return(mk_status)
    })

    MD2_filter = lapply( TIDX, function(y) 
    {
        ts_name    = y
        no_ts_name = CINFO[which(CINFO$tissue != ts_name), "tissue"] %>% unique()
        ds_group_name    = CINFO[which(CINFO$tissue == ts_name), "disease"] %>% unique()
        no_ds_group_name = CINFO[which(CINFO$tissue != ts_name), "disease"] %>% unique()

        mk = MD2[[ts_name]]
        if( nrow(mk) > 0 )
        {
            mk_status = data.frame()
            for( k in 1:nrow(mk) )
            {
                mks = data.frame(
                    gene = mk$gene[k],
                    Matched.Tissue = sum(MD2_stat[which(MD2_stat$gene == mk$gene[k]), ts_name]),
                    NonMatched.Tissue = sum(MD2_stat[which(MD2_stat$gene == mk$gene[k]), no_ts_name]),
                    TissueMatched.DiseaseGroup = sum(MD1_stat[which(MD1_stat$gene == mk$gene[k]), ds_group_name]),
                    TissueNonMatched.DiseaseGroup = sum(MD1_stat[which(MD1_stat$gene == mk$gene[k]), no_ds_group_name])
                )
                mk_status = rbind(mk_status, mks)
            }
        }else{
            mk_status = data.frmae()
        }
        message(y)
        return(mk_status)
    })

    save( MD1, MD2, MD1_stat, MD2_stat, MD1_filter, MD2_filter, file="TissueSpecificMarkers_Candidates.Rdata" )  
    
    names(MD1_filter) = gsub("High grade ovarian serous adenocarcinoma", "Ovarian serous AD", names(MD1_filter))
    names(MD1_filter) = gsub("Invasive ductal carcinoma, not otherwise specified", "Invasive ductal carcinoma", names(MD1_filter))  
    names(MD1_filter) = gsub("adenocarcinoma", "AD", names(MD1_filter))
    names(MD1_filter) = gsub("carcinoma", "CA", names(MD1_filter))

    MD1_filter_rev = lapply(MD1_filter, function(y)
    {
        # y %>% filter(
        #     NonMatched.Disease == 0,
        #     Matched.Tissue == 1,
        #     NonMatched.Tissue  == 0
        # )
        y[which(
            y$NonMatched.Disease == 0 &
            y$Matched.Tissue == 1 &
            y$NonMatched.Tissue  == 0 ) , ]
        
    })
    MD2_filter_rev = lapply(MD2_filter, function(y)
    {
        y %>% filter(
            NonMatched.Tissue == 0,
            TissueMatched.DiseaseGroup == 1,
            TissueNonMatched.DiseaseGroup == 0
        )
    })
    openxlsx::write.xlsx(MD1_filter_rev, "DiseaseLevel_Markers_Candidates.xlsx")
    
    MD2_filter_rev2 = lapply(TIDX, function(y)
    {
        mk = MD2_filter_rev[[y]]
        dsg = CINFO[which(CINFO$tissue == y), "disease"] %>% unique()

        dsg = gsub("High grade ovarian serous adenocarcinoma", "Ovarian serous AD", dsg)
        dsg = gsub("Invasive ductal carcinoma, not otherwise specified", "Invasive ductal carcinoma", dsg)  
        dsg = gsub("adenocarcinoma", "AD", dsg)
        dsg = gsub("carcinoma", "CA", dsg)
        mk_ds = MD1_filter_rev[dsg]

        mk$DS = sapply(mk$gene, function(z) 
        {
            paste(names(which(lapply( mk_ds, function(w) length(which(w$gene == z)) ) == 1)), collapse=";")
        })
        return(mk)      
    })
    openxlsx::write.xlsx(MD2_filter_rev2, "TissueLevel_Markers_Candidates.xlsx")

#--------------------------------------------------------------------------------------------------#

    markerSummary = data.frame(rbind(
        c("Bladder", "Bladder carcinoma", "BLCA" ),
        c("Soft_Tissue", "Ewing sarcoma", "EW-SARC"),
        c("Soft_Tissue", "Extrarenal rhabdoid tumor", "RH-SARC"),
        c("Brain", "Astrocytoma", "LGG"),
        c("Brain", "Glioblastoma", "GBM"),
        c("Breast", "Breast adenocarcinoma", "BRCA"),
        c("Breast", "Breast carcinoma", "BRCA"),
        c("Breast", "Ductal breast carcinoma", "BRCA",
        c("Breast", "Invasive ductal carcinoma, not otherwise specified", "BRCA"),
        c("Colorectal", "Cecum adenocarcinoma", "COAD"),
        c("Colorectal", "Rectal adenocarcinoma", "COAD"),
        c("Colorectal", "Colon adenocarcinoma", "COAD"),
        c("Esophagus",  "Esophageal squamous cell carcinoma", "ESCA"),
        c("Kidney",     "Clear cell renal cell carcinoma", "KIRC"),
        c("Kidney", "Renal cell carcinoma", "KIRC"),
        c("Liver", "Hepatocellular carcinoma", "LIHC"),
        c("Lung", "Large cell lung carcinoma", "LULC"),
        c("Lung", "Lung adenocarcinoma", "LUAD"),
        c("Lung", "Non-small cell lung carcinoma", "NSCLC"),
        c("Lung", "Squamous cell lung carcinoma", "LUSC"),
        c("Lung", "Small cell lung carcinoma", "SCLC"),
        c("Ovaray", "High grade ovarian serous adenocarcinoma", "OV"),
        c("Pancreas", "Pancreatic ductal adenocarcinoma", "PAAD"),
        c("Pancreas", "Pancreatic adenocarcinoma", "PAAD"),
        c("Peripheral Nerve", "Neuroblastoma", "NBL"),
        c("Prostate", "Prostate carcinoma", "PRAD"),
        c("Skin", "Melanoma", "SKCM"),
        c("Skin", "Cutaneous melanoma", "SKCM"),
        c("Stomach", "Gastric adenocarcinoma", "STAD"),
        c("Uterus", "Endometrial adenocarcinoma", "UCEC")
    ))
    colnames(markerSummary) = c("Tissue","Disease","Code")


    SCLC_MK = MD1[["Small cell lung carcinoma"]] %>% filter(
        gene %in% MD2_filter_rev2[["LungS"]]$gene
    ) %>% filter( mean_ratio > 5 ) %>% .$gene
    MD2_filter_rev2[["LungS"]] = MD2_filter_rev2[["LungS"]] %>% filter( gene %in% SCLC_MK )

    NBL_MK = MD1[["Neuroblastoma"]] %>% filter(
        gene %in% MD2_filter_rev2[["Peripheral Nerve"]]$gene
    ) %>% filter( mean_ratio > 2 ) %>% .$gene
    MD2_filter_rev2[["Peripheral Nerve"]] = MD2_filter_rev2[["Peripheral Nerve"]] %>% filter( gene %in% NBL_MK )

    COAD_MK = MD1[["Colon adenocarcinoma"]] %>% filter(
        gene %in% MD2_filter_rev2[["Colorectal"]]$gene
    ) %>% filter( mean_ratio > 2 ) %>% .$gene
    MD2_filter_rev2[["Colorectal"]] = MD2_filter_rev2[["Colorectal"]] %>% filter( gene %in% COAD_MK )

    MD2_filter_rev3 = lapply( MD2_filter_rev2, function(y)
    {
        DS2 = markerSummary$Disease
        DS2 = gsub("High grade ovarian serous adenocarcinoma", "Ovarian serous AD", DS2)
        DS2 = gsub("Invasive ductal carcinoma, not otherwise specified", "Invasive ductal carcinoma", DS2)  
        DS2 = gsub("adenocarcinoma", "AD", DS2)
        DS2 = gsub("carcinoma", "CA", DS2)

        y[which(y$DS %in% DS2), ]
    })

    ccle_tissue_markers = lapply(MD2_filter_rev3, function(y) y$gene )
    ccle_tissue_markers = ccle_tissue_markers[names(ccle_tissue_markers) != c("Prostate")]
    ccle_tissue_markers_genes = unlist(ccle_tissue_markers)

    save( MD1_filter_rev, MD2_filter_rev, MD2_filter_rev2, markerSummary, MD2_filter_rev3, ccle_tissue_markers, ccle_tissue_markers_genes,
     file="CCLE_CLUSTERING_MARKERS_v1.Rdata")


#---| UMAP

    TISSUE_COLORS1 = c(
        "Bladder carcinoma"                        = "#00b2b2", 
        "Glioblastoma"                             = "#E78AC3", 
        "Astrocytoma"                              = "#FCCDE5", 
        "Ductal breast carcinoma"                  = "#FFBF00",
        "Invasive ductal carcinoma, not otherwise specified"="#FFBF00",
        "Colon adenocarcinoma"                     = "#563e00", 
        "Esophageal squamous cell carcinoma"       = "#8DD3C7", 
        "Renal cell carcinoma"                     = "#FB8072", 
        "Hepatocellular carcinoma"                 = "#E41A1C", 
        "Lung adenocarcinoma"                      = "#0040FF", 
        "Small cell lung carcinoma"                = "#006EAF", 
        "High grade ovarian serous adenocarcinoma" = "#000059", 
        "Pancreatic adenocarcinoma"                = "#66C2A5", 
        "Pancreatic ductal adenocarcinoma"         = "#66C2A5", 
        "Neuroblastoma"                            = "#8DA0CB", 
        "Melanoma"                                 = "#4DAF4A", 
        "Cutaneous melanoma"                       = "#4DAF4A", 
        "Ewing sarcoma"                            = "#323299",
        "Extrarenal rhabdoid tumor"                = "#5900b2", 
        "Gastric adenocarcinoma"                   = "#FF7F0E", 
        "Endometrial adenocarcinoma"               = "#9932ff"
    )
    TISSUE_COLORS2 = c(
        "GBM"     = "#E78AC3", 
        "LGG"     = "#FCCDE5", 
        "BRCA"    = "#FFBF00",
        "ESCA"    = "#295914", 
        "KIRC"    = "#ac7c00", 
        "LIHC"    = "#E41A1C", 
        "STAD"    = "#FF7F0E", 
        "COAD"    = "#563e00", 
        "LUAD"    = "#0040FF", 
        "SCLC"    = "#1BA1E2", 
        "PAAD"    = "#66C2A5", 
        "SKCM"    = "#4DAF4A", 
        "EW-SARC" = "#5900b2",
        "RH-SARC" = "#9932ff", 
        "NBL"     = "#8DA0CB", 
        "OV"      = "#323299", 
        "UCEC"    = "#A6D854"
    )
    Draw.CCLE.UMAP.Clustering <- function( ExprMatrix, CellInfo, PointSize=8, 
        UMAP_METRIC = c("cosine","pearson","pearson2"),
        NEIGHBORS=NA, TISSUE_COLORS, ZSCORE=TRUE )
    {
        require(ggplot2)
        require(umap)
        require(dplyr)
        require(threejs)
        # tissue filtering
        message(" >> tissue filtering...")
        exprM = ExprMatrix[, CellInfo$cellid]

        # gene variance
        message(" >> gene filtering...")
        # gene filtering and expression values transformation
        if( ZSCORE )
        {
            mad_chk = apply(exprM, 1, mad)
            inputM = exprM[ names(which(mad_chk != 0)), ]
            INPUT = t(apply(inputM, 1, function(z) (z-median(z))/mad(z)))

        }else{
            INPUT = exprM
        }
        # umap config
        umapConfig              = umap.defaults
        umapConfig$n_components = 2
        if( is.na(NEIGHBORS) )
        { umapConfig$n_neighbors = as.integer(ncol(INPUT)^(1/2)) }else{ umapConfig$n_neighbors = NEIGHBORS } 
        umapConfig$metric       = UMAP_METRIC
        umapConfig$random_state = 123
        umapConfig$transform_state = 123
        umapConfig$knn_repeats = 3
        # run umap
        message(" >> Run UMAP 2D...")
        umapResult = umap(t(INPUT), config=umapConfig)
        umapPlotData = as.data.frame(umapResult$layout[colnames(INPUT),])
        umapPlotData = data.frame(
            cellid=rownames(umapPlotData),
            X=umapPlotData[,1],
            Y=umapPlotData[,2],
            CellInfo[match(rownames(umapPlotData), CellInfo$cellid), c("cell_name","tissue","code")]
        )
        # plots
        message(" >> create plot2D...")
        umapPlots2D <- ggplot(umapPlotData, aes(x=X, y=Y, fill=code, text=cell_name)) + geom_point(col="white", pch=21, size=PointSize, alpha=0.9, stroke=0.5) +
            scale_fill_manual(values=TISSUE_COLORS) +
            labs( x="", y="", fill="Cancer Type") +
            theme(
                axis.line.x  = element_blank(), 
                axis.line.y  = element_blank(),
                axis.text.x  = element_blank(),
                axis.text.y  = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour = "#a3a3a3", fill=NA, size=1),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'), 
                plot.title=element_text(face="bold", size=14),
                plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                legend.key.size = unit(3, 'mm')
            )
        #
        message(" >> Run UMAP 3D...")
        umapConfig              = umap.defaults
        umapConfig$n_components = 3
        if( is.na(NEIGHBORS) )
        { umapConfig$n_neighbors = as.integer(ncol(INPUT)^(1/2)) }else{ umapConfig$n_neighbors = NEIGHBORS } 
        umapConfig$metric       = UMAP_METRIC
        umapConfig$random_state = 123
        umapConfig$transform_state = 123
        umapConfig$knn_repeats = 3
        umapResult2 = umap(t(INPUT), config=umapConfig)
        umapPlotData2 = as.data.frame(umapResult2$layout[colnames(INPUT),])
        umapPlotData2 = data.frame(
            cellid=rownames(umapPlotData2),
            X=umapPlotData2[,1],
            Y=umapPlotData2[,2],
            Z=umapPlotData2[,3],
            CellInfo[match(rownames(umapPlotData2), CellInfo$cellid), c("cell_name","tissue","disease")]
        )
        message(" >> create plot3D...")
        umapPlots3D = scatterplot3js( x = umapPlotData2$X, y = umapPlotData2$Y, z = umapPlotData2$Z, color = umapPlotData2$COLORS, labels = umapPlotData2$disease, axis=TRUE, num.ticks = c(0,0,0), stroke = '#e5e5e5' )

        ResultValues = list(
            data = umapResult,
            plot2D = umapPlots2D,
            plot3D = umapPlots3D
        )
        message(" >> DONE.")
        return(ResultValues)
    }

    cinfo_umap = CINFO %>% filter( disease %in% names(TISSUE_COLORS1) ) 
    cinfo_umap$code = markerSummary[match(cinfo_umap$disease, markerSummary$Disease), "Code"]

    library(plotly)
    res_umap_zexpr_cos = Draw.CCLE.UMAP.Clustering(
        ExprMatrix = as.matrix(SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, ]),
        CellInfo   = cinfo_umap,
        PointSize  = 6,
        UMAP_METRIC = "cosine",
        NEIGHBORS  = NA, 
        TISSUE_COLORS = TISSUE_COLORS2, 
        ZSCORE = TRUE
    )
    res_umap_zexpr_ps = Draw.CCLE.UMAP.Clustering(
        ExprMatrix = as.matrix(SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, ]),
        CellInfo   = cinfo_umap,
        PointSize  = 6,
        UMAP_METRIC = "pearson",
        NEIGHBORS  = NA, 
        TISSUE_COLORS = TISSUE_COLORS2, 
        ZSCORE = TRUE
    )
    res_umap_zexpr_ps2 = Draw.CCLE.UMAP.Clustering(
        ExprMatrix = as.matrix(SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, ]),
        CellInfo   = cinfo_umap,
        PointSize  = 6,
        UMAP_METRIC = "pearson2",
        NEIGHBORS  = NA, 
        TISSUE_COLORS = TISSUE_COLORS2, 
        ZSCORE = TRUE
    )

    ggplotly(res_umap_zexpr_cos$plot2D)
    ggplotly(res_umap_zexpr_ps$plot2D)
    ggplotly(res_umap_zexpr_ps2$plot2D)

#---| NMF
    
    nmf_res = nmf( x=SOLID_TUMOR_TPM_FILTER2, rank=18, method="brunet", nrun=30, seed=123456)

    NMF_INPUT2 = SOLID_TUMOR_TPM_FILTER2[, cinfo_umap$cellid]
    nmf_res2 = nmf( x=NMF_INPUT2, rank=18, method="brunet", nrun=30, seed=123456)

    NMF_INPUT3 = SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, cinfo_umap$cellid]
    nmf_res3 = nmf( x=NMF_INPUT3, rank=18, method="brunet", nrun=30, seed=123456)


    nmf2_consensus = consensusmap(nmf_res2)
    nmf3_consensus = consensusmap(nmf_res3)

    str(nmf2_consensus$Rowv) # h 0.757
    str(nmf3_consensus$Rowv) # h 0.628

    cinfoCluster = cinfo_umap

    for( H in c(0.8, 0.7, 0.6 ) )
    {
        cell_cluster = dendextend::cutree(nmf2_consensus$Rowv, h=H)
        cellClusterMember = lapply( 1:length(unique(cell_cluster)), function(y) names(cell_cluster[cell_cluster==y]) )
        names(cellClusterMember) = 1:length(cellClusterMember)
        cellClusterMemberRev = ldply(lapply(cellClusterMember, function(z) cinfo_umap %>% dplyr::filter( cellid %in% z) ), .id="cluster")
        
        cinfoCluster$nmfCluster = cellClusterMemberRev[match(cinfoCluster$cellid, cellClusterMemberRev$cellid), "cluster"]
        cn = sprintf("nmf2_h%s", H)
        colnames(cinfoCluster)[ncol(cinfoCluster)] = cn
    }
    for( H in c(0.7, 0.6, 0.5) )
    {
        cell_cluster = dendextend::cutree(nmf3_consensus$Rowv, h=H)
        cellClusterMember = lapply( 1:length(unique(cell_cluster)), function(y) names(cell_cluster[cell_cluster==y]) )
        names(cellClusterMember) = 1:length(cellClusterMember)
        cellClusterMemberRev = ldply(lapply(cellClusterMember, function(z) cinfo_umap %>% dplyr::filter( cellid %in% z) ), .id="cluster")
        
        cinfoCluster$nmfCluster = cellClusterMemberRev[match(cinfoCluster$cellid, cellClusterMemberRev$cellid), "cluster"]
        cn = sprintf("nmf3_h%s", H)
        colnames(cinfoCluster)[ncol(cinfoCluster)] = cn
    }

    png("NMF.All.Genes.Consensus.png", width=10, height=10, units="in",res=200  )
    nmf2_consensus = consensusmap(nmf_res2)
    dev.off()
    png("NMF.Marker.Genes.Consensus.png", width=10, height=10, units="in",res=200  )
    nmf3_consensus = consensusmap(nmf_res3)
    dev.off()

    nmf.cluster.id.cinfo = cinfo_umap

    for( H in c(0.8, 0.7, 0.6 ) )
    {
        cell_cluster = dendextend::cutree(nmf2_consensus$Rowv, h=H)
        cellClusterMember = lapply( 1:length(unique(cell_cluster)), function(y) names(cell_cluster[cell_cluster==y]) )
        names(cellClusterMember) = 1:length(cellClusterMember)
        cellClusterMemberRev = ldply(lapply(cellClusterMember, function(z) cinfo_umap %>% dplyr::filter( cellid %in% z) ), .id="cluster")
        
        nmf.cluster.id.cinfo$nmfCluster = cellClusterMemberRev[match(nmf.cluster.id.cinfo$cellid, cellClusterMemberRev$cellid), "cluster"]
        cn = sprintf("nmf2_h%s", H)
        colnames(nmf.cluster.id.cinfo)[ncol(nmf.cluster.id.cinfo)] = cn
    }
    for( H in c(0.7, 0.6, 0.5) )
    {
        cell_cluster = dendextend::cutree(nmf3_consensus$Rowv, h=H)
        cellClusterMember = lapply( 1:length(unique(cell_cluster)), function(y) names(cell_cluster[cell_cluster==y]) )
        names(cellClusterMember) = 1:length(cellClusterMember)
        cellClusterMemberRev = ldply(lapply(cellClusterMember, function(z) cinfo_umap %>% dplyr::filter( cellid %in% z) ), .id="cluster")
        
        nmf.cluster.id.cinfo$nmfCluster = cellClusterMemberRev[match(nmf.cluster.id.cinfo$cellid, cellClusterMemberRev$cellid), "cluster"]
        cn = sprintf("nmf3_h%s", H)
        colnames(nmf.cluster.id.cinfo)[ncol(nmf.cluster.id.cinfo)] = cn
    }


#---| ConsensusClusterPlus
    mad_chk2 = apply(SOLID_TUMOR_TPM_FILTER2[, cinfo_umap$cellid], 1, function(y) mad(y))
    rmRows2 = which(mad_chk2 == 0 )
    if( length(rmRows2) > 0 )
    {
        ccp_mat2 = SOLID_TUMOR_TPM_FILTER2[, cinfo_umap$cellid][-rmRows2, ]
    }
    CCP_INPUT2 = t(apply(ccp_mat2, 1, function(z) (z-median(z))/mad(z)))
    #CCP_INPUT2_rev = sweep(ccp_mat2, 1, apply(ccp_mat2,1,median))

    ccp_res2 = ConsensusClusterPlus(CCP_INPUT2, maxK=20 ,reps=100, pItem=0.9, pFeature=1, title="CCP_All_Genes", clusterAlg="hc", distance="pearson",seed=123456,plot="png")
    ccp_res2_rev = ConsensusClusterPlus(as.matrix(CCP_INPUT2_rev), maxK=20 ,reps=100, pItem=0.9, pFeature=1, title="CCLE_ExprAll2", clusterAlg="hc", distance="pearson",seed=123456,plot="png")

    mad_chk3 = apply(SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, cinfo_umap$cellid], 1, function(y) mad(y))
    rmRows3 = which(mad_chk3 == 0)
    if( length(rmRows3) > 0 )
    {
        ccp_mat3 = SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, cinfo_umap$cellid][-rmRows3, ]
    }
    CCP_INPUT3 = t(apply(ccp_mat3, 1, function(z) (z-median(z))/mad(z)))
    #CCP_INPUT3_rev = sweep(ccp_mat3, 1, apply(ccp_mat3,1,median))

    ccp_res3 = ConsensusClusterPlus(CCP_INPUT3, maxK=20 ,reps=100,pItem=0.9, pFeature=1,title="CCLE_ExprMark", clusterAlg="hc", distance="pearson",seed=123456,plot="png")
    ccp_res3_rev = ConsensusClusterPlus(as.matrix(CCP_INPUT3_rev), maxK=20 ,reps=100,pItem=0.9, pFeature=1,title="CCLE_ExprMark2", clusterAlg="hc", distance="pearson",seed=123456,plot="png")

    # 

    for( CLUSTERS_N in 17:19 )
    {
        cluster_membership = ccp_res2[[ CLUSTERS_N ]]$consensusClass
        #cinfoCluster = cinfoCluster %>% mutate( cluster =  cluster_membership[ cellid ] ) 
        cinfoCluster$cluster = cluster_membership[ cinfoCluster$cellid ]
        cn = sprintf("ccp2_k%s", CLUSTERS_N)
        colnames(cinfoCluster)[ncol(cinfoCluster)] = cn
        message(CLUSTERS_N)
    }
    for( CLUSTERS_N in 17:19 )
    {
        cluster_membership = ccp_res3[[ CLUSTERS_N ]]$consensusClass
        #cinfoCluster = cinfoCluster %>% mutate( cluster =  cluster_membership[ cellid ] ) 
        cinfoCluster$cluster = cluster_membership[ cinfoCluster$cellid ]
        cn = sprintf("ccp3_k%s", CLUSTERS_N)
        colnames(cinfoCluster)[ncol(cinfoCluster)] = cn
        message(CLUSTERS_N)
    }
    for( CLUSTERS_N in 17:19 )
    {
        cluster_membership = ccp_res2_rev[[ CLUSTERS_N ]]$consensusClass
        #cinfoCluster = cinfoCluster %>% mutate( cluster =  cluster_membership[ cellid ] ) 
        cinfoCluster$cluster = cluster_membership[ cinfoCluster$cellid ]
        cn = sprintf("ccp2rev_k%s", CLUSTERS_N)
        colnames(cinfoCluster)[ncol(cinfoCluster)] = cn
        message(CLUSTERS_N)
    }
    for( CLUSTERS_N in 17:19 )
    {
        cluster_membership = ccp_res3_rev[[ CLUSTERS_N ]]$consensusClass
        #cinfoCluster = cinfoCluster %>% mutate( cluster =  cluster_membership[ cellid ] ) 
        cinfoCluster$cluster = cluster_membership[ cinfoCluster$cellid ]
        cn = sprintf("ccp3rev_k%s", CLUSTERS_N)
        colnames(cinfoCluster)[ncol(cinfoCluster)] = cn
        message(CLUSTERS_N)
    }

    openxlsx::write.xlsx(cinfoCluster, "CCLE_NMF_and_ConsensusClusterPlus_Result.xlsx")


    ccp_res[[2]]["consensusTree"]
    icl = calcICL(ccp_res,title="CCLE",plot="png")

    #cinfoCluster = cinfo_umap

    for( CLUSTERS_N in 10:18 )
    {
        cluster_membership = ccp_res[[ CLUSTERS_N ]]$consensusClass
        #cinfoCluster = cinfoCluster %>% mutate( cluster =  cluster_membership[ cellid ] ) 
        cinfoCluster$cluster = cluster_membership[ cinfoCluster$cellid ]
        cn = sprintf("ccp_k%s", CLUSTERS_N)
        colnames(cinfoCluster)[ncol(cinfoCluster)] = cn
        message(CLUSTERS_N)
    }
#--------------------------------------------------------------------------------------------------#


#---| FILTERED CELL RE-CLUSTERING CHECK
    clist01 = openxlsx::read.xlsx("CCLE_NMF_and_ConsensusClusterPlus_Result_rev.xlsx", sheet=2) # NMF2
    clist02 = openxlsx::read.xlsx("CCLE_NMF_and_ConsensusClusterPlus_Result_rev.xlsx", sheet=3) # NMF3
    clist03 = openxlsx::read.xlsx("CCLE_NMF_and_ConsensusClusterPlus_Result_rev.xlsx", sheet=4) # CCP3

    clist04 = rbind(clist02, clist03) %>% unique()

    clist04_rev = clist04 %>% filter( cell_name %nin% c("SW48","HPAF-II","SUIT-2","AsPC-1","TCC-PAN2","AGS", "SW1990","PaTu 8988s","CADO-ES1" ) )

    cinfo_umap2 = CINFO %>% filter( cellid %in% clist04_rev$cellid ) 
    cinfo_umap2$code = markerSummary[match(cinfo_umap2$disease, markerSummary$Disease), "Code"]

    res_umap2 = Draw.CCLE.UMAP.Clustering(
        ExprMatrix = as.matrix(SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, cinfo_umap2$cellid]),
        CellInfo   = cinfo_umap2,
        PointSize  = 5,
        UMAP_METRIC = "cosine",
        NEIGHBORS  = NA, 
        TISSUE_COLORS = TISSUE_COLORS2, 
        ZSCORE = TRUE
    )

    ggplotly(res_umap2$plot2D)

    ggsave(res_umap2$plot2D, file="UMAP.CCLE.Clustering.png", width=10, height=10, unit='in', dpi=200, type='cairo')


    res_umap2$plot3D

#--------------------------------------------------------------------------------------------------#

#---| FINAL RESULT OBJECTS |-----------------------------------------------------------------------#

    # NMF 
        NMF.input.cinfo            = cinfo_umap
        NMF.all.genes.input.mat    = NMF_INPUT2
        NMF.marker.genes.input.mat = NMF_INPUT3
        NMF.all.genes.results      = nmf_res2 # nmf2
        NMF.marker.genes.results   = nmf_res3 # nmf3
        nmf2_consensus             = consensusmap(nmf_res2)
        nmf3_consensus             = consensusmap(nmf_res3)
        NMF.dendrogram.height      = c( "NMF2" = 0.757, "NMF3" = 0.628 )
        nmf.cluster.id.cinfo       # cluster id labeled result


        save(NMF.input.cinfo, NMF.all.genes.input.mat, NMF.marker.genes.input.mat, NMF.all.genes.results,  NMF.marker.genes.results, nmf2_consensus, nmf3_consensus, NMF.dendrogram.height, nmf.cluster.id.cinfo, file = "NMF.Clustering.Data.Rdata" )
    #

    # ConsensusClustering
        CCP.input.cinfo            = cinfo_umap
        CCP.all.genes.input.mat    = CCP_INPUT2
        CCP.marker.genes.input.mat = CCP_INPUT3

        CCP.all.genes.results = ConsensusClusterPlus(CCP.all.genes.input.mat, maxK=20 ,reps=100, pItem=0.9, pFeature=1, title="CCP_All_Genes", clusterAlg="hc", distance="pearson",seed=123456, plot="png")
        CCP.marker.genes.results = ConsensusClusterPlus(CCP.marker.genes.input.mat, maxK=20 ,reps=100, pItem=0.9, pFeature=1, title="CCP_Marker_Genes", clusterAlg="hc", distance="pearson",seed=123456, plot="png")

        CCP.all.genes.icl    = calcICL(CCP.all.genes.results,    title="CCP_All_Genes",    plot="png")
        CCP.marker.genes.icl = calcICL(CCP.marker.genes.results, title="CCP_Marker_Genes", plot="png")

        ccp.cluster.id.cinfo = cinfo_umap
        for( CLUSTERS_N in 17:19 )
        {
            cluster_membership = CCP.all.genes.results[[ CLUSTERS_N ]]$consensusClass
            ccp.cluster.id.cinfo$cluster = cluster_membership[ ccp.cluster.id.cinfo$cellid ]
            cn = sprintf("ccp2_k%s", CLUSTERS_N)
            colnames(ccp.cluster.id.cinfo)[ncol(ccp.cluster.id.cinfo)] = cn
            message(CLUSTERS_N)
        }
        for( CLUSTERS_N in 17:19 )
        {
            cluster_membership = CCP.marker.genes.results[[ CLUSTERS_N ]]$consensusClass
            ccp.cluster.id.cinfo$cluster = cluster_membership[ ccp.cluster.id.cinfo$cellid ]
            cn = sprintf("ccp3_k%s", CLUSTERS_N)
            colnames(ccp.cluster.id.cinfo)[ncol(ccp.cluster.id.cinfo)] = cn
            message(CLUSTERS_N)
        }
        CCP.cluster.K = c("CCP2"=18, "CCP3"=18)

        save( CCP.input.cinfo, CCP.all.genes.input.mat, CCP.marker.genes.input.mat, CCP.all.genes.results, CCP.marker.genes.results, CCP.all.genes.icl, CCP.marker.genes.icl, ccp.cluster.id.cinfo, CCP.cluster.K, file = "CCP.ConsensusClusterPlus.Data.Rdata")
    #

    # Cell Line List
    ccle.expr.clustering.cinfo = clist04_rev

    cinfo.summary = ccle.expr.clustering.cinfo %>% group_by(tissue,disease,code) %>% reframe( cells = length(cellid) ) %>% as.data.frame() %>%
        arrange( tissue, desc(cells))

    openxlsx::write.xlsx(
        list(
            Cell_Info = ccle.expr.clustering.cinfo, 
            Cell_Type_Summary = cinfo.summary
        ),
        "CCLE.ExprClustering.Result_Cell.Info.xlsx"
    )

    TISSUE_COLORS = c(
        "GBM"     = "#E78AC3", 
        "LGG"     = "#FCCDE5", 
        "BRCA"    = "#FFBF00",
        "ESCA"    = "#295914", 
        "KIRC"    = "#ac7c00", 
        "LIHC"    = "#b25a4e", 
        "STAD"    = "#FF7F0E", 
        "COAD"    = "#563e00", 
        "LUAD"    = "#0040FF", 
        "SCLC"    = "#1BA1E2", 
        "PAAD"    = "#66C2A5", 
        "SKCM"    = "#4DAF4A", 
        "EW-SARC" = "#5900b2",
        "RH-SARC" = "#9932ff", 
        "NBL"     = "#8DA0CB", 
        "OV"      = "#323299", 
        "UCEC"    = "#A6D854"
    )
   
    # Markers
    ccle_tissue_markers
    ccle.tissue.expr.marker.genes = ccle_tissue_markers_genes

    save(ccle_tissue_markers, ccle.tissue.expr.marker.genes, file="CCLE.ExprClustering_GeneMarkers.Rdata")

    # Expr Matrix
    ccle.clustering.cell.lines.expr.mat   = SOLID_TUMOR_TPM_FILTER2[ccle_tissue_markers_genes, ccle.expr.clustering.cinfo$cellid]
    ccle.clustering.cell.lines.expr.z.mat = t(apply(ccle.clustering.cell.lines.expr.mat, 1, function(z) (z-median(z))/mad(z)))

    cell.line.median = apply(ccle.clustering.cell.lines.expr.mat, 1, function(w) median(w))
    cell.line.mad    = apply(ccle.clustering.cell.lines.expr.mat, 1, function(w) mad(w))

    cell.line.z.median = apply(ccle.clustering.cell.lines.expr.z.mat, 1, function(w) median(w))
    cell.line.z.mad    = apply(ccle.clustering.cell.lines.expr.z.mat, 1, function(w) mad(w))

    ccle.clustering.cell.lines.expr.stat = list(
        tpm_median   = cell.line.median,
        tpm_mad      = cell.line.mad,
        z_tpm_median = cell.line.z.median,
        z_tpm_mad    = cell.line.z.mad
    )

    save( ccle.clustering.cell.lines.expr.mat, ccle.clustering.cell.lines.expr.z.mat, ccle.clustering.cell.lines.expr.stat, file = "CCLE.ExprClustering.Result_CellLineExprMatrix.Rdata")

    # UMAP Clustering Functions
    Draw.CCLE.UMAP.Clustering <- function( ExprMatrix, CellInfo, PointSize=8, 
        UMAP_METRIC = c("cosine","pearson","pearson2"),
        NEIGHBORS=NA, TISSUE_COLORS, ZSCORE=TRUE )
    {
        require(ggplot2)
        require(umap)
        require(dplyr)
        require(threejs)
        # tissue filtering
        message(" >> tissue filtering...")
        exprM = ExprMatrix[, CellInfo$cellid]

        # gene variance
        message(" >> gene filtering...")
        # gene filtering and expression values transformation
        if( ZSCORE )
        {
            mad_chk = apply(exprM, 1, mad)
            inputM = exprM[ names(which(mad_chk != 0)), ]
            INPUT = t(apply(inputM, 1, function(z) (z-median(z))/mad(z)))

        }else{
            INPUT = exprM
        }
        # umap config
        umapConfig              = umap.defaults
        umapConfig$n_components = 2
        if( is.na(NEIGHBORS) )
        { umapConfig$n_neighbors = as.integer(ncol(INPUT)^(1/2)) }else{ umapConfig$n_neighbors = NEIGHBORS } 
        umapConfig$metric       = UMAP_METRIC
        umapConfig$random_state = 123
        umapConfig$transform_state = 123
        umapConfig$knn_repeats = 3
        # run umap
        message(" >> Run UMAP 2D...")
        umapResult = umap(t(INPUT), config=umapConfig)
        umapPlotData = as.data.frame(umapResult$layout[colnames(INPUT),])
        umapPlotData = data.frame(
            cellid=rownames(umapPlotData),
            X=umapPlotData[,1],
            Y=umapPlotData[,2],
            CellInfo[match(rownames(umapPlotData), CellInfo$cellid), c("cell_name","tissue","code")]
        )
        # plots
        message(" >> create plot2D...")
        umapPlots2D <- ggplot(umapPlotData, aes(x=X, y=Y, fill=code, text=cell_name)) + geom_point(col="white", pch=21, size=PointSize, alpha=0.9, stroke=0.5) +
            scale_fill_manual(values=TISSUE_COLORS) +
            labs( x="", y="", fill="Cancer Type") +
            theme(
                axis.line.x  = element_blank(), 
                axis.line.y  = element_blank(),
                axis.text.x  = element_blank(),
                axis.text.y  = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour = "#a3a3a3", fill=NA, size=1),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'), 
                plot.title=element_text(face="bold", size=14),
                plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                legend.key.size = unit(3, 'mm')
            )
        #
        message(" >> Run UMAP 3D...")
        umapConfig              = umap.defaults
        umapConfig$n_components = 3
        if( is.na(NEIGHBORS) )
        { umapConfig$n_neighbors = as.integer(ncol(INPUT)^(1/2)) }else{ umapConfig$n_neighbors = NEIGHBORS } 
        umapConfig$metric       = UMAP_METRIC
        umapConfig$random_state = 123
        umapConfig$transform_state = 123
        umapConfig$knn_repeats = 3
        umapResult2 = umap(t(INPUT), config=umapConfig)
        umapPlotData2 = as.data.frame(umapResult2$layout[colnames(INPUT),])
        umapPlotData2 = data.frame(
            cellid=rownames(umapPlotData2),
            X=umapPlotData2[,1],
            Y=umapPlotData2[,2],
            Z=umapPlotData2[,3],
            CellInfo[match(rownames(umapPlotData2), CellInfo$cellid), c("cell_name","tissue","disease")]
        )
        message(" >> create plot3D...")
        umapPlots3D = scatterplot3js( x = umapPlotData2$X, y = umapPlotData2$Y, z = umapPlotData2$Z, color = umapPlotData2$COLORS, labels = umapPlotData2$disease, axis=TRUE, num.ticks = c(0,0,0), stroke = '#e5e5e5' )

        ResultValues = list(
            data = umapResult,
            plot2D = umapPlots2D,
            plot3D = umapPlots3D
        )
        message(" >> DONE.")
        return(ResultValues)
    }

    
    # ExprHeatmap
    t(scale(t(ccle.clustering.cell.lines.expr.mat)))[1:5,1:5]
    library(circlize)

    col_fun = colorRamp2(c(-2, 0, 2), c("#14426c", "#ffffff", "#7f271b"))

    TissueOrders = c("SCLC","LUAD","SKCM","GBM","LGG","COAD","NBL","KIRC","ESCA","EW-SARC","RH-SARC","BRCA","PAAD","OV","LIHC","STAD","UCEC", "BLCA")
    GeneOrders   = c("LungS","Lung","Skin","Brain","Colorectal","Peripheral Nerve","Kidney","Esophagus","Soft_Tissue","Breast","Pancreas","Ovary","Liver","Stomach","Uterus","Bladder")

    cinfo_heatmap = data.frame()
    for( k in 1:17)
    {
        cinfo_heatmap = rbind(
            cinfo_heatmap,
            ccle.expr.clustering.cinfo %>% filter( code == TissueOrders[k] ) %>% arrange( NMF3 )
        )
    }
    genes_heatmap = unlist(ccle_tissue_markers[GeneOrders])


    annotCells <- HeatmapAnnotation(
        Tissue = cinfo_heatmap$code,
        col    = list( Tissue = TISSUE_COLORS ),
        gap    = unit(1.2, "mm"),
        annotation_name_gp      = gpar(fontsize = 8, bold=T), 
        annotation_name_side    = 'left', 
        annotation_legend_param = list( title_gp=gpar(fontsize=12, bold=T), labels_gp=gpar(fontsize=9) )        
    )

    png("CCLE.ExprClustering_GeneMarkersHeatmap.png", width=12,height=12,units="in",res=200)
    Heatmap(
        as.matrix( t(scale(t(ccle.clustering.cell.lines.expr.mat[genes_heatmap,cinfo_heatmap$cellid])))),
        col = col_fun,
        name = "TPM",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        top_annotation = annotCells,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 3)
    )
    dev.off()

    
    # Object For Analysis
    
    ccle.expr.clustering.cinfo = cinfo_heatmap

    ccle.clustering.cell.lines.expr.mat = ccle.clustering.cell.lines.expr.mat[genes_heatmap,cinfo_heatmap$cellid]

    ccle.clustering.cell.lines.expr.z.mat = ccle.clustering.cell.lines.expr.z.mat[genes_heatmap,cinfo_heatmap$cellid]

    ccle.expr = SOLID_TUMOR_TPM_FILTER
    class.marker.genes = genes_heatmap

    save( ccle.expr, ccle.expr.clustering.cinfo, TISSUE_COLORS, Draw.CCLE.UMAP.Clustering, ccle.clustering.cell.lines.expr.mat, ccle.clustering.cell.lines.expr.z.mat, ccle.clustering.cell.lines.expr.stat, ccle.tissue.expr.marker.genes, file = "CCLE_ExprClustering_AnalysisData.Rdata")

    save( ccle.expr, ccle.expr.clustering.cinfo, TISSUE_COLORS, Draw.CCLE.UMAP.Clustering, ccle.clustering.cell.lines.expr.mat, ccle.clustering.cell.lines.expr.z.mat, ccle.clustering.cell.lines.expr.stat, ccle.tissue.expr.marker.genes, file = "~/myDB/rds/CCLE_ExprClustering_AnalysisData.Rdata")


    # scores (I)
    indexing = function(x){ sapply(unique(x), function(z) list(z)) }

    cidx = indexing(ccle.expr.clustering.cinfo$cellid)

    ccle.clustering.markers.tpm.list.cell = lapply( cidx, function(cx) 
    {
        as.matrix(ccle.clustering.cell.lines.expr.mat)[,cx]
    })
    ccle.clustering.markers.z.tpm.list.cell = lapply( cidx, function(cx) 
    {
        as.matrix(ccle.clustering.cell.lines.expr.z.mat)[,cx]
    })

    tsidx = indexing( TissueOrders[1:17] )

    ccle.clustering.markers.tpm.list.tissue.median = lapply( tsidx, function(tx)
    {
        ts_cell = ccle.expr.clustering.cinfo %>% filter( code == tx ) %>% .$cellid
        v = apply(as.matrix(ccle.clustering.cell.lines.expr.mat)[, ts_cell], 1, median)
        return(v)
    })
    ccle.clustering.markers.tpm.list.tissue.mean = lapply( tsidx, function(tx)
    {
        ts_cell = ccle.expr.clustering.cinfo %>% filter( code == tx ) %>% .$cellid
        v = apply(as.matrix(ccle.clustering.cell.lines.expr.mat)[, ts_cell], 1, mean)
        return(v)
    })
    ccle.clustering.markers.z.tpm.list.tissue.median = lapply( tsidx, function(tx)
    {
        ts_cell = ccle.expr.clustering.cinfo %>% filter( code == tx ) %>% .$cellid
        v = apply(as.matrix(ccle.clustering.cell.lines.expr.z.mat)[, ts_cell], 1, median)
        return(v)
    })
    ccle.clustering.markers.z.tpm.list.tissue.mean = lapply( tsidx, function(tx)
    {
        ts_cell = ccle.expr.clustering.cinfo %>% filter( code == tx ) %>% .$cellid
        v = apply(as.matrix(ccle.clustering.cell.lines.expr.z.mat)[, ts_cell], 1, mean)
        return(v)
    })

    save( class.marker.genes, ccle.expr.clustering.cinfo, TISSUE_COLORS, Draw.CCLE.UMAP.Clustering, ccle.clustering.cell.lines.expr.mat, ccle.clustering.cell.lines.expr.stat, ccle.tissue.expr.marker.genes, file = "CCLE_ExprClustering_AnalysisData.Rdata")

    save( class.marker.genes, ccle.expr.clustering.cinfo, TISSUE_COLORS, Draw.CCLE.UMAP.Clustering, ccle.clustering.cell.lines.expr.mat, ccle.clustering.cell.lines.expr.stat, ccle.tissue.expr.marker.genes, file = "~/myDB/rds/CCLE_ExprClustering_AnalysisData.Rdata")




    TEST_CELL_ID = "CCL_1830"
    
    TEST_EXPR_VALUES = as.matrix(TpmMat)[class.marker.genes ,TEST_CELL_ID]
    

    CLUSTER.TEST <- function( 
        CCLE_EXPR_MATRIX = ccle.clustering.cell.lines.expr.mat,
        CCLE_MEDIAN      = ccle.clustering.cell.lines.expr.stat$tpm_median, 
        CCLE_MAD         = ccle.clustering.cell.lines.expr.stat$tpm_mad, 
        CCLE_INFO        = ccle.expr.clustering.cinfo, 
        MARKER_GENES     = class.marker.genes,
        TEST_EXPR_VALUES, 
        PointSize      = 5, 
        UMAP_METRIC    = c("cosine","pearson","pearson2"),
        CORES          = 25,
        SCORE_CELL_TOP = NA,
        TISSUE_COLORS  = TISSUE_COLORS,
        NEIGHBORS      = NA, 
        ZSCORE         = TRUE
    )
    {
        require(ggplot2)
        require(umap)
        require(dplyr)
        require(plyr)
        require(lsa)
        require(parallel)
        require(threejs)

        indexing = function(x){ sapply(unique(x), function(z) list(z)) }
        # tissue filtering
        message(" >> create matrix...")
        exprM = as.matrix(data.frame(CCLE_EXPR_MATRIX, TEST=TEST_EXPR_VALUES[rownames(CCLE_EXPR_MATRIX)]))

        # cell info & color
        message(" >> update cell-info & colors...")
        CellInfo = rbind(
            CCLE_INFO[,c("cellid","cell_name","tissue","code")],
            c("TEST","Unknown","Unknown","Unknown")
        )
        TISSUE_COLORS_TEST = c(TISSUE_COLORS, "Unknown" = "#E41A1C")

        # transform
        if( ZSCORE )
        {
            message(" >> transforming expr-matrix...")
            mad_chk = apply(exprM, 1, mad)
            inputM = exprM[ names(which(mad_chk != 0)), ]
            INPUT = apply(inputM, 2, function(y) (y-CCLE_MEDIAN[rownames(inputM)])/CCLE_MAD[rownames(inputM)] )  
        }else{
            INPUT = exprM
        }
        # umap config
        umapConfig              = umap.defaults
        umapConfig$n_components = 2
        if( is.na(NEIGHBORS) )
        { umapConfig$n_neighbors = as.integer(ncol(INPUT)^(1/2)) }else{ umapConfig$n_neighbors = NEIGHBORS } 
        umapConfig$metric       = UMAP_METRIC
        umapConfig$random_state = 1234
        umapConfig$transform_state = 1234
        umapConfig$knn_repeats = 3
        # run umap
        message(" >> Run UMAP ...")
        umapResult = umap(t(INPUT), config=umapConfig)
        umapPlotData = as.data.frame(umapResult$layout[colnames(INPUT),])
        umapPlotData = data.frame(
            cellid=rownames(umapPlotData),
            X=umapPlotData[,1],
            Y=umapPlotData[,2],
            CellInfo[match(rownames(umapPlotData), CellInfo$cellid), c("cell_name","tissue","code")]
        )
        # plots
        message(" >> create UMAP Plot ...")
        umapPlots2D <- ggplot(umapPlotData, aes(x=X, y=Y, fill=code, text=cell_name)) + geom_point(col="white", pch=21, size=PointSize, alpha=0.9, stroke=0.4) +
            scale_fill_manual(values=TISSUE_COLORS_TEST) +
            labs( x="", y="", fill="Cancer Type") +
            theme(
                axis.line.x  = element_blank(), 
                axis.line.y  = element_blank(),
                axis.text.x  = element_blank(),
                axis.text.y  = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour = "#a3a3a3", fill=NA, linewidth=1),
                panel.background = element_rect(fill='white'), 
                panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted'), 
                plot.title=element_text(face="bold", size=14),
                plot.margin  = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
                legend.key.size = unit(3, 'mm')
            )
        # 
        # umap3 config
        umapConfig              = umap.defaults
        umapConfig$n_components = 3
        if( is.na(NEIGHBORS) )
        { umapConfig$n_neighbors = as.integer(ncol(INPUT)^(1/2)) }else{ umapConfig$n_neighbors = NEIGHBORS } 
        umapConfig$metric       = UMAP_METRIC
        umapConfig$random_state = 12345
        umapConfig$transform_state = 12345
        umapConfig$knn_repeats = 3
        # run umap
        message(" >> Run UMAP ...")
        umapResult2 = umap(t(INPUT), config=umapConfig)
        umapPlotData2 = as.data.frame(umapResult2$layout[colnames(INPUT),])
        umapPlotData2 = data.frame(
            cellid=rownames(umapPlotData2),
            X=umapPlotData2[,1],
            Y=umapPlotData2[,2],
            Z=umapPlotData2[,3],
            CellInfo[match(rownames(umapPlotData2), CellInfo$cellid), c("cell_name","tissue","code")]
        )
        umapPlotData2$COLORS = TISSUE_COLORS_TEST[umapPlotData2$code]
        umapPlots3D = scatterplot3js( x = umapPlotData2$X, y = umapPlotData2$Y, z = umapPlotData2$Z, color = umapPlotData2$COLORS, labels = umapPlotData2$code, axis=TRUE, num.ticks = c(0,0,0), stroke = '#e5e5e5' )
        # centroid distance
        centroid_table = umapPlotData2 %>% group_by(code) %>% reframe( cent_x = mean(X), cent_y = mean(Y), cent_z = mean(Z) )
        
        for( k in 1:18 )
        {
            centroid_table[k, "cosine_dist"] = 
                lsa::cosine(
                    as.numeric(c(centroid_table[k,2],  centroid_table[k,3],  centroid_table[k,4]) ), 
                    as.numeric(c(centroid_table[18,2], centroid_table[18,3], centroid_table[18,4]))
                )
            centroid_table[k, "euclidean_dist"] = 
                dist(rbind(
                    as.numeric(c(centroid_table[k,2],  centroid_table[k,3],  centroid_table[k,4]) ), 
                    as.numeric(c(centroid_table[18,2], centroid_table[18,3], centroid_table[18,4]))
                ))
        }
        centroid_table = centroid_table %>% filter( code != "Unknown" ) %>% 
            arrange(desc(cosine_dist)) %>% mutate( cosine_dist_rank = 1:nrow(.)) %>%
            arrange(euclidean_dist) %>% mutate( euclidean_dist_rank = 1:nrow(.))

        # scores
        message(" >> Calculate Tissue Assign Score ...")
        cidx = indexing(CCLE_INFO$cellid)
        tidx = indexing(CCLE_INFO$code)

        ccleExprValueListCell   = lapply( cidx, function(cx) { as.matrix(INPUT)[,cx] } )
        ccleExprValueListTissue = lapply( tidx, function(tx) { apply(as.matrix(INPUT)[,CCLE_INFO[which(CCLE_INFO$code == tx), "cellid"]], 1, median) } ) 
        testExprValue     = as.matrix(INPUT)[,"TEST"]

        message(" >> Calculate Tissue Assign Score : cell-based...")
        scoreCell = ldply(mclapply(ccleExprValueListCell, function(cevl)
        {
            data.frame(
                rho = cor( cevl[names(testExprValue)], testExprValue, method = "spearman"),
                tau = cor( cevl[names(testExprValue)], testExprValue, method = "kendall"),
                cos = lsa::cosine( cevl[names(testExprValue)], testExprValue )[1,1]
            )
        }, mc.cores = CORES ), .id="cellid") %>% mutate( cellid = as.character(cellid) ) %>%
            left_join(., CCLE_INFO , by='cellid')
        
        if( is.na(SCORE_CELL_TOP) ){ TopCutOff = nrow(scoreCell) }else{ TopCutOff = SCORE_CELL_TOP }
        scoreCellSummary = scoreCell %>% arrange( desc(cos) ) %>% head(TopCutOff) %>% 
            group_by(code) %>% reframe( RHO=median(rho), TAU=median(tau), COS=median(cos) ) %>% arrange(desc(COS))

        message(" >> Calculate Tissue Assign Score : tissue-based...")
        scoreTissue = ldply(mclapply(ccleExprValueListTissue, function(cevl)
        {
            data.frame(
                rho = cor( cevl[names(testExprValue)], testExprValue, method = "spearman"),
                tau = cor( cevl[names(testExprValue)], testExprValue, method = "kendall"),
                cos = lsa::cosine( cevl[names(testExprValue)], testExprValue )[1,1]
            )
        }, mc.cores = 20 ), .id="code") %>% arrange(desc(rho))

        message(" >> Summarize All Results ...")
        ResultValues = list(
            plot2D = umapPlots2D,
            plot3D = umapPlots3D,
            score_tissue = scoreTissue,
            score_cell = scoreCellSummary,
            score_cell_raw = scoreCell,
            umap_centroid = centroid_table
        )
        message(" >> DONE.")
        return(ResultValues)
    }



    save( class.marker.genes, ccle.expr.clustering.cinfo, TISSUE_COLORS, Draw.CCLE.UMAP.Clustering, ccle.clustering.cell.lines.expr.mat, ccle.clustering.cell.lines.expr.stat, ccle.tissue.expr.marker.genes, CLUSTER.TEST, file = "CCLE_ExprClustering_AnalysisData.Rdata")

    save( class.marker.genes, ccle.expr.clustering.cinfo, TISSUE_COLORS, Draw.CCLE.UMAP.Clustering, ccle.clustering.cell.lines.expr.mat, ccle.clustering.cell.lines.expr.stat, ccle.tissue.expr.marker.genes, CLUSTER.TEST, file = "~/myDB/rds/CCLE_ExprClustering_AnalysisData.Rdata")



### TEST 

ccle_cinfo = openxlsx::read.xlsx("/storage/home/kangsm/myDB/ccle_rnaseq/CCLE_Cell_Info_rev.xlsx")

load("/storage/home/kangsm/myDB/ccle_rnaseq/CCLE_RNAseq_v2_TPM_Matrix.Rdata")
load( "~/myDB/rds/CCLE_ExprClustering_AnalysisData.Rdata" )

#------------------------------------------------------------------------------#

    cellTypeDisease = data.frame(rbind(
        c("Bladder", "Bladder carcinoma", "BLCA" ),
        c("Soft_Tissue", "Ewing sarcoma", "EW-SARC"),
        c("Soft_Tissue", "Extrarenal rhabdoid tumor", "RH-SARC"),
        c("Brain", "Astrocytoma", "LGG"),
        c("Brain", "Glioblastoma", "GBM"),
        c("Breast", "Breast adenocarcinoma", "BRCA"),
        c("Breast", "Breast carcinoma", "BRCA"),
        c("Breast", "Ductal breast carcinoma", "BRCA"),
        c("Breast", "Invasive ductal carcinoma, not otherwise specified", "BRCA"),
        c("Colorectal", "Cecum adenocarcinoma", "COAD"),
        c("Colorectal", "Rectal adenocarcinoma", "COAD"),
        c("Colorectal", "Colon adenocarcinoma", "COAD"),
        c("Esophagus",  "Esophageal squamous cell carcinoma", "ESCA"),
        c("Kidney",     "Clear cell renal cell carcinoma", "KIRC"),
        c("Kidney", "Renal cell carcinoma", "KIRC"),
        c("Liver", "Hepatocellular carcinoma", "LIHC"),
        c("Lung", "Large cell lung carcinoma", "LULC"),
        c("Lung", "Lung adenocarcinoma", "LUAD"),
        c("Lung", "Non-small cell lung carcinoma", "NSCLC"),
        c("Lung", "Squamous cell lung carcinoma", "LUSC"),
        c("Lung", "Small cell lung carcinoma", "SCLC"),
        c("Ovaray", "High grade ovarian serous adenocarcinoma", "OV"),
        c("Pancreas", "Pancreatic ductal adenocarcinoma", "PAAD"),
        c("Pancreas", "Pancreatic adenocarcinoma", "PAAD"),
        c("Peripheral Nerve", "Neuroblastoma", "NBL"),
        c("Prostate", "Prostate carcinoma", "PRAD"),
        c("Skin", "Melanoma", "SKCM"),
        c("Skin", "Cutaneous melanoma", "SKCM"),
        c("Stomach", "Gastric adenocarcinoma", "STAD"),
        c("Uterus", "Endometrial adenocarcinoma", "UCEC")
    ))
    colnames(cellTypeDisease) = c("Tissue","Disease","Code")
    ccle_cinfo$code = NA
    for( k in 1:nrow(cellTypeDisease))
    {
        ccle_cinfo[which(ccle_cinfo$tissue == cellTypeDisease[k,1] & ccle_cinfo$disease == cellTypeDisease[k,2]), "code"] = cellTypeDisease[k,3]
    }

    testSet = ccle_cinfo %>% filter( !is.na(code) ) %>% filter( code %in% c(names(TISSUE_COLORS))) %>% 
        filter( cellid %nin% colnames(ccle.clustering.cell.lines.expr.mat)) # 173

    testCellIdList = testSet$cellid

    TestResultSummary = data.frame()
    #TestResultData    = list()

    for( TEST_CELL_ID in testCellIdList )
    {
        RunNo = which(testCellIdList == TEST_CELL_ID)

        message(sprintf(">>> %s / 173 : Testing...", RunNo ))
        TEST_EXPR_VALUES = as.matrix(TpmMat)[class.marker.genes ,TEST_CELL_ID]
        cat(">>> analysis ")
        test_res = CLUSTER.TEST( 
            CCLE_EXPR_MATRIX = ccle.clustering.cell.lines.expr.mat,
            CCLE_MEDIAN      = ccle.clustering.cell.lines.expr.stat$tpm_median, 
            CCLE_MAD         = ccle.clustering.cell.lines.expr.stat$tpm_mad, 
            CCLE_INFO        = ccle.expr.clustering.cinfo, 
            MARKER_GENES     = class.marker.genes,
            TEST_EXPR_VALUES, 
            PointSize      = 5, 
            UMAP_METRIC    = "cosine",
            CORES          = 25,
            SCORE_CELL_TOP = NA,
            TISSUE_COLORS  = TISSUE_COLORS,
            NEIGHBORS      = NA, 
            ZSCORE         = TRUE
        )
        cat(" --->>> create figure ")
        ggsave(test_res$plot2D, file = sprintf("TEST/CCLE_CLUSTERING_TEST_%s.png", TEST_CELL_ID), width=8, height=8, unit='in', dpi=150, type='cairo')
        cat(" --->>> summary results ")
        rank_res = data.frame(
            sample_id      = TEST_CELL_ID,
            sample_name    = testSet %>% filter( cellid == TEST_CELL_ID ) %>% .$cell_name,
            sample_tissue  = testSet %>% filter( cellid == TEST_CELL_ID ) %>% .$code,
            sample_disease = testSet %>% filter( cellid == TEST_CELL_ID ) %>% .$disease,
            res_TS_COSINE = test_res$score_tissue %>% arrange(desc(cos)) %>% .$code %>% head(1),
            res_CL_COSINE = test_res$score_cell %>% arrange(desc(COS)) %>% .$code %>% head(1),
            res_UMAP_CENT_EUCLIDEAN = test_res$umap_centroid %>% arrange(euclidean_dist_rank) %>% .$code %>% head(1),
            res_UMAP_CENT_COSINE    = test_res$umap_centroid %>% arrange(cosine_dist_rank) %>% .$code %>% head(1)
        )
        TestResultSummary = rbind(TestResultSummary, rank_res)      
        cat(sprintf(" --->>> %s / 173 : DONE.\n", RunNo ))
    }


    colnames(TestResultSummary)[6] = 'cell_base_result'


    openxlsx::write.xlsx(TestResultSummary, "TEST/CCLE_NonClustedCells_TEST_RESULT.xlsx")

    TestResultSummary$sample_tissue = factor(TestResultSummary$sample_tissue, levels=names(TISSUE_COLORS))
    TestResultSummary$res_TS_COSINE = factor(TestResultSummary$res_TS_COSINE, levels=names(TISSUE_COLORS))
    TestResultSummary$res_CL_COSINE = factor(TestResultSummary$res_CL_COSINE, levels=names(TISSUE_COLORS))
    TestResultSummary$res_UMAP_CENT_EUCLIDEAN = factor(TestResultSummary$res_UMAP_CENT_EUCLIDEAN, levels=names(TISSUE_COLORS))
    TestResultSummary$res_UMAP_CENT_COSINE = factor(TestResultSummary$res_UMAP_CENT_COSINE, levels=names(TISSUE_COLORS))


    list1 = TestResultSummary[,c("sample_tissue","sample_disease")] %>% unique()
    list2 = TestResultSummary$sample_tissue %>% unique()
    list3 = TestResultSummary$sample_disease %>% unique()

    EVAL_RES1 = EVAL_RES2 = EVAL_RES3 = data.frame()
    for( k in 1:nrow(list1))
    {
        eval_res = data.frame(
            sample_tissue = list1[k,1],
            sample_disease = list1[k,2],
            sample_N = TestResultSummary %>% filter( sample_tissue == list1[k,1], sample_disease == list1[k,2]) %>% nrow(.),
            TS_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_tissue == list1[k,1], sample_disease == list1[k,2]), 
                truth=sample_tissue, estimate=res_TS_COSINE
            ) %>% .$.estimate %>% head(1),
            CL_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_tissue == list1[k,1], sample_disease == list1[k,2]), 
                truth=sample_tissue, estimate=res_CL_COSINE
            )%>% .$.estimate %>% head(1),
            UMAP_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_tissue == list1[k,1], sample_disease == list1[k,2]), 
                truth=sample_tissue, estimate=res_UMAP_CENT_EUCLIDEAN
            )%>% .$.estimate %>% head(1)
        )   
        EVAL_RES1 = rbind(EVAL_RES1, eval_res)
    }

    for( k in 1:length(list2))
    {
        eval_res = data.frame(
            sample_tissue = list2[k],
            TS_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_tissue == list2[k]), 
                truth=sample_tissue, estimate=res_TS_COSINE
            ) %>% .$.estimate %>% head(1),
            CL_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_tissue == list2[k]), 
                truth=sample_tissue, estimate=res_CL_COSINE
            )%>% .$.estimate %>% head(1),
            UMAP_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_tissue == list2[k]), 
                truth=sample_tissue, estimate=res_UMAP_CENT_EUCLIDEAN
            )%>% .$.estimate %>% head(1)
        )   
        EVAL_RES2 = rbind(EVAL_RES2, eval_res)
    }

    for( k in 1:length(list3))
    {
        eval_res = data.frame(
            sample_disease = list3[k],
            TS_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_disease == list3[k]), 
                truth=sample_tissue, estimate=res_TS_COSINE
            ) %>% .$.estimate %>% head(1),
            CL_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_disease == list3[k]), 
                truth=sample_tissue, estimate=res_CL_COSINE
            )%>% .$.estimate %>% head(1),
            UMAP_AUC = yardstick::accuracy(data=TestResultSummary %>% filter( sample_disease == list3[k]), 
                truth=sample_tissue, estimate=res_UMAP_CENT_EUCLIDEAN
            )%>% .$.estimate %>% head(1)
        )   
        EVAL_RES3 = rbind(EVAL_RES3, eval_res)
    }

    openxlsx::write.xlsx(EVAL_RES1, "TEST/CCLE_NonClustedCells_TEST_RESULT_extra.xlsx")


    yardstick::accuracy(data=TestResultSummary, 
        truth=sample_tissue, estimate=res_UMAP_CENT_COSINE
    )

#------------------------------------------------------------------------------#

    testSet2 = ccle_cinfo %>% filter( code %in% c("LUSC","LULC","NSCLC")) 
    testCellIdList2 = testSet2$cellid

    LungTestRes = data.frame()
    #TestResultData    = list()

    for( TEST_CELL_ID in testCellIdList2 )
    {
        RunNo = which(testCellIdList2 == TEST_CELL_ID)

        message(sprintf(">>> %s / 51 : Testing...", RunNo ))
        TEST_EXPR_VALUES = as.matrix(TpmMat)[class.marker.genes ,TEST_CELL_ID]
        cat(">>> analysis ")
        test_res = CLUSTER.TEST( 
            CCLE_EXPR_MATRIX = ccle.clustering.cell.lines.expr.mat,
            CCLE_MEDIAN      = ccle.clustering.cell.lines.expr.stat$tpm_median, 
            CCLE_MAD         = ccle.clustering.cell.lines.expr.stat$tpm_mad, 
            CCLE_INFO        = ccle.expr.clustering.cinfo, 
            MARKER_GENES     = class.marker.genes,
            TEST_EXPR_VALUES, 
            PointSize      = 5, 
            UMAP_METRIC    = "cosine",
            CORES          = 25,
            SCORE_CELL_TOP = NA,
            TISSUE_COLORS  = TISSUE_COLORS,
            NEIGHBORS      = NA, 
            ZSCORE         = TRUE
        )
        cat(" --->>> create figure ")
        ggsave(test_res$plot2D, file = sprintf("TEST/LUNG_CANCER_TEST_%s.png", TEST_CELL_ID), width=8, height=8, unit='in', dpi=150, type='cairo')
        cat(" --->>> summary results ")
        rank_res = data.frame(
            sample_id      = TEST_CELL_ID,
            sample_name    = testSet2 %>% filter( cellid == TEST_CELL_ID ) %>% .$cell_name,
            sample_tissue  = testSet2 %>% filter( cellid == TEST_CELL_ID ) %>% .$code,
            sample_disease = testSet2 %>% filter( cellid == TEST_CELL_ID ) %>% .$disease,
            res_TS_COSINE = test_res$score_tissue %>% arrange(desc(cos)) %>% .$code %>% head(1),
            res_CL_COSINE = test_res$score_cell %>% arrange(desc(COS)) %>% .$code %>% head(1),
            res_UMAP_CENT_EUCLIDEAN = test_res$umap_centroid %>% arrange(euclidean_dist_rank) %>% .$code %>% head(1),
            res_UMAP_CENT_COSINE    = test_res$umap_centroid %>% arrange(cosine_dist_rank) %>% .$code %>% head(1)
        )
        LungTestRes = rbind(LungTestRes, rank_res)      
        cat(sprintf(" --->>> %s / 51 : DONE.\n", RunNo ))
    }
    LungTestRes2 = LungTestRes
    LungTestRes2$sample_tissue = "LUNG"
    LungTestRes2$res_TS_COSINE = gsub("LUAD", "LUNG", LungTestRes2$res_TS_COSINE)
    LungTestRes2$res_TS_COSINE = gsub("SCLC", "LUNG", LungTestRes2$res_TS_COSINE)
    LungTestRes2$res_CL_COSINE = gsub("LUAD", "LUNG", LungTestRes2$res_CL_COSINE)
    LungTestRes2$res_CL_COSINE = gsub("SCLC", "LUNG", LungTestRes2$res_CL_COSINE)
    LungTestRes2$res_UMAP_CENT_EUCLIDEAN = gsub("LUAD", "LUNG", LungTestRes2$res_UMAP_CENT_EUCLIDEAN)
    LungTestRes2$res_UMAP_CENT_EUCLIDEAN = gsub("SCLC", "LUNG", LungTestRes2$res_UMAP_CENT_EUCLIDEAN)
    LungTestRes2$sample_tissue = factor(LungTestRes2$sample_tissue, levels=c(names(TISSUE_COLORS), "LUNG"))
    LungTestRes2$res_TS_COSINE = factor(LungTestRes2$res_TS_COSINE, levels=c(names(TISSUE_COLORS), "LUNG"))
    LungTestRes2$res_CL_COSINE = factor(LungTestRes2$res_CL_COSINE, levels=c(names(TISSUE_COLORS), "LUNG"))
    LungTestRes2$res_UMAP_CENT_EUCLIDEAN = factor(LungTestRes2$res_UMAP_CENT_EUCLIDEAN, levels=c(names(TISSUE_COLORS), "LUNG"))

    yardstick::accuracy(data=LungTestRes2, 
        truth=sample_tissue, estimate=res_TS_COSINE
    )
    yardstick::accuracy(data=LungTestRes2, 
        truth=sample_tissue, estimate=res_CL_COSINE
    )
    yardstick::accuracy(data=LungTestRes2, 
        truth=sample_tissue, estimate=res_UMAP_CENT_EUCLIDEAN
    )

    openxlsx::write.xlsx(LungTestRes, "TEST/CCLE_LUNG_CANCER_OTHER_TYPE_TEST_RESULT.xlsx")


#------------------------------------------------------------------------------#
    test_res$score_tissue
    test_res$score_cell
    test_res$plot2D
    plotly::ggplotly(test_res$plot2D)

    GTF = "/storage/home/kangsm/myDB/gtf/hg19/gencode.v22.annotation.gtf.gz"
    gtf_gr = rtracklayer::import(GTF)
    gtf = as.data.frame(gtf_gr)

    ID_MAP = gtf[,c("seqnames","gene_id","gene_name","transcript_id")] %>% unique()
    ID_MAP = ID_MAP %>% filter( !is.na(transcript_id))



    library(data.table)
    library(htmlwidgets)
    library(threejs)
    geo_brca_dir = "/storage/home/kangsm/analysis/gcx.rnaseq.based.clustering/GEO_TESTSET/GSE235167_BRCA"
    data_list    = list.files(geo_brca_dir, pattern = "sf.gz")
    
    GEO_BRCA_PDX_RESULT = data.frame()

    ST="BRCA_PDX"
    for( k in 1:length(data_list))
    {
        exprData = as.data.frame(fread(sprintf("%s/%s", geo_brca_dir, data_list[k] )))

        exprData = exprData %>% left_join(., ID_MAP %>% dplyr::rename(Name=transcript_id), by="Name" )
        exprGene = exprData %>% group_by(gene_id) %>% reframe( TPM = sum(TPM), READS=round(sum(NumReads),0) )
        exprGene$ens_id = sapply(exprGene$gene_id, function(y) unlist(strsplit(y, "\\."))[1])
        exprGene$gene = ginfo[match(exprGene$ens_id, ginfo$ens_geneid), "symbol"]

        exprGeneValue = exprGene$TPM
        names(exprGeneValue) = exprGene$gene

        TEST_EXPR_VALUES = exprGeneValue[class.marker.genes]
        TEST_EXPR_VALUES[is.na(TEST_EXPR_VALUES)] = 0

        TEST_ID = gsub("\\.trimmed.quant.sf.gz$", "", data_list[k])

        test_res = CLUSTER.TEST( 
            CCLE_EXPR_MATRIX = ccle.clustering.cell.lines.expr.mat,
            CCLE_MEDIAN      = ccle.clustering.cell.lines.expr.stat$tpm_median, 
            CCLE_MAD         = ccle.clustering.cell.lines.expr.stat$tpm_mad, 
            CCLE_INFO        = ccle.expr.clustering.cinfo, 
            MARKER_GENES     = class.marker.genes,
            TEST_EXPR_VALUES, 
            PointSize      = 5, 
            UMAP_METRIC    = "cosine",
            CORES          = 25,
            SCORE_CELL_TOP = NA,
            TISSUE_COLORS  = TISSUE_COLORS,
            NEIGHBORS      = NA, 
            ZSCORE         = TRUE
        )
        ggsave(test_res$plot2D, file = sprintf("TEST/GEO_TEST_%s.png", TEST_ID), width=8, height=8, unit='in', dpi=150, type='cairo')

        rank_res = data.frame(
            sample_id      = TEST_ID,
            sample_tissue  = ST,
            res_TS_COSINE = test_res$score_tissue %>% arrange(desc(cos)) %>% .$code %>% head(1),
            res_CL_COSINE = test_res$score_cell %>% arrange(desc(COS)) %>% .$code %>% head(1),
            res_UMAP_CENT_EUCLIDEAN = test_res$umap_centroid %>% arrange(euclidean_dist_rank) %>% .$code %>% head(1),
            res_UMAP_CENT_COSINE    = test_res$umap_centroid %>% arrange(cosine_dist_rank) %>% .$code %>% head(1)
        )
        
        GEO_BRCA_PDX_RESULT = rbind(GEO_BRCA_PDX_RESULT, rank_res)
    }

    openxlsx::write.xlsx(GEO_BRCA_PDX_RESULT, "TEST/BRCA_PDX_TEST_RESULT.xlsx")

    GEO_BRCA_PDX_RESULT$sample_tissue = "BRCA"
    GEO_BRCA_PDX_RESULT$sample_tissue = factor(GEO_BRCA_PDX_RESULT$sample_tissue, levels=names(TISSUE_COLORS))
    GEO_BRCA_PDX_RESULT$res_TS_SCORE = factor(GEO_BRCA_PDX_RESULT$res_TS_SCORE, levels=names(TISSUE_COLORS))
    GEO_BRCA_PDX_RESULT$res_CL_SCORE = factor(GEO_BRCA_PDX_RESULT$res_CL_SCORE, levels=names(TISSUE_COLORS))
    GEO_BRCA_PDX_RESULT$res_UMAP_CENT_EUCLIDEAN = factor(GEO_BRCA_PDX_RESULT$res_UMAP_CENT_EUCLIDEAN, levels=names(TISSUE_COLORS))
    GEO_BRCA_PDX_RESULT$res_UMAP_CENT_COSINE = factor(GEO_BRCA_PDX_RESULT$res_UMAP_CENT_COSINE, levels=names(TISSUE_COLORS))

    yardstick::accuracy(data=GEO_BRCA_PDX_RESULT, truth=sample_tissue, estimate=res_UMAP_CENT_COSINE)

    geo_brca$SAMPLE.TISSUE = "BRCA"
    geo_brca$SAMPLE.TISSUE = factor(geo_brca$SAMPLE.TISSUE, levels=names(TISSUE_COLORS))
    geo_brca$'TISSUE-based.SCORE' = factor(geo_brca$'TISSUE-based.SCORE', levels=names(TISSUE_COLORS))
    geo_brca$'CELL-based.SCORE' = factor(geo_brca$'CELL-based.SCORE', levels=names(TISSUE_COLORS))
    geo_brca$UMAP.CENTROID.DISTANCE = factor(geo_brca$UMAP.CENTROID.DISTANCE, levels=names(TISSUE_COLORS))

    yardstick::accuracy(data=geo_brca %>% filter( BRCA.TYPE == "ER+"), truth=SAMPLE.TISSUE, estimate='TISSUE-based.SCORE')
    yardstick::accuracy(data=geo_brca %>% filter( BRCA.TYPE == "ER+"), truth=SAMPLE.TISSUE, estimate='CELL-based.SCORE')
    yardstick::accuracy(data=geo_brca %>% filter( BRCA.TYPE == "ER+"), truth=SAMPLE.TISSUE, estimate='UMAP.CENTROID.DISTANCE')

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

    # LIHC PDX
    lihc_pdx_fpkm = as.data.frame(fread("/storage/home/kangsm/analysis/gcx.rnaseq.based.clustering/GEO_TESTSET/GSE130160_LUNG_PDX/GSE130160_TR_5253.genes.fpkm_table3.txt.gz"))

    pdx_id = c("TR_5253_001", "TR_5253_002","TR_5253_003","TR_5253_004","TR_5253_005","TR_5253_006")

    rownames(lihc_pdx_fpkm) = lihc_pdx_fpkm$transcript_id
    pdx_fpkm = lihc_pdx_fpkm[, pdx_id]

    # TPM= apply(RPKM, 2, function(x) x / sum(as.numeric(x)) * 10^6) %>% as.data.frame()
    pdx_tpm = apply(pdx_fpkm, 2, function(x) x / sum(as.numeric(x)) * 10^6)
    
    pdx_tpm = data.frame(pdx_tpm, tx_id=rownames(pdx_tpm), gene= lihc_pdx_fpkm[match(rownames(pdx_tpm), lihc_pdx_fpkm$transcript_id), "gene_name"] )

    pdx_tpm_rev = pdx_tpm %>% group_by( gene ) %>%
        reframe( 
            pdx01 = sum(TR_5253_001),
            pdx02 = sum(TR_5253_002),
            pdx03 = sum(TR_5253_003),
            pdx04 = sum(TR_5253_004),
            pdx05 = sum(TR_5253_005),
            pdx06 = sum(TR_5253_006)
        )
    pdx_tpm_rev = as.data.frame(pdx_tpm_rev)
    rownames(pdx_tpm_rev) = pdx_tpm_rev$gene

    PDX_TPM = as.matrix(pdx_tpm_rev[,-1])

    ST="LIHC_PDX"
    pdx_id_list = c("pdx01","pdx02","pdx03","pdx04","pdx05")

    GEO_LIHC_PDX_RESULT = data.frame()

    for( k in 1:length(pdx_id_list))
    {
        selCol = which(colnames(pdx_tpm_rev) == pdx_id_list[k] )
        exprGeneValue = PDX_TPM[, selCol]

        # exprGeneValue = exprGene$TPM
        # names(exprGeneValue) = exprGene$gene

        TEST_EXPR_VALUES = exprGeneValue[class.marker.genes]
        names(TEST_EXPR_VALUES) = class.marker.genes
        TEST_EXPR_VALUES[is.na(TEST_EXPR_VALUES)] = 0

        TEST_ID = pdx_id_list[k]

        test_res = CLUSTER.TEST( 
            CCLE_EXPR_MATRIX = ccle.clustering.cell.lines.expr.mat,
            CCLE_MEDIAN      = ccle.clustering.cell.lines.expr.stat$tpm_median, 
            CCLE_MAD         = ccle.clustering.cell.lines.expr.stat$tpm_mad, 
            CCLE_INFO        = ccle.expr.clustering.cinfo, 
            MARKER_GENES     = class.marker.genes,
            TEST_EXPR_VALUES, 
            PointSize      = 5, 
            UMAP_METRIC    = "cosine",
            CORES          = 25,
            SCORE_CELL_TOP = NA,
            TISSUE_COLORS  = TISSUE_COLORS,
            NEIGHBORS      = NA, 
            ZSCORE         = TRUE
        )

        #ggsave(test_res$plot2D, file = sprintf("TEST/LIHC_PDX_TEST_%s.png", TEST_ID), width=8, height=8, unit='in', dpi=150, type='cairo')

        rank_res = data.frame(
            sample_id      = TEST_ID,
            sample_tissue  = ST,
            res_TS_COSINE = test_res$score_tissue %>% arrange(desc(cos)) %>% .$code %>% head(1),
            res_CL_COSINE = test_res$score_cell %>% arrange(desc(COS)) %>% .$code %>% head(1),
            res_UMAP_CENT_EUCLIDEAN = test_res$umap_centroid %>% arrange(euclidean_dist_rank) %>% .$code %>% head(1),
            res_UMAP_CENT_COSINE    = test_res$umap_centroid %>% arrange(cosine_dist_rank) %>% .$code %>% head(1)
        )
        
        GEO_LIHC_PDX_RESULT = rbind(GEO_LIHC_PDX_RESULT, rank_res)
    }

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
    # LUNG PDX
    ID_MAP2 = ID_MAP[,c(2,3)] %>% unique()

    lung_pdx_list = c("GSE69405_GSM1702587_LUAD_PDX","GSE69405_GSM1702587_LungCellLine","LUNG_CANCER_PDX_01","LUNG_CANCER_PDX_02","LUNG_CANCER_PDX_03","LUNG_CANCER_PDX_04")

    GEO_LUNG_PDX_RESULT = data.frame()
    ST="LUNG_PDX"

    for( k in 1:length(lung_pdx_list))
    {
        TEST_ID = lung_pdx_list[k]

        rsem_res_file = sprintf("/storage/home/kangsm/analysis/gcx.rnaseq.based.clustering/GEO_TESTSET/LUNG_PDX/%s/quant/rsem/%s.genes.results", TEST_ID, TEST_ID)
        rsem_res = read.delim(rsem_res_file)
        
        rsem_res$gene = ID_MAP2[match(rsem_res$gene_id, ID_MAP2$gene_id), "gene_name"]

        exprGeneValue = rsem_res$TPM
        names(exprGeneValue) = rsem_res$gene

        TEST_EXPR_VALUES = exprGeneValue[class.marker.genes]
        TEST_EXPR_VALUES[is.na(TEST_EXPR_VALUES)] = 0

        test_res = CLUSTER.TEST( 
            CCLE_EXPR_MATRIX = ccle.clustering.cell.lines.expr.mat,
            CCLE_MEDIAN      = ccle.clustering.cell.lines.expr.stat$tpm_median, 
            CCLE_MAD         = ccle.clustering.cell.lines.expr.stat$tpm_mad, 
            CCLE_INFO        = ccle.expr.clustering.cinfo, 
            MARKER_GENES     = class.marker.genes,
            TEST_EXPR_VALUES, 
            PointSize      = 5, 
            UMAP_METRIC    = "cosine",
            CORES          = 25,
            SCORE_CELL_TOP = NA,
            TISSUE_COLORS  = TISSUE_COLORS,
            NEIGHBORS      = NA, 
            ZSCORE         = TRUE
        )

        # ggsave(test_res$plot2D, file = sprintf("TEST/GEO_LUNG_PDX_TEST_%s.png", TEST_ID), width=8, height=8, unit='in', dpi=150, type='cairo')

        # rank_res = data.frame(
        #     sample_id      = TEST_ID,
        #     sample_tissue  = ST,
        #     res_TS_COSINE = test_res$score_tissue %>% arrange(desc(cos)) %>% .$code %>% head(1),
        #     res_CL_COSINE = test_res$score_cell %>% arrange(desc(COS)) %>% .$code %>% head(1),
        #     res_UMAP_CENT_EUCLIDEAN = test_res$umap_centroid %>% arrange(euclidean_dist_rank) %>% .$code %>% head(1),
        #     res_UMAP_CENT_COSINE    = test_res$umap_centroid %>% arrange(cosine_dist_rank) %>% .$code %>% head(1)
        # )
        
        # GEO_LUNG_PDX_RESULT = rbind(GEO_LUNG_PDX_RESULT, rank_res)
        test_res$score_tissue %>% arrange(desc(cos))
        test_res$score_cell %>% arrange(desc(COS))
        test_res$umap_centroid %>% arrange(euclidean_dist_rank)

    }

    openxlsx::write.xlsx(GEO_LUNG_PDX_RESULT, "TEST/GEO_LUNG_PDX_TEST_RESULT.xlsx")

















