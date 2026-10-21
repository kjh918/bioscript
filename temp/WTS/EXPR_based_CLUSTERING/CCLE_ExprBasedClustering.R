
#---| PACKAGES |-----------------------------------------------------------------------------------3
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("umap"))
    suppressPackageStartupMessages(library("openxlsx"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("scale"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    ccle_data_dir = "/storage/home/kangsm/myDB/ccle_rnaseq"
    res_dir       = "/storage/home/kangsm/analysis/gcx.rnaseq.based.clustering"

#---| LOAD CELL LINE INFO |------------------------------------------------------------------------#

    cinfo <- read.xlsx(sprintf("%s/Cell_Lines_Info_Master_v1.0_20220321.xlsx", ccle_data_dir))
    load(sprintf("%s/CCLE_RNAseq_v2_TPM_Matrix.Rdata", ccle_data_dir))          # TpmMat , REF 

    ccle_cinfo <- cinfo %>% filter( cellid %in% colnames(TpmMat) )
    ccle_cinfo <- ccle_cinfo %>% arrange( tissue, disease_class3, disease )

    ccle_cinfo_stat = ccle_cinfo %>% group_by( tissue, disease_class3, disease ) %>% reframe( Cells = length(cellid) )


    write.xlsx(list(CINFO=ccle_cinfo, STAT=ccle_cinfo_stat), sprintf("%s/CCLE_Cell_Info.xlsx", ccle_data_dir))

    cinfo = read.xlsx(sprintf("%s/CCLE_Cell_Info_rev.xlsx", ccle_data_dir))

    cinfo_stat = cinfo %>% group_by(tissue, disease) %>% reframe( cells = length(cellid) ) %>% 
        arrange( tissue, desc(cells) )

#--------------------------------------------------------------------------------------------------#

#---| COLOR CODE |---------------------------------------------------------------------------------#
##
    #             tissue                                            disease cells
    # 1           Bladder                                  Bladder carcinoma    25    #00b2b2
    # 2             Blood                                Plasma cell myeloma    28    #B85450
    # 3             Blood                       Adult acute myeloid leukemia    17    #F8CECC
    # 4             Brain                                       Glioblastoma    38    #E78AC3
    # 5             Brain                                        Astrocytoma    13    #E78AC3
    # 6            Breast                            Ductal breast carcinoma    16    #82B366
    # 7            Breast                                   Breast carcinoma    15    #82B366
    # 8            Breast Invasive ductal carcinoma, not otherwise specified    14    #82B366
    # 9            Breast                              Breast adenocarcinoma     9    #82B366
    # 10       Colorectal                               Colon adenocarcinoma    34    #D79B00
    # 11       Colorectal                               Cecum adenocarcinoma     8    #D79B00
    # 12       Colorectal                              Rectal adenocarcinoma     7    #D79B00
    # 13        Esophagus                 Esophageal squamous cell carcinoma    23    #8DD3C7
    # 14           Kidney                               Renal cell carcinoma    21    #FB8072
    # 15           Kidney                    Clear cell renal cell carcinoma     8    #FB8072
    # 16            Liver                           Hepatocellular carcinoma    22    #DAE8FC
    # 17             Lung                                Lung adenocarcinoma    72    #6C8EBF
    # 18             Lung                          Small cell lung carcinoma    48    #006EAF
    # 19             Lung                       Squamous cell lung carcinoma    26    #6C8EBF
    # 20             Lung                          Large cell lung carcinoma    14    #6C8EBF
    # 21             Lung                      Non-small cell lung carcinoma    11    #6C8EBF
    # 22            Lymph                      Diffuse large B-cell lymphoma    16    #db9e6c
    # 23            Lymph                     Anaplastic large cell lymphoma    10    #e7be9d
    # 24            Lymph                                   Hodgkin lymphoma    10    #f3dece
    # 25            Ovary           High grade ovarian serous adenocarcinoma    19    #9673A6
    # 26         Pancreas                          Pancreatic adenocarcinoma    18    #66C2A5
    # 27         Pancreas                   Pancreatic ductal adenocarcinoma    15    #66C2A5
    # 28 Peripheral Nerve                                      Neuroblastoma    18    #8DA0CB
    # 29         Prostate                                 Prostate carcinoma     6    #1eb0a1
    # 30             Skin                                           Melanoma    38    #FFED6F
    # 31             Skin                                 Cutaneous melanoma    11    #FFED6F
    # 32      Soft_Tissue                                      Ewing sarcoma    13    #59b200
    # 33      Soft_Tissue                          Extrarenal rhabdoid tumor    10    #323299
    # 34      Soft_Tissue                                       Osteosarcoma    10    #999999
    # 35          Stomach                             Gastric adenocarcinoma    13    #E1D5E7
    # 36           Uterus                         Endometrial adenocarcinoma    15    #FFE6CC
##
TISSUE_COLORS = c(
    "Bladder carcinoma"                        = "#00b2b2", "Plasma cell myeloma"           = "#B85450", 
    "Adult acute myeloid leukemia"             = "#F8CECC", "Glioblastoma"                  = "#E78AC3", 
    "Astrocytoma"                              = "#E78AC3", "Ductal breast carcinoma"       = "#82B366",
    "Breast carcinoma"                         = "#82B366", "Invasive ductal carcinoma, not otherwise specified"="#82B366",
    "Breast adenocarcinoma"                    = "#82B366", "Colon adenocarcinoma"          = "#D79B00", 
    "Cecum adenocarcinoma"                     = "#D79B00", "Rectal adenocarcinoma"         = "#D79B00", 
    "Esophageal squamous cell carcinoma"       = "#8DD3C7", "Renal cell carcinoma"          = "#FB8072", 
    "Clear cell renal cell carcinoma"          = "#FB8072", "Hepatocellular carcinoma"      = "#DAE8FC", 
    "Lung adenocarcinoma"                      = "#6C8EBF", "Small cell lung carcinoma"     = "#006EAF", 
    "Squamous cell lung carcinoma"             = "#6C8EBF", "Large cell lung carcinoma"     = "#6C8EBF", 
    "Non-small cell lung carcinoma"            = "#6C8EBF", "Diffuse large B-cell lymphoma" = "#db9e6c", 
    "Anaplastic large cell lymphoma"           = "#e7be9d", "Hodgkin lymphoma"              = "#f3dece", 
    "High grade ovarian serous adenocarcinoma" = "#9673A6", "Pancreatic adenocarcinoma"     = "#66C2A5", 
    "Pancreatic ductal adenocarcinoma"         = "#66C2A5", "Neuroblastoma"                 = "#8DA0CB", 
    "Prostate carcinoma"                       = "#1eb0a1", "Melanoma"                      = "#FFED6F", 
    "Cutaneous melanoma"                       = "#FFED6F", "Ewing sarcoma"                 = "#59b200",
    "Extrarenal rhabdoid tumor"                = "#323299", "Osteosarcoma"                  = "#999999", 
    "Gastric adenocarcinoma"                   = "#E1D5E7", "Endometrial adenocarcinoma"    = "#FFE6CC"
)
cinfo$COLORS = TISSUE_COLORS[ cinfo$disease ]

#--------------------------------------------------------------------------------------------------#

#---| CCLE Only Clustering |-----------------------------------------------------------------------#

    load(sprintf("%s/CCLE_RNAseq_v2_TPM_Matrix.Rdata", ccle_data_dir)) # TpmMat
    
    CINFO = cinfo %>% filter( tissue %nin% c("Blood","Lymph"))
    SOLID_TUMOR_TPM = TpmMat[,CINFO$cellid]
    
    SOLID_TUMOR_TPM_FILTER = SOLID_TUMOR_TPM[which(apply(SOLID_TUMOR_TPM, 1, function(y) mad(y) ) != 0), ]
    SOLID_TUMOR_TPM_Z = t(apply(SOLID_TUMOR_TPM_FILTER, 1, function(z) (z-median(z))/mad(z)))

    indexing = function(x){ sapply(unique(x), function(z) list(z)) }
    gidx = indexing(rownames(SOLID_TUMOR_TPM_FILTER))    
    
    TissueName = "Lung adenocarcinoma"

    TIDX = indexing(CINFO$disease)

    TISSUE_MARKER_LIST = lapply( TIDX, function(TissueName) 
    {
        message( sprintf("TISSUE = %s", TissueName))
        cid1 = CINFO %>% filter( disease == TissueName ) %>% .$cellid
        cid2 = CINFO %>% filter( disease != TissueName ) %>% .$cellid
        
        tissue_markers = mclapply( gidx,  function(y) 
        {
            gv1 = SOLID_TUMOR_TPM_FILTER[y, cid1]
            gv2 = SOLID_TUMOR_TPM_FILTER[y, cid2]
            res = data.frame( 
                T_test = t.test(gv1,gv2)$p.value,
                W_test = wilcox.test(gv1,gv2)$p.value,
                M1 = median(gv1),
                M2 = median(gv2)
            )
            return(res)
        }, mc.cores= 60)
        tissue_markers = ldply(tissue_markers, .id="gene") %>% mutate( gene = as.character(gene) ) %>% 
            mutate( median_gap = M2 - M1 )
    
        MARKER_GENE = unique(rbind(
            tissue_markers %>% filter( W_test < 0.001 , T_test < 0.01 ) %>% filter( abs(median_gap) >= 50 ),
            tissue_markers %>% filter( T_test < 0.001 , W_test < 0.01, M1 < 0 ) %>% filter( abs(median_gap) > 1 ),
            tissue_markers %>% filter( T_test < 0.001 , W_test < 0.01, M2 < 0 ) %>% filter( abs(median_gap) > 1 )
        )) %>% .$gene

        return(MARKER_GENE)
        #return(tissue_markers[1,])
    })    

    GENES = rownames(SOLID_TUMOR_TPM_Z)

    MarkerCheck = mclapply(indexing(GENES), function(y)
    {
        chk = lapply(TISSUE_MARKER_LIST, function(z) ifelse(y %in% z, 1, 0))

        chk2 = data.frame(matrix(unlist(chk), nrow=1))
        colnames(chk2) = names(chk)
        #chk2$GENE = y
        return(chk2)
    }, mc.cores=60)
    MarkerCheck = ldply(MarkerCheck, .id="GENE")

    UniqMarker = as.character(MarkerCheck[which(apply(MarkerCheck[,-1], 1, function(y) sum(y)) == 1 ), ]$GENE)

    TISSUE_MARKER_LIST2 = lapply(TISSUE_MARKER_LIST, function(qq) qq[ qq %in% UniqMarker])

    TISSUE_MARKER_LIST2[c("Colon adenocarcinoma","Hepatocellular carcinoma","Lung adenocarcinoma","")]


    MARKER="PON3"
    par(mfrow=c(1,2))
    boxplot(list( SOLID_TUMOR_TPM_Z[MARKER, cid1], SOLID_TUMOR_TPM_Z[MARKER, cid2]) )
    boxplot(list( SOLID_TUMOR_TPM_FILTER[MARKER, cid1], SOLID_TUMOR_TPM_FILTER[MARKER, cid2]) )


    load("/storage/home/kangsm/myDB/ProteinAtlas/HPA_TissueSpecificExpressionMarkers.Rdata") #tis_marker_list
    tis_marker_list_rev = lapply(tis_marker_list, function(ts)
    {
        ts %>% filter( RNA.tissue.specificity == "Tissue enriched") %>% .$Gene
    })

    TISSUES = c("Lung","Bladder","Brain","Breast","Colorectal","Kidney","Liver","Ovary","Pancreas","Prostate","Skin","Stomach","Uterus")
    EXPR_MATRIX = TpmMat[, CINFO$cellid]

    VARIANCE_CUTOFF = 0.6 ## --> means use top 40% variance genes
 
    TISSUE_COLORS = c(
        "Bladder carcinoma"                        = "#00b2b2", 
        "Plasma cell myeloma"                      = "#B85450", 
        "Adult acute myeloid leukemia"             = "#F8CECC", 
        "Glioblastoma"                             = "#FFED6F", 
        "Astrocytoma"                              = "#E78AC3", 
        "Ductal breast carcinoma"                  = "#FFBF00",
        "Breast carcinoma"                         = "#FFBF00", 
        "Invasive ductal carcinoma, not otherwise specified"="#82B366",
        "Breast adenocarcinoma"                    = "#82B366", 
        "Colon adenocarcinoma"                     = "#FF7F0E", 
        "Esophageal squamous cell carcinoma"       = "#8DD3C7", 
        "Renal cell carcinoma"                     = "#FB8072", 
        "Clear cell renal cell carcinoma"          = "#FB8072", 
        "Hepatocellular carcinoma"                 = "#E41A1C", 
        "Lung adenocarcinoma"                      = "#1BA1E2", 
        "Small cell lung carcinoma"                = "#006EAF", 
        "Squamous cell lung carcinoma"             = "#8DA0CB", 
        "Large cell lung carcinoma"                = "#306493", 
        "Non-small cell lung carcinoma"            = "#8DD3C7", 
        "Diffuse large B-cell lymphoma"            = "#db9e6c", 
        "Anaplastic large cell lymphoma"           = "#e7be9d", 
        "Hodgkin lymphoma"                         = "#f3dece", 
        "High grade ovarian serous adenocarcinoma" = "#000059", 
        "Pancreatic adenocarcinoma"                = "#66C2A5", 
        "Pancreatic ductal adenocarcinoma"         = "#66C2A5", 
        "Neuroblastoma"                            = "#8DA0CB", 
        "Prostate carcinoma"                       = "#FFBF00", 
        "Melanoma"                                 = "#4DAF4A", 
        "Cutaneous melanoma"                       = "#4DAF4A", 
        "Ewing sarcoma"                            = "#59b200",
        "Extrarenal rhabdoid tumor"                = "#323299", 
        "Osteosarcoma"                             = "#999999", 
        "Gastric adenocarcinoma"                   = "#5900b2", 
        "Endometrial adenocarcinoma"               = "#FFE6CC"
    )

    c("Bladder","Blood","Soft_Tissue","Brain","Breast","Colorectal","Esophagus","Kidney","Liver","Lung","Lymph"
      "Ovary","Pancreas","Peripheral Nerve","Prostate","Skin","Stomach","Uterus")

    
    Draw.CCLE.UMAP.Clustering <- function( ExprMatrix, SelectTissues, VarianceCutoff, CellInfo, PointSize=8, NEIGHBORS=NA, TISSUE_COLORS, ZSCORE=TRUE )
    {
        require(ggplot2)
        require(umap)
        require(dplyr)
        require(threejs)
        # tissue filtering
        message(" >> tissue filtering...")
        if( all(is.na(SelectTissues)) )
        { 
            cell.info = CellInfo  
        }else{
            cell.info = CellInfo %>% filter( tissue %in% SelectTissues )
        }
        exprM = ExprMatrix[, cell.info$cellid]
        # gene variance
        message(" >> gene filtering...")
        exprM.mad = apply(exprM, 1, function(y) mad(y) )
        exprM.var = exprM[which(exprM.mad != 0 ), ]
        exprM.mad.nonZero = rev(sort(exprM.mad[ exprM.mad != 0 ]))
        exprM.mad.filtered = exprM.mad.nonZero[ exprM.mad.nonZero >= quantile(exprM.mad.nonZero, VarianceCutoff) ]
        # gene filtering and expression values transformation
        if( ZSCORE )
        {
            INPUT = t(apply(exprM.var[names(exprM.mad.filtered), ], 1, function(z) (z-median(z))/mad(z)))
        }else{
            INPUT = exprM.var[names(exprM.mad.filtered), ]
        }
        # umap config
        umapConfig              = umap.defaults
        umapConfig$n_components = 2
        if( is.na(NEIGHBORS) )
        { umapConfig$n_neighbors = as.integer(ncol(INPUT)^(1/2)) }else{ umapConfig$n_neighbors = NEIGHBORS } 
        umapConfig$metric       = "cosine"
        # run umap
        message(" >> Run UMAP 2D...")
        umapResult = umap(t(INPUT), config=umapConfig)
        umapPlotData = as.data.frame(umapResult$layout[colnames(INPUT),])
        umapPlotData = data.frame(
            cellid=rownames(umapPlotData),
            X=umapPlotData[,1],
            Y=umapPlotData[,2],
            cell.info[match(rownames(umapPlotData), cell.info$cellid), c("cell_name","tissue","disease","COLORS")]
        )
        # plots
        message(" >> create plot2D...")
        umapPlots2D <- ggplot(umapPlotData, aes(x=X, y=Y, fill=disease, text=cell_name)) + geom_point(col="white", pch=21, size=PointSize, alpha=0.9, stroke=0.5) +
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
        umapConfig$metric       = "cosine"
        umapResult2 = umap(t(INPUT), config=umapConfig)
        umapPlotData2 = as.data.frame(umapResult2$layout[colnames(INPUT),])
        umapPlotData2 = data.frame(
            cellid=rownames(umapPlotData2),
            X=umapPlotData2[,1],
            Y=umapPlotData2[,2],
            Z=umapPlotData2[,3],
            cell.info[match(rownames(umapPlotData2), cell.info$cellid), c("cell_name","tissue","disease","COLORS")]
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

    DS = c(
        "Colon adenocarcinoma",
        "Renal cell carcinoma",
        "Breast carcinoma",
        "Hepatocellular carcinoma",
        "Lung adenocarcinoma",
        #"Squamous cell lung carcinoma",
        "Small cell lung carcinoma",
        "High grade ovarian serous adenocarcinoma",
        "Pancreatic adenocarcinoma",
        "Melanoma",
        "Gastric adenocarcinoma"
        #"Prostate carcinoma" 
    )
    cinfo2 = CINFO %>% filter( disease %in% DS )

    rmList = c("SLR 20","NUGC-2","GCIY","HCC1833","NCI-H196","SW1271","SBC-5","NCI-H841","DMS 454","DMS 53","DU4475")
    cinfo2 = cinfo2 %>% filter( cell_name %nin% rmList )


    MARKERSET1 = Reduce(union, TISSUE_MARKER_LIST2[DS])
    MARKERSET2 = Reduce(union, tis_marker_list_rev[c(1,2,4,5,6,7,8,9,10,12)])

    ccle.umap = Draw.CCLE.UMAP.Clustering(
        ExprMatrix     = SOLID_TUMOR_TPM[unique(c(MARKERSET1,MARKERSET2)), ],
        SelectTissues  = NA,
        VarianceCutoff = 0.2,
        CellInfo       = cinfo2,
        PointSize      = 5,
        NEIGHBORS      = NA,
        TISSUE_COLORS  = TISSUE_COLORS,
        ZSCORE = TRUE
    )
    #ccle.umap$plot2D
    ggplotly(ccle.umap$plot2D)


    rmList =c("SLR 20")


    ccle.umap$plot2D

    ccle.umap$plot


            umapPlotData = data.frame(
            cellid=rownames(umapPlotData),
            X=umapPlotData[,1],
            Y=umapPlotData[,2],
            Z=umapPlotData[,3],
            cinfo[match(rownames(umapPlotData), cinfo$cellid), c("cell_name","tissue","disease","COLORS")]
        )


        scatterplot3js( x = umapPlotData$X, y = umapPlotData$Y, z = umapPlotData$Z, color = umapPlotData$COLORS, labels = umapPlotData$disease, axis=TRUE, num.ticks = c(0,0,0), stroke = '#e5e5e5' )

