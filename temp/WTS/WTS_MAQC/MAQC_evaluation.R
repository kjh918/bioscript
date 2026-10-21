
setwd("/data/wts/WTS_MAQC/performance")

library(plyr)
library(dplyr)
library(Hmisc)
library(openxlsx)

## benchmark data
bmk_cq <- read.xlsx("RT_PCR_Cq.xlsx")

## gcx maqc data
gcx_data <- readRDS("/data/wts/RDS/RDS_DegAnalysis/WTS_MAQC.AID001.DEG.analysis.result.table.rds")
gcx_tpm  <- data.frame(
    MAQC_A_R1=readRDS("/data/wts/RDS/RDS_TPM/MAQC_A_rep1.TPM.rds"),
    MAQC_A_R2=readRDS("/data/wts/RDS/RDS_TPM/MAQC_A_rep2.TPM.rds"),
    MAQC_B_R1=readRDS("/data/wts/RDS/RDS_TPM/MAQC_B_rep1.TPM.rds"),
    MAQC_B_R2=readRDS("/data/wts/RDS/RDS_TPM/MAQC_B_rep2.TPM.rds")
)
rownames(gcx_tpm) <- sapply(rownames(gcx_tpm), function(id) unlist(strsplit(id, "\\."))[1])


###
res <- bmk_cq %>% dplyr::rename( gene_id=1 )
res <- data.frame(
    res,
    gcx_tpm[res$gene_id, ]
)
res <- res %>% mutate(GCX_MAQC_A=(MAQC_A_rep1+MAQC_A_rep2)/2, GCX_MAQC_B=(MAQC_B_rep1+MAQC_B_rep2)/2)
res <- res[which(!is.na(res$GCX_MAQC_A)), ]
res <- res[which(!is.na(res$GCX_MAQC_B)), ]
res$keep <- apply(res[,c("GCX_MAQC_A","GCX_MAQC_B")], 1, function(z) all(z > 0.5) )
table(res$keep)

RES <- res %>% filter( keep ) %>% mutate( GCX.MAQCA = GCX_MAQC_A, GCX.MAQCB = GCX_MAQC_B) %>% dplyr::select(c(1,2,3,10,11,12))

RES <- RES %>% mutate( ref.FC = log2(MAQCA/MAQCB), gcx.FC = log2(GCX.MAQCB/GCX.MAQCA) )

which(is.nan(RES$ref.FC))

RES <- RES[which(!is.nan(RES$ref.FC)),  ]
RES <- RES[which(!is.infinite(RES$ref.FC)),  ]

cor( x=RES$ref.FC, y=RES$gcx.FC, method="spearman")

summary(lm( gcx.FC ~ ref.FC, data = RES))

plot( x=RES$ref.FC, y=RES$gcx.FC, xlim=c(-1,1))

qpcr <- ggplot( RES, aes(x=ref.FC, y=gcx.FC) ) +
    geom_hline( yintercept=0, color="#2d2d2d", alpha=0.8, linewidth=0.7) +
    geom_vline( xintercept=0, color="#2d2d2d", alpha=0.8, linewidth=0.7) +
    geom_point( pch=21, size=3, fill="#D79B00", alpha=0.7, color="#FFFFFF" ) + 
    coord_cartesian( xlim=c(-1.2, 1.2)) +
    labs( x = "MAQC-B/MAQC-A qRT-PCR FoldChange", y="MAQC-B/MAQC-A RNAseq FoldChange") +
    scale_x_continuous( breaks = c(0), label = c("0") ) +
    scale_y_continuous( breaks = c(0), label = c("0") ) +
    theme(
        panel.background = element_rect(fill = '#FFFFFF'),
        panel.grid.major = element_line(color = '#cccccc', linetype = 'dotted'),
        axis.line    = element_blank(),
        axis.ticks   = element_blank(),
        axis.text    = element_text(family="Helvetica", face="bold", color="#2d2d2d", size=16, hjust=0.5, vjust=0.5),
        axis.title   = element_text(family="Helvetica", face="bold", color="#2d2d2d", size=18, hjust=0.5, vjust=0.5),
        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )

ggsave(qpcr, file = "MAQC_FoldChange_Correlation_Scatter_green.png", width=6, height=6, unit='in', dpi=150)

cnv %>% filter( sample_info %in% c("hCRC","hLC","hRCC","hGC") ) %>%
    ggplot(aes(x=sample_origin, y=V2, color=sample_origin) ) + 
    geom_boxplot() +
    geom_point() + 
    geom_line(aes(group=sample_group)) +
    facet_wrap(~ sample_info)


res <- res %>% mutate(
    ref.FC1 = MAQCB/MAQCA,
    ref.FC2 = log2(res$MAQCB+1) - log2(res$MAQCA+1),
    gcx.FC1 = GCX_MAQC_B/GCX_MAQC_A,
    gcx.FC2 = log2(res$GCX_MAQC_B+1) - log2(res$GCX_MAQC_A+1)
)

plot(res$ref.FC2, y=res$gcx.FC2)

cor(res$ref.FC2, y=res$gcx.FC2, na.rm=T)
plot(res$ref.FC2, y=res$gcx.FC2)

res$ref.FC <- round(log2(res$MAQCB+1) - log2(res$MAQCA+1),6) 
res$gcx.FC <- round(log2(res$GCX_MAQC_B+1) - log2(res$GCX_MAQC_A+1),6) 

res$ref.FC1 <- MAQCB/MAQCA

round(log2(res$MAQCB+1) - log2(res$MAQCA+1),6) 
res$gcx.FC <- round(log2(res$GCX_MAQC_B+1) - log2(res$GCX_MAQC_A+1),6) 


plot( x = res$ref.FC, y = res$gcx.FC )


res$ref.FC <- apply(bmk_cq[,2:3], 1, function(y)
{
    if( y[1] == 0 & y[2] > 0 ){ 
        fc = round(y[2]/1 ,6)
    }else if( y[2] == 0 & y[1] > 0 ){ 
        fc = round(1 / y[1], 6)
    }else if( y[1] == 0 & y[2] == 0 ){
        fc = 0
    }else{
        fc = round(y[2] / y[1], 6)
    }
    return(fc)
})


%>% mutate( ref_FC = (MAQCB+0.1)/(MAQCA+0.1))

bmk_cq <- bmk_cq %>% mutate( fc = MAQCB/MAQCA )



bmk_cq$gcx_logfc <- 2^(gcx_data[match( bmk_cq[,1], gcx_data$ens_geneid), "log2_foldchange" ])