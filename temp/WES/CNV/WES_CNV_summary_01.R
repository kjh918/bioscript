
### DATASET 1

# Paired Analysis Only, 
# ==> ADTEx_RES, CNVkitPairedRES, EXCAVATOR2PairedRES, cnMOPSPairedRefRes, cnMOPSPairedCnmopsRes

K = 'PON01_05'
cnvP1 = draw.CNVplot(CNV_DATA = ADTEx_RES[[K]], CHR = "all", byCHR = FALSE, TITLE="ADTEx : PAIRED",      MAX_Y=4,  baseZero = TRUE)
cnvP2 = draw.CNVplot(CNV_DATA = CNVkitPairedRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : PAIRED",          MAX_Y=6,  baseZero = TRUE)
cnvP3 = draw.CNVplot(CNV_DATA = EXCAVATOR2_RES[[K]], CHR = "all", byCHR = FALSE, TITLE="EXCAVATOR2 : PAIRED", MAX_Y=4,  baseZero = TRUE)
cnvP4 = draw.CNVplot(CNV_DATA = cnMOPSPairedRefRes[[K]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : PAIRED_REF", MAX_Y=4,  baseZero = TRUE)
cnvP5 = draw.CNVplot(CNV_DATA = cnMOPSPairedCnmopsRes[[K]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : PAIRED", MAX_Y=4,  baseZero = TRUE)

CNVPLOT = cowplot::plot_grid(cnvP1, cnvP2, cnvP3, cnvP4, cnvP5, ncol=1)
CNVPLOT
 
# 9 x 2.5
ggsave(CNVPLOT, file = sprintf('%s_paired_res_CNVplot.png',K), width=10, height=12.5, unit='in', dpi=150, type='cairo')



CNVPLOT = cowplot::plot_grid(cnvP1, cnvP2, cnvP3, cnvP4, ncol=1)
ggsave(CNVPLOT, file = sprintf('%s_paired_res_CNVplot.png',K), width=10, height=10, unit='in', dpi=150, type='cairo')

# Paired - PON 비교 
# CNVkitPairedRES, CNVkitPONRES, EXCAVATOR2PairedRES, EXCAVATOR2PoNRES, XHMM_RES, cnMOPSSinglePoNRefRes

K = 'PON01_05'
cnvP1 = draw.CNVplot(CNV_DATA = CNVkitPairedRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : PAIRED",      MAX_Y=6,  baseZero = TRUE)
cnvP2 = draw.CNVplot(CNV_DATA = CNVkitPONRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : PoN",          MAX_Y=6,  baseZero = TRUE)
cnvP3 = draw.CNVplot(CNV_DATA = EXCAVATOR2PairedRES[[K]], CHR = "all", byCHR = FALSE, TITLE="EXCAVATOR2 : PAIRED", MAX_Y=4,  baseZero = TRUE)
cnvP4 = draw.CNVplot(CNV_DATA = EXCAVATOR2PoNRES[[K]], CHR = "all", byCHR = FALSE, TITLE="EXCAVATOR2 : PoN", MAX_Y=4,  baseZero = TRUE)
cnvP5 = draw.CNVplot(CNV_DATA = XHMM_RES[[K]], CHR = "all", byCHR = FALSE, TITLE="XHMM : PoN", MAX_Y=NA, baseZero = FALSE)
cnvP6 = draw.CNVplot(CNV_DATA = cnMOPSSinglePoNRefRes[[sprintf("%s_T",K)]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : PoN_REF", MAX_Y=4,  baseZero = TRUE)

CNVPLOT = cowplot::plot_grid(cnvP1, cnvP2, cnvP3, cnvP4, cnvP5, cnvP6, ncol=1)
CNVPLOT
ggsave(CNVPLOT, file = sprintf('%s_paired_n_PoN_res_CNVplot.png',K), width=10, height=15, unit='in', dpi=150, type='cairo')

# CNVkit Mode별 차이
K = 'PON01_02'
for( K in c("PON01_01","PON01_02","PON01_03","PON01_04","PON01_05") )
{
cnvP1 = draw.CNVplot(CNV_DATA = CNVkitPairedRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : PAIRED",      MAX_Y=6,  baseZero = TRUE)
cnvP2 = draw.CNVplot(CNV_DATA = CNVkitPONRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : PoN",          MAX_Y=6,  baseZero = TRUE)
cnvP3 = draw.CNVplot(CNV_DATA = CNVkitBatchRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : T-only Batch",      MAX_Y=8,  baseZero = TRUE)
cnvP4 = draw.CNVplot(CNV_DATA = CNVkitSingleRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : T-only Single",          MAX_Y=8,  baseZero = TRUE)

CNVPLOT = cowplot::plot_grid(cnvP1, cnvP2, cnvP3, cnvP4, ncol=1)
CNVPLOT
ggsave(CNVPLOT, file = sprintf('%s_CNVkit_diff.modes_CNVplot.png',K), width=10, height=10, unit='in', dpi=150, type='cairo')
}


### DATASET 2
SID = c("NA12878_01","NA12878_02")
sidx <- INDEXING(SID)

  CNVKIT_BASE_DIR <- "/storage/home/kangsm/analysis/cnv_test/cnvkit"
  CNVkitPONRES = mclapply( sidx, function(y) {
    cnvkitRegionFile = sprintf("%s/pooled_normal/%s/%s.analysis.ready.cnr",     CNVKIT_BASE_DIR, y,y)
    cnvkitCallFile   = sprintf("%s/pooled_normal/%s/%s.analysis.ready.call.cns",CNVKIT_BASE_DIR, y,y)
    cnvkitRes = Read.CNVkit.Result( cnvkitRegionFile, cnvkitCallFile )
    return(cnvkitRes)
  }, mc.cores = length(sidx))

  CNVkitBatchRES = mclapply( sidx, function(y) {
    cnvkitRegionFile = sprintf("%s/batch/%s.analysis.ready.cnr",     CNVKIT_BASE_DIR, y)
    cnvkitCallFile   = sprintf("%s/batch/%s.analysis.ready.call.cns",CNVKIT_BASE_DIR, y)
    cnvkitRes = Read.CNVkit.Result( cnvkitRegionFile, cnvkitCallFile )
    return(cnvkitRes)
  }, mc.cores = length(sidx))

  CNVkitSingleRES = mclapply( sidx, function(y) {
    cnvkitRegionFile = sprintf("%s/single/%s/%s.analysis.ready.cnr",     CNVKIT_BASE_DIR, y,y)
    cnvkitCallFile   = sprintf("%s/single/%s/%s.analysis.ready.call.cns",CNVKIT_BASE_DIR, y,y)
    cnvkitRes = Read.CNVkit.Result( cnvkitRegionFile, cnvkitCallFile )
    return(cnvkitRes)
  }, mc.cores = length(sidx))

  EXCAVATOR2_BASE_DIR <- "/storage/home/kangsm/analysis/cnv_test/excavator2/NA12878_w10k"
  EXCAVATOR2PoNRES = mclapply( sidx, function(y) {
    xcvtrRegionFile = sprintf("%s/result.pooling.%s/Results/%s/HSLMResults_%s.txt", EXCAVATOR2_BASE_DIR, y, y, y )
    xcvtrCallFile   = sprintf("%s/result.pooling.%s/Results/%s/FastCallResults_%s.txt", EXCAVATOR2_BASE_DIR, y, y, y )
    excvator2Res = Read.EXCAVATOR2.Result( xcvtrRegionFile, xcvtrCallFile )
    return(excvator2Res)
  }, mc.cores = length(sidx))

  sampleIdList = c("23WES01_02","23WES01_03")
  names(sampleIdList) = SID
  XHMM_RES = mclapply( sidx, function(y) {
    xhmmZscoreFile = sprintf("%s/%s_pooled.normal/pcaNormalized.zScores.RD.txt", XHMM_BASE_DIR,y)
    xhmCallFile    = sprintf("%s/%s_pooled.normal/CNV.analysisResult.xcnv", XHMM_BASE_DIR,y)
    xhmmRes = Read.XHMM.Result( xhmmZscoreFile, xhmCallFile, sampleID = sampleIdList[y])
    return(xhmmRes)
  }, mc.cores = length(sidx))


K = 'NA12878_02'
cnvP1 = draw.CNVplot(CNV_DATA = CNVkitPONRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : PoN",      MAX_Y=6,  baseZero = TRUE)
cnvP2 = draw.CNVplot(CNV_DATA = CNVkitBatchRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : T-only Batch",          MAX_Y=8,  baseZero = TRUE)
cnvP3 = draw.CNVplot(CNV_DATA = CNVkitSingleRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : T-only Single",          MAX_Y=8,  baseZero = TRUE)
cnvP4 = draw.CNVplot(CNV_DATA = EXCAVATOR2PoNRES[[K]], CHR = "all", byCHR = FALSE, TITLE="EXCAVATOR2 : PoN", MAX_Y=4,  baseZero = TRUE)
cnvP5 = draw.CNVplot(CNV_DATA = XHMM_RES[[K]], CHR = "all", byCHR = FALSE, TITLE="XHMM : PoN", MAX_Y=NA, baseZero = FALSE)

CNVPLOT = cowplot::plot_grid(cnvP1, cnvP2, cnvP3, cnvP4, cnvP5, ncol=1)
CNVPLOT
ggsave(CNVPLOT, file = sprintf('%s_PoN_res_CNVplot.png',K), width=10, height=12.5, unit='in', dpi=150, type='cairo')



