
library(Cairo)
grDevices::X11.options(type='cairo')
options(device='x11')
options(stringsAsFactors=FALSE) 
library(plyr)
library(dplyr)
library(ggplot2)
library(CopyNumberPlots)
INDEXING = function(X) { sapply(unique(X), function(y) list(y)) }

### FUNCTIONS ###------------------------------------------------------------------------------------------------------------------------------------#
# Draw CNV Plot -----------------------------------------------------------------------------------#
  draw.CNVplot = function( CNV_DATA, CHR="all", byCHR=FALSE, TITLE=NULL, MAX_Y=NA, baseZero = TRUE, refVer='hg19') 
  {
    if( CHR == 'all' ){ Chromosome = c(seq(1,22,1), "X","Y") }else{ Chromosome = CHR }
    if( refVer == "hg38" ){ CNV_DATA$chr = gsub("^chr", "", CNV_DATA$chr) }
    if( CHR == "all" ){ 
      chr_label = Chromosome 
      chr_lbl_size = 4
      X_lbl = "chromosome"
      segSize = 0.2
    }else{ 
      chr_label = paste0("chr",Chromosome) 
      chr_lbl_size = 6
      X_lbl = ""
      segSize = 0.5
    }
    if( byCHR && CHR %in% "all") 
    {
      inputDat = CNV_DATA %>% filter( chr %in% Chromosome ) %>% arrange(as.numeric(X)) %>% group_by(chr) %>% mutate( index=1:length(chr)) %>%
        mutate(chr = factor(chr, levels=c(seq(1,22,1), "X","Y")))
      chrLabel = paste0("chr", c(seq(1,22,1), "X","Y"))
      names(chrLabel) = c(seq(1,22,1), "X","Y")
    }else{
      inputDat = CNV_DATA %>% filter( chr %in% Chromosome ) %>% arrange(as.numeric(X)) %>% mutate(index=1:length(chr)) %>% 
        mutate(chr = factor(chr, levels=c(seq(1,22,1), "X","Y")))
      chr_label = intersect(chr_label, inputDat$chr)
    }
    annot_X = inputDat %>% group_by(as.numeric(chr2)) %>% summarise( chr_start=head(index, 1)) %>% .$chr_start
    if( is.na(MAX_Y) ){ maxY = max(abs(inputDat$segment)) }else{ maxY = MAX_Y }
    if( maxY < 4 ){ maxY = 4 }
    if( baseZero ) { minY = 0 }else{ minY = -maxY }

    cnvPlot = ggplot(inputDat) + 
      geom_point(  aes( x = index, y = segment), color="#D9D9D9", size=segSize, alpha=0.5) + 
      geom_segment(aes( x = index, y = cnv, xend = index+1, yend = cnv), color="#931b3b", linewidth=1) +
      labs( x = X_lbl, y = "CNV" ) +
      theme(
        axis.text.x = element_blank(),
        panel.background = element_rect(fill='white'), 
        panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth=0.4),
        axis.line.y  = element_line(colour = "#6f6d6d"),
        axis.text.y  = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 5)),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 5)),
        axis.ticks.x = element_blank(),
        plot.title   = element_text(face="bold", size=15)
      )
    if( byCHR && CHR == "all" )
    { 
      cnvPlot + facet_wrap( vars(chr), scales = "free_x", labeller = as_labeller(chrLabel[intersect(names(chrLabel), unique(inputDat$chr))]) ) + ggtitle(TITLE)
    }else{
      cnvPlot + coord_cartesian(xlim=c(1,nrow(inputDat)+1), ylim = c(minY,maxY) ) + 
        geom_vline( xintercept = annot_X, color="#69aadb", linetype='dashed', linewidth=0.2, alpha=0.7) +
        annotate("text", x = annot_X, y = minY, label = chr_label, hjust=-0.25, vjust=-0.25, size=chr_lbl_size) + ggtitle(TITLE)
    }
  }
# READ ADTEx DATA and CONVERT ---------------------------------------------------------------------#
  Read.ADTEx.Result = function( ADTExResgionFile, refVer=c("hg19","hg38") ) 
  {
    ADTExSeg = read.delim( ADTExResgionFile )
    if( refVer == "hg38" ){ ADTExSeg$chr = gsub("chr", "", ADTExSeg$chr) }
    ADTExSeg$chr2 = ADTExSeg$chr ; ADTExSeg$chr2 = gsub("X", "23", ADTExSeg$chr2) ; ADTExSeg$chr2 = gsub("Y", "24", ADTExSeg$chr2)
    ADTExSeg = ADTExSeg %>% dplyr::rename( start = exon_start, end = exon_end )
    maxW = max(nchar(ADTExSeg$start))
    ADTExSeg$start = formatC(ADTExSeg$start, digits=0, width=maxW, format='f', flag='0')
    ADTExSeg$end   = formatC(ADTExSeg$end,   digits=0, width=maxW, format='f', flag='0')
    ADTExSeg2 = ADTExSeg %>% mutate( segment = 2^(rationormalized_after_smoothing), cnv = 2^(seg_mean), X = paste(chr2, start, sep="."), X2 = paste(chr2, end, sep=".") ) %>%
      mutate( X = as.numeric(X), X2 = as.numeric(X2))
    return(ADTExSeg2)
  }
# READ CNVkit DATA and CONVERT --------------------------------------------------------------------#
  Read.CNVkit.Result = function( cnvkitRegionFile, cnvkitCallFile, refVer=c("hg19","hg38") )
  {
      resRegion = read.delim(cnvkitRegionFile)
      resCall   = read.delim(cnvkitCallFile)
      resRegion$cnv = NA
      for( k in 1:nrow(resCall) ){ resRegion[which( resRegion$chromosome == resCall[k,"chromosome"] & resRegion$start >= resCall[k,"start"] & resRegion$end <= resCall[k, "end"] ), "cnv" ] = resCall[k, "cn"] }
      cnvSeg = resRegion %>% filter( gene != "Antitarget" ) %>% select(c("chromosome","start","end","log2","cnv")) %>% dplyr::rename(chr=chromosome)
      maxW = max(nchar(cnvSeg$start))
      if( refVer == "hg38" ){ cnvSeg$chr = gsub("chr", "", cnvSeg$chr) }
      cnvSeg$chr2 = cnvSeg$chr
      cnvSeg$chr2 = gsub("X", "23", cnvSeg$chr2) ; cnvSeg$chr2 = gsub("Y", "24", cnvSeg$chr2)
      cnvSeg$start = formatC(cnvSeg$start, digits=0, width=maxW, format='f', flag='0')
      cnvSeg$end   = formatC(cnvSeg$end,   digits=0, width=maxW, format='f', flag='0')
      cnvSeg2 = cnvSeg %>% mutate( segment = 2^(log2)*2 , X = paste(chr2, start, sep="."), X2 = paste(chr2, end, sep=".") ) %>%
        mutate( X = as.numeric(X), X2 = as.numeric(X2)) 
      return(cnvSeg2)
  }
# READ EXCAVATOR2 DATA and CONVERT ----------------------------------------------------------------#
  Read.EXCAVATOR2.Result = function( xcvtrRegionFile, xcvtrCallFile, refVer=c("hg19","hg38") )
  {
    resRegion = read.delim( xcvtrRegionFile )
    resCall   = read.delim( xcvtrCallFile )
    resRegion$cnv = NA
    resRegion = resRegion %>% dplyr::rename( chr=Chromosome, start=Start, end=End, log2=SegMean )
    for( k in 1:nrow(resCall) ){ resRegion[which( resRegion$chr == resCall[k,"Chromosome"] & resRegion$start >= resCall[k,"Start"] & resRegion$end <= resCall[k, "End"] ), "cnv" ] = resCall[k, "CN"] }
    resRegion[which(is.na(resRegion$cnv)), 'cnv'] = 2
    cnvSeg = resRegion %>% select(c("chr","start","end","log2","cnv"))
    maxW = max(nchar(cnvSeg$start))
    if( refVer == "hg38" ){ cnvSeg$chr = gsub("chr", "", cnvSeg$chr) }
    cnvSeg$chr2 = cnvSeg$chr
    cnvSeg$chr2 = gsub("X", "23", cnvSeg$chr2) ; cnvSeg$chr2 = gsub("Y", "24", cnvSeg$chr2)
    cnvSeg$start = formatC(cnvSeg$start, digits=0, width=maxW, format='f', flag='0')
    cnvSeg$end   = formatC(cnvSeg$end,   digits=0, width=maxW, format='f', flag='0')
    cnvSeg2 = cnvSeg %>% mutate( segment = 2^(log2)*2 , X = paste(chr2, start, sep="."), X2 = paste(chr2, end, sep=".") ) %>%
      mutate( X = as.numeric(X), X2 = as.numeric(X2))
  }
# READ XHMM DATA and CONVERT ----------------------------------------------------------------------#
  Read.XHMM.RD = function( XHMM_FILE, what = double() )
  {
    READ_MAT_FILE_PREFIX = ""
      READ_MAT_FILE_SUFFIX = paste0(" ", XHMM_FILE)
      if (regexpr("\\.gz$", XHMM_FILE) > 0) {
          READ_MAT_FILE_PREFIX = paste("gzip -cd ", XHMM_FILE, " | ", sep="")
          READ_MAT_FILE_SUFFIX = ""
      }

    getDims = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' '{if (NR>1) print NF}'", READ_MAT_FILE_SUFFIX, " | sort | uniq -c | awk '{print $1,$2}'", sep="")
    dimsVec = system(getDims, intern=TRUE)
    if (length(dimsVec) != 1) { stop(paste("Cannot read jagged matrix: ", XHMM_FILE, sep="")) }
    
    rows_cols = as.numeric(strsplit(dimsVec, "\\s+")[[1]])
    rows = rows_cols[1]
    cols = rows_cols[2] - 1

    n = rows * cols
    readMat = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {$1=\"\"; gsub(\"^\"OFS\"+\", \"\", $0); print $0;}'", READ_MAT_FILE_SUFFIX, sep="")

    if (log2(n) < 31) {
      con <- pipe(readMat)
      m = matrix(scan(con, what=what, n=n, skip=1, quiet=TRUE, sep="\t"), rows, cols, byrow = TRUE)
      close(con)
    }
    else {
      m = data.frame()
      for (r in 1:rows) {
        con <- pipe(readMat)
        m = rbind(m, scan(con, what=what, n=cols, skip=(1+(r-1)), quiet=TRUE, sep="\t"))
        close(con)
      }
    }

    readColNames = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {if (NR == 1) {print $0;", ifelse(READ_MAT_FILE_PREFIX == "", "exit", ""), "}}'", READ_MAT_FILE_SUFFIX, sep="")
    con <- pipe(readColNames)
          COLNAMES = scan(con, what="", quiet=TRUE, sep="\t")
          if (length(COLNAMES) > cols) {
              COLNAMES = COLNAMES[-1]
          }
          if (length(COLNAMES) != cols) {
      stop(paste("Number of header column names [", length(COLNAMES), "] does not match number of columns [", cols, "] in matrix: ", XHMM_FILE, sep=""))            
          }
    colnames(m) = COLNAMES
    close(con)

    readRowNames = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {if (NR > 1) print $1}'", READ_MAT_FILE_SUFFIX, sep="")
    con <- pipe(readRowNames)
    rownames(m) = scan(con, what="", quiet=TRUE, sep="\t")
    close(con)

    return(m)
  }
  Read.XHMM.Result = function( xhmmZscoreFile, xhmCallFile, sampleID, refVer=c("hg19","hg38") )
  {
    xhmmZscore     = Read.XHMM.RD( xhmmZscoreFile , what=double() ) 
    xhmmResZ       = as.data.frame(t(xhmmZscore))
    xhmmResZ$chr   = sapply(rownames(xhmmResZ), function(y) unlist(strsplit(y, "\\:"))[1])
    xhmmResZ$start = as.numeric(sapply(rownames(xhmmResZ), function(y) unlist(strsplit(unlist(strsplit(y, "\\:"))[2], "\\-"))[1]))
    xhmmResZ$end   = as.numeric(sapply(rownames(xhmmResZ), function(y) unlist(strsplit(unlist(strsplit(y, "\\:"))[2], "\\-"))[2]))
    SAMPLE_ID = sampleID
    xhmmCall = read.delim(xhmCallFile)
    xhmmCall = xhmmCall %>% filter(SAMPLE == SAMPLE_ID) 
    xhmmCall$START = as.numeric(sapply(xhmmCall$INTERVAL, function(y) unlist(strsplit(unlist(strsplit(y, "\\:"))[2], "\\-"))[1]  ))
    xhmmCall$END   = as.numeric(sapply(xhmmCall$INTERVAL, function(y) unlist(strsplit(unlist(strsplit(y, "\\:"))[2], "\\-"))[2]  ))
    xhmmResZ$cnv = NA
    for( k in 1:nrow(xhmmCall) ) {
      xhmmResZ[which(xhmmResZ$chr == xhmmCall[k,"CHR"] & xhmmResZ$start >= xhmmCall[k,"START"] & xhmmResZ$end <= xhmmCall[k,"END"]), "cnv" ] = xhmmCall[k, "CNV"]
    }
    xhmmRes = xhmmResZ %>% dplyr::rename(segment=which(colnames(xhmmResZ)==SAMPLE_ID) ) %>% select(c("chr","start","end","segment","cnv"))
    xhmmRes[which(is.na(xhmmRes$cnv)), "cnv"]   = 0
    xhmmRes[which(xhmmRes$cnv == "DEL"), "cnv"] = -3
    xhmmRes[which(xhmmRes$cnv == "DUP"), "cnv"] = 3
    xhmmRes$cnv = as.numeric(xhmmRes$cnv)

    maxW = max(nchar(xhmmRes$start))
    xhmmRes$chr2 = xhmmRes$chr
    xhmmRes$chr2 = gsub("X", "23", xhmmRes$chr2) ; xhmmRes$chr2 = gsub("Y", "24", xhmmRes$chr2)
    if( refVer == "hg38" ){ cnvSeg$chr2 = gsub("chr", "", cnvSeg$chr2) }
    xhmmRes$start = formatC(xhmmRes$start, digits=0, width=maxW, format='f', flag='0')
    xhmmRes$end   = formatC(xhmmRes$end,   digits=0, width=maxW, format='f', flag='0')
    xhmmRes2 = xhmmRes %>% mutate( X = paste(chr2, start, sep="."), X2 = paste(chr2, end, sep=".") ) %>%
      mutate( X = as.numeric(X), X2 = as.numeric(X2))
    return(xhmmRes2)
  }
# READ cn.MOPS DATA and CONVERT
  Read.cnMOPS.Result = function( cnMOPS_RESULT_LIST )
  {
    resSeg  = as.data.frame(individualCall( cnMOPS_RESULT_LIST ))
    resCall = as.data.frame(  segmentation( cnMOPS_RESULT_LIST ))
    resSeg$cnv = NA
    for( k in 1:nrow(resCall) ) {
      resSeg[which( resSeg$seqnames == resCall[k,"seqnames"] & resSeg$start >= resCall[k,"start"] & resSeg$end <= resCall[k, "end"] ), "cnv" ] = resCall[k, "CN"]  
    }
    cnmopsRes = resSeg %>% dplyr::rename(chr=1, log2=6) %>% select(c("chr","start","end","log2","cnv")) %>% mutate( cnv = as.numeric(gsub("CN","",cnv)) )
    maxW = max(nchar(cnmopsRes$start))
    cnmopsRes$chr2 = cnmopsRes$chr
    cnmopsRes$chr2 = gsub("X", "23", cnmopsRes$chr2) ; cnmopsRes$chr2 = gsub("Y", "24", cnmopsRes$chr2)
    cnmopsRes$start = formatC(cnmopsRes$start, digits=0, width=maxW, format='f', flag='0')
    cnmopsRes$end   = formatC(cnmopsRes$end,   digits=0, width=maxW, format='f', flag='0')
    cnmopsRes2 = cnmopsRes %>% mutate( segment = 2^(log2)*2, X = paste(chr2, start, sep="."), X2 = paste(chr2, end, sep=".") ) %>% mutate( X = as.numeric(X), X2 = as.numeric(X2)) 
    return(cnmopsRes2)
  }
  cnMOPS.Result2List = function( cnMOPS_BATCH_RESULT )
  {
    resSeg  = as.data.frame(individualCall( cnMOPS_BATCH_RESULT ))
    resCall = as.data.frame(  segmentation( cnMOPS_BATCH_RESULT ))
    SampleList = unique(resCall$sampleName)
    sidx = INDEXING( SampleList )
    ResList = mclapply( sidx, function(y) 
    {
      seg.y  = resSeg[,c(1:5, which(colnames(resSeg)==y))]
      call.y = resCall %>% filter( sampleName == y ) 
      seg.y$cnv = NA
      for( k in 1:nrow(call.y) ) {
        seg.y[which( seg.y$seqnames == call.y[k,"seqnames"] & seg.y$start >= call.y[k,"start"] & seg.y$end <= call.y[k, "end"] ), "cnv" ] = call.y[k, "CN"]  
      }
      res.y = seg.y %>% dplyr::rename(chr=1, log2=6) %>% select(c("chr","start","end","log2","cnv")) %>% mutate( cnv = as.numeric(gsub("CN","",cnv)) )
      maxW = max(nchar(res.y$start))
      res.y$chr2 = res.y$chr
      res.y$chr2 = gsub("X", "23", res.y$chr2) ; res.y$chr2 = gsub("Y", "24", res.y$chr2)
      res.y$start = formatC(res.y$start, digits=0, width=maxW, format='f', flag='0')
      res.y$end   = formatC(res.y$end,   digits=0, width=maxW, format='f', flag='0')
      res2.y = res.y %>% mutate( segment = 2^(log2)*2, X = paste(chr2, start, sep="."), X2 = paste(chr2, end, sep=".") ) %>% mutate( X = as.numeric(X), X2 = as.numeric(X2)) 
      return(res2.y)
    }, mc.cores = length(sidx))
    return(ResList)
  }

#--------------------------------------------------------------------------------------------------#

### Sample IDs
SID = c("PON01_01","PON01_02","PON01_03","PON01_04","PON01_05")
sidx <- INDEXING(SID)

### ADTEx ###
  ADTEX_BASE_DIR <- "/storage/home/kangsm/analysis/cnv_test/adtex"
  ADTEx_RES = mclapply(sidx, function(y) {
    ADTExResgionFile = sprintf("%s/%s/cnv.result", ADTEX_BASE_DIR, y)
    ADTExRes = Read.ADTEx.Result( ADTExResgionFile )
    return(ADTExRes)
  }, mc.cores = length(sidx) )

  #  draw.CNVplot(CNV_DATA = ADTEx_RES[[1]], CHR = "all", byCHR = FALSE, TITLE="ADTEx")

### CNVkit ###
  CNVKIT_BASE_DIR <- "/storage/home/kangsm/analysis/cnv_test/cnvkit"

  CNVkitPairedRES = mclapply( sidx, function(y) {
    cnvkitRegionFile = sprintf("%s/paired/%s_T.analysis.ready.cnr",     CNVKIT_BASE_DIR, y)
    cnvkitCallFile   = sprintf("%s/paired/%s_T.analysis.ready.call.cns",CNVKIT_BASE_DIR, y)
    cnvkitRes = Read.CNVkit.Result( cnvkitRegionFile, cnvkitCallFile )
    return(cnvkitRes)
  }, mc.cores = length(sidx))

  CNVkitPONRES = mclapply( sidx, function(y) {
    cnvkitRegionFile = sprintf("%s/pooled_normal/%s/%s_T.analysis.ready.cnr",     CNVKIT_BASE_DIR, y,y)
    cnvkitCallFile   = sprintf("%s/pooled_normal/%s/%s_T.analysis.ready.call.cns",CNVKIT_BASE_DIR, y,y)
    cnvkitRes = Read.CNVkit.Result( cnvkitRegionFile, cnvkitCallFile )
    return(cnvkitRes)
  }, mc.cores = length(sidx))

  CNVkitBatchRES = mclapply( sidx, function(y) {
    cnvkitRegionFile = sprintf("%s/batch/%s_T.analysis.ready.cnr",     CNVKIT_BASE_DIR, y)
    cnvkitCallFile   = sprintf("%s/batch/%s_T.analysis.ready.call.cns",CNVKIT_BASE_DIR, y)
    cnvkitRes = Read.CNVkit.Result( cnvkitRegionFile, cnvkitCallFile )
    return(cnvkitRes)
  }, mc.cores = length(sidx))

  CNVkitSingleRES = mclapply( sidx, function(y) {
    cnvkitRegionFile = sprintf("%s/single/%s/%s_T.analysis.ready.cnr",     CNVKIT_BASE_DIR, y,y)
    cnvkitCallFile   = sprintf("%s/single/%s/%s_T.analysis.ready.call.cns",CNVKIT_BASE_DIR, y,y)
    cnvkitRes = Read.CNVkit.Result( cnvkitRegionFile, cnvkitCallFile )
    return(cnvkitRes)
  }, mc.cores = length(sidx))



  #  draw.CNVplot(CNV_DATA = CNVkitRES[[2]], CHR = "all", byCHR = FALSE, TITLE="CNVkit")

### EXCAVATOR2 ###
  EXCAVATOR2_BASE_DIR <- "/storage/home/kangsm/analysis/cnv_test/excavator2/PON01_w10k"
  EXCAVATOR2PoNRES = mclapply( sidx, function(y) {
    xcvtrRegionFile = sprintf("%s/result.pooling.%s_T/Results/%s_T/HSLMResults_%s_T.txt", EXCAVATOR2_BASE_DIR, y, y, y )
    xcvtrCallFile   = sprintf("%s/result.pooling.%s_T/Results/%s_T/FastCallResults_%s_T.txt", EXCAVATOR2_BASE_DIR, y, y, y )
    excvator2Res = Read.EXCAVATOR2.Result( xcvtrRegionFile, xcvtrCallFile )
    return(excvator2Res)
  }, mc.cores = length(sidx))

  EXCAVATOR2PairedRES = mclapply( sidx, function(y) {
    xcvtrRegionFile = sprintf("%s/result.paired/Results/%s_T/HSLMResults_%s_T.txt", EXCAVATOR2_BASE_DIR, y, y, y )
    xcvtrCallFile   = sprintf("%s/result.paired/Results/%s_T/FastCallResults_%s_T.txt", EXCAVATOR2_BASE_DIR, y, y, y )
    excvator2Res = Read.EXCAVATOR2.Result( xcvtrRegionFile, xcvtrCallFile )
    return(excvator2Res)
  }, mc.cores = length(sidx))

  #  draw.CNVplot(CNV_DATA = EXCAVATOR2_RES[[1]], CHR = "all", byCHR = FALSE, TITLE="EXCAVATOR2", MAX_Y=NA)

### XHMM ###
  XHMM_BASE_DIR = "/storage/home/kangsm/analysis/cnv_test/xhmm"
  sampleIdList = c("coad.tumor.01","coad.tumor.02","coad.tumor.03","coad.tumor.04","coad.tumor.05")
  names(sampleIdList) = SID
  XHMM_RES = mclapply( sidx, function(y) {
    xhmmZscoreFile = sprintf("%s/%s_pooled.normal/pcaNormalized.zScores.RD.txt", XHMM_BASE_DIR,y)
    xhmCallFile    = sprintf("%s/%s_pooled.normal/CNV.analysisResult.xcnv", XHMM_BASE_DIR,y)
    xhmmRes = Read.XHMM.Result( xhmmZscoreFile, xhmCallFile, sampleID = sampleIdList[y])
    return(xhmmRes)
  }, mc.cores = length(sidx))

  #  draw.CNVplot(CNV_DATA = XHMM_RES[[2]], CHR = "all", byCHR = FALSE, TITLE="XHMM", MAX_Y=NA, baseZero = FALSE)


### cn.MOPS
  library(cn.mops)
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/NA12878.batch.exomeCNMOPS.Rdata"           )
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/PanelofNormal_segmentMatrix.Rdata"         )
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/PON01_TumorOnlyBatch_exomeCNMOPS.Rdata"    )
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/PON01.NTpaired.exomeCNMOPS.list.Rdata"     )
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/PON01.NTpaired.referenceCNMOPS.list.Rdata" )
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/PON01.Pooled.All.exomeCNMOPS.Rdata"        )
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/singleTumor.PoN.exomeCNMOPS.list.Rdata"    )
  load("/storage/home/kangsm/analysis/cnv_test/cn.mops/singleTumor.PoN.referenceCNMOPS.list.Rdata")



  # Paired reference-method --> PON01.NTpaired.referenceCNMOPS.list
  cnMOPSPairedRefRes = mclapply( PON01.NTpaired.referenceCNMOPS.list, function(y) { Read.cnMOPS.Result(y) }, mc.cores = 5 )
  # Paired exome.cn.mops-method --> PON01.NTpaired.exomeCNMOPS.list
  cnMOPSPairedCnmopsRes = mclapply( PON01.NTpaired.exomeCNMOPS.list[-3], function(y) { Read.cnMOPS.Result(y) }, mc.cores = 5 )
  # Tumor + PON + reference-method --> 
  cnMOPSSinglePoNRefRes = mclapply( singleTumor.PoN.referenceCNMOPS.list, function(y) { Read.cnMOPS.Result(y) }, mc.cores = 5 )
  # Tumor + PON + exome.cn.mops-method --> 
  cnMOPSSinglePoNCnmopsRes = mclapply( singleTumor.PoN.exomeCNMOPS.list, function(y) { Read.cnMOPS.Result(y) }, mc.cores = 5 )
  # T-only batch
  cnMOPSTumorBatchCnmopsRes = cnMOPS.Result2List( PON01_TumorOnlyBatch_exomeCNMOPS )
  # Pooled All 
  cnMOPSPoolAllCnmopsRes = cnMOPS.Result2List( PON01.Pooled.All.exomeCNMOPS )

  K=3
  cnvP1 = draw.CNVplot(CNV_DATA = cnMOPSPairedRefRes[[K]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : Paired REF", MAX_Y=4,  baseZero = TRUE)
  #cnvP2 = draw.CNVplot(CNV_DATA = cnMOPSPairedCnmopsRes[[K-1]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : Paired MOPS", MAX_Y=4,  baseZero = TRUE)
  cnvP3 = draw.CNVplot(CNV_DATA = cnMOPSSinglePoNRefRes[[K]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : Single REF", MAX_Y=4,  baseZero = TRUE)
  cnvP4 = draw.CNVplot(CNV_DATA = cnMOPSSinglePoNCnmopsRes[[K]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : Single MOPS", MAX_Y=4,  baseZero = TRUE)
  cnvP5 = draw.CNVplot(CNV_DATA = cnMOPSTumorBatchCnmopsRes[[K]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : Batch", MAX_Y=4,  baseZero = TRUE)
  cnvP6 = draw.CNVplot(CNV_DATA = cnMOPSPoolAllCnmopsRes[[K]], CHR = "all", byCHR = FALSE, TITLE="cn.MOPS : Pool", MAX_Y=4,  baseZero = TRUE)
  cowplot::plot_grid(cnvP1, cnvP3, cnvP4, cnvP5, cnvP6, ncol=1)

  K=3
  cnvP1 = draw.CNVplot(CNV_DATA = cnMOPSPairedRefRes[[K]], CHR = c(17,18,19), byCHR = FALSE, TITLE="cn.MOPS : Paired REF", MAX_Y=4,  baseZero = TRUE)
  #cnvP2 = draw.CNVplot(CNV_DATA = cnMOPSPairedCnmopsRes[[K-1]], CHR = c(17,18,19), byCHR = FALSE, TITLE="cn.MOPS : Paired MOPS", MAX_Y=4,  baseZero = TRUE)
  cnvP3 = draw.CNVplot(CNV_DATA = cnMOPSSinglePoNRefRes[[K]], CHR = c(17,18,19), byCHR = FALSE, TITLE="cn.MOPS : Single REF", MAX_Y=4,  baseZero = TRUE)
  cnvP4 = draw.CNVplot(CNV_DATA = cnMOPSSinglePoNCnmopsRes[[K]], CHR = c(17,18,19), byCHR = FALSE, TITLE="cn.MOPS : Single MOPS", MAX_Y=4,  baseZero = TRUE)
  cnvP5 = draw.CNVplot(CNV_DATA = cnMOPSTumorBatchCnmopsRes[[K]], CHR = c(17,18,19), byCHR = FALSE, TITLE="cn.MOPS : Batch", MAX_Y=4,  baseZero = TRUE)
  cnvP6 = draw.CNVplot(CNV_DATA = cnMOPSPoolAllCnmopsRes[[K]], CHR = c(17,18,19), byCHR = FALSE, TITLE="cn.MOPS : Pool", MAX_Y=4,  baseZero = TRUE)
  cowplot::plot_grid(cnvP1, cnvP3, cnvP4, cnvP5, cnvP6, ncol=1)

### MERGED PLOTS ### 
  K=5
  cnvP1 = draw.CNVplot(CNV_DATA = ADTEx_RES[[K]],      CHR = "all", byCHR = FALSE, TITLE="ADTEx",      MAX_Y=4,  baseZero = TRUE)
  cnvP2 = draw.CNVplot(CNV_DATA = CNVkitRES[[K]],      CHR = "all", byCHR = FALSE, TITLE="CNVkit",     MAX_Y=6,  baseZero = TRUE)
  cnvP3 = draw.CNVplot(CNV_DATA = EXCAVATOR2_RES[[K]], CHR = "all", byCHR = FALSE, TITLE="EXCAVATOR2", MAX_Y=4,  baseZero = TRUE)
  cnvP4 = draw.CNVplot(CNV_DATA = XHMM_RES[[K]],       CHR = "all", byCHR = FALSE, TITLE="XHMM",       MAX_Y=NA, baseZero = FALSE)
  cowplot::plot_grid(cnvP1, cnvP2, cnvP3, cnvP4, ncol=1)
 
  K=3
  cnvP1 = draw.CNVplot(CNV_DATA = CNVkitPairedRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : PAIRED",          MAX_Y=6,  baseZero = TRUE)
  cnvP2 = draw.CNVplot(CNV_DATA = CNVkitPONRES[[K]],    CHR = "all", byCHR = FALSE, TITLE="CNVkit : PANEL OF NORMAL", MAX_Y=6,  baseZero = TRUE)
  cnvP3 = draw.CNVplot(CNV_DATA = CNVkitBatchRES[[K]],  CHR = "all", byCHR = FALSE, TITLE="CNVkit : T-ONLY BATCH",    MAX_Y=20,  baseZero = TRUE)
  cnvP4 = draw.CNVplot(CNV_DATA = CNVkitSingleRES[[K]], CHR = "all", byCHR = FALSE, TITLE="CNVkit : T-ONLY SINGLE",   MAX_Y=20,  baseZero = TRUE)
  cowplot::plot_grid(cnvP1, cnvP2, cnvP3, cnvP4, ncol=1)
 




