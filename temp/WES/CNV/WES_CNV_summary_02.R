
library(Cairo)
grDevices::X11.options(type='cairo')
options(device='x11')
options(stringsAsFactors=FALSE) 
library(plyr)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(CopyNumberPlots)
INDEXING = function(X) { sapply(unique(X), function(y) list(y)) }

setwd("/storage/home/kangsm/analysis/cnv_test/wes_cnv_results_02_new")

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
    ADTExSeg$chr2 = ADTExSeg$chr ; ADTExSeg$chr2 = gsub("X", "23", ADTExSeg$chr2) ; ADTExSeg$chr2 = gsub("Y", "24", ADTExSeg$chr2)
    if( refVer == "hg38" ){ cnvSeg$chr2 = gsub("chr", "", cnvSeg$chr2) }
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
      cnvSeg$chr2 = cnvSeg$chr
      cnvSeg$chr2 = gsub("X", "23", cnvSeg$chr2) ; cnvSeg$chr2 = gsub("Y", "24", cnvSeg$chr2)
      if( refVer == "hg38" ){ cnvSeg$chr2 = gsub("chr", "", cnvSeg$chr2) }
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
    cnvSeg$chr2 = cnvSeg$chr
    cnvSeg$chr2 = gsub("X", "23", cnvSeg$chr2) ; cnvSeg$chr2 = gsub("Y", "24", cnvSeg$chr2)
    if( refVer == "hg38" ){ cnvSeg$chr2 = gsub("chr", "", cnvSeg$chr2) }
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

#--------------------------------------------------------------------------------------------------#


### HCC1395
BASE_DIR = "/storage/home/kangsm/analysis/cnv_test/wes_cnv_results_02_new"
sampleID = "HCC1395"
RES_DIR = sprintf("%s/%s", BASE_DIR, sampleID)

resList = list.files(RES_DIR)
resList = resList[ resList %nin% c("XHMM_coverage.files", "excavator2_run")]

setwd(RES_DIR)
ADTEx.Paired.Res = Read.ADTEx.Result( ADTExResgionFile = "ADTEx_paired/cnv.result", refVer="hg38" )
CNVkit.Paired.Res = Read.CNVkit.Result(
    cnvkitRegionFile = sprintf("CNVkit_paired/%s_T.analysis.ready.cnr",      sampleID),
    cnvkitCallFile   = sprintf("CNVkit_paired/%s_T.analysis.ready.call.cns", sampleID), 
    refVer="hg38" 
)
CNVkit.PoN.Res = Read.CNVkit.Result(
    cnvkitRegionFile = sprintf("CNVkit_pon/%s_T.analysis.ready.cnr",      sampleID),
    cnvkitCallFile   = sprintf("CNVkit_pon/%s_T.analysis.ready.call.cns", sampleID), 
    refVer="hg38" 
)
EXCAVATOR2.Paired.Res = Read.EXCAVATOR2.Result(
    xcvtrRegionFile = sprintf("EXCAVATOR2_paired/Results/%s/HSLMResults_%s.txt", sampleID, sampleID), 
    xcvtrCallFile   = sprintf("EXCAVATOR2_paired/Results/%s/FastCallResults_%s.txt", sampleID, sampleID), 
    refVer = "hg38"
)
EXCAVATOR2.PoN.Res = Read.EXCAVATOR2.Result(
    xcvtrRegionFile = sprintf("EXCAVATOR2_pon/Results/%s/HSLMResults_%s.txt", sampleID, sampleID), 
    xcvtrCallFile   = sprintf("EXCAVATOR2_pon/Results/%s/FastCallResults_%s.txt", sampleID, sampleID), 
    refVer = "hg38"
)

cnvPlot1 = draw.CNVplot( ADTEx.Paired.Res, CHR="all", byCHR=FALSE, TITLE="ADTEx (PAIRED)", MAX_Y=6, baseZero = TRUE, refVer='hg38')
cnvPlot2 = draw.CNVplot( CNVkit.Paired.Res, CHR="all", byCHR=FALSE, TITLE="CNVkit (PAIRED)", MAX_Y=6, baseZero = TRUE, refVer='hg38')
cnvPlot3 = draw.CNVplot( CNVkit.PoN.Res, CHR="all", byCHR=FALSE, TITLE="CNVkit (PoN)", MAX_Y=6, baseZero = TRUE, refVer='hg38')
cnvPlot4 = draw.CNVplot( EXCAVATOR2.PoN.Res, CHR="all", byCHR=FALSE, TITLE="EXCAVATOR2 (PoN)", MAX_Y=6, baseZero = TRUE, refVer='hg38')

cowplot::plot_grid(cnvPlot1, cnvPlot2, cnvPlot3, cnvPlot4, ncol=1 )


ANS = read.delim("meta")

ADTEx.Call = read.delim(sprintf("ADTEx_paired/%s.cnv", sampleID))
CNVkit.Call = read.delim(sprintf("CNVkit_paired/%s_T.analysis.ready.cnr", sampleID))

  Read.CNVkit.Result2 = function( cnvkitRegionFile, cnvkitCallFile, answerFile, refVer=c("hg19","hg38") )
  {
      resRegion = read.delim(cnvkitRegionFile)
      resCall   = read.delim(cnvkitCallFile)
      ans       = read.delim(answerFile)
      resRegion$cnv = NA
      resRegion$cnv2 = NA
      resRegion$ans_log2 = NA
      for( k in 1:nrow(resCall) ){ resRegion[which( resRegion$chromosome == resCall[k,"chromosome"] & resRegion$start >= resCall[k,"start"] & resRegion$end <= resCall[k, "end"] ), "cnv" ] = resCall[k, "cn"] }
      for( k in 1:nrow(ans) ){ resRegion[which( resRegion$chromosome == ans[k,"chromosome"] & resRegion$start >= ans[k,"start"] & resRegion$end <= ans[k, "end"] ), "cnv2" ] = ans[k, "cn"] }
      for( k in 1:nrow(ans) ){ resRegion[which( resRegion$chromosome == ans[k,"chromosome"] & resRegion$start >= ans[k,"start"] & resRegion$end <= ans[k, "end"] ), "ans_log2" ] = ans[k, "log2"] }
      cnvSeg = resRegion %>% filter( gene != "Antitarget" ) %>% select(c("chromosome","start","end","log2","cnv","cnv2", "ans_log2")) %>% dplyr::rename(chr=chromosome)
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

CNVkit.Paired.Res2 = Read.CNVkit.Result2(
    cnvkitRegionFile = sprintf("CNVkit_paired/%s_T.analysis.ready.cnr",      sampleID),
    cnvkitCallFile   = sprintf("CNVkit_paired/%s_T.analysis.ready.call.cns", sampleID), 
    answerFile       = "/storage/home/kangsm/analysis/cnv_test/wes_cnv_results_02_new/meta/seqc_hcc1395_cnvkit_trueset.txt",
    refVer="hg38" 
)

CNVkit.PoN.Res2 = Read.CNVkit.Result2(
    cnvkitRegionFile = sprintf("CNVkit_pon/%s_T.analysis.ready.cnr",      sampleID),
    cnvkitCallFile   = sprintf("CNVkit_pon/%s_T.analysis.ready.call.cns", sampleID), 
    answerFile       = "/storage/home/kangsm/analysis/cnv_test/wes_cnv_results_02_new/meta/seqc_hcc1395_cnvkit_trueset.txt",
    refVer="hg38" 
)

CNVkit.PoN.Res2 = Read.CNVkit.Result2(
    cnvkitRegionFile = sprintf("CNVkit_pon/%s_T.analysis.ready.cnr",      sampleID),
    cnvkitCallFile   = sprintf("CNVkit_pon/%s_T.analysis.ready.call.cns", sampleID), 
    answerFile       = "/storage/home/kangsm/analysis/cnv_test/wes_cnv_results_02_new/meta/seqc_hcc1395_cnvkit_trueset.txt",
    refVer="hg38" 
)

  draw.CNVplot2 = function( CNV_DATA, CHR="all", byCHR=FALSE, TITLE=NULL, MAX_Y=NA, baseZero = TRUE, refVer='hg19') 
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
      geom_segment(aes( x = index, y = cnv2, xend = index+1, yend = cnv2), color="#1a5387", linewidth=1) +
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

draw.CNVplot2( CNVkit.Paired.Res2, CHR="all", byCHR=FALSE, TITLE="CNVkit (PAIRED)", MAX_Y=8, baseZero = TRUE, refVer='hg38')

draw.CNVplot2( CNVkit.PoN.Res2, CHR="all", byCHR=FALSE, TITLE="CNVkit (PoN)", MAX_Y=8, baseZero = TRUE, refVer='hg38')


chrom='8'
cnvPlot1=draw.CNVplot2( CNVkit.Paired.Res2, CHR=chrom, byCHR=FALSE, TITLE="CNVkit (PAIRED)", MAX_Y=12, baseZero = TRUE, refVer='hg38')
cnvPlot2=draw.CNVplot( CNVkit.Paired.Res2, CHR=chrom, byCHR=FALSE, TITLE="CNVkit (PAIRED)", MAX_Y=12, baseZero = TRUE, refVer='hg38')
cnvPlot3=draw.CNVplot2( CNVkit.PoN.Res2, CHR=chrom, byCHR=FALSE, TITLE="CNVkit (PoN)", MAX_Y=12, baseZero = TRUE, refVer='hg38')
cnvPlot4=draw.CNVplot( ADTEx.Paired.Res,  CHR = chrom, byCHR = FALSE, TITLE="ADTEx",      MAX_Y=12,  baseZero = TRUE, refVer='hg38')
#cnvPlot5 = draw.CNVplot( EXCAVATOR2.Paired.Res, CHR=chrom, byCHR=FALSE, TITLE="EXCAVATOR2 (PoN)", MAX_Y=6, baseZero = TRUE, refVer='hg38')
#cnvPlot6 = draw.CNVplot( EXCAVATOR2.PoN.Res, CHR=chrom, byCHR=FALSE, TITLE="EXCAVATOR2 (PoN)", MAX_Y=6, baseZero = TRUE, refVer='hg38')
cowplot::plot_grid(cnvPlot1, cnvPlot2, cnvPlot3, cnvPlot4, ncol=2 )



