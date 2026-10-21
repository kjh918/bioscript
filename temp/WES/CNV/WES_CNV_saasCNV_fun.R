

vcf2txt2 <- function(vcf.file, normal.col=10, tumor.col=11, MQ.cutoff=30)
{
	inputFile <- file(vcf.file, "r")
    ## skip some definition lines
    while(length(lines <- readLines(inputFile, n = 1, warn = FALSE)) > 0) {
        ##message(lines)
        if (length(grep("\\#CHROM",lines))==1) break
    }
    header <- sub("\\#","",lines)
    header <- strsplit(header,"\t")[[1]]
    header[c(normal.col, tumor.col)] <- c("Normal","Tumor")

    ## read in data table
    vcf <- read.table(inputFile,header=FALSE,sep="\t",as.is=TRUE,comment.char="") ##na.strings="."
    names(vcf) <- header
    vcf <- vcf[,c(1:9, normal.col, tumor.col)]
    close(inputFile)

    ## filter1: chr1-22,X,Y and PASS
    chrs <- paste0("chr",c(1:22,"X","Y"))
    idx <- with(vcf, CHROM %in% chrs & FILTER=="PASS")
    vcf <- vcf[idx,]

    ## filter2: only have 1 ALT allele
    idx <- grep(",",vcf$ALT)
    if (length(idx)>=1) vcf <- vcf[-idx,]

    ## find REF and ALT read depth
    ## normal
    normalRD = apply(vcf[,c("FORMAT", "Normal")], 1, function(y) {
        vcfValuesN = unlist(strsplit(y[2], ":"))
        if( length(vcfValuesN) == 1 )
        {
            res = data.frame( Normal.GT="./.", Normal.REF.DP=0, Normal.ALT.DP=0 )
        }else{
            valueNames = unlist(strsplit(y[1], ":"))
            names(vcfValuesN) = valueNames
            val.GT = vcfValuesN["GT"]
            val.AD = unlist(strsplit(vcfValuesN["AD"], ","))
            val.AD.REF = val.AD[1]
            val.AD.ALT = val.AD[2]
            res = data.frame( Normal.GT=val.GT, Normal.REF.DP=as.numeric(val.AD.REF), Normal.ALT.DP=as.numeric(val.AD.ALT) )
        }
        return(res)
    })
    normalRD = do.call(rbind,normalRD)

    tumorRD = apply(vcf[,c("FORMAT", "Tumor")], 1, function(y) {
        vcfValuesN = unlist(strsplit(y[2], ":"))
        if( length(vcfValuesN) == 1 )
        {
            res = data.frame( Tumor.GT="./.", Tumor.REF.DP=0, Tumor.ALT.DP=0 )
        }else{
            valueNames = unlist(strsplit(y[1], ":"))
            names(vcfValuesN) = valueNames
            val.GT = vcfValuesN["GT"]
            val.AD = unlist(strsplit(vcfValuesN["AD"], ","))
            val.AD.REF = val.AD[1]
            val.AD.ALT = val.AD[2]
            res = data.frame( Tumor.GT=val.GT, Tumor.REF.DP=as.numeric(val.AD.REF), Tumor.ALT.DP=as.numeric(val.AD.ALT) )
        }
        return(res)
    })
    tumorRD = do.call(rbind,tumorRD)
    vcf = cbind(vcf, normalRD, tumorRD)

    ## reduce # columns
    vcf <- vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "INFO", 
                  "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]

    ## filter3: normal or tumor GT=="./."
    idx <- which(vcf$Normal.GT=="./." | vcf$Tumor.GT=="./.")
    if (length(idx)>=1) vcf <- vcf[-idx,]

    ## in INFO find MQ
    # info <- strsplit(vcf$INFO,";")
    # extract.MQ <- function(x) {
    #     idx.MQ <- grep("^MQ\\=",x)
    #     if (length(idx.MQ)==1) return(x[idx.MQ]) else return(NA)
    #     }
    # info1 <- sapply(info,FUN=extract.MQ)
    # vcf$MQ <- as.numeric(sub("^MQ\\=","",info1))

    ## filter4: MQ
    # vcf <- vcf[vcf$MQ>=MQ.cutoff,]

    ## reduce # columns
    # vcf <- vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "MQ", 
    #               "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]

    ## remove duplicate
    a <- vcf[,c("CHROM","POS")]
    idx <- which(duplicated(a))
    if (length(idx)>=1) vcf <- vcf[-idx,]
    
    return(vcf)
}

cnv.data <- function(vcf, min.chr.probe=100, verbose=FALSE) 
{
	## initialize
	idx <- which(vcf$Normal.GT=="0/1")
    vcf <- vcf[idx,]
    normal.total.coverage <- sum(vcf$Normal.REF.DP + vcf$Normal.ALT.DP, na.rm=TRUE)
    tumor.total.coverage  <- sum(vcf$Tumor.REF.DP + vcf$Tumor.ALT.DP, na.rm=TRUE)

    ## count #probes in each chr
    chr.probe.count <- table(vcf$CHROM)
    chrs <- names(chr.probe.count)[chr.probe.count >= min.chr.probe]
    chrs <- sub("^chr","",chrs)
    all.chrs <- c(1:22,"X","Y")
    chrs <- all.chrs[all.chrs %in% chrs] ## this is to make the order of chr as what we usually see 1-22, X, Y

    data <- NULL
    for (chr in chrs) {
	    ## data processing
	    vcf.chr <- vcf[vcf$CHROM==paste0("chr",chr),]
	    vcf.chr <- vcf.chr[order(vcf.chr$POS),]
	    data.chr <- cnv.data.chr(chr=vcf.chr$CHROM, position=vcf.chr$POS, 
	                             normal.ref.rd=vcf.chr$Normal.REF.DP, normal.alt.rd=vcf.chr$Normal.ALT.DP, 
	                             tumor.ref.rd=vcf.chr$Tumor.REF.DP, tumor.alt.rd=vcf.chr$Tumor.ALT.DP,
	                             normal.total.coverage=normal.total.coverage, tumor.total.coverage=tumor.total.coverage)
    
        ## collect information
        data <- rbind(data, data.chr)
    
        if (verbose) cat("finish data transformation of chr =",chr,"\n")
    } ## for chr

    return(data)
}

NGS.CNV <-
function (vcf, output.dir, sample.id,
          do.GC.adjust = FALSE,
          gc.file = system.file("extdata","GC_1kb_hg19.txt.gz",package="saasCNV"),
          min.chr.probe=100,
          min.snps=10,
          joint.segmentation.pvalue.cutoff=1e-4,
          max.chpts=30,
          do.merge=TRUE, use.null.data=TRUE, num.perm=1000, maxL=NULL, 
          merge.pvalue.cutoff=0.05,
          do.cnvcall.on.merge=TRUE, 
          cnvcall.pvalue.cutoff=0.05,
          do.plot=TRUE, cex=0.3, ref.num.probe=NULL,
          do.gene.anno=FALSE,
          gene.anno.file=NULL,
          seed=NULL,
          verbose=TRUE)
{	
	## create output dir
	if (!file.exists(output.dir)) dir.create(output.dir)
	
	## data process
	seq.data <- cnv.data(vcf=vcf, min.chr.probe=min.chr.probe, verbose=verbose) 
	
	if (do.GC.adjust) {
		gc <- read.delim(file = gc.file, as.is=TRUE)
        seq.data <- GC.adjust(data = seq.data, gc = gc, maxNumDataPoints = 10000)
	}
	
	data.dir <- paste0(output.dir,"/mid_data")
    if (!file.exists(data.dir)) dir.create(data.dir)
    save(list="seq.data", file=file.path(data.dir,"seq.data.RData"))
	
	## Step I: joint segmentation
	seq.segs <- joint.segmentation(data=seq.data, min.snps=min.snps, 
	                               global.pval.cutoff=joint.segmentation.pvalue.cutoff, 
	                               max.chpts=max.chpts, verbose=verbose)
    write.table(seq.segs, file=file.path(data.dir,"seq.segs.txt"), 
                quote=FALSE, row.names=FALSE, sep="\t")
    
    ## segmentation diagnosis plot
    if (do.plot) {
        plot.dir <- paste0(output.dir,"/mid_seg_plot")
        if (!file.exists(plot.dir)) dir.create(plot.dir)
        chrs <- sub("^chr","",unique(seq.segs$chr))
        
        for (chr in chrs) {
            png(filename=file.path(plot.dir, paste0("seg_chr",chr,".png")), width=240*5, height=240*4)
            diagnosis.seg.plot.chr(data=seq.data, segs=seq.segs, sample.id=sample.id, chr=chr, cex=cex)
            dev.off()
        }
    } ## do.plot
    
    ## other diagnosis plot
    if (do.plot) {
        diagnosis.plot.dir <- file.path(output.dir,"mid_diagnosis_plot")
        if (!file.exists(diagnosis.plot.dir)) dir.create(diagnosis.plot.dir)
        chrs <- sub("^chr","",unique(seq.segs$chr))
        
        png(filename=file.path(diagnosis.plot.dir,"diag_QQ_plot.png"), width=240*3, height=240*2)
        diagnosis.QQ.plot(data=seq.data, chrs=chrs)
        dev.off()    	
    } ##do.plot


    ## Step II (optional): merge segments not far apart in log2ratio and log2mBAF
    if (do.merge) {
    	seq.segs.merge <- merging.segments(data=seq.data, segs.stat=seq.segs, 
    	                                 use.null.data=use.null.data,
                                         N=num.perm, maxL=maxL, 
                                         merge.pvalue.cutoff=merge.pvalue.cutoff, 
                                         seed=seed,
                                         verbose=verbose)

        write.table(seq.segs.merge, file=file.path(data.dir,"seq.segs.merge.txt"), 
                    quote=FALSE, row.names=FALSE, sep="\t")
        
        if (do.plot) {
            plot.dir <- paste0(output.dir,"/mid_seg_merge_plot")
            if (!file.exists(plot.dir)) dir.create(plot.dir)
            chrs <- sub("^chr","",unique(seq.segs.merge$chr))
        
            for (chr in chrs) {
                png(filename=file.path(plot.dir, paste0("merge_seg_chr",chr,".png")), width=240*5, height=240*4)
                diagnosis.seg.plot.chr(data=seq.data, segs=seq.segs.merge, sample.id=sample.id, chr=chr, cex=cex)
                dev.off()
            }
        } ## do.plot
    } ## do.merge


    ## Step III: call CNV
    cat("CNV calling ...\n")
    if (do.merge==TRUE & do.cnvcall.on.merge==TRUE) {
        seq.cnv <- cnv.call(data=seq.data, sample.id=sample.id, segs.stat=seq.segs.merge, 
                            maxL=maxL, N=num.perm,
                            pvalue.cutoff=cnvcall.pvalue.cutoff,
                            seed=seed)
    } else {
    	seq.cnv <- cnv.call(data=seq.data, sample.id=sample.id, segs.stat=seq.segs, 
    	                    maxL=maxL, N=num.perm,
                            pvalue.cutoff=cnvcall.pvalue.cutoff,
                            seed=seed)

    }
    res.dir <- paste0(output.dir,"/mid_res")
    if (!file.exists(res.dir)) dir.create(res.dir)
    write.table(seq.cnv, file=file.path(res.dir,"seq.cnv.txt"), 
                quote=FALSE, row.names=FALSE, sep="\t")


    ## Step IV (option): gene annotation
    if (do.gene.anno) {
    	cat("gene annotation ...\n")
        gene <- read.delim(file=gene.anno.file, as.is=TRUE, comment.char="")
        seq.cnv.anno <- reannotate.CNV.res(res=seq.cnv, gene=gene, only.CNV=TRUE)
        write.table(seq.cnv.anno, file=file.path(res.dir,"seq.cnv.anno.txt"),
                    quote=FALSE, row.names=FALSE, sep="\t")
    }

    ## visualization
    ## genome-wide plot
    png(filename=file.path(res.dir,"cnv_gw_plot.png"), width=240*6, height=240*4)
    genome.wide.plot(data=seq.data, segs=seq.cnv, 
                     sample.id=sample.id, 
                     chrs=sub("^chr","",unique(seq.cnv$chr)), 
                     cex=cex)
    dev.off()

    ## diagnosis cluster plot
    pdf(file=file.path(res.dir,"cnv_cluster_plot.pdf"), width=8, height=8)
    diagnosis.cluster.plot(segs=seq.cnv, 
                           chrs=sub("^chr","",unique(seq.cnv$chr)), 
                           min.snps=min.snps, max.cex=3, 
                           ref.num.probe=ref.num.probe)
    dev.off()
}

