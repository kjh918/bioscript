


### cn.MOPS ###------------------------------------------------------------------------------------#
library(cn.mops)
library(plyr)
library(dplyr)
library(parallel)

setwd("/storage/home/kangsm/analysis/cnv_test/cn.mops")

NORMAL_BAM_DIR = "/storage/home/kangsm/analysis/cnv_test/normal_bams"
TUMOR_BAM_DIR  = "/storage/home/kangsm/analysis/cnv_test/tumor_bams"

twistBed   = "/storage/home/kangsm/myDB/bed/ADTEx_hg19/ADTEx_twist.exome.2.0.target.bed"
twistSeg   = read.table(twistBed, sep="\t", as.is=TRUE)
twistRange = GRanges(twistSeg[,1], IRanges(twistSeg[,2], twistSeg[,3]))

agilentBed   = "/storage/home/kangsm/myDB/bed/ADTEx_hg19/ADTEx_SureSelect.HumanAllExon.v6.extended.target.bed"
agilentSeg   = read.table(agilentBed, sep="\t", as.is=TRUE)
agilentRange = GRanges(agilentSeg[,1], IRanges(agilentSeg[,2], agilentSeg[,3]))

# 1) TUMOR ONLY SMAPLES : exomecn.mops
    tumorOnlySamples = c("PON01_01_T","PON01_02_T","PON01_03_T","PON01_04_T","PON01_05_T")

    tumorSamplesBamFiles = c()
    for( i in 1:length(tumorOnlySamples) )
    {
        bamFile = list.files( TUMOR_BAM_DIR, pattern=tumorOnlySamples[i] )
        bamFile = bamFile[ grep("bam$", bamFile) ]
        tumorSamplesBamFiles = c(tumorSamplesBamFiles, bamFile)
    }
    tumorSamplesBamFiles = paste(TUMOR_BAM_DIR, tumorSamplesBamFiles, sep="/")

    gRangeTumorOnlyGroupBams = getSegmentReadCountsFromBAM(
        BAMFiles = tumorSamplesBamFiles,
        GR = agilentRange,
        sampleNames = tumorOnlySamples,
        parallel = length(tumorOnlySamples)
    )
    TumorOnlySampleGroupExomeCNMOPS = exomecn.mops(gRangeTumorOnlyGroupBams)
    PON01_TumorOnlyBatch_exomeCNMOPS = calcIntegerCopyNumbers(TumorOnlySampleGroupExomeCNMOPS)

    save(PON01_TumorOnlyBatch_exomeCNMOPS, file="PON01_TumorOnlyBatch_exomeCNMOPS.Rdata")

## PoN Samples
    PoNsamples = c("PON01_01_N","PON01_02_N","PON01_03_N","PON01_04_N","PON01_05_N")
    PoNsamplesBamFiles = c()
    for( i in 1:length(PoNsamples) )
    {
        bamFile = list.files( NORMAL_BAM_DIR, pattern=PoNsamples[i] )
        bamFile = bamFile[ grep("bam$", bamFile) ]
        PoNsamplesBamFiles = c(PoNsamplesBamFiles, bamFile)
    }
    PoNsamplesBamFiles = paste(NORMAL_BAM_DIR, PoNsamplesBamFiles, sep="/")

    gRangePoNbams = getSegmentReadCountsFromBAM(
        BAMFiles    = PoNsamplesBamFiles,
        GR          = agilentRange,
        sampleNames = PoNsamples,
        parallel    = length(PoNsamples)
    )
    PoN_segmentMatrix = as.data.frame(gRangePoNbams)
    
    save(PoN_segmentMatrix, file="PanelofNormal_segmentMatrix.Rdata")

    #==>>  apply(PoNsegmentMatrix[,PoNsamples], 1, median)

# 2) 1 Tumor Samples + PoN Samples

    idx = sapply(tumorOnlySamples, function(y) list(y))
    
    singleTumor.PoN.referenceCNMOPS.list = mclapply( idx, function(ts) {
        gRangeSingleTumorBam = getSegmentReadCountsFromBAM(
            BAMFiles = sprintf("%s/%s.analysis.ready.bam", TUMOR_BAM_DIR, ts),
            GR = agilentRange,
            sampleNames = ts
        )
        gRangeSingleTumorBam$PON_MEDIAN = apply(PoN_segmentMatrix[,PoNsamples], 1, median)

        singleTumorRefCNMOPS = referencecn.mops( 
            cases    = gRangeSingleTumorBam[,1], 
            controls = gRangeSingleTumorBam[,2],
            classes  = c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN16","CN32","CN64","CN128"),
            I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64), segAlgorithm = "DNAcopy"
        )
        singleTumorRefCNMOPS = calcIntegerCopyNumbers(singleTumorRefCNMOPS)
        return(singleTumorRefCNMOPS)
    }, mc.cores = length(idx) )

    save(singleTumor.PoN.referenceCNMOPS.list, file="singleTumor.PoN.referenceCNMOPS.list.Rdata")

# 3) 1 Tumor + 5 Normal exomecn.mops
    normalBam = list.files( NORMAL_BAM_DIR )
    normalBam = normalSamples = normalBam[ grep("bam$", normalBam) ]
    normalBam = paste(NORMAL_BAM_DIR, normalBam, sep="/")
    normalSamples = gsub("\\.analysis.ready.bam$","", normalSamples)

    singleTumor.PoN.exomeCNMOPS.list = list()
    for( ts in tumorOnlySamples )
    {
        tumorBam = sprintf("%s/%s.analysis.ready.bam", TUMOR_BAM_DIR, ts)
        BamList  = c(tumorBam, normalBam)
        gRangeSingleTumorPoNbams = getSegmentReadCountsFromBAM(
            BAMFiles = BamList,
            GR = agilentRange,
            sampleNames = c(ts,normalSamples),
            parallel = length(c(ts,normalSamples))
        )
        singleTumorPoNexomeCNMOPS = exomecn.mops(gRangeSingleTumorPoNbams)
        singleTumorPoNexomeCNMOPS = calcIntegerCopyNumbers(singleTumorPoNexomeCNMOPS)
        
        singleTumor.PoN.exomeCNMOPS.list = c(singleTumor.PoN.exomeCNMOPS.list, list(singleTumorPoNexomeCNMOPS))
        names(singleTumor.PoN.exomeCNMOPS.list)[length(singleTumor.PoN.exomeCNMOPS.list)] = ts
        print(ts)
    }
    save(singleTumor.PoN.exomeCNMOPS.list, file="singleTumor.PoN.exomeCNMOPS.list.Rdata")
    
# 4) N-T paired exomecn.mops
    ID.list = c("PON01_01","PON01_02","PON01_03","PON01_04","PON01_05")
    idx2 = sapply(ID.list, function(y) list(y))

    NTpairedBAM.list = list()
    for( sn in ID.list )
    {
        bam_N = list.files( NORMAL_BAM_DIR, pattern=sprintf("%s_N", sn) )
        bam_T = list.files( TUMOR_BAM_DIR,  pattern=sprintf("%s_T", sn) )
        bam_N = bam_N[ grep("bam$", bam_N) ]
        bam_T = bam_T[ grep("bam$", bam_T) ]
        pairedSampleBamFiles = c(
            paste(NORMAL_BAM_DIR, bam_N, sep="/"), 
            paste(TUMOR_BAM_DIR,  bam_T, sep="/")
        )
        gRangeNTpairedBams = getSegmentReadCountsFromBAM(
            BAMFiles = pairedSampleBamFiles,
            GR = agilentRange,
            sampleNames = c(sprintf("%s_N", sn),sprintf("%s_T", sn)),
            parallel = 2
        )
        NTpairedBAM.list = c(NTpairedBAM.list, list(gRangeNTpairedBams))
        names(NTpairedBAM.list)[length(NTpairedBAM.list)] = sn
    }

    PON01.NTpaired.exomeCNMOPS.list = mclapply( NTpairedBAM.list, function(nt) {
        excnmops = exomecn.mops(nt)
        if( length(excnmops@cnvs@seqnames) == 0 )
        {
            return("No CNVs")
        }else{
            excnmops = calcIntegerCopyNumbers(excnmops)
            return(excnmops)
        }
    }, mc.cores = length(ID.list))

    PON01.NTpaired.referenceCNMOPS.list = mclapply( NTpairedBAM.list, function(nt) {
        NTrefcnmops = referencecn.mops( 
            cases    = nt[,2], 
            controls = nt[,1],
            classes  = c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN16","CN32","CN64","CN128"),
            I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64), segAlgorithm = "DNAcopy"
        )
        NTrefcnmops = calcIntegerCopyNumbers(NTrefcnmops)
        return(NTrefcnmops)
    }, mc.cores = length(ID.list) )

    save(PON01.NTpaired.exomeCNMOPS.list, file="PON01.NTpaired.exomeCNMOPS.list.Rdata")
    save(PON01.NTpaired.referenceCNMOPS.list, file="PON01.NTpaired.referenceCNMOPS.list.Rdata")

#  5) N-T Pair Multiple Samples
    tumorSamples  = c("PON01_01_T","PON01_02_T","PON01_03_T","PON01_04_T","PON01_05_T")
    normalSamples = c("PON01_01_N","PON01_02_N","PON01_03_N","PON01_04_N","PON01_05_N")
    groupNormalTumorBamFiles = c()
    for( i in 1:length(normalSamples) )
    {
        bamFile = list.files( NORMAL_BAM_DIR, pattern=normalSamples[i] )
        bamFile = bamFile[ grep("bam$", bamFile) ]
        bamFile = paste(NORMAL_BAM_DIR, bamFile, sep="/")
        groupNormalTumorBamFiles = c(groupNormalTumorBamFiles, bamFile)
    }
    for( i in 1:length(tumorSamples) )
    {
        bamFile = list.files( TUMOR_BAM_DIR, pattern=tumorSamples[i] )
        bamFile = bamFile[ grep("bam$", bamFile) ]
        bamFile = paste(TUMOR_BAM_DIR, bamFile, sep="/")
        groupNormalTumorBamFiles = c(groupNormalTumorBamFiles, bamFile)
    }
    gRangeNormTumorGroupBams = getSegmentReadCountsFromBAM(
        BAMFiles = groupNormalTumorBamFiles,
        GR = agilentRange,
        sampleNames = c(normalSamples,tumorSamples),
        parallel = length(c(normalSamples,tumorSamples))
    )
    PON01.Pooled.All.exomeCNMOPS = exomecn.mops(gRangeNormTumorGroupBams)
    PON01.Pooled.All.exomeCNMOPS = calcIntegerCopyNumbers(PON01.Pooled.All.exomeCNMOPS)

    save(PON01.Pooled.All.exomeCNMOPS, file="PON01.Pooled.All.exomeCNMOPS.Rdata")

#  6) NA12878 batch
    NA12878_BAMs = c(
        "/storage/home/kangsm/analysis/cnv_test/tumor_bams/NA12878_01.analysis.ready.bam",
        "/storage/home/kangsm/analysis/cnv_test/tumor_bams/NA12878_02.analysis.ready.bam"
    )

    gRangeNA12878.Bams = getSegmentReadCountsFromBAM(
        BAMFiles = NA12878_BAMs,
        GR = twistRange,
        sampleNames = c("NA12878_01","NA12878_02"),
        parallel = 2
    )
    NA12878.batch.exomeCNMOPS = exomecn.mops(gRangeNA12878.Bams)
    NA12878.batch.exomeCNMOPS = calcIntegerCopyNumbers(NA12878.batch.exomeCNMOPS)

    save(NA12878.batch.exomeCNMOPS, file="NA12878.batch.exomeCNMOPS.Rdata")

#

InputCNV = singleTumor.PoN.referenceCNMOPS.list[[1]]    
    
table(as.data.frame(cnvs(InputCNV))[,c("sampleName","CN")])

segplot(TumorOnlySampleGroupExomeCNMOPS, seqnames=6)

SEG = as.data.frame(segmentation(TumorOnlySampleGroupExomeCNMOPS))
CNVR = as.data.frame(cnvr(normTumorGroupExomeCNMOPS))

###

#--------------------------------------------------------------------------------------------------#

as.data.frame(cnvs(NA12878.batch.exomeCNMOPS))



table(CNV.tumorGroup[,c("sampleName","CN")])
table(CNV.NT.cnmops[,c("sampleName","CN")])
table(CNV.NT.ref[,c("sampleName","CN")])
table(CNV.NT.ref.mean[,c("sampleName","CN")])
table(CNV.allSample[,c("sampleName","CN")])


RES_DIR="/storage/home/kangsm/analysis/cnv_test/saascnv"

as.data.frame(cnvs(NA12878.batch.exomeCNMOPS)) %>% head
as.data.frame(cnvr(NA12878.batch.exomeCNMOPS)) %>% head
as.data.frame(segmentation(NA12878.batch.exomeCNMOPS)) %>% head

as.data.frame(integerCopyNumber(NA12878.batch.exomeCNMOPS)) %>% head


#--------------------------------------------------------------------------------------------------#
library(GenVisR)
library(stringr)

CN =  c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64)
names(CN) = c("CN0","CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN16","CN32","CN64","CN128")



a <- function(x){
    # read data and set column names
    data <- read.delim(x, header=FALSE)
    colnames(data) <- c("chromosome", "start", "end", "probes", "segmean")

    # get the sample name from the file path
    sampleName <- str_extract(x, "H_OM.+cc2")
    sampleName <- gsub(".cc2", "", sampleName)
    data$sample <- sampleName

    # return the data
    return(data)
}
cnData <- lapply(files, a)
cnData <- do.call("rbind", cnData)

cnFreq(cnData, genome="hg19")
cnFreq(cnData, genome="hg19", CN_low_cutoff = 0, CN_high_cutoff = .1)
library(ggplot2)
layer1 <- geom_vline(xintercept=c(39709170))
cnFreq(cnData, genome="hg19", plotChr="chr17", plotLayer=layer1)

InputCNV = singleTumor.PoN.referenceCNMOPS.list[[1]]    
as.data.frame(segmentation(InputCNV)) %>% 
    select(c("seqnames","start","end", "width","mean", "sampleName","CN")) %>%
    dplyr::rename(chromosome=1, pobe=width, sample=sampleName) %>%
    mutate( chromosome = gsub("^","chr", chromosome)) -> cnData
cnData$segmean = CN[ cnData$CN ]
cnData$segmean = cnData$segmean-1


cnFreq(A, genome="hg19", CN_low_cutoff = 0.5, CN_high_cutoff = 1.5)


as.data.frame(segmentation(InputCNV)) %>% select(c(1,2,3,6)) %>% 
    dplyr::rename(chromosome=1, start=2, end=3, CN=4) %>%
    mutate( chromosome = gsub("^","chr", chromosome), sample="PON01_01_T" ) -> cnData
cnData$segmean = CN[ cnData$CN ]

cnFreq(cnData, genome="hg19", plotChr="chr9", CN_low_cutoff = 0.7, CN_high_cutoff = 1.3) 


### RESULT 정리

# 1) ADTEx 
# METHOD : N-T Paired 
RES_BASE_DIR="/storage/home/kangsm/analysis/cnv_test/adtex"

seqID="PON01_01"
cnvCall = read.delim(sprintf("%s/%s/%s.cnv", RES_BASE_DIR, smeaneqID, seqID))
cnvRes  = read.delim(sprintf("%s/%s/cnv.result", RES_BASE_DIR, seqID))
cnvRes2 = cnvRes %>% select(c("chr","cnv","CNV_start","CNV_end","seg_mean")) %>% unique()


