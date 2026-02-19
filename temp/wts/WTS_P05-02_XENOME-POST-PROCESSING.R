
### Packages 
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library(Hmisc)
        library(plyr)
        library(dplyr)
        library(setwidth)
        library(parallel)
        library(optparse)
    })
###
#---| OPTPARSE |---------------------------------------------------------------#
    option_list = list(
        make_option(c("--baseDir"), action="store", default=NULL, type="character", help="Base Analysis Dir"),
        make_option(c("--batchID"), action="store", default=NULL, type="character", help="Analysis batch-id")
        make_option(c("--seqID"), action="store", default=NULL, type="character", help="Sample-id"),
        make_option(c("--threads"), action="store", default=NULL, type="character", help="parallel processing cores")
    )
    ARGS = parse_args(OptionParser(option_list=option_list))
    BaseDir <- ARGS$baseDir
    BatchID <- ARGS$batchID
    SeqID   <- ARGS$seqID
    Threads <- as.numeric(ARGS$threads)
#------------------------------------------------------------------------------#
#---| Manual Params
    # BaseDir <- "/data/wts"
    # BatchID <- "WTS_25_02"
    # SeqID   <- "WTS_25_02_03"
    # Threads <- 10
#------------------------------------------------------------------------------#
#---| Index |
    FastqGroupIndex = list(
        input     = ".trimmed_R1.fastq",
        human     = "_hg38_1.fastq",
        mouse     = "_mm10_1.fastq",
        both      = "_both_1.fastq",
        ambiguous = "_ambiguous_1.fastq",
        neither   = "_neither_1.fastq"
    )

    ReadCountList <- mclapply(FastqGroupIndex, function(SuffixTag){
        RdCounts <- system(sprintf("cat %s/%s/%s/fastq/%s%s | wc -l", BaseDir, BatchID, SeqID, SeqID, SuffixTag), intern=T)
        RD_CNT <- data.frame(ReadCounts=RdCounts)
        return(RD_CNT)
    }, mc.cores=6 )

    ReadCount <- ldply(ReadCountList, .id="FastqGroup")
    InputRD <- ReadCount %>% filter( FastqGroup == "input" ) %>% .$ReadCounts %>% as.numeric()
    ReadCount <- ReadCount %>% mutate( GroupRatio = round(as.numeric(ReadCounts)/InputRD, 3))

    write.table( ReadCount, 
        sprintf("%s/%s/%s/qcfiles/%s.Fastq.Xenome.Split.Stats.tsv",BaseDir, BatchID, SeqID, SeqID),
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    write.table( ReadCount, 
        sprintf("%s/%s/%s/fastq/%s.Fastq.Xenome.Split.Stats.tsv",BaseDir, BatchID, SeqID, SeqID),
        quote=F, col.names=T, row.names=F, sep="\t"
    )

    FastqDir <- sprintf("%s/%s/%s/fastq", BaseDir, BatchID, SeqID)
    system(sprintf("bgzip -@ %s -c %s/%s_hg38_1.fastq > %s/%s.human.trimmed_R1.fastq.gz", Threads, FastqDir, SeqID, FastqDir, SeqID))
    system(sprintf("bgzip -@ %s -c %s/%s_hg38_2.fastq > %s/%s.human.trimmed_R2.fastq.gz", Threads, FastqDir, SeqID, FastqDir, SeqID))
    system(sprintf("bgzip -@ %s -c %s/%s_mm10_1.fastq > %s/%s.mouse.trimmed_R1.fastq.gz", Threads, FastqDir, SeqID, FastqDir, SeqID))
    system(sprintf("bgzip -@ %s -c %s/%s_mm10_2.fastq > %s/%s.mouse.trimmed_R2.fastq.gz", Threads, FastqDir, SeqID, FastqDir, SeqID))
    system(sprintf("rm -rf %s/*.fastq", FastqDir))












