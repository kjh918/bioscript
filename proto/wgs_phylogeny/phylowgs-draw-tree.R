

###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        library("jsonlite")
        library("data.table")
    })
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
###
### 
    BaseDir="/data/wgs"
    BatchID="WGS-26-01"
    AnalysisList="/data/wgs/WGS-26-01/meta/analysis.list"
    
    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)
    WorkDir    <- sprintf("%s/%s/analysis/phylowgs", BaseDir, BatchID)
    system(sprintf("mkdir -p %s", WorkDir))

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###
###
    trees <- fromJSON(gzfile(sprintf("%s/%s.summ.json.gz", WorkDir,BatchID )))
    muts  <- fromJSON(gzfile(sprintf("%s/%s.muts.json.gz", WorkDir,BatchID )))

    
    tree_list <- trees$trees
    best_tree <- tree_list[[ which.max(lapply( tree_list, function(x) x$llh )) ]]

