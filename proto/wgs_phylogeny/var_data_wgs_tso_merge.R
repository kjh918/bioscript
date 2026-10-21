
###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        library("GenomicRanges")
        library("data.table")
    })
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
###
### 
    BaseDir="/data/wgs"
    BatchID="WGS-26-01"
    AnalysisList="/data/wgs/WGS-26-01/meta/analysis.list"
    
    VarDataDir <- sprintf("%s/%s/analysis/variant-sites-allele-count", BaseDir, BatchID)

    sinfo <- read.table(AnalysisList, header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    sinfo$SampleID <- make.names(sinfo$SeqID)
    tumor_ids      <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
    tumor_samples  <- sinfo %>% filter( Group != "normal" ) %>% .$SampleID
###
### var-data
    var_data_wgs <- read.delim(sprintf("%s/%s.tumors.variants.vaf.tsv", VarDataDir, BatchID)) %>% 
        mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) 

    var_data_tso <- read.delim(sprintf("%s/%s.tumors.tso.variants.vaf.tsv", VarDataDir, BatchID)) %>% 
        mutate( vkey=paste(CHROM,POS,REF,ALT, sep="_") ) 


    var_data_merged <- var_data_wgs

#    for( tid in tumor_ids ){
    for( j in 1:nrow(var_data_tso) ){
        wgs_matched <- which(var_data_merged$vkey == var_data_tso[j, "vkey"])
        if( length(wgs_matched) == 0 ){
            var_data_merged <- rbind(var_data_merged, var_data_tso[j, ])
        }else{
            if( var_data_tso[ j, make.names(tumor_ids[1]) ] != "." ){
                var_data_merged[wgs_matched, make.names(tumor_ids[1])] <- var_data_tso[ j, make.names(tumor_ids[1]) ]
            } 
            if( var_data_tso[ j, make.names(tumor_ids[2]) ] != "." ){
                var_data_merged[wgs_matched, make.names(tumor_ids[2])] <- var_data_tso[ j, make.names(tumor_ids[1]) ]
            }
        }
    }

    write.table(var_data_merged, sprintf("%s/%s.wgs.tso.merged.variants.vaf.tsv", VarDataDir, BatchID), quote=F, col.names=T, row.names=F, sep="\t")


