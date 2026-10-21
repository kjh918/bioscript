
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        suppressWarnings(library("dndscv"))
    })

    BaseDir="/data/wgs"
    BatchID="WGS-26-01"
    AnalysisList="/data/wgs/WGS-26-01/meta/analysis.list"
    WorkDir="/data/wgs/WGS-26-01/analysis/dndscv"

    sinfo <- read.table(sprintf("%s",AnalysisList), header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    tumor_samples <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID

    var_data = read.delim(sprintf("%s/%s/analysis/%s.tumors.variants.vaf.tsv", BaseDir, BatchID, BatchID))

    VAR_DATA = data.frame()
    for( tid in tumor_samples ){
        cn = which(colnames(var_data) == make.names(tid))

        vd <- var_data[which(var_data[,cn] != "."), c("CHROM","POS","REF","ALT")] %>% mutate(sampleID=tid) %>% dplyr::rename(
            chr=CHROM, pos=POS, ref=REF, mut=ALT
        ) %>% filter( nchar(ref) == 1, nchar(mut) == 1)
        VAR_DATA = rbind(VAR_DATA, vd)
    } 
    VAR_DATA <- VAR_DATA[,c("sampleID","chr","pos","ref","mut")]


    VAR_DATA2 = VAR_DATA %>% mutate( chr = paste0("chr", chr))

    dndout <- dndscv( 
        mutations=VAR_DATA,
        max_muts_per_gene_per_sample = 10,
        mingenecovs = 50
    )

    saveRDS(dndout, file=sprintf("%s/%s.dndscv.output.RDS", WorkDir, BatchID))

