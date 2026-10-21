

options(stringsAsFactors=FALSE) 
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("lubridate"))

#---| INPUT ARGS
    option_list = list(
        make_option(c("--sampleListFile"), action="store", default=NA, type="character", help="$SAMPLE_ID_FILE"),
        make_option(c("--dataDir"),        action="store", default=NA, type="character", help="$DATA_DIR"),
        make_option(c("--ngsLib"),         action="store", default=NA, type="character", help="$RGLB"),
    )
    ARGS = parse_args(OptionParser(option_list=option_list))
    SampleListFile = ARGS$sampleListFile
    DataDir        = ARGS$dataDir
    NgsLib         = ARGS$ngsLib
  #---/
    SampleList = read.delim(SampleListFile, header=FALSE)[,1]
    DatasetId  = tail(unlist(strsplit(DataDir, "\\/")),1)
    MultiqcResultDir = sprintf("%s/qc_files/qc_bam_data", DataDir)
    ResultSaveDir    = sprintf("%s/Results", DataDir)

#---| QC-REPORT DATA RE-FORMAT
    align.stat = read.delim(sprintf('%s/multiqc_picard_AlignmentSummaryMetrics.txt', MultiqcResultDir))
    dup.stat   = read.delim(sprintf('%s/multiqc_picard_dups.txt',                    MultiqcResultDir))
    hs.stat    = read.delim(sprintf('%s/multiqc_picard_HsMetrics.txt',               MultiqcResultDir))
    is.stat    = read.delim(sprintf('%s/multiqc_picard_insertSize.txt',              MultiqcResultDir))
    oxo.stat   = read.delim(sprintf('%s/multiqc_picard_OxoGMetrics.txt',             MultiqcResultDir))
#---|
    
bamQC = data.frame(Sample=align.stat$Sample, SAMPLE_ID=align.stat$Sample, BAM_STATUS=align.stat$Sample ) 
bamQC$BAM_STATUS = sapply(bamQC$BAM_STATUS, function(bq) paste(unlist(strsplit(bq, "\\."))[-1], collapse="." ))
bamQC$SAMPLE_ID  = sapply(bamQC$SAMPLE_ID,  function(si) unlist(strsplit(si, "\\."))[1])
bamQC$BAM_STATUS = gsub("sorted.bam","S",bamQC$BAM_STATUS)
bamQC$BAM_STATUS = gsub("sorted.dup.marked.bam","SD",bamQC$BAM_STATUS)
bamQC$BAM_STATUS = gsub("sorted.dup.marked.recal.bam","SDR",bamQC$BAM_STATUS)

bamQC = merge(bamQC, align.stat[, colnames(align.stat) %nin% c("SAMPLE","LIBRARY","READ_GROUP")], by='Sample', all.x=T)
bamQC = merge(bamQC, dup.stat[,colnames(dup.stat) %nin% 'LIBRARY'], by='Sample', all.x=T)
bamQC = merge(bamQC, hs.stat[, colnames(hs.stat) %nin% c("TOTAL_READS","PF_READS","PCT_PF_READS","SAMPLE","LIBRARY","READ_GROUP")], by='Sample', all.x=T)
bamQC = merge(bamQC, is.stat[,colnames(is.stat) %nin% c('Sample','SAMPLE','LIBRARY','READ_GROUP')] %>% rename( Sample = SAMPLE_NAME ), by='Sample', all.x=T)
bamQC = merge(bamQC, oxo.stat %>% filter( CONTEXT == 'CCG' ) %>%  select(c('Sample','OXIDATION_ERROR_RATE')) %>% 
  mutate( Sample = gsub("_CCG$", "", Sample)) %>% rename(CCG_OXIDATION_ERROR_RATE=OXIDATION_ERROR_RATE),  by='Sample', all.x=T)
bamQC$LIBRARY = NgsLib
bamQC$DATASET_ID = DatasetId
bamQC = bamQC[, colnames(bamQC) != "Sample"]    

#---| IMPORT INTO DB    
  #---/ Import Into DB
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='GCX_WES' )
    preExistClear = dbGetQuery(dbCon, sprintf("DELETE FROM wes_bam_stats WHERE DATASET_ID = '%s'", DatasetId) )
    dbWriteTable(dbCon, name='wes_bam_stats' , value = bamQC, row.names=FALSE, append=TRUE )
    dbDisconnect(dbCon)
  #---/
#---/  







