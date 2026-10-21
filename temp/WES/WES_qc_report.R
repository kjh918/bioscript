

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
        make_option(c("--dataDir"),        action="store", default=NA, type="character", help="$DATA_DIR")
    )
    ARGS = parse_args(OptionParser(option_list=option_list))
    SampleListFile = ARGS$sampleListFile
    DataDir        = ARGS$dataDir
  #---/
    SampleList = read.delim(SampleListFile, header=FALSE)[,1]
    DatasetId  = tail(unlist(strsplit(DataDir, "\\/")),1)
    MultiqcResultDir = sprintf("%s/qc_files/qc_report_data", DataDir)
    ResultSaveDir    = sprintf("%s/Results", DataDir)

#---| QC-REPORT DATA RE-FORMAT
    general.stat = read.delim(sprintf('%s/multiqc_general_stats.txt',                  MultiqcResultDir), header=T)
    align.stat   = read.delim(sprintf('%s/multiqc_picard_AlignmentSummaryMetrics.txt', MultiqcResultDir))
    dup.stat     = read.delim(sprintf('%s/multiqc_picard_dups.txt',                    MultiqcResultDir))
    hs.stat      = read.delim(sprintf('%s/multiqc_picard_HsMetrics.txt',               MultiqcResultDir))
    is.stat      = read.delim(sprintf('%s/multiqc_picard_insertSize.txt',              MultiqcResultDir))
  #---/ General Stats
    GS = general.stat %>% filter( Sample %in% SampleList) %>% 
        select(c('Sample','fastp_mqc.generalstats.fastp.after_filtering_q30_rate','fastp_mqc.generalstats.fastp.after_filtering_gc_content')) %>% 
        rename( Q30_Rate=fastp_mqc.generalstats.fastp.after_filtering_q30_rate, GC_Rate=fastp_mqc.generalstats.fastp.after_filtering_gc_content) 
  # Picard AlignmentSummary Metric
    AS = align.stat %>% filter( Sample %in% SampleList) %>% 
        select(c('Sample','PF_READS','PF_READS_ALIGNED','PCT_PF_READS_ALIGNED')) %>% 
        rename( TOTAL_READS=PF_READS, ALIGNED_READS=PF_READS_ALIGNED, ALIGNED_READS_Rate=PCT_PF_READS_ALIGNED )
  #---/ Picard DuplicateMark Metric
    DS = dup.stat %>% filter( Sample %in% SampleList) %>% select(c('Sample','PERCENT_DUPLICATION')) %>% rename( DUPLICATION_Rate=PERCENT_DUPLICATION )
  #---/ Picard HS Metric
    HS = hs.stat %>% filter( Sample %in% SampleList) %>%
        select(c('Sample','PF_BASES','PCT_USABLE_BASES_ON_TARGET','MEAN_TARGET_COVERAGE',
                 'PCT_TARGET_BASES_1X','PCT_TARGET_BASES_10X','PCT_TARGET_BASES_30X','PCT_TARGET_BASES_50X','PCT_TARGET_BASES_100X','PCT_TARGET_BASES_250X')) %>% 
        rename( TOTAL_BASES=PF_BASES, ON_TARGET_Rate=PCT_USABLE_BASES_ON_TARGET, MEAN_TARGET_DEPTH=MEAN_TARGET_COVERAGE, 
                TARGET_COVERAGE_1X=PCT_TARGET_BASES_1X, TARGET_COVERAGE_10X=PCT_TARGET_BASES_10X, TARGET_COVERAGE_30X=PCT_TARGET_BASES_30X, 
                TARGET_COVERAGE_50X=PCT_TARGET_BASES_50X, TARGET_COVERAGE_100X=PCT_TARGET_BASES_100X, TARGET_COVERAGE_250X=PCT_TARGET_BASES_250X )
  #---/ Picard InsertSize Metric
    IS = is.stat %>% filter( SAMPLE_NAME %in% SampleList) %>% 
        select(c('SAMPLE_NAME','MEDIAN_INSERT_SIZE')) %>% rename( Sample=SAMPLE_NAME, FRAGMENT_SIZE=MEDIAN_INSERT_SIZE )
  #---/ MERGE AND REORDER COLUMNS & DIGIT ADJUSTMENT
    REPORT = merge( AS, GS, by=1)
    REPORT = merge(REPORT, DS, by=1)
    REPORT = merge(REPORT, HS, by=1)
    REPORT = merge(REPORT, IS, by=1)
    REPORT = REPORT[match(SampleList, REPORT$Sample),
        c('Sample','TOTAL_BASES','TOTAL_READS','ALIGNED_READS','ALIGNED_READS_Rate','Q30_Rate','GC_Rate','DUPLICATION_Rate','FRAGMENT_SIZE','MEAN_TARGET_DEPTH',
          'ON_TARGET_Rate','TARGET_COVERAGE_1X','TARGET_COVERAGE_10X','TARGET_COVERAGE_30X','TARGET_COVERAGE_50X','TARGET_COVERAGE_100X','TARGET_COVERAGE_250X')]
    for( k in c(5,6,7,8,11,12,13,14,15,16,17) ){ REPORT[,k] = round( REPORT[,k]*100, digit=2 ) }
    REPORT$MEAN_TARGET_DEPTH = round(REPORT$MEAN_TARGET_DEPTH, digit=2)
#---/    
    
#---| SAVE AS XLSX FILE & IMPORT INTO DB    
  #---/ Save File   
    openxlsx::write.xlsx(REPORT, sprintf('%s/%s_WES_QC_Report.xlsx', ResultSaveDir, DatasetId ), rowNames=F, overwrite=T)
  #---/ Import Into DB
    DbReport = REPORT %>% mutate( DatasetID = DatasetId, ReportDate = today() )
    dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.35', user='gcx', port=3306, password='gencurix!!', db='GCX_WES' )
    
    
    dbWriteTable(dbCon, name='wes_qc_report' , value = DbReport, row.names=FALSE, append=TRUE )
    dbDisconnect(dbCon)
  #---/
#---/  







