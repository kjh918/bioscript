

    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
    })

    option_list = list( 
        make_option(c("--WORK_DIR"), action="store", default=NULL, type="character", help="lichee analysis dir"),   
        make_option(c("--ANALYSIS_LIST"), action="store", default=NULL, type="character", help="analysis list file (sample info file)"),
        make_option(c("--VAF_TABLE"), action="store", default=NULL, type="character", help="vaf table file ( from vcf )")
    )
    #---------------------------------------------------------------------------  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #---------------------------------------------------------------------------
    work_dir           <- ARGS$WORK_DIR
    analysis_list_file <- ARGS$ANALYSIS_LIST
    vcf_vaf_table_file <- ARGS$VAF_TABLE


    # BaseDir            <- "/data/wgs"
    # BatchID            <- "WGS-26-01"
    # work_dir           <- "/data/wgs/WGS-26-01/analysis/phylowgs"
    # analysis_list_file <- "/data/wgs/WGS-26-01/meta/analysis.list"
    
    # vcf_vaf_table_file <- "WGS-26-01.tumors.merged.sage.vaf.tsv"
    
    sinfo <- read.table(sprintf("%s",analysis_list_file), header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    tumor_samples <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID

    for( tid in tumor_samples ){
        
        purple_dir = sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)
        
        cnv_data = read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",purple_dir,tid))
        qc_data   = read.delim(sprintf("%s/%s.purple.purity.tsv",purple_dir,tid))

        cnv_res <- cnv_data %>% mutate( minor_cn=round(minorAlleleCopyNumber, 0), major_cn=round(majorAlleleCopyNumber, 0)) %>% 
            mutate( copy_number = minor_cn+major_cn, cellular_prevalence=qc_data$purity ) %>% 
            dplyr::select(
                c("chromosome","start","end","copy_number","minor_cn","major_cn","cellular_prevalence")
            ) 
        write.table(cnv_res,
            sprintf("%s/%s.cnv.data.tsv", work_dir, tid),
            quote=F, col.names=T, row.names=F, sep="\t"
        )
    }
###











