

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
    lichee_dir         <- ARGS$WORK_DIR
    analysis_list_file <- ARGS$ANALYSIS_LIST
    vcf_vaf_table_file <- ARGS$VAF_TABLE


    # lichee_dir         <- "/data/wgs/WGS-26-01/analysis/lichee"
    # analysis_list_file <- "analysis.list"
    # vcf_vaf_table_file <- "WGS-26-01.tumors.merged.sage.vaf.tsv"
    
    sinfo <- read.table(sprintf("%s/%s", lichee_dir, analysis_list_file), header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    tumor_samples <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID

    vcf_data <- read.delim(sprintf("%s/%s", lichee_dir, vcf_vaf_table_file))
    for( tid in tumor_samples ){

        id <- make.names(tid)
        vcf_data[, id] <- sapply(vcf_data[, id], function(ad) {
            if( ad == "." ){ 
                VAF <- 0
            }else{
                rd  <- unlist(strsplit(ad, ","))
                VAF <- round(as.numeric(rd[2]) / ( as.numeric(rd[1]) + as.numeric(rd[2]) ), 6)
            }
            return(VAF)
        })
    }
    
    for( ts in tumor_samples ){
        #ts = tumor_samples[1]
        colnames(vcf_data)[which(colnames(vcf_data) == make.names(ts))] <- sinfo[which(sinfo$SeqID == ts), "Group"]
    }

    output_colnames <- c("#chrom","pos","DESC","normal", c(sinfo %>% filter( Group !="normal" ) %>% arrange(Order) %>% .$Group))

    output_data <- vcf_data %>% mutate(DESC=paste0(CHROM,"_",POS), normal=0) %>% 
        dplyr::rename( '#chrom'=CHROM, 'pos'=POS ) %>% 
        dplyr::select(all_of(output_colnames)) 

    write.table(output_data, "lichee_analysis_input.txt", quote=F, col.names=T, row.names=F, sep="\t")


    










