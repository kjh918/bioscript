
# MODIFICATION SCRIPT OF MUTECT2 FILTERED VCF TO PureCN Panel-of-Normal INPUT VCF FORMAT

#---| PACKAGES |-----------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("dplyr"))
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--SEQ_ID"),  action="store", default=NA, type="character", help="$SEQ_ID"),
        make_option(c("--VCF_DIR"), action="store", default=NA, type="character", help="$VCF_DIR")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    SEQ_ID  <- ARGS$SEQ_ID
    VCF_DIR <- ARGS$VCF_DIR
#--------------------------------------------------------------------------------------------------#

#---| FILES AND FOLDERS |--------------------------------------------------------------------------#
    VCF_FILE     <- sprintf("%s.mutect2.filtered.vcf", SEQ_ID)
    PON_VCF_FILE <- sprintf("%s.purecn.pon.vcf",       SEQ_ID)
    #--------------------------------------------------------------------------#
    if( ! file.exists(sprintf("%s/%s", VCF_DIR, VCF_FILE)) )
    { message("\n ** NO INPUT VCF FILE. STOPPED.") ; quit( save="no", status=0 ) }
#--------------------------------------------------------------------------------------------------#

#---| READ VCF AND REFORMAT |----------------------------------------------------------------------#
    read_vcf <- read.vcfR(sprintf("%s/%s", VCF_DIR, VCF_FILE))
    vcf      <- cbind(as.data.frame(getFIX(read_vcf)), extract.gt( read_vcf, element="AD"))
    vcf_pf   <- vcf[which(vcf$FILTER=="PASS"),]
    #--------------------------------------------------------------------------#
    vcf_pf$INFO   <- "."
    vcf_pf$ID     <- "."
    vcf_pf$QUAL   <- "."
    vcf_pf$FORMAT <- "AD"
    #--------------------------------------------------------------------------#
    vcf_out <- vcf_pf %>% select(c(1,2,3,4,5,6,7,9,10,8)) %>% dplyr::rename("#CHROM"=CHROM) 
#--------------------------------------------------------------------------------------------------#

#---| WRITE NEW VCF |------------------------------------------------------------------------------#
    system(sprintf("grep ^##fileformat %s/%s > %s/%s",    VCF_DIR, VCF_FILE, VCF_DIR, PON_VCF_FILE))
    system(sprintf("grep ^##FORMAT %s/%s >> %s/%s",       VCF_DIR, VCF_FILE, VCF_DIR, PON_VCF_FILE))
    system(sprintf("grep ^##INFO %s/%s >>%s/%s",          VCF_DIR, VCF_FILE, VCF_DIR, PON_VCF_FILE))
    system(sprintf("grep ^##contig %s/%s >> %s/%s",       VCF_DIR, VCF_FILE, VCF_DIR, PON_VCF_FILE))
    system(sprintf("grep ^##tumor_sample %s/%s >> %s/%s", VCF_DIR, VCF_FILE, VCF_DIR, PON_VCF_FILE))
    #--------------------------------------------------------------------------#
    write.table(vcf_out, sprintf("%s/%s", VCF_DIR, PON_VCF_FILE), quote=F, col.names=T, row.names=F, sep="\t", append=T)
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
    message("PureCN Panel-of-Normal Input VCF Generated.")
#--------------------------------------------------------------------------------------------------#

