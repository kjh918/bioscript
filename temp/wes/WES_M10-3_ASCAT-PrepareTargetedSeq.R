

#---| PACKAGES |------------------------------------------------------------------------------------
    options(stringsAsFactors=FALSE)
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("ASCAT"))
#---------------------------------------------------------------------------------------------------
#---| ARGUMENTS |-----------------------------------------------------------------------------------
    option_list = list(
        make_option(c("--BASE_DIR"        ), action="store", default=NULL,  type="character", help="ngs data base folder. REQUIRED."),
        make_option(c("--NGS_LIB"         ), action="store", default=NULL,  type="character", help="NGS library. REQUIRED."),
        make_option(c("--REFERENCE_GENOME"), action="store", default=NULL,  type="character", help="reference genome version. default = hg19 "),
        make_option(c("--FORCE"           ), action="store", default=FALSE, type="logical",   help="RE-create pon data forcely")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR         = ARGS$BASE_DIR 
    NGS_LIB          = ARGS$NGS_LIB 
    REFERENCE_GENOME = ARGS$REFERENCE_GENOME  
    FORCE            = ARGS$FORCE  
#---------------------------------------------------------------------------------------------------
#---| MANUAL PARAMS |-------------------------------------------------------------------------------
    # BASE_DIR         = "/storage/home/kangsm/analysis/gcx.MRD/05.VALIDATION/VSID12_Liver.Cancer"
    # NGS_LIB          = "seqcap.ez.exome.v3.utr"
    # REFERENCE_GENOME = "hg19"
    # FORCE            = FALSE
    #---------------------------------------------------------------------------
    ## REQUIRED : ascat_ts_normal_sample.list -->> ${BASE_DIR}/cnv_pon/normal_sample.list
    ## REQUIRED : ascat_ts_worksheet.tsv 
    ## >> Patient_ID      Normal_ID       Normal_file     Gender
#--------------------------------------------------------------------------------------------------- 
#---| CHECK PARAMS |--------------------------------------------------------------------------------
    if(is.null(BASE_DIR)){ stop("No 'BASE_DIR' found. REQUIRED.") }
    if(is.null(NGS_LIB) ){ stop("No 'NGS_LIB' found. REQUIRED.")  }
    if(is.null(REFERENCE_GENOME)){ REFERENCE_GENOME = "hg19" }else{ REFERENCE_GENOME = REFERENCE_GENOME }
#--------------------------------------------------------------------------------------------------- 
#---| FOLDERS and RUN PARAMS |----------------------------------------------------------------------
    # ascat targeted-seq mode panel-of-normals result data folder
    PON_BASE_DIR       = "/storage/references_and_index/cnv/ascat/targeted_seq"
    ASCAT_CNV_PON_DIR  = sprintf("%s/%s", PON_BASE_DIR, NGS_LIB)
    if( dir.exists(ASCAT_CNV_PON_DIR) ){ 
        if( !FORCE )
        {
            stop("ASCAT Targeted-Seq Panel-of-Normals for Selected NGS LIBRARY is ALREADY EXISTS. If want to create again, Set 'FORCE=TRUE' ")
        }
    }else{
        system(sprintf("mkdir -p %s", ASCAT_CNV_PON_DIR)) 
    }
    # reference genome 
    if( REFERENCE_GENOME == "hg19" )
    {      ALLELE_PREFIX = "/storage/references_and_index/cnv/ascat/G1000_allelesAll_hg19/G1000_alleles_hg19_chr"
    }else{ ALLELE_PREFIX = "/storage/references_and_index/cnv/ascat/G1000_allelesAll_hg38/G1000_alleles_hg38_chr" }
    # bed file
    BED_FILE = sprintf("/storage/references_and_index/%s/bed/%s/%s_%s.target.bed", REFERENCE_GENOME, NGS_LIB, REFERENCE_GENOME, NGS_LIB)
    if( !file.exists(BED_FILE) ){ 
        stop("No 'BED_FILE' found. NEED FOR RUNNING.")
    }else{
        # chromosome ranges
        CHROM = unique(read.delim(BED_FILE, header=FALSE)[,1])
        CHROM = CHROM[ CHROM %nin%  c("X","Y","MT") ]
    }
    # normal sample list
    NormalSampleListFile = sprintf("%s/cnv_pon/normal_sample.list", BASE_DIR)
    if( !file.exists(NormalSampleListFile) ){ 
        stop("No Normal-Sample-List found. NEED FOR RUNNING.")
    }else{
        NORMAL_SAMPLE_LIST = read.table(NormalSampleListFile, header=FALSE,sep=" ")[,1:2]
    }
    # ascat input worksheet file
    PON_WORKSHEET  = sprintf("%s/cnv_pon/ascat_ts_workshhet.tsv", BASE_DIR)
    WorksheetTable = NORMAL_SAMPLE_LIST
    colnames(WorksheetTable) = c("Patient_ID","Normal_ID")
    WorksheetTable = WorksheetTable %>% mutate( 
        Normal_file = paste0( BASE_DIR, "/", Patient_ID, "/", Normal_ID, "/bam/", Normal_ID, ".analysisReady.bam"),
        Gender      = "XX"
    )
    write.table(WorksheetTable, PON_WORKSHEET, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
#---------------------------------------------------------------------------------------------------
#---| RUN ASCAT TARGETED-SEQ MODE PANEL-OF-NORMALS |------------------------------------------------
    ascat.prepareTargetedSeq(
        Worksheet         = PON_WORKSHEET,
        Workdir           = ASCAT_CNV_PON_DIR,
        alleles.prefix    = ALLELE_PREFIX,
        BED_file          = BED_FILE,
        allelecounter_exe = "/storage/apps/alleleCount-4.2.1/bin/alleleCounter",
        genomeVersion     = REFERENCE_GENOME,
        nthreads          = length(CHROM),
        minCounts         = 10,
        is_chr_based      = FALSE,
        chrom_names       = CHROM,
        min_base_qual     = 20,
        min_map_qual      = 20,
        plotQC            = TRUE
    )





