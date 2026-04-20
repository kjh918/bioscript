
#---| metaDataFile |-------------------------------------------------------------------------------# 
# BAM        : path to bam file of sample.  An index file is required as well.
# VCF        : path to vcf file of sample.  This should include both somatic and germline variants. 
# NAME       : unique identifier of the sample.
# INDIVIDUAL : unique identifier of the indiviudal the sample comes from.  
#              This information is used pair up matched normal samples for somatic mutation filtering
#              This column is also used to group samples that should be compared to each other
# NORMAL     : Should be YES if the sample is normal, NO otherwise. A normal sample is assumed to have no somatic mutations. 
# TIMEPOINT  : Mostly used as label in plots together with the sample. Can be Diagnosis, Relapse, resistant or similar labels. 
#              Does not influence the analyss otherwise.
#--------------------------------------------------------------------------------------------------#

#---| PACKAGES |-----------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("superFreq"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Cairo"))
grDevices::X11.options(type='cairo')
options(device='x11')
options(bitmapType='cairo')
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR        <- "/data/wes"
    GENOME_ASSEMBLY <- "hg19"
    WES_LIB_KIT     <- "twist.exome.2.0"
    TUMOR_SEQ_ID    <- ""
    NORMAL_SEQ_ID   <- ""
    ORG_SEQ_ID      <- ""
    MODE            <- "exome"
    THREADS         <- 15
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
      make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER"),   
      make_option(c("--tumorSeqID"), action="store", default="", type="character", help="TUMOR_SEQ_ID"),
      make_option(c("--normalSeqID"), action="store", default="", type="character", help="NORMAL_SEQ_ID"),
      make_option(c("--organoidSeqID"), action="store", default="", type="character", help="ORG_SEQ_ID"),
      make_option(c("--baseDir"), action="store", default=NA, type="character", help="BASE_DIR"),
      make_option(c("--assembly"), action="store", default="hg19", type="character", help="GENOME_ASSEMBLY"),
      make_option(c("--wesLibKit"), action="store", default="twist.exome.2.0", type="character", help="WES_LIB_KIT"),
      make_option(c("--mode"), action="store", default="exome", type="character", help="run mode : exome"),
      make_option(c("--threads"), action="store", default=15, type="double", help="THREADS"),
      make_option(c("--nType"), action="store", default=15, type="character", help="Type of Normal Sample"),
      make_option(c("--setID"), action="store", default=15, type="character", help="result folder name"),
      make_option(c("--runDir"), action="store", default=NULL, type="character", help="result base folder")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    SEQ_FOLDER        <- ARGS$seqFolder
    TUMOR_SEQ_ID      <- ARGS$tumorSeqID
    NORMAL_SEQ_ID     <- ARGS$normalSeqID
    ORG_SEQ_ID        <- ARGS$organoidSeqID
    BASE_DIR          <- ARGS$baseDir
    GENOME_ASSEMBLY   <- ARGS$assembly
    WES_LIB_KIT       <- ARGS$wesLibKit
    THREADS           <- as.numeric(ARGS$threads)
    NORMAL_TYPE       <- ARGS$nType
    SET_ID            <- ARGS$setID
    RUN_BASE_DIR      <- ARGS$runDir
#--------------------------------------------------------------------------------------------------#

#---| FUNCTIONS |----------------------------------------------------------------------------------#
    Calculate.superFreq.Ploidy = function( SegRes )
    {
        require(dplyr)
        require(Hmisc)
        SegRes = SegRes %>% filter( chr %nin% c("X","Y","chrX","chrY") ) %>% mutate( W = x2 - x1 )
        PLOIDY = weighted.mean( 2*2^(SegRes$M), SegRes$W)
        return(PLOIDY)
    }
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR = sprintf("%s/%s", BASE_DIR, SEQ_FOLDER)
    #--------------------------------------------------------------------------#
    if(is.null(RUN_BASE_DIR)){
        RUN_DIR = sprintf("%s/%s/cnv/superfreq", DATA_DIR, TUMOR_SEQ_ID)
        system(sprintf("mkdir -p %s", RUN_DIR))
    }else{
        RUN_DIR = RUN_BASE_DIR
        system(sprintf("mkdir -p %s", RUN_DIR))
    }
#--------------------------------------------------------------------------------------------------#

#---| VCF TAG |------------------------------------------------------------------------------------#
    VCF_TAG="mutect2.keep.germline.filtered"
#--------------------------------------------------------------------------------------------------#

#---| SET 'SET_ID' |----------------------------------------------------------------------------#
    # if( SET_ID %in% c("tonly","Tonly") )
    # {
    #     SET_ID = sprintf("%s.%s", TUMOR_SEQ_ID, "Tonly")
    # }else{
    #     if( NORMAL_TYPE == "ORG" )
    #     {
    #         SET_ID = sprintf("%s.%s", TUMOR_SEQ_ID, "ORG.NT")
    #     }else{
    #         SET_ID = sprintf("%s.%s", TUMOR_SEQ_ID, "NT")
    #     }
    # }
#--------------------------------------------------------------------------------------------------#

#---| META DATA CREATE |---------------------------------------------------------------------------#
    if( TUMOR_SEQ_ID  != "" )
    {
        TUMOR_BAM  = sprintf("%s/%s/bam/%s.analysisReady.bam",    DATA_DIR, TUMOR_SEQ_ID, TUMOR_SEQ_ID)
        TUMOR_VCF  = sprintf("%s/%s/vcf/%s.%s.vcf", DATA_DIR, TUMOR_SEQ_ID, TUMOR_SEQ_ID, VCF_TAG)
    }else{
        TUMOR_BAM = TUMOR_VCF = ""
    }
    #--------------------------------------------------------------------------# 
    if( NORMAL_SEQ_ID != "" )
    { 
        NORMAL_BAM = sprintf("%s/%s/bam/%s.analysisReady.bam",    DATA_DIR, NORMAL_SEQ_ID, NORMAL_SEQ_ID)
        NORMAL_VCF = sprintf("%s/%s/vcf/%s.%s.vcf", DATA_DIR, NORMAL_SEQ_ID, NORMAL_SEQ_ID, VCF_TAG)
    }else{
        NORMAL_BAM = NORMAL_VCF = ""
    }
    #--------------------------------------------------------------------------# 
    if( ORG_SEQ_ID != "" )
    { 
        ORG_BAM = sprintf("%s/%s/bam/%s.analysisReady.bam",    DATA_DIR, ORG_SEQ_ID, ORG_SEQ_ID)
        ORG_VCF = sprintf("%s/%s/vcf/%s.%s.vcf", DATA_DIR, ORG_SEQ_ID, ORG_SEQ_ID, VCF_TAG)
    }else{
        ORG_BAM = ORG_VCF = ""
    }
    #--------------------------------------------------------------------------#
    BAM_LIST = c(NORMAL_BAM, TUMOR_BAM, ORG_BAM) 
    VCF_LIST = c(NORMAL_VCF, TUMOR_VCF, ORG_VCF)    
    #--------------------------------------------------------------------------#
    metaInfo = data.frame(
        BAM        = BAM_LIST,
        VCF        = VCF_LIST,
        NAME       = c(NORMAL_SEQ_ID, TUMOR_SEQ_ID, ORG_SEQ_ID),
        INDIVIDUAL = SET_ID,
        NORMAL     = c("YES","NO","NO"),
        TIMEPOINT  = c("1","2","3")
    )
    metaInfo = metaInfo[which(metaInfo$NAME != ""), ]
    #--------------------------------------------------------------------------#
    if( TUMOR_SEQ_ID == "" & NORMAL_SEQ_ID == "" & ORG_SEQ_ID == "" )
    {
        message("No INPUT SEQ_ID !!")
        quit( save="no", status = 0 )
    }
    #--------------------------------------------------------------------------#
    if( nrow(metaInfo) == 1 & metaInfo[1, "NORMAL"] == "YES" ) { metaInfo[1, "NORMAL"] = "NO" }
    #--------------------------------------------------------------------------#
    write.table(metaInfo, sprintf("%s/%s_metaData.tsv", RUN_DIR, SET_ID), quote=F, col.names=T, row.names=F, sep="\t")
#--------------------------------------------------------------------------------------------------#

#---| RUN PARAMETERS |-----------------------------------------------------------------------------#
    if( GENOME_ASSEMBLY == "hg38" ) {
        GENOME_FASTA = "/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta"
    }else{
        GENOME_FASTA = "/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta"
    }
    #--------------------------------------------------------------------------#
    SUPERFREQ_RESOURCE_DIR = "/storage/references_and_index/cnv/superfreq"
#--------------------------------------------------------------------------------------------------#


#---| RUN SUPERFREQ |------------------------------------------------------------------------------#

    RUN_superFeq = superFreq( 
        metaDataFile      = sprintf("%s/%s_metaData.tsv", RUN_DIR, SET_ID), 
        captureRegions    = "/storage/references_and_index/cnv/superfreq/captureRegions/hg19exons.bed",
        normalDirectory   = sprintf("%s/%s_%s_pon_bam", SUPERFREQ_RESOURCE_DIR, GENOME_ASSEMBLY, WES_LIB_KIT ),
        Rdirectory        = RUN_DIR, 
        plotDirectory     = sprintf("%s/%s", RUN_DIR, SET_ID),
        reference         = GENOME_FASTA, 
        genome            = GENOME_ASSEMBLY,
        cpus              = THREADS, 
        mode              = MODE,
        resourceDirectory = sprintf("%s/superFreqResources", SUPERFREQ_RESOURCE_DIR),
        forceRedo         = forceRedoNothing() #forceRedoNothing() forceRedoEverything()
    )

#--------------------------------------------------------------------------------------------------#

#---| SUPERFREQ SUMMARY FOR CNV ANALYSIS |---------------------------------------------------------#
    PLOIDY_DATA_DIR = sprintf("%s/%s/%s/data", RUN_DIR, SET_ID, SET_ID)
    #--------------------------------------------------------------------------#
    SuperFreqRes <- data.frame(NAME=metaInfo$NAME, PLOIDY=NA, PLOIDY_INT=NA)
    for( n in 1:nrow(SuperFreqRes) )
    {
        SegRes <- read.delim(sprintf("%s/CNAsegments_%s.tsv", PLOIDY_DATA_DIR, SuperFreqRes[n,"NAME"]))
        ploidy <- Calculate.superFreq.Ploidy( SegRes )
        SuperFreqRes[n, "PLOIDY"]     = round(as.numeric(ploidy), 2)
        SuperFreqRes[n, "PLOIDY_INT"] = round(as.numeric(ploidy))
    }
    #--------------------------------------------------------------------------#
    load(sprintf("%s/%s/fit.Rdata", RUN_DIR, SET_ID)) # fit
    sex_res          <- fit$fit$sex
    names(sex_res)   <- gsub("\\-normal", "", names(sex_res))
    SuperFreqRes$SEX <- sex_res[ SuperFreqRes$NAME ]
    #--------------------------------------------------------------------------#
    load(sprintf("%s/%s/stories.Rdata", RUN_DIR, SET_ID)) # stories
    PR                  <- stories$stories[[SET_ID]]$clusters$cloneStories %>% filter( call != "germline" ) 
    purity              <- apply(PR$stories, 2, function(y) max(y))
    SuperFreqRes$PURITY <- round(purity[ SuperFreqRes$NAME ], 2)
    #--------------------------------------------------------------------------#
    write.table(SuperFreqRes, sprintf("%s/%s.superFreq.summary.txt", RUN_DIR, SET_ID), quote=F, col.names=F, row.names=F, sep=" ")
    write.table(SuperFreqRes, sprintf("%s/%s.superFreq.summary.tsv", RUN_DIR, SET_ID), quote=F, col.names=T, row.names=F, sep="\t")

#--------------------------------------------------------------------------------------------------#



# xs = chrToX(fCsExon$annotation$Chr, fCsExon$annotation$Start)

# weight = fCsExon$counts[,NORMAL_SEQ_ID]

# weightInRegion = function(x1, x2) sum(weight[xs > x1 & xs < x2])
# weightsInRegion = function(xs1, xs2) sapply(1:length(xs1), function(i) weightInRegion(xs1[i], xs2[i]))


# weight = fCsExon$counts[,TUMOR_SEQ_ID]

# cnvCancer = clusters[[ TUMOR_SEQ_ID ]]$clusters
# cancerPloidy = 2*sum(weightsInRegion(cnvCancer$x1,cnvCancer$x2)*2^cnvCancer$M)/sum(weight)




# extractPloidy = function(truthData) {
#   xs = chrToX(truthData$fCsExon$annotation$Chr, truthData$fCsExon$annotation$Start)
#   weight = truthData$fCsExon$counts[,'normal1']
#   weightInRegion = function(x1, x2) sum(weight[xs > x1 & xs < x2])
#   weightsInRegion = function(xs1, xs2) sapply(1:length(xs1), function(i) weightInRegion(xs1[i], xs2[i]))
  
#   cnvCancer = truthData$clusters$cancer$clusters
#   cancerPloidy = 2*sum(weightsInRegion(cnvCancer$x1,cnvCancer$x2)*2^cnvCancer$M)/sum(weight)
#   cnvNormal1 = truthData$clusters$normal1$clusters
#   normal1Ploidy = 2*sum(weightsInRegion(cnvNormal1$x1,cnvNormal1$x2)*2^cnvNormal1$M)/sum(weight)
#   cnvNormal2 = truthData$clusters$normal2$clusters
#   normal2Ploidy = 2*sum(weightsInRegion(cnvNormal2$x1,cnvNormal2$x2)*2^cnvNormal2$M)/sum(weight)
#   return(list('normal1'=normal1Ploidy, 'normal2'=normal2Ploidy, 'cancer'=cancerPloidy))
# }

