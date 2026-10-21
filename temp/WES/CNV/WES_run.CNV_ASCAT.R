

#---| PACKAGES |-----------------------------------------------------------------------------------#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("Cairo"))
grDevices::X11.options(type='cairo')
options(device='x11')
options(stringsAsFactors=FALSE) 
options(bitmapType='cairo')
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR        <- "/data/wes"
    GENOME_ASSEMBLY <- "hg19"
    WES_LIB_KIT     <- "twist.exome.2.0"
    TUMOR_SEQ_ID    <- ""
    NORMAL_SEQ_ID   <- ""
    THREADS         <- 15
    RUN_MODE        <- "tonly"
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--baseDir"), action="store", default=NA, type="character", help="BASE_DIR"),
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--tumorSeqID"), action="store", default="", type="character", help="TUMOR_SEQ_ID"),
        make_option(c("--normalSeqID"), action="store", default="", type="character", help="NORMAL_SEQ_ID"),
        make_option(c("--bed"), action="store", default=NA, type="character", help="SET_ID"),
        make_option(c("--assembly"), action="store", default="hg19", type="character", help="GENOME_ASSEMBLY"),
        make_option(c("--wesLibKit"), action="store", default="twist.exome.2.0", type="character", help="WES_LIB_KIT"),
        make_option(c("--threads"), action="store", default=15, type="double", help="THREADS"),
        make_option(c("--runMode"), action="store", default="tonly", type="character", help="RUN_MODE")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR        <- ARGS$baseDir
    SEQ_FOLDER      <- ARGS$seqFolder
    TUMOR_SEQ_ID    <- ARGS$tumorSeqID
    NORMAL_SEQ_ID   <- ARGS$normalSeqID
    BED_FILE        <- ARGS$bed
    GENOME_ASSEMBLY <- ARGS$assembly
    WES_LIB_KIT     <- ARGS$wesLibKit
    THREADS         <- as.numeric(ARGS$threads)
    RUN_MODE        <- ARGS$runMode
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR = sprintf("%s/%s", BASE_DIR, SEQ_FOLDER)
    #--------------------------------------------------------------------------#
    ASCAT_RES_DIR = sprintf("%s/%s/cnv/ascat", DATA_DIR, TUMOR_SEQ_ID)
    if( !dir.exists(ASCAT_RES_DIR) ){ system(sprintf("mkdir -p %s", ASCAT_RES_DIR)) }
    #--------------------------------------------------------------------------#
    

#--------------------------------------------------------------------------------------------------#
    TUMOR_BAM_FILE  <- sprintf( "%s/%s/bam/%s.analysisReady.bam", DATA_DIR, TUMOR_SEQ_ID,  TUMOR_SEQ_ID )
    NORMAL_BAM_FILE <- sprintf( "%s/%s/bam/%s.analysisReady.bam", DATA_DIR, NORMAL_SEQ_ID, NORMAL_SEQ_ID)
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
    if( GENOME_ASSEMBLY == "hg38" )
    {
        GENOME_FASTA     <- "/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta"
        ASCAT_REF_ALLELE <- "/storage/home/kangsm/myDB/cnv_ascat/hg38/wes/G1000_allelesAll_hg38/G1000_alleles_hg38_chr"
        ASCAT_REF_LOCI   <- "/storage/home/kangsm/myDB/cnv_ascat/hg38/wes/G1000_lociAll_hg38/G1000_loci_hg38_chr"
        ASCAT_REF_GC     <- "/storage/home/kangsm/myDB/cnv_ascat/hg38/wes/GC_G1000_hg38.txt"
        ASCAT_REF_RT     <- "/storage/home/kangsm/myDB/cnv_ascat/hg38/wes/RT_G1000_hg38.txt"
    }else{
        GENOME_FASTA     <- "/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta"
        ASCAT_REF_ALLELE <- "/storage/home/kangsm/myDB/cnv_ascat/hg19/wes/G1000_allelesAll_hg19/G1000_alleles_hg19_chr"
        ASCAT_REF_LOCI   <- "/storage/home/kangsm/myDB/cnv_ascat/hg19/wes/G1000_lociAll_hg19/G1000_loci_hg19_chr"
        ASCAT_REF_GC     <- "/storage/home/kangsm/myDB/cnv_ascat/hg19/wes/GC_G1000_hg19.txt"
        ASCAT_REF_RT     <- "/storage/home/kangsm/myDB/cnv_ascat/hg19/wes/RT_G1000_hg19.txt"
    }
#--------------------------------------------------------------------------------------------------#

#---| PREPARE HTS DATA |---------------------------------------------------------------------------#
    if( RUN_MODE == "nt" ) 
    {
        TUMOR_LOGR_FILE  <- sprintf("%s/%s.NT.Tumor.LogR.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID)
        TUMOR_BAF_FILE   <- sprintf("%s/%s.NT.Tumor.BAF.txt",  ASCAT_RES_DIR, TUMOR_SEQ_ID)
        NORMAL_LOGR_FILE <- sprintf("%s/%s.NT.Germline.LogR.txt", ASCAT_RES_DIR, NORMAL_SEQ_ID)
        NORMAL_BAF_FILE  <- sprintf("%s/%s.NT.Germline.BAF.txt",  ASCAT_RES_DIR, NORMAL_SEQ_ID)

        ascat.prepareHTS(
            tumourseqfile     = TUMOR_BAM_FILE,
            normalseqfile     = NORMAL_BAM_FILE,
            tumourname        = TUMOR_SEQ_ID,
            normalname        = NORMAL_SEQ_ID,
            BED_file          = BED_FILE,
            allelecounter_exe = "/storage/apps/alleleCount-4.2.1/bin/alleleCounter",
            alleles.prefix    = ASCAT_REF_ALLELE,
            loci.prefix       = ASCAT_REF_LOCI,
            gender            = "XX",
            genomeVersion     = GENOME_ASSEMBLY,
            nthreads          = THREADS,
            tumourLogR_file   = TUMOR_LOGR_FILE,
            tumourBAF_file    = TUMOR_BAF_FILE,
            normalLogR_file   = NORMAL_LOGR_FILE,
            normalBAF_file    = NORMAL_BAF_FILE,
            ref.fasta         = GENOME_FASTA
        )
    }else{
        TUMOR_LOGR_FILE  <- sprintf("%s/%s.NT.Tumor.LogR.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID)
        TUMOR_BAF_FILE   <- sprintf("%s/%s.NT.Tumor.BAF.txt",  ASCAT_RES_DIR, TUMOR_SEQ_ID)

        ascat.prepareHTS(
            tumourseqfile     = TUMOR_BAM_FILE,
            tumourname        = TUMOR_SEQ_ID,
            BED_file          = BED_FILE,
            allelecounter_exe = "/storage/apps/alleleCount-4.2.1/bin/alleleCounter",
            alleles.prefix    = ASCAT_REF_ALLELE,
            loci.prefix       = ASCAT_REF_LOCI,
            gender            = "XX",
            genomeVersion     = GENOME_ASSEMBLY,
            nthreads          = THREADS,
            tumourLogR_file   = TUMOR_LOGR_FILE,
            tumourBAF_file    = TUMOR_BAF_FILE,
            ref.fasta         = GENOME_FASTA
        )
    }
#--------------------------------------------------------------------------------------------------#

#---| LOADING ASCAT DATA |-------------------------------------------------------------------------#
    if( RUN_MODE == "nt" ) 
    {
        ascat.bc <- ascat.loadData(
            Tumor_LogR_file    = TUMOR_LOGR_FILE, 
            Tumor_BAF_file     = TUMOR_BAF_FILE, 
            Germline_LogR_file = NORMAL_LOGR_FILE, 
            Germline_BAF_file  = NORMAL_BAF_FILE, 
            gender             = 'XX', 
            genomeVersion      = GENOME_ASSEMBLY, 
            isTargetedSeq      = TRUE
        )
    }else{
        ascat.bc <- ascat.loadData(
            Tumor_LogR_file    = TUMOR_LOGR_FILE, 
            Tumor_BAF_file     = TUMOR_BAF_FILE, 
            gender             = rep('XX',100), 
            genomeVersion      = GENOME_ASSEMBLY, 
            isTargetedSeq      = TRUE
        )
    }
#--------------------------------------------------------------------------------------------------#

#---| CREATE PLOTS OF BEFORE-CORRECTION |----------------------------------------------------------#
    ascat.plotRawData( 
        ascat.bc, 
        img.prefix    = sprintf("%s/PreCorrection.", ASCAT_RES_DIR),
        logr.y_values = c(-3, 3) 
    )
#--------------------------------------------------------------------------------------------------#

#---| CORRECTION |---------------------------------------------------------------------------------# 
    ascat.bc <- ascat.correctLogR( 
        ascat.bc, 
        GCcontentfile    = ASCAT_REF_GC, 
        replictimingfile = ASCAT_REF_RT
    )
#--------------------------------------------------------------------------------------------------#

#---| CREATE PLOTS OF AFTER-CORRECTION |----------------------------------------------------------#
    ascat.plotRawData(
        ascat.bc, 
        img.prefix    = sprintf("%s/PostCorrection.", ASCAT_RES_DIR),
        logr.y_values = c(-3, 3) 
    )
#--------------------------------------------------------------------------------------------------#

#---| SEGMENTATION |-------------------------------------------------------------------------------#
    if( RUN_MODE == "nt" ) 
    {
        ascat.bc <- ascat.aspcf(
            ascat.bc, 
            penalty=25
        )
    }else{
        GG <- ascat.predictGermlineGenotypes( ascat.bc )
        #----------------------------------------------------------------------#
        ascat.bc <- ascat.aspcf(
            ascat.bc, 
            penalty  = 25,
            ascat.gg = GG
        )
    }
#--------------------------------------------------------------------------------------------------#

#---| SEGMENT PLOTS |------------------------------------------------------------------------------#
    ascat.plotSegmentedData( ascat.bc )
#--------------------------------------------------------------------------------------------------#

#---| RUN ASCAT ANALYSIS |-------------------------------------------------------------------------#
    ascat.output = ascat.runAscat(
        ascat.bc, 
        gamma          = 1, 
        y_limit        = 9
        write_segments = TRUE,
        pdfPlot        = TRUE
    )
#--------------------------------------------------------------------------------------------------#

#---| ASCAT QC RESULTS |---------------------------------------------------------------------------#
    ascat.qc <- ascat.metrics(
        ascat.bc,
        ascat.output
    )
#--------------------------------------------------------------------------------------------------#

#---| SAVE RESULTS |-------------------------------------------------------------------------------#
    ASCAT_RESULT_RDATA <- sprintf("%s/%s.ASCAT.results.Rdata", ASCAT_RES_DIR, TUMOR_SEQ_ID)
    #--------------------------------------------------------------------------#
    save( ascat.bc, ascat.output, ascat.qc, file = ASCAT_RESULT_RDATA )
    #--------------------------------------------------------------------------#
    write.table(
        ascat.qc, 
        sprintf("%s/%s.ASCAT.qc.result.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID ), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    #--------------------------------------------------------------------------#
    qc.simple            <- data.frame(ID=rownames(ascat.qc), ascat.qc[,c("purity", "sex", "ploidy")], ploidy_int="") 
    qc.simple$sex        <- sapply( qc.simple$sex, function(sx) ifelse( sx == "XX", "female", "male" ) )
    qc.simple$ploidy_int <- round( qc.simple$ploidy+0.01, 0 )
    #--------------------------------------------------------------------------#
    write.table(
        qc.simple, 
        sprintf("%s/%s.ASCAT.purity.sex.ploidy.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID ), 
        quote=F, col.names=F, row.names=F, sep=" "
    )
#--------------------------------------------------------------------------------------------------#



