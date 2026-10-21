
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("RMySQL")
        library("optparse")
        library("ggplot2")
        library("superFreq")
        library("txdbmaker")
    })
    # scipen error patch
    {
        ga  <- getAnywhere("annotateSomaticQs")
        idx <- which(ga$where == "namespace:superFreq")
        f   <- ga$objs[[idx]]
        body_txt <- deparse(body(f))
        body_txt <- gsub(
            "options\\(scipen = sp\\)",
            "options(scipen = 999)",
            body_txt
        )
        body(f) <- parse(text = paste(body_txt, collapse = "\n"))[[1]]
        ns <- asNamespace("superFreq")
        unlockBinding("annotateSomaticQs", ns)
        assign("annotateSomaticQs", f, envir = ns)
        lockBinding("annotateSomaticQs", ns)
    }
    # makeTxDb patch
    {
        ns <- asNamespace("superFreq")
        h  <- "annotateQ"
        f <- get(h, envir = ns)
        txt <- deparse(body(f))
        txt <- gsub(
            "GenomicFeatures::makeTxDb",
            "txdbmaker::makeTxDb",
            txt,
            fixed = TRUE
        )
        body(f) <- parse(text = paste(txt, collapse = "\n"))[[1]]
        unlockBinding(h, ns)
        assign(h, f, envir = ns)
        lockBinding(h, ns)
    }
###
    BaseDir <- "/data/wgs"
    BatchID <- "WGS-26-01"
    AnalysisList <- "analysis.list"

    
###
    WorkDir <- sprintf("%s/%s/analysis/superfreq", BaseDir, BatchID)
    if( !dir.exists(WorkDir) ){ system(sprintf("mkdir -p %s/results", WorkDir)) }
    ainfo <- read.table(sprintf("%s/%s/meta/%s", BaseDir, BatchID, AnalysisList), header=F, sep=" ") 
    colnames(ainfo) <- c("SeqID","Group","Order")

    GENOME_FASTA           <- "/storage/references_and_index/hg19/fasta/hg19.primary_decoy_HLA_no.ALT.fa"
    SUPERFREQ_RESOURCE_DIR <- "/storage/references_and_index/cnv/superfreq/superFreqResources"
    NORMAL_REFERENCE_DIR   <- "/storage/references_and_index/cnv/superfreq/reference_normal_hg19_wgs"
    CAPTURE_REGION         <- "/storage/references_and_index/cnv/superfreq/captureRegions/hg19genome.10k.bed"

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

    bam_file_list <- sapply( ainfo$SeqID, function(sid) sprintf("%s/%s/%s/bam/%s.analysisReady.bam", BaseDir, BatchID, sid, sid))
    vcf_file_list <- sapply( ainfo$SeqID, function(sid) sprintf("%s/%s/%s/vcf/mutect2_tumor_only/%s.mutect2.Tonly.filtered.only.vcf", BaseDir, BatchID, sid, sid))
    name_list     <- ainfo$SeqID
    normal_tag    <- sapply(ainfo$Group, function(g) ifelse(g == "normal", "YES", "NO"))

    meta_data <- data.frame(
        BAM        = bam_file_list,
        VCF        = vcf_file_list,
        NAME       = name_list,
        INDIVIDUAL = BatchID,
        NORMAL     = normal_tag,
        TIMEPOINT  = c("T1","T2","N")
    )
    meta_data = meta_data[order(meta_data$TIMEPOINT), ]
    write.table(meta_data, sprintf("%s/%s.superfreq.meta.data.tsv", WorkDir, BatchID), row.names=F, col.names=T, sep="\t")

    RUN_superFeq = superFreq( 
        metaDataFile      = sprintf("%s/%s.superfreq.meta.data.tsv", WorkDir, BatchID), 
        Rdirectory        = WorkDir, 
        plotDirectory     = WorkDir,
        reference         = GENOME_FASTA, 
        genome            = "hg19",
        cpus              = 24, 
        mode              = "genome",
        captureRegions    = CAPTURE_REGION,
        normalDirectory   = NORMAL_REFERENCE_DIR,
        resourceDirectory = SUPERFREQ_RESOURCE_DIR,
        forceRedo         = forceRedoNothing() #forceRedoNothing() forceRedoEverything()
    )



