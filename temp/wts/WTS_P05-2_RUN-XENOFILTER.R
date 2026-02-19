### Packages ---------------------------------------------------------------------------------------    
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("BiocParallel"))
    suppressPackageStartupMessages(library("XenofilteR"))
###
### Option Arguments -------------------------------------------------------------------------------
    option_list = list(
        make_option(c("--baseDir"), action="store", default=NULL, type="character", help="Human-reference aligned BAM file (absolute path)"),
        make_option(c("--batchID"), action="store", default=NULL, type="character", help="Human-reference aligned BAM file (absolute path)"),
        make_option(c("--seqID"), action="store", default=NULL, type="character", help="Human-reference aligned BAM file (absolute path)"),
        make_option(c("--threads"), action="store", default=NULL, type="character", help="Human-reference aligned BAM file (absolute path)")
    )
    ARGS = parse_args(OptionParser(option_list=option_list))
    BaseDir <- ARGS$baseDir
    BatchID <- ARGS$batchID
    SeqID   <- ARGS$seqID
    Threads <- ARGS$threads
###
### Param Check ------------------------------------------------------------------------------------
    if( is.null(BaseDir) ){ stop("No 'BaseDir'.") }
    if( is.null(BatchID) ){ stop("No 'BatchID'.") }
    if( is.null(SeqID)   ){ stop("No 'SeqID'.")   }
    if( is.null(Threads) ){ Threads  <- 4 }else{ Threads <- as.numeric(Threads) }
###
### Manual Params (debug or manual run) ------------------------------------------------------------
    # BaseDir = "/data/wts"
    # BatchID = "WTS_25_02"
    # SeqID   = "WTS_25_02_03"
    # Threads = 10
###
### Parallel Processing Cores ----------------------------------------------------------------------
    RunCores <- SnowParam(workers=Threads, type="SOCK")
###
### Data & Folders ---------------------------------------------------------------------------------
    BamDir   <- sprintf("%s/%s/%s/bam", BaseDir, BatchID, SeqID)
    FastqDir <- sprintf("%s/%s/%s/fastq", BaseDir, BatchID, SeqID)
    InputHumanBamFile <- sprintf("%s/%s.human.Aligned.out.Sorted.bam", BamDir, SeqID)
    InputMouseBamFile <- sprintf("%s/%s.mouse.Aligned.out.Sorted.bam", BamDir, SeqID)
    OutputPrefixHuman <- sprintf("%s.human", SeqID)
    OutputPrefixMouse <- sprintf("%s.mouse", SeqID)
###
### Sample List ------------------------------------------------------------------------------------
    SampleList_Run1 <- data.frame( Graft = InputHumanBamFile, Host = InputMouseBamFile )
    SampleList_Run2 <- data.frame( Graft = InputMouseBamFile, Host = InputHumanBamFile )
###
### RUN1 - Human BAM Generation --------------------------------------------------------------------
    message(sprintf("XenofilteR Run-1 : Human BAM Generation"))
    XenofilteR(
        sample.list        = SampleList_Run1,
        destination.folder = BamDir,
        bp.param           = RunCores,
        output.names       = OutputPrefixHuman,
        MM_threshold       = 6,
        Unmapped_penalty   = 8
    )
    XenofilterHumanBam <- sprintf("%s/Filtered_bams/%s_Filtered.bam", BamDir, OutputPrefixHuman)
    system(sprintf("mv %s %s/%s.human.xenofilter.bam", XenofilterHumanBam, BamDir, SeqID))
    system(sprintf("mv %s.bai %s/%s.human.xenofilter.bam.bai", XenofilterHumanBam, BamDir, SeqID))
    system(sprintf("rm -rf  %s/Filtered_bams", BamDir))
    message(sprintf("XenofilteR Run-1 Done."))
###
### RUN2 - Mouse BAM Generation --------------------------------------------------------------------
    message(sprintf("XenofilteR Run-2 : MOuse BAM Generation"))
    XenofilteR(
        sample.list        = SampleList_Run2,
        destination.folder = BamDir,
        bp.param           = RunCores,
        output.names       = OutputPrefixMouse,
        MM_threshold       = 6,
        Unmapped_penalty   = 8
    )
    XenofilterMouseBam <- sprintf("%s/Filtered_bams/%s_Filtered.bam", BamDir, OutputPrefixMouse)
    system(sprintf("mv %s %s/%s.mouse.xenofilter.bam", XenofilterMouseBam, BamDir, SeqID))
    system(sprintf("mv %s.bai %s/%s.mouse.xenofilter.bam.bai", XenofilterMouseBam, BamDir, SeqID))
    system(sprintf("rm -rf  %s/Filtered_bams", BamDir))
    message(sprintf("XenofilteR Run-2 Done."))
###
### Post-Processing (Bam-to-Fastq) -----------------------------------------------------------------
    # Human
    XfInitialBamHuman <- sprintf("%s/%s.human.xenofilter.bam", BamDir, SeqID)
    XfPairedBamHuman  <- sprintf("%s/%s.human.xenofilter.paired.bam", BamDir, SeqID)
    XfSortedBamHuman  <- sprintf("%s/%s.human.xenofilter.paired.sorted.bam", BamDir, SeqID)
    FastqPrefixHuman  <- sprintf("%s/%s.human.trimmed", FastqDir, SeqID)
    system(sprintf("samtools view -@ %s -b -f 2 %s > %s", Threads, XfInitialBamHuman, XfPairedBamHuman))
    system(sprintf("samtools sort -@ %s -n -o %s %s", Threads, XfSortedBamHuman, XfPairedBamHuman))       
    system(sprintf("samtools fastq -@ %s -n -1 %s_R1.fastq.gz -2 %s_R2.fastq.gz -0 /dev/null -s /dev/null %s", 
        Threads, FastqPrefixHuman, FastqPrefixHuman, XfSortedBamHuman
    )) 

    # Mouse
    XfInitialBamMouse <- sprintf("%s/%s.mouse.xenofilter.bam", BamDir, SeqID)
    XfPairedBamMouse  <- sprintf("%s/%s.mouse.xenofilter.paired.bam", BamDir, SeqID)
    XfSortedBamMouse  <- sprintf("%s/%s.mouse.xenofilter.paired.sorted.bam", BamDir, SeqID)
    FastqPrefixMouse  <- sprintf("%s/%s.mouse.trimmed", FastqDir, SeqID)
    system(sprintf("samtools view -@ %s -b -f 2 %s > %s", Threads, XfInitialBamMouse, XfPairedBamMouse))
    system(sprintf("samtools sort -@ %s -n -o %s %s", Threads, XfSortedBamMouse, XfPairedBamMouse))       
    system(sprintf("samtools fastq -@ %s -n -1 %s_R1.fastq.gz -2 %s_R2.fastq.gz -0 /dev/null -s /dev/null %s", 
        Threads, FastqPrefixMouse, FastqPrefixMouse, XfSortedBamMouse
    )) 
###

### Final message ----------------------------------------------------------------------------------
    message(sprintf("XenofilteR PROCESSING FINISHED."))
###=================================================================================================