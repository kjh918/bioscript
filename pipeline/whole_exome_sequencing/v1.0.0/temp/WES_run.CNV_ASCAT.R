

#---| PACKAGES |-----------------------------------------------------------------------------------#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("parallel"))
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
    NORMAL_TYPE     <- "TS"
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--baseDir"), action="store", default=NA, type="character", help="BASE_DIR"),
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--tumorSeqID"), action="store", default="", type="character", help="TUMOR_SEQ_ID"),
        make_option(c("--normalSeqID"), action="store", default="", type="character", help="NORMAL_SEQ_ID"),
        make_option(c("--wesLibKit"), action="store", default=NA, type="character", help="BED_FILE"),
        make_option(c("--assembly"), action="store", default="hg19", type="character", help="GENOME_ASSEMBLY"),
        make_option(c("--threads"), action="store", default=15, type="double", help="THREADS"),
        make_option(c("--nType"), action="store", default="TS", type="character", help="MATCHED NORMAL TYPE. 'TS' or 'ORG'"),
        make_option(c("--superfreqSetID"), action="store", default="TumorOnly", type="character", help="superFreq set-ID")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR         <- ARGS$baseDir
    SEQ_FOLDER       <- ARGS$seqFolder
    TUMOR_SEQ_ID     <- ARGS$tumorSeqID
    NORMAL_SEQ_ID    <- ARGS$normalSeqID
    WES_LIB_KIT      <- ARGS$wesLibKit
    GENOME_ASSEMBLY  <- ARGS$assembly
    THREADS          <- as.numeric(ARGS$threads)
    NORMAL_TYPE      <- ARGS$nType
    SUPERFREQ_SET_ID <- ARGS$superfreqSetID
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR <- sprintf("%s/%s", BASE_DIR, SEQ_FOLDER)    
    LOG_DIR  <- sprintf("%s/%s/log", DATA_DIR, TUMOR_SEQ_ID)
    #--------------------------------------------------------------------------#
    LOG_FILE <- sprintf("%s/%s.ASCAT.%s.log", LOG_DIR, TUMOR_SEQ_ID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
    #--------------------------------------------------------------------------#
    ASCAT_BASE_DIR <- sprintf("%s/%s/cnv/ascat", DATA_DIR, TUMOR_SEQ_ID)
    if( !dir.exists(ASCAT_BASE_DIR) ){ system(sprintf("mkdir -p %s", ASCAT_BASE_DIR)) }
    #--------------------------------------------------------------------------#
    ASCAT_AF_DIR  <- sprintf("%s/alleleFreqs", ASCAT_BASE_DIR)
    if( !dir.exists(ASCAT_AF_DIR) ){ system(sprintf("mkdir -p %s", ASCAT_AF_DIR)) }
#--------------------------------------------------------------------------------------------------#

#---| ASCAT REFERENCES |---------------------------------------------------------------------------#
    if( GENOME_ASSEMBLY == "hg38" )
    {
        GENOME_FASTA     <- "/storage/references_and_index/hg38/fasta/Homo_sapiens_assembly38.fasta"
        REF_FLAT         <- "/storage/references_and_index/cnv/ascat/hg38.refFlat.txt"
        ASCAT_REF_ALLELE <- "/storage/references_and_index/cnv/ascat/G1000_allelesAll_hg38/G1000_alleles_hg38_chr"
        ASCAT_REF_LOCI   <- "/storage/references_and_index/cnv/ascat/G1000_lociAll_hg38/G1000_loci_hg38_chr"
        ASCAT_REF_GC     <- "/storage/references_and_index/cnv/ascat/GC_G1000_hg38.txt"
        ASCAT_REF_RT     <- "/storage/references_and_index/cnv/ascat/RT_G1000_hg38.txt"
    }else{
        GENOME_FASTA     <- "/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta"
        REF_FLAT         <- "/storage/references_and_index/cnv/ascat/hg19.refFlat.txt"
        ASCAT_REF_ALLELE <- "/storage/references_and_index/cnv/ascat/G1000_allelesAll_hg19/G1000_alleles_hg19_chr"
        ASCAT_REF_LOCI   <- "/storage/references_and_index/cnv/ascat/G1000_lociAll_hg19/G1000_loci_hg19_chr"
        ASCAT_REF_GC     <- "/storage/references_and_index/cnv/ascat/GC_G1000_hg19.txt"
        ASCAT_REF_RT     <- "/storage/references_and_index/cnv/ascat/RT_G1000_hg19.txt"
    }
    #--------------------------------------------------------------------------#
    BED_FILE <- sprintf("/storage/references_and_index/%s/bed/%s/%s_%s.target.bed", GENOME_ASSEMBLY, WES_LIB_KIT, GENOME_ASSEMBLY, WES_LIB_KIT)
    GENE_BED <- sprintf("/storage/references_and_index/%s/bed/%s/%s_%s.probe.cnvkit.gene.bed", GENOME_ASSEMBLY, WES_LIB_KIT, GENOME_ASSEMBLY, WES_LIB_KIT)
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES ASCAT: Purity,Ploidy,Sex,CNV log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE       : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER     : %s", SEQ_FOLDER), LOG_FILE, append=T) 
    write(sprintf(" TUMOR SEQ ID   : %s", TUMOR_SEQ_ID), LOG_FILE, append=T) 
    write(sprintf(" NORMAL SEQ ID  : %s", NORMAL_SEQ_ID), LOG_FILE, append=T) 
    write(sprintf(" NORMAL TYPE    : %s", NORMAL_TYPE), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

#==================================================================================================# 
    write(sprintf("%s | [ Preapre HTS-DATA ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ Preapre HTS-DATA ] START.")
#==================================================================================================# 

#---| RESULT TAG & FOLDER |------------------------------------------------------------------------#
    if( NORMAL_TYPE == "ORG" ){ RES_TAG = "ORG.NT" }else{ RES_TAG = "NT" }
    ASCAT_RES_DIR <- sprintf("%s/%s.%s.result", ASCAT_BASE_DIR, TUMOR_SEQ_ID, RES_TAG)
    if( ! dir.exists(ASCAT_RES_DIR) ){ system(sprintf("mkdir -p %s", ASCAT_RES_DIR)) }
#--------------------------------------------------------------------------------------------------#

#---| MOVE TO ALLELE-FREQ DIR |-------------------------------------------------------------------------#
    setwd(ASCAT_AF_DIR)
#--------------------------------------------------------------------------------------------------#

#---| SUPERFREQ RESULT FOR SEX |-------------------------------------------------------------------#
    SUPERFREQ_RES <- sprintf("%s/%s/cnv/superfreq/%s.%s.superFreq.summary.tsv", DATA_DIR, TUMOR_SEQ_ID, TUMOR_SEQ_ID, SUPERFREQ_SET_ID)   
    if( !file.exists(SUPERFREQ_RES) )
    { 
        SUPERFREQ_RES <- sprintf("%s/%s/cnv/superfreq/%s.%s.superFreq.summary.tsv", DATA_DIR, TUMOR_SEQ_ID, TUMOR_SEQ_ID, SUPERFREQ_SET_ID) 
        if( !file.exists(SUPERFREQ_RES) )
        {
            message("RUN superFreq First to determine sex of input")
            quit(save="no", status=0)
        }
        superFreqRes <- read.delim(SUPERFREQ_RES)
        superFreqRes <- superFreqRes[which(superFreqRes$NAME == TUMOR_SEQ_ID ), ]
    }else{
        superFreqRes <- read.delim(SUPERFREQ_RES)
    }
    INPUT_SEX    <- ifelse( superFreqRes$SEX == "male", "XY", "XX" )
    
#--------------------------------------------------------------------------------------------------#

#---| PREPARE HTS DATA |---------------------------------------------------------------------------#
        TUMOR_LOGR_FILE  <- sprintf("%s/%s.%s.Tumor.LogR.txt",    ASCAT_RES_DIR, TUMOR_SEQ_ID,  RES_TAG )
        TUMOR_BAF_FILE   <- sprintf("%s/%s.%s.Tumor.BAF.txt",     ASCAT_RES_DIR, TUMOR_SEQ_ID,  RES_TAG )
        NORMAL_LOGR_FILE <- sprintf("%s/%s.%s.Germline.LogR.txt", ASCAT_RES_DIR, NORMAL_SEQ_ID, RES_TAG )
        NORMAL_BAF_FILE  <- sprintf("%s/%s.%s.Germline.BAF.txt",  ASCAT_RES_DIR, NORMAL_SEQ_ID, RES_TAG )
        #----------------------------------------------------------------------#
        ascat.prepareHTS(
            tumourseqfile     = sprintf( "%s/%s/bam/%s.analysisReady.bam", DATA_DIR, TUMOR_SEQ_ID,  TUMOR_SEQ_ID ),
            normalseqfile     = sprintf( "%s/%s/bam/%s.analysisReady.bam", DATA_DIR, NORMAL_SEQ_ID, NORMAL_SEQ_ID),
            tumourname        = TUMOR_SEQ_ID,
            normalname        = NORMAL_SEQ_ID,
            BED_file          = BED_FILE,
            allelecounter_exe = "/storage/apps/alleleCount-4.2.1/bin/alleleCounter",
            alleles.prefix    = ASCAT_REF_ALLELE,
            loci.prefix       = ASCAT_REF_LOCI,
            genomeVersion     = GENOME_ASSEMBLY,
            nthreads          = THREADS,
            gender            = INPUT_SEX,
            tumourLogR_file   = TUMOR_LOGR_FILE,
            tumourBAF_file    = TUMOR_BAF_FILE,
            normalLogR_file   = NORMAL_LOGR_FILE,
            normalBAF_file    = NORMAL_BAF_FILE,
            #ref.fasta         = GENOME_FASTA,
            min_base_qual     = 20,
            min_map_qual      = 20
        )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | [ Preapre HTS-DATA ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("<<<<< [ Preapre HTS-DATA ] FINISHED.")
    write(sprintf("%s | [ Load HTS-DATA ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ Load HTS-DATA ] START.")
#==================================================================================================# 

#---| MOVE TO RESULT DIR |-------------------------------------------------------------------------#
    setwd(ASCAT_RES_DIR)
#--------------------------------------------------------------------------------------------------#

#---| LOADING ASCAT DATA |-------------------------------------------------------------------------#
    ascat.bc <- ascat.loadData(
        Tumor_LogR_file    = TUMOR_LOGR_FILE, 
        Tumor_BAF_file     = TUMOR_BAF_FILE, 
        Germline_LogR_file = NORMAL_LOGR_FILE, 
        Germline_BAF_file  = NORMAL_BAF_FILE, 
        genomeVersion      = GENOME_ASSEMBLY, 
        gender             = INPUT_SEX,
        isTargetedSeq      = TRUE
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | [ Load HTS-DATA ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("<<<<< [ Load HTS-DATA ] FINISHED.")
    write(sprintf("%s | [ CORRECTION ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s |   - Create pre-correction plots Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ CORRECTION ] START.")
#==================================================================================================# 

#---| CREATE PLOTS OF BEFORE-CORRECTION |----------------------------------------------------------#
    ascat.plotRawData( 
        ascat.bc, 
        img.dir       = ASCAT_RES_DIR,
        img.prefix    = sprintf("pre-correction.%s.", RES_TAG),
        logr.y_values = c(-3, 3) 
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s |   - Pre-correction plots created.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s |   - Perform correction Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
#==================================================================================================# 

#---| CORRECTION |---------------------------------------------------------------------------------# 
    ascat.bc <- ascat.correctLogR( 
        ascat.bc, 
        GCcontentfile    = ASCAT_REF_GC, 
        replictimingfile = ASCAT_REF_RT
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s |   - Perform correction Finished.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s |   - Create post-correction plots Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
#==================================================================================================# 

#---| CREATE PLOTS OF AFTER-CORRECTION |-----------------------------------------------------------#
    ascat.plotRawData(
        ascat.bc, 
        img.dir       = ASCAT_RES_DIR,
        img.prefix    = sprintf("post-correction.%s.", RES_TAG),
        logr.y_values = c(-3, 3) 
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s |   - Post-correction plots created.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s | [ CORRECTION ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("<<<<< [ CORRECTION ] FINISHED.")
    write(sprintf("%s | [ SEGMENTATION ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s |   - Perform segmentation Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ SEGMENTATION ] START.")
#==================================================================================================# 

#---| SEGMENTATION |-------------------------------------------------------------------------------#
    ascat.bc <- ascat.aspcf( ascat.bc, penalty = 70 )  # <- create BAF.PCFed.txt & LogR.PCFed.txt
    #--------------------------------------------------------------------------#
    # rename_baf_cmd  <- sprintf("mv %s/%s.BAF.PCFed.txt %s/%s.%s.BAF.PCFed.txt",   ASCAT_RES_DIR, TUMOR_SEQ_ID, ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
    # rename_logr_cmd <- sprintf("mv %s/%s.LogR.PCFed.txt %s/%s.%s.LogR.PCFed.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID, ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
    # system(rename_baf_cmd)
    # system(rename_logr_cmd)
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s |   - Perform segmentation Finished.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s |   - Create segment-plots Start.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
#==================================================================================================# 

#---| SEGMENT PLOTS |------------------------------------------------------------------------------#
    ascat.plotSegmentedData( ascat.bc )
    #--------------------------------------------------------------------------#
    # rename_aspcf_cmd  <- sprintf("mv %s/%s.ASPCF.png %s/%s.%s.ASPCF.png", ASCAT_RES_DIR, TUMOR_SEQ_ID, ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
    # system(rename_aspcf_cmd)
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s |   - Segment-plots created.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s | [ SEGMENTATION ] FINISHED", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("<<<<< [ SEGMENTATION ] FINISHED.")
    write(sprintf("%s | [ ASCAT-ANALYSIS ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ ASCAT-ANALYSIS ] START.")
#==================================================================================================# 

#---| RUN ASCAT ANALYSIS |-------------------------------------------------------------------------#
    ascat.output = ascat.runAscat(
        ascat.bc, 
        gamma          = 1, 
        y_limit        = 9,
        write_segments = TRUE,
        pdfPlot        = TRUE
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | [ ASCAT-ANALYSIS ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("<<<<< [ ASCAT-ANALYSIS ] FINISHED.")
    write(sprintf("%s | [ CREATE QC and SAVE ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ CREATE QC and SAVE ] START.")
#==================================================================================================# 

#---| ASCAT QC RESULTS |---------------------------------------------------------------------------#
    ascat.qc <- ascat.metrics(
        ascat.bc,
        ascat.output
    )
#--------------------------------------------------------------------------------------------------#

#---| SAVE RESULTS |-------------------------------------------------------------------------------#
    ASCAT_RESULT_RDATA <- sprintf("%s/%s.%s.ASCAT.results.Rdata", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
    #--------------------------------------------------------------------------#
    save( ascat.bc, ascat.output, ascat.qc, file = ASCAT_RESULT_RDATA )
    #--------------------------------------------------------------------------#
    write.table(
        ascat.qc, 
        sprintf("%s/%s.%s.ASCAT.qc.result.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG ), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    #--------------------------------------------------------------------------#
    qc.simple            <- data.frame(ID=rownames(ascat.qc), ascat.qc[,c("purity", "sex", "ploidy")], ploidy_int="") 
    qc.simple$sex        <- sapply( qc.simple$sex, function(sx) ifelse( sx == "XX", "female", "male" ) )
    qc.simple$ploidy_int <- round( qc.simple$ploidy )
    #--------------------------------------------------------------------------#
    write.table(
        qc.simple, sprintf("%s/%s.%s.ASCAT.summary.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG ), 
        quote=F, col.names=F, row.names=F, sep=" "
    )
    write.table(
        qc.simple, sprintf("%s/%s.%s.ASCAT.summary.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG ), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    #--------------------------------------------------------------------------#
    write.table(
        qc.simple, sprintf("%s/%s.%s.ASCAT.summary.txt", ASCAT_BASE_DIR, TUMOR_SEQ_ID, RES_TAG ), 
        quote=F, col.names=F, row.names=F, sep=" "
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | [ CREATE QC and SAVE ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("<<<<< [ CREATE QC and SAVE ] FINISHED.")
    write(sprintf("%s | [ Post-Analysis of ASCAT ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ ASCAT POST-ANALYSIS ] START.")
#==================================================================================================# 

#---| ADDITIONAL PACKAGES |------------------------------------------------------------------------#
    library(plyr)
    library(dplyr)
    library(Hmisc)
    library(ggplot2)
#--------------------------------------------------------------------------------------------------#

#---| FUNCTIONS FOR REFORMAT-RESULTS |-------------------------------------------------------------#
    #' @description add genomic region coordicates 
    #' @param ASCAT_LOGR       ASCAT logRatio results table
    #' @param CHR_LENGTH_TABLE genomic region coordinates table
    #' @export
        AddChrCoords <- function( ASCAT_LOGR, CHR_LENGTH_TABLE ) 
        {
            ASCAT_LOGR$position2 <- apply( ASCAT_LOGR[,c("Chromosome","Position")], 1, function(Z) 
            { as.numeric(Z[2]) + as.numeric(CHR_LENGTH_TABLE[which(CHR_LENGTH_TABLE$chr == Z[1]), "add_X"]) })
            return(ASCAT_LOGR)
        }
    #----------------------------------------------------------------------------------------------#
    #' @description merge ASCAT segmentation result and genomic region coordinates to logRatio result
    #' @param ASCAT_LOGR       ASCAT logRatio results
    #' @param ASCAT_SEGMENT    ASCAT segmentation results 
    #' @param CHR_LENGTH_TABLE genomic region coordinates table
    #' @export
        MergeSegmentsToLogRatio <- function( ASCAT_LOGR, ASCAT_SEGMENT, CHR_LENGTH_TABLE, PLOIDY )
        {
            ASCAT_LOGR$CN = NA
            
            for( k in 1:nrow(ASCAT_SEGMENT) )
            {
                ASCAT_LOGR[which(
                    ASCAT_LOGR$Position   >= ASCAT_SEGMENT[k, "startpos"] &
                    ASCAT_LOGR$Position   <= ASCAT_SEGMENT[k, "endpos"  ] &
                    ASCAT_LOGR$Chromosome == ASCAT_SEGMENT[k, "chr"     ] 
                ), "CN" ] <- ( ASCAT_SEGMENT[k, "nMajor"] + ASCAT_SEGMENT[k, "nMinor"] )
            }

            ASCAT_LOGR$CNA_ABS <- sapply( ASCAT_LOGR$CN, function(n) 
            {
                if( n == 0         ){ cna1 <- "DELETION"      }
                if( n == 1         ){ cna1 <- "LOSS"          }
                if( n == 2         ){ cna1 <- "NEUTRAL"       }
                if( n  > 2 & n < 8 ){ cna1 <- "GAIN"          }
                if( n >= 8         ){ cna1 <- "AMPLIFICATION" }
                return(cna1)
            })

            ASCAT_LOGR$CNA_PLOIDY <- sapply( ASCAT_LOGR$CN, function(n) 
            {
                if( n == 0                     ){ cna2 <- "DELETION"      }
                if( n >= 1 & n < PLOIDY        ){ cna2 <- "LOSS"          }
                if( n == PLOIDY                ){ cna2 <- "NEUTRAL"       }
                if( n  > PLOIDY & n < PLOIDY*4 ){ cna2 <- "GAIN"          }
                if( n >= PLOIDY*4              ){ cna2 <- "AMPLIFICATION" }
                return(cna2)
            })

            ASCAT_LOGR <- AddChrCoords( ASCAT_LOGR = ASCAT_LOGR, CHR_LENGTH_TABLE = CHR_LENGTH_TABLE )

            return(ASCAT_LOGR)
        }  
    #----------------------------------------------------------------------------------------------#
    #' @description indexing list generation
    #' @export
        indexing = function(x){ sapply(unique(x), function(z) list(z)) }    
    #----------------------------------------------------------------------------------------------#
    #' @description make gene-level copy number alteration results from ASCAT result
    #' @param AscatRes ASCAT result
    #' @param GeneBed  WES library BED file including gene-annotation 
    #' @param Threads  multicore N-cpu
    #' @param Ploidy   ploidy as integer
    #' @export
        AscatGeneLevelCna <- function( AscatRes = ascat_res, GeneBed = gene_bed, Threads = THREADS, Ploidy = INT_PLOIDY  )
        {
            bedGenes <- mclapply( indexing(AscatRes$X), function(asct)
            {
                loc = unlist(strsplit(asct, "_"))
                regionGenes <- GeneBed[which(
                    GeneBed$start <= as.numeric(loc[2]) &
                    GeneBed$end   >= as.numeric(loc[2]) & 
                    GeneBed$chr   == as.character(loc[1])
                ), "gene"] 
                RG = paste(unique(unlist(strsplit(regionGenes, ","))), collapse=",")
                return(RG)
            }, mc.cores = Threads)
            
            CNA_RES = AscatRes %>% 
                left_join(., (ldply(bedGenes, .id = "X") %>% 
                mutate( X = as.character(X) )), by="X") %>%
                dplyr::rename( GENE = V1 )

            CNA_RES_GENE = unique(CNA_RES[,c("Chromosome","CN","CNA_ABS","CNA_PLOIDY","GENE")])

            GENE_CNA = ldply(mclapply( indexing(1:nrow(CNA_RES_GENE)), function(n) 
            {   
                genecn <- CNA_RES_GENE[n, ]
                genelist <- unlist(strsplit(genecn[, "GENE"], ","))
                if( length(genelist) > 1 )
                {
                    GCN_rev = data.frame( genecn[,c("Chromosome","CN","CNA_ABS","CNA_PLOIDY")], GENE=genelist )
                }else{
                    GCN_rev = genecn
                }
                return(GCN_rev)
            }, mc.cores = Threads)) %>% 
                filter( GENE %in% cnv_ref_genes$gene ) %>% unique() 
            
            GIDX <- indexing( GENE_CNA$GENE )

            GENE_CNA_rev <- mclapply( GIDX, function(gn) 
            {
                gcn <- GENE_CNA %>% filter( GENE == gn ) 

                if( nrow(gcn) > 1 )
                {   
                    colOrders = colnames(gcn)
                    copies = gcn$CN

                    maxAbsCopy = ifelse( length(copies[ copies >= 2 ]) > 0, max(copies[ copies >= 2 ]), NA )
                    minAbsCopy = ifelse( length(copies[ copies  < 2 ]) > 0, min(copies[ copies  < 2 ]), NA )

                    if( is.na(maxAbsCopy) )
                    { 
                        maxAbsCNA = NA 
                    }else{
                        if( maxAbsCopy == 2 ){ maxAbsCNA = "NEUTRAL" }else{ maxAbsCNA = ifelse( maxAbsCopy >= 8, "AMPLIFICATION", "GAIN" ) }
                    }

                    if( is.na(minAbsCopy) )
                    {
                        minAbsCNA = NA
                    }else{
                        minAbsCNA = ifelse( minAbsCopy == 0 , "DELETION", "LOSS" )
                    }

                    maxPldCopy = ifelse( length(copies[ copies >= Ploidy ]) > 0, max(copies[ copies >= Ploidy ]), NA )
                    minPldCopy = ifelse( length(copies[ copies  < Ploidy ]) > 0, min(copies[ copies  < Ploidy ]), NA )

                    if( is.na(maxPldCopy) )
                    { 
                        maxPldCNA = NA
                    }else{
                        if( maxPldCopy == Ploidy)
                        { maxPldCNA = "NEUTRAL" }else{ maxPldCNA = ifelse( maxPldCopy >= Ploidy*4, "AMPLIFICATION", "GAIN") }
                    }

                    if( is.na(minPldCopy) )
                    {
                        minPldCNA = NA
                    }else{
                        minPldCNA = ifelse( minPldCopy == 0 , "DELETION", "LOSS" ) 
                    }

                    abs_cna_res = c(maxAbsCNA, minAbsCNA)
                    if( length(abs_cna_res[!is.na(abs_cna_res)]) < 2 ){ abs_cna_res = abs_cna_res[!is.na(abs_cna_res)] }

                    ploidy_cna_res = c(maxPldCNA, minPldCNA)
                    if( length(ploidy_cna_res[!is.na(ploidy_cna_res)]) < 2 ){ ploidy_cna_res = ploidy_cna_res[!is.na(ploidy_cna_res)] }

                    abs_cna_res = paste(abs_cna_res, collapse=)

                    abs_cn = c(maxAbsCopy, minAbsCopy)
                    abs_cn = abs_cn[!is.na(abs_cn)]

                    gcn2 = data.frame(
                        unique(gcn[,c("Chromosome","GENE")]),
                        CN           = abs_cn,
                        CNA_ABS = abs_cna_res,
                        CNA_PLOIDY   = ploidy_cna_res
                    )

                    gcn = gcn2[ which(!is.na(gcn2$CN)), colOrders ]
                }

                return(gcn)

            }, mc.cores = Threads)

            GENE_CNA_RESULTS = ldply(GENE_CNA_rev)[,-1] 
            
            return(GENE_CNA_RESULTS)
        }
    #----------------------------------------------------------------------------------------------#
    #' @description draw CNA class distribution bar plot
    #' @param GENE_CNA_DATA gene-level ASCAT CNA results
    #' @param GENE_TYPE     gene group type. 'gcx': gcx cancer genes, 'oc' : oncogenes, 'ts' : tumor suppressor genes
    #' @param CNA_TYPE      CNA-class. 'abs': absolute CN-based, 'pld': ploidy CN-based
    #' @export
        AscatCNABarplot <- function( GENE_CNA_DATA, GENE_TYPE, CNA_TYPE )
        {
            if( GENE_TYPE == "gcx" )
            {
                if( CNA_TYPE ==  "abs" )
                {
                    GENE_CNA_DATA <- GENE_CNA_DATA %>% dplyr::rename( CNA = CNA_ABS, GROUP = GCX_CancerGenes ) 
                    TITLE         <- "Abs CNA: GCX-CancerGenes"
                }else{
                    GENE_CNA_DATA <- GENE_CNA_DATA %>% dplyr::rename( CNA = CNA_PLOIDY, GROUP = GCX_CancerGenes ) 
                    TITLE         <- "Ploidy CNA: GCX-CancerGenes"
                }
                X.max = 203
                             
            }
            if( GENE_TYPE == "oc" )
            {
                if( CNA_TYPE ==  "abs" )
                {
                    GENE_CNA_DATA <- GENE_CNA_DATA %>% dplyr::rename( CNA = CNA_ABS, GROUP = OncoGenes ) 
                    TITLE         <- "Abs CNA: Oncogenes"
                }else{
                    GENE_CNA_DATA <- GENE_CNA_DATA %>% dplyr::rename( CNA = CNA_PLOIDY, GROUP = OncoGenes ) 
                    TITLE         <- "Ploidy CNA: GCX-CancerGenes"
                }   
                X.max = 418          
            }
            if( GENE_TYPE == "ts" )
            {
                if( CNA_TYPE ==  "abs" )
                {
                    GENE_CNA_DATA <- GENE_CNA_DATA %>% dplyr::rename( CNA = CNA_ABS, GROUP = TumorSuppressorGenes ) 
                    TITLE         <- "Abs CNA: TumorSuppressor Genes"
                }else{
                    GENE_CNA_DATA <- GENE_CNA_DATA %>% dplyr::rename( CNA = CNA_PLOIDY, GROUP = TumorSuppressorGenes ) 
                    TITLE         <- "Ploidy CNA: TumorSuppressor Genes"
                }     
                X.max = 365                
            }
            #------------------------------------------------------------------#
            alter_levels = c("LOSS","GAIN","DELETION","AMPLIFICATION")
            alter_colors = c("DELETION"="#1a5387", "LOSS"="#8ca9c3", "NEUTRAL"="#4c4c4c", "GAIN"="#cf9890", "AMPLIFICATION"="#9f3122")
            #------------------------------------------------------------------#
            pdata <- GENE_CNA_DATA %>% filter( GROUP == 1 ) %>% group_by( CNA ) %>% reframe( counts = length(unique(GENE)) ) 
            missing_cna <- alter_levels[ alter_levels %nin% pdata$CNA ]
            if( length(missing_cna) > 0 )
            {
                pdata <- rbind(
                    pdata, 
                    data.frame( CNA = alter_levels[ alter_levels %nin% pdata$CNA ], counts = 0 )
                )
            }
            cnabarp <- pdata %>%
                mutate( CNA = factor(CNA, levels = alter_levels) ) %>%
                ggplot() + 
                    coord_cartesian( xlim = c(0, X.max) ) +
                    geom_bar( aes( x = counts, y = CNA, fill = CNA ), stat="identity", position = "dodge", width = 0.8) +
                    geom_text( aes( x = counts, y = CNA, label = counts ), hjust = -1, size = 3) +
                    scale_fill_manual( values = alter_colors ) +
                    geom_vline( xintercept = 0, colour = "#6f6d6d") +
                    labs( x = "", y = "CNA") +
                    theme(
                        panel.background = element_rect(fill='white'), 
                        panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth = 0.5),
                        plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                        plot.title   = element_text(face = "bold", size = 11),
                        #axis.line.y  = element_line(colour = "#6f6d6d"),
                        axis.line.x  = element_line(colour = "#6f6d6d"),
                        axis.text.y  = element_text(face = "bold", size = 10),
                        axis.text.x  = element_text(face = "bold", size = 10),
                        axis.title.x = element_text(face = "bold", size = 10, margin = margin(t = 5)),
                        axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 5)),
                        axis.ticks.y = element_blank(), 
                        legend.position="none"
                    ) +
                    ggtitle( TITLE )
            return(cnabarp)
        }
    #----------------------------------------------------------------------------------------------#
    #' @description draw CNV plot
    #' @param SEGMENT          ASCAT segment, logRatio, baf result table
    #' @param CHR_LENGTH_TABLE genomic region coordinates table
    #' @param CHR              chromosome for plotting
    #' @param GENOME_ASSEMBLY  genome assembly version
    #' @param CNA_TYPE         CNA-class
    #' @param Show_BAF         add BAF to plot or not. default = FALSE
    #' @export
        AscatCNVPlot <- function( SEGMENT , CHR_LENGTH_TABLE, CHR = "all", GENOME_ASSEMBLY, CNA_TYPE, Show_BAF = FALSE )
        {
            # CNV COLORS
            CNA.COLORS <- c("DELETION"="#14426c", "LOSS"="#6C8EBF", "NEUTRAL"="#4c4c4c", "GAIN"="#B85450", "AMPLIFICATION"="#7f271b")
            #
            SEGMENT <- SEGMENT %>% dplyr::rename( logR = 4 )
            if( CNA_TYPE == "abs" ){
                SEGMENT <- SEGMENT %>% dplyr::rename( CNA = CNA_ABS )
            }else{
                SEGMENT <- SEGMENT %>% dplyr::rename( CNA = CNA_PLOIDY )
            }
            
            log2Outlier = 3
            SEGMENT[which(SEGMENT$logR >  log2Outlier), "logR"] <-  log2Outlier
            SEGMENT[which(SEGMENT$logR < -log2Outlier), "logR"] <- -log2Outlier
            
            # PLOT PARAMETERS
            if( CHR == "all"     ){ CHROMOSOME = unique(SEGMENT$Chromosome) }else{ CHROMOSOME = CHR }
            if( CHR == "all"     ){ PS    = 0.1 }else{ PS    = 0.3 }
            if( CHR == "all"     ){ chrFS = 5   }else{ chrFS = 12  }
            if( CHR == "all"     ){ ChromosomeTitle = "All Chromosomes" }else{ ChromosomeTitle = CHR }
            if( CNA_TYPE == "abs" ){ PlotType = "Abs CNA" }else{ PlotType = "Ploidy CNA" }

            if( Show_BAF )
            {
                breakPoint.Y    <- c(-1,-0,1,1.5,2)
                limit.Y         <- c(-1.2, 2.2)
                label.Y         <- rep("", length(breakPoint.Y))
            }else{
                breakPoint.Y    <- c(-1,-0,1)
                limit.Y         <- c(-1.2, 1.2)
                label.Y         <- rep("", length(breakPoint.Y))
            }
            breakPoint.X    <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$add_X
            label.X         <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$chr
            chrEndPoint.X   <- chrLengthTable %>% filter( chr %in% CHROMOSOME ) %>% .$cum_chr_length
            chrStartPoint.X <- ifelse( CHR == "all", 0, chrLengthTable %>% filter( chr %in% CHR ) %>% .$add_X %>% head(1))
            
            if( GENOME_ASSEMBLY == "hg19" & CHR != "all" )
            { label.X <- paste0("chr", label.X ) ; ChromosomeTitle <- paste0("chr", ChromosomeTitle ) }

            # CNV PLOT
            cnvPlot <- SEGMENT %>% filter( Chromosome %in% CHROMOSOME ) %>% ggplot() +
                coord_cartesian( ylim = limit.Y ) +
                theme(
                    panel.background = element_rect(fill='white'),
                    axis.line.y      = element_line(colour = "#7a7a7a"),
                    axis.ticks.x     = element_blank(),
                    axis.text.x      = element_text(hjust = 0, angle = 30, vjust=1, size=chrFS ),
                    legend.position  = "none"
                ) +
                labs( x = "Chromosome", y = "", color = "CNA" ) +
                geom_hline( yintercept = 0, col = "#7a7a7a", alpha = 0.5 ) +
                geom_point( aes( x = position2, y = logR/3 , color = CNA ), size = PS, alpha = 0.1 ) +
                scale_color_manual(values = CNA.COLORS)

            if( Show_BAF )
            {
                cnvPlot <- cnvPlot + 
                    geom_hline( yintercept = c(-1,1,1.5,2), color = "#b7b7b7", alpha = 0.2, linetype = 'dashed' ) +
                    geom_point( aes( x = position2, y = BAF+1), size = 0.1, alpha = 0.6, color = "#cccccc")
            }else{
                cnvPlot <- cnvPlot + 
                    geom_hline( yintercept = c(-1,1), color = "#b7b7b7", alpha = 0.2, linetype = 'dashed' )
            }
                
            cnvPlot <- cnvPlot  +
                scale_y_continuous( breaks = breakPoint.Y, labels = label.Y ) +
                scale_x_continuous( breaks = breakPoint.X, labels = label.X ) +
                geom_vline( xintercept = chrStartPoint.X, col = "#7a7a7a", alpha = 0.5, linewidth = 0.2, linetype = 'dashed') +
                geom_vline( xintercept = chrEndPoint.X,   col = "#7a7a7a", alpha = 0.5, linewidth = 0.2, linetype = 'dashed' ) +
                ggtitle( paste0( PlotType, " : ", ChromosomeTitle) )

            return(cnvPlot)
        }
#--------------------------------------------------------------------------------------------------#

#---| PRE-DEFINED DATA |---------------------------------------------------------------------------#
    # CHR LENGTH
    load("/storage/home/kangsm/myDB/rds/chr_length.Rdata")                     # CHR_LIST,hg19_chr_length,hg38_chr_length"
    #--------------------------------------------------------------------------#
    # CNV REF GENES
    load("/storage/home/kangsm/myDB/rds/cnv_refGenes_list_hg19_hg38.Rdata")    # CNV_REF_GENES
    # CANCER GENES
    load("/storage/home/kangsm/myDB/rds/PresetGenes.Rdata")                    # c11_preset_genes, gcx_preset_genes
    #--------------------------------------------------------------------------#
    if( GENOME_ASSEMBLY == "hg38" ){ chrLengthTable <- hg38_chr_length    }else{ chrLengthTable <- hg19_chr_length }
    if( GENOME_ASSEMBLY == "hg38" ){ cnv_ref_genes  <- CNV_REF_GENES$hg38 }else{ cnv_ref_genes  <- CNV_REF_GENES$hg19 }
    #--------------------------------------------------------------------------#
    GeneBedFile <- sprintf("/storage/references_and_index/%s/bed/%s/hg19_twist.exome.2.0.probe.cnvkit.gene.bed", GENOME_ASSEMBLY, WES_LIB_KIT)
#--------------------------------------------------------------------------------------------------#

#---| GENE LEVEL CNA |-----------------------------------------------------------------------------#
    ascat_logr    <- read.delim(sprintf("%s/%s.%s.Tumor.LogR.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG))
    ascat_segment <- read.delim(sprintf("%s/%s.segments.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID))
    ascat_call    <- read.delim(sprintf("%s/%s.LogR.PCFed.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID), header=FALSE) %>% dplyr::rename( "X" = 1, "mean_logR" = 2 )
    ascat_baf     <- read.delim(sprintf("%s/%s.%s.Tumor.BAF.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG))
    INT_PLOIDY    <- qc.simple$ploidy_int    
    if( qc.simple$ploidy > 1 & qc.simple$ploidy < 2.5 ){ INT_PLOIDY = 2 }
    #--------------------------------------------------------------------------#
    ascat_res <- MergeSegmentsToLogRatio( 
        ASCAT_LOGR       = ascat_logr, 
        ASCAT_SEGMENT    = ascat_segment, 
        CHR_LENGTH_TABLE = chrLengthTable, 
        PLOIDY           = INT_PLOIDY 
    ) %>% left_join(., ascat_call, by="X")
    #--------------------------------------------------------------------------#
    gene_bed  <- read.delim(GENE_BED, header=F) %>% dplyr::rename( chr=1, start=2, end=3, gene=4) 
    #--------------------------------------------------------------------------#
    ascat.gene.level.cna <- AscatGeneLevelCna( 
        AscatRes = ascat_res, 
        GeneBed  = gene_bed, 
        Threads  = THREADS,
        Ploidy   = INT_PLOIDY 
    )
    #--------------------------------------------------------------------------#
    ascat.gene.level.cna$OncoGenes = ascat.gene.level.cna$TumorSuppressorGenes = ascat.gene.level.cna$GCX_CancerGenes = 0
    ascat.gene.level.cna[which(ascat.gene.level.cna$GENE %in% gcx_preset_genes$oncokb_oncogene     ), "OncoGenes"            ] = 1
    ascat.gene.level.cna[which(ascat.gene.level.cna$GENE %in% gcx_preset_genes$oncokb_tsgene       ), "TumorSuppressorGenes" ] = 1
    ascat.gene.level.cna[which(ascat.gene.level.cna$GENE %in% gcx_preset_genes$gcx_cancer_essential), "GCX_CancerGenes"      ] = 1
    #--------------------------------------------------------------------------#
    ascat.GCNA1 <- ascat.gene.level.cna
    ascat.GCNA2 <- ascat.gene.level.cna %>% 
        filter( CNA_ABS != "NEUTRAL" | CNA_PLOIDY != "NEUTRAL" )
    ascat.GCNA3 <- ascat.gene.level.cna %>% 
        filter( CNA_ABS != "NEUTRAL" | CNA_PLOIDY != "NEUTRAL" ) %>%
        filter( GCX_CancerGenes != 0 | TumorSuppressorGenes != 0 | OncoGenes != 0 ) 
    #--------------------------------------------------------------------------#
    write.table( ascat.GCNA1, 
        sprintf("%s/%s.%s.ASCAT.gene.level.CNA.raw.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG),
        quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
    )
    write.table( ascat.GCNA2, 
        sprintf("%s/%s.%s.ASCAT.gene.level.CNA.altered.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG),
        quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
    )
    write.table( ascat.GCNA3, 
        sprintf("%s/%s.%s.ASCAT.gene.level.CNA.altered.cancer.genes.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG),
        quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
    )
#--------------------------------------------------------------------------------------------------#

#---| DISTRIBUTION FIGURES |-----------------------------------------------------------------------#
    Gene_CN_stats <- as.data.frame( table(ascat.GCNA1$CN) ) %>% dplyr::rename( CN = 1, GENES = 2) %>% 
        mutate( GENES_pct = round(GENES/sum(GENES)*100, 1) , ID = TUMOR_SEQ_ID )
    #--------------------------------------------------------------------------#
    CN.DistributePlot <- ggplot(Gene_CN_stats) +
        geom_line( aes( x = CN, y = GENES_pct, group = ID), linewidth = 1, color = "#82B366", alpha = 0.5 ) +
        geom_point( aes( x = CN, y = GENES_pct), size = 5, color = "#B85450" ) +
        geom_text(aes(x=CN, y=GENES_pct, label=GENES), vjust=-1.5) +
        coord_cartesian( ylim = c(0,100) ) +
        labs( x = "Estimated CN", y = "Genes (%)" ) +
        theme(
            panel.background = element_rect(fill='white'), 
            panel.grid.major = element_line(color = '#dcdbdb', linetype = 'dotted', linewidth = 0.5),
            plot.margin  = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            plot.title   = element_text(face = "bold", size = 10),
            axis.line.x  = element_line(colour = "#6f6d6d"),
            axis.line.y  = element_line(colour = "#6f6d6d"),
            axis.text.x  = element_text(face = "bold", size = 12),
            axis.text.y  = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 5)),
            axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
            axis.ticks.x = element_blank()
        ) + 
        ggtitle( "Distribution of Estimated CopyNumbers" )
    #--------------------------------------------------------------------------#
    ggsave(CN.DistributePlot, width=5, height=5, unit='in', dpi=200, type='cairo',
        file = sprintf("%s/%s.%s.ASCAT.estimated.CN.distribution.png", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG )  )
#--------------------------------------------------------------------------------------------------#

#---| CNA-CLASS DISTRIBUTION BARPLOT |-------------------------------------------------------------#
    CNA.abs.barplots <- cowplot::plot_grid(
        AscatCNABarplot( GENE_CNA_DATA = ascat.GCNA3, GENE_TYPE = "gcx", CNA_TYPE = "abs" ),
        AscatCNABarplot( GENE_CNA_DATA = ascat.GCNA3, GENE_TYPE = "oc", CNA_TYPE = "abs" ),
        AscatCNABarplot( GENE_CNA_DATA = ascat.GCNA3, GENE_TYPE = "ts", CNA_TYPE = "abs" ),
        nrow = 1
    )
    CNA.ploidy.barplots <- cowplot::plot_grid(
        AscatCNABarplot( GENE_CNA_DATA = ascat.GCNA3, GENE_TYPE = "gcx", CNA_TYPE = "pld" ),
        AscatCNABarplot( GENE_CNA_DATA = ascat.GCNA3, GENE_TYPE = "oc", CNA_TYPE = "pld" ),
        AscatCNABarplot( GENE_CNA_DATA = ascat.GCNA3, GENE_TYPE = "ts", CNA_TYPE = "pld" ),
        nrow = 1
    )
    ggsave(CNA.abs.barplots, width=14, height=2.5, unit='in', dpi=150, type='cairo',
        file = sprintf("%s/%s.%s.ASCAT.abs.CNA.class.distribution.png", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG )  )
    ggsave(CNA.ploidy.barplots, width=14, height=2.5, unit='in', dpi=150, type='cairo',
        file = sprintf("%s/%s.%s.ASCAT.ploidy.CNA.class.distribution.png", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG )  )
#--------------------------------------------------------------------------------------------------#

#---| CNV PLTOS |----------------------------------------------------------------------------------#
    ascat_res$BAF <- ascat_baf[match( ascat_res$X, ascat_baf[,1]), 4]
    #--------------------------------------------------------------------------#
    write.table(ascat_res, 
        sprintf("%s/%s.%s.ASCAT.segments.CNA.with.baf.raw.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG),
        quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t"
    )
    #--------------------------------------------------------------------------#
    cnvPlot_dir <- sprintf("%s/cnvplots", ASCAT_RES_DIR)
    if( !dir.exists(cnvPlot_dir) ){ system(sprintf("mkdir -p %s", cnvPlot_dir)) }
    #--------------------------------------------------------------------------#
    CHR_LIST <- c(1:22, "X", "Y")
    if( GENOME_ASSEMBLY == "hg38" ){ CHR_LIST <- paste0("chr", CHR_LIST) }
    
    for( CHR in c("all", CHR_LIST) )
    {
        if( CHR == "all" ){ FW = 12 ; FH = 3 }else{ FW = 8 ; FH = 3 }

        chr.cnvplot <- AscatCNVPlot( 
            SEGMENT          = ascat_res , 
            CHR_LENGTH_TABLE = chrLengthTable, 
            CHR              = CHR, 
            GENOME_ASSEMBLY  = GENOME_ASSEMBLY, 
            CNA_TYPE         = "abs", 
            Show_BAF         = TRUE
        )
        ggsave( chr.cnvplot,  width=FW, height=FH, unit='in', dpi=150, type='cairo',
            file = sprintf("%s/%s.%s.chr.%s.CNV.plot.with.BAF.png", cnvPlot_dir, TUMOR_SEQ_ID, RES_TAG, CHR )  
        )

        chr.cnvplot2 <- AscatCNVPlot( 
            SEGMENT          = ascat_res , 
            CHR_LENGTH_TABLE = chrLengthTable, 
            CHR              = CHR, 
            GENOME_ASSEMBLY  = GENOME_ASSEMBLY, 
            CNA_TYPE         = "abs", 
            Show_BAF         = FALSE
        )
        ggsave( chr.cnvplot2,  width=FW, height=FH, unit='in', dpi=150, type='cairo',
            file = sprintf("%s/%s.%s.chr.%s.CNV.plot.no.BAF.png", cnvPlot_dir, TUMOR_SEQ_ID, RES_TAG, CHR )  
        )

    }
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#
    write(sprintf("%s | [ Post-Analysis of ASCAT ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ ASCAT POST-ANALYSIS ] FINISHED.")
    write("   ", LOG_FILE, append=T)
    write(sprintf("%s | WES ASCAT: Purity,Ploidy,Sex,CNV DONE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
    write("   ", LOG_FILE, append=T)
#==================================================================================================# 
