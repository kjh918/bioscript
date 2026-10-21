
### Pakages ----------------------------------------------------------------------------------------
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("setwidth")
        library("optparse")
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("RMySQL")
        library("GenomicRanges")
        library("openxlsx")
    })
#
### Parse options ----------------------------------------------------------------------------------
    option_list = list(
        make_option(c("--baseDir"), action="store", default=NULL, type="character", help="base working folder [REQUIRED]"),
        make_option(c("--batchID"), action="store", default=NULL, type="character", help="analysis batch-id folder  [REQUIRED]"),
        make_option(c("--seqID"),   action="store", default=NULL, type="character", help="sample-id  [REQUIRED]"),
        make_option(c("--genomeAssembly"), action="store", default=NULL, type="character", help="reference genome version. default = 'hg38'"),
        make_option(c("--circleMap"), action="store", default=TRUE, type="logical", help="run for circleMap result. default = TRUE"),
        make_option(c("--circleFinder"), action="store", default=TRUE, type="logical", help="run for circle_finder result. default = TRUE")
    )
    ARGS = parse_args(OptionParser(option_list=option_list))
    #...........................................................................
    BaseDir         <- ARGS$baseDir
    BatchID         <- ARGS$batchID
    SeqID           <- ARGS$seqID
    GenomeAssembly  <- ARGS$genomeAssembly
    RunCircleMap    <- ARGS$circleMap
    RunCircleFinder <- ARGS$circleFinder
#
### Manual Params ----------------------------------------------------------------------------------
    # BaseDir <- "/data/wgs"
    # BatchID <- "cancer_cell_lines"
    # SeqID   <- "SRR22509543_HeLa"
    # GenomeAssembly <- "hg38"
#
### Param Check ------------------------------------------------------------------------------------
    if(is.null(BaseDir)        ){  stop("No 'baseDir'. REQUIRED.") }
    if(is.null(BatchID)        ){  stop("No 'batchID'. REQUIRED.") }
    if(is.null(SeqID)          ){  stop("No 'seqID'. REQUIRED.")   }
    if(is.null(GenomeAssembly) ){  GenomeAssembly  <- "hg38" }
    if(is.null(RunCircleMap)   ){  RunCircleMap    <- TRUE   }
    if(is.null(RunCircleFinder)){  RunCircleFinder <- TRUE   }
#
### Functions --------------------------------------------------------------------------------------
    indexing <- function(X) { sapply(unique(X), function(y) list(y)) }
    #
    GeneExonCovByRegion <- function( GENE, CHROM, START, END, GENOME ){
        if( GENOME == 'hg19' ){ DbTable <- "hg19_gencode_v19_all_transcripts_seq" }else{ DbTable <- "hg38_gencode_v38_all_transcripts_seq" }
        dbCon = dbConnect(dbDriver("MySQL"), host='192.168.0.34', user='gcx', port=3306, password='gencurix!!', db='public_database' )
        on.exit(dbDisconnect(dbCon))

        mane_tx <- dbGetQuery(dbCon, sprintf("SELECT gene,ensembl FROM MANE_selection WHERE gene = '%s'", GENE))
        if( nrow(mane_tx) == 0 ){ 
            GeneData <- dbGetQuery(dbCon, 
                sprintf("SELECT gene_name,transcript_id,chr,exon_start,exon_end,exon_number,gene_start,gene_end,strand FROM %s WHERE gene_name = '%s'", DbTable, GENE )
            )
            MaxExonTx <- GeneData %>% group_by(transcript_id) %>% reframe(max_exon=max(exon_number))
            TX <- MaxExonTx[which.max(MaxExonTx$max_exon), "transcript_id"]
            GeneData  <- GeneData %>% filter( transcript_id == TX )
        }else{
            GeneData <- dbGetQuery(dbCon, 
                sprintf("SELECT gene_name,transcript_id,chr,exon_start,exon_end,exon_number,gene_start,gene_end,strand FROM %s WHERE transcript_id = '%s'", DbTable, mane_tx$ensembl )
            )
        }
        RegionGeneData <- GeneData %>% filter( exon_start <= END, exon_end >= START )
        if( nrow(RegionGeneData) == 0 ){
            GCOV <- "" 
        }else{
            if( all(GeneData$exon_number %in% RegionGeneData$exon_number) ){
                GCOV <- GENE
            }else{ 
                GCOV <- paste0(GENE, " (ex", min(RegionGeneData$exon_number), "-", max(RegionGeneData$exon_number),")")
            }
        }
        return(GCOV)
    }
#
### Resources
    GeneGR_List <- readRDS("/storage/home/kangsm/myDB/rds/Genes.GRanges.List.By.GTF.RDS")
    if( GenomeAssembly == "hg19" ){ GenesGR <- GeneGR_List$hg19 }else{ GenesGR <- GeneGR_List$hg38 }
    ChromIndex <- c(paste0("chr", 1:22), "chrX","chrY","chrM")
#

### CIRCLE-MAP ###----------------------------------------------------------------------------------
    if( RunCircleMap )
    {
        ## Folders 
        CircleMapDataDir <- sprintf("%s/%s/%s/ecDNA/circlemap", BaseDir, BatchID, SeqID)

        ## Read Result
        # chrom               chromosome
        # start postiion      start position (0-based)
        # end position        end position (0-based)
        # discordants         the number of discordant read pairs supporting the detected circular DNA
        # split reads         the number of split reads supporting the detected circular DNA
        # circle score        50-200 decent score, >200 high qual 
        # mean coverage       mean base coverage
        # cov sd              standard deviation of the base coverage
        # cov increase pos    > 0.33
        # cov continuity      > 0.33
        
        ResCols <- c("Chrom","StartPos","EndPos","Discordants","SupportSplitReads","CircleScore","CovMean","CovSd","CovIncreaseStart","CovIncreaseEnd","CovContinuity")
        cm_res  <- read.delim(sprintf("%s/%s.ecDNA.circlemap.result.bed", CircleMapDataDir, SeqID), header=FALSE)
        colnames(cm_res) <- ResCols
        cm_res <- cm_res %>% mutate(Size_kb=round((EndPos-StartPos+1)/1000, 3))

        ## Filtering CutOff reference : PubMed: PMID 41273082 (CircleBase V2)

        cm_res_filtered <- cm_res %>% filter(
            Chrom %in% ChromIndex,
            Discordants >= 3, SupportSplitReads >= 3, CircleScore > 50, CovIncreaseStart > 0.33 , CovIncreaseEnd > 0.33
        ) %>% mutate(
            CovMean=round(CovMean,2), CovSd=round(CovSd,2), CovIncreaseStart=round(CovIncreaseStart,2), CovIncreaseEnd=round(CovIncreaseEnd,2), CovContinuity=round(CovContinuity,3)
        )
        rownames(cm_res_filtered) <- NULL

        CircleMap_GR <- GRanges(cm_res_filtered$Chrom, IRanges(start=cm_res_filtered$StartPos, end=cm_res_filtered$EndPos))
        CircleMap_OvGene <- findOverlaps( CircleMap_GR, GenesGR )

        CircleMapGeneMatch <- CircleMap_OvGene %>% as.data.frame() %>% 
            mutate( cm_res_filtered[ .$queryHits, c("Chrom","StartPos","EndPos")] ) %>% 
            mutate( Gene = GenesGR[ .$subjectHits ]$gene )

        CircleMapGeneMatch$GeneExonCov <- apply(CircleMapGeneMatch, 1, function(gm){
            GeneExonCovByRegion( GENE=gm[6], CHROM=gm[3], START=as.numeric(gm[4]), END=as.numeric(gm[5]), GENOME=GenomeAssembly )
        } )

        CircleMapGeneMatchRev <- CircleMapGeneMatch %>% filter( GeneExonCov != "" ) %>% group_by(queryHits) %>% reframe(
            NumberOfGene = length(GeneExonCov),
            GeneExonCoverage = paste0( sort(GeneExonCov), collapse="; " )
        ) %>% as.data.frame()

        cm_res_filtered <- cm_res_filtered %>% mutate(Size_kb = round((EndPos-StartPos+1)/1000, 3) ) %>% mutate( NumberOfGene=0, GeneExonCoverage = "" )
        for( i in 1:nrow(CircleMapGeneMatchRev) ){
            cm_res_filtered[CircleMapGeneMatchRev[i, "queryHits"], "NumberOfGene"    ] <- CircleMapGeneMatchRev[i, "NumberOfGene"    ]
            cm_res_filtered[CircleMapGeneMatchRev[i, "queryHits"], "GeneExonCoverage"] <- CircleMapGeneMatchRev[i, "GeneExonCoverage"]
        }
        
        write.table( cm_res_filtered,
            sprintf("%s/%s.ecDNA.circlemap.gene.annotated.result.tsv",CircleMapDataDir, SeqID),
            quote=F, col.names=T, row.names=F, sep="\t"
        )
        write.xlsx( cm_res_filtered,
            sprintf("%s/%s.ecDNA.circlemap.gene.annotated.result.xlsx",CircleMapDataDir, SeqID),
            rowNames=F, overwrite=T
        )
    }
###

### CIRCLE_FINDER ###-------------------------------------------------------------------------------
    if( RunCircleFinder )
    {
        ## Folders 
        CicrleFinderDataDir <- sprintf("%s/%s/%s/ecDNA/circlefinder", BaseDir, BatchID, SeqID)

        ## Read Result
        cf_res <- read.delim(sprintf("%s/%s-%s.microDNA-JT.txt",CicrleFinderDataDir, SeqID, GenomeAssembly), header=FALSE)
        colnames(cf_res) <- c("Chrom","StartPos","EndPos","SplitReads")

        cf_res_filtered <- cf_res %>% mutate(Size_kb=round((EndPos-StartPos+1)/1000, 2)) %>% filter( Chrom %in% ChromIndex, SplitReads >= 5, Size_kb <= 10 )

        circleFinder_GR <- GRanges(cf_res_filtered$Chrom, IRanges(start=cf_res_filtered$StartPos, end=cf_res_filtered$EndPos))
        circleFinder_OvGene <- findOverlaps( circleFinder_GR, GenesGR )
    
        CircleFinderGeneMatch <- circleFinder_OvGene %>% as.data.frame() %>% 
            mutate( cf_res_filtered[ .$queryHits, c("Chrom","StartPos","EndPos")] ) %>% 
            mutate( Gene = GenesGR[ .$subjectHits ]$gene )

        CircleFinderGeneMatch$GeneExonCov <- apply(CircleFinderGeneMatch, 1, function(gm){
            GeneExonCovByRegion( GENE=gm[6], CHROM=gm[3], START=as.numeric(gm[4]), END=as.numeric(gm[5]), GENOME=GenomeAssembly )
        } )
        CircleFinderGeneMatchRev <- CircleFinderGeneMatch %>% filter( GeneExonCov != "" ) %>% group_by(queryHits) %>% reframe(
            NumberOfGene = length(GeneExonCov),
            GeneExonCoverage = paste0( sort(GeneExonCov), collapse="; " )
        ) %>% as.data.frame()

        cf_res_filtered <- cf_res_filtered %>% mutate(Size_kb = round((EndPos-StartPos+1)/1000, 1) ) %>% mutate( NumberOfGene=0, GeneExonCoverage = "" )
        for( i in 1:nrow(CircleFinderGeneMatchRev) ){
            cf_res_filtered[CircleFinderGeneMatchRev[i, "queryHits"], "NumberOfGene"    ] <- CircleFinderGeneMatchRev[i, "NumberOfGene"    ]
            cf_res_filtered[CircleFinderGeneMatchRev[i, "queryHits"], "GeneExonCoverage"] <- CircleFinderGeneMatchRev[i, "GeneExonCoverage"]
        }

        write.table( cf_res_filtered,
            sprintf("%s/%s.ecDNA.circlefinder.gene.annotated.result.tsv",CicrleFinderDataDir, SeqID),
            quote=F, col.names=T, row.names=F, sep="\t"
        )
        write.xlsx( cf_res_filtered,
            sprintf("%s/%s.ecDNA.circlefinder.gene.annotated.result.xlsx",CicrleFinderDataDir, SeqID),
            rowNames=F, overwrite=T
        )
    }
###


