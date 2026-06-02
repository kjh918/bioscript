

#---| PACKAGES |-----------------------------------------------------------------------------------#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("Cairo"))
grDevices::X11.options(type='cairo')
options(device='x11')
options(stringsAsFactors=FALSE) 
options(bitmapType='cairo')
#--------------------------------------------------------------------------------------------------#

#---| Set Defualt Values |-------------------------------------------------------------------------#
    BASE_DIR          = "/data/wes"
    WES_LIB_KIT       = "twist.exome.2.0"
    TUMOR_SEQ_ID      = ""
    THREADS           = 10
    NORMAL_TYPE       = "TS"
    INPUT_SEX         = "XX"
    TARGETED_SEQ      = TRUE
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--baseDir"), action="store", default=NA, type="character", help="BASE_DIR"),
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="SEQ_FOLDER"),   
        make_option(c("--tumorSeqID"), action="store", default="", type="character", help="TUMOR_SEQ_ID"),
        make_option(c("--wesLibKit"), action="store", default=NA, type="character", help="BED_FILE"),
        make_option(c("--threads"), action="store", default=15, type="double", help="THREADS"),
        make_option(c("--nType"), action="store", default="TS", type="character", help="MATCHED NORMAL TYPE. 'TS' or 'ORG'")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR          = ARGS$baseDir
    SEQ_FOLDER        = ARGS$seqFolder
    TUMOR_SEQ_ID      = ARGS$tumorSeqID
    WES_LIB_KIT       = ARGS$wesLibKit
    THREADS           = as.numeric(ARGS$threads)
    NORMAL_TYPE       = ARGS$nType
#--------------------------------------------------------------------------------------------------#

#---| MODIFIED FUNCTIONS |-------------------------------------------------------------------------#
    ascat.prepareHTS_rev = function (tumourseqfile, normalseqfile = NA, tumourname, normalname = NA,
        allelecounter_exe, alleles.prefix, loci.prefix, gender, genomeVersion, 
        nthreads = 1, tumourLogR_file = NA, tumourBAF_file = NA,
        normalLogR_file = NA, normalBAF_file = NA, minCounts = 10,
        BED_file = NA, probloci_file = NA, chrom_names = c(1:22), min_base_qual = 20, min_map_qual = 35, additional_allelecounter_flags = NA,
        skip_allele_counting_tumour = FALSE, skip_allele_counting_normal = FALSE)
    {
        suppressPackageStartupMessages(library("foreach"))
        suppressPackageStartupMessages(library("doParallel"))
        suppressPackageStartupMessages(library("parallel"))
        registerDoParallel(cores = nthreads)
        if (is.na(tumourLogR_file))
            tumourLogR_file = paste0(tumourname, "_tumourLogR.txt")
        if (is.na(tumourBAF_file))
            tumourBAF_file = paste0(tumourname, "_tumourBAF.txt")
        if (is.na(normalLogR_file))
            normalLogR_file = paste0(tumourname, "_normalLogR.txt")
        if (is.na(normalBAF_file))
            normalBAF_file = paste0(tumourname, "_normalBAF.txt")
        if (is.na(normalseqfile)) {
            skip_allele_counting_normal = TRUE
            tumour_only_mode = TRUE
            print("No normal bam file provided, running in Tumour-only mode")
        }else{
            tumour_only_mode = FALSE
        }
        if (!skip_allele_counting_tumour) {
            foreach(CHR = chrom_names) %dopar% {
                CHR = get("CHR")
                ascat.getAlleleCounts(seq.file = tumourseqfile, output.file = paste0(tumourname,
                    "_alleleFrequencies_chr", CHR, ".txt"), loci.file = paste0(loci.prefix,
                    CHR, ".txt"), min.base.qual = min_base_qual,
                    min.map.qual = min_map_qual, allelecounter.exe = allelecounter_exe,
                    additional_allelecounter_flags = additional_allelecounter_flags)
            }
        }
        if (!skip_allele_counting_normal) {
            foreach(CHR = chrom_names) %dopar% {
                CHR = get("CHR")
                ascat.getAlleleCounts(seq.file = normalseqfile, output.file = paste0(normalname,
                    "_alleleFrequencies_chr", CHR, ".txt"), loci.file = paste0(loci.prefix,
                    CHR, ".txt"), min.base.qual = min_base_qual,
                    min.map.qual = min_map_qual, allelecounter.exe = allelecounter_exe,
                    additional_allelecounter_flags = additional_allelecounter_flags)
            }
        }
        ascat.getBAFsAndLogRs_v2(samplename = tumourname, tumourAlleleCountsFile.prefix = paste0(tumourname,
            "_alleleFrequencies_chr"), normalAlleleCountsFile.prefix = paste0(normalname,
            "_alleleFrequencies_chr"), tumourLogR_file = tumourLogR_file,
            tumourBAF_file = tumourBAF_file, normalLogR_file = normalLogR_file,
            normalBAF_file = normalBAF_file, alleles.prefix = alleles.prefix,
            gender = gender, genomeVersion = genomeVersion, chrom_names = chrom_names,
            minCounts = minCounts, BED_file = BED_file, probloci_file = probloci_file,
            tumour_only_mode = tumour_only_mode)
        ascat.synchroniseFiles_v2(samplename = tumourname, tumourLogR_file = tumourLogR_file,
            tumourBAF_file = tumourBAF_file, normalLogR_file = normalLogR_file,
            normalBAF_file = normalBAF_file)
    }
    #---------------------------------------------------------------------------
    ascat.synchroniseFiles_v2 = function (samplename, tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file)
    {
        FILES = lapply(c(tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file), function(x) {
            if (!file.exists(x))
            {
                return(NA)
            }else{
                tmp = data.frame(fread(x, sep = "\t", showProgress = FALSE, header = TRUE, na.strings = c("-Inf", "Inf", "NA","NaN", "", "-")), 
                    row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
                colnames(tmp) = c("Chromosome", "Position", samplename)
                tmp = tmp[!is.na(tmp[, 3]), ]
                return(tmp)
            }
        })
        names(FILES) = c(tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file)
        FILES = FILES[!is.na(FILES)]
        IDs   = Reduce(intersect, lapply(FILES, rownames))
        FILES = lapply(FILES, function(x) x[rownames(x) %in% IDs, ])
        rm(IDs)
        for( i in 2:length(FILES) )
        {
            FILES[[i]] = FILES[[i]][rownames(FILES[[1]]), ]
        }
        for (i in 1:length(FILES)) {
            write.table(FILES[[i]], file = names(FILES)[i], sep = "\t",
                quote = FALSE, row.names = TRUE, col.names = NA)
        }
        rm(i)
    }
    ascat.getBAFsAndLogRs_v2 = function (samplename, tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix,
    tumourLogR_file, tumourBAF_file, normalLogR_file, normalBAF_file,
    alleles.prefix, gender, genomeVersion, chrom_names = c(1:22), minCounts = 10, BED_file = NA, probloci_file = NA,
    tumour_only_mode = FALSE)
    {
        #set.seed(seed)
        stopifnot(gender %in% c("XX", "XY"))
        stopifnot(genomeVersion %in% c("hg19", "hg38", "CHM13"))
        tumour_input_data = readAlleleCountFiles(tumourAlleleCountsFile.prefix,
            ".txt", chrom_names, 1)
        if (tumour_only_mode) {
            normal_input_data = tumour_input_data
        }else{
            normal_input_data = readAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", chrom_names, minCounts)
        }
        allele_data  = readAllelesFiles(alleles.prefix, ".txt", chrom_names)
        matched_data = Reduce(intersect, list(rownames(tumour_input_data), rownames(normal_input_data), rownames(allele_data)))
        tumour_input_data = tumour_input_data[rownames(tumour_input_data) %in% matched_data, ]
        normal_input_data = normal_input_data[rownames(normal_input_data) %in% matched_data, ]
        allele_data       = allele_data[rownames(allele_data) %in% matched_data, ]
        rm(matched_data)
        if (!is.na(probloci_file)) {
            stopifnot(file.exists(probloci_file) && file.info(probloci_file)$size > 0)
            probloci = data.frame(fread(probloci_file, sep = "\t", showProgress = FALSE, header = TRUE), stringsAsFactors = FALSE)
            probloci = paste0(gsub("^chr", "", probloci[, 1]), "_", probloci[, 2])
            probloci = which(rownames(tumour_input_data) %in% probloci)
            if (length(probloci) > 0) 
            {
                tumour_input_data = tumour_input_data[-probloci, ]
                normal_input_data = normal_input_data[-probloci, ]
                allele_data = allele_data[-probloci, ]
            }else{
                warning("The probloci did not remove any SNPs, it might be worth checking the data.")
            }
            rm(probloci)
        }
        stopifnot(isTRUE(all.equal(allele_data[, 1], tumour_input_data[, 1])) && isTRUE(all.equal(allele_data[, 1], normal_input_data[, 1])))
        stopifnot(isTRUE(all.equal(allele_data[, 2], tumour_input_data[, 2])) && isTRUE(all.equal(allele_data[, 2], normal_input_data[, 2])))
        tumour_input_data = tumour_input_data[, 3:6]
        normal_input_data = normal_input_data[, 3:6]
        if (!is.na(BED_file)) 
        {
            stopifnot(file.exists(BED_file) && file.info(BED_file)$size > 0)
            BED = read.table(BED_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1:3]
            colnames(BED) = c("chr", "start", "end")
            BED$chr       = gsub("^chr", "", BED$chr)
            BED$start     = BED$start + 1
            BED           = BED[BED$chr %in% chrom_names, ]
            if (nrow(BED) == 0) { stop("Major issue with BED file, please double-check its content") }
            suppressPackageStartupMessages(library("GenomicRanges"))
            suppressPackageStartupMessages(library("IRanges"))
            overlaps = findOverlaps(
                GRanges(seqnames = BED$chr, ranges = IRanges(start = BED$start, end = BED$end)), 
                GRanges(seqnames = allele_data$chromosome, ranges = IRanges(start = allele_data$position, end = allele_data$position))
            )
            if (length(overlaps) > 0) {
                tumour_input_data = tumour_input_data[unique(overlaps@to), ]
                normal_input_data = normal_input_data[unique(overlaps@to), ]
                allele_data = allele_data[unique(overlaps@to), ]
            }else{
                print(head(allele_data))
                print(head(BED))
                stop("The overlap between the BED file and loci is empty. Data must be checked!")
            }
            rm(BED, overlaps)
        }
        len = nrow(allele_data)
        tumour_input_data$REF = tumour_input_data[cbind(1:len, allele_data[, 3])]
        tumour_input_data$ALT = tumour_input_data[cbind(1:len, allele_data[, 4])]
        normal_input_data$REF = normal_input_data[cbind(1:len, allele_data[, 3])]
        normal_input_data$ALT = normal_input_data[cbind(1:len, allele_data[, 4])]
        if (tumour_only_mode) {
            TO_KEEP = which(tumour_input_data$REF + tumour_input_data$ALT >= 1)
        }else{
            TO_KEEP = which(tumour_input_data$REF + tumour_input_data$ALT >= 1 & normal_input_data$REF + normal_input_data$ALT >= minCounts)
        }
        stopifnot(length(TO_KEEP) > 0)
        allele_data = allele_data[TO_KEEP, ]
        tumour_input_data = tumour_input_data[TO_KEEP, ]
        normal_input_data = normal_input_data[TO_KEEP, ]
        rm(TO_KEEP)
        len         = nrow(allele_data)
        mutCount1   = tumour_input_data$REF
        mutCount2   = tumour_input_data$ALT
        totalTumour = mutCount1 + mutCount2
        normCount1  = normal_input_data$REF
        normCount2  = normal_input_data$ALT
        totalNormal = normCount1 + normCount2
        rm(tumour_input_data, normal_input_data)
        normalBAF  = vector(length = len, mode = "numeric")
        tumourBAF  = vector(length = len, mode = "numeric")
        normalLogR = vector(length = len, mode = "numeric")
        tumourLogR = vector(length = len, mode = "numeric")
        normalBAF_unmirrored = normCount2/totalNormal
        tumourBAF_unmirrored = mutCount2/totalTumour

        germline.BAF_unmirrored = data.frame(Chromosome = allele_data$chromosome, Position = allele_data$position, baf = normalBAF_unmirrored,
            ID = rownames(allele_data), row.names = 4, stringsAsFactors = FALSE)
        tumor.BAF_unmirrored = data.frame(Chromosome = allele_data$chromosome, Position = allele_data$position, baf = tumourBAF_unmirrored,
            ID = rownames(allele_data), row.names = 4, stringsAsFactors = FALSE)
        colnames(tumor.BAF_unmirrored)[3]    = samplename
        colnames(germline.BAF_unmirrored)[3] = samplename
        write.table(tumor.BAF_unmirrored, file = gsub("\\.txt$", "_rawBAF.txt", tumourBAF_file), row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
        if (!tumour_only_mode) 
        {
            write.table(germline.BAF_unmirrored, file = gsub("\\.txt$", "_rawBAF.txt", normalBAF_file), row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
        }
        rm(normalBAF_unmirrored, tumourBAF_unmirrored, germline.BAF_unmirrored, tumor.BAF_unmirrored)
        selector = round(runif(len))
        normalBAF[which(selector == 0)] = normCount1[which(selector == 0)]/totalNormal[which(selector == 0)]
        normalBAF[which(selector == 1)] = normCount2[which(selector == 1)]/totalNormal[which(selector == 1)]
        tumourBAF[which(selector == 0)] = mutCount1[which(selector == 0)]/totalTumour[which(selector == 0)]
        tumourBAF[which(selector == 1)] = mutCount2[which(selector == 1)]/totalTumour[which(selector == 1)]
        rm(selector)
        if (tumour_only_mode) 
        {
            tumourLogR = log2(totalTumour/median(totalTumour, na.rm = TRUE))
        }else{
            tumourLogR = totalTumour/totalNormal
            tumourLogR = log2(tumourLogR/mean(tumourLogR, na.rm = TRUE))
            if (gender == "XY") {
                if (genomeVersion == "hg19") { nonPAR = c(2699521, 154931043) 
                }else if (genomeVersion == "hg38") { nonPAR = c(2781480, 155701382)
                }else if (genomeVersion == "CHM13") { nonPAR = c(2394411, 153925834)
                }
                nonPAR = which(allele_data$chromosome %in% c("X","chrX") & allele_data$position >= nonPAR[1] & allele_data$position <= nonPAR[2])
                tumourLogR[nonPAR] = tumourLogR[nonPAR] - 1
            }
        }
        tumor.LogR = data.frame(Chromosome = allele_data$chromosome, Position = allele_data$position, logr = tumourLogR, ID = rownames(allele_data),
            row.names = 4, stringsAsFactors = FALSE)
        tumor.BAF = data.frame(Chromosome = allele_data$chromosome, Position = allele_data$position, baf = tumourBAF, ID = rownames(allele_data),
            row.names = 4, stringsAsFactors = FALSE)
        germline.LogR = data.frame(Chromosome = allele_data$chromosome, Position = allele_data$position, logr = normalLogR, ID = rownames(allele_data),
            row.names = 4, stringsAsFactors = FALSE)
        germline.BAF = data.frame(Chromosome = allele_data$chromosome, Position = allele_data$position, baf = normalBAF, ID = rownames(allele_data),
            row.names = 4, stringsAsFactors = FALSE)
        colnames(tumor.LogR)[3]    = samplename
        colnames(tumor.BAF)[3]     = samplename
        colnames(germline.LogR)[3] = samplename
        colnames(germline.BAF)[3]  = samplename
        write.table(tumor.LogR, file = tumourLogR_file, row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
        write.table(tumor.BAF, file  = tumourBAF_file, row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
        if (!tumour_only_mode) 
        {
            write.table(germline.LogR, file = normalLogR_file, row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
            write.table(germline.BAF, file = normalBAF_file, row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
        }
    }
    readAlleleCountFiles = function(prefix, suffix, chrom_names, minCounts, keep_chr_string=FALSE) 
    {
        files=paste0(prefix, chrom_names, suffix)
        files=files[sapply(files, function(x) file.exists(x) && file.info(x)$size>0)]
        stopifnot(length(files)>0)
        data=do.call(rbind, lapply(files, function(x) {
            tmp=data.frame(fread(x, sep="\t", showProgress=FALSE, header=TRUE), stringsAsFactors=FALSE)
            tmp=tmp[tmp[, 7]>=minCounts, ]
            if (nrow(tmp)>0) {
            if (!keep_chr_string) tmp[, 1]=gsub("^chr", "", tmp[, 1])
            rownames(tmp)=paste0(tmp[, 1], "_", tmp[, 2])
            }
            return(tmp)
        }))
        stopifnot(nrow(data)>0)
        return(data)
    }
    readAllelesFiles=function(prefix, suffix, chrom_names, add_chr_string=FALSE) 
    {
        files=paste0(prefix, chrom_names, suffix)
        files=files[sapply(files, function(x) file.exists(x) && file.info(x)$size>0)]
        stopifnot(length(files)>0)
        data=do.call(rbind, lapply(files, function(x) {
            tmp=data.frame(fread(x, sep="\t", showProgress=FALSE, header=TRUE))
            tmp=tmp[!is.na(tmp[, 2] & !is.na(tmp[, 3])), ]
            tmp=tmp[!duplicated(tmp[, 1]), ]
            tmp$chromosome=gsub(paste0(prefix, "(", paste(chrom_names, collapse="|"), ")", suffix), "\\1", x)
            if (add_chr_string) tmp$chromosome=paste0("chr", tmp$chromosome)
            tmp=tmp[, c(4, 1:3)]
            rownames(tmp)=paste0(tmp[, 1], "_", tmp[, 2])
            return(tmp)
        }))
        stopifnot(nrow(data)>0)
        return(data)
    }
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR = sprintf("%s/%s", BASE_DIR, SEQ_FOLDER)    
    LOG_DIR  = sprintf("%s/%s/log", DATA_DIR, TUMOR_SEQ_ID)
    #--------------------------------------------------------------------------#
    LOG_FILE = sprintf("%s/%s.ASCAT.Tonly.%s.log", LOG_DIR, TUMOR_SEQ_ID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
    #--------------------------------------------------------------------------#
    ASCAT_BASE_DIR = sprintf("%s/%s/cnv/ascat_tonly", DATA_DIR, TUMOR_SEQ_ID)
    if( !dir.exists(ASCAT_BASE_DIR) ){ system(sprintf("mkdir -p %s", ASCAT_BASE_DIR)) }
    system(sprintf("rm -r %s/*", ASCAT_BASE_DIR))
    #--------------------------------------------------------------------------#
    ASCAT_AF_DIR  = sprintf("%s/alleleFreqs", ASCAT_BASE_DIR)
    if( !dir.exists(ASCAT_AF_DIR) ){ system(sprintf("mkdir -p %s", ASCAT_AF_DIR)) }
#--------------------------------------------------------------------------------------------------#

#---| ASCAT REFERENCES |---------------------------------------------------------------------------#
    GENOME_FASTA      = "/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta"
    REF_FLAT          = "/storage/references_and_index/cnv/ascat/hg19.refFlat.txt"
    #ASCAT_REF_ALLELE = "/storage/references_and_index/cnv/ascat/G1000_allelesAll_hg19/G1000_alleles_hg19_chr"
    #ASCAT_REF_LOCI   = "/storage/references_and_index/cnv/ascat/G1000_lociAll_hg19/G1000_loci_hg19_chr"
    ASCAT_REF_GC      = "/storage/references_and_index/cnv/ascat/GC_G1000_hg19.txt"
    ASCAT_REF_RT      = "/storage/references_and_index/cnv/ascat/RT_G1000_hg19.txt"
    ASCAT_TS_ALLELE   = sprintf("/storage/references_and_index/cnv/ascat/targeted_seq/%s/alleleData/Cleaned/alleleData_chr", WES_LIB_KIT)
    ASCAT_TS_LOCI     = sprintf("/storage/references_and_index/cnv/ascat/targeted_seq/%s/alleleData/Cleaned/loci_chr", WES_LIB_KIT)
    #--------------------------------------------------------------------------#
    BED_FILE = sprintf("/storage/references_and_index/hg19/bed/%s/hg19_%s.target.bed", WES_LIB_KIT, WES_LIB_KIT)
    #GENE_BED = sprintf("/storage/references_and_index/%s/bed/%s/%s_%s.probe.cnvkit.gene.bed", "hg19", WES_LIB_KIT, "hg19", WES_LIB_KIT)
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES ASCAT CNV Tumor-Only Mode log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE       : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER     : %s", SEQ_FOLDER), LOG_FILE, append=T) 
    write(sprintf(" TUMOR SEQ ID   : %s", TUMOR_SEQ_ID), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

#==================================================================================================# 
    write(sprintf("%s | [ Preapre HTS-DATA ] START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message(">>>>> [ Preapre HTS-DATA ] START.")
#==================================================================================================# 

#---| RESULT TAG & FOLDER |------------------------------------------------------------------------#
    if( NORMAL_TYPE == "ORG" ){ RES_TAG = "ORG.TumorOnly" }else{ RES_TAG = "TumorOnly" }
    ASCAT_RES_DIR = sprintf("%s/%s.%s.result", ASCAT_BASE_DIR, TUMOR_SEQ_ID, RES_TAG)
    if( ! dir.exists(ASCAT_RES_DIR) ){ system(sprintf("mkdir -p %s", ASCAT_RES_DIR)) }
#--------------------------------------------------------------------------------------------------#

#---| MOVE TO ALLELE-FREQ DIR |-------------------------------------------------------------------------#
    setwd(ASCAT_AF_DIR)
#--------------------------------------------------------------------------------------------------#

#---| SUPERFREQ RESULT FOR SEX |-------------------------------------------------------------------#
    # if( is.null(INPUT_SEX) )
    # {
    #     SUPERFREQ_RES = sprintf("%s/%s/cnv/superfreq/%s.%s.superFreq.summary.tsv", DATA_DIR, TUMOR_SEQ_ID, TUMOR_SEQ_ID, SUPERFREQ_SET_ID)   
    #     if( !file.exists(SUPERFREQ_RES) )
    #     { 
    #         message("RUN superFreq First to determine sex of input OR assign sample's sex manually")
    #         quit(save="no", status=0)
    #     }
    #     superFreqRes = read.delim(SUPERFREQ_RES)
    #     superFreqRes = superFreqRes[which(superFreqRes$NAME == TUMOR_SEQ_ID ), ]
    #     INPUT_SEX    = ifelse( superFreqRes$SEX == "male", "XY", "XX" )
    # }
#--------------------------------------------------------------------------------------------------#

#---| PREPARE HTS DATA |---------------------------------------------------------------------------#
    source_normal_logR_file = sprintf("/storage/references_and_index/cnv/ascat/targeted_seq/%s/Normal_data/Normal_LogR.txt", WES_LIB_KIT)
    source_normal_baf_file  = sprintf("/storage/references_and_index/cnv/ascat/targeted_seq/%s/Normal_data/Normal_BAF.txt", WES_LIB_KIT)
    NORMAL_LOGR_FILE = sprintf("%s/Normal_LogR.txt", ASCAT_RES_DIR)
    NORMAL_BAF_FILE  = sprintf("%s/Normal_BAF.txt", ASCAT_RES_DIR)
    TUMOR_LOGR_FILE  = sprintf("%s/%s.%s.Tumor.LogR.txt",    ASCAT_RES_DIR, TUMOR_SEQ_ID,  RES_TAG )
    TUMOR_BAF_FILE   = sprintf("%s/%s.%s.Tumor.BAF.txt",     ASCAT_RES_DIR, TUMOR_SEQ_ID,  RES_TAG )
    system(sprintf("cp %s %s", source_normal_logR_file, NORMAL_LOGR_FILE))
    system(sprintf("cp %s %s", source_normal_baf_file, NORMAL_BAF_FILE))
    #----------------------------------------------------------------------#
    ascat.prepareHTS_rev(
        tumourseqfile     = sprintf( "%s/%s/bam/%s.analysisReady.bam", DATA_DIR, TUMOR_SEQ_ID, TUMOR_SEQ_ID ),
        tumourname        = TUMOR_SEQ_ID,
        BED_file          = BED_FILE,
        allelecounter_exe = "/storage/apps/alleleCount-4.2.1/bin/alleleCounter",
        alleles.prefix    = ASCAT_TS_ALLELE,
        loci.prefix       = ASCAT_TS_LOCI,
        genomeVersion     = "hg19",
        nthreads          = THREADS,
        gender            = INPUT_SEX,
        chrom_names       = C(1:22),
        tumourLogR_file   = TUMOR_LOGR_FILE,
        tumourBAF_file    = TUMOR_BAF_FILE,
        normalLogR_file   = NORMAL_LOGR_FILE,
        normalBAF_file    = NORMAL_BAF_FILE,
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
    ascat.bc = ascat.loadData(
        Tumor_LogR_file    = TUMOR_LOGR_FILE, 
        Tumor_BAF_file     = TUMOR_BAF_FILE, 
        Germline_LogR_file = NORMAL_LOGR_FILE, 
        Germline_BAF_file  = NORMAL_BAF_FILE, 
        genomeVersion      = "hg19", 
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
    ascat.bc = ascat.correctLogR( 
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
    ascat.bc = ascat.aspcf( ascat.bc, penalty = 70 )  # = create BAF.PCFed.txt & LogR.PCFed.txt    
    #--------------------------------------------------------------------------#
    # rename_baf_cmd  = sprintf("mv %s/%s.BAF.PCFed.txt %s/%s.%s.BAF.PCFed.txt",   ASCAT_RES_DIR, TUMOR_SEQ_ID, ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
    # rename_logr_cmd = sprintf("mv %s/%s.LogR.PCFed.txt %s/%s.%s.LogR.PCFed.txt", ASCAT_RES_DIR, TUMOR_SEQ_ID, ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
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
    # rename_aspcf_cmd  = sprintf("mv %s/%s.ASPCF.png %s/%s.%s.ASPCF.png", ASCAT_RES_DIR, TUMOR_SEQ_ID, ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
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
    ascat.qc = ascat.metrics(
        ascat.bc,
        ascat.output
    )
#--------------------------------------------------------------------------------------------------#

#---| SAVE RESULTS |-------------------------------------------------------------------------------#
    ASCAT_RESULT_RDATA = sprintf("%s/%s.%s.ASCAT.results.Rdata", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG)
    #--------------------------------------------------------------------------#
    save( ascat.bc, ascat.output, ascat.qc, file = ASCAT_RESULT_RDATA )
    #--------------------------------------------------------------------------#
    write.table(
        ascat.qc, 
        sprintf("%s/%s.%s.ASCAT.qc.result.tsv", ASCAT_RES_DIR, TUMOR_SEQ_ID, RES_TAG ), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    #--------------------------------------------------------------------------#
    qc.simple            = data.frame(ID=rownames(ascat.qc), ascat.qc[,c("purity", "sex", "ploidy")], ploidy_int="") 
    qc.simple$sex        = sapply( qc.simple$sex, function(sx) ifelse( sx == "XX", "female", "male" ) )
    qc.simple$ploidy_int = round( qc.simple$ploidy )
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
        qc.simple, sprintf("%s/%s.%s.ASCAT.summary.tsv", ASCAT_BASE_DIR, TUMOR_SEQ_ID, RES_TAG ), 
        quote=F, col.names=T, row.names=F, sep="\t"
    )
    write.table(
        qc.simple, sprintf("%s/%s.%s.ASCAT.summary.txt", ASCAT_BASE_DIR, TUMOR_SEQ_ID, RES_TAG ), 
        quote=F, col.names=F, row.names=F, sep=" "
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | [ CREATE QC and SAVE ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("<<<<< [ CREATE QC and SAVE ] FINISHED.")
    write(sprintf("%s | [ ASCAT CNV ANALYSIS ] FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
    write("   ", LOG_FILE, append=T)
    message(">>>>> [ ASCAT CNV ANALYSIS ] FINISHED.")
#==================================================================================================# 

