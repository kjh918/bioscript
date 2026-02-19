
# script version : v0.1

#---| PACKAGES |-----------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("lubridate"))
#--------------------------------------------------------------------------------------------------# 

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--seqFolder")           , action="store", default=NA,      type="character", help="$SEQ_FOLDER" ),
        make_option(c("--seqID")               , action="store", default=NA,      type="character", help="$SEQ_ID" ),
        make_option(c("--baseDir")             , action="store", default="/data/wes", type="character", help="Base Work Directory" ),
        make_option(c("--cov")                 , action="store", default=NA,      type="character", help="Alfred QC Coverage File, required only when useReferenceMethod = FALSE" ),
        make_option(c("--bed")                 , action="store", default=NA,      type="character", help="Target BED file, required only when useReferenceMethod = FALSE" ),
        make_option(c("--nType")               , action="store", default="NT",      type="character", help="type of normal sample" ),
        make_option(c("--DP")                  , action="store", default=25,      type="integer",   help="Variant Cut-Off: Tumor Read Depth, required only when useReferenceMethod = FALSE" ), 
        make_option(c("--tumorALT")            , action="store", default=3,       type="integer",   help="Variant Cut-Off: Minimum Tumor-Allele Count, required only when useReferenceMethod = FALSE" ), 
        make_option(c("--vafSNV")              , action="store", default=0.02,    type="double",    help="Variant Cut-Off: SNV VAF, required only when useReferenceMethod = FALSE" ),
        make_option(c("--vafINDEL")            , action="store", default=0.05,    type="double",    help="Variant Cut-Off: INDEL VAF, required only when useReferenceMethod = FALSE" ),
        make_option(c("--popAF")               , action="store", default=0.01,    type="double",    help="Variant Cut-Off: Population AF, required only when useReferenceMethod = FALSE" ),
        make_option(c("--useReferenceMethod")  , action="store", default=TRUE,    type="logical",   help="use reference method or not" ), 
        make_option(c("--useReferenceSize")    , action="store", default=TRUE,    type="logical",   help="use reference exonic size or not" ), 
        make_option(c("--excludeSilent")       , action="store", default=TRUE,    type="logical",   help="exclude silent variants or not, required only when useReferenceMethod = FALSE" ), 
        make_option(c("--excludeClinvar")      , action="store", default=FALSE,   type="logical",   help="exclude clinvar's pathogenic variants or not, required only when useReferenceMethod = FALSE" ), 
        make_option(c("--excludeAlphaMissense"), action="store", default=FALSE,   type="logical",   help="exclude alpha-missense's pathogenic variants or not, required only when useReferenceMethod = FALSE" ), 
        make_option(c("--variantCallMode")     , action="store", default="tonly", type="character", help="variant calling mode" ),
        make_option(c("--importDB")            , action="store", default=FALSE,   type="logical",   help="Import Result into DB or Not" )
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    seqFolder                         <- ARGS$seqFolder
    seqID                             <- ARGS$seqID
    BASE_DIR                          <- ARGS$baseDir 
    COVERAGE_FILE                     <- ARGS$cov
    TARGET_BED_FILE                   <- ARGS$bed
    NORMAL_TYPE                       <- ARGS$nType
    gcx_DP                            <- ARGS$DP
    gcx_TUMOR_ALT                     <- ARGS$tumorALT
    gcx_SNV_VAF                       <- ARGS$vafSNV
    gcx_INDEL_VAF                     <- ARGS$vafINDEL
    gcx_POPAF                         <- ARGS$popAF
    useReferenceFilter                <- ARGS$useReferenceMethod
    useReferenceSize                  <- ARGS$useReferenceSize
    gcxExcludeSilent                  <- ARGS$excludeSilent
    gcxExcludeClinvarPathogenic       <- ARGS$excludeClinvar
    gcxExcludeAlphaMissensePathogenic <- ARGS$excludeAlphaMissense
    variantCallMode                   <- ARGS$variantCallMode
    importDB                          <- ARGS$importDB
#--------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTIONS |-----------------------------------------------------------------------#
    source("/data/wes/params/ruo_wes_db.R")
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    TMB_RES_DIR <- sprintf("%s/%s/%s/tmb", BASE_DIR, seqFolder, seqID)
    LOG_DIR     <- sprintf("%s/%s/%s/log", BASE_DIR, seqFolder, seqID)
    if( ! dir.exists(TMB_RES_DIR) ) { system(sprintf("mkdir -p %s", TMB_RES_DIR)) }
    if( ! dir.exists(LOG_DIR)     ) { system(sprintf("mkdir -p %s", LOG_DIR))     }
#--------------------------------------------------------------------------------------------------#

#---| LOG FILE |-----------------------------------------------------------------------------------#
    LOG_FILE <- sprintf("%s/%s.TMB.calculation.%s.log", LOG_DIR, seqID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
#--------------------------------------------------------------------------------------------------#

#---| RUN METHOD |---------------------------------------------------------------------------------#
    RUN_METHOD <- ifelse( useReferenceFilter, "TMB Harmonization Project", "Gencurix")
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES TMB Calculation log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE   : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")),      LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER : %s", seqFolder),  LOG_FILE, append=T) 
    write(sprintf(" SEQ ID     : %s", seqID),      LOG_FILE, append=T) 
    write(sprintf(" METHOD     : %s", RUN_METHOD), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

#==================================================================================================# 
    write(sprintf("%s | TMB CALCULATION START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("TMB CALCULATION START.")
#==================================================================================================# 

#---| EXON REGION DICISION |-----------------------------------------------------------------------#
    if( useReferenceSize )
    {
        #---| USE TMB HARMONIZATION PROJECT EXON LENGTH |------------------------------------------#
        LIB_LENGTH <- 32102474  # Friends of Cancer Research TMB Harmonization Poject Phase I reference 
        ##
        #==========================================================================================# 
            write(sprintf("%s | USE FIXED Exonic Region Size", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("USE FIXED Exonic Region Size.")
        #==========================================================================================# 
    }else{
        #==========================================================================================#
            write(sprintf("%s | CALCULATE Exonic Region Size START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("Calculating Exonic Region Size START.")
        #==========================================================================================#

        #---| MANUALLY EXON LENGTH CALCULATION |---------------------------------------------------#
        ##
        COV <- read.delim(COVERAGE_FILE) %>% mutate( key = paste(Chr,Start,End,sep="_"))
        #----------------------------------------------------------------------#
        BED           <- read.delim(TARGET_BED_FILE, header=F)
        colnames(BED) <- c("chr","start","end","gene_raw")
        BED$gene      <- sapply(BED$gene_raw, function(gr) unlist(strsplit(gr, ";"))[1])
        BED           <- BED[,c("chr","start","end","gene")] %>% mutate( key = paste(chr,start,end,sep="_"))
        #----------------------------------------------------------------------#
        covRes <- COV %>% 
            left_join(., BED[,c("key","gene")], by="key") %>%
            filter( gene != "", AvgCov >= gcx_DP ) %>%
            group_by(Chr, Start ) %>%
            summarise( End = max(End), AvgCov = mean(AvgCov) ) %>%
            mutate( LibLength =  End - Start )
        #----------------------------------------------------------------------#
        LIB_LENGTH <- sum(covRes$LibLength) # Exon Region Length (bp) coverage >= cutoff_DP
        ##
        #==========================================================================================#
            write(sprintf("%s | CALCULATE Exonic Region Size FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("Calculating Exonic Region Size FINISHED.")
        #==========================================================================================#
    }
#--------------------------------------------------------------------------------------------------#

#---| VARIANT DATA LOADING |-----------------------------------------------------------------------#
    if( variantCallMode == "tonly" )
    {
        MAF_TAG = "mutect2.bias.filtered"
    }else if( variantCallMode == "tonlygermlinefiltered" ){  MAF_TAG = "mutect2.bias.filtered.germline.filtered"
    }else{
        MAF_TAG = ifelse( NORMAL_TYPE == "ORG", "mutect2.ORG.NT.bias.filtered", "mutect2.NT.bias.filtered" )
    }
    #--------------------------------------------------------------------------#    
    MAF_FILE <- sprintf("%s/%s/%s/vcf/%s.%s.somatic.variants.only.maf", BASE_DIR, seqFolder, seqID, seqID, MAF_TAG)
    if( ! file.exists(MAF_FILE) )
    {
        message("No MAF FILE for TMB Calculation. Analysis Stopped.")
        write(sprintf("%s | No MAF FILE for TMB Calculation. Analysis Stopped.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        quit( save = "no", status = 0 )
    }else{
        MF <- read.delim(MAF_FILE) %>% mutate( VAF = round(TUMOR_ALT_COUNT/(TUMOR_REF_COUNT+TUMOR_ALT_COUNT), 4) )
    }
    ##
#--------------------------------------------------------------------------------------------------#

    #==============================================================================================# 
        write(sprintf("%s | VARIANT FILTERING FOR TMB START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        write(sprintf("%s | Variant Filtering Method : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S"), RUN_METHOD), LOG_FILE, append=T)
        message("VARIANT FILTERING FOR TMB START.")
    #==============================================================================================# 

#---| VARIANT FILTER |-----------------------------------------------------------------------------#
    ##
    if( useReferenceFilter )
    {
        #---| USE TMB HARMONIZATION PROJECT METHOD |-----------------------------------------------#
        IncludeVariantClass <- c(
            "Missense_Mutation", "Nonsense_Mutation",
            "In_Frame_Del"     , "In_Frame_Ins",
            "Frame_Shift_Del"  , "Frame_Shift_Ins"
        )
        MF_tmb <- MF %>% 
            filter( Variant_Classification %in% IncludeVariantClass ) %>% 
            filter( FILTER          == "PASS" ) %>%
            filter( TUMOR_DEPTH     >= 25     ) %>%
            filter( TUMOR_ALT_COUNT >= 3      ) %>%
            filter( VAF             >= 0.05   )
        ##
        #==========================================================================================#
            write(sprintf("%s |  - VARIANT CLASS   : Non-Synonymous Only (6 Classes)", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            write(sprintf("%s |  - TUMOR DEPTH     : 25",   format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            write(sprintf("%s |  - TUMOR ALT DEPTH : 3",    format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            write(sprintf("%s |  - TUMOR VAF       : 0.05", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            write(sprintf("%s | VARIANT FILTERING FOR TMB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("Calculating Exonic Region Size FINISHED.")
        #==========================================================================================#
    }else{
        #---| MANUALLY VARIANT FILTERING |---------------------------------------------------------#
        ExcludeVariantClass <- c("3'Flank","3'UTR","5'Flank","5'UTR","Intron","RNA","IGR","Splice_Region")
        ## 
        MF_gcx <- MF %>%
            filter( Variant_Classification %nin% ExcludeVariantClass ) %>%
            filter( FILTER == "PASS" )
        ##
        if( gcxExcludeSilent )
        { 
            MF_gcx <- MF_gcx %>% filter( Variant_Classification != "Silent" )
            #======================================================================================#
                 write(sprintf("%s |  - VARIANT CLASS : Non-Synonymous Only", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            #======================================================================================#
        }else{
            #======================================================================================#
                 write(sprintf("%s |  - VARIANT CLASS : Non-Synonymous and Synonymous All", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            #======================================================================================#
        }        
        ##
        if( gcxExcludeClinvarPathogenic )
        { 
            MF_gcx <- MF_gcx %>% filter( CLINVAR_SIG != "Pathogenic" ) 
            #======================================================================================#
                 write(sprintf("%s |  - EXCLUDE CLINVAR PATHOGENIC : Yes ", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            #======================================================================================#
        }else{
            #======================================================================================#
                 write(sprintf("%s |  - EXCLUDE CLINVAR PATHOGENIC : No ", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            #======================================================================================#
        }
        ##
        if( gcxExcludeAlphaMissensePathogenic )
        { 
            MF_gcx <- MF_gcx %>% filter( ALPHA_MISSENSE_CLASS != "pathogenic" )
            #======================================================================================#
                 write(sprintf("%s |  - EXCLUDE ALPHA_MISSENSE PATHOGENIC : Yes ", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            #======================================================================================#
        }else{
            #======================================================================================#
                 write(sprintf("%s |  - EXCLUDE ALPHA_MISSENSE PATHOGENIC : No ", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            #======================================================================================#
        }
        ##
        MF_tmb <- rbind(
            MF_gcx %>%     
                filter( Variant_Type %nin% c("INS","DEL") ) %>%
                filter( TUMOR_DEPTH     >= gcx_DP         ) %>% 
                filter( TUMOR_ALT_COUNT >= gcx_TUMOR_ALT  ) %>%
                filter( VAF             >= gcx_SNV_VAF    ),
            MF_gcx %>%     
                filter( Variant_Type %in% c("INS","DEL")  ) %>%
                filter( TUMOR_DEPTH     >= gcx_DP         ) %>% 
                filter( TUMOR_ALT_COUNT >= gcx_TUMOR_ALT  ) %>%
                filter( VAF             >= gcx_INDEL_VAF  )
        ) %>%
            filter( G1K_AF <= gcx_POPAF         | is.na(G1K_AF)         ) %>%
            filter( G1K_EAS_AF <= gcx_POPAF     | is.na(G1K_EAS_AF)     ) %>%
            filter( gnomADe_AF <= gcx_POPAF     | is.na(gnomADe_AF)     ) %>%
            filter( gnomADe_EAS_AF <= gcx_POPAF | is.na(gnomADe_EAS_AF) ) 
        ##
        exonRegionCheck <- apply( MF_tmb[,c("Chromosome","Start_Position","End_Position")], 1, function(V) { 
            covRes %>% filter( Chr == V[1], Start <= as.numeric(V[2]), End >= as.numeric(V[2]) ) %>% nrow() 
        })
        rmRows <- which(exonRegionCheck == 0)
        if( length(rmRows) > 0 ) { MF_tmb <- MF_tmb[-rmRows, ] }
        #==========================================================================================#
            write(sprintf("%s |  - TUMOR DEPTH     : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S"), gcx_DP), LOG_FILE, append=T)
            write(sprintf("%s |  - TUMOR ALT DEPTH : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S"), gcx_TUMOR_ALT), LOG_FILE, append=T)
            write(sprintf("%s |  - TUMOR VAF       : %s (SNV), %s (INDEL)", format(now() ,format = "%Y-%m-%d %H:%M:%S"), gcx_SNV_VAF, gcx_INDEL_VAF), LOG_FILE, append=T)
            write(sprintf("%s |  - POPULATION AF   : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S"), gcx_POPAF), LOG_FILE, append=T)
            write(sprintf("%s | VARIANT FILTERING FOR TMB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("Calculating Exonic Region Size FINISHED.")
        #==========================================================================================#
    }
    ##
    rownames(MF_tmb) <- NULL
#--------------------------------------------------------------------------------------------------#

#---| TMB CALCULATION |----------------------------------------------------------------------------#
    TMB <- round(nrow(MF_tmb)/(LIB_LENGTH/1000000), 2)                        
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#
    write(sprintf("%s | TMB CALCULATED. ", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
#==================================================================================================#

#---| RESULT SUMMARY & SAVE |----------------------------------------------------------------------#
    if( variantCallMode == "tonly" ) 
    { 
        VarCallMethod = "Tumor-only"
    }else if( variantCallMode == "tonlygermlinefiltered" ){
        VarCallMethod = "TonlyGermlineFiltered"
    }else{
        VarCallMethod = ifelse( NORMAL_TYPE == "ORG", "ORG.NT-paired", "NT-paired")
    }
    #--------------------------------------------------------------------------#
    if( useReferenceFilter )
    {
        TMB.varFilterMethod   <- "TMB Harmonization Project"
        TMB.IncludeSilent     <- "No"
        TMB.pethogenicClinvar <- "No"
        TMB.pethogenicAM      <- "No"
        TMB.dp                <- 25
        TMB.alt               <- 3
        TMB.popaf             <- "No"
        TMB.vaf               <- 0.05        
    }else{
        TMB.varFilterMethod   <- "Gencurix"
        TMB.IncludeSilent     <- ifelse( gcxExcludeSilent,                  "No",  "Yes" )
        TMB.pethogenicClinvar <- ifelse( gcxExcludeClinvarPathogenic,       "Yes", "No"  )
        TMB.pethogenicAM      <- ifelse( gcxExcludeAlphaMissensePathogenic, "Yes", "No"  )
        TMB.dp                <- gcx_DP
        TMB.alt               <- gcx_TUMOR_ALT
        TMB.popaf             <- 0.01
        TMB.vaf               <- sprintf("SNV=%s;INDEL=%s",gcx_SNV_VAF,gcx_INDEL_VAF)
    }
    #--------------------------------------------------------------------------#
    TMB_RESULTS <- data.frame(
        seq_folder                 = seqFolder,
        seq_id                     = seqID,
        TMB                        = paste0(TMB, " Mut/Mb"),
        TMB_value                  = TMB,
        variants                   = nrow(MF_tmb),
        exon_region_Mb             = (LIB_LENGTH/1000000),
        variant_call_mode          = VarCallMethod,
        calculation_method         = TMB.varFilterMethod,
        include_silent_vars        = TMB.IncludeSilent,
        exclude_pathogenic_clinvar = TMB.pethogenicClinvar,
        exclude_pathogenic_am      = TMB.pethogenicAM,
        tumor_depth_cutoff         = TMB.dp,
        tumor_alt_count_cutoff     = TMB.alt,
        pop_af_cutoff              = TMB.popaf,
        vaf_cutoff                 = TMB.vaf
    )
    #--------------------------------------------------------------------------#
    write.table( TMB_RESULTS,
        sprintf("%s/%s.TMB.calculation.result.txt", TMB_RES_DIR, seqID ),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#
    write(sprintf("%s | TMB RESULT SAVED as FILE. ", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
#==================================================================================================#

#---| RESULT IMPORT INTO DATABASE -----------------------------------------------------------------#
    if( importDB )
    {
        #==========================================================================================# 
            write(sprintf("%s | RESULTS IMPORT INTO DB START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("RESULTS IMPORT INTO DB START.")
        #==========================================================================================# 
       
        dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
        preDataClear <- dbGetQuery(dbCon, 
            sprintf("DELETE FROM tmb WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", seqFolder, seqID, VarCallMethod)
        )
        writeNewData <- dbWriteTable(dbCon, name="tmb", value=TMB_RESULTS, row.names=FALSE, append=TRUE)

        dbDisconnect(dbCon)

        #==========================================================================================# 
            write(sprintf("%s | RESULTS IMPORT INTO DB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("RESULTS IMPORT INTO DB FINISHED.")
        #==========================================================================================#
    }
#----------------------------------------------------------------------------------------------------------------------------------------------------#
 






