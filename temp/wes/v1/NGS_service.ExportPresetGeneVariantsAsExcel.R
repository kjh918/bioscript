
#---| PACKAGES |-----------------------------------------------------------------------------------#
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("RMySQL"))
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--BASE_DIR"      ), action="store", default=NA, type="character", help="BaseDir"         ),
        make_option(c("--SEQ_FOLDER"    ), action="store", default=NA, type="character", help="SeqFolder"       ),
        make_option(c("--EXCLUDE_SEQ_ID"), action="store", default=NA, type="character", help="SeqID to exclude"),
        make_option(c("--INCLUDE_SEQ_ID"), action="store", default=NA, type="character", help="SeqID to include")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR       <- ARGS$BASE_DIR
    SEQ_FOLDER     <- ARGS$SEQ_FOLDER
    EXCLUDE_SEQ_ID <- ARGS$EXCLUDE_SEQ_ID
    INCLUDE_SEQ_ID <- ARGS$INCLUDE_SEQ_ID
#--------------------------------------------------------------------------------------------------#

#---| DATA LOAD |----------------------------------------------------------------------------------#
    message(sprintf(">>>>> SEQ_FOLDER = %s", SEQ_FOLDER))
    #--------------------------------------------------------------------------#
    dbCon     = dbConnect(dbDriver("MySQL"), host='192.168.0.34', port=3306, user='gcx', password='gencurix!!', db = 'gcx_ngs_service' )
    orderInfo = dbGetQuery(dbCon, sprintf("SELECT * FROM seqfolder_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    sinfo     = dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_folder = '%s'", SEQ_FOLDER))
    if( !is.na(INCLUDE_SEQ_ID) )
    {
        include_samples   <- unlist(strsplit(INCLUDE_SEQ_ID, ","))
        QueryIncludeSeqID <- paste(paste0("'", include_samples, "'"), collapse=",")
        sinfo_extra <- dbGetQuery(dbCon, sprintf("SELECT * FROM seqid_info WHERE seq_id IN (%s)", QueryIncludeSeqID))
        sinfo <- rbind(sinfo, sinfo_extra)
    }
    cinfo = dbGetQuery(dbCon, sprintf("SELECT seqfolder_info.seq_folder,facility_info.facility_name,seqfolder_info.ngs_order_id FROM seqfolder_info INNER JOIN facility_info ON seqfolder_info.facility_id = facility_info.facility_id WHERE seq_folder = '%s'", SEQ_FOLDER))
    AnalysisSeqID = paste(paste0("'", sinfo$seq_id, "'"), collapse=",")
    var_pf        = dbGetQuery(dbCon, sprintf("SELECT * FROM gcx_wes.variants_somatic WHERE seq_id IN (%s)", AnalysisSeqID))
    dbDisconnect(dbCon)
    #--------------------------------------------------------------------------#
    varSinfo <- sinfo %>% filter( seq_id %in%  var_pf$seq_id )
    #--------------------------------------------------------------------------#
    orderID <- cinfo$ngs_order_id
    #--------------------------------------------------------------------------#
    if( all(!is.na(EXCLUDE_SEQ_ID)) )
    {
        ExcludeSeqID <- unlist(strsplit(EXCLUDE_SEQ_ID, ","))
        varSinfo     <- varSinfo %>% filter( seq_id %nin% ExcludeSeqID )
    }
    #--------------------------------------------------------------------------#
    load("/storage/home/kangsm/myDB/rds/PresetGenes.Rdata")
#--------------------------------------------------------------------------------------------------#

#---| FUNCTIONS |----------------------------------------------------------------------------------#
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
    #--------------------------------------------------------------------------#
    PRESET_GENE_VARIANTS_TABLE2 <- function( PF_VARS, CANCER_TISSUE, GENE_LIST, SAMPLE_INFO, PRESET_NAME, VARIANT_CALL_METHOD )
    {
        seqIdList <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% .$seq_id

        cancerType <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% .$sample_info %>% unique()

        if( VARIANT_CALL_METHOD == "TS_N"  ){ VCM = "Matched Normal (Tissue)"   }
        if( VARIANT_CALL_METHOD == "ORG_N" ){ VCM = "Matched Normal (Organoid)" }
        if( VARIANT_CALL_METHOD == "Tonly" ){ VCM = "Tumor-only" }
        if( VARIANT_CALL_METHOD == "TonlyGermlineFiltered" ){ VCM = "Tumor-only Germline Filtered" }

        presetGeneVarsTable <- PF_VARS %>% 
            filter( variant_call_mode == VARIANT_CALL_METHOD, seq_id %in% seqIdList, HGNC_SYMBOL %in% GENE_LIST ) %>%
            mutate( SampleName = SAMPLE_INFO[match(seq_id, SAMPLE_INFO$seq_id), c("sample_name")], varDNA=paste(Refseq_TRANSCRIPT,HGVSc, sep=":")) %>%
            arrange( seq_id, Chromosome, Start_Position) %>% 
            select(c("SampleName","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) %>%
            dplyr::rename(
                'Gene'            = HGNC_SYMBOL,
                'DNA Change'      = varDNA, 
                'PROTEIN Change'  = HGVSp_Short,
                'Depth'           = DEPTH, 
                'Tumor Alt Count' = TUMOR_ALT_COUNT
            ) 
        return(presetGeneVarsTable)
    }
#--------------------------------------------------------------------------------------------------#

#---| PRESET GENES BY CANCER TYPE |----------------------------------------------------------------#
    cancerTypeList <- lapply(indexing(varSinfo$sample_tissue), function(st) 
    {
        genelist = ifelse( st %in% names(c11_preset_genes), st, "pan_cancers" )
        genelist = c(genelist, "mutation_panel")
        return(genelist)
    })
    message(">>>>>")
    message(">>>>> Detected Tissue-Info and Gene-Presets")
    for( k in 1:length(cancerTypeList) )
    {
        message(sprintf("TISSUE / GENE-PRESET : %s / %s , %s", names(cancerTypeList)[k], cancerTypeList[[k]][1], cancerTypeList[[k]][2]))
    }
    message(">>>>>")
#--------------------------------------------------------------------------------------------------#

#---|
    variantCallModeIndex <- unique(var_pf$variant_call_mode)
    message(">>>>> Variant Call Modes")
    for( k in 1:length(variantCallModeIndex) )
    { message(sprintf(">>>>> %s. %s", k, variantCallModeIndex[k])) }
    message(">>>>>")
    #--------------------------------------------------------------------------#
    analysisVisGroup    <- lapply( indexing(variantCallModeIndex), function(vcm){ 
        var_pf %>% filter( variant_call_mode == vcm ) %>% .$seq_id %>% unique() }
    ) 
    ##
    analysisVisGroupTag <- c(
        "TS_N"  = "Matched.normal.TISSUE",
        "ORG_N" = "Matched.normal.ORGANOID",
        "Tonly" = "Tumor.only",
        "TonlyGermlineFiltered" = "Tumor.only.Germline.Filtered"
    )
    #--------------------------------------------------------------------------#
    tissueCode              <- unique(varSinfo[,c("sample_tissue","sample_info")])
    noTsName <- which(is.na(tissueCode$sample_tissue))
    for( ntn in noTsName ){ tissueCode[ntn, 1] = tissueCode[ntn ,2] }
    SampleTissueCode        <- tissueCode$sample_info
    names(SampleTissueCode) <- tissueCode$sample_tissue
    #--------------------------------------------------------------------------#
    AnalysisGroups <- lapply(indexing( varSinfo %>% filter(sample_group != "") %>% .$sample_group ), function(sg) 
    {
        var_pf %>% 
            filter( seq_id %in% (varSinfo %>% filter( sample_group == sg ) %>% .$seq_id) ) %>% 
            .$variant_call_mode %>% unique() 
    })
    ##
    MatchedAnalysisGroup <- AnalysisGroups[ which(lapply(AnalysisGroups, function(ag) all(ag == "Tonly") ) == FALSE) ]
    TonlyAnalysisGroup   <- AnalysisGroups[ which(lapply(AnalysisGroups, function(ag) all(ag == "Tonly") ) == TRUE ) ]
    ##
    TonlyAnalysisSamples <- lapply(indexing(varSinfo %>% filter(sample_group == "") %>% .$sample_info), function(sf) 
    { varSinfo %>% filter( sample_group == "", sample_info == sf ) %>% .$seq_id %>% unique() })
#--------------------------------------------------------------------------------------------------#
  
  
#---| 
    for( VARIANT_CALL_METHOD in variantCallModeIndex )
    {   
        message(sprintf(">>>>> Generating Variant List for Mode '%s' ", analysisVisGroupTag[VARIANT_CALL_METHOD]) )
        #----------------------------------------------------------------------#
        vidx <- unlist(lapply( 1:length(cancerTypeList), function(i) paste(names(cancerTypeList)[i], cancerTypeList[[i]], sep=";") ))
        #----------------------------------------------------------------------#
        VAR_LIST <- lapply( indexing(vidx), function(vx)
        {
            PSG = unlist(strsplit(vx,";"))[2]
            CTS = unlist(strsplit(vx,";"))[1]

            matchMethod.matchTissue.sinfo = varSinfo %>% filter( 
                seq_id        %in% analysisVisGroup[[VARIANT_CALL_METHOD]], 
                sample_tissue %in% CTS 
            )
            matchVARS <- var_pf %>% filter( 
                variant_call_mode  ==  VARIANT_CALL_METHOD,
                seq_id            %in% matchMethod.matchTissue.sinfo$seq_id, 
                HGNC_SYMBOL       %in% c11_preset_genes[[ PSG ]] 
            ) 
            PSG_VAR_TABLE <- PRESET_GENE_VARIANTS_TABLE2(
                PF_VARS       = var_pf,
                CANCER_TISSUE = CTS,
                GENE_LIST     = c11_preset_genes[[ PSG ]],
                SAMPLE_INFO   = matchMethod.matchTissue.sinfo,
                PRESET_NAME   = PSG,
                VARIANT_CALL_METHOD = VARIANT_CALL_METHOD
            )
            return(PSG_VAR_TABLE)
        })
        #----------------------------------------------------------------------#
        for(k in 1:length(SampleTissueCode))
        {
            names(VAR_LIST) = gsub( names(SampleTissueCode)[k], SampleTissueCode[k], names(VAR_LIST) )
        }
        names(VAR_LIST) = gsub(";", "_", names(VAR_LIST))
        #----------------------------------------------------------------------#
        exportFileName <- sprintf("%s/%s/meta/%s.Preset.Genes.Variants.%s.xlsx", BASE_DIR, SEQ_FOLDER, orderID, analysisVisGroupTag[VARIANT_CALL_METHOD] )
        #----------------------------------------------------------------------#
        openxlsx::write.xlsx( VAR_LIST, exportFileName, rowNames=F, overwrite=T)
        message(sprintf(">>>>> Generating Variant List for Mode '%s' ==>> Finished and Result Saved.", analysisVisGroupTag[VARIANT_CALL_METHOD]) )
    }
    #--------------------------------------------------------------------------#
    message(">>>>>")
    message("Exporting Preset-Gene-Variant-Result-Table as Excel-file has Done.")
    message("     ")
#--------------------------------------------------------------------------------------------------#




    # source("/storage/home/kangsm/shinyWeb/resources/Fun_QC_Report.R")
    # source("/storage/home/kangsm/shinyWeb/resources/Fun_VARIANTS.R")

    # PRESET_GENE_VARIANTS_TABLE2 <- function( PF_VARS, CANCER_TISSUE, GENE_LIST, SAMPLE_INFO, PRESET_NAME, VARIANT_CALL_METHOD )
    # {
    #     seqIdList <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% .$seq_id

    #     cancerType <- SAMPLE_INFO %>% filter( sample_tissue == CANCER_TISSUE ) %>% .$sample_info %>% unique()

    #     if( VARIANT_CALL_METHOD == "TS_N"  ){ VCM = "Matched Normal (Tissue)"   }
    #     if( VARIANT_CALL_METHOD == "ORG_N" ){ VCM = "Matched Normal (Organoid)" }
    #     if( VARIANT_CALL_METHOD == "Tonly" ){ VCM = "Tumor-only" }

    #     presetGeneVarsTable <- PF_VARS %>% 
    #         filter( variant_call_mode == VARIANT_CALL_METHOD, seq_id %in% seqIdList, HGNC_SYMBOL %in% GENE_LIST ) %>%
    #         mutate( SampleName = SAMPLE_INFO[match(seq_id, SAMPLE_INFO$seq_id), c("sample_name")], varDNA=paste(Refseq_TRANSCRIPT,HGVSc, sep=":")) %>%
    #         arrange( seq_id, Chromosome, Start_Position) %>% 
    #         select(c("SampleName","HGNC_SYMBOL", "varDNA","HGVSp_Short","DEPTH","TUMOR_ALT_COUNT","VAF")) %>%
    #         dplyr::rename(
    #             'Gene'            = HGNC_SYMBOL,
    #             'DNA Change'      = varDNA, 
    #             'PROTEIN Change'  = HGVSp_Short,
    #             'Depth'           = DEPTH, 
    #             'Tumor Alt Count' = TUMOR_ALT_COUNT
    #         ) 
    #     return(presetGeneVarsTable)
    # }
    # indexing = function(X) { sapply(unique(X), function(y) list(y)) }


    # varSinfo <- sinfo[ which(sinfo$seq_id %in% var_pf$seq_id), ]






