

#   script version v1.1
#   update 2023-12-08
#   summary of change : merge default-vep, refseq-vep (remove gatk-funcotator)

#---| PACKAGES |-----------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages({
        library("optparse")
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("RMySQL")
        library("maftools")
        library("lubridate")
    })
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--BASE_DIR"), action="store", default=NULL, type="character", help="analysis base folder"),
        make_option(c("--BATCH_ID"), action="store", default=NULL, type="character", help="batch-id"),
        make_option(c("--SEQ_ID"), action="store", default=NULL, type="character", help="sample-id"),
        make_option(c("--ANNOTATE_TARGET_VCF_TYPE"), action="store", default=NULL, type="character", help="annotate target VCF type. 'NT_TS', 'NT_ORG', or 'Tonly'"),
        make_option(c("--CUTOFF_DEPTH"), action="store", default="10", type="character", help="somatic variant filter minimum depth"),
        make_option(c("--CUTOFF_ALT_READS"), action="store", default="2", type="character", help="somatic variant filter minimum alt-reads count"),
        make_option(c("--CUTOFF_POPAF"), action="store", default="0.01", type="character", help="somatic variant filter maximum population allele frequency"),
        make_option(c("--DB_IMPORT") , action="store", default=TRUE, type="logical", help="write filtered somatic variants maf into database")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BaseDir                     <- ARGS$BASE_DIR
    BatchID                     <- ARGS$BATCH_ID
    SeqID                       <- ARGS$SEQ_ID
    AnnotateTargetVcfType       <- ARGS$ANNOTATE_TARGET_VCF_TYPE
    WriteSomaticVarinatsOnlyMAF <- ARGS$WRITE_SOMATIC_VARIANT_ONLY_MAF
    CutOff_Depth                <- ARGS$CUTOFF_DEPTH
    CutOff_AltReadsCount        <- ARGS$CUTOFF_ALT_READS
    CutOff_PopAF                <- ARGS$CUTOFF_POPAF
    DatabaseImport              <- ARGS$DB_IMPORT
#--------------------------------------------------------------------------------------------------#

#---| Check Params and Set Default |----------------------------------------------------------------
    if( is.null(BaseDir) ){ BaseDir <- "/data/wes" }
    if( is.null(BatchID) ){ stop("No Batch-ID. REQUIRED.") }
    if( is.null(SeqID)   ){ stop("No Seq-ID. REQUIRED.") }
    if( is.null(AnnotateTargetVcfType)   ){ stop("No Annotate-Target-VCF-Type. REQUIRED.") }
    if( is.null(CutOff_Depth) ){ CutOff_Depth <- 10 }else{ CutOff_Depth <- as.numeric(CutOff_Depth) }
    if( is.null(CutOff_AltReadsCount) ){ CutOff_AltReadsCount <- 2 }else{ CutOff_AltReadsCount <- as.numeric(CutOff_AltReadsCount) }
    if( is.null(CutOff_PopAF) ){ CutOff_PopAF <- 0.01 }else{ CutOff_PopAF <- as.numeric(CutOff_PopAF) }
    if( is.null(DatabaseImport) ){ DatabaseImport <- FALSE }
#---------------------------------------------------------------------------------------------------

#---| Manual Params |-------------------------------------------------------------------------------
    # BaseDir               = "/data/wes"
    # BatchID               = "WES_25_09-WES_25_10_revised"
    # SeqID                 = "WES_25_09_09"
    # AnnotateTargetVcfType = "NT_TS"
    # CutOff_Depth          = 10
    # CutOff_AltReadsCount  = 2
    # CutOff_PopAF          = 0.01
    # DatabaseImport        = FALSE
#---------------------------------------------------------------------------------------------------

#==================================================================================================#    
    LogDir <- sprintf("%s/%s/%s/log", BaseDir, BatchID, SeqID)
    RunLog <- sprintf("%s/%s.maf.post.processing.%s.log", LogDir, SeqID, today())
    if( ! file.exists(RunLog) ){ system(sprintf("touch %s", RunLog)) }

    write(" ", RunLog, append=T)
    write("WES MAF Post-Processing log", RunLog, append=T)
    write("----------------------------------------------------------------------", RunLog, append=T)
    write(sprintf(" RUN DATE    %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), RunLog, append=T) 
    write(sprintf(" BATCH ID    %s", BatchID), RunLog, append=T) 
    write(sprintf(" SEQ ID      %s", SeqID), RunLog, append=T) 
    write("----------------------------------------------------------------------", RunLog, append=T)
    write(" ", RunLog, append=T)
#==================================================================================================# 

#---| INPUT MAF |-----------------------------------------------------------------------------------

    VcfDir <- sprintf("%s/%s/%s/vcf", BaseDir, BatchID, SeqID)
    
    if( AnnotateTargetVcfType == "NT_TS" ){
        VcfPrefix <- sprintf("%s.NT.somaticseq", SeqID)
        VariantCallMethodTag <- "TS_N"
        VariantCallMode <- ""
    }else if( AnnotateTargetVcfType == "NT_ORG" ){
        VcfPrefix <- sprintf("%s.NT.ORG.somaticseq", SeqID)
        VariantCallMethodTag <- "ORG_N"
    }else{
        VcfPrefix <- sprintf("%s.Tonly.somaticseq", SeqID)
        VariantCallMethodTag <- "Tonly"
    }
    
    FileList     <- list.files(VcfDir, pattern=VcfPrefix)
    NoPurecn_CHK <- FileList[ grep("no.purecn.vcf$", FileList) ]
    Purecn_CHK   <- FileList[ grep("purecn.vcf$", FileList)    ]
    if( length(Purecn_CHK) > 0 ){
        if( length(NoPurecn_CHK) > 0 ){
            InputMAF         <- sprintf("%s/%s.no.purecn.vep.norm.maf",  VcfDir, VcfPrefix)
            IsPurecnAddedMaf <- FALSE
        }else{
            InputMAF         <- sprintf("%s/%s.purecn.vep.norm.maf",  VcfDir, VcfPrefix)
            IsPurecnAddedMaf <- TRUE
        }
    }

    message(sprintf("--- INPUT MAF FILE ---"))
    message(sprintf("-> %s", InputMAF))
    message(" ")
#---------------------------------------------------------------------------------------------------

#==================================================================================================# 
    write(sprintf(" VARIANT CALL MODE   %s", VariantCallMethodTag), RunLog, append=T) 
    write(sprintf(" INPUT MAF FILE      %s", InputMAF), RunLog, append=T) 
    write("----------------------------------------------------------------------", RunLog, append=T)
    write(" ", RunLog, append=T)
#==================================================================================================# 

#---| LOAD DB DATA |-------------------------------------------------------------------------------#
    write(sprintf("[ %s ] Load Resources.", Sys.time()), RunLog, append=T)
    message(sprintf(" >> Load Resources..."))
    source("/data/wes/params/ruo_wes_db.R")
    dbCon          <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_resource )
    ginfo          <- dbGetQuery(dbCon, "SELECT * FROM gene_info_hgnc")
    ens2hgnc       <- dbGetQuery(dbCon, "SELECT * FROM hg19_ens2hgnc")
    kova           <- dbGetQuery(dbCon, "SELECT * FROM kova_hg19")
    essentialGenes <- dbGetQuery(dbCon, "SELECT * FROM ngs_panel_essential_genes WHERE esgd_id = '1'")
    dbDisconnect(dbCon)
    dbCon      <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
    addon_data <- dbGetQuery(dbCon, "SELECT * FROM custom_addon_data")
    dbDisconnect(dbCon)
#--------------------------------------------------------------------------------------------------#

#---| MAF FILEDS |---------------------------------------------------------------------------------#
    
    # MAF required columns
    RequiredFields <- c(
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', "Strand","Reference_Allele", 
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode'
    )
    # Additional columns with Rename Targets
    RenameFields <- c(
        'NCBI_Build'='REFERENCE_GENOME',
        'FILTER'='FILTER',
        'Existing_variation'='dbSNP',
        't_depth'='TUMOR_DEPTH',
        't_ref_count'='TUMOR_REF_COUNT',
        't_alt_count'='TUMOR_ALT_COUNT',
        'n_depth'='MATCHED_NORMAL_DEPTH',
        'n_ref_count'='MATCHED_NORMAL_REF_COUNT',
        'n_alt_count'='MATCHED_NORMAL_ALT_COUNT',
        'VAF'='VAF',
        'Consequence'='CONSEQUENCE',
        'SYMBOL'='GENE_SYMBOL',
        'Transcript_ID'='ENSEMBL_TRANSCIRPT_ID',
        'Gene'='ENSEMBL_GENE_ID',
        'Refseq_Feature'='REFSEQ_TRANSCRIPT_ID',
        'Refseq_Gene'='ENTREZ',
        'Refseq_HGVSc'='REFSEQ_HGVSc',
        'Refseq_HGVSp'='REFSEQ_HGVSp',
        'HGVSc'='HGVSc',
        'HGVSp'='HGVSp',
        'HGVSp_Short'='HGVSp_Short',
        'HGVSg'='HGVSg',
        'HGVS_OFFSET'='HGVS_OFFSET',
        'EXON'='EXON',
        'INTRON'='INTRON',
        'cDNA_position'='cDNA_POSITION',
        'CDS_position'='CDS_POSITION',
        'Amino_acids'='AMINO_ACID_CHANGE',
        'Codons'='CODON_CHANGE',
        'Protein_position'='PROTEIN_POSITION',
        'CANONICAL'='CANONICAL',
        'AF'='AF',
        'EAS_AF'='EAS_AF',
        'gnomADe_AF'='GNOMADe_AF',
        'gnomADe_EAS_AF'='GNOMADe_EAS_AF',
        'MAX_AF'='MAX_AF',
        'MAX_AF_POPS'='MAX_AF_POPS',
        'SIFT'='SIFT',
        'PolyPhen'='POLYPHEN',
        'GENE_PHENO'='GENE_PHENO',
        'COSMIC'='COSMIC_ID',
        'COSMIC_HGVSC'='COSMIC_HGVSc',
        'COSMIC_LEGACY_ID'='COSMIC_LEGACY_ID',
        'COSMIC_TIER'='COSMIC_TIER',
        'ClinVar_CLNSIG'='CLINVAR_SIGINIFICANCE',
        'ClinVar_CLNDN'='CLINVAR_DISEASE',
        'ClinVar_ORIGIN'='CLINVAR_ORIGIN',
        'am_class'='ALPHA_MISSENSE_CLASS',
        'am_pathogenicity'='ALPHA_MISSENSE_PATHOGENICITY',
        'PureCN.PS'='PUCRECN_POSTERIOR_SOMATIC',
        'PureCN.ML.C'='PURECN_CN',
        'PureCN.CS'='PURECN_SUBCLONAL_CN',
        'PureCN.ML.LOH'='PURECN_LOH',
        'PureCN.GHOMOZYGOUS'='PURECN_GERMLINE_HOMOZYGOUS',
        'PureCN.LR'='PURECN_LOG_RATIO',
        'PureCN.CF'='PURECN_CELLULAR_FRACTION',
        'IMPACT'='IMPACT',
        'SYMBOL_SOURCE'='SYMBOL_SOURCE',
        'HGNC_ID'='HGNC_ID',
        'BIOTYPE'='BIOTYPE',
        'CCDS'='CCDS',
        'SWISSPROT'='SWISSPROT',
        'UNIPARC'='UNIPARC',
        'DISTANCE'='DISTANCE',
        'all_effects'='ALL_EFFECTS',
        'DOMAINS'='DOMAINS'  , 
        'POPAF'='POPAF',
        'PON'='PON',
        'TLOD'='TLOD',
        'VKEY'='VKEY'
    )
    # initial selected columns
    SelectFields <- c(RequiredFields,RenameFields)
    names(SelectFields) <- NULL
    SelectFields <- c(SelectFields[1:69], "PURECN_FLAG", SelectFields[70:length(SelectFields)])
    # exclude variant classes in somatic variants decision
    RemoveClasses <- c("Silent","Intron","5'Flank","3'Flank","Targeted_Region","IGR","3'UTR","5'UTR","RNA")
    # KOVA modification
    kova <- kova %>% mutate( VKEY = paste0(Chr, ":", Pos, "_", Ref, ">", Alt) ) %>% dplyr::rename(KOVA_AF=AF)
    # NGS panel test essential genes
    essGenes <- essentialGenes$hgnc_symbol
    # HLA blacklist genes
    blacklistGenes <- unlist(strsplit(addon_data[which(addon_data$tag == "hla_genes"), "data_value"],";"))
    # final selected columns
    OutputMafFields <- c(       
        "Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
        "Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","Variant_Classification","Variant_Type","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode",
        "REFERENCE_GENOME","FILTER","dbSNP","TUMOR_DEPTH","TUMOR_REF_COUNT","TUMOR_ALT_COUNT",
        "MATCHED_NORMAL_DEPTH","MATCHED_NORMAL_REF_COUNT","MATCHED_NORMAL_ALT_COUNT",
        "VAF","CONSEQUENCE","HGVSc","HGVSp","HGVSp_Short","HGVSg","HGVS_OFFSET","REFSEQ_HGVSc","REFSEQ_HGVSp","EXON","INTRON",
        "cDNA_POSITION","CDS_POSITION","AMINO_ACID_CHANGE","CODON_CHANGE","PROTEIN_POSITION",
        "GENE_SYMBOL","HGNC_NAME","HGNC_LOCUS_TYPE","HGNC_ID","ENTREZ","ENSEMBL_GENE_ID","ENSEMBL_TRANSCIRPT_ID","REFSEQ_TRANSCRIPT_ID","CANONICAL",
        "LOCATION","SYMBOL_SOURCE","BIOTYPE","SWISSPROT","UNIPARC",
        "AF","EAS_AF","GNOMADe_AF","GNOMADe_EAS_AF","MAX_AF","MAX_AF_POPS","KOVA_AF",
        "COSMIC_ID","COSMIC_LEGACY_ID","COSMIC_TIER","COSMIC_HGVSc",
        "PUCRECN_POSTERIOR_SOMATIC","PURECN_CN","PURECN_SUBCLONAL_CN","PURECN_LOH","PURECN_GERMLINE_HOMOZYGOUS","PURECN_LOG_RATIO",
        "PURECN_CELLULAR_FRACTION","PURECN_FLAG",
        "CLINVAR_SIGINIFICANCE","CLINVAR_DISEASE","CLINVAR_ORIGIN","ALPHA_MISSENSE_CLASS","ALPHA_MISSENSE_PATHOGENICITY",
        "IMPACT","SIFT","POLYPHEN","GENE_PHENO","CCDS","DISTANCE","ALL_EFFECTS","DOMAINS","POPAF","PON","TLOD","VKEY"
    )

#--------------------------------------------------------------------------------------------------#

#---| MAF PROCESSING |------------------------------------------------------------------------------
    write(sprintf("[ %s ] MAF-PROCESSING Start...", Sys.time()), RunLog, append=T)
    message(sprintf(" >> MAF-PROCESSING Start..."))

    NumVars <- system(sprintf("cat %s | wc -l", InputMAF), intern=T) %>% as.numeric() -1 

    if( NumVars <=  0 ){
        write(sprintf("[ %s ] NO VARIANTS.", Sys.time()), RunLog, append=T)
        message(sprintf(" >> NO VARIANTS."))
        quit(save="no", status=0)
    }

    ## read maf as text-table
    INITIAL_MAF <- read.delim(InputMAF, skip=1, header=TRUE) %>% mutate(
        VKEY = paste0(Chromosome, ":", Start_Position, "_", Reference_Allele, ">", Tumor_Seq_Allele2),
        VAF  = round(t_alt_count/(t_alt_count+t_ref_count),3)
    )

    write(sprintf("[ %s ] %s variants found in MAF.", Sys.time(), nrow(INITIAL_MAF)), RunLog, append=T)
    message(sprintf(" >> %s variants found in MAF.", nrow(INITIAL_MAF)))
    
    UNFILTERED_MAF <- INITIAL_MAF %>% 
        dplyr::rename(!!!setNames(names(RenameFields), RenameFields)) %>% 
        dplyr::mutate( PURECN_FLAG="" ) %>%
        dplyr::select(all_of(SelectFields)) 

    UNFILTERED_MAF$PURECN_FLAG <- sapply( UNFILTERED_MAF$PUCRECN_POSTERIOR_SOMATIC, function(pps) {
        if( is.na(pps) ){ "Not_Evaluated" 
        }else if( as.numeric(pps) > 0.9 ){ "PureCN_High_Confidence"
        }else{ "PureCN_Less_Somatic"}
    })
    UNFILTERED_MAF$HGNC_ID <- ens2hgnc[match(UNFILTERED_MAF$ENSEMBL_GENE_ID, ens2hgnc$ens_geneid), "hgnc_id"]
    UNFILTERED_MAF <- UNFILTERED_MAF %>% 
        left_join(., kova[,c("VKEY","KOVA_AF")], by="VKEY" ) %>% 
        left_join(., 
            ginfo[,c("hgnc_id","symbol","name","locus_type","location")] %>% unique() %>% 
                dplyr::rename(HGNC_ID=hgnc_id, HGNC_SYMBOL=symbol, HGNC_NAME=name, HGNC_LOCUS_TYPE=locus_type, LOCATION=location), 
            by="HGNC_ID" 
        )

    FILTERED_MAF <- UNFILTERED_MAF %>% 
        filter( FILTER %in% c("PASS","LowQual") ) %>% 
        filter( TUMOR_REF_COUNT+TUMOR_ALT_COUNT >= CutOff_Depth, TUMOR_ALT_COUNT >= CutOff_AltReadsCount ) 

    EssentialPanelGeneMaf <- UNFILTERED_MAF %>%  filter( FILTER %in% c("PASS","LowQual"), HGNC_SYMBOL %in% essGenes )
    if( nrow(EssentialPanelGeneMaf) > 0 ){
        FILTERED_MAF <- unique(rbind(FILTERED_MAF, EssentialPanelGeneMaf))
    }

    FILTERED_MAF = FILTERED_MAF[,OutputMafFields]

    write(sprintf("[ %s ] Filtered-Raw-MAF created.", Sys.time()), RunLog, append=T)
    message(sprintf(" >> Filtered-Raw-MAF created."))

#---------------------------------------------------------------------------------------------------

#---| SOMATIC VARIANTS MAF : RESEARCH & REPORT |----------------------------------------------------

    SOMATIC_RESEARCH_MAF <- rbind(
        FILTERED_MAF %>% 
            filter( Variant_Classification %nin% RemoveClasses ) %>%
            filter( MAX_AF < CutOff_PopAF | is.na(MAX_AF) ) %>% 
            filter( KOVA_AF < CutOff_PopAF | is.na(KOVA_AF) ) %>% 
            filter( Variant_Type %nin% c("INS","DEL") ) %>% 
            filter( VAF >= 0.02 ),
        FILTERED_MAF %>% 
            filter( Variant_Classification %nin% RemoveClasses ) %>%
            filter( MAX_AF < CutOff_PopAF | is.na(MAX_AF) ) %>% 
            filter( KOVA_AF < CutOff_PopAF | is.na(KOVA_AF) ) %>% 
            filter( Variant_Type %in% c("INS","DEL") ) %>% 
            filter( VAF >= 0.05 )
    ) %>% as.data.frame() 

    SOMATIC_REPORT_MAF <- SOMATIC_RESEARCH_MAF %>% filter( FILTER == "PASS" )

    write(sprintf("[ %s ] Somatic-Variants-MAF created.", Sys.time()), RunLog, append=T)
    message(sprintf(" >> Somatic-Variants-MAF created."))

#---------------------------------------------------------------------------------------------------

#---| VARIANT-CLASS STATS |-------------------------------------------------------------------------

    RawVarSyn    <- FILTERED_MAF %>% filter( Variant_Classification %in%  RemoveClasses )
    RawVarNonSyn <- FILTERED_MAF %>% filter( Variant_Classification %nin% RemoveClasses )

    RawVarSyn_Stats <- RawVarSyn %>% group_by(Variant_Classification, Variant_Type) %>% 
        reframe(stats=length(HGVSg)) %>% as.data.frame() %>% 
        mutate( seq_folder=BatchID, seq_id=SeqID, variant_group='synonymous' )

    RawVarNonSyn_Stats <- RawVarNonSyn %>% group_by(Variant_Classification, Variant_Type) %>% 
        reframe(stats=length(HGVSg)) %>% as.data.frame() %>% 
        mutate( seq_folder=BatchID, seq_id=SeqID, variant_group='non_synonymous' )

    SomaticVar_Stats <- SOMATIC_REPORT_MAF %>% group_by(Variant_Classification, Variant_Type) %>% 
        reframe(stats=length(HGVSg)) %>% as.data.frame() %>% 
        mutate( seq_folder=BatchID, seq_id=SeqID, variant_group='filtered_somatic_only' )

    VARIANT_CLASS_STATS <- rbind( RawVarSyn_Stats, RawVarNonSyn_Stats, SomaticVar_Stats )

    write(sprintf("[ %s ] Variant-Class-Stats created.", Sys.time()), RunLog, append=T)
    message(sprintf(" >> Variant-Class-Stats created."))

#---------------------------------------------------------------------------------------------------

#---| WRITE RESULTS |-------------------------------------------------------------------------------

    Output_FilteredRaw_MAF     <- gsub("norm.maf","norm.filtered.raw.maf"    , InputMAF)
    Output_SomaticResearch_MAF <- gsub("norm.maf","norm.somatic.research.maf", InputMAF)
    Output_SomaticReport_MAF   <- gsub("norm.maf","norm.somatic.report.maf"  , InputMAF)
    Output_VarClassStats       <- gsub("norm.maf","variant.class.summary.tsv", InputMAF)

    # raw-maf 
    write.table( FILTERED_MAF, Output_FilteredRaw_MAF, quote=F, col.names=T, row.names=F, sep="\t")
   
    # somatic-research-maf 
    write.table( SOMATIC_RESEARCH_MAF, Output_SomaticResearch_MAF, quote=F, col.names=T, row.names=F, sep="\t")
 
    # somatic-report-maf 
    write.table( SOMATIC_REPORT_MAF, Output_SomaticReport_MAF, quote=F, col.names=T, row.names=F, sep="\t")

    # variant-class
    write.table( VARIANT_CLASS_STATS, Output_VarClassStats, quote=F, col.names=T, row.names=F, sep="\t")

    write(sprintf("[ %s ] Results saved as files.", Sys.time()), RunLog, append=T)
    message(sprintf(" >> Results saved as files."))

#--------------------------------------------------------------------------------------------------#

#---| IMPORT INTO DDATABASE |-----------------------------------------------------------------------

    if( DatabaseImport ){

        # soamtic-vars db-table
        DB_IMPORT_MAF <- SOMATIC_REPORT_MAF %>% mutate(
            seq_folder        = BatchID, 
            seq_id            = SeqID,
            variant_call_mode = VariantCallMethodTag
        )
        # variant-class-stats db-table
        DB_IMPORT_STATS <- VARIANT_CLASS_STATS %>% mutate(
            seq_folder        = BatchID, 
            seq_id            = SeqID,
            variant_call_mode = VariantCallMethodTag
        )

        
        dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )

        # delete pre-exists records         
        ClearStats <- dbGetQuery(dbCon, 
            sprintf("DELETE FROM variants_summary WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", BatchID, SeqID, VariantCallMethodTag)
        )
        ClearVars <- dbGetQuery(dbCon, 
            sprintf("DELETE FROM variants_somatic WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", BatchID, SeqID, VariantCallMethodTag)
        )
        # write new=records
        WriteStats <- dbWriteTable(dbCon, name="variants_summary", value=DB_IMPORT_STATS, row.names=FALSE, append=TRUE)
        WriteVars  <- dbWriteTable(dbCon, name="variants_somatic", value=DB_IMPORT_MAF,   row.names=FALSE, append=TRUE)

        dbDisconnect(dbCon)

        write(sprintf("[ %s ] Results imported into database.", Sys.time()), RunLog, append=T)
        message(sprintf(" Results imported into database."))

    }else{

        write(sprintf("[ %s ] Database import skipped.", Sys.time()), RunLog, append=T)
        message(sprintf(" Database import skipped."))

    }

#---------------------------------------------------------------------------------------------------

#==================================================================================================#
    write("   ", RunLog, append=T)
    write(sprintf("[ %s ] MAF-PROCESSING DONE.", Sys.time()), RunLog, append=T)
    write("----------------------------------------------------------------------", RunLog, append=TRUE ) 
    write("   ", RunLog, append=T)
#==================================================================================================# 


