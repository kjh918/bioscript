#   script version v1.4
#   update 2026-06-04
#   summary of change : Fix 'Names must be unique' error by removing pre-existing target columns before DB join/rename

#---| PACKAGES |-----------------------------------------------------------------------------------#
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("lubridate"))
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--seqFolder"), action="store", type="character", help="$SEQ_FOLDER"),
        make_option(c("--seqID")    , action="store", type="character", help="$SEQ_ID"),
        make_option(c("--baseDir")  , action="store", type="character", help="$BASE_DIR"),
        make_option(c("--nType")    , action="store", type="character", help="type of normal sample"),
        make_option(c("--writeSomaticVarinatsOnlyMAF") , action="store", type="logical", help="write another maf file including only filtered somatic variants"),
        make_option(c("--DP")       , action="store", type="double", help="variant filter: DEPTH "),
        make_option(c("--ALT")      , action="store", type="double", help="variant filter: ALT DEPTH "),
        make_option(c("--POPAF")    , action="store", type="double", help="variant filter: POP AF "),
        make_option(c("--writeDB")  , action="store", type="logical", help="write filtered variants maf into database"),
        make_option(c("--variantCallMode"), action="store", type="character", help="variant calling run mode"),
        make_option(c("--GERMLINE") , action="store", type="logical", help="Germline variant MAF or not.")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    
    if(is.null(ARGS$seqFolder)) stop("[ERROR] --seqFolder argument is missing.")
    if(is.null(ARGS$seqID)) stop("[ERROR] --seqID argument is missing.")
    if(is.null(ARGS$baseDir)) stop("[ERROR] --baseDir argument is missing.")
    if(is.null(ARGS$nType)) stop("[ERROR] --nType argument is missing.")
    if(is.null(ARGS$writeSomaticVarinatsOnlyMAF)) stop("[ERROR] --writeSomaticVarinatsOnlyMAF argument is missing.")
    if(is.null(ARGS$DP)) stop("[ERROR] --DP argument is missing.")
    if(is.null(ARGS$ALT)) stop("[ERROR] --ALT argument is missing.")
    if(is.null(ARGS$POPAF)) stop("[ERROR] --POPAF argument is missing.")
    if(is.null(ARGS$writeDB)) stop("[ERROR] --writeDB argument is missing.")
    if(is.null(ARGS$variantCallMode)) stop("[ERROR] --variantCallMode argument is missing.")
    if(is.null(ARGS$GERMLINE)) stop("[ERROR] --GERMLINE argument is missing.")
    #--------------------------------------------------------------------------#
    SEQ_FOLDER                  <- ARGS$seqFolder
    SEQ_ID                      <- ARGS$seqID
    BASE_DIR                    <- ARGS$baseDir
    NORMAL_TYPE                 <- ARGS$nType
    writeSomaticVarinatsOnlyMAF <- ARGS$writeSomaticVarinatsOnlyMAF
    CutOff_DP                   <- ARGS$DP
    CutOff_ALT_DP               <- ARGS$ALT
    CutOff_POP_AF               <- ARGS$POPAF
    writeDB                     <- ARGS$writeDB
    variantCallRunMode          <- ARGS$variantCallMode
    is_gemlineMAF               <- ARGS$GERMLINE
#--------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTIONS |-----------------------------------------------------------------------#
    source("/data/wes/params/ruo_wes_db.R")
#--------------------------------------------------------------------------------------------------#

#---| LOAD DB DATA |-------------------------------------------------------------------------------#
    dbCon    <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_resource )
    ginfo    <- dbGetQuery(dbCon, "SELECT * FROM gene_info_hgnc")
    ens2hgnc <- dbGetQuery(dbCon, "SELECT * FROM hg19_ens2hgnc")
    kova     <- dbGetQuery(dbCon, "SELECT * FROM kova_hg19")
    essentialGenes <- dbGetQuery(dbCon, "SELECT * FROM ngs_panel_essential_genes WHERE esgd_id = '1'")
    dbDisconnect(dbCon)

    dbCon      <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
    addon_data <- dbGetQuery(dbCon, "SELECT * FROM custom_addon_data")
    dbDisconnect(dbCon)
    
    kova <- kova %>% mutate( VKEY = paste0(Chr, ":", Pos, "_", Ref, ">", Alt) ) %>% dplyr::rename(KOVA_AF=AF)
    essGenes <- essentialGenes$hgnc_symbol
    blacklistGenes <- unlist(strsplit(addon_data[which(addon_data$tag == "hla_genes"), "data_value"],";"))
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS & LOG |------------------------------------------------------------------------------#
    DATA_DIR <- sprintf("%s/%s", BASE_DIR,  SEQ_FOLDER)
    VCF_DIR  <- sprintf("%s/%s/vcf", DATA_DIR, SEQ_ID)
    LOG_DIR  <- sprintf("%s/%s/log", DATA_DIR, SEQ_ID)

    LOG_FILE <- sprintf("%s/%s.maf.post-processing.%s.log", LOG_DIR, SEQ_ID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
#--------------------------------------------------------------------------------------------------#

#---| MAF FIELDS |---------------------------------------------------------------------------------#
    RequiredFields <- c(
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', "Strand","Reference_Allele", 
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode'
    )
    
    RenameFields <- c(
        'NCBI_Build'='REFERENCE_GENOME', 'FILTER'='FILTER', 'Existing_variation'='dbSNP',
        't_depth'='TUMOR_DEPTH', 't_ref_count'='TUMOR_REF_COUNT', 't_alt_count'='TUMOR_ALT_COUNT',
        'n_depth'='MATCHED_NORMAL_DEPTH', 'n_ref_count'='MATCHED_NORMAL_REF_COUNT', 'n_alt_count'='MATCHED_NORMAL_ALT_COUNT',
        'VAF'='VAF', 'Consequence'='CONSEQUENCE', 'SYMBOL'='GENE_SYMBOL', 'Transcript_ID'='ENSEMBL_TRANSCIRPT_ID',
        'Gene'='ENSEMBL_GENE_ID', 'Refseq_Feature'='REFSEQ_TRANSCRIPT_ID', 'Refseq_Gene'='ENTREZ',
        'Refseq_HGVSc'='REFSEQ_HGVSc', 'Refseq_HGVSp'='REFSEQ_HGVSp', 'HGVSc'='HGVSc', 'HGVSp'='HGVSp',
        'HGVSp_Short'='HGVSp_Short', 'HGVSg'='HGVSg', 'HGVS_OFFSET'='HGVS_OFFSET', 'EXON'='EXON', 'INTRON'='INTRON',
        'cDNA_position'='cDNA_POSITION', 'CDS_position'='CDS_POSITION', 'Amino_acids'='AMINO_ACID_CHANGE',
        'Codons'='CODON_CHANGE', 'Protein_position'='PROTEIN_POSITION', 'CANONICAL'='CANONICAL',
        'AF'='AF', 'EAS_AF'='EAS_AF', 'gnomADe_AF'='GNOMADe_AF', 'gnomADe_EAS_AF'='GNOMADe_EAS_AF',
        'MAX_AF'='MAX_AF', 'MAX_AF_POPS'='MAX_AF_POPS', 'SIFT'='SIFT', 'PolyPhen'='POLYPHEN',
        'GENE_PHENO'='GENE_PHENO', 'COSMIC'='COSMIC_ID', 'COSMIC_HGVSC'='COSMIC_HGVSc', 'COSMIC_LEGACY_ID'='COSMIC_LEGACY_ID',
        'COSMIC_TIER'='COSMIC_TIER', 'ClinVar_CLNSIG'='CLINVAR_SIGINIFICANCE', 'ClinVar_CLNDN'='CLINVAR_DISEASE',
        'ClinVar_ORIGIN'='CLINVAR_ORIGIN', 'am_class'='ALPHA_MISSENSE_CLASS', 'am_pathogenicity'='ALPHA_MISSENSE_PATHOGENICITY',
        'PureCN.PS'='PUCRECN_POSTERIOR_SOMATIC', 'PureCN.ML.C'='PURECN_CN', 'PureCN.CS'='PURECN_SUBCLONAL_CN',
        'PureCN.ML.LOH'='PURECN_LOH', 'PureCN.GHOMOZYGOUS'='PURECN_GERMLINE_HOMOZYGOUS', 'PureCN.LR'='PURECN_LOG_RATIO',
        'PureCN.CF'='PURECN_CELLULAR_FRACTION', 'IMPACT'='IMPACT', 'SYMBOL_SOURCE'='SYMBOL_SOURCE', 'HGNC_ID'='HGNC_ID',
        'BIOTYPE'='BIOTYPE', 'CCDS'='CCDS', 'SWISSPROT'='SWISSPROT', 'UNIPARC'='UNIPARC', 'DISTANCE'='DISTANCE',
        'all_effects'='ALL_EFFECTS', 'DOMAINS'='DOMAINS', 'POPAF'='POPAF', 'PON'='PON', 'TLOD'='TLOD', 'VKEY'='VKEY'
    )
    
    SelectFields <- c(RequiredFields, RenameFields)
    names(SelectFields) <- NULL
    SelectFields <- c(SelectFields[1:69], "PURECN_FLAG", SelectFields[70:length(SelectFields)])
    
    RemoveClasses <- c("Silent","Intron","5'Flank","3'Flank","Targeted_Region","IGR","3'UTR","5'UTR","RNA")
    
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

#---| REFORMAT FUNCTION |--------------------------------------------------------------------------#
    MergeMAF <- function( refseq_maf, ens_maf, RequiredFields, RenameFields, OutputMafFields, geneInfo, kovaDB )
    {
        if (!"all_effects" %in% colnames(refseq_maf)) stop("[ERROR] 'all_effects' column is missing in refseq_maf.")
        if (!"all_effects" %in% colnames(ens_maf)) stop("[ERROR] 'all_effects' column is missing in ens_maf.")

        print(1)
        refseq_maf$all_effects = sapply(refseq_maf$all_effects, function(y) { 
            if(is.na(y) || as.character(y) == "") return("")
            paste(sapply(unlist(strsplit(as.character(y),";")), function(z) tail(unlist(strsplit(z, ",")), 1)), collapse=";") 
        })
        ens_maf$all_effects = sapply(ens_maf$all_effects, function(y) { 
            if(is.na(y) || as.character(y) == "") return("")
            paste(sapply(unlist(strsplit(as.character(y),";")), function(z) tail(unlist(strsplit(z, ",")), 1)), collapse=";") 
        })
        
        refseq_maf <- refseq_maf %>% mutate(VKEY = paste0(Chromosome, ":", Start_Position, "_", Reference_Allele, ">", Tumor_Seq_Allele2))
        ens_maf <- ens_maf %>% mutate(VKEY = paste0(Chromosome, ":", Start_Position, "_", Reference_Allele, ">", Tumor_Seq_Allele2))

        rename_vector <- setNames(names(RenameFields), unname(RenameFields))
        
        reformatMAF <- refseq_maf %>% 
            left_join(ens_maf %>% select(VKEY, any_of(names(RenameFields))), by='VKEY', suffix = c("", "_ens")) %>%
            dplyr::rename(any_of(rename_vector))
        
        # [MODIFIED] HGNC_ID, ENTREZ 등 DB 정보와 충돌할 수 있는 컬럼을 사전에 제거하여 Uniqueness 오류 원천 차단
        reformatMAF <- reformatMAF %>% 
            select(-any_of(c("HGNC_ID", "HGNC_STATUS", "HGNC_SYMBOL", "HGNC_NAME", "LOCUS_GROUP", "HGNC_LOCUS_TYPE", "LOCATION", "ENTREZ", "UCSC_ID", "UNIPROT", "OMIM", "ORPHANET"))) %>%
            mutate( geneInfo[match( ens2hgnc[match(ENSEMBL_GENE_ID, ens2hgnc$ens_geneid), "hgnc_id"], geneInfo$input_id ), c(4,2,5,6,7,8,9,10,12,14,16,19)] ) %>%
            dplyr::rename( 
                HGNC_ID=hgnc_id, HGNC_STATUS=status, HGNC_SYMBOL=symbol, HGNC_NAME=name, LOCUS_GROUP=locus_group, 
                HGNC_LOCUS_TYPE=locus_type, LOCATION=location, ENTREZ=entrez, UCSC_ID=ucsc, UNIPROT=uniprot, OMIM=omim, ORPHANET=orphanet
            ) %>%
            mutate( KOVA_AF = kovaDB[match( VKEY, kovaDB$VKEY), "KOVA_AF"] ) %>%
            mutate( PURECN_FLAG = NA ) 
            
        final_maf <- reformatMAF %>% select(any_of(OutputMafFields))
        return(final_maf)
    }
#--------------------------------------------------------------------------------------------------#
    
#==================================================================================================# 
    write(" ", LOG_FILE, append=T)
    write("WES MAF Post-Processing log", LOG_FILE, append=T)
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(sprintf(" RUN DATE   : %s", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T) 
    write(sprintf(" SEQ FOLDER : %s", SEQ_FOLDER), LOG_FILE, append=T) 
    write(sprintf(" SEQ ID     : %s", SEQ_ID), LOG_FILE, append=T) 
    write("----------------------------------------------------------------------", LOG_FILE, append=T)
    write(" ", LOG_FILE, append=T)
#==================================================================================================# 

    write(sprintf("%s | Merge MAF START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Merge MAF START.")

#---| MAF FILE PREFIX |----------------------------------------------------------------------------#
    if( is_gemlineMAF ) {
        MAF_TAG = "deepvariant.germline.variant.filtered"
    } else {
        if( variantCallRunMode %in% c("tonly","Tonly") ) {
            MAF_TAG = "mutect2.bias.filtered"
        } else {
            MAF_TAG = ifelse( NORMAL_TYPE == "ORG", "mutect2.ORG.NT.bias.filtered", "mutect2.NT.bias.filtered")
        }
    }
#--------------------------------------------------------------------------------------------------#

#---| MERGE MAF |----------------------------------------------------------------------------------#
    refseq_file <- sprintf("%s/%s.%s.vep.refseq.maf",  VCF_DIR, SEQ_ID, MAF_TAG )
    ens_file    <- sprintf("%s/%s.%s.vep.ensembl.maf", VCF_DIR, SEQ_ID, MAF_TAG )
    
    if(!file.exists(refseq_file)) stop(sprintf("[ERROR] Refseq MAF file not found: %s", refseq_file))
    if(!file.exists(ens_file)) stop(sprintf("[ERROR] Ensembl MAF file not found: %s", ens_file))

    refseq_maf <- read.delim(refseq_file, skip=1, header=T)
    ens_maf    <- read.delim(ens_file, skip=1, header=T)
    
    annot_maf <- MergeMAF( refseq_maf, ens_maf, RequiredFields, RenameFields, OutputMafFields, geneInfo=ginfo, kovaDB=kova )
    
    write.table( annot_maf, sprintf("%s/%s.%s.vep.annotated.maf", VCF_DIR, SEQ_ID, MAF_TAG), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
#--------------------------------------------------------------------------------------------------#

    write(sprintf("%s | Merge MAF FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s | Variants Class Summary START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Merge MAF FINISHED.")
    message("Variants Class Summary START.")

#---| VARIANT CLASS SUMMARY |----------------------------------------------------------------------#
    mf <- read.maf(sprintf("%s/%s.%s.vep.annotated.maf", VCF_DIR, SEQ_ID, MAF_TAG))
    mf_nonsyn <- mf@data %>% filter(Variant_Classification %nin% RemoveClasses)
    mf_syn    <- mf@data %>% filter(Variant_Classification %in% RemoveClasses)
    
    VariantClassSummary <- rbind(
        mf_nonsyn %>% 
            group_by(Variant_Classification, Variant_Type) %>% 
            summarise( stats=n(), .groups = 'drop' ) %>%
            mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID, variant_group='non_synonymous' ),
        mf_syn %>% 
            group_by(Variant_Classification, Variant_Type) %>% 
            summarise( stats=n(), .groups = 'drop' ) %>%
            mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID, variant_group='synonymous' )
    )
#--------------------------------------------------------------------------------------------------#

    write(sprintf("%s | Variants Class Summary FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Variants Class Summary FINISHED.")

#---| SOMATIC VARIANTS ONLY RESULTS |--------------------------------------------------------------#
    if( writeSomaticVarinatsOnlyMAF ) {
        write(sprintf("%s | SomaticVariantsOnly MAF Generate START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message("writeSomaticVarinatsOnlyMAF = TRUE, SomaticVariantsOnly MAF Generate START.")

        if( is_gemlineMAF ) {
            write.table( mf_nonsyn, sprintf("%s/%s.%s.non_synonymous.maf", VCF_DIR, SEQ_ID, MAF_TAG), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
            write.table( mf_syn, sprintf("%s/%s.%s.synonymous.maf", VCF_DIR, SEQ_ID, MAF_TAG), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
            write.table( VariantClassSummary, sprintf("%s/%s.%s.maf.variant.class.summary.txt", VCF_DIR, SEQ_ID, MAF_TAG), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
        } else {
            mf_nonsyn_filtered <- mf@data %>% 
                filter( FILTER=="PASS" ) %>%
                filter( TUMOR_DEPTH >= CutOff_DP, TUMOR_ALT_COUNT >= CutOff_ALT_DP ) %>%
                filter( is.na(AF)             | AF             < CutOff_POP_AF ) %>%
                filter( is.na(EAS_AF)         | EAS_AF         < CutOff_POP_AF ) %>%
                filter( is.na(GNOMADe_AF)     | GNOMADe_AF     < CutOff_POP_AF ) %>%
                filter( is.na(GNOMADe_EAS_AF) | GNOMADe_EAS_AF < CutOff_POP_AF ) %>%
                mutate( VAF = round(TUMOR_ALT_COUNT/TUMOR_DEPTH, 3) )
            
            mf_nonsyn_filtered <- unique(
                rbind(
                    mf_nonsyn_filtered %>% filter(Variant_Type %nin% c("INS","DEL"), VAF >= 0.02 ),  
                    mf_nonsyn_filtered %>% filter(Variant_Type %in%  c("INS","DEL"), VAF >= 0.05 ),
                    mf@data %>% filter( FILTER=="PASS", GENE_SYMBOL %in% essGenes ) %>% mutate(VAF = round(TUMOR_ALT_COUNT/TUMOR_DEPTH, 3) )
                )
            ) %>% 
            filter( GENE_SYMBOL %nin% blacklistGenes ) %>%
            filter( Variant_Classification %nin% RemoveClasses )
            
            VCS <- mf_nonsyn_filtered %>% 
                group_by(Variant_Classification, Variant_Type) %>% 
                summarise( stats=n(), .groups = 'drop' ) %>%
                mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID, variant_group='filtered_somatic_only' )
            
            VariantClassSummary <- rbind(VariantClassSummary, VCS)
            
            write.table( mf_nonsyn_filtered, sprintf("%s/%s.%s.somatic.variants.only.maf", VCF_DIR, SEQ_ID, MAF_TAG), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
            write.table( VariantClassSummary, sprintf("%s/%s.%s.maf.variant.class.summary.txt", VCF_DIR, SEQ_ID, MAF_TAG), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
            
            write(sprintf("%s | SomaticVariantsOnly MAF Generate FINISHED", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("SomaticVariantsOnly MAF Generate FINISHED")
        }
    } else {
        write.table( VariantClassSummary, sprintf("%s/%s.%s.maf.variant.class.summary.txt", VCF_DIR, SEQ_ID, MAF_TAG), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t" )
        write(sprintf("%s | SomaticVariantsOnly MAF Generation SKIPPED", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message("SomaticVariantsOnly MAF Generation SKIPPED")
    }
#--------------------------------------------------------------------------------------------------#

#---| IMPORT INTO DATABSE |------------------------------------------------------------------------#
if( writeDB ) {
    write(sprintf("%s | RESULTS IMPORT INTO DB START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("RESULTS IMPORT INTO DB START.")

    if( MAF_TAG == "mutect2.bias.filtered"        ){ VAR_CALL_MODE = "Tonly" }
    if( MAF_TAG == "mutect2.NT.bias.filtered"     ){ VAR_CALL_MODE = "TS_N"  }
    if( MAF_TAG == "mutect2.ORG.NT.bias.filtered" ){ VAR_CALL_MODE = "ORG_N" }

    dbCon <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
    preDataClear <- dbGetQuery(dbCon, sprintf("DELETE FROM variants_summary WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", SEQ_FOLDER, SEQ_ID, VAR_CALL_MODE))
    writeData <- dbWriteTable(dbCon, name="variants_summary", value=data.frame(variant_call_mode=VAR_CALL_MODE, VariantClassSummary), row.names=FALSE, append=TRUE)
    
    if( writeSomaticVarinatsOnlyMAF ) {
        preDataClear <- dbGetQuery(dbCon, sprintf("DELETE FROM variants_somatic WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", SEQ_FOLDER, SEQ_ID, VAR_CALL_MODE))
        target_somatic_data <- if(is_gemlineMAF) mf_nonsyn else mf_nonsyn_filtered
        writeData <- dbWriteTable(dbCon, name="variants_somatic", value=data.frame(variant_call_mode=VAR_CALL_MODE, target_somatic_data), row.names=FALSE, append=TRUE)
        write(sprintf("%s | filtered somatic variants IMPORTED INTO DB", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    }
    
    dbDisconnect(dbCon)
    write(sprintf("%s | RESULTS IMPORT INTO DB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("RESULTS IMPORT INTO DB FINISHED.")
} else {
    write(sprintf("%s | RESULTS IMPORT INTO DB SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("RESULTS IMPORT INTO DB SKIPPED.")
}
#--------------------------------------------------------------------------------------------------#

    write("   ", LOG_FILE, append=T)
    write(sprintf("%s | WES MAF Post-Processing DONE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE )