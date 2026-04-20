

#   script version v1.1
#   update 2023-12-08
#   summary of change : merge default-vep, refseq-vep (remove gatk-funcotator)

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
        make_option(c("--seqFolder"), action="store", default=NA, type="character", help="$SEQ_FOLDER"),
        make_option(c("--seqID")    , action="store", default=NA, type="character", help="$SEQ_ID"),
        make_option(c("--baseDir")  , action="store", default=NA, type="character", help="$BASE_DIR"),
        make_option(c("--nType"), action="store", default="TS", type="character", help="type of normal sample"),
        make_option(c("--writeSomaticVarinatsOnlyMAF") , action="store", default=TRUE, type="logical", 
            help="write another maf file including only filtered somatic variants"),
        make_option(c("--DP"), action="store", default=10, type="double", help="variant filter: DEPTH "),
        make_option(c("--ALT"), action="store", default=2, type="double", help="variant filter: ALT DEPTH "),
        make_option(c("--POPAF"), action="store", default=0.01, type="double", help="variant filter: POP AF "),
        make_option(c("--writeDB") , action="store", default=TRUE, type="logical", 
            help="write filtered variants maf into database"),
        make_option(c("--variantCallMode"), action="store", default="tonly", type="character", help="variant calling run mode"),
        make_option(c("--GERMLINE"), action="store", default=FALSE, type="logical", help="Germline variant MAF or not. default = FALSE ( = somatic )")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
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
    dbDisconnect(dbCon)
#--------------------------------------------------------------------------------------------------#

#---| FOLDERS |------------------------------------------------------------------------------------#
    DATA_DIR <- sprintf("%s/%s", BASE_DIR,  SEQ_FOLDER)
    VCF_DIR  <- sprintf("%s/%s/vcf", DATA_DIR, SEQ_ID)
    LOG_DIR  <- sprintf("%s/%s/log", DATA_DIR, SEQ_ID)
#--------------------------------------------------------------------------------------------------#

#---| LOG FILE |-----------------------------------------------------------------------------------#
    LOG_FILE <- sprintf("%s/%s.maf.post-processing.%s.log", LOG_DIR, SEQ_ID, today())
    if( ! file.exists(LOG_FILE) ){ system(sprintf("touch %s", LOG_FILE)) }
#--------------------------------------------------------------------------------------------------#

#---| MAF FILEDS |---------------------------------------------------------------------------------#
    requireFields = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', "Reference_Allele", 
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode'
    )
    cols_refseqMaf = c("NCBI_Build","Strand","Transcript_ID",
        "Gene","Exon_Number","Consequence","ENSP","VARIANT_CLASS","CANONICAL","cDNA_position",
        "CDS_position","Protein_position","Amino_acids","Codons","t_depth","t_ref_count","t_alt_count",
        "n_depth","n_ref_count","n_alt_count",'DP','AF','AFR_AF','AMR_AF','ASN_AF','EAS_AF','EUR_AF',
        'SAS_AF','AA_AF','EA_AF',"BIOTYPE","dbSNP_RS",'SIFT','PolyPhen','HGVS_OFFSET','PHENO','am_class',
        'am_pathogenicity',"CADD_PHRED","CADD_RAW","ClinVar_CLNSIG","ClinVar_CLNDN","ClinVar_ORIGIN",
        "gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF",
        "gnomADe_NFE_AF","gnomADe_OTH_AF","gnomADe_SAS_AF","MAX_AF","MAX_AF_POPS","TRANSCRIPTION_FACTORS",
        "COSMIC","COSMIC_HGVSC","COSMIC_LEGACY_ID","COSMIC_TIER","GERMQ","MBQ","MFRL","MMQ","MPOS",
        "POPAF","PON","ROQ","TLOD","PUBMED",'flanking_bps','FILTER','pArtifact','artiStatus','vcf_pos','vkey'
    )    
    cols_ensMaf = c("Transcript_ID","Gene","HGVSc","HGVSp","HGVSp_Short","HGVSg","ENSP","cDNA_position","Exon_Number","DISTANCE",
        "CANONICAL","CCDS","UNIPARC","DOMAINS","vkey"
    )
#--------------------------------------------------------------------------------------------------#

#---| REFORMAT FUNCTION |--------------------------------------------------------------------------#
#' @description combine refseq-MAF and ensembl-MAF to integrated-MAF
#' @param refseq_maf     RefSeq annotated maf table1
#' @param ens_maf        Ensembl annotated maf table
#' @param requireFields  Required Common Columns
#' @param cols_refseqMaf RefSeq MAF Columns
#' @param cols_ensMaf    Ensembl MAF Columns 
#' @export
    MergeMAF <- function( refseq_maf, ens_maf, requireFields, cols_refseqMaf, cols_ensMaf, geneInfo )
    {
        refseq_maf$REFSEQ_ALL_TRANSCRIPTS = sapply(refseq_maf$all_effects, function(y)
        { paste(sapply(unlist(strsplit(y,";")), function(z) tail(unlist(strsplit(z, ",")), 1)),collapse=";") })
        ###
        MAF1 = refseq_maf %>% select(any_of(c(requireFields,cols_refseqMaf,"REFSEQ_ALL_TRANSCRIPTS"))) %>% 
            dplyr::rename( 
                REF_GEOME = NCBI_Build, 
                Refseq_TRANSCRIPT=Transcript_ID, 
                Refseq_GENE=Gene, 
                Refseq_EXON=Exon_Number, 
                Refseq_PROTEIN=ENSP, 
                Refseq_CANONICAL=CANONICAL, 
                Refseq_cDNA_POS=cDNA_position,
                CDS_POS=CDS_position, 
                PROTEIN_POS=Protein_position, 
                AMINO_ACIDS=Amino_acids, 
                CODONS=Codons, 
                TUMOR_DEPTH=t_depth, 
                TUMOR_REF_COUNT=t_ref_count, 
                TUMOR_ALT_COUNT=t_alt_count,
                NORMAL_DEPTH=n_depth, 
                NORMAL_REF_COUNT=n_ref_count, 
                NORMAL_ALT_COUNT=n_alt_count, 
                DEPTH=DP, 
                G1K_AF=AF, 
                G1K_AFR_AF=AFR_AF, 
                G1K_AMR_AF=AMR_AF, 
                G1K_ASN_AF=ASN_AF, 
                G1K_EAS_AF=EAS_AF, 
                G1K_EUR_AF=EUR_AF, 
                G1K_SAS_AF=SAS_AF, 
                G1K_AA_AF=AA_AF, 
                G1K_EA_AF=EA_AF, 
                POLYPHEN=PolyPhen,
                ALPHA_MISSENSE_CLASS=am_class, 
                ALPHA_MISSENSE_PATHGENICITY=am_pathogenicity, 
                CLINVAR_SIG=ClinVar_CLNSIG, 
                CLINVAR_DISEASE=ClinVar_CLNDN, 
                CLINVAR_ORIGIN=ClinVar_ORIGIN,
                COSMIC_HGVSc=COSMIC_HGVSC, 
                COSMIC_CANCER_GENE_CENSUS_TIER=COSMIC_TIER,
                FLANKING_BPS=flanking_bps, 
                VCF_POS=vcf_pos,
                BIAS_PVAL=pArtifact,
                BIAS_STATUS=artiStatus
            )
        
        ens_maf$ENS_ALL_TRANSCRIPTS = sapply(ens_maf$all_effects, function(y)
        { paste(sapply(unlist(strsplit(y,";")), function(z) tail(unlist(strsplit(z, ",")), 1)),collapse=";") })
        ###
        MAF2 = ens_maf %>% select(any_of(c(cols_ensMaf,"ENS_ALL_TRANSCRIPTS"))) %>% 
            dplyr::rename( 
                Ensembl_TRANSCRIPT=Transcript_ID, 
                Ensembl_GENE=Gene, 
                Ensembl_PROTEIN=ENSP, 
                Ensembl_cDNA_POS=cDNA_position, 
                Ensembl_EXON=Exon_Number, 
                Ensembl_CANONICAL=CANONICAL 
            )
        ###
        reformatMAF <- MAF1 %>% left_join(., MAF2, by='vkey') 
        ###
        reformatMAF = reformatMAF %>% 
            mutate( geneInfo[match( ens2hgnc[match(Ensembl_GENE, ens2hgnc$ens_geneid), "hgnc_id"], geneInfo$input_id ), c(4,2,5,6,7,8,9,10,12,14,16,19)] ) %>%
            dplyr::rename( 
                HGNC_ID=hgnc_id, 
                HGNC_STATUS=status, 
                HGNC_SYMBOL=symbol, 
                HGNC_NAME=name, 
                LOCUS_GROUP=locus_group, 
                LOCUS_TYPE=locus_type,
                CHR_LOCATION=location, 
                ENTREZ=entrez, 
                UCSC_ID=ucsc, 
                UNIPROT=uniprot, 
                OMIM=omim, 
                ORPHANET=orphanet
            ) %>%
            mutate( seq_folder = SEQ_FOLDER, seq_id = SEQ_ID ) %>%
            mutate( KOVA_AF = kova[match( vkey, kova$vkey), "AF"] )
        ### 
        return(reformatMAF)
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

#==================================================================================================# 
    write(sprintf("%s | Merge MAF START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Merge MAF START.")
#==================================================================================================# 

#---| MAF FILE PREFIX |----------------------------------------------------------------------------#
    if( is_gemlineMAF )
    {
        MAF_TAG = "deepvariant.germline.variant.filtered"
    }else{
        if( variantCallRunMode %in% c("tonly","Tonly") )
        {
            MAF_TAG = "mutect2.bias.filtered"
        }else{
            MAF_TAG = ifelse( NORMAL_TYPE == "ORG", "mutect2.ORG.NT.bias.filtered", "mutect2.NT.bias.filtered")
        }
    }
#--------------------------------------------------------------------------------------------------#

#---| MERGE MAF |----------------------------------------------------------------------------------#
    ##
    refseq_maf <- read.delim(sprintf("%s/%s.%s.vep.refseq.maf",  VCF_DIR, SEQ_ID, MAF_TAG ), skip=1, header=T) %>%
        mutate(vkey=paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep=":"))
    ens_maf    <- read.delim(sprintf("%s/%s.%s.vep.ensembl.maf", VCF_DIR, SEQ_ID, MAF_TAG ), skip=1, header=T) %>%
        mutate(vkey=paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep=":"))
    ##
    annot_maf <- MergeMAF( refseq_maf, ens_maf, requireFields, cols_refseqMaf, cols_ensMaf, geneInfo=ginfo )
    ##
    write.table( annot_maf, 
        sprintf("%s/%s.%s.vep.annotated.maf", VCF_DIR, SEQ_ID, MAF_TAG),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | Merge MAF FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write(sprintf("%s | Variants Class Summary START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Merge MAF FINISHED.")
    message("Variants Class Summary START.")
#==================================================================================================# 

#---| VARIANT CLASS SUMMARY |----------------------------------------------------------------------#
    mf <- read.maf(sprintf("%s/%s.%s.vep.annotated.maf", VCF_DIR, SEQ_ID, MAF_TAG))
    mf_nonsyn <- mf@data 
    mf_syn    <- mf@maf.silent 
    ##
    VariantClassSummary <- rbind(
        mf_nonsyn %>% 
            group_by(Variant_Classification, Variant_Type, VARIANT_CLASS) %>% 
            summarise( stats=length(vkey) ) %>%
            mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID, variant_group='non_synonymous' ),
        mf_syn %>% 
            group_by(Variant_Classification, Variant_Type, VARIANT_CLASS) %>% 
            summarise( stats=length(vkey) ) %>%
            mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID, variant_group='synonymous' )
    )
#--------------------------------------------------------------------------------------------------#

#==================================================================================================# 
    write(sprintf("%s | Variants Class Summary FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    message("Variants Class Summary FINISHED.")
#==================================================================================================# 

#---| SOMATIC VARIANTS ONLY RESULTS |--------------------------------------------------------------#
    if( writeSomaticVarinatsOnlyMAF )
    {
        #==========================================================================================# 
            write(sprintf("%s | SomaticVariantsOnly MAF Generate START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("writeSomaticVarinatsOnlyMAF = TRUE, SomaticVariantsOnly MAF Generate START.")
        #==========================================================================================# 

        #----------------------------------------------------------------------#
        dbCon      <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
        addon_data <- dbGetQuery(dbCon, "SELECT * FROM custom_addon_data")
        dbDisconnect(dbCon)
        dbCon          <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_resource )
        essentialGenes <- dbGetQuery(dbCon, "SELECT * FROM ngs_panel_essential_genes WHERE esgd_id = '1'")
        dbDisconnect(dbCon)
        #----------------------------------------------------------------------#
        kova           <- kova %>% mutate(vkey = paste( Chr,Pos,Ref,Alt, sep=":") )
        essGenes       <- essentialGenes$hgnc_symbol
        blacklistGenes <- unlist(strsplit(addon_data[which(addon_data$tag == "hla_genes"), "data_value"],";"))
        #----------------------------------------------------------------------#
        if( is_gemlineMAF )
        {
            mf <- read.maf(sprintf("%s/%s.%s.vep.annotated.maf", VCF_DIR, SEQ_ID, MAF_TAG))
            mf_nonsyn <- mf@data
            mf_syn    <- mf@maf.silent
            #------------------------------------------------------------------#
            write.table( mf_nonsyn, 
                sprintf("%s/%s.%s.non_synonymous.maf", VCF_DIR, SEQ_ID, MAF_TAG),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
            )
            #------------------------------------------------------------------#
            write.table( mf_syn, 
                sprintf("%s/%s.%s.synonymous.maf", VCF_DIR, SEQ_ID, MAF_TAG),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
            )
            #------------------------------------------------------------------#
            write.table( VariantClassSummary, 
                sprintf("%s/%s.%s.maf.variant.class.summary.txt", VCF_DIR, SEQ_ID, MAF_TAG),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
            )

        }else{
            mf <- read.maf(sprintf("%s/%s.%s.vep.annotated.maf", VCF_DIR, SEQ_ID, MAF_TAG))
            ##
            mf_nonsyn <- mf@data %>% 
                filter( FILTER=="PASS" , LOCUS_GROUP == "protein-coding gene" ) %>%
                filter( DEPTH >= CutOff_DP, TUMOR_ALT_COUNT >= CutOff_ALT_DP ) %>%
                filter( is.na(G1K_AF)         | G1K_AF     < CutOff_POP_AF ) %>%
                filter( is.na(G1K_EAS_AF)     | G1K_EAS_AF < CutOff_POP_AF ) %>%
                filter( is.na(G1K_EAS_AF)     | G1K_EAS_AF < CutOff_POP_AF ) %>%
                filter( is.na(gnomADe_AF)     | gnomADe_AF < CutOff_POP_AF ) %>%
                filter( is.na(gnomADe_EAS_AF) | gnomADe_EAS_AF < CutOff_POP_AF ) %>%
                mutate(VAF = round(TUMOR_ALT_COUNT/DEPTH, 3) )
            ##
            mf_nonsyn <- unique(
                rbind(
                    mf_nonsyn %>% filter(Variant_Type %nin% c("INS","DEL"), VAF >= 0.02 ),  
                    mf_nonsyn %>% filter(Variant_Type %in%  c("INS","DEL"), VAF >= 0.05 ),
                    mf@data %>% filter( FILTER=="PASS", HGNC_SYMBOL %in% essGenes ) %>% mutate(VAF = round(TUMOR_ALT_COUNT/DEPTH, 3) )
                )
            ) %>% 
                filter( HGNC_SYMBOL %nin% blacklistGenes )    
            #------------------------------------------------------------------#
            VCS <- mf_nonsyn %>% 
                group_by(Variant_Classification, Variant_Type, VARIANT_CLASS) %>% 
                summarise( stats=length(vkey) ) %>%
                mutate( seq_folder=SEQ_FOLDER, seq_id=SEQ_ID, variant_group='filtered_somatic_only' )
            ##
            VariantClassSummary <- rbind(VariantClassSummary, VCS)
            #------------------------------------------------------------------#
            write.table( mf_nonsyn, 
                sprintf("%s/%s.%s.somatic.variants.only.maf", VCF_DIR, SEQ_ID, MAF_TAG),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
            )
            write.table( VariantClassSummary, 
                sprintf("%s/%s.%s.maf.variant.class.summary.txt", VCF_DIR, SEQ_ID, MAF_TAG),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
            )
            #======================================================================================# 
                write(sprintf("%s | SomaticVariantsOnly MAF Generate FINISHED", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
                write(sprintf("%s | SomaticVariantsOnly MAF File SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
                write(sprintf("%s | VariantClassSummary File SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
                message("SomaticVariantsOnly MAF Generate FINISHED")
                message("SomaticVariantsOnly MAF & VariantClassSummary Table Saved as File.")
            #======================================================================================# 
        }
    }else{
        write.table( VariantClassSummary, 
            sprintf("%s/%s.%s.maf.variant.class.summary.txt", VCF_DIR, SEQ_ID, MAF_TAG),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t"
        )
        #==========================================================================================# 
            write(sprintf("%s | SomaticVariantsOnly MAF Generation SKIPPED", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            write(sprintf("%s | VariantClassSummary File SAVED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            message("SomaticVariantsOnly MAF Generation SKIPPED")
            message("Only VariantClassSummary Table Saved as File.")
        #==========================================================================================# 
    }
#--------------------------------------------------------------------------------------------------#

#---| IMPORT INTO DATABSE |------------------------------------------------------------------------#
if( writeDB )
{
    #==============================================================================================# 
        write(sprintf("%s | RESULTS IMPORT INTO DB START.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message("RESULTS IMPORT INTO DB START.")
    #==============================================================================================# 

    if( MAF_TAG == "mutect2.bias.filtered"        ){ VAR_CALL_MODE = "Tonly" }
    if( MAF_TAG == "mutect2.NT.bias.filtered"     ){ VAR_CALL_MODE = "TS_N"  }
    if( MAF_TAG == "mutect2.ORG.NT.bias.filtered" ){ VAR_CALL_MODE = "ORG_N" }


    dbCon        <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
    preDataClear <- dbGetQuery(dbCon, 
        sprintf("DELETE FROM variants_summary WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", SEQ_FOLDER, SEQ_ID, VAR_CALL_MODE)
    )
    writeData <- dbWriteTable(dbCon, name="variants_summary", value=data.frame(variant_call_mode=VAR_CALL_MODE, VariantClassSummary), row.names=FALSE, append=TRUE)
    #--------------------------------------------------------------------------#
    if( writeSomaticVarinatsOnlyMAF )
    {
        preDataClear <- dbGetQuery(dbCon, 
            sprintf("DELETE FROM variants_somatic WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", SEQ_FOLDER, SEQ_ID, VAR_CALL_MODE)
        )
        writeData <- dbWriteTable(dbCon, name="variants_somatic", value=data.frame(variant_call_mode=VAR_CALL_MODE, mf_nonsyn), row.names=FALSE, append=TRUE)

        #==========================================================================================# 
            write(sprintf("%s | filtered somatic variants IMPORTED INTO DB", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
            write(sprintf("%s | variant class summary IMPORTED INTO DB", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        #==========================================================================================#
    }else{
        #==========================================================================================# 
            write(sprintf("%s | variant class summary IMPORTED INTO DB", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        #==========================================================================================#
    }
    #--------------------------------------------------------------------------#
    dbDisconnect(dbCon)
    ##
    #==============================================================================================# 
        write(sprintf("%s | RESULTS IMPORT INTO DB FINISHED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message("RESULTS IMPORT INTO DB FINISHED.")
    #==============================================================================================#
}else{
    #==============================================================================================# 
        write(sprintf("%s | RESULTS IMPORT INTO DB SKIPPED.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
        message("RESULTS IMPORT INTO DB SKIPPED.")
    #==============================================================================================#
}
#--------------------------------------------------------------------------------------------------#

#==================================================================================================#
    write("   ", LOG_FILE, append=T)
    write(sprintf("%s | WES MAF Post-Processing DONE.", format(now() ,format = "%Y-%m-%d %H:%M:%S")), LOG_FILE, append=T)
    write("----------------------------------------------------------------------", file = LOG_FILE, append = TRUE ) 
    write("   ", LOG_FILE, append=T)
#==================================================================================================# 
