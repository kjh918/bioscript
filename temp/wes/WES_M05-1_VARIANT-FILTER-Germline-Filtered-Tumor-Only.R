
#---| PACKAGES |-----------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)   
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("GenomicRanges"))
    suppressPackageStartupMessages(library("plyranges"))
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list( 
        make_option(c("--BASE_DIR"), action="store", default=NULL, type="character", help="analysis base folder"),
        make_option(c("--BATCH_ID"), action="store", default=NULL, type="character", help="analysis batch-id"),
        make_option(c("--TONLY_MAF_SEQ_ID"), action="store", default=NULL, type="character", help="input tumor-only call variants maf sample id"),
        make_option(c("--GERMLINE_MAF_SEQ_IDS"), action="store", default=NULL, type="character", help="deepvariant germline variants maf sample ids list"),   
        make_option(c("--DB_IMPORT"), action="store", default=FALSE, type="logical", help="database import")
    )
    #--------------------------------------------------------------------------#  
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR             <- ARGS$BASE_DIR
    BATCH_ID             <- ARGS$BATCH_ID
    TONLY_MAF_SEQ_ID     <- ARGS$TONLY_MAF_SEQ_ID
    GERMLINE_MAF_SEQ_IDS <- ARGS$GERMLINE_MAF_SEQ_IDS
    DB_IMPORT            <- ARGS$DB_IMPORT
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
    # BASE_DIR = "/data/wes"
    # BATCH_ID = "WES_25_05"
    # TONLY_MAF_SEQ_ID = "WES_25_05_01"
    # GERMLINE_MAF_SEQ_IDS = "WES_25_05_01"
    # DB_IMPORT=TRUE
#--------------------------------------------------------------------------------------------------#

    Filter.TSO.Blacklist = function( InputMAF, TSO_BLACKLIST_GR )
    {
        MAF_GR = GRanges(
            InputMAF$Chromosome,
            IRanges( start=InputMAF$Start_Position, end=InputMAF$End_Position ),
            vkey = InputMAF$vkey
        )
        RmTargets = findOverlaps(MAF_GR, TSO_BLACKLIST_GR)
        if( length(RmTargets) > 0 )
        {
            RmVkeys = unique(MAF_GR[queryHits(RmTargets)]$"vkey")
            tso_filtered_maf = InputMAF %>% filter( vkey %nin% RmVkeys )
        }else{
            tso_filtered_maf = InputMAF
        }
        return(tso_filtered_maf)
    }


    load("/storage/home/kangsm/myDB/rds/hg19_gcx.MRD.blacklists_gr.Rdata") # TSO_BLACKLIST_GR,ENCODE_BLACKLIST_GR,MAPPABILITY_EXCLUDE_GR,DNASE_HYPERSENSE_GR,FAIRESEQ_PEAKS_GR



#---| FILTER VARINATS : GERMLINE VARIANTS |--------------------------------------------------------#
    if(is.null(TONLY_MAF_SEQ_ID)){ stop(">> No Input-Tonly-MAF. REQUIRED.") }
    if(is.null(GERMLINE_MAF_SEQ_IDS)){ stop(">> No Germline Variant Sample-IDs. REQUIRED.") }
    message(">> create common germline variants index...")
    GermlineVars_SeqIDs  = unlist(strsplit(GERMLINE_MAF_SEQ_IDS, ","))
    GermlineVariantIndex = Reduce( 
        intersect,
        lapply( GermlineVars_SeqIDs, function(sid) 
        {
            GV = read.delim(sprintf("%s/%s/%s/vcf/%s.deepvariant.germline.variant.filtered.vep.annotated.maf", BASE_DIR, BATCH_ID, sid, sid))
            GV = GV %>% mutate( vkey=paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep=":") )
            return(GV$vkey)
        })
    )
    message(">> read tumor-only somatic variants...")
    Tonly_Maf_File = sprintf("%s/%s/%s/vcf/%s.mutect2.bias.filtered.somatic.variants.only.maf", BASE_DIR, BATCH_ID, TONLY_MAF_SEQ_ID, TONLY_MAF_SEQ_ID)
    TonlyVars = read.delim(Tonly_Maf_File) %>% 
        mutate( vkey=paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep=":") )
    
    message(">> filter germline variants...")
    Tonly_GV_Filtered_Vars = TonlyVars %>% filter( vkey %nin% GermlineVariantIndex )
#--------------------------------------------------------------------------------------------------#

#---| FILTER VARINATS : TSO Blacklist |------------------------------------------------------------#    
    message(">> filter TSO500 balcklist...")
    TSO_Tonly_GV_Filtered_Vars = Filter.TSO.Blacklist(InputMAF=Tonly_GV_Filtered_Vars, TSO_BLACKLIST_GR=TSO_BLACKLIST_GR)
#--------------------------------------------------------------------------------------------------#

#---| SAVE NEW MAF |-------------------------------------------------------------------------------#    
    message(">> save germline filtered tumor-only called somatic variants as MAF file...")
    output_maf_name = sprintf("%s/%s/%s/vcf/%s.mutect2.bias.filtered.germline.filtered.somatic.variants.only.maf", BASE_DIR, BATCH_ID, TONLY_MAF_SEQ_ID, TONLY_MAF_SEQ_ID )  
    write.table( Tonly_GV_Filtered_Vars, output_maf_name, quote=F, col.names=T, row.names=F, sep="\t" )

    message(">> create variant-class summary table...")
    VariantClassSummary = Tonly_GV_Filtered_Vars %>% 
        group_by(Variant_Classification, Variant_Type, VARIANT_CLASS) %>% 
        summarise( stats=length(vkey) ) %>%
        mutate( seq_folder=BATCH_ID, seq_id=TONLY_MAF_SEQ_ID, variant_group='filtered_somatic_only' ) 
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------# 
    if( DB_IMPORT )
    {
        message(">> import results into database...")
        source("/data/wes/params/ruo_wes_db.R")
        VAR_CALL_MODE = "TonlyGermlineFiltered"
        #----------------------------------------------------------------------#
        dbCon        <- dbConnect(dbDriver("MySQL"), host=db_host, user=db_user, port=db_port, password=db_pw, db=db_name_data )
        preDataClear1 <- dbGetQuery(dbCon, 
            sprintf("DELETE FROM variants_summary WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", BATCH_ID, TONLY_MAF_SEQ_ID, VAR_CALL_MODE)
        )
        writeData1 <- dbWriteTable(dbCon, name="variants_summary", value=data.frame(variant_call_mode=VAR_CALL_MODE, VariantClassSummary), row.names=FALSE, append=TRUE)
        
        preDataClear2 <- dbGetQuery(dbCon, 
            sprintf("DELETE FROM variants_somatic WHERE seq_folder = '%s' AND seq_id = '%s' AND variant_call_mode = '%s'", BATCH_ID, TONLY_MAF_SEQ_ID, VAR_CALL_MODE)
        )
        writeData2 <- dbWriteTable(dbCon, name="variants_somatic", value=data.frame(variant_call_mode=VAR_CALL_MODE, Tonly_GV_Filtered_Vars), row.names=FALSE, append=TRUE)
        dbDisconnect(dbCon)
    }
#--------------------------------------------------------------------------------------------------#
    message(">> FINISHED.")

  

