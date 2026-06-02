
# SCRIPT VERSION : v1.0
# GUIDELINE ID   : GMC-NGS-B02-A03
# DATE           : 2024-03-11
# AUTHOR         : kangsm@gencurix.com

#---| PACKAGES |-----------------------------------------------------------------------------------#
    options(stringsAsFactors=FALSE)
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("lubridate"))
    suppressPackageStartupMessages(library("parallel"))
#--------------------------------------------------------------------------------------------------#

#---| DEFAULT VALUE |------------------------------------------------------------------------------#
    WES_RunGermlineVariant = FALSE
#--------------------------------------------------------------------------------------------------#

#---| ARGUMENTS |----------------------------------------------------------------------------------#
    option_list = list(
        make_option(c("--BASE_DIR"        ), action="store", default=NULL,  type="character", help="ngs base folder"),
        make_option(c("--BATCH_ID"        ), action="store", default=NULL,  type="character", help="batch id ( = seq folder)"),
        make_option(c("--CLIENT_ID"       ), action="store", default=NULL,  type="character", help="client id"),
        make_option(c("--SAMPLE_INFO_FILE"), action="store", default=NULL,  type="character", help="sample info file"),
        make_option(c("--ORDER_ID"        ), action="store", default=NULL,  type="character", help="order id"),
        make_option(c("--SEQ_TYPE"        ), action="store", default=NULL,  type="character", help="sequencing sample type. default = DNA"),
        make_option(c("--THREADS"         ), action="store", default=NULL,  type="character", help="threads for sge scripts"),
        make_option(c("--DP_RD_FROM_FASTQ"), action="store", default=NULL,  type="logical",   help="get depth and read-length info from fastq file"),
        make_option(c("--DP"              ), action="store", default=NULL,  type="character", help="depth. required when DP_RD_FROM_FASTQ = FALSE"),
        make_option(c("--RD"              ), action="store", default=NULL,  type="character", help="read-length. required when DP_RD_FROM_FASTQ = FALSE"),
        make_option(c("--WES_NGS_LIB"     ), action="store", default=NULL,  type="character", help="NGS_LIBRARY for sge scripts"),
        make_option(c("--WES_RunGermlineVariant"), action="store", default=FALSE, type="logical",   help="run germline variant (processing)"),
        make_option(c("--WTS_Stranded"          ), action="store", default=FALSE, type="logical",   help="WTS standed library or not. default = FALSE "),
        make_option(c("--WTS_PairedSeq"         ), action="store", default=TRUE,  type="logical",   help="WTS seq. type. default = 'PE'"),
        make_option(c("--WTS_ReferenceGenome"   ), action="store", default=NULL,  type="character", help="WTS reference genome. default = 'hg19'")
    )
    #--------------------------------------------------------------------------#
    ARGS = parse_args(OptionParser(option_list=option_list))
    #--------------------------------------------------------------------------#
    BASE_DIR         = ARGS$BASE_DIR 
    BATCH_ID         = ARGS$BATCH_ID
    CLIENT_ID        = ARGS$CLIENT_ID
    SAMPLE_INFO_FILE = ARGS$SAMPLE_INFO_FILE
    ORDER_ID         = ARGS$ORDER_ID
    SEQ_TYPE         = ARGS$SEQ_TYPE 
    THREADS          = ARGS$THREADS 
    DP_RD_FROM_FASTQ = ARGS$DP_RD_FROM_FASTQ
    DP               = ARGS$WES_DP 
    RD               = ARGS$WES_RD 
    NGS_LIB          = ARGS$WES_NGS_LIB 
    WES_RunGermlineVariant = ARGS$WES_RunGermlineVariant  
    WTS_Stranded        = ARGS$WTS_Stranded  
    WTS_PairedSeq       = ARGS$WTS_PairedSeq  
    WTS_ReferenceGenome = ARGS$WTS_ReferenceGenome  
#--------------------------------------------------------------------------------------------------#
#---| MANUAL PARAMS |------------------------------------------------------------------------------#
    # BASE_DIR         = "/data/wts"
    # BATCH_ID         = "WTS_25_02" 
    # CLIENT_ID        = 11
    # SAMPLE_INFO_FILE = "WTS_25_02_seqid_info.xlsx"
    # ORDER_ID         = "GCX-C11-251215"
    # SEQ_TYPE         = "RNA" 
    # WES_RunGermlineVariant = TRUE  
    # THREADS          = 15
    # NGS_LIB          = "twist.exome.2.0"
    # DP_RD_FROM_FASTQ = TRUE
    # DP               = NA
    # RD               = 150
    # WTS_Stranded     = FALSE
    # WTS_PairedSeq    = TRUE
    # WTS_ReferenceGenome = "hg19"
#--------------------------------------------------------------------------------------------------#   
#---| CHECK PARAMS |-------------------------------------------------------------------------------#
    if(is.null(BASE_DIR) ){ stop("No 'BASE_DIR' found. REQURIED.")  }
    if(is.null(BATCH_ID) ){ stop("No 'BATCH_ID' found. REQURIED.")  }
    if(is.null(CLIENT_ID)){ stop("No 'CLIENT_ID' found. REQURIED.") }
    if(is.null(SAMPLE_INFO_FILE)){ stop("No 'SAMPLE_INFO_FILE' found. REQURIED.") }
    if(is.null(ORDER_ID)){ stop("No 'ORDER_ID' found. REQURIED.") }
    #---------------------------------------------------------------------------
    if(is.null(SEQ_TYPE)){ SEQ_TYPE = "DNA" }else{ SEQ_TYPE = SEQ_TYPE }
    if(is.null(THREADS) ){ THREADS = 15 }else{ THREADS = as.numeric(THREADS) }
    #---------------------------------------------------------------------------
    if( SEQ_TYPE == "DNA" ){ if(is.null(NGS_LIB)){ stop("No 'NGS_LIB' found. REQUIRED When SeqType = DNA'") } }
    if( SEQ_TYPE == "RNA" ){ 
        if(is.null(WTS_ReferenceGenome) ){ WTS_ReferenceGenome = "hg19" }else{ WTS_ReferenceGenome = WTS_ReferenceGenome } 
        if( WTS_PairedSeq ){ 
            wts_seq_type = "PE"
            wts_paired   = "true"
        }else{
            wts_seq_type = "SE"
            wts_paired   = "false"
        }
        if( WTS_Stranded ){ wts_stranded = "true" }else{ wts_stranded = "false" }
    }
#--------------------------------------------------------------------------------------------------#   
#---| FUNCTION |-----------------------------------------------------------------------------------#
    Get.Fastq.Info = function(fastq_file, bed=NULL)
    {
        suppressPackageStartupMessages(library("MODIS"))
        if( is.null(bed) )
        {
            exomeSize = 32.102474
        }else{
            bed_file = sprintf("/storage/references_and_index/hg19/bed/%s/hg19_%s.target.bed", bed, bed)
            bed_info = read.delim(bed_file, header=FALSE)
            exomeSize = sum(bed_info[,3] - bed_info[,2] + 1)/1000000
            if( exomeSize >= 32.102474 ){ exomeSize = 32.102474 }
        }

        fastq_szie_mb = fileSize( fastq_file, units="MB")
        fastq_szie_gb = fileSize( fastq_file, units="GB")
        fastq_dp      = fastq_szie_mb*1.5/exomeSize
        if( fastq_dp < 100 ){ fastq_dp = round(fastq_dp, -1) }else{ fastq_dp = round(fastq_dp, -2) }

        R100 = system(sprintf("/storage/apps/bin/samtools view %s | head -100", fastq_file), intern=T) 
        fastq_len = mean(sapply( R100, function(rd) nchar(unlist(strsplit(rd, "\t"))[10]) ))

        res = data.frame(
            Fastq.File=fastq_file,
            Fastq.Size.GB=round(fastq_szie_gb, 1),
            Fastq.Size.MB=round(fastq_szie_mb, 2),
            Fastq.Depth=fastq_dp,
            Fastq.Length=round(fastq_len, -1)
        )
        return(res)
    }
    #---------------------------------------------------------------------------
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
#--------------------------------------------------------------------------------------------------#

#---| DATABASE CONNECTION |------------------------------------------------------------------------#
    if( SEQ_TYPE == "RNA" )
    { 
        source("/data/wts/params/ruo_wts_db.R") 
    }else{ 
        source("/data/wes/params/ruo_wes_db.R") 
        DB_HOST      = db_host
        DB_ID        = db_user
        DB_PORT      = db_port
        DB_PW        = db_pw
        DB_NAME_INFO = db_name_info
    }
#--------------------------------------------------------------------------------------------------#

#---| READ SAMPLE-INFO FILE (SEQ_ID INFO FILE) |---------------------------------------------------#
    META_DIR <- sprintf("%s/%s/meta", BASE_DIR, BATCH_ID)
    setwd(META_DIR)
    #--------------------------------------------------------------------------#
    fileExtension <- tail(unlist(strsplit(SAMPLE_INFO_FILE, "\\.")), 1)
    if( fileExtension %in% c("xlsx", "xls") ){ fileType <- "XLSX" }else{ fileType <- "TEXT" }
    #--------------------------------------------------------------------------#
    if( fileType == "XLSX" ){ 
        sinfo <- openxlsx::read.xlsx(sprintf("%s/%s", META_DIR, SAMPLE_INFO_FILE))
    }else{ 
        sinfo <- read.delim(sprintf("%s/%s", META_DIR, SAMPLE_INFO_FILE), header=T, check.names=FALSE) 
    }
#--------------------------------------------------------------------------------------------------#

#---| CLIENT INFO |--------------------------------------------------------------------------------#
    CLIENT_ID <- CLIENT_ID
    #--------------------------------------------------------------------------#
    dbCon       <- dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_INFO )
    client_info <- dbGetQuery(dbCon, sprintf("SELECT * FROM facility_info WHERE facility_id = '%s'", CLIENT_ID))
    dbDisconnect(dbCon)
    if(nrow(client_info) == 0){ stop(">> NO CLIENT INFO in DATABASE. REQUIRED.") }
    #--------------------------------------------------------------------------#
    if(       tail(unlist(strsplit(BASE_DIR, "\\/")), 1) == "wes" ){ PANEL <- "WES" 
    }else if( tail(unlist(strsplit(BASE_DIR, "\\/")), 1) == "wts" ){ PANEL <- "WTS" 
    }else if( tail(unlist(strsplit(BASE_DIR, "\\/")), 1) == "wgs" ){ PANEL <- "WGS" 
    }else { PANEL = "" }
#--------------------------------------------------------------------------------------------------#

#---| SEQ_FOLDER INFO |----------------------------------------------------------------------------#
    message(">> Batch-ID Info Table create and import into database ....")
    batch_info = data.frame(
        seq_folder        = BATCH_ID,
        facility_id       = CLIENT_ID,
        date              = today(),
        client_name       = client_info$facility_name,
        panel             = PANEL,
        ngs_order_id      = ORDER_ID,
        ext_hdd_id        = NA,
        ext_hdd_send_date = NA
    )
    # import into database -----------------------------------------------------
    dbCon = dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_INFO)
    Check = dbGetQuery(dbCon, sprintf("SELECT * FROM seqfolder_info WHERE ngs_order_id = '%s'", ORDER_ID))
    if( nrow(Check) == 0 )
    {
        #clear1 = dbGetQuery(dbCon, sprintf("DELETE FROM seqfolder_info WHERE ngs_order_id = '%s'", ORDER_ID))
        write1 = dbWriteTable(dbCon, name="seqfolder_info", value=batch_info, row.names=FALSE, append=TRUE)
    }
    dbDisconnect(dbCon)
#--------------------------------------------------------------------------------------------------#

#---| SEQ-ID INFO IMPORT |-------------------------------------------------------------------------#
    message(">> Seq-ID Info import into database ....")
    dbCon  = dbConnect(dbDriver("MySQL"), host=DB_HOST, user=DB_ID, port=DB_PORT, password=DB_PW, db=DB_NAME_INFO )
    clear2 = dbGetQuery(dbCon, sprintf("DELETE FROM seqid_info WHERE seq_folder = '%s'", BATCH_ID))
    write2 = dbWriteTable(dbCon, name="seqid_info", value=sinfo, row.names=FALSE, append=TRUE)
    dbDisconnect(dbCon)
#--------------------------------------------------------------------------------------------------#

#---| WES PROCESSING |-----------------------------------------------------------------------------#
    if( PANEL == "WES" )
    {
        if( DP_RD_FROM_FASTQ )
        {
            message(">> Get Fastq Info ....")
            FASTQ_DIR = sprintf("%s/fastq", BASE_DIR)
            #---------------------------------------------------------------------------
            raw_fastq_info = data.frame()
            for( i in 1:nrow(sinfo) )
            {
                fqinfo = data.frame(rbind(
                    Get.Fastq.Info(fastq_file=sprintf("%s/%s_R1.fastq.gz", FASTQ_DIR, sinfo[i, "seq_id"])) %>% mutate( SeqID = sinfo[i, "seq_id"] ),
                    Get.Fastq.Info(fastq_file=sprintf("%s/%s_R2.fastq.gz", FASTQ_DIR, sinfo[i, "seq_id"])) %>% mutate( SeqID = sinfo[i, "seq_id"] )
                ))
                raw_fastq_info = rbind(raw_fastq_info, fqinfo)
            }
            fastq_info = raw_fastq_info %>% group_by( SeqID ) %>% reframe( DP = min(Fastq.Depth), LEN = min(Fastq.Length) ) %>% as.data.frame()
            #---------------------------------------------------------------------------
            write.table( raw_fastq_info, sprintf("%s/raw_fastq_info.tsv", META_DIR), quote=F, col.names=T, row.names=F, sep="\t")
            #------------------------------------------------------------------------------------------#
            # SGE SAMPLE PROCESSING SAMPLE LIST
            message(">> SGE processing scripts create ....")
            # SEQ_ID SEQ_DEPTH FLOWCELL_ID READ_LENGTH #
            sample_list_file = data.frame()
            for( i in 1:nrow(sinfo) )
            {
                sample_list_file = rbind(
                    sample_list_file,
                    data.frame(
                        SEQ_ID      = sinfo[i, "seq_id"], 
                        SEQ_DEPTH   = fastq_info[which(fastq_info$SeqID == sinfo[i, "seq_id"]), "DP"],
                        FLOWCELL_ID = sinfo[i, "flowcell_id"],
                        READ_LENGTH = fastq_info[which(fastq_info$SeqID == sinfo[i, "seq_id"]), "LEN"]
                    )
                )
            }
        }else{
            if(is.null(DP) ){ stop("No 'DP' found. REQURIED WHEN 'DP_RD_FROM_FASTQ' = FALSE")  }
            if(is.null(RD) ){ stop("No 'RD' found. REQURIED WHEN 'DP_RD_FROM_FASTQ' = FALSE")  }
            sample_list_file = data.frame(
                SEQ_ID      = sinfo$seq_id,
                SEQ_DEPTH   = DP,
                FLOWCELL_ID = sinfo$flowcell_id,
                READ_LENGTH = RD
            )
        }
        #---------------------------------------------------------------------------
        write.table(sample_list_file, sprintf("%s/sample.list.tsv", META_DIR), quote=F, col.names=T, row.names=F, sep="\t")
        write.table(sample_list_file, sprintf("%s/sample.list", META_DIR), quote=F, col.names=F, row.names=F, sep=" ")
        #---------------------------------------------------------------------------
        default_sge_processing_script = "/storage/home/kangsm/myScripts/Default_Scripts/DEFAULT_wes_SeqID_Processing_TumorOnly_Analysis.sh"
        sge_processing_run_script     = sprintf("%s/%s_wes_SeqID_Processing_TumorOnly_Analysis.sh", META_DIR, BATCH_ID)
        system(sprintf("cp %s %s", default_sge_processing_script, sge_processing_run_script))
        #---------------------------------------------------------------------------
        sge_base_dir   = gsub("/","\\\\/", BASE_DIR)
        sge_job_name   = sprintf("%s", BATCH_ID)
        sge_total_task = nrow(sample_list_file)
        sge_job        = ifelse( nrow(sample_list_file) > 8, 8, nrow(sample_list_file ) )
        sge_threads    = THREADS
        sge_ngs_lib    = NGS_LIB
        sge_error_log  = sprintf("%s/%s/meta/analysis_log/%s.SGE.SeqID.Processing.%s.err", BASE_DIR, BATCH_ID, BATCH_ID, today() )        
        sge_out_log    = sprintf("%s/%s/meta/analysis_log/%s.SGE.SeqID.Processing.%s.out", BASE_DIR, BATCH_ID, BATCH_ID, today() )
        sge_error_log  = gsub("/","\\\\/", sge_error_log)
        sge_out_log    = gsub("/","\\\\/", sge_out_log)
        #---------------------------------------------------------------------------
        run_germline_vars = ifelse( WES_RunGermlineVariant, "ON", "OFF" )
        ascat_ts_normal   = sprintf("/storage/references_and_index/cnv/ascat/targeted_seq/%s/Normal_data/Normal_BAF.txt", NGS_LIB)
        if( file.exists(ascat_ts_normal) ){ run_ascat = "ON" }else{ run_ascat = "OFF" }
        purecn_normal     = sprintf("/storage/references_and_index/cnv/purecn/normalDB/hg19_%s/normalDB_hg19_%s_hg19.rds", NGS_LIB, NGS_LIB)
        if( file.exists(purecn_normal) ){ run_purecn = "ON" }else{ run_purecn = "OFF" }
        if( dir.exists("/storage/references_and_index/cnv/superfreq/hg19_twist.exome.2.0_pon_bam/bam") )
        { run_superfreq = "ON" }else{ run_superfreq = "OFF" }
        #---------------------------------------------------------------------------
        system(sprintf("sed -i 's/JOB_NAME_HERE/%s/' %s",   sge_job_name,   sge_processing_run_script))
        system(sprintf("sed -i 's/THREADS_N_HERE/%s/g' %s", sge_threads,    sge_processing_run_script))
        system(sprintf("sed -i 's/TOTAL_TASK_HERE/%s/' %s", sge_total_task, sge_processing_run_script))
        system(sprintf("sed -i 's/JOB_N_HERE/%s/' %s",      sge_job,        sge_processing_run_script))
        system(sprintf("sed -i 's/ERROR_LOG_HERE/%s/' %s",  sge_error_log,  sge_processing_run_script))
        system(sprintf("sed -i 's/OUT_LOG_HERE/%s/' %s",    sge_out_log,    sge_processing_run_script))
        system(sprintf("sed -i 's/BASE_DIR_HERE/%s/' %s",   sge_base_dir,   sge_processing_run_script))
        system(sprintf("sed -i 's/BATCH_ID_HERE/%s/' %s",   BATCH_ID,       sge_processing_run_script))
        system(sprintf("sed -i 's/NGS_LIB_HERE/%s/' %s",    NGS_LIB,        sge_processing_run_script))
        #-----------------------------------------------------------------------
        system(sprintf("sed -i 's/GERMLINE_VARS_RUN_OPTION/%s/g' %s", run_germline_vars, sge_processing_run_script))
        system(sprintf("sed -i 's/ASCAT_CNV_RUN_OPTION/%s/' %s",      run_ascat,         sge_processing_run_script))
        system(sprintf("sed -i 's/PURECN_CNV_RUN_OPTION/%s/' %s",     run_purecn,        sge_processing_run_script))
        system(sprintf("sed -i 's/SUPERFREQ_CNV_RUN_OPTION/%s/' %s",  run_superfreq,     sge_processing_run_script))
        #-----------------------------------------------------------------------
        vscode_dir = list.files("/storage/home/kangsm/myScripts/Analysis", pattern=BATCH_ID)
        system(sprintf("ln -Tsf %s /storage/home/kangsm/myScripts/Analysis/%s/%s_wes_SeqID_Processing_TumorOnly_Analysis.sh", sge_processing_run_script, vscode_dir, BATCH_ID))
        # MATCHED ANALYSIS
        # NORMAL_SEQ_ID TUMOR_SEQ_ID SEQ_DEPTH NORMAL_TYPE #
        sinfo_groups    = table(sinfo$sample_group)
        analysis_groups = names(sinfo_groups[sinfo_groups >= 2])
        matched_sinfo   = sinfo %>% filter( sample_group %in% analysis_groups )
        matched_normal_group = matched_sinfo %>% filter( matched_normal %in% c(2,3) ) %>% .$sample_group %>% unique()
        #-----------------------------------------------------------------------
        if( length(matched_normal_group) > 0 )
        {
            message(">> create matched-analysis scripts...")
            matched_list = data.frame()
            for( mag in matched_normal_group )
            {
                normal_ts_list  = sinfo %>% filter( matched_normal == 3, sample_group == mag ) %>% .$seq_id %>% unique()
                normal_org_list = sinfo %>% filter( matched_normal == 2, sample_group == mag ) %>% .$seq_id %>% unique()
                tumor_list      = sinfo %>% filter( matched_normal %nin% c(2,3), sample_group == mag ) %>% .$seq_id %>% unique()

                if( length(normal_ts_list) > 0 )
                {
                    for( j in 1:length(normal_ts_list) )
                    {
                        matched_list = rbind( 
                            matched_list,
                            data.frame( N_ID=normal_ts_list[j], T_ID=tumor_list, DP=NA, TYPE="TS" )
                        )
                        if( length(normal_org_list) > 0 )
                        {
                            matched_list = rbind( 
                                matched_list,
                                data.frame( N_ID=normal_ts_list[j], T_ID=normal_org_list, DP=NA, TYPE="TS" )
                            )
                        }
                    }
                }
                if( length(normal_org_list) > 0 )
                {
                    for( j in 1:length(normal_ts_list) )
                    {
                        matched_list = rbind( 
                            matched_list,
                            data.frame( N_ID=normal_org_list[j], T_ID=tumor_list, DP=NA, TYPE="ORG" )
                        )
                    }
                }
            }
            matched_list$DP = fastq_info[match(matched_list$T_ID,fastq_info$SeqID), "DP"]
            matched_list = matched_list %>% filter( !is.na(N_ID), !is.na(T_ID) ) %>% unique()
            #-----------------------------------------------------------------------
            write.table(matched_list, sprintf("%s/matched_analysis.list.tsv", META_DIR), quote=F, col.names=T, row.names=F, sep="\t")
            write.table(matched_list, sprintf("%s/matched_analysis.list", META_DIR),     quote=F, col.names=F, row.names=F, sep=" " )
            #-----------------------------------------------------------------------
            default_sge_analysis_script = "/storage/home/kangsm/myScripts/Default_Scripts/DEFAULT_wes_SampleGroup_MatchedNormal_Analysis.sh"
            sge_analysis_run_script     = sprintf("%s/%s_wes_SampleGroup_MatchedNormal_Analysis.sh", META_DIR, BATCH_ID)
            system(sprintf("cp %s %s", default_sge_analysis_script, sge_analysis_run_script))
            #-----------------------------------------------------------------------
            nt.sge_base_dir   = gsub("/","\\\\/", BASE_DIR)
            nt.sge_job_name   = sprintf("NT_%s", BATCH_ID)
            nt.sge_total_task = nrow(matched_list)
            nt.sge_job        = ifelse( nrow(matched_list) > 8, 8, nrow(matched_list ) )
            nt.sge_threads    = THREADS
            nt.sge_ngs_lib    = NGS_LIB
            nt.sge_error_log  = sprintf("%s/%s/meta/analysis_log/%s.SGE.Matched.Analysis.%s.err", BASE_DIR, BATCH_ID, BATCH_ID, today() )        
            nt.sge_out_log    = sprintf("%s/%s/meta/analysis_log/%s.SGE.Matched.Analysis.%s.out", BASE_DIR, BATCH_ID, BATCH_ID, today() )
            nt.sge_error_log  = gsub("/","\\\\/", nt.sge_error_log)
            nt.sge_out_log    = gsub("/","\\\\/", nt.sge_out_log)
            #-----------------------------------------------------------------------
            system(sprintf("sed -i 's/JOB_NAME_HERE/%s/' %s",   nt.sge_job_name,   sge_analysis_run_script))
            system(sprintf("sed -i 's/THREADS_N_HERE/%s/g' %s", nt.sge_threads,    sge_analysis_run_script))
            system(sprintf("sed -i 's/TOTAL_TASK_HERE/%s/' %s", nt.sge_total_task, sge_analysis_run_script))
            system(sprintf("sed -i 's/JOB_N_HERE/%s/' %s",      nt.sge_job,        sge_analysis_run_script))
            system(sprintf("sed -i 's/ERROR_LOG_HERE/%s/' %s",  nt.sge_error_log,  sge_analysis_run_script))
            system(sprintf("sed -i 's/OUT_LOG_HERE/%s/' %s",    nt.sge_out_log,    sge_analysis_run_script))
            system(sprintf("sed -i 's/BASE_DIR_HERE/%s/' %s",   nt.sge_base_dir,   sge_analysis_run_script))
            system(sprintf("sed -i 's/BATCH_ID_HERE/%s/' %s",   BATCH_ID,          sge_analysis_run_script))
            system(sprintf("sed -i 's/NGS_LIB_HERE/%s/' %s",    NGS_LIB,           sge_analysis_run_script))
            #-----------------------------------------------------------------------
            vscode_dir = list.files("/storage/home/kangsm/myScripts/Analysis", pattern=BATCH_ID)
            system(sprintf("ln -Tsf %s /storage/home/kangsm/myScripts/Analysis/%s/%s_wes_SampleGroup_MatchedNormal_Analysis.sh", sge_analysis_run_script, vscode_dir, BATCH_ID))
            #-----------------------------------------------------------------------
            matched_analysis_script_created = TRUE
        }else{
            matched_analysis_script_created = FALSE
            message(">> NO matched-analysis samples.")
        }
        #-----------------------------------------------------------------------
        wes_batch_info = data.frame( INFO = c(
            sprintf("[ Batch ID    ] : %s", BATCH_ID), 
            sprintf("[ Samples     ] : %s", nrow(sample_list_file)),
            sprintf("[ NGS Library ] : %s", NGS_LIB)
        ))
        write.table(wes_batch_info, sprintf("%s/wes_info.tsv", META_DIR), quote=F, col.names=T, row.names=F, sep="\t")
        #---------------------------------------------------------------------------
        message(">> READY TO RUN WES PROCESSING and ANALYSIS.")
    }
#--------------------------------------------------------------------------------------------------#
    if( PANEL == "WTS" )
    {
        # sample.list column
        # SEQ_ID SEQ_TYPE PAIRED THREADS STRAND GENOME_ASSEMBLY #
        sample_list = data.frame(
            SEQ_ID          = sinfo$seq_id,
            SEQ_TYPE        = wts_seq_type,
            PAIRED          = wts_paired,
            THREADS         = THREADS,
            STRAND          = wts_stranded,
            GENOME_ASSEMBLY = WTS_ReferenceGenome
        )
        # write sample.list
        write.table(sample_list, sprintf("%s/sample.list.tsv", META_DIR), quote=F, col.names=T, row.names=F, sep="\t")
        write.table(sample_list, sprintf("%s/sample.list", META_DIR), quote=F, col.names=F, row.names=F, sep=" ")
        # script copy
        default_wts_processing_script = "/storage/home/kangsm/myScripts/Default_Scripts/DEFAULT_wts_SeqID_Processing.sh"
        sge_processing_run_script = sprintf("%s/WTS_run.SGE.SeqID.Processing.sh", META_DIR )
        system(sprintf("cp %s %s", default_wts_processing_script, sge_processing_run_script))
        #-----------------------------------------------------------------------
        sge_base_dir   = gsub("/","\\\\/", BASE_DIR)
        sge_job_name   = sprintf("%s", BATCH_ID)
        sge_total_task = nrow(sample_list_file)
        sge_job        = ifelse( nrow(sample_list_file) > 8, 8, nrow(sample_list_file ) )
        sge_threads    = THREADS
        sge_ngs_lib    = NGS_LIB
        sge_error_log  = sprintf("%s/%s/meta/analysis_log/%s.SGE.SeqID.Processing.%s.err", BASE_DIR, BATCH_ID, BATCH_ID, today() )        
        sge_out_log    = sprintf("%s/%s/meta/analysis_log/%s.SGE.SeqID.Processing.%s.out", BASE_DIR, BATCH_ID, BATCH_ID, today() )
        sge_error_log  = gsub("/","\\\\/", sge_error_log)
        sge_out_log    = gsub("/","\\\\/", sge_out_log)
        #---------------------------------------------------------------------------
        run_germline_vars = ifelse( WES_RunGermlineVariant, "ON", "OFF" )
        ascat_ts_normal   = sprintf("/storage/references_and_index/cnv/ascat/targeted_seq/%s/Normal_data/Normal_BAF.txt", NGS_LIB)
        if( file.exists(ascat_ts_normal) ){ run_ascat = "ON" }else{ run_ascat = "OFF" }
        purecn_normal     = sprintf("/storage/references_and_index/cnv/purecn/normalDB/hg19_%s/normalDB_hg19_%s_hg19.rds", NGS_LIB, NGS_LIB)
        if( file.exists(purecn_normal) ){ run_purecn = "ON" }else{ run_purecn = "OFF" }
        if( dir.exists("/storage/references_and_index/cnv/superfreq/hg19_twist.exome.2.0_pon_bam/bam") )
        { run_superfreq = "ON" }else{ run_superfreq = "OFF" }
        #---------------------------------------------------------------------------
        system(sprintf("sed -i 's/JOB_NAME_HERE/%s/' %s",   sge_job_name,   sge_processing_run_script))
        system(sprintf("sed -i 's/THREADS_N/%s/g' %s", sge_threads,    sge_processing_run_script))
        system(sprintf("sed -i 's/1-TOTAL_TAKS_N/%s/' %s", sge_total_task, sge_processing_run_script))
        system(sprintf("sed -i 's/TASK_N/%s/' %s",      sge_job,        sge_processing_run_script))
        system(sprintf("sed -i 's/ERROR_LOG/%s/' %s",  sge_error_log,  sge_processing_run_script))
        system(sprintf("sed -i 's/OUTPUT_LOG/%s/' %s",    sge_out_log,    sge_processing_run_script))
        #system(sprintf("sed -i 's/BASE_DIR_HERE/%s/' %s",   sge_base_dir,   sge_processing_run_script))
        system(sprintf("sed -i 's/SEQ_FOLDER_ID/%s/' %s",   BATCH_ID,       sge_processing_run_script))
        #system(sprintf("sed -i 's/NGS_LIB_HERE/%s/' %s",    NGS_LIB,        sge_processing_run_script))
        #-----------------------------------------------------------------------
        #system(sprintf("sed -i 's/GERMLINE_VARS_RUN_OPTION/%s/g' %s", run_germline_vars, sge_processing_run_script))
        #system(sprintf("sed -i 's/ASCAT_CNV_RUN_OPTION/%s/' %s",      run_ascat,         sge_processing_run_script))
        #system(sprintf("sed -i 's/PURECN_CNV_RUN_OPTION/%s/' %s",     run_purecn,        sge_processing_run_script))
        #system(sprintf("sed -i 's/SUPERFREQ_CNV_RUN_OPTION/%s/' %s",  run_superfreq,     sge_processing_run_script))

    
    }

###