# ============================================================
# SCRIPT  : ngs_batch_prep.R
# VERSION : v2.0
# DATE    : 2024-03-11
# AUTHOR  : kangsm@gencurix.com
#
# PURPOSE : Per-sample SGE shell script generation for WES/WTS
#           - Reads sample info file (xlsx or tsv)
#           - Registers batch & sample metadata to MySQL DB
#           - Generates one qsub-ready .sh per sample in meta/scripts/
# ============================================================

# ── 0. Packages ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(lubridate)
    library(RMySQL)
    library(openxlsx)
    library(Hmisc)       # %nin%
})

options(stringsAsFactors = FALSE)


# ── 1. CLI Arguments ─────────────────────────────────────────────────────────
option_list <- list(
    make_option("--BASE_DIR",            type = "character", help = "NGS base directory (e.g. /data/wes)"),
    make_option("--BATCH_ID",            type = "character", help = "Batch ID (= sequencing folder name)"),
    make_option("--CLIENT_ID",           type = "character", help = "Client (facility) ID"),
    make_option("--SAMPLE_INFO_FILE",    type = "character", help = "Sample info file name (xlsx or tsv)"),
    make_option("--ORDER_ID",            type = "character", help = "NGS order ID"),
    make_option("--SEQ_TYPE",            type = "character", default = "DNA",  help = "DNA | RNA  [default: DNA]"),
    make_option("--THREADS",             type = "integer",   default = 15,     help = "Threads per job [default: 15]"),
    make_option("--DP_RD_FROM_FASTQ",    type = "logical",   default = FALSE,  help = "Infer depth/read-length from FASTQ"),
    make_option("--DP",                  type = "character", default = NULL,   help = "Depth (required when DP_RD_FROM_FASTQ = FALSE)"),
    make_option("--RD",                  type = "character", default = NULL,   help = "Read length (required when DP_RD_FROM_FASTQ = FALSE)"),
    make_option("--WES_NGS_LIB",         type = "character", default = NULL,   help = "NGS library key (WES only)"),
    make_option("--WES_RunGermlineVariant", type = "logical", default = FALSE, help = "Run germline variant calling (WES only)"),
    make_option("--WTS_Stranded",        type = "logical",   default = FALSE,  help = "Stranded library (WTS only)"),
    make_option("--WTS_PairedSeq",       type = "logical",   default = TRUE,   help = "Paired-end sequencing (WTS only)"),
    make_option("--WTS_ReferenceGenome", type = "character", default = "hg19", help = "Reference genome (WTS only) [default: hg19]")
)

args <- parse_args(OptionParser(option_list = option_list))


# ── 2. Parameter Validation & Defaults ───────────────────────────────────────
required_args <- c("BASE_DIR", "BATCH_ID", "CLIENT_ID", "SAMPLE_INFO_FILE", "ORDER_ID")
for (arg in required_args) {
    if (is.null(args[[arg]])) stop(sprintf("Missing required argument: --%s", arg))
}

BASE_DIR         <- args$BASE_DIR
BATCH_ID         <- args$BATCH_ID
CLIENT_ID        <- args$CLIENT_ID
SAMPLE_INFO_FILE <- args$SAMPLE_INFO_FILE
ORDER_ID         <- args$ORDER_ID
SEQ_TYPE         <- args$SEQ_TYPE
THREADS          <- args$THREADS
DP_RD_FROM_FASTQ <- args$DP_RD_FROM_FASTQ
DP               <- args$DP
RD               <- args$RD
NGS_LIB          <- args$WES_NGS_LIB
RUN_GERMLINE     <- args$WES_RunGermlineVariant
WTS_STRANDED     <- args$WTS_Stranded
WTS_PAIRED       <- args$WTS_PairedSeq
WTS_GENOME       <- args$WTS_ReferenceGenome

# Panel detection from BASE_DIR tail
PANEL <- switch(
    basename(BASE_DIR),
    wes = "WES",
    wts = "WTS",
    wgs = "WGS",
    ""
)

# WES-specific checks
if (PANEL == "WES") {
    if (is.null(NGS_LIB)) stop("--WES_NGS_LIB is required when SEQ_TYPE = DNA")
    if (!DP_RD_FROM_FASTQ) {
        if (is.null(DP)) stop("--DP is required when --DP_RD_FROM_FASTQ = FALSE")
        if (is.null(RD)) stop("--RD is required when --DP_RD_FROM_FASTQ = FALSE")
    }
}

# WTS-specific derived values
if (PANEL == "WTS") {
    wts_seq_type <- ifelse(WTS_PAIRED, "PE", "SE")
    wts_paired   <- ifelse(WTS_PAIRED, "true", "false")
    wts_stranded <- ifelse(WTS_STRANDED, "true", "false")
}


# ── 3. Helper Functions ───────────────────────────────────────────────────────

#' Connect to MySQL, run a query, then disconnect
#'
#' Wraps every statement in a tryCatch so FK / syntax errors surface clearly.
#' @param query  SQL string
#' @param db_cfg Named list: host, user, port, password, dbname
#' @return data.frame of results (empty data.frame for non-SELECT statements)
db_query <- function(query, db_cfg) {
    con <- dbConnect(dbDriver("MySQL"),
                     host     = db_cfg$host,
                     user     = db_cfg$user,
                     port     = db_cfg$port,
                     password = db_cfg$password,
                     db       = db_cfg$dbname)
    on.exit(dbDisconnect(con))
    tryCatch(
        dbGetQuery(con, query),
        error = function(e) stop(sprintf("[DB ERROR] %s\nQuery: %s", conditionMessage(e), query))
    )
}

# ── Child tables that reference seqid_info (FK ON DELETE NO ACTION) ──────────
# Add any additional child tables here to keep the delete order maintainable.
SEQID_INFO_CHILD_TABLES <- c(
    "hla_type"
    # e.g. "tmb_result", "msi_result"   ← append as needed
)

db_write <- function(table_name, df, db_cfg) {
    con <- dbConnect(dbDriver("MySQL"),
                     host     = db_cfg$host,
                     user     = db_cfg$user,
                     port     = db_cfg$port,
                     password = db_cfg$password,
                     db       = db_cfg$dbname)
    on.exit(dbDisconnect(con))
    dbWriteTable(con, name = table_name, value = df, row.names = FALSE, append = TRUE)
}

#' Read FASTQ size/depth/read-length for one sample
get_fastq_info <- function(seq_id, fastq_dir, bed = NULL) {
    suppressPackageStartupMessages(library(MODIS))

    exome_size_mb <- if (is.null(bed)) {
        32.102474
    } else {
        bed_path <- sprintf("/storage/references_and_index/hg19/bed/%s/hg19_%s.target.bed", bed, bed)
        bed_tbl  <- read.delim(bed_path, header = FALSE)
        min(sum(bed_tbl[, 3] - bed_tbl[, 2] + 1) / 1e6, 32.102474)
    }

    read_fastq <- function(r_tag) {
        fq_path    <- sprintf("%s/%s_%s.fastq.gz", fastq_dir, seq_id, r_tag)
        size_mb    <- fileSize(fq_path, units = "MB")
        size_gb    <- fileSize(fq_path, units = "GB")
        depth      <- size_mb * 1.5 / exome_size_mb
        depth      <- if (depth < 100) round(depth, -1) else round(depth, -2)
        raw_reads  <- system(sprintf("samtools view %s | head -100", fq_path), intern = TRUE)
        read_len   <- mean(sapply(raw_reads, function(r) nchar(strsplit(r, "\t")[[1]][10])))
        data.frame(
            SeqID        = seq_id,
            Tag          = r_tag,
            Fastq.File   = fq_path,
            Fastq.Size.GB = round(size_gb, 1),
            Fastq.Size.MB = round(size_mb, 2),
            Fastq.Depth  = depth,
            Fastq.Length = round(read_len, -1)
        )
    }
    rbind(read_fastq("R1"), read_fastq("R2"))
}

#' Escape forward slashes for sed replacement strings
escape_for_sed <- function(x) gsub("/", "\\\\/", x)

#' Write text lines to a file
write_lines_to_file <- function(lines, path) {
    writeLines(lines, con = path)
    invisible(path)
}


# ── 4. Read Sample Info File ──────────────────────────────────────────────────
META_DIR    <- file.path(BASE_DIR, BATCH_ID, "meta")
SCRIPT_DIR  <- file.path(META_DIR, "scripts")   # per-sample sh output directory
LOG_DIR     <- file.path(META_DIR, "analysis_log")

dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,    recursive = TRUE, showWarnings = FALSE)

ext      <- tolower(tools::file_ext(SAMPLE_INFO_FILE))
sinfo_path <- file.path(META_DIR, SAMPLE_INFO_FILE)

sinfo <- if (ext %in% c("xlsx", "xls")) {
    openxlsx::read.xlsx(sinfo_path)
} else {
    read.delim(sinfo_path, header = TRUE, check.names = FALSE)
}

message(sprintf(">> Loaded %d samples from '%s'", nrow(sinfo), SAMPLE_INFO_FILE))


# ── 5. Database Config ────────────────────────────────────────────────────────
db_source <- if (PANEL == "WTS") "/data/wts/params/ruo_wts_db.R" else "/data/wes/params/ruo_wes_db.R"
source(db_source)

db_cfg <- list(
    host     = db_host,
    user     = db_user,
    port     = db_port,
    password = db_pw,
    dbname   = db_name_info
)


# ── 6. Client Info ────────────────────────────────────────────────────────────
client_info <- db_query(
    sprintf("SELECT * FROM facility_info WHERE facility_id = '%s'", CLIENT_ID),
    db_cfg
)
if (nrow(client_info) == 0) stop(sprintf("No client found for facility_id = %s", CLIENT_ID))


# ── 7. Register Batch Info ────────────────────────────────────────────────────
message(">> Registering batch info to DB ...")

batch_info <- data.frame(
    seq_folder        = BATCH_ID,
    facility_id       = CLIENT_ID,
    date              = today(),
    client_name       = client_info$facility_name,
    panel             = PANEL,
    ngs_order_id      = ORDER_ID,
    ext_hdd_id        = NA,
    ext_hdd_send_date = NA
)

existing_batch <- db_query(
    sprintf("SELECT * FROM seqfolder_info WHERE ngs_order_id = '%s'", ORDER_ID),
    db_cfg
)
if (nrow(existing_batch) == 0) db_write("seqfolder_info", batch_info, db_cfg)


# ── 8. Register Sample Info ───────────────────────────────────────────────────
message(">> Registering sample info to DB ...")

#' Delete rows from seqid_info safely, honouring FK constraints.
#'
#' Child tables that hold a FK → seqid_info.seq_id (ON DELETE NO ACTION) must
#' be cleared first; otherwise MySQL raises ER_ROW_IS_REFERENCED_2.
#'
#' Strategy:
#'   1. Collect seq_ids about to be deleted.
#'   2. Delete matching rows from every known child table (gcx_wes schema).
#'   3. Delete from seqid_info (parent).
#'
#' @param batch_id  BATCH_ID / seq_folder value to scope the delete
#' @param db_cfg    DB connection config list
delete_seqid_info_safe <- function(batch_id, db_cfg) {

    # 1. Identify seq_ids in scope
    seq_ids_df <- db_query(
        sprintf("SELECT seq_id FROM gcx_ngs_service.seqid_info WHERE seq_folder = '%s'", batch_id),
        db_cfg
    )

    if (nrow(seq_ids_df) == 0) {
        message("   No existing rows for this batch — skipping delete.")
        return(invisible(NULL))
    }

    seq_id_list <- paste0("'", seq_ids_df$seq_id, "'", collapse = ", ")

    # 2. Delete from child tables first (FK children → parent order)
    for (child_table in SEQID_INFO_CHILD_TABLES) {
        child_query <- sprintf(
            "DELETE FROM gcx_wes.%s WHERE seq_id IN (%s)",
            child_table, seq_id_list
        )
        tryCatch({
            db_query(child_query, db_cfg)
            message(sprintf("   Cleared child table: %s", child_table))
        }, error = function(e) {
            # Table may not exist in every environment — warn but don't abort
            warning(sprintf("   Could not clear '%s': %s", child_table, conditionMessage(e)))
        })
    }

    # 3. Delete from parent
    db_query(
        sprintf("DELETE FROM gcx_ngs_service.seqid_info WHERE seq_folder = '%s'", batch_id),
        db_cfg
    )
    message(sprintf("   Deleted %d row(s) from seqid_info.", nrow(seq_ids_df)))
}

delete_seqid_info_safe(BATCH_ID, db_cfg)
db_write("seqid_info", sinfo, db_cfg)


# ── 9. CNV Tool Availability (WES) ───────────────────────────────────────────
check_cnv_tools <- function(ngs_lib) {
    ascat_path   <- sprintf("/storage/references_and_index/cnv/ascat/targeted_seq/%s/Normal_data/Normal_BAF.txt", ngs_lib)
    purecn_path  <- sprintf("/storage/references_and_index/cnv/purecn/normalDB/hg19_%s/normalDB_hg19_%s_hg19.rds", ngs_lib, ngs_lib)
    superfreq_dir <- "/storage/references_and_index/cnv/superfreq/hg19_twist.exome.2.0_pon_bam/bam"
    list(
        ascat     = ifelse(file.exists(ascat_path),       "ON", "OFF"),
        purecn    = ifelse(file.exists(purecn_path),      "ON", "OFF"),
        superfreq = ifelse(dir.exists(superfreq_dir),     "ON", "OFF")
    )
}


# ── 10. WES: Per-Sample Shell Script Generation ───────────────────────────────
generate_wes_sample_scripts <- function(sinfo, fastq_info = NULL) {

    default_template <- "/storage/home/kangsm/myScripts/Default_Scripts/DEFAULT_wes_SeqID_Processing_TumorOnly_Analysis.sh"
    cnv_flags        <- check_cnv_tools(NGS_LIB)
    germline_flag    <- ifelse(RUN_GERMLINE, "ON", "OFF")

    for (i in seq_len(nrow(sinfo))) {
        seq_id      <- sinfo[i, "seq_id"]
        flowcell_id <- sinfo[i, "flowcell_id"]

        # Depth / read-length
        if (!is.null(fastq_info)) {
            sample_dp  <- fastq_info[fastq_info$SeqID == seq_id, "DP"]
            sample_len <- fastq_info[fastq_info$SeqID == seq_id, "LEN"]
        } else {
            sample_dp  <- DP
            sample_len <- RD
        }

        # Paths
        out_sh    <- file.path(SCRIPT_DIR, sprintf("%s_wes_processing.sh", seq_id))
        err_log   <- file.path(LOG_DIR, sprintf("%s.SGE.Processing.%s.err", seq_id, today()))
        out_log   <- file.path(LOG_DIR, sprintf("%s.SGE.Processing.%s.out", seq_id, today()))

        # Copy template
        file.copy(default_template, out_sh, overwrite = TRUE)

        # Build sed substitution map
        subs <- list(
            JOB_NAME_HERE          = seq_id,
            THREADS_N_HERE         = as.character(THREADS),
            TOTAL_TASK_HERE        = "1",
            JOB_N_HERE             = "1",
            ERROR_LOG_HERE         = escape_for_sed(err_log),
            OUT_LOG_HERE           = escape_for_sed(out_log),
            BASE_DIR_HERE          = escape_for_sed(BASE_DIR),
            BATCH_ID_HERE          = BATCH_ID,
            NGS_LIB_HERE           = NGS_LIB,
            SEQ_ID_HERE            = seq_id,
            FLOWCELL_ID_HERE       = flowcell_id,
            SEQ_DEPTH_HERE         = as.character(sample_dp),
            READ_LENGTH_HERE       = as.character(sample_len),
            GERMLINE_VARS_RUN_OPTION = germline_flag,
            ASCAT_CNV_RUN_OPTION   = cnv_flags$ascat,
            PURECN_CNV_RUN_OPTION  = cnv_flags$purecn,
            SUPERFREQ_CNV_RUN_OPTION = cnv_flags$superfreq
        )

        apply_sed_substitutions(out_sh, subs)
        message(sprintf("   [WES] Script created: %s", basename(out_sh)))
    }
}


# ── 11. WTS: Per-Sample Shell Script Generation ───────────────────────────────
generate_wts_sample_scripts <- function(sinfo) {

    default_template <- "/storage/home/kangsm/myScripts/Default_Scripts/DEFAULT_wts_SeqID_Processing.sh"

    for (i in seq_len(nrow(sinfo))) {
        seq_id <- sinfo[i, "seq_id"]

        out_sh  <- file.path(SCRIPT_DIR, sprintf("%s_wts_processing.sh", seq_id))
        err_log <- file.path(LOG_DIR, sprintf("%s.SGE.Processing.%s.err", seq_id, today()))
        out_log <- file.path(LOG_DIR, sprintf("%s.SGE.Processing.%s.out", seq_id, today()))

        file.copy(default_template, out_sh, overwrite = TRUE)

        subs <- list(
            JOB_NAME_HERE    = seq_id,
            THREADS_N        = as.character(THREADS),
            `1-TOTAL_TAKS_N` = "1",
            TASK_N           = "1",
            ERROR_LOG        = escape_for_sed(err_log),
            OUTPUT_LOG       = escape_for_sed(out_log),
            SEQ_FOLDER_ID    = BATCH_ID,
            SEQ_ID_HERE      = seq_id,
            SEQ_TYPE_HERE    = wts_seq_type,
            PAIRED_HERE      = wts_paired,
            STRAND_HERE      = wts_stranded,
            GENOME_HERE      = WTS_GENOME
        )

        apply_sed_substitutions(out_sh, subs)
        message(sprintf("   [WTS] Script created: %s", basename(out_sh)))
    }
}


#' Apply a named list of sed substitutions to a file (in-place)
apply_sed_substitutions <- function(file_path, subs) {
    for (pattern in names(subs)) {
        replacement <- subs[[pattern]]
        cmd <- sprintf("sed -i 's/%s/%s/g' %s", pattern, replacement, file_path)
        system(cmd)
    }
}


# ── 12. Dispatch by Panel ─────────────────────────────────────────────────────
if (PANEL == "WES") {

    message(">> [WES] Preparing sample list ...")

    # --- Fastq info -----------------------------------------------------------
    if (DP_RD_FROM_FASTQ) {
        message(">> [WES] Extracting depth & read-length from FASTQ ...")
        FASTQ_DIR     <- file.path(BASE_DIR, "fastq")
        raw_fq_info   <- do.call(rbind, lapply(sinfo$seq_id, get_fastq_info, fastq_dir = FASTQ_DIR))
        fastq_info    <- raw_fq_info %>%
                         group_by(SeqID) %>%
                         reframe(DP = min(Fastq.Depth), LEN = min(Fastq.Length)) %>%
                         as.data.frame()
        write.table(raw_fq_info, file.path(META_DIR, "raw_fastq_info.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    } else {
        fastq_info <- NULL
    }

    # --- sample.list (for reference) ------------------------------------------
    sample_list <- data.frame(
        SEQ_ID      = sinfo$seq_id,
        SEQ_DEPTH   = if (!is.null(fastq_info)) fastq_info$DP[match(sinfo$seq_id, fastq_info$SeqID)] else DP,
        FLOWCELL_ID = sinfo$flowcell_id,
        READ_LENGTH = if (!is.null(fastq_info)) fastq_info$LEN[match(sinfo$seq_id, fastq_info$SeqID)] else RD
    )
    write.table(sample_list, file.path(META_DIR, "sample.list.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(sample_list, file.path(META_DIR, "sample.list"),     quote = FALSE, sep = " ",  row.names = FALSE, col.names = FALSE)

    # --- Per-sample scripts ---------------------------------------------------
    message(">> [WES] Generating per-sample SGE scripts ...")
    generate_wes_sample_scripts(sinfo, fastq_info)

    # --- Matched-normal analysis list -----------------------------------------
    sinfo_groups        <- table(sinfo$sample_group)
    matched_groups      <- names(sinfo_groups[sinfo_groups >= 2])
    matched_sinfo       <- sinfo %>% filter(sample_group %in% matched_groups)
    normal_groups       <- matched_sinfo %>% filter(matched_normal %in% c(2, 3)) %>% pull(sample_group) %>% unique()

    if (length(normal_groups) > 0) {
        message(">> [WES] Building matched-normal analysis list ...")

        matched_list <- do.call(rbind, lapply(normal_groups, function(grp) {
            ts_normals  <- sinfo %>% filter(matched_normal == 3, sample_group == grp) %>% pull(seq_id)
            org_normals <- sinfo %>% filter(matched_normal == 2, sample_group == grp) %>% pull(seq_id)
            tumors      <- sinfo %>% filter(matched_normal %nin% c(2, 3), sample_group == grp) %>% pull(seq_id)

            rows <- data.frame()
            if (length(ts_normals) > 0) {
                rows <- rbind(rows, expand.grid(N_ID = ts_normals, T_ID = tumors,      TYPE = "TS",  stringsAsFactors = FALSE))
                if (length(org_normals) > 0)
                    rows <- rbind(rows, expand.grid(N_ID = ts_normals, T_ID = org_normals, TYPE = "TS",  stringsAsFactors = FALSE))
            }
            if (length(org_normals) > 0)
                rows <- rbind(rows, expand.grid(N_ID = org_normals, T_ID = tumors,     TYPE = "ORG", stringsAsFactors = FALSE))
            rows
        }))

        if (!is.null(fastq_info)) {
            matched_list$DP <- fastq_info$DP[match(matched_list$T_ID, fastq_info$SeqID)]
        } else {
            matched_list$DP <- DP
        }

        matched_list <- matched_list %>% filter(!is.na(N_ID), !is.na(T_ID)) %>% distinct()
        write.table(matched_list, file.path(META_DIR, "matched_analysis.list.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
        write.table(matched_list, file.path(META_DIR, "matched_analysis.list"),     quote = FALSE, sep = " ",  row.names = FALSE, col.names = FALSE)
        message(sprintf(">> [WES] Matched-normal pairs written: %d", nrow(matched_list)))
    } else {
        message(">> [WES] No matched-normal samples found.")
    }

    # --- Summary --------------------------------------------------------------
    summary_info <- data.frame(INFO = c(
        sprintf("[ Batch ID    ] : %s", BATCH_ID),
        sprintf("[ Samples     ] : %d", nrow(sample_list)),
        sprintf("[ NGS Library ] : %s", NGS_LIB),
        sprintf("[ Script Dir  ] : %s", SCRIPT_DIR)
    ))
    write.table(summary_info, file.path(META_DIR, "wes_batch_info.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    message(">> [WES] All scripts ready. Submit with:")
    message(sprintf("   for f in %s/*_wes_processing.sh; do qsub \"$f\"; done", SCRIPT_DIR))
}


if (PANEL == "WTS") {

    message(">> [WTS] Generating per-sample SGE scripts ...")

    # --- sample.list (for reference) ------------------------------------------
    sample_list <- data.frame(
        SEQ_ID          = sinfo$seq_id,
        SEQ_TYPE        = wts_seq_type,
        PAIRED          = wts_paired,
        THREADS         = THREADS,
        STRAND          = wts_stranded,
        GENOME_ASSEMBLY = WTS_GENOME
    )
    write.table(sample_list, file.path(META_DIR, "sample.list.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(sample_list, file.path(META_DIR, "sample.list"),     quote = FALSE, sep = " ",  row.names = FALSE, col.names = FALSE)

    # --- Per-sample scripts ---------------------------------------------------
    generate_wts_sample_scripts(sinfo)

    message(">> [WTS] All scripts ready. Submit with:")
    message(sprintf("   for f in %s/*_wts_processing.sh; do qsub \"$f\"; done", SCRIPT_DIR))
}

message(">> Done.")