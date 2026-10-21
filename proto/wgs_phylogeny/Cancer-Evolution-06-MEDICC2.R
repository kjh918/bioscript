
###
    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
        library("GenomicRanges")
        library("data.table")
    })
    indexing = function(X) { sapply(unique(X), function(y) list(y)) }
###
###
    merge_same_cn_segments <- function(df) {
        dt <- as.data.table(df)
        setorder(dt, sample_id, chrom, start, end)

        dt[, grp := cumsum(
            is.na(shift(end)) |
            start != shift(end) + 1 |
            cn_a != shift(cn_a) |
            cn_b != shift(cn_b)
        ), by = .(sample_id, chrom)]

        out <- dt[, .(
            start = min(start),
            end = max(end),
            cn_a = first(cn_a),
            cn_b = first(cn_b)
        ), by = .(sample_id, chrom, grp)]

        out[, grp := NULL]
        setorder(out, sample_id, chrom, start)

        return(out)
    }

    make_common_ranges <- function(gr_list) {
        all_seq <- unique(unlist(lapply(gr_list, function(gr) as.character(seqnames(gr)))))

        common <- lapply(all_seq, function(chr) {
            bps <- sort(unique(unlist(lapply(gr_list, function(gr) {
                x <- gr[as.character(seqnames(gr)) == chr]
                c(start(x), end(x) + 1L)
            }))))

            bps <- bps[!is.na(bps)]

            if (length(bps) < 2) return(NULL)

            GRanges(
                seqnames = chr,
                ranges = IRanges(
                    start = bps[-length(bps)],
                    end   = bps[-1] - 1L
                )
            )
        })

        do.call(c, common)
    }

    remap_to_common <- function(common, gr, metadata_cols = colnames(mcols(gr))) {
        hits <- findOverlaps(common, gr, type = "within")

        if (length(hits) == 0) {
            stop("No overlaps found.")
        }

        hit_df <- data.frame(
            q = queryHits(hits),
            s = subjectHits(hits),
            width = width(pintersect(common[queryHits(hits)], gr[subjectHits(hits)]))
        )

        # common segment가 여러 원본 segment와 겹치면 가장 overlap 큰 것 선택
        hit_df <- hit_df[order(hit_df$q, -hit_df$width), ]
        hit_df <- hit_df[!duplicated(hit_df$q), ]

        out <- common

        for (col in metadata_cols) {
            vals <- rep(NA, length(common))
            vals[hit_df$q] <- mcols(gr)[[col]][hit_df$s]
            mcols(out)[[col]] <- vals
        }

        out
    }
###
###
    BaseDir            <- "/data/wgs"
    BatchID            <- "WGS-26-01"
    work_dir           <- "/data/wgs/WGS-26-01/analysis/medicc2"
    analysis_list_file <- "/data/wgs/WGS-26-01/meta/analysis.list"
    sinfo <- read.table(sprintf("%s",analysis_list_file), header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    tumor_samples <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID
###
###
    cnv_data_list = lapply(  indexing(tumor_samples), function(tid) {
        purple_dir = sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)       
        cnv_data = read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",purple_dir,tid))        

        cnv_res <- cnv_data %>% mutate( 
            cn_a=floor(majorAlleleCopyNumber+0.5), 
            cn_b=floor(minorAlleleCopyNumber+0.5) 
        ) %>% mutate(
            sample_id=tid
        ) %>% dplyr::select( 
            c("sample_id","chromosome","start","end","cn_a","cn_b") 
        ) %>% dplyr::rename( 
            chrom=chromosome 
        ) %>% mutate(
            temp_max = pmax(cn_a, cn_b),
            temp_min = pmin(cn_a, cn_b),
            cn_a = temp_max,
            cn_b = temp_min
        ) %>% dplyr::select(-temp_max, -temp_min)

        return(cnv_res)
    })
###
###
    cnv_data_list_rev = lapply( cnv_data_list, function(y) merge_same_cn_segments(y) )
    cnv_data_gr       = lapply( cnv_data_list_rev, function(y) {
        gr <- GRanges( y$chrom, IRanges(start=y$start, end=y$end), cn_a=y$cn_a, cn_b=y$cn_b)   
        seqlevels(gr) <- c(paste0(1:22), "X","Y")
        gr <- sort(gr)
        return(gr)
    })
    common_gr <- suppressWarnings(make_common_ranges(cnv_data_gr))
###
###
    medicc2_input <- data.frame()

    for( j in 1:length(cnv_data_gr) ){
        
        remap_data <- remap_to_common(common_gr, cnv_data_gr[[j]] ) %>% 
            as.data.frame() %>% 
            mutate( sample_id=names(cnv_data_gr)[j] ) %>% 
            dplyr::rename(chrom=seqnames) %>% 
            dplyr::select(c("sample_id","chrom","start","end","cn_a","cn_b"))
         
        medicc2_input <- rbind(medicc2_input, remap_data)
    }
###
###
    write.table(medicc2_input,
        sprintf("%s/%s.medicc2.input.tsv", work_dir, BatchID),
        quote=F, col.names=T, row.names=F, sep="\t"
    )
###
###   
    ## medicc2 run
    RunMedicc2 <- "singularity exec -B /storage,/data /storage/images/medicc2.sif medicc2"
    RUN_CMD    <- sprintf("%s --input-type t --events --input-allele-columns cn_a,cn_b --normal-name diploid %s %s", 
        RunMedicc2, sprintf("%s/%s.medicc2.input.tsv", work_dir, BatchID), work_dir
    )
    system(RUN_CMD)
###
###
    # singularity exec -B /storage,/data /storage/images/medicc2.sif medicc2 \
    #     --input-type t \
    #     "/data/wgs/WGS-26-01/analysis/medicc2/WGS-26-01_02.medicc2.input.data.tsv" \
    #     "/data/wgs/WGS-26-01/analysis/medicc2" 



