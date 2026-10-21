
## R

    options(stringsAsFactors=FALSE) 
    suppressPackageStartupMessages({
        library("Hmisc")
        library("plyr")
        library("dplyr")
        library("setwidth")
        library("optparse")
    })
    BaseDir            <- "/data/wgs"
    BatchID            <- "WGS-26-01"
    work_dir           <- "/data/wgs/WGS-26-01/analysis/medicc2"
    analysis_list_file <- "/data/wgs/WGS-26-01/meta/analysis.list"
    sinfo <- read.table(sprintf("%s",analysis_list_file), header=F, sep=" ") %>% dplyr::rename(SeqID=1, Group=2, Order=3)
    tumor_samples <- sinfo %>% filter( Group != "normal" ) %>% .$SeqID

    cnv_data_list = lapply(  tumor_samples, function(tid) {
        purple_dir = sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)       
        cnv_data = read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",purple_dir,tid))        

        cnv_res <- cnv_data %>% mutate( cn_a=round(majorAlleleCopyNumber, 0), cn_b=round(minorAlleleCopyNumber, 0) ) %>% 
            mutate(sample_id=tid) %>% 
            dplyr::select( c("sample_id","chromosome","start","end","cn_a","cn_b") ) %>% 
            dplyr::rename( chrom=chromosome )
        return(cnv_res)
    })



    CNV_DATA = data.frame()

    for( tid in tumor_samples ){
        
        purple_dir = sprintf("%s/%s/%s/hmftools/purple", BaseDir, BatchID, tid)       
        cnv_data = read.delim(sprintf("%s/%s.purple.cnv.somatic.tsv",purple_dir,tid))        

        cnv_res <- cnv_data %>% mutate( cn_a=round(majorAlleleCopyNumber, 0), cn_b=round(minorAlleleCopyNumber, 0) ) %>% 
            mutate(sample_id=tid) %>% 
            dplyr::select( c("sample_id","chromosome","start","end","cn_a","cn_b") ) %>% 
            dplyr::rename( chrom=chromosome )

        CNV_DATA = rbind(CNV_DATA, cnv_res)
    }

    cn_gr <- GRanges(
        seqnames = CNV_DATA$chrom,
        ranges = IRanges(start = CNV_DATA$start + 1, end = CNV_DATA$end),
        sample = CNV_DATA$sample_id,
        cn_a = CNV_DATA$cn_a,
        cn_b = CNV_DATA$cn_b
    )

    bp_dt <- data.table::data.table(CNV_DATA)[, .(
        pos = c(start, end)
    ), by = chrom]

    bp_dt <- unique(bp_dt)
    data.table::setorder(bp_dt, chrom, pos)

common_dt <- bp_dt[, {
  x <- sort(unique(pos))
  if (length(x) < 2) {
    NULL
  } else {
    data.table(
      start = x[-length(x)],
      end = x[-1]
    )
  }
}, by = chrom]

common_dt <- common_dt[end > start]
setorder(common_dt, chrom, start, end)

common_gr <- GRanges(
  seqnames = common_dt$chrom,
  ranges = IRanges(start = common_dt$start + 1, end = common_dt$end)
)

message("Number of common segments: ", length(common_gr))

map_one_sample <- function(sample_name) {
  message("Mapping sample: ", sample_name)

  s_gr <- cn_gr[mcols(cn_gr)$sample == sample_name]

  hits <- findOverlaps(common_gr, s_gr, type = "within")

  if (length(hits) != length(common_gr)) {
    missing_idx <- setdiff(seq_along(common_gr), queryHits(hits))
    stop(sprintf(
      "Some common segments are not covered by sample %s. Example: %s:%s-%s",
      sample_name,
      as.character(seqnames(common_gr[missing_idx[1]])),
      start(common_gr[missing_idx[1]]) - 1,
      end(common_gr[missing_idx[1]])
    ))
  }

  # If multiple hits occur, choose the source segment with largest overlap.
  hit_dt <- data.table(
    q = queryHits(hits),
    s = subjectHits(hits),
    width = width(pintersect(common_gr[queryHits(hits)], s_gr[subjectHits(hits)]))
  )

  setorder(hit_dt, q, -width)
  hit_dt <- hit_dt[, .SD[1], by = q]

  data.table(
    sample = sample_name,
    chrom = as.character(seqnames(common_gr[hit_dt$q])),
    start = start(common_gr[hit_dt$q]) - 1,
    end = end(common_gr[hit_dt$q]),
    major_cn = mcols(s_gr)$major_cn[hit_dt$s],
    minor_cn = mcols(s_gr)$minor_cn[hit_dt$s]
  )
}


    write.table(CNV_DATA,
        sprintf("%s/%s.medicc2.input.data.tsv", work_dir, tid),
        quote=F, col.names=T, row.names=F, sep="\t"
    )

    
    ## medicc2 run

    singularity exec -B /storage,/data /storage/images/medicc2.sif medicc2 \
        --input-type t \
        "/data/wgs/WGS-26-01/analysis/medicc2/WGS-26-01_02.medicc2.input.data.tsv" \
        "/data/wgs/WGS-26-01/analysis/medicc2" 

