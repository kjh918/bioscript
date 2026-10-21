
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(VariantAnnotation)
})

BatchID="WGS-26-01"
SeqID="WGS-26-01_01"

# SNV
var_list <- lapply( c("WGS-26-01_01","WGS-26-01_02"), function(SeqID) {
    vcf  <- readVcf(sprintf("/data/wgs/RUO_Service/WGS-26-01/%s/hmftools/purple/%s.purple.somatic.vcf.gz", SeqID, SeqID))
    rr   <- rowRanges(vcf)
    geno <- geno(vcf)
    ad <- geno$AD[, SeqID, drop = FALSE]
    ad_vec <- ad[, 1]
    ref_count <- vapply(ad_vec, function(x) {
        if (length(x) < 2 || any(is.na(x[1:2]))) return(NA)
        as.numeric(x[1])
    }, numeric(1))
    alt_count <- vapply(ad_vec, function(x) {
        if (length(x) < 2 || any(is.na(x[1:2]))) return(NA)
        as.numeric(x[2])
    }, numeric(1))
    total_count <- ref_count + alt_count
    var_data    <- data.table(
        sample = SeqID,
        chr = as.character(seqnames(rr)),
        pos = start(rr),
        ref = as.character(ref(vcf)),
        alt = sapply(alt(vcf), function(x) as.character(x[1])),
        R   = alt_count,
        X   = total_count,
        VAF_calc = alt_count / total_count
    )
    var_data <- var_data[nchar(ref) == 1 & nchar(alt) == 1]
    var_data <- var_data[X >= 10 & R >= 2]
    return(var_data)
})

# CNV
cnv_list <- lapply( c("WGS-26-01_01","WGS-26-01_02"), function(SeqID) {
    cnv  <- fread(sprintf("/data/wgs/RUO_Service/WGS-26-01/%s/hmftools/purple/%s.purple.cnv.somatic.tsv", SeqID, SeqID))
    qc   <- fread(sprintf("/data/wgs/RUO_Service/WGS-26-01/%s/hmftools/purple/%s.purple.purity.tsv", SeqID, SeqID))
    out <- data.table(
        sample = SeqID,
        chr    = as.character(cnv[["chromosome"]]),
        start  = as.integer(cnv[["start"]]),
        end    = as.integer(cnv[["end"]]),
        major_cn = as.numeric(cnv[["majorAlleleCopyNumber"]]),
        minor_cn = as.numeric(cnv[["minorAlleleCopyNumber"]])
    )
    out[, len := end - start + 1]
    out[, total_cn := major_cn + minor_cn]
    out[, ploidy := qc$ploidy]
    out <- out[len >= 1e6 & total_cn <= qc$ploidy * 3 ]
    return(out)
})

##
#
all_snv <- rbindlist(var_list) %>% mutate(mut_id=paste(chr,pos,ref,alt,sep=":"))
all_cnv <- rbindlist(cnv_list) %>% mutate(cna_id=paste(chr,start,end,sep=":"))
sample_list <- c("WGS-26-01_01","WGS-26-01_02")

all_snv <- all_snv[ chr %in% as.character(1:22) ]
all_cnv <- all_cnv[ chr %in% as.character(1:22) ]

mut_ids <- sort(unique(all_snv$mut_id)) 
R <- matrix(0, nrow = length(mut_ids), ncol = length(sample_list), dimnames = list(mut_ids, sample_list))
X <- matrix(0, nrow = length(mut_ids), ncol = length(sample_list), dimnames = list(mut_ids, sample_list))
for (s in sample_list) {
    z <- all_snv[sample == s]
    R[z$mut_id, s] <- z$R
    X[z$mut_id, s] <- z$X
}

#
cna_ids <- sort(unique(all_cnv$cna_id))
WM <- matrix(1, nrow = length(cna_ids), ncol = length(sample_list), dimnames = list(cna_ids, sample_list))
Wm <- matrix(1, nrow = length(cna_ids), ncol = length(sample_list), dimnames = list(cna_ids, sample_list))
for (s in sample_list) {
    z <- all_cnv[sample == s]
    WM[z$cna_id, s] <- z$major_cn
    Wm[z$cna_id, s] <- z$minor_cn
}

#
mut_dt <- unique(all_snv[, .(mut_id, chr, pos)], by = "mut_id")
mut_dt <- mut_dt[match(mut_ids, mut_id)]

cna_dt <- unique(all_cnv[, .(cna_id, chr, start, end)], by = "cna_id")
cna_dt <- cna_dt[match(cna_ids, cna_id)]

gr_mut <- GRanges(mut_dt$chr, IRanges(start = mut_dt$pos,   end = mut_dt$pos), mut_id = mut_dt$mut_id )
gr_cna <- GRanges(cna_dt$chr, IRanges(start = cna_dt$start, end = cna_dt$end), cna_id = cna_dt$cna_id )

hits  <- findOverlaps(gr_mut, gr_cna)
n_sna <- nrow(mut_dt)
n_cna <- nrow(cna_dt)
canopy_Y <- matrix(0, nrow=n_sna, ncol=n_cna+1, dimnames = list(mut_dt$mut_id, c("noCNA", cna_dt$cna_id)))  
qh <- S4Vectors::queryHits(hits)      # SNA index
sh <- S4Vectors::subjectHits(hits)    # CNA-region index
if (length(hits) > 0) { canopy_Y[cbind(qh, sh + 1L)] <- 1 }
no_overlap <- rowSums(canopy_Y[, -1, drop = FALSE]) == 0
canopy_Y[no_overlap, "noCNA"] <- 1L
canopy_C <- diag(n_cna)
rownames(canopy_C) <- cna_dt$cna_id
colnames(canopy_C) <- cna_dt$cna_id


cnv_meta <- copy(all_cnv)
cnv_meta[, rel_cn := total_cn / ploidy]
cnv_meta[, loh_like := minor_cn < 0.5]
cnv_meta[, high_amp_like := rel_cn >= 3]
cnv_meta[, broad_5mb := len >= 5e6]
cnv_meta[, broad_10mb := len >= 1e7]
cnv_meta[, delta_from_ploidy := abs(total_cn - ploidy)]
cnv_meta[, allele_imbalance := abs(major_cn - minor_cn)]
cnv_meta[, score := 3.0 * loh_like + 2.0 * broad_10mb + 1.0 * broad_5mb + delta_from_ploidy + 0.5 * allele_imbalance ]
cna_core <- cnv_meta[ 
    !high_amp_like & ( 
        broad_10mb | (broad_5mb & loh_like) | (broad_5mb & delta_from_ploidy >= 0.75)
    ) 
][order(-score)]
cna_core_balanced <- cna_core[ , .SD[order(-loh_like, -broad_10mb, -score)][1:min(.N, 20)], by = sample ]
keep_cna <- unique(cna_core_balanced$cna_id)
keep_cna <- keep_cna[keep_cna %in% rownames(WM)]


noCNA_col <- colnames(canopy_Y)[1]
WM2 <- WM[keep_cna, , drop = FALSE]
Wm2 <- Wm[keep_cna, , drop = FALSE]
Y_tmp <- canopy_Y[, c(noCNA_col, keep_cna), drop = FALSE]
C2 <- diag(length(keep_cna))
rownames(C2) <- keep_cna
colnames(C2) <- keep_cna


auto_mut_ids <- unique(all_snv$mut_id)
auto_mut_ids <- auto_mut_ids[auto_mut_ids %in% rownames(R)]

R_auto <- R[auto_mut_ids, , drop = FALSE]
X_auto <- X[auto_mut_ids, , drop = FALSE]
Y_auto <- Y_tmp[auto_mut_ids, , drop = FALSE]

vaf <- R_auto / X_auto
vaf[!is.finite(vaf)] <- NA

sna_meta <- data.table(
  mut_id = rownames(R_auto),
  depth_min = apply(X_auto, 1, min, na.rm = TRUE),
  depth_mean = rowMeans(X_auto, na.rm = TRUE),
  alt_max = apply(R_auto, 1, max, na.rm = TRUE),
  vaf_mean = rowMeans(vaf, na.rm = TRUE),
  vaf_max = apply(vaf, 1, max, na.rm = TRUE),
  vaf_min = apply(vaf, 1, min, na.rm = TRUE)
)

sna_meta[, vaf_diff := abs(vaf[, 1] - vaf[, 2])]
sna_meta[, cna_overlap := rowSums(Y_auto[, -1, drop = FALSE])]
sna_meta[, region := ifelse(cna_overlap > 0, "CNA_overlap", "noCNA")]

sna_pool <- sna_meta[
  depth_min >= 30 &
    alt_max >= 5 &
    vaf_max >= 0.03 &
    vaf_max <= 0.90
]

sna_pool[, .N, by = region]

sna_pool[, score :=
  2.0 * vaf_diff +
    0.5 * abs(vaf_mean - 0.25) +
    0.2 * log10(depth_mean + 1)
]

keep_sna_noCNA <- sna_pool[
  region == "noCNA"
][order(-score)]$mut_id

keep_sna_CNA <- sna_pool[
  region == "CNA_overlap"
][order(-score)]$mut_id

keep_sna_noCNA <- head(keep_sna_noCNA, 100)
keep_sna_CNA <- head(keep_sna_CNA, 50)

keep_mut <- unique(c(keep_sna_noCNA, keep_sna_CNA))
length(keep_mut)

R2 <- R_auto[keep_mut, , drop = FALSE]
X2 <- X_auto[keep_mut, , drop = FALSE]
Y2 <- Y_auto[keep_mut, , drop = FALSE]

##
library(Canopy)
sampchain_test <- canopy.sample(
    R = R2,
    X = X2,
    WM = WM2,
    Wm = Wm2,
    epsilonM = 0.1,
    epsilonm = 0.1,
    C = C2,
    Y = Y2,
    K = 2:5,
    numchain = 4,
    max.simrun = 2000,
    min.simrun = 500,
    writeskip = 100,
    projectname = "WGS-26-01_canopy_smoke"
)
bic_test <- canopy.BIC(
  sampchain = sampchain_test,
  projectname = "WGS-26-01_canopy_smoke",
  K = 2:5,
  numchain = 4,
  burnin = 1,
  thin = 1
)
optK <- c(2:5)[which.min(bic_test)]
post_test <- canopy.post(
    sampchain = sampchain_test,
    projectname = "WGS-26-01_canopy_smoke",
    K = 2:5,
    numchain = 4,
    burnin = 1,
    thin = 1,
    optK = optK,
    C = C2
)
tree_test <- canopy.output(
  post = post_test,
  config.i = 1,
  C = C2
)
canopy.plottree(
  tree = tree_test,
  pdf = TRUE,
  pdf.name = "WGS-26-01_canopy_smoke_tree.pdf"
)

sampchain <- canopy.sample(
    R = R,
    X = X,
    WM = WM,
    Wm = Wm,
    epsilonM = 0.05,
    epsilonm = 0.05,
    C = canopy_C,
    Y = canopy_Y,
    K = 3:10,
    numchain = 20,
    max.simrun = 50000,
    min.simrun = 10000,
    writeskip = 200,
    projectname = "WGS-26-01_hmftools.canopy"
)


    canopy_input_data_list <- list(
        R = R,
        X = X,
        WM = WM,
        Wm = Wm,
        epsilonM = epsilon,
        epsilonm = epsilon,
        Y = Y,
        C = NULL,
        sna.name = mut_ids,
        cna.name = cna_ids,
        snv_table = mut_dt,
        cna_table = cna_dt
    )