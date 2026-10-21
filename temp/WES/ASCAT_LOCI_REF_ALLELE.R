

    indexing = function(X) { sapply(unique(X), function(y) list(y)) }

ref_dir = "/storage/references_and_index/cnv/ascat/G1000_lociAll_hg19"

CHRS = paste0("chr", c(1:22, "X"))

ref_allele_query = function(chr, pos)
{
    Q = sprintf("samtools faidx /storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta %s:%s-%s", chr, pos, pos)
    query_allele = system(Q, intern=T)[2]
    return(query_allele)
}

    loci = read.delim(sprintf("%s/G1000_loci_hg19_%s.txt", ref_dir, CHRS[1]), header=FALSE) %>%
        dplyr::rename(chr=1, pos=2) %>% mutate(index=paste(chr, pos, sep="_"))
    
    idx =indexing(loci$index)
    




ref_allele_query(1, 10642)
