

library(data.table)
library(plyr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)

mut_hek293 = fread("HEK293_mutations.csv")
mut_mcf7   = fread("MCF7_mutations.csv")
genes <- unique(c(mut_hek293$Gene,mut_mcf7$Gene))

gtf_file <- sprintf("/storage/references_and_index/RNAseq/%s/genes/genes.gtf", "hg38_ercc.plus")
gtf <- rtracklayer::import(gtf_file) %>% as.data.frame() %>% data.table()

target_genes <- gtf %>% filter( gene_name %in% genes, type == "exon")
target_genes$start = target_genes$start-50
target_genes$end   = target_genes$end+50

write.table(target_genes[,c("seqnames","start","end")], "mutation_gene_regions_raw.bed", quote=F, col.names=F, row.names=F, sep="\t")


cat mutation_gene_regions_raw.bed | sort -k1,1 -k2,2n |
    bedtools merge -i - -c 3 -o distinct |
    bedtools slop -i - -g /storage/references_and_index/RNAseq/hg38_ercc.plus/fasta/genome.chrom.sizes -b 50 |
    bedtools merge -i - -c 4 -o distinct \
    > mutation_gene_regions.bed

mutation_gene_regions.bed


SEQ_DICT=/storage/references_and_index/RNAseq/hg38_ercc.plus/fasta/genome.dict
java -jar /storage/apps/bin/picard.jar BedToIntervalList \
    --INPUT /data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/mutation_gene_regions.bed \
    --OUTPUT /data/scRNAseq/02.Public-Datasets/004.GSE151334.PRJNA635600.cell-line-smartseq2/meta/mutation_gene_regions.interval_list \
    --SEQUENCE_DICTIONARY ${SEQ_DICT}


