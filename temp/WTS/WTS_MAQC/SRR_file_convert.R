

FastqDir <- "/storage/home/kangsm/analysis/gcx.rnaseq.maqc/fastq"
MetaDir  <- "/storage/home/kangsm/analysis/gcx.rnaseq.maqc/meta"

setwd(FastqDir)

data_list <- openxlsx::read.xlsx(sprintf("%s/srr_list.xlsx", MetaDir))

maqc_ids <- c("MAQC_II_A_rep1","MAQC_II_A_rep2","MAQC_II_B_rep1","MAQC_II_B_rep2")
maqc_ids <- c("MAQC_II_A_rep3","MAQC_II_B_rep3")
for( ID in maqc_ids )
{
    message(sprintf("|---> ID = %s", ID) )
    srr_files_list <- data_list %>% filter( SAMPLE_ID == ID )
    srr_fq_r1 <- paste(paste0(srr_files_list$SRR, "_1.fastq"), collapse=" ")
    srr_fq_r2 <- paste(paste0(srr_files_list$SRR, "_2.fastq"), collapse=" ")
    message("|---> R1 merging...")
    system(sprintf("cat %s  > %s_R1.fastq", srr_fq_r1, ID))
    message("|---> R2 merging...")
    system(sprintf("cat %s > %s_R2.fastq", srr_fq_r2, ID))
    message("|---> R1 compressing...")
    system(sprintf("bgzip --threads 20 %s_R1.fastq", ID))
    message("|---> R2 compressing...")
    system(sprintf("bgzip --threads 20 %s_R2.fastq", ID))
}

