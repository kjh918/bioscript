
library(scarHRD)
library(data.table)

## https://github.com/sztup/scarHRD

SEQZ_RES = "/storage/home/kangsm/analysis/gcx.MRD/05.VALIDATION/VSID02_Colon.Cancer.Organoid/CRC02/CRC02_TUMOR_P100/cnv/sequenza/CRC02_TUMOR_P100_bin50.seqz.gz"

hrd = scar_score(SEQZ_RES, reference = "grch37", seqz=T, chr.in.names=F)

seqz_res = fread(SEQZ_RES)
SEQZ_SEGMENT = "/storage/home/kangsm/analysis/gcx.MRD/05.VALIDATION/VSID02_Colon.Cancer.Organoid/CRC02/CRC02_TUMOR_P100/cnv/sequenza/CRC02_TUMOR_P100_segments.txt"
seqz_segment = fread(SEQZ_SEGMENT)
# BRCA1
cnv_brca1 = seqz_segment %>% filter( chromosome == 17 , start.pos <= 41196311, end.pos >= 41277500 )
#seqz_res %>% filter( chromosome == 17 , position >= 41196311, position <= 41277500 )
# BRCA2
cnv_brca2 = seqz_segment %>% filter( chromosome == 13 , start.pos <= 32889610, end.pos >= 32973805  )
#seqz_res %>% filter( chromosome == 13 , position >= 32889610, position <= 32973805 )
BRCA_CNV_RES = rbind(
    cnv_brca1 %>% mutate( Gene = "BRCA1", Gene_Region = "chr17:41196311-41277500"),
    cnv_brca2 %>% mutate( Gene = "BRCA2", Gene_Region = "chr13:32889610-32973805")
) %>% dplyr::select(c("Gene","Gene_Region","chromosome","start.pos","end.pos","Bf","depth.ratio","CNt","A","B")) %>% dplyr::rename(
    Chrom=chromosome, Segment_Start=start.pos, Segment_End=end.pos, B_allele_Freq=Bf, logRatio=depth.ratio, Total_CN= CNt, CN_A=A, CN_B=B
) %>% mutate( B_allele_Freq = round(B_allele_Freq, 4), logRatio = round(log2(logRatio),4) )

### 
VARS = read.delim("/storage/home/kangsm/analysis/gcx.MRD/05.VALIDATION/VSID02_Colon.Cancer.Organoid/CRC02/CRC02_TUMOR_P100/vcf/CRC02_TUMOR_P100.mutect2.NT.bias.filtered.vep.annotated.maf" )

VARS_BRCA = VARS %>% filter( Hugo_Symbol %in% c("BRCA1","BRCA2")) %>% filter( BIAS_STATUS != "artifact") %>% mutate( VAF = round(TUMOR_ALT_COUNT/(TUMOR_ALT_COUNT+TUMOR_REF_COUNT), 3) )
VARS_BRCA = VARS_BRCA[,c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","TUMOR_REF_COUNT","TUMOR_ALT_COUNT","Variant_Classification","Consequence","VAF",
    "HGVSc","HGVSp","HGVSp_Short","HGVSg","FILTER","BIAS_STATUS","Refseq_TRANSCRIPT","Ensembl_TRANSCRIPT","Refseq_EXON","Ensembl_EXON")] %>% dplyr::rename(Mutect2_FILTER=FILTER, SOBDetector_BiasStatus=BIAS_STATUS)

BRCA_VARIANTS_RES = rbind(
    VARS_BRCA %>% filter( Hugo_Symbol == "BRCA1" ),
    VARS_BRCA %>% filter( Hugo_Symbol == "BRCA2" )
)



