

WORK_DIR=/storage/home/kangsm/analysis/cnv_test/adtex_saascnv
GENOME_FASTA=/storage/references_and_index/hg19/fasta/human_g1k_v37_decoy.fasta

RUN_GATK="singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk"

VCF_FILE=/storage/home/kangsm/analysis/cnv_test/vcfs/PON01_NT_01.merged.vcf.gz

N_ID="coad.normal.01"
T_ID="coad.tumor.01"


${RUN_GATK} SelectVariants \
    --reference ${GENOME_FASTA} \
    --variant ${VCF_FILE} \
    -sn ${N_ID} -sn ${T_ID} \
    --select-type-to-include SNP --restrict-alleles-to BIALLELIC \
    -select "vc.getGenotype('$N_ID').getDP() > 9 && vc.getGenotype('$T_ID').getDP() > 19" \
    --output ${WORK_DIR}/snp_trimmed.vcf

## R
library(saasCNV)
library(plyr)
library(dplyr)
library(strinr)
vcf.file="snp_trimmed.vcf"





