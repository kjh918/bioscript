

# SRR donwload

### 1st 
```
FASTQ_DIR=/storage/home/kangsm/analysis/gcx.rnaseq.maqc
SRR_LIST=$(cat /storage/home/kangsm/analysis/gcx.rnaseq.maqc/download_srr.list | awk -F " " '{print $1}')
parallel -j 12 -k "fastq-dump --split-3 --outdir ${FASTQ_DIR} {}" ::: ${SRR_LIST}
```

### 2nd 
```
FASTQ_DIR=/storage/home/kangsm/analysis/gcx.rnaseq.maqc/fastq
SRR_LIST=$(cat /storage/home/kangsm/analysis/gcx.rnaseq.maqc/meta/download_srr2.list | awk -F " " '{print $1}')
parallel -j 32 -k "fastq-dump --split-3 --outdir ${FASTQ_DIR} {}" ::: ${SRR_LIST}
```

### 3rd
```
FASTQ_DIR=/storage/home/kangsm/analysis/gcx.rnaseq.maqc/fastq
SRR_LIST=$(cat /storage/home/kangsm/analysis/gcx.rnaseq.maqc/meta/download_srr3.list | awk -F " " '{print $1}')
parallel -j 32 -k "fastq-dump --split-3 --outdir ${FASTQ_DIR} {}" ::: ${SRR_LIST}
```

## SINGLE END PROCESSING
SE_SAMPLE_N=$(cat /storage/home/kangsm/analysis/gcx.rnaseq.maqc/single-end.list | wc -l)
for N in seq 1 $SE_SAMPLE_N
do
    SRR_INFO=$(sed -n -e "${N} p" /storage/home/kangsm/analysis/gcx.rnaseq.maqc/single-end.list)
    SRR=$(echo ${SRR_INFO} | awk -F " " '{print $1}')
    FQ_NAME=$(echo ${SRR_INFO} | awk -F " " '{print $3}')
done

for SRR in ${SRR_LIST}
do
    mv ${FASTQ_DIR}/${SRR}_1.fastq ${FASTQ_DIR}/${SRR}_R1.fastq 
    mv ${FASTQ_DIR}/${SRR}_2.fastq ${FASTQ_DIR}/${SRR}_R2.fastq 
    bgzip --threads 15 ${FASTQ_DIR}/${SRR}_R1.fastq 
    bgzip --threads 15 ${FASTQ_DIR}/${SRR}_R2.fastq 
done

# MAQC-A, MAQC-B
    cat SRR3670989.fastq SRR3670990.fastq SRR3670991.fastq SRR3670992.fastq > MAQC_B_rep2.fastq
    cat SRR3670985.fastq SRR3670986.fastq SRR3670987.fastq SRR3670988.fastq > MAQC_B_rep1.fastq
    cat SRR3670981.fastq SRR3670982.fastq SRR3670983.fastq SRR3670984.fastq > MAQC_A_rep2.fastq
    cat SRR3670977.fastq SRR3670978.fastq SRR3670979.fastq SRR3670980.fastq > MAQC_A_rep1.fastq

    bgzip --threads 15 MAQC_A_rep1.fastq
    bgzip --threads 15 MAQC_A_rep2.fastq
    bgzip --threads 15 MAQC_B_rep1.fastq
    bgzip --threads 15 MAQC_B_rep2.fastq

# MAQC-C
    cat SRR896663_1.fastq SRR896665_1.fastq SRR896667_1.fastq SRR896669_1.fastq SRR896671_1.fastq SRR896673_1.fastq SRR896675_1.fastq SRR896677_1.fastq > MAQC_C_F1_rep1_R1.fastq
    cat SRR896663_2.fastq SRR896665_2.fastq SRR896667_2.fastq SRR896669_2.fastq SRR896671_2.fastq SRR896673_2.fastq SRR896675_2.fastq SRR896677_2.fastq > MAQC_C_F1_rep1_R2.fastq
    cat SRR896664_1.fastq SRR896666_1.fastq SRR896668_1.fastq SRR896670_1.fastq SRR896672_1.fastq SRR896674_1.fastq SRR896676_1.fastq SRR896678_1.fastq > MAQC_C_F2_rep1_R1.fastq
    cat SRR896664_2.fastq SRR896666_2.fastq SRR896668_2.fastq SRR896670_2.fastq SRR896672_2.fastq SRR896674_2.fastq SRR896676_2.fastq SRR896678_2.fastq > MAQC_C_F2_rep1_R2.fastq
    cat SRR896679_1.fastq SRR896681_1.fastq SRR896683_1.fastq SRR896685_1.fastq SRR896687_1.fastq SRR896689_1.fastq SRR896691_1.fastq SRR896693_1.fastq > MAQC_C_F1_rep2_R1.fastq
    cat SRR896679_2.fastq SRR896681_2.fastq SRR896683_2.fastq SRR896685_2.fastq SRR896687_2.fastq SRR896689_2.fastq SRR896691_2.fastq SRR896693_2.fastq > MAQC_C_F1_rep2_R2.fastq
    cat SRR896680_1.fastq SRR896682_1.fastq SRR896684_1.fastq SRR896686_1.fastq SRR896688_1.fastq SRR896690_1.fastq SRR896692_1.fastq SRR896694_1.fastq > MAQC_C_F2_rep2_R1.fastq
    cat SRR896680_2.fastq SRR896682_2.fastq SRR896684_2.fastq SRR896686_2.fastq SRR896688_2.fastq SRR896690_2.fastq SRR896692_2.fastq SRR896694_2.fastq > MAQC_C_F2_rep2_R2.fastq

    bgzip --threads 15 MAQC_C_F1_rep1_R1.fastq 
    bgzip --threads 15 MAQC_C_F1_rep1_R2.fastq
    bgzip --threads 15 MAQC_C_F2_rep1_R1.fastq
    bgzip --threads 15 MAQC_C_F2_rep1_R2.fastq
    bgzip --threads 15 MAQC_C_F1_rep2_R1.fastq
    bgzip --threads 15 MAQC_C_F1_rep2_R2.fastq
    bgzip --threads 15 MAQC_C_F2_rep2_R1.fastq
    bgzip --threads 15 MAQC_C_F2_rep2_R2.fastq



