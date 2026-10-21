
#---/ OptiType Test Run (2023.09.23)

# USE '23WES02' DATA 

FQ_RAW_DIR=/BiO/DATA/WES/23WES02/fastq_merged
SAMPLE_LIST=$(cat /BiO/DATA/WES/23WES02/sample.list2)

BASE_DIR=/BiO/DATA/WES/23WES02/extra_analysis/HLA
HLA_FQ_DIR=${BASE_DIR}/fastq_hla
HLA_RES_DIR=${BASE_DIR}/optitype_results
HLA_FASTA=/BiO/Install/optitype-1.3.3/data/hla_reference_dna.fasta

start pm 6:40
#---/ HLA-Fasta Selection
    parallel -j 6 -k "razers3 -i 95 -m 1 -dr 0 \
        -o ${HLA_FQ_DIR}/{}_HLA_R1.bam ${HLA_FASTA} ${FQ_RAW_DIR}/{}_R1.fastq.gz" ::: ${SAMPLE_LIST}

    parallel -j 6 -k "razers3 -i 95 -m 1 -dr 0 \
        -o ${HLA_FQ_DIR}/{}_HLA_R2.bam ${HLA_FASTA} ${FQ_RAW_DIR}/{}_R2.fastq.gz" ::: ${SAMPLE_LIST}
#---/ BAM to Fastq
    parallel -j 6 -k "/BiO/SWTools/samtools bam2fq ${HLA_FQ_DIR}/{}_HLA_R1.bam > ${HLA_FQ_DIR}/{}_HLA_R1.fastq" ::: ${SAMPLE_LIST}
    parallel -j 6 -k "/BiO/SWTools/samtools bam2fq ${HLA_FQ_DIR}/{}_HLA_R2.bam > ${HLA_FQ_DIR}/{}_HLA_R2.fastq" ::: ${SAMPLE_LIST}

#---/ Run Optitype
    parallel -j 6 -k "python /BiO/Install/miniconda3/envs/HLA/bin/OptiTypePipeline.py \
        --input ${HLA_FQ_DIR}/{}_HLA_R1.fastq ${HLA_FQ_DIR}/{}_HLA_R2.fastq \
        --dna \
        --outdir ${HLA_RES_DIR} \
        --prefix {}" ::: ${SAMPLE_LIST}


#---/ HLA-LA Test Run (2023.09.24)

SAMPLE_LIST=$(cat /BiO/DATA/WES/23WES02/sample.list2)
BAM_DIR=/BiO/DATA/WES/23WES02/bam_files

parallel -j 6 -k "/BiO/Install/HLA-LA-1.0.3/src/HLA-LA.pl \
    --BAM ${BAM_DIR}/{}.sorted.dup.marked.recal.bam \
    --graph PRG_MHC_GRCh38_withIMGT \
    --sampleID {} \
    --maxThreads 25 \
    --workingDir /BiO/DATA/WES/23WES02/extra_analysis/HLA/hla-la_results" ::: ${SAMPLE_LIST}





