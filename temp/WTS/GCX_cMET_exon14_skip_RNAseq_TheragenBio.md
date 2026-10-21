
# GCX cMET exon14 Skipping Sample RNAseq

### 1. 데이터 다운로드 ( from 테라젠바이오 FTP )

<pre style="font-family:Consolas">
# Truseq RNA Exome kit Raw data&기본분석 결과
wget -r ftp://"GCRX_GMR":"genCRX_00736!"@ftp.theragenbio.com/TBD231186_18685_20231020

# SMARTer picomamalian kit 기본분석 결과
wget -r ftp://"GCRX_GMR":"genCRX_00736!"@ftp.theragenbio.com/TBD231186_18832_20231020

# SMARTer rawdata
wget -r ftp://"GCRX_GMR":"genCRX_00736!"@ftp.theragenbio.com/TBD231186_18832_20231013
</pre>

---

## 2. 데이터 프로세싱
### 데이터 프로세싱 개요
- FFPE 샘플인것을 감안하여 fastq 및 aligned BAM QC를 우선함.
- Genome Alignment 와 `7번 염색체`만 reference로 사용하는 방식을 모두 진행
- hg19 먼저 프로세싱하고 hg38도 수행

> `cMET 유전자 위치` : 7번 염색체 ( 7q21-q31 )  
---
#### A. chr7 레퍼런스 파일 생성
1. Fasta File   
<pre style="font-family:Consolas">
# hg19 chr7 fasta
samtools faidx /storage/references_and_index/gcx.hg19/gcx.hg19.fasta 1 > /data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.fasta

# hg38 chr7 fasta
samtools faidx /storage/references_and_index/gcx.hg38/gcx.hg38.fa chr1 > /data/wts/gcx-cmet-exon14/ref/gcx.hg38.chr7.fasta
</pre>

2. GFT File 
<pre style="font-family:Consolas">
# hg19 chr7 GTF
java -jar /storage/home/kangsm/swtools/jvarkit/dist/gtfsplitter.jar \
    -m contig \
    -o /data/wts/gcx-cmet-exon14/ref/hg19_gtf_split \
    /storage/references_and_index/gcx.hg19/gcx.hg19.gtf 

# hg38 chr7 GTF
java -jar /storage/home/kangsm/swtools/jvarkit/dist/gtfsplitter.jar \
    -m contig \
    -o /data/wts/gcx-cmet-exon14/ref/hg38_gtf_split \
    /storage/references_and_index/gcx.hg38/gcx.hg38.gtf 
</pre>

3. RSEM-STAR Index  
<pre style="font-family:Consolas">
# chr7 rsem-index generation
# hg19
/storage/apps/rsem-1.3.3/bin/rsem-prepare-reference \
    --num-threads 20 \
    --gtf /data/wts/gcx-cmet-exon14/ref/hg19_gtf_split/7.gtf \
	--star --star-path /storage/apps/STAR-2.7.11a/Linux_x86_64 \
    --star-sjdboverhang 150 \
   	/data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.fasta \
    /data/wts/gcx-cmet-exon14/ref/rsem-index-hg19-chr7/gcx.hg19.chr7
#hg38
/storage/apps/rsem-1.3.3/bin/rsem-prepare-reference \
    --num-threads 20 \
    --gtf /data/wts/gcx-cmet-exon14/ref/hg38_gtf_split/chr7.gtf \
	--star --star-path /storage/apps/STAR-2.7.11a/Linux_x86_64 \
    --star-sjdboverhang 150 \
   	/data/wts/gcx-cmet-exon14/ref/gcx.hg38.chr7.fasta \
    /data/wts/gcx-cmet-exon14/ref/rsem-index-hg38-chr7/gcx.hg38.chr7
</pre>
---

#### B. 데이터 정보 & 데이터 위치 및 이름 재설정

1. 데이터 정보 
 
| 테라젠 폴더이름 |  생성날짜 | 라이브러리 | 데이터 사이즈<br>(압축파일) | Read Length | Note |
|:---:|:---:|:---:|---:|:---:|:---|
|TBD231186_18832_20231013|  2023-10-13 | SMARTer<br>picomamalian kit | 2.6 GB x 2 | PE-150 | R2 = sense strand |
|TBD231186_18832_20231020|  2023-10-19 | SMARTer<br>picomamalian kit | 2.6 GB x 2 | PE-150 |R2 = sense strand |
|TBD231186_18685_20231020|  2023-09-22 | Illumina<br>TruSeq RNA Exome kit | 2.4 GB x 2 | PE-150 | R2 = sense strand |
- 모두 Paired-End 150bp (`PE-150`) 방식 시퀀싱
- 두 라이브러리 모두 stranded-RNAseq 방식임 >> `Read 2 = sense strand` ( = mRNA starnd )
  


---
- fastq 파일들을  ${SEQ_ID}.fastq.gz 로 변경

- SEQ_FOLDER=`/data/wts/gcx-cmet-exon14`

| 테라젠 폴더이름 | SEQ_ID | 생성날짜 | Library | Data Size | Read Length | Note |
|:---:|:---:|:---:|:---:|---:|---:|:---|
|TBD231186_18832_20231013| gcx_cmet_smarter_1013 | 2023-10-13 | SMARTer<br>picomamalian kit | 2.6G x 2 | PE-150 | R2 = sense strand |
|TBD231186_18832_20231020| gcx_cmet_smarter_1020 | 2023-10-19 | SMARTer<br>picomamalian kit | 2.6G x 2 | PE-150 |R2 = sense strand |
|TBD231186_18685_20231020| gcx_cmet_truseq       | 2023-09-22 | TruSeq RNA Exome kit | 2.4G x 2 | PE-150 | R2 = sense strand |


---






---
---
### fastQC : raw-fastq 품질 검사
1. fastQC Run
<pre style="font-family:Consolas">
SEQ_FOLDER=/data/wts/gcx-cmet-exon14
for SEQ_ID in gcx_cmet_truseq gcx_cmet_smarter_1020 gcx_cmet_smarter_1013
do
    /storage/apps/bin/fastqc --noextract -t 2 \
        -o ${SEQ_FOLDER}/${SEQ_ID}/qcfiles \
        ${SEQ_FOLDER}/${SEQ_ID}/fastq/*.fastq.gz &
done
</pre>



---

### fastQC 결과
### A. TruSeq Reads
<img src="../img/gcx_cmet_truseq_r1_fastqc.png" width="400px" height="200px" title="gcx_cmet_truseq_R1"/>   
<img src="../img/gcx_cmet_truseq_r2_fastqc.png" width="400px" height="200px" title="gcx_cmet_truseq_R2"/>   
---

### B. SMARTer Picomammalian
<img src="../img/gcx_cmet_smarter_R1.png" width="400px" height="200px" title="gcx_cmet_smarter_R1"/>  
<img src="../img/gcx_cmet_smarter_R2.png" width="400px" height="200px" title="gcx_cmet_smarter_R2"/>   
---


| Sample | Read Length | GC (%) | Total Sequences |
|---:|:---:|:---:|---:|
| TruSeq-R1 | 151 | 48 | 32,226,503 |
| TruSeq-R2 | 151 | 49 | 32,226,503 |
| SMARTer-R1 | 151 | 60 | 36,834,596 |
| SMARTer-R2 | 151 | 61 | 36,834,596 |
---

- `gcx_cmet_smarter_1013`와 `gcx_cmet_smarter_1020`은 `IDENTICAL` 샘플   
    
    ==> 하나만 진행 = `gcx_cmet_smarter_1020`

> SEQ_ID 1 = `gcx_cmet_truseq`  
> SEQ_ID 2 = `gcx_cmet_smarter` <- gcx_cmet_smarter_1020에서 이름 변경

---

### fastq trimming & filtering 
<pre style="font-family:Consolas">
SEQ_IDS="gcx_cmet_truseq gcx_cmet_smarter"
parallel -j 4 -k "/storage/apps/bin/fastp \
    --in1 /data/wts/gcx-cmet-exon14/{}/fastq/{}_R1.fastq.gz \
    --in2 /data/wts/gcx-cmet-exon14/{}/fastq/{}_R2.fastq.gz \
    --html /data/wts/gcx-cmet-exon14/{}/qcfiles/{}.fastp.html \
    --json /data/wts/gcx-cmet-exon14/{}/qcfiles/{}.fastp.json \
    --length_required 50 \
    --qualified_quality_phred 20 \
    --out1 /data/wts/gcx-cmet-exon14/{}/fastq/{}.trimmed_R1.fastq.gz \
    --out2 /data/wts/gcx-cmet-exon14/{}/fastq/{}.trimmed_R2.fastq.gz" ::: ${SEQ_IDS}
</pre>

---
## STAR ALIGNENT

### (1) hg19 + 2-pass basic + genome ( rsem-index-star-index )

<pre style="font-family:Consolas">
SEQ_IDS="gcx_cmet_truseq gcx_cmet_smarter"

parallel -j 2 -k "/storage/apps/bin/STAR \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --twopassMode Basic \
    --genomeLoad NoSharedMemory \
    --runThreadN 40 \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --readQualityScoreBase 20 \
    --outFilterMismatchNmax 30 \
    --outFilterMismatchNoverLmax 0.1 \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped None\
    --genomeDir /storage/references_and_index/gcx.hg19/rsem-index-star-pe150 \
    --readFilesIn /data/wts/gcx-cmet-exon14/{}/fastq/{}.trimmed_R1.fastq.gz /data/wts/gcx-cmet-exon14/{}/fastq/{}.trimmed_R2.fastq.gz \
    --outFileNamePrefix /data/wts/gcx-cmet-exon14/{}/bam/{}. \
    --sjdbGTFfile /storage/references_and_index/gcx.hg19/gcx.hg19.gtf " ::: ${SEQ_IDS}
</pre>

    
### (2) hg19 + 2-pass basic + chr7 only
<pre style="font-family:Consolas">
SEQ_IDS="gcx_cmet_truseq gcx_cmet_smarter"

parallel -j 2 -k "/storage/apps/bin/STAR \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --twopassMode Basic \
    --genomeLoad NoSharedMemory \
    --runThreadN 40 \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --readQualityScoreBase 20 \
    --outFilterMismatchNmax 30 \
    --outFilterMismatchNoverLmax 0.1 \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped None\
    --genomeDir /data/wts/gcx-cmet-exon14/ref/rsem-index-hg19-chr7 \
    --readFilesIn /data/wts/gcx-cmet-exon14/{}/fastq/{}.trimmed_R1.fastq.gz /data/wts/gcx-cmet-exon14/{}/fastq/{}.trimmed_R2.fastq.gz \
    --outFileNamePrefix /data/wts/gcx-cmet-exon14/{}/bam/{}_chr7. \
    --sjdbGTFfile /data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.gtf " ::: ${SEQ_IDS}
</pre>

### (3) Sorting BAM Generation
<pre style="font-family:Consolas">
SEQ_IDS="gcx_cmet_truseq gcx_cmet_smarter"
parallel -j 2 -k "/storage/apps/bin/samtools sort -@ 25 \
    -o  /data/wts/gcx-cmet-exon14/{}/bam/{}.Aligned.out.Sorted.bam \
    /data/wts/gcx-cmet-exon14/{}/bam/{}.Aligned.out.bam" ::: ${SEQ_IDS}

/storage/apps/bin/samtools sort -@ 25 \
    -o  /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/bam/gcx_cmet_truseq_chr7.Aligned.out.Sorted.bam \
    /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/bam/gcx_cmet_truseq_chr7.Aligned.out.bam
/storage/apps/bin/samtools sort -@ 25 \
    -o  /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/bam/gcx_cmet_smarter_chr7.Aligned.out.Sorted.bam \
    /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/bam/gcx_cmet_smarter_chr7.Aligned.out.bam

</pre>
gcx_cmet_truseq_chr7.Aligned.out.bam
### (4) Picard Collect RNAseq Metric
<pre style="font-family:Consolas">
# Genome
java -Xmx16G -jar /storage/apps/bin/gatk-local.jar CollectRnaSeqMetrics \
    --INPUT /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/bam/gcx_cmet_truseq.Aligned.out.Sorted.bam \
    --OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_truseq.RNA_Metrics \
    --REF_FLAT /data/wts/gcx-cmet-exon14/ref/gcx.hg19.refFlat.txt \
    --STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
    --CHART_OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_truseq.picard.rnaseq.metric.pdf \
    --MINIMUM_LENGTH 100 \
    --RIBOSOMAL_INTERVALS /data/wts/gcx-cmet-exon14/ref/gcx.hghe19.rRNA.insterval_list \
    --REFERENCE_SEQUENCE /storage/references_and_index/gcx.hg19/gcx.hg19.fasta

java -Xmx16G -jar /storage/apps/bin/gatk-local.jar CollectRnaSeqMetrics \
    --INPUT /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/bam/gcx_cmet_smarter.Aligned.out.Sorted.bam \
    --OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_smarter.RNA_Metrics \
    --REF_FLAT /data/wts/gcx-cmet-exon14/ref/gcx.hg19.refFlat.txt \
    --STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
    --CHART_OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_smarter.picard.rnaseq.metric.pdf \
    --MINIMUM_LENGTH 100 \
    --RIBOSOMAL_INTERVALS /data/wts/gcx-cmet-exon14/ref/gcx.hg19.rRNA.insterval_list \
    --REFERENCE_SEQUENCE /storage/references_and_index/gcx.hg19/gcx.hg19.fasta

# chr7
java -Xmx16G -jar /storage/apps/bin/gatk-local.jar CollectRnaSeqMetrics \
    --INPUT /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/bam/gcx_cmet_truseq_chr7.Aligned.out.Sorted.bam \
    --OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_truseq_chr7.RNA_Metrics \
    --REF_FLAT /data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.refFlat.txt \
    --STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
    --CHART_OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_truseq_chr7.picard.rnaseq.metric.pdf \
    --MINIMUM_LENGTH 100 \
    --REFERENCE_SEQUENCE /data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.fasta

java -Xmx16G -jar /storage/apps/bin/gatk-local.jar CollectRnaSeqMetrics \
    --INPUT /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/bam/gcx_cmet_smarter_chr7.Aligned.out.Sorted.bam \
    --OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_smarter_chr7.RNA_Metrics \
    --REF_FLAT /data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.refFlat.txt \
    --STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
    --CHART_OUTPUT /data/wts/gcx-cmet-exon14/meta/gcx_cmet_smarter_chr7.picard.rnaseq.metric.pdf \
    --MINIMUM_LENGTH 100 \
    --REFERENCE_SEQUENCE /data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.fasta
</pre>

### RNA-SeQC
<pre style="font-family:Consolas">
# Genome

/storage/apps/rnaseqc-2.4.2/rnaseqc \
    --sample=gcx_cmet_truseq.rna-seqc \
    --fasta=/storage/references_and_index/gcx.hg19/gcx.hg19.fasta \
    --base-mismatch=30 \
    --gene-length=100 \
    --coverage \
    --detection-threshold=3 \
    /storage/references_and_index/gcx.hg19/gcx.hg19.gtf \
    /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/bam/gcx_cmet_truseq.Aligned.out.Sorted.bam \
    /data/wts/gcx-cmet-exon14/meta

/storage/apps/rnaseqc-2.4.2/rnaseqc \
    --sample=gcx_cmet_smarter.rna-seqc \
    --fasta=/storage/references_and_index/gcx.hg19/gcx.hg19.fasta \
    --base-mismatch=30 \
    --gene-length=100 \
    --coverage \
    --detection-threshold=3 \
    /storage/references_and_index/gcx.hg19/gcx.hg19.gtf \
    /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/bam/gcx_cmet_smarter.Aligned.out.Sorted.bam \
    /data/wts/gcx-cmet-exon14/meta &


SEQ_IDS="gcx_cmet_truseq gcx_cmet_smarter"
parallel -j 2 -k "/storage/apps/rnaseqc-2.4.2/rnaseqc \
    --sample={} \
    --fasta=/storage/references_and_index/gcx.hg19/gcx.hg19.fasta \
    --base-mismatch=30 \
    --gene-length=100 \
    --coverage \
    --detection-threshold=3 \
    /storage/references_and_index/gcx.hg19/gcx.hg19.gtf \
    /data/wts/gcx-cmet-exon14/{}/bam/{}.Aligned.out.Sorted.bam \
    /data/wts/gcx-cmet-exon14/meta/RNA-SeQC" ::: ${SEQ_IDS}

# chr7
SEQ_IDS="gcx_cmet_truseq_chr7 gcx_cmet_smarter_chr7"
parallel -j 2 -k "/storage/apps/rnaseqc-2.4.2/rnaseqc \
    --sample={} \
    --fasta=/data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.fasta \
    --base-mismatch=30 \
    --gene-length=100 \
    --coverage \
    --detection-threshold=3 \
    /data/wts/gcx-cmet-exon14/ref/gcx.hg19.chr7.gtf \
    /data/wts/gcx-cmet-exon14/{}/bam/{}.Aligned.out.bam \
    /data/wts/gcx-cmet-exon14/meta/RNA-SeQC" ::: ${SEQ_IDS}
</pre>

## RSEM
/storage/apps/rsem-1.3.3/bin/rsem-calculate-expression \
    --paired-end --no-bam-output --estimate-rspd \
    --strandedness reverse \
    --num-threads 30 \
    --alignments /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/bam/gcx_cmet_truseq.Aligned.toTranscriptome.out.bam \
    /storage/references_and_index/gcx.hg19/rsem-index-star-pe150/gcx.hg19 \
    /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/quant/gcx_cmet_truseq
      
/storage/apps/rsem-1.3.3/bin/rsem-calculate-expression \
    --paired-end --no-bam-output --estimate-rspd \
    --strandedness reverse \
    --num-threads 30 \
    --alignments /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/bam/gcx_cmet_smarter.Aligned.toTranscriptome.out.bam \
    /storage/references_and_index/gcx.hg19/rsem-index-star-pe150/gcx.hg19 \
    /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/quant/gcx_cmet_smarter
      
cat /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/quant/gcx_cmet_smarter.genes.results | grep ENSG00000105976
cat /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/quant/gcx_cmet_truseq.genes.results | grep ENSG00000105976

MET = ENSG00000105976

BRAF = ENSG00000157764
cat /data/wts/gcx-cmet-exon14/gcx_cmet_smarter/quant/gcx_cmet_smarter.genes.results | grep ENSG00000157764
cat /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/quant/gcx_cmet_truseq.genes.results | grep ENSG00000157764

- MultiQC Summarize
<pre style="font-family:Consolas">
mkdir -p /data/wts/gcx-cmet-exon14/meta/qc-data

cp /data/wts/gcx-cmet-exon14/gcx_cmet_truseq/qcfiles/*fastqc.zip /data/wts/gcx-cmet-exon14/meta/qc-data/
cp /data/wts/gcx-cmet-exon14/gcx_cmet_smarter_1020/qcfiles/*fastqc.zip /data/wts/gcx-cmet-exon14/meta/qc-data/
cp /data/wts/gcx-cmet-exon14/gcx_cmet_smarter_1013/qcfiles/*fastqc.zip /data/wts/gcx-cmet-exon14/meta/qc-data/

source /storage/apps/miniconda3/bin/activate wes
multiqc --force \
    --filename gcx-cMet-RNAseq-QC \
    --outdir /data/wts/gcx-cmet-exon14/meta \
    --data-dir /data/wts/gcx-cmet-exon14/meta/qc-data
</pre>

samtools sort -@25 -o NCI-H596_D1KNF.hg38_Aligned.out.Sorted.bam NCI-H596_D1KNF.Aligned.out.bam &
samtools sort -@25 -o NCI-H596_D1KYP.hg38_Aligned.out.Sorted.bam NCI-H596_D1KYP.Aligned.out.bam &

samtools index -@25 NCI-H596_D1KNF.hg38_Aligned.out.Sorted.bam &
samtools index -@25 NCI-H596_D1KYP.hg38_Aligned.out.Sorted.bam &
samtools index -@25 NCI-H596_C1EHH.hg38_Aligned.out.Sorted.bam & 
