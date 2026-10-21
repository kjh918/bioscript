
# 1. DATA 개요

| 테라젠 폴더이름 |  생성날짜 | 라이브러리 | Type |Seq.Type | Note |
|:---:|:---:|:---:|:---:|:---:|:---|
|TBD231186_18832_20231013|  2023-10-13 | SMARTer<br>picomamalian kit | Total RNA | PE-150 | stranded seq |
|TBD231186_18832_20231020|  2023-10-19 | SMARTer<br>picomamalian kit | Total RNA | PE-150 | stranded seq |
|TBD231186_18685_20231020|  2023-09-22 | Illumina<br>TruSeq RNA Exome kit | mRNA Exon | PE-150 | stranded seq |


- `TBD231186_18832_20231020`과 `TBD231186_18832_20231020`은 같은 샘플임

> - Sample ID assign
>   - TBD231186_18832_20231013 = <b>gcx_cmet_smarter</b>
>   - TBD231186_18685_20231020 = <b>gcx_cmet_truseq</b>

---

# 2. fastq file QC

### A. gcx_cmet_truseq
<p align="center">
<img src="../img/gcx_cmet_truseq_r1_fastqc.png" width="400px" height="250px" title="gcx_cmet_truseq_R1"/>   
<img src="../img/gcx_cmet_truseq_r2_fastqc.png" width="400px" height="250px" title="gcx_cmet_truseq_R2"/>   
</p>
---

### B. gcx_cmet_smarter
<p align="center">
<img src="../img/gcx_cmet_smarter_R1.png" width="400px" height="250px" title="gcx_cmet_smarter_R1"/>  
<img src="../img/gcx_cmet_smarter_R2.png" width="400px" height="250px" title="gcx_cmet_smarter_R2"/>   
</p>
---

<br>
<p align="center">

| Sample | Read Length | GC (%) | Total Sequences | Total Sequences(M)
|:---:|:---:|:---:|---:|---:|
| gcx_cmet_truseq_R1 | 151 | 48 | 32,226,503 | 32.2 M |
| gcx_cmet_truseq_R2 | 151 | 49 | 32,226,503 | 32.3 M |
| gcx_cmet_smarter_R1 | 151 | 60 | 36,834,596 | 36.8 M |
| gcx_cmet_smarter_R2 | 151 | 61 | 36,834,596 | 36.8 M |
</p>
<br>

# 3. fastq file trimming & filtering
- Filtering tools & parameters
  
| Ints. | swtools | quality cut-off | minimum length cut-off |
|---:|:---:|:---:|:---:|
| 테라젠   | cutadapt | 20 | 50 |
| 젠큐릭스 | fastp | 20 | 50 |

---
<br>
<p align="center">
<img src="../img/gcx_cmet_fastp_filter.png" width="600px" height="350px" title="gcx_cmet_smarter_R1"/>  
</p>

---

<br>

# 4. Alignment
- 젠큐릭스 : `STAR` --twopass Basic 모드

| Paramter | 테라젠 | 젠큐릭스 | 기본값 |
| ---:|:---:|:---:|:---:|
| genome ver. | GRCh38.p13 | GRCh37 (hg19) | - |
| STAR ver. | 2.7.1a | 2.7.11a | - |
| STAR Index sjdbOverhang | - | 150 | 99 |
| outFilterType | BySJout | 기본값 | BySJout |
| outFilterMultimapNmax | 20 | 기본값 | 20 |
| alignSJoverhangMin | 8 | 기본값 | 5 |
| alignSJDBoverhangMin  | 1 | 기본값 | 1 | 
| `outFilterMismatchNmax` | `999` | `30` | `10` |
| `outFilterMismatchNoverLmax`  | `0.04` | `0.1` | `0.3` |
| alignIntronMin | 20 | 기본값 | 20 |
| alignIntronMax | 1000000 | 기본값 | 1000000 |
| alignMatesGapMax | 1000000 | 기본값 | 1000000 |
---
<br> 

# 5. Alignment QC 
- `RNA-SeQC` (v2.4.2)
- picard `CollectRnaSeqMetrics` (v3.1.0)
  
### 5.1 alignmnet qc(1)

| SAMPLE ID        | Mapped Reads | Uniq Reads | Sliced<br>(Alternative.Alignments) |  InsertSize | 
| ---:|---:|---:|---:|---:|
| 테라젠 SMARTer    |  62,436,900 |24,389,070 | 4,197,674 | `20.98` |  
| gcx_cmet_smarter |  57,940,802 | 57,940,802 | 35,999,548 | -  |
| 테라젠 TruSeq     |  57,921,210 | 45,805,890 | 45,805,890 | 169.01 |
| gcx_cmet_truseq  |  58,475,931 | 58,475,931 | 20,919,143 | -  |

---
---
<br>

- `Too Low Proper-Sized Inserts`

<p align="center">
<img src="../img/gcx_cmet_fastp_insertsize.png" width="600px" height="350px" title="gcx_cmet_smarter_R1"/>  
</p>
---

### 5.2 Alignment qc(2)

| Parameter |  gcx_cmet_smarter | gcx_cmet_truseq |
| ---:|---:|---:|
| `Intergenic.Rate` (%) | `57.1` | 3.0|
| Intronic.Rate (%) | 22.7 | 2.7 |
| `Intragenic.Rate` (%) | `40.0` | `84.0` |
| `Exonic.Rate` (%) | 16.9 | `81.4` | 
| `High.Quality.Exonic.Rate` (%) | 29.3 | `85.5` |
| rRNA.rate (%) | 2.1 | 2.3 |
| Intergenic.Reads | 33,063,385 | 1,761,396 |
| Intronic.Reads   | 13,143,591 | 1,575,122 |
| `Intragenic.Reads` | 22,914,165 | `49,163,880` |
| `Median.of.Avg.Transcript.Coverage` | 1.67 | `0.57` |
| PCT_USABLE_BASES (%) | 12.4 | 76.1 |
| PCT_CODING_BASES (%) | 8.4 | 66.8 |
| PCT_INTRONIC_BASES (%) | 23.9 | 7.0 |
| PCT_INTERGENIC_BASES (%) | 63.1 | 16.0 |
| PCT_MRNA_BASES (%) | 12.9 | 77.0 | 
---
<br>
<p align="center">
<img src="../img/gcx_cmet_tx_type_ratio.png" width="700px" height="350px" title="gcx_cmet_smarter_R1"/>  
</p>
<br>

---
<br>
<p align="center">
<img src="../img/gcx_met_base_distribution.png" width="700px" height="350px" align="center" title="gcx_cmet_smarter_R1"/>  
</p>
---

## 6. Quantitation

| Library| MET counts | MET TPM | total reads on genes |
| ---:|:---:|:---:|---:|
| gcx_cmet_smarter | 140 | 17.28 |  3,036,608 |
| gcx_cmet_truseq  | 920 | 53.99 | 20,411,729 |


---

## 6. IGV view

### NSCLC Cell Lines 
- NCI-H596 : MET exon14 skipping (+)
- A549 : MET exon14 skipping (-)
<p align="center">
<img src="../img/cell_line_IGV.png" width="1000px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

| chromosome | start | end | NCI-H596 | A549 | exon no. |
|:---:|---:|---:|---:|---:|---:|
| chr7 | 116763049 | 116763049 | 253 | 392 | #11 |
| chr7 | 116763049 | 116763268 | 544 | 921 | #11 |
| chr7 | 116769665 | 116769791 | 312 | 586 | #12 |
| chr7 | 116771497 | 116771654 | 274 | 608 | `#13` |
| chr7 | 116771848 | 116771989 | `106` | 585 | `#14` |
| chr7 | 116774880 | 116775072 | 333 | 675 | `#15` |
| chr7 | 116774880 | 116775111 | 388 | 752 | `#15` |
| chr7 | 116777388 | 116777469 | 324 | 462 | #16 |
| chr7 | 116778775 | 116778957 | 458 | 691 | #17 |

---

## GCX Processing : MET

<p align="center">
<img src="../img/gcx_MET_IGV.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

| chromosome | start | end | TruSeq | SMARTer | exon no. |
|:---:|---:|---:|---:|---:|---:|
| 7 | 116403103 | 116403134 | 0 | 10 | #11 | 
| 7 | 116403103 | 116403322 | 0 | 26 | #11 |
| 7 | 116409698 | 116409845 | 214 | 24 | #12 |
| 7 | 116409698 | 116409963 | 214 | 24 | #12 |
| 7 | 116409719 | 116409845 | 214 | 24 | #12 |
| 7 | 116411551 | 116411708 | `0` | `0`  | `#13` |
| 7 | 116411708 | 116412043 | `0` | `13` | `#14` |
| 7 | 116414934 | 116415126 | `0` | `8`  | `#15` |
| 7 | 116414934 | 116415165 | `0` | `8`  | `#15` |
| 7 | 116417442 | 116417523 | 146 | 47 | #16 |
| 7 | 116418829 | 116419011 | 124 | 22 | #17 |

---
<br>

## GCX Processing : Other Genes

### BRAF
<p align="center">
<img src="../img/gcx_BRAF_igv_01.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

### EGFR
<p align="center">
<img src="../img/gcx_EGFR_igv_01.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

### ERBB2
<p align="center">
<img src="../img/gcx_ERBB2_igv_01.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

### SMAD4 
<p align="center">
<img src="../img/gcx_SMAD4_igv_01.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

### KRAS
<p align="center">
<img src="../img/gcx_KRAS_igv_01.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

### TP53
<p align="center">
<img src="../img/gcx_TP53_igv_01.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>

---

### 테라젠 MET exon14
<p align="center">
<img src="../img/theragen_met_exon14.png" width="800px" height="450px" title="gcx_cmet_smarter_R1"/>  
</p>



---

## Conclusion

#### 1. 이 두 샘플의 결과로는 MET exon14 skipping 에대한 `판단을 할 수 없다`.
    --->  SMARTer 라이브러리 데이터는 활용이 어려운 상태이며,
          TruSeq 라이브러리는 유전자별 exon 커버리지 차이가 심해서 skipping에 대한 판단을 내릴 수 없다.

#### 2. Library 제작 과정에서 intact 한 라이브러리가 제작 되지 않았을 가능성

#### 3. FFPE --> RNA 과정의 문제일 가능성

 
         
---


## Question

 1. SMARTer Library 샘플의 Alignment QC에서 Insert Size = ~20 bp 로 레포트에 나와있는데 이 샘플의 매핑 결과가 의미가 있다고 판단한 근거는 무엇인지?
  
 2. RNA-SeQC 또는 Picard 같이 잘 알려진 툴로 별도의 RNAseq alignment QC 체크를 하신 것이 있는지?

 3. STAR mapping parameter에서 'outFilterMismatchNmax' 값을 999로 설정한 이유가 있는지?

 4. IGV 상의 exon14 부분 매핑에서 대부분이 미스매치로 보이는데 이것이 올바른 매핑 결과로 간주할수 있는것인지?

 5. Exon 지역에대한 비율, 실제 리드수가 모두 현저히 낮은데 SMARTer 라이브러리 결과를 채택하신 근거가 있는지?

 6. (1) 질문과 관련하여 라이브러리 제작후 수행하셨을 tape station peak 결과를 확인할수 있을까요?
