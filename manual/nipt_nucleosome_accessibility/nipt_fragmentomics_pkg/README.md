# nipt_fragmentomics

NIPT cfDNA fragmentomics 분석 파이프라인.  
short/long fragment 분리 기반으로 CNV, Fetal Fraction, BAF, marker WPS 를 통합 계산합니다.

---

## 목차

1. [개념 및 설계 원칙](#개념-및-설계-원칙)
2. [파이프라인 흐름](#파이프라인-흐름)
3. [디렉토리 구조](#디렉토리-구조)
4. [입력 파일](#입력-파일)
5. [출력 파일](#출력-파일)
6. [실행 방법](#실행-방법)
7. [파일별 기능 설명](#파일별-기능-설명)
8. [파라미터 레퍼런스](#파라미터-레퍼런스)
9. [알려진 제한사항](#알려진-제한사항)

---

## 개념 및 설계 원칙

### Fragment 분류

| 구분 | 길이 | 특성 |
|---|---|---|
| short | ≤ 150 bp | fetal-enriched (태아 DNA 비율 높음) |
| long  | ≥ 151 bp | maternal-enriched (모체 DNA 비율 높음) |

모든 지표(count, WPS, BAF)를 short/long 각각 따로 계산하여  
태아 신호와 모체 신호를 구분합니다.

### Fragment 귀속 방식

- **midpoint 기반 단일 귀속**: fragment 의 중간점이 속한 bin/marker 에만 집계
- 중복 집계 없음
- paired-end: read1 만 처리 + seen qname set 으로 중복 방지

### WPS (Window Protection Score)

Snyder et al. 2016 정의:
```
WPS(i, k) = N_spanning(i) - N_endpoints(i)
  N_spanning : [i-k/2, i+k/2] 를 완전히 가로지르는 fragment 수
  N_endpoints: 윈도우 내에 5'/3' 말단이 있는 fragment 수
```

| 모드 | fragment 길이 | window k | 해석 |
|---|---|---|---|
| L-WPS | 120–180 bp | 120 bp | 뉴클레오솜 보호 |
| S-WPS | 35–80 bp  | 16 bp  | TF binding footprint |

bin 단위(100kb)에서는 spanning/endpoints 가 상쇄되어 median ≈ 0 → **WPS 는 marker BED 영역 단위로 계산합니다.**

### GC 보정 방식

- CPM = count / library_total × 1e6  
  (library_total: 전체 autosome quality-pass bin 합산 — subset 아님)
- LOWESS fit on GC ∈ [0.25, 0.75] 구간
- correction = clip(fit − baseline, −0.35, +0.35)

### LOO (Leave-One-Out) 정규화

CNV 계산 시 각 염색체를 정규화할 때 자기 자신을 제외한 나머지 상염색체 median 을 baseline 으로 사용합니다.  
→ 거대 이수성(chr21 trisomy 등)이 전체 baseline 을 왜곡하는 것을 방지.

---

## 파이프라인 흐름

```
BAM
 │
 ├─ Step 1: bin_extractor
 │    100kb grid (--bin-bed) 기준으로 BAM 스캔
 │    → bins_raw.parquet
 │      short_count, long_count, breadth, gc, mappability
 │
 ├─ Step 2: gc_corrector
 │    LOWESS GC 보정 + mappability 마스킹
 │    → bins_corrected.parquet
 │      log2_corrected_short/long, mappability_pass
 │
 ├─ Step 3: marker_extractor  [--marker-bed 있을 때]
 │    marker BED 영역별로 BAM 재스캔 (크기 무관, 100bp 도 OK)
 │    → marker_stats.parquet
 │      short/long count, ratio, WPS_L/S, breadth
 │
 ├─ Step 4: fetal_fraction
 │    SeqFF (short/long ratio) + Y-chr (남아 판정)
 │    → fetal_fraction.json
 │
 ├─ Step 5: cnv_caller
 │    LOO 정규화 → MAD z-score → CBS 세그멘테이션 → CNV call
 │    → cnv_calls.parquet
 │
 └─ Step 6: baf_calculator  [--vcf 있을 때]
      population SNP site 기반 BAF 계산 (short/long/combined)
      → bins_baf.parquet → cnv_baf.parquet (cnv_calls 와 join)
```

---

## 디렉토리 구조

```
nipt_fragmentomics/
├── run.py                          ← CLI entry point
├── setup.py
├── README.md
└── nipt_fragmentomics/
    ├── core/
    │   ├── constants.py            ← 전역 파라미터 (WPS_PARAMS, FNAME, ...)
    │   └── schema.py               ← FragmentScore dataclass
    ├── steps/
    │   ├── bin_extractor.py        ← Step 1: BAM → bins_raw.parquet
    │   ├── gc_corrector.py         ← Step 2: GC/Mappability 보정
    │   ├── marker_extractor.py     ← Step 3: marker BED → marker_stats.parquet
    │   ├── fetal_fraction.py       ← Step 4: FF 추정
    │   ├── cnv_caller.py           ← Step 5: LOO + z-score + CBS
    │   └── baf_calculator.py       ← Step 6: VCF 기반 BAF
    ├── scripts/
    │   ├── wps_compute.py          ← (독립 실행) marker WPS 계산
    │   └── wps_profile_plot.py     ← (독립 실행) WPS 프로파일 시각화
    ├── pipeline.py                 ← Step 1~6 오케스트레이터
    └── viz/
        ├── cnv_track.py            ← CNV + BAF genome track
        ├── qc_dashboard.py         ← GC bias / FF / breadth QC
        └── wps_profile.py          ← WPS 시각화 (marker_stats 기반)
```

---

## 입력 파일

| 옵션 | 설명 | 필수 |
|---|---|---|
| `--bam` | 인덱스된 BAM (.bai 필요) | ✅ |
| `--out-dir` | 결과 저장 디렉토리 | ✅ |
| `--bin-bed` | 100kb bin BED.gz (없으면 자동 생성) | 권장 |
| `--marker-bed` | marker BED (아래 포맷 참조) | Step 3 |
| `--fasta` | 참조 FASTA (.fai 필요, GC 계산용) | 권장 |
| `--bw` | Mappability bigWig (hg38.100mer.bw 등) | 권장 |
| `--vcf` | population SNP VCF.gz (.tbi 필요, BAF용) | Step 6 |

### marker BED 포맷

```
# header 없음, 탭 구분
chrom   start   end   marker_id             marker_type
chr1    104896  105048  chr1:104896-105048  cCREs
chr1    138866  139134  chr1:138866-139134  promoter
```

- 5열 고정 (marker_type 없으면 "marker" 로 자동 설정)
- 영역 크기 제한 없음 (100bp 짜리 promoter 도 OK)

---

## 출력 파일

```
{out_dir}/
├── bins_raw.parquet          Step 1 출력
├── bins_corrected.parquet    Step 2 출력
├── marker_stats.parquet      Step 3 출력 (--marker-bed 있을 때)
├── fetal_fraction.json       Step 4 출력
├── cnv_calls.parquet         Step 5 출력
├── bins_baf.parquet          Step 6 출력 (--vcf 있을 때)
├── cnv_baf.parquet           Step 6 출력 (cnv_calls + BAF join)
├── run_manifest.json         파라미터 + 소요 시간 기록
├── pipeline.log              전체 로그
└── viz/
    ├── qc_dashboard.html     GC bias / FF / breadth QC
    ├── cnv_track.html        CNV + BAF genome-wide track
    └── marker_wps.html       marker WPS 시각화
```

### marker_stats.parquet 컬럼

| 컬럼 | 타입 | 설명 |
|---|---|---|
| chrom, start, end | str/int | 좌표 |
| marker_id | str | marker 식별자 |
| marker_type | str | marker 그룹 (cCREs, promoter 등) |
| short_count | int32 | ≤150bp fragment 수 |
| long_count | int32 | ≥151bp fragment 수 |
| total_count | int32 | 전체 fragment 수 |
| short_ratio | float32 | short/total (fetal fraction proxy) |
| short_wps_L | float32 | short L-WPS median (뉴클레오솜 보호) |
| long_wps_L | float32 | long L-WPS median |
| short_wps_S | float32 | short S-WPS median (TF binding) |
| long_wps_S | float32 | long S-WPS median |
| short_breadth | float32 | short coverage 비율 |
| long_breadth | float32 | long coverage 비율 |
| gc | float32 | 해당 영역 GC 비율 |
| mappability | float32 | 해당 영역 mappability |

### cnv_baf.parquet 컬럼 (cnv_calls + BAF join)

| 컬럼 | 설명 |
|---|---|
| short/long_copy_number | LOO copy number (정상=2.0) |
| short/long_log2_norm | log2(CN/2), 정상=0.0 |
| short/long_zscore | MAD z-score |
| short/long_cnv_call | gain / loss / normal |
| segment_id | CBS 세그먼트 ID |
| sex_call | XX / XY |
| combined_baf_median | short+long 합산 BAF (정상≈0.5) |
| short_baf_median | short BAF |
| long_baf_median | long BAF |

---

## 실행 방법

### 기본 (CNV + FF)

```bash
python -m run \
    --bam  sample.bam \
    --out-dir ./results/SID001 \
    --bin-bed bins_100kb.bed.gz \
    --fasta hg38.fa \
    --bw hg38.100mer.bw \
    --jobs 8 --sample-id SID001
```

### marker WPS 포함

```bash
python -m run \
    --bam  sample.bam \
    --out-dir ./results/SID001 \
    --bin-bed bins_100kb.bed.gz \
    --marker-bed cCREs.bed \
    --fasta hg38.fa \
    --bw hg38.100mer.bw \
    --jobs 8 --sample-id SID001 --resume
```

### 전체 (CNV + FF + marker WPS + BAF)

```bash
python -m run \
    --bam  sample.bam \
    --out-dir ./results/SID001 \
    --bin-bed bins_100kb.bed.gz \
    --marker-bed cCREs.bed \
    --fasta hg38.fa \
    --bw hg38.100mer.bw \
    --vcf population_snp.vcf.gz \
    --jobs 8 --sample-id SID001 --resume
```

### 재실행 (중간 파일 있으면 건너뜀)

```bash
python -m run ... --resume
```

---

## 파일별 기능 설명

---

### `core/constants.py`

전역 파라미터 관리. 알고리즘 파라미터를 한 곳에서 관리합니다.

```python
SHORT_MAX = 150          # short/long 경계
WPS_PARAMS = {
    "L": {"frag_min": 120, "frag_max": 180, "window": 120},
    "S": {"frag_min":  35, "frag_max":  80, "window":  16},
}
FNAME = {
    "bins_raw":      "bins_raw.parquet",
    "marker_stats":  "marker_stats.parquet",
    "cnv_baf":       "cnv_baf.parquet",
    ...
}
```

---

### `core/schema.py`

`FragmentScore` dataclass. BAM read → fragment 좌표/길이/is_short 를 추출합니다.

```python
@dataclass
class FragmentScore:
    frag_start: int
    frag_end:   int
    frag_len:   int
    midpoint:   int
    is_short:   bool    # frag_len <= SHORT_MAX

    @classmethod
    def from_read(cls, read) -> Optional["FragmentScore"]: ...
```

---

### `steps/bin_extractor.py`

**BAM → bins_raw.parquet**

- 염색체별 1회 순차 스캔 (ProcessPoolExecutor, 염색체 단위 병렬)
- midpoint 기반 단일 bin 귀속 (bisect, O(log n))
- BED 없으면 BAM 헤더에서 자동 grid 생성
- BED 에 gc/mappability 컬럼 있으면 FASTA/bigWig 재계산 생략

**주요 함수**

| 함수 | 설명 |
|---|---|
| `run()` | 공개 진입점. BED → 병렬 스캔 → parquet 저장 |
| `_scan_chrom()` | 단일 염색체 스캔 worker (모듈 레벨, pickle 가능) |
| `_BinAccumulator` | bin 단위 count/breadth 누적기 |

---

### `steps/gc_corrector.py`

**bins_raw → bins_corrected**

- `mappability_pass` 마스킹: gc=0, mappability < 0.75, 성염색체 제외
- CPM = count / **전체** autosome quality-pass bin 합산 × 1e6
- LOWESS GC fit (frac=0.5, GC 범위 0.25~0.75 구간만)
- correction = clip(fit − baseline, −0.35, +0.35)
- 불량 bin → log2_corrected = NaN

**주의**: fit_mask subset 합산이 아닌 전체 library_total 사용 (subset 사용 시 ~17% 과대평가)

---

### `steps/marker_extractor.py`

**BAM + marker BED → marker_stats.parquet**

bin_extractor 와 동일한 `_MarkerAccumulator` 로직 사용. 차이점:

- 입력이 100kb grid 가 아닌 marker BED (크기 무관)
- 같은 염색체 marker 를 1회 fetch 로 일괄 처리 (속도 최적화)
- WPS (L/S 모드 × short/long) 를 marker 영역 배열로 직접 계산
- `to_profile()`: bp 단위 WPS/coverage 배열 (center 기준 상대 위치)

**주요 함수**

| 함수 | 설명 |
|---|---|
| `run()` | marker BED 로드 → 염색체별 병렬 스캔 → parquet 저장 |
| `load_marker_bed()` | BED 파싱, header 없음 5열 포맷 |
| `_MarkerAccumulator` | marker 단위 누적기 |
| `_MarkerAccumulator.to_row()` | WPS median 단일값 dict 반환 |
| `_MarkerAccumulator.to_profile()` | bp 단위 WPS 배열 반환 (center 기준) |
| `_scan_chrom_markers()` | 염색체 내 모든 marker 1회 fetch 처리 |

---

### `steps/fetal_fraction.py`

**bins_corrected → fetal_fraction.json**

**방법 1 — SeqFF**: `FF = α + β × (Σshort_CPM / Σlong_CPM)`  
(Larsen et al. 2017, autosome mappability-pass bin 만 사용)

**방법 2 — Y-chr**: chrY mean count / autosome median  
Y fraction > threshold → 남아(XY), `FF ≈ 2 × y_fraction`

**consensus_ff**: SeqFF 와 Y-chr FF 의 가중 평균

---

### `steps/cnv_caller.py`

**bins_corrected → cnv_calls.parquet**

1. **LOO 정규화** (`loo_normalize()`)
   - quality_mask = mappability_pass & gc>0 & isfinite 적용
   - 각 염색체: baseline = 나머지 상염색체 quality-pass bin median
   - copy_number = 2.0 × linear / baseline (정상=2.0)
   - 불량 bin → NaN

2. **MAD z-score** (`_robust_zscore()`)
   - 기준: autosome quality-pass bin median/MAD
   - MAD=0 이면 std fallback

3. **성염색체 z-score 보정**
   - FF 있으면 chrX/chrY 기대값 보정 (fetal XX/XY 비율 반영)

4. **CBS 세그멘테이션** (`_segment_chrom()`)
   - 재귀 Welch t-test, short log2_norm 기준

---

### `steps/baf_calculator.py`

**BAM + VCF → bins_baf.parquet**

- VCF 에서 population AF 0.2~0.8 범위 SNP site 추출
- BAM pileup으로 site 별 ref/alt depth 수집 (short/long 분리)
- bin 단위 BAF median/std/MAD/n_sites 집계
- `merge_into_cnv()`: cnv_calls.parquet 와 left join → cnv_baf.parquet

**BAF 해석**
- 정상 heterozygous: BAF ≈ 0.5
- 이수성 gain: BAF 편향 (태아 FF 낮으면 신호 약함)

---

### `scripts/wps_compute.py`

**독립 실행 가능. marker BED 기반 WPS 계산 (원저자 방식)**

```bash
python -m nipt_fragmentomics.scripts.wps_compute \
    --bam sample.bam \
    --out-prefix ./results/SID001/wps/SID001 \
    --marker-bed cCREs.bed \
    --mode L --frag long \
    --extend 1000 --jobs 8
```

- 같은 염색체 marker 를 1회 BAM fetch 로 처리
- `adjusted_WPS = (WPS − windowMedian) / Coverage × 100`
  (원저자 shendurelab/cfDNA normalizeWPSwigs.py 방식)
- 출력: `{prefix}.wps_{mode}_{frag}.npy`
  ```python
  {
    "wps_norm": {marker_id: float32 array [2*extend+1]},
    "mode": "L",
    "frag": "long",
    "extend": 1000,
  }
  ```

---

### `scripts/wps_profile_plot.py`

**독립 실행 가능. wps_compute npy → WPS 평균 프로파일 HTML**

```bash
python -m nipt_fragmentomics.scripts.wps_profile_plot \
    --npy SID001.wps_L_long.npy SID001.wps_L_short.npy \
    --labels "Long (maternal)" "Short (fetal)" \
    --marker-bed cCREs.bed \
    --extend 2000 \
    --out profile.html
```

- x축: marker center 기준 ±extend bp
- y축: 상하위 1% 제외 trimmed median (그룹 간 공유 고정)
- SEM 음영 표시
- `--group-col marker_type` 으로 marker_type 별 비교

---

### `viz/cnv_track.py`

**cnv_calls/cnv_baf.parquet → CNV genome track HTML**

- Row 1: Short log₂ norm (LOO), gain=red / loss=blue
- Row 2: Long log₂ norm
- Row 3: BAF (combined/short/long), y=0.3~0.7 고정
- Row 4: VAF (있을 때)
- **y축 대칭 고정**: |min| = |max| (양수/음수 동일 범위)

---

### `viz/qc_dashboard.py`

**bins_corrected + fetal_fraction → QC 대시보드 HTML**

- GC bias curve (LOWESS fit 오버레이)
- Mappability 분포
- short/long breadth 분포
- short_count vs long_count scatter
- Fetal Fraction 요약

---

### `viz/wps_profile.py`

**marker_stats.parquet → WPS 시각화**

**`plot_marker_wps()`**: marker 단위 WPS 요약 산점도 + box

**`plot_fragment_track()`**: 특정 region WPS/Endpoint/Coverage 3-track

**`plot_fragment_track_compare()`**: short vs long 6-track 비교

**y축 정책**
- short/long 동일 지표끼리 공유 고정
- 상하위 1% trimmed range + pad=0.12

---

## 파라미터 레퍼런스

### run.py 전체 옵션

```
참조 파일
  --bin-bed         CNV bin BED (없으면 --bin-size 로 자동 생성)
  --marker-bed      marker BED (Step 3, WPS + count 통계)
  --fasta           GC 계산 참조 FASTA
  --bw              Mappability bigWig
  --vcf             Population SNP VCF (Step 6, BAF 계산)

Bin 파라미터
  --bin-size        100000  자동 grid bin 크기 (bp)
  --min-mapq        20
  --min-baseq       20
  --min-mappability 0.75

CNV 파라미터
  --zscore-gain     3.0   gain 판정 z-score
  --zscore-loss    -3.0   loss 판정 z-score

BAF 파라미터
  --baf-af-min      0.2   population AF 하한
  --baf-af-max      0.8   population AF 상한
  --baf-min-depth   5     site 당 최소 depth

실행 옵션
  --jobs            4     병렬 프로세스 수
  --resume          중간 파일 있으면 해당 Step 건너뜀
  --no-viz          시각화 생략
  --sample-id       로그/그래프 제목에 표시될 샘플 ID
  --log-level       INFO
```

---

## 알려진 제한사항

1. **marker WPS 프로파일 그래프 (PENDING)**  
   marker_stats.parquet 는 WPS median 단일값만 저장.  
   TSS 기준 position별 WPS 변화 프로파일(x=±extend, y=WPS)은  
   `wps_compute.py` 를 별도 실행 후 `wps_profile_plot.py` 로 시각화.

2. **BAF 는 population VCF 필요**  
   내부 SNP calling 없음. KOVA/gnomAD 기반 population SNP VCF 준비 필요.

3. **SeqFF 계수 검증 필요**  
   constants.py 의 SEQFF_ALPHA/BETA 는 문헌 근사값.  
   기관 정상 코호트로 재보정 권장.

4. **CBS 세그멘테이션 속도**  
   t-test 기반 재귀 분할 → bin 수 많으면 느림.  
   bin 수 > 5000 이면 SEG_MIN_BINS 값 조정 고려.
