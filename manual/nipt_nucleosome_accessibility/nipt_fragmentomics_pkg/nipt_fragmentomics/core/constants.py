"""
전역 상수 및 파라미터 기본값.
알고리즘 파라미터는 모두 여기서 관리합니다.
"""

# ── Fragment 길이 경계 ────────────────────────────────────────────────
SHORT_MAX: int = 150   # ≤ 150 bp : fetal-enriched short fragment
LONG_MIN:  int = 151   # ≥ 151 bp : maternal-enriched long fragment

# ── WPS 파라미터 (Snyder et al. 2016) ────────────────────────────────
WPS_PARAMS: dict = {
    "L": {"frag_min": 120, "frag_max": 180, "window": 120},  # 뉴클레오솜
    "S": {"frag_min":  35, "frag_max":  80, "window":  16},  # TF footprint
}

# ── Bin 기본값 ───────────────────────────────────────────────────────
DEFAULT_BIN_SIZE:   int   = 100_000   # 100 kb
MIN_MAPQ:           int   = 20
MIN_BASEQ:          int   = 20
MIN_MAPPABILITY:    float = 0.75

# ── GC 보정 ─────────────────────────────────────────────────────────
GC_RANGE        = (0.25, 0.75)
GC_LOWESS_FRAC  = 0.5
GC_PSEUDOCOUNT  = 0.25
GC_CLIP         = 0.35

# ── Fetal Fraction ───────────────────────────────────────────────────
# SeqFF 선형 계수 (Larsen et al. 2017 근사 — 검증 후 교체)
SEQFF_ALPHA:          float = -0.186
SEQFF_BETA:           float =  1.658
Y_MALE_THRESHOLD:     float =  1.10e-4   # ChrY% = Y_reads/total_reads 기준
                                          # FMY hotspot 제거 후 최적값
                                          # Kim et al. AUC=0.996, 민감도 99.53%

# ── CNV 판정 ────────────────────────────────────────────────────────
ZSCORE_GAIN:    float =  3.0
ZSCORE_LOSS:    float = -3.0
SEG_MIN_BINS:   int   =  4
SEG_ALPHA:      float =  0.01

# ── 출력 파일명 ──────────────────────────────────────────────────────
FNAME = {
    "bins_raw":       "bins_raw.parquet",
    "bins_corrected": "bins_corrected.parquet",
    "marker_stats":   "marker_stats.parquet",
    "fetal_fraction": "fetal_fraction.json",
    "bins_baf":       "bins_baf.parquet",
    "cnv_calls":      "cnv_calls.parquet",
    "cnv_baf":        "cnv_baf.parquet",
    "marker_profiles":"marker_stats_profiles.npy",
    "manifest":       "run_manifest.json",
}