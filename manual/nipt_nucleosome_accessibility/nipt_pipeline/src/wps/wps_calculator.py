"""
wps_calculator.py
=================
Window Protection Score (WPS) 계산 모듈.
Snyder et al. 2016 정의 기반.

L-WPS : 120~180 bp 단편, k=120 bp 윈도우  → 뉴클레오솜 위치 / TSS 발현량
S-WPS :  35~ 80 bp 단편, k= 16 bp 윈도우  → TF 결합 부위 (CTCF 등)

처리 파이프라인
--------------
1.  calc_wps_region()          BAM → bp-단위 raw WPS 배열
2.  normalize_wps_local()      1000 bp 블록 median 차감 (local baseline 제거)
3.  smooth_wps_savgol()        Savitzky-Golay 필터 (peak 형상 보존 denoising)
4.  call_wps_peaks()           연속 상승 구간 → 뉴클레오솜 dyad 좌표
5.  calc_wps_marker_batch()    marker BED 배치 처리 (위 1~4 통합)
6.  aggregate_wps_to_marker()  marker 당 WPS 배열 → 요약 지표
7.  run_gene_tss_wps()         유전자 TSS L-WPS 래퍼
8.  run_enhancer_wps()         Enhancer S-WPS 래퍼
9.  normalize_wps_scores()     샘플 간 Z-score / quantile 정규화

CLI 사용
--------
python wps_calculator.py \\
    --bam sample.bam \\
    --bed markers.bed \\
    --mode L \\
    --marker-type gene \\
    --out-prefix ./results/sample
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import pysam
from scipy.signal import savgol_filter
from scipy.ndimage import label as nd_label
from typing import Literal, Optional
from utils import log

# ──────────────────────────────────────────────
# 파라미터 정의
# ──────────────────────────────────────────────
WPS_PARAMS = {
    "L": {
        "frag_min":    120,
        "frag_max":    180,
        "window_size": 120,
    },
    "S": {
        "frag_min":    35,
        "frag_max":    80,
        "window_size": 16,
    },
}

WPSMode = Literal["L", "S"]


# ──────────────────────────────────────────────
# 핵심: 단일 region WPS 계산
# ──────────────────────────────────────────────
def calc_wps_region(
    bam_handle: pysam.AlignmentFile,
    chrom: str,
    region_start: int,      # 0-based, inclusive
    region_end:   int,      # 0-based, exclusive
    mode: WPSMode = "L",
    min_mapq: int = 20,
) -> np.ndarray:
    """
    region 내 모든 bp 위치 i에 대해 WPS(i, k) 를 계산합니다.

    WPS(i, k) = N_spanning - N_endpoints

    N_spanning  : 윈도우 [i - k//2, i + k//2] 를 완전히 가로지르는 단편 수
    N_endpoints : 윈도우 내부에 5'/3' 말단이 하나라도 떨어지는 단편 수

    Returns
    -------
    wps : np.ndarray (int32), shape=(region_end - region_start,)
          index 0 == region_start
    """
    params  = WPS_PARAMS[mode]
    frag_min    = params["frag_min"]
    frag_max    = params["frag_max"]
    k           = params["window_size"]
    half_k      = k // 2

    region_len = region_end - region_start
    spanning   = np.zeros(region_len, dtype=np.int32)   # N_spanning 누적
    endpoints  = np.zeros(region_len, dtype=np.int32)   # N_endpoints 누적

    # fetch padding: 윈도우가 region 경계를 넘을 수 있으므로 half_k 만큼 확장
    fetch_start = max(0, region_start - half_k)
    fetch_end   = region_end + half_k

    for read in bam_handle.fetch(chrom, fetch_start, fetch_end):
        if (read.is_unmapped or read.is_duplicate or
                read.is_secondary or read.is_supplementary):
            continue
        if read.mapping_quality < min_mapq:
            continue

        # fragment 좌표 (0-based)
        if read.is_paired and read.template_length != 0:
            frag_start = read.reference_start
            frag_end   = frag_start + abs(read.template_length)  # exclusive
        else:
            frag_start = read.reference_start
            frag_end   = read.reference_end

        frag_len = frag_end - frag_start
        if not (frag_min <= frag_len <= frag_max):
            continue

        # ── N_spanning 누적 ────────────────────────────────────────────
        # 윈도우 [i-half_k, i+half_k+1) 을 완전히 가로지르려면:
        #   frag_start < i - half_k  AND  frag_end > i + half_k
        # → i 범위: (frag_start + half_k, frag_end - half_k)  (정수)
        # region 좌표로 변환 후 clamp
        span_lo = frag_start + half_k + 1   # i > frag_start + half_k  이므로 +1
        span_hi = frag_end   - half_k       # i < frag_end - half_k    이므로 exclusive

        lo = max(span_lo, region_start) - region_start
        hi = min(span_hi, region_end)   - region_start
        if lo < hi:
            spanning[lo:hi] += 1

        # ── N_endpoints 누적 ──────────────────────────────────────────
        # 5' 말단(frag_start)이 윈도우 [i-half_k, i+half_k+1) 안에 있는 i 범위:
        #   frag_start >= i - half_k  →  i <= frag_start + half_k
        #   frag_start <  i + half_k  →  i >  frag_start - half_k
        # → i in (frag_start - half_k, frag_start + half_k + 1)
        def _add_endpoint(pos):
            ep_lo = pos - half_k
            ep_hi = pos + half_k + 1
            lo_ = max(ep_lo, region_start) - region_start
            hi_ = min(ep_hi, region_end)   - region_start
            if lo_ < hi_:
                endpoints[lo_:hi_] += 1

        _add_endpoint(frag_start)   # 5' 말단
        if frag_end != frag_start:  # 중복 방지 (길이 0 방어)
            _add_endpoint(frag_end - 1)  # 3' 말단 (inclusive 위치)

    return spanning - endpoints


# ──────────────────────────────────────────────
# Step 2. 로컬 정규화 (1000 bp 블록 median 차감)
# ──────────────────────────────────────────────
def normalize_wps_local(
    wps: np.ndarray,
    block_size: int = 1000,
) -> np.ndarray:
    """
    1000 bp 비중첩 블록 단위 median을 차감해 로컬 커버리지 편향을 제거합니다.

    각 블록 내 모든 bp 위치에서 해당 블록의 median WPS를 빼므로
    GC bias / 국소 depth 변동에 의한 baseline drift 를 억제합니다.

    Parameters
    ----------
    wps        : raw WPS 배열 (int32 or float)
    block_size : 블록 크기 (bp), 기본값 1000

    Returns
    -------
    norm_wps : float32 배열, 원본과 동일 길이
    """
    arr  = wps.astype(np.float32)
    n    = len(arr)
    out  = arr.copy()

    for start in range(0, n, block_size):
        end = min(start + block_size, n)
        out[start:end] -= np.median(arr[start:end])

    return out


# ──────────────────────────────────────────────
# Step 3. Savitzky-Golay 필터 (peak 형상 보존 smoothing)
# ──────────────────────────────────────────────
def smooth_wps_savgol(
    wps: np.ndarray,
    window_length: int = 21,
    polyorder: int = 2,
) -> np.ndarray:
    """
    Savitzky-Golay 필터로 고주파 톱니 노이즈를 제거합니다.

    단순 이동 평균과 달리 국소 다항식 회귀(least squares) 기반이므로
    뉴클레오솜 피크의 height / width 를 수학적으로 보존합니다.

    Parameters
    ----------
    wps           : normalize_wps_local() 출력
    window_length : 슬라이딩 윈도우 크기 (홀수여야 함), 기본값 21
                    L-WPS 뉴클레오솜 스케일: 21~101 권장
                    S-WPS TF footprint 스케일: 11~21 권장
    polyorder     : 다항식 차수 (window_length > polyorder), 기본값 2

    Returns
    -------
    smoothed : float32 배열
    """
    # window_length는 홀수 + 신호 길이보다 작아야 함
    wl = min(window_length, len(wps))
    if wl % 2 == 0:
        wl -= 1
    if wl <= polyorder:
        return wps.astype(np.float32)

    return savgol_filter(wps, window_length=wl, polyorder=polyorder).astype(np.float32)


# ──────────────────────────────────────────────
# Step 4. Peak Calling (뉴클레오솜 dyad 좌표 추출)
# ──────────────────────────────────────────────
def call_wps_peaks(
    wps_smoothed: np.ndarray,
    region_start: int,
    chrom: str,
    min_peak_score: float = 0.0,
    min_peak_width: int   = 50,
    max_peak_width: int   = 250,
) -> pd.DataFrame:
    """
    S-G 필터 적용 후 WPS 배열에서 뉴클레오솜 dyad 위치를 추출합니다.

    알고리즘
    --------
    1. WPS > min_peak_score 인 연속 구간을 후보 peak region 으로 탐지
       (scipy.ndimage.label 이용)
    2. 구간 길이가 [min_peak_width, max_peak_width] 범위인 것만 통과
    3. 각 통과 구간에서 WPS 최대값 위치 → dyad (뉴클레오솜 중심) 좌표

    Parameters
    ----------
    wps_smoothed   : smooth_wps_savgol() 출력
    region_start   : fetch 시작 좌표 (0-based), 절대 좌표 변환에 사용
    chrom          : 염색체명
    min_peak_score : 연속 상승 기준 threshold (기본 0.0 = 양수 구간)
    min_peak_width : 유효 peak 최소 폭 (bp)
    max_peak_width : 유효 peak 최대 폭 (bp), 뉴클레오솜 ~147 bp 고려

    Returns
    -------
    DataFrame 컬럼:
      chrom, dyad_pos (0-based), peak_start, peak_end,
      peak_width, peak_score (최대 WPS), peak_mean_score
    """
    above       = (wps_smoothed > min_peak_score).astype(int)
    labeled, n_features = nd_label(above)

    records = []
    for label_id in range(1, n_features + 1):
        idx = np.where(labeled == label_id)[0]
        width = len(idx)
        if not (min_peak_width <= width <= max_peak_width):
            continue

        peak_scores  = wps_smoothed[idx]
        local_max    = int(idx[np.argmax(peak_scores)])
        dyad_abs     = region_start + local_max

        records.append({
            "chrom":           chrom,
            "dyad_pos":        dyad_abs,
            "peak_start":      region_start + int(idx[0]),
            "peak_end":        region_start + int(idx[-1]) + 1,
            "peak_width":      width,
            "peak_score":      float(np.max(peak_scores)),
            "peak_mean_score": float(np.mean(peak_scores)),
        })

    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=["chrom", "dyad_pos", "peak_start", "peak_end",
                 "peak_width", "peak_score", "peak_mean_score"]
    )


# ──────────────────────────────────────────────
# marker 단위 WPS 집계 (정규화 + S-G + peak 통합)
# ──────────────────────────────────────────────
def aggregate_wps_to_marker(
    wps_array: np.ndarray,
    region_start: int,
    chrom: str,
    flanking: int = 2000,
    sg_window: int = 21,
    sg_polyorder: int = 2,
    block_size: int = 1000,
    min_peak_score: float = 0.0,
) -> dict:
    """
    단일 marker 의 raw WPS 배열을 정규화 → S-G smoothing → peak calling 후
    요약 지표를 반환합니다.

    Returns
    -------
    dict:
      wps_mean         : 정규화 후 전체 평균
      wps_center_mean  : 중심 ±500 bp 평균 (TSS 활성 proxy)
      wps_auc          : 정규화 후 AUC
      wps_peak         : S-G smoothed 최대값
      wps_trough       : S-G smoothed 최소값 (NDR 깊이)
      wps_oscillation  : peak - trough
      n_nucleosomes    : 탐지된 뉴클레오솜(dyad) 수
      smoothed         : S-G 적용 배열 (np.ndarray)
      peaks_df         : peak 좌표 DataFrame
    """
    # Step 2: 로컬 정규화
    norm = normalize_wps_local(wps_array, block_size=block_size)

    # Step 3: S-G smoothing
    smoothed = smooth_wps_savgol(norm, window_length=sg_window, polyorder=sg_polyorder)

    # Step 4: Peak calling
    peaks_df = call_wps_peaks(
        smoothed, region_start=region_start, chrom=chrom,
        min_peak_score=min_peak_score,
    )

    n   = len(smoothed)
    center     = n // 2
    half_inner = min(500, flanking)
    c_lo = max(0, center - half_inner)
    c_hi = min(n, center + half_inner)

    return {
        "wps_mean":        float(np.mean(norm)),
        "wps_center_mean": float(np.mean(norm[c_lo:c_hi])),
        "wps_auc":         float(np.trapz(norm)),
        "wps_peak":        float(np.max(smoothed)),
        "wps_trough":      float(np.min(smoothed)),
        "wps_oscillation": float(np.max(smoothed) - np.min(smoothed)),
        "n_nucleosomes":   len(peaks_df),
        "smoothed":        smoothed,
        "peaks_df":        peaks_df,
    }


# ──────────────────────────────────────────────
# marker BED 배치 처리
# ──────────────────────────────────────────────
def calc_wps_marker_batch(
    bam_path: str,
    marker_df: pd.DataFrame,
    mode: WPSMode = "L",
    flanking: int = 2000,
    smoothing_window: int = 100,
    min_mapq: int = 20,
    strand_aware: bool = True,
) -> pd.DataFrame:
    """
    marker_df 각 행에 대해 WPS 를 계산하고 요약 지표를 반환합니다.

    marker_df 필수 컬럼
    -------------------
    chrom       : str  (e.g. "chr1")
    center      : int  (유전자 TSS 또는 enhancer 중심, 0-based)
    strand      : str  ("+" | "-")  strand_aware=True 일 때 필요
    marker_id   : str  (유전자명, enhancer ID 등)

    Returns
    -------
    DataFrame : marker_id 당 요약 지표 + wps_profile (list 컬럼)
    """
    required = {"chrom", "center", "marker_id"}
    missing  = required - set(marker_df.columns)
    if missing:
        raise ValueError(f"marker_df에 필수 컬럼 누락: {missing}")

    records    = []
    all_peaks  = []

    with pysam.AlignmentFile(bam_path, "rb") as bam_handle:
        for row in marker_df.itertuples(index=False):
            chrom     = row.chrom
            center    = int(row.center)
            marker_id = row.marker_id
            strand    = getattr(row, "strand", "+")

            reg_start = max(0, center - flanking)
            reg_end   = center + flanking + 1

            try:
                wps = calc_wps_region(
                    bam_handle, chrom, reg_start, reg_end,
                    mode=mode, min_mapq=min_mapq,
                )
            except Exception as e:
                log(f"[WPS] {marker_id} ({chrom}:{reg_start}-{reg_end}) 실패: {e}")
                continue

            # minus strand → 반전 (TSS 기준 프로파일 정렬)
            if strand_aware and strand == "-":
                wps = wps[::-1]

            summary = aggregate_wps_to_marker(
                wps,
                region_start=reg_start,
                chrom=chrom,
                flanking=flanking,
                sg_window=smoothing_window if smoothing_window % 2 == 1 else smoothing_window + 1,
            )

            records.append({
                "marker_id":       marker_id,
                "chrom":           chrom,
                "center":          center,
                "strand":          strand,
                "mode":            mode,
                "wps_mean":        summary["wps_mean"],
                "wps_center_mean": summary["wps_center_mean"],
                "wps_auc":         summary["wps_auc"],
                "wps_peak":        summary["wps_peak"],
                "wps_trough":      summary["wps_trough"],
                "wps_oscillation": summary["wps_oscillation"],
                "n_nucleosomes":   summary["n_nucleosomes"],
                "wps_profile":     summary["smoothed"].tolist(),
            })

            # peak 좌표에 marker 정보 추가
            if not summary["peaks_df"].empty:
                pk = summary["peaks_df"].copy()
                pk["marker_id"] = marker_id
                all_peaks.append(pk)

    result_df = pd.DataFrame(records)
    peaks_df  = pd.concat(all_peaks, ignore_index=True) if all_peaks else pd.DataFrame()

    log(f"[WPS] {mode}-WPS 완료: {len(result_df)} markers, "
        f"{len(peaks_df)} peaks detected")
    return result_df, peaks_df


# ──────────────────────────────────────────────
# 유전자 / Enhancer marker 적용 전략 래퍼
# ──────────────────────────────────────────────
def run_gene_tss_wps(
    bam_path: str,
    gene_bed_df: pd.DataFrame,
    flanking: int = 5000,
    smoothing_window: int = 100,
    min_mapq: int = 20,
) -> pd.DataFrame:
    """
    유전자 TSS 기반 L-WPS 분석.

    gene_bed_df 컬럼: chrom, tss (0-based), strand, gene_id / gene_name
    → marker_df 포맷으로 변환 후 calc_wps_marker_batch 호출

    해석:
      wps_center_mean 높음 → TSS 주변 뉴클레오솜 배제(NDR) → 활성 전사
      wps_oscillation 큼   → 뉴클레오솜 위치 규칙성 강함 (positioned nucleosome)
    """
    marker_df = gene_bed_df.rename(columns={"tss": "center", "gene_id": "marker_id"}).copy()
    if "marker_id" not in marker_df.columns and "gene_name" in marker_df.columns:
        marker_df = marker_df.rename(columns={"gene_name": "marker_id"})

    return calc_wps_marker_batch(
        bam_path, marker_df,
        mode="L",
        flanking=flanking,
        smoothing_window=smoothing_window,
        min_mapq=min_mapq,
        strand_aware=True,
    )


def run_enhancer_wps(
    bam_path: str,
    enhancer_bed_df: pd.DataFrame,
    flanking: int = 1000,
    smoothing_window: int = 21,
    min_mapq: int = 20,
) -> tuple:
    """
    Enhancer 중심 S-WPS 분석 (TF footprint 탐지용).

    enhancer_bed_df 컬럼: chrom, start, end, enhancer_id, strand (optional)
    center = (start + end) // 2 자동 계산

    해석:
      wps_trough 깊음      → NDR 내 절단 빈번 → enhancer 활성화 상태
      wps_center_mean 낮음 → TF 결합 부재 또는 closed chromatin

    Returns
    -------
    (summary_df, peaks_df)
    """
    df = enhancer_bed_df.copy()
    if "center" not in df.columns:
        df["center"] = ((df["start"] + df["end"]) // 2).astype(int)
    if "marker_id" not in df.columns:
        df["marker_id"] = df["enhancer_id"]
    if "strand" not in df.columns:
        df["strand"] = "+"

    return calc_wps_marker_batch(
        bam_path, df,
        mode="S",
        flanking=flanking,
        smoothing_window=smoothing_window,
        min_mapq=min_mapq,
        strand_aware=False,
    )


# ──────────────────────────────────────────────
# 샘플 간 정규화 (Z-score 또는 quantile)
# ──────────────────────────────────────────────
def normalize_wps_scores(
    df: pd.DataFrame,
    score_col: str = "wps_center_mean",
    method: Literal["zscore", "quantile"] = "zscore",
) -> pd.DataFrame:
    """
    여러 샘플 비교를 위한 marker-level WPS 정규화.

    zscore  : (x - mean) / std  — 샘플 내 상대적 순위
    quantile: 0-1 사이 rank 변환 — 분포 형태 무관 비교
    """
    df = df.copy()
    vals = df[score_col].values.astype(float)

    if method == "zscore":
        std = vals.std()
        df[f"{score_col}_norm"] = (vals - vals.mean()) / (std if std > 0 else 1.0)
    elif method == "quantile":
        from scipy.stats import rankdata
        df[f"{score_col}_norm"] = rankdata(vals) / len(vals)

    return df


# ──────────────────────────────────────────────
# BED 로더 (공통)
# ──────────────────────────────────────────────
def _load_bed(bed_path: str, marker_type: str) -> pd.DataFrame:
    """
    BED 파일을 읽어 marker_df 포맷으로 반환합니다.

    gene BED 예상 컬럼  : chrom, tss, strand, gene_id  (또는 gene_name)
    enhancer BED 예상 컬럼: chrom, start, end, enhancer_id  (strand 선택)
    """
    df = pd.read_csv(bed_path, sep="\t", comment="#", header=0)

    if marker_type == "gene":
        # tss 컬럼이 없으면 start 를 tss 로 사용
        if "tss" not in df.columns:
            if "start" in df.columns:
                df["tss"] = df["start"]
            else:
                raise ValueError("gene BED: 'tss' 또는 'start' 컬럼이 필요합니다.")
        if "gene_id" not in df.columns and "gene_name" in df.columns:
            df["gene_id"] = df["gene_name"]
        if "strand" not in df.columns:
            df["strand"] = "+"
        df["center"]    = df["tss"].astype(int)
        df["marker_id"] = df["gene_id"].astype(str)

    elif marker_type == "enhancer":
        if "center" not in df.columns:
            df["center"] = ((df["start"].astype(int) + df["end"].astype(int)) // 2)
        if "enhancer_id" not in df.columns:
            df["enhancer_id"] = (df["chrom"].astype(str) + ":" +
                                  df["start"].astype(str) + "-" +
                                  df["end"].astype(str))
        df["marker_id"] = df["enhancer_id"].astype(str)
        if "strand" not in df.columns:
            df["strand"] = "+"
    else:
        raise ValueError(f"marker_type은 'gene' 또는 'enhancer' 여야 합니다: {marker_type}")

    return df


# ──────────────────────────────────────────────
# CLI main
# ──────────────────────────────────────────────
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="wps_calculator.py",
        description="WPS 계산 파이프라인 (Snyder et al. 2016 기반)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam",          required=True,  help="입력 BAM 파일 경로 (인덱스 필요)")
    p.add_argument("--bed",          required=True,  help="marker BED 파일 경로")
    p.add_argument("--out-prefix",   required=True,  help="출력 파일 prefix (디렉토리 포함)")
    p.add_argument("--mode",         choices=["L", "S"], default="L",
                   help="WPS 모드: L=뉴클레오솜(120-180bp), S=TF footprint(35-80bp)")
    p.add_argument("--marker-type",  choices=["gene", "enhancer"], default="gene",
                   help="marker 종류: gene(TSS L-WPS) / enhancer(S-WPS)")
    p.add_argument("--flanking",     type=int, default=None,
                   help="marker 중심 ± flanking bp (미지정 시 mode별 자동: gene=5000, enhancer=1000)")
    p.add_argument("--sg-window",    type=int, default=21,
                   help="Savitzky-Golay 윈도우 크기 (홀수)")
    p.add_argument("--sg-polyorder", type=int, default=2,
                   help="Savitzky-Golay 다항식 차수")
    p.add_argument("--block-size",   type=int, default=1000,
                   help="로컬 정규화 블록 크기 (bp)")
    p.add_argument("--min-mapq",     type=int, default=20,
                   help="최소 매핑 품질 (MAPQ)")
    p.add_argument("--norm-method",  choices=["zscore", "quantile"], default="zscore",
                   help="샘플 내 marker 점수 정규화 방법")
    p.add_argument("--no-peaks",     action="store_true",
                   help="Peak calling 결과 저장 생략")
    return p


def main():
    parser = _build_parser()
    args   = parser.parse_args()

    # ── 입력 검증 ────────────────────────────────────────────────────
    for path in [args.bam, args.bed]:
        if not os.path.exists(path):
            print(f"[ERROR] 파일을 찾을 수 없습니다: {path}", file=sys.stderr)
            sys.exit(1)

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # ── flanking 기본값 ──────────────────────────────────────────────
    if args.flanking is None:
        args.flanking = 5000 if args.marker_type == "gene" else 1000

    log(f"[WPS] 시작 | mode={args.mode} | marker={args.marker_type} "
        f"| flanking=±{args.flanking} | sg_window={args.sg_window}")

    # ── BED 로드 ─────────────────────────────────────────────────────
    marker_df = _load_bed(args.bed, args.marker_type)
    log(f"[WPS] marker 수: {len(marker_df)}")

    # ── WPS 계산 ─────────────────────────────────────────────────────
    strand_aware = (args.marker_type == "gene")
    summary_df, peaks_df = calc_wps_marker_batch(
        bam_path      = args.bam,
        marker_df     = marker_df,
        mode          = args.mode,
        flanking      = args.flanking,
        smoothing_window = args.sg_window,
        min_mapq      = args.min_mapq,
        strand_aware  = strand_aware,
    )

    # ── 샘플 내 정규화 ───────────────────────────────────────────────
    if not summary_df.empty:
        summary_df = normalize_wps_scores(
            summary_df, score_col="wps_center_mean", method=args.norm_method
        )

    # ── 출력 저장 ────────────────────────────────────────────────────
    summary_path = f"{args.out_prefix}.wps_summary.tsv"
    summary_df.drop(columns=["wps_profile"], errors="ignore").to_csv(
        summary_path, sep="\t", index=False
    )
    log(f"[WPS] summary 저장: {summary_path}")

    if not args.no_peaks and not peaks_df.empty:
        peaks_path = f"{args.out_prefix}.wps_peaks.tsv"
        peaks_df.to_csv(peaks_path, sep="\t", index=False)
        log(f"[WPS] peaks 저장: {peaks_path} ({len(peaks_df)} peaks)")

    # wps_profile (smoothed 배열) 은 npy 로 저장 (downstream 재사용 가능)
    if not summary_df.empty and "wps_profile" in summary_df.columns:
        profiles = np.array(summary_df["wps_profile"].tolist(), dtype=np.float32)
        npy_path = f"{args.out_prefix}.wps_profiles.npy"
        np.save(npy_path, profiles)
        log(f"[WPS] profiles 저장: {npy_path} shape={profiles.shape}")

    log("[WPS] 완료")


if __name__ == "__main__":
    main()