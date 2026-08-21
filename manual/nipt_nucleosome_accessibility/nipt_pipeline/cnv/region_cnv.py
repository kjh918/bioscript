"""
cnv/region_cnv.py
=================
PELT (Pruned Exact Linear Time) changepoint detection 기반
region-level CNV call

논문 근거 (Chen et al. 2022, Front. Genet.)
------------------------------------------
- bin size 100kb, Z-score threshold ±3
- Z < -3 연속 region → microdeletion
- Z > +3 연속 region → microduplication
- Sliding window 대신 PELT로 changepoint 탐지 → region 경계 정확도 향상

흐름
----
normalized.tsv (bin-level Z-score)
  → chromosome별 PELT segmentation
  → segment mean Z-score 계산
  → size / Z-score threshold 필터
  → region CNV call TSV 출력

출력 컬럼
---------
chrom, start, end, size_mb, n_bins,
mean_zscore_total, mean_zscore_short, mean_zscore_long,
call (GAIN/LOSS/NORMAL), ff_corrected_zscore
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from .utils import read_tsv, write_tsv, STANDARD_CHROMS

log = logging.getLogger(__name__)

# ── 상수 ──────────────────────────────────────────────────────────────
CALL_GAIN   = "GAIN"
CALL_LOSS   = "LOSS"
CALL_NORMAL = "NORMAL"

FIELDNAMES = [
    "chrom", "start", "end", "size_mb", "n_bins",
    "mean_zscore_total", "mean_zscore_short", "mean_zscore_long",
    "call_total", "call_short", "call_long",
    "ff_corrected_zscore_total",
]


# ── PELT 구현 ─────────────────────────────────────────────────────────
def _pelt_changepoints(
    signal:   np.ndarray,
    penalty:  float,
    min_size: int = 2,
) -> list[int]:
    """
    PELT (Pruned Exact Linear Time) changepoint detection.
    비용함수: 구간 내 잔차제곱합 (RSS) — 정규분포 가정.

    Parameters
    ----------
    signal  : 1-D float array (bin Z-score)
    penalty : BIC-based penalty = log(n) * sigma^2
              높을수록 changepoint 수 감소
    min_size: 최소 segment 크기 (bin 수)

    Returns
    -------
    changepoint 인덱스 리스트 (0-based, exclusive end)
    마지막 원소는 항상 len(signal)
    """
    n = len(signal)
    if n < min_size * 2:
        return [n]

    # 누적합 / 누적 제곱합 (RSS 빠른 계산용)
    cs  = np.concatenate([[0.0], np.cumsum(signal)])
    cs2 = np.concatenate([[0.0], np.cumsum(signal ** 2)])

    def rss(s: int, e: int) -> float:
        """구간 [s, e) 의 RSS"""
        length = e - s
        if length <= 0:
            return 0.0
        sum_x  = cs[e]  - cs[s]
        sum_x2 = cs2[e] - cs2[s]
        return sum_x2 - (sum_x ** 2) / length

    # DP
    F   = np.full(n + 1, np.inf)
    F[0] = -penalty
    last = np.zeros(n + 1, dtype=int)  # 이전 changepoint 추적

    # admissible set
    admissible = [0]

    for t in range(min_size, n + 1):
        candidates = [
            (F[s] + rss(s, t) + penalty, s)
            for s in admissible
            if t - s >= min_size
        ]
        if not candidates:
            F[t]    = F[t - 1]
            last[t] = last[t - 1]
            continue

        best_cost, best_s = min(candidates, key=lambda x: x[0])
        F[t]    = best_cost
        last[t] = best_s

        # PELT 가지치기: F[s] + penalty >= F[t] 이면 s 제거
        admissible = [
            s for s in admissible
            if F[s] + penalty <= F[t]
        ]
        admissible.append(t)

    # changepoint 역추적
    cps: list[int] = []
    t = n
    while t > 0:
        cps.append(t)
        t = last[t]
    cps.reverse()
    return cps


# ── segment → region 변환 ─────────────────────────────────────────────
@dataclass
class CnvRegion:
    chrom:                  str
    start:                  int
    end:                    int
    n_bins:                 int
    mean_zscore_total:      float
    mean_zscore_short:      float
    mean_zscore_long:       float
    call_total:             str
    call_short:             str
    call_long:              str
    ff_corrected_zscore:    float  # FF 보정 Z-score (total 기준)

    @property
    def size_mb(self) -> float:
        return round((self.end - self.start) / 1e6, 4)

    def to_dict(self) -> dict:
        return {
            "chrom":                    self.chrom,
            "start":                    self.start,
            "end":                      self.end,
            "size_mb":                  self.size_mb,
            "n_bins":                   self.n_bins,
            "mean_zscore_total":        f"{self.mean_zscore_total:.4f}",
            "mean_zscore_short":        f"{self.mean_zscore_short:.4f}",
            "mean_zscore_long":         f"{self.mean_zscore_long:.4f}",
            "call_total":               self.call_total,
            "call_short":               self.call_short,
            "call_long":                self.call_long,
            "ff_corrected_zscore_total": f"{self.ff_corrected_zscore:.4f}",
        }


def _call(z: float, gain_thr: float, loss_thr: float) -> str:
    if z >= gain_thr:
        return CALL_GAIN
    if z <= loss_thr:
        return CALL_LOSS
    return CALL_NORMAL


def _penalty(signal: np.ndarray, multiplier: float = 1.0) -> float:
    """
    BIC-based penalty: log(n) × σ²
    multiplier로 sensitivity 조정 (높을수록 changepoint 감소)
    """
    n     = len(signal)
    sigma = float(np.std(signal, ddof=1)) if n > 1 else 1.0
    return multiplier * np.log(n) * (sigma ** 2)


def _segment_chrom(
    chrom:       str,
    bins:        list[dict],
    gain_thr:    float,
    loss_thr:    float,
    min_size_mb: float,
    bin_size:    int,
    penalty_mul: float,
    ff:          Optional[float],
    zscore_col:  str,
) -> list[CnvRegion]:
    """
    단일 염색체 bin 리스트 → PELT segmentation → CnvRegion 리스트
    """
    if not bins:
        return []

    zvals = np.array([
        float(r[zscore_col]) if r.get(zscore_col, "NA") != "NA" else 0.0
        for r in bins
    ], dtype=np.float64)

    pen    = _penalty(zvals, penalty_mul)
    cps    = _pelt_changepoints(zvals, pen, min_size=1)

    regions: list[CnvRegion] = []
    prev    = 0

    for cp in cps:
        seg_bins = bins[prev:cp]
        if not seg_bins:
            prev = cp
            continue

        seg_z_total = np.array([
            float(r["zscore_total"]) if r.get("zscore_total", "NA") != "NA" else 0.0
            for r in seg_bins
        ])
        seg_z_short = np.array([
            float(r["zscore_short"]) if r.get("zscore_short", "NA") != "NA" else 0.0
            for r in seg_bins
        ])
        seg_z_long  = np.array([
            float(r["zscore_long"])  if r.get("zscore_long",  "NA") != "NA" else 0.0
            for r in seg_bins
        ])

        mean_total = float(np.mean(seg_z_total))
        mean_short = float(np.mean(seg_z_short))
        mean_long  = float(np.mean(seg_z_long))

        # call
        call_total = _call(mean_total, gain_thr, loss_thr)
        call_short = _call(mean_short, gain_thr, loss_thr)
        call_long  = _call(mean_long,  gain_thr, loss_thr)

        # region 좌표
        start = int(seg_bins[0]["bin_start"])
        end   = int(seg_bins[-1]["bin_end"])
        size_bp = end - start

        # FF 보정 Z-score: z_corrected = z_raw / (FF / 2)
        # FF=10% → 태아 기여 비율 = 0.1, diploid 기준 signal = 0.5 × FF
        ff_corrected = mean_total
        if ff and ff > 0:
            ff_corrected = mean_total / (ff / 200.0)  # FF는 % 단위

        # 크기 필터: NORMAL region은 min_size_mb 무관하게 skip 안 함
        # GAIN/LOSS는 min_size_mb 이상만 보고
        size_mb = size_bp / 1e6
        if call_total in (CALL_GAIN, CALL_LOSS) and size_mb < min_size_mb:
            prev = cp
            continue

        regions.append(CnvRegion(
            chrom                = chrom,
            start                = start,
            end                  = end,
            n_bins               = len(seg_bins),
            mean_zscore_total    = round(mean_total, 4),
            mean_zscore_short    = round(mean_short, 4),
            mean_zscore_long     = round(mean_long,  4),
            call_total           = call_total,
            call_short           = call_short,
            call_long            = call_long,
            ff_corrected_zscore  = round(ff_corrected, 4),
        ))
        prev = cp

    return regions


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def run_region_cnv(
    norm_tsv_path:          str,
    output_path:            str,
    zscore_gain_threshold:  float = 3.0,
    zscore_loss_threshold:  float = -3.0,
    min_size_mb:            float = 1.0,   # call 보고 최소 region 크기 (Mb)
    bin_size:               int   = 100_000,
    penalty_multiplier:     float = 2.0,   # PELT penalty 배수 (높을수록 보수적)
    ff_percent:             Optional[float] = None,  # fetal fraction (%)
    target_chroms:          Optional[list[str]] = None,
    report_normal:          bool  = False,  # NORMAL segment도 출력할지
) -> list[CnvRegion]:
    """
    normalized TSV → PELT segmentation → region-level CNV call TSV

    Parameters
    ----------
    norm_tsv_path          : run_normalize 출력 TSV
    output_path            : region CNV call TSV 저장 경로
    zscore_gain_threshold  : GAIN Z-score threshold
    zscore_loss_threshold  : LOSS Z-score threshold
    min_size_mb            : GAIN/LOSS call 최소 region 크기 (Mb)
                             microdeletion: 1Mb, chromosome aneuploidy: 10Mb
    bin_size               : bin 크기 (default 100kb)
    penalty_multiplier     : PELT penalty 배수
                             1.0 = BIC, 2.0 = 보수적, 0.5 = 민감
    ff_percent             : fetal fraction (%), FF 보정 Z-score 계산용
    target_chroms          : 분석할 염색체 (None이면 standard 전체)
    report_normal          : True면 NORMAL segment도 출력

    Returns
    -------
    CnvRegion 리스트
    """
    rows = read_tsv(norm_tsv_path)
    if not rows:
        raise ValueError(f"빈 TSV: {norm_tsv_path}")

    # chrom별로 bins 그룹화 (순서 유지)
    target = target_chroms or STANDARD_CHROMS
    chrom_bins: dict[str, list[dict]] = {c: [] for c in target}

    for r in rows:
        chrom = r.get("chrom", "")
        if chrom in chrom_bins:
            chrom_bins[chrom].append(r)

    all_regions: list[CnvRegion] = []

    for chrom in target:
        bins = chrom_bins[chrom]
        if not bins:
            log.debug("%s: bin 없음 — skip", chrom)
            continue

        regions = _segment_chrom(
            chrom       = chrom,
            bins        = bins,
            gain_thr    = zscore_gain_threshold,
            loss_thr    = zscore_loss_threshold,
            min_size_mb = min_size_mb,
            bin_size    = bin_size,
            penalty_mul = penalty_multiplier,
            ff          = ff_percent,
            zscore_col  = "zscore_total",
        )

        n_abnormal = sum(1 for r in regions if r.call_total != CALL_NORMAL)
        log.info(
            "  %s: %d segments, %d GAIN/LOSS (≥%.1fMb)",
            chrom, len(regions), n_abnormal, min_size_mb,
        )
        all_regions.extend(regions)

    # 출력: NORMAL 제외 여부
    out_regions = all_regions if report_normal else [
        r for r in all_regions if r.call_total != CALL_NORMAL
    ]

    rows_out = [r.to_dict() for r in out_regions]
    write_tsv(rows_out, output_path, FIELDNAMES)
    log.info(
        "region CNV 완료: %d regions (GAIN/LOSS) → %s",
        len(out_regions), output_path,
    )
    return all_regions