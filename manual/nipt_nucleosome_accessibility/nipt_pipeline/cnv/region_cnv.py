"""
cnv/region_cnv.py
=================
ruptures PELT (l2 모델) 기반 region-level CNV segmentation + call

판정 기준
---------
FF가 있을 때: ff_corrected_zscore 기준
  ff_corrected_z = median_z / (FF% / 200)
  이유: FF=10% 샘플에서 deletion signal = FF/2 = 5%
        raw Z-score 기준으로는 threshold 미달이어도
        FF로 나누면 실제 copy number 변화에 비례한 Z가 됨

FF가 없을 때 (female fetus): raw Z-score 기준

출력 컬럼
---------
chrom, start, end, size_mb, n_bins,
median_zscore_total, median_zscore_short, median_zscore_long,
median_zscore_short_fraction,
call_total, call_short, call_long, call_short_fraction,
ff_corrected_zscore_total, call_ff_corrected
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from .utils import read_tsv, write_tsv, STANDARD_CHROMS

log = logging.getLogger(__name__)

CALL_GAIN   = "GAIN"
CALL_LOSS   = "LOSS"
CALL_NORMAL = "NORMAL"

FIELDNAMES = [
    "chrom", "start", "end", "size_mb", "n_bins",
    "median_zscore_total", "median_zscore_short", "median_zscore_long",
    "median_zscore_short_fraction",
    "call_total", "call_short", "call_long", "call_short_fraction",
    "ff_corrected_zscore_total", "call_ff_corrected",
]


# ── PELT ─────────────────────────────────────────────────────────────
def _pelt_ruptures(
    signal:   np.ndarray,
    penalty:  float = 3.0,
    min_size: int   = 3,
    jump:     int   = 3,
) -> list[int]:
    """ruptures PELT l2 모델 — O(n), 빠름"""
    try:
        import ruptures as rpt
        model = rpt.Pelt(model="l2", min_size=min_size, jump=jump)
        model.fit(signal.reshape(-1, 1))
        return model.predict(pen=penalty)
    except ImportError:
        log.warning("ruptures 없음 → RSS fallback")
        return _pelt_rss_fallback(signal, penalty, min_size)
    except Exception as e:
        log.warning("ruptures 오류 (%s) → RSS fallback", e)
        return _pelt_rss_fallback(signal, penalty, min_size)


def _pelt_rss_fallback(signal: np.ndarray, penalty: float, min_size: int) -> list[int]:
    n   = len(signal)
    cs  = np.concatenate([[0.0], np.cumsum(signal)])
    cs2 = np.concatenate([[0.0], np.cumsum(signal ** 2)])

    def rss(s: int, e: int) -> float:
        length = e - s
        if length <= 0:
            return 0.0
        return cs2[e] - cs2[s] - (cs[e] - cs[s]) ** 2 / length

    F    = np.full(n + 1, np.inf)
    F[0] = -penalty
    last = np.zeros(n + 1, dtype=int)
    adm  = [0]

    for t in range(min_size, n + 1):
        cands = [(F[s] + rss(s, t) + penalty, s)
                 for s in adm if t - s >= min_size]
        if not cands:
            F[t] = F[t - 1]; last[t] = last[t - 1]
            continue
        best_cost, best_s = min(cands)
        F[t] = best_cost; last[t] = best_s
        adm  = [s for s in adm if F[s] + penalty <= F[t]]
        adm.append(t)

    cps: list[int] = []
    t = n
    while t > 0:
        cps.append(t); t = last[t]
    cps.reverse()
    return cps


# ── Segment ───────────────────────────────────────────────────────────
@dataclass
class Segment:
    chrom:                   str
    start:                   int
    end:                     int
    n_bins:                  int
    start_idx:               int
    end_idx:                 int
    median_z_total:          float
    median_z_short:          float
    median_z_long:           float
    median_z_short_fraction: float
    call_total:              str
    call_short:              str
    call_long:               str
    call_short_fraction:     str
    ff_corrected:            float
    call_ff_corrected:       str

    @property
    def size_mb(self) -> float:
        return round((self.end - self.start) / 1e6, 4)

    def to_dict(self) -> dict:
        return {
            "chrom":                        self.chrom,
            "start":                        self.start,
            "end":                          self.end,
            "size_mb":                      self.size_mb,
            "n_bins":                       self.n_bins,
            "median_zscore_total":          f"{self.median_z_total:.4f}",
            "median_zscore_short":          f"{self.median_z_short:.4f}",
            "median_zscore_long":           f"{self.median_z_long:.4f}",
            "median_zscore_short_fraction": f"{self.median_z_short_fraction:.4f}",
            "call_total":                   self.call_total,
            "call_short":                   self.call_short,
            "call_long":                    self.call_long,
            "call_short_fraction":          self.call_short_fraction,
            "ff_corrected_zscore_total":    f"{self.ff_corrected:.4f}",
            "call_ff_corrected":            self.call_ff_corrected,
        }


def _call(z: float, gain_thr: float, loss_thr: float) -> str:
    if z >= gain_thr: return CALL_GAIN
    if z <= loss_thr: return CALL_LOSS
    return CALL_NORMAL


def _ff_corrected_zscore(median_z: float, ff_percent: Optional[float]) -> tuple[float, bool]:
    """
    FF 보정 Z-score 계산.

    원리:
      태아 heterozygous deletion/duplication의 경우
      expected signal = FF / 2 (copy number 변화 × fetal fraction)
      ff_corrected_z = median_z / (FF% / 200)

    Returns
    -------
    (ff_corrected_z, ff_used)
      ff_used: True면 FF 보정 적용, False면 raw Z 그대로
    """
    if ff_percent and ff_percent > 0:
        expected_signal = ff_percent / 200.0  # FF%/2/100
        if expected_signal > 0:
            return median_z / expected_signal, True
    return median_z, False


def _make_segments(
    chrom:       str,
    bins:        list[dict],
    bkps:        list[int],
    gain_thr:    float,
    loss_thr:    float,
    min_size_mb: float,
    ff:          Optional[float],
) -> list[Segment]:
    segs: list[Segment] = []
    prev = 0

    for bkp in bkps:
        seg_bins = bins[prev:bkp]
        if not seg_bins:
            prev = bkp; continue

        def _z(key: str) -> np.ndarray:
            return np.array([
                float(r[key]) if r.get(key, "NA") != "NA" else 0.0
                for r in seg_bins
            ])

        z_total          = _z("zscore_total")
        z_short          = _z("zscore_short")
        z_long           = _z("zscore_long")
        z_short_fraction = _z("zscore_short_fraction")

        med_total          = float(np.median(z_total))
        med_short          = float(np.median(z_short))
        med_long           = float(np.median(z_long))
        med_short_fraction = float(np.median(z_short_fraction))

        # raw Z-score 기반 call (reference)
        call_total          = _call(med_total,          gain_thr, loss_thr)
        call_short          = _call(med_short,          gain_thr, loss_thr)
        call_long           = _call(med_long,           gain_thr, loss_thr)
        call_short_fraction = _call(med_short_fraction, gain_thr, loss_thr)

        # FF 보정 Z-score 기반 call (primary)
        ff_z, ff_used = _ff_corrected_zscore(med_total, ff)
        call_ff = _call(ff_z, gain_thr, loss_thr)

        start   = int(seg_bins[0]["bin_start"])
        end     = int(seg_bins[-1]["bin_end"])
        size_mb = (end - start) / 1e6

        # primary call: FF 보정이 적용됐으면 ff_corrected call,
        #               아니면 raw total call 사용
        primary_call = call_ff if ff_used else call_total

        # 크기 필터: primary call이 GAIN/LOSS이고 min_size_mb 미달이면 skip
        if primary_call != CALL_NORMAL and size_mb < min_size_mb:
            prev = bkp; continue

        segs.append(Segment(
            chrom                   = chrom,
            start                   = start,
            end                     = end,
            n_bins                  = len(seg_bins),
            start_idx               = prev,
            end_idx                 = bkp,
            median_z_total          = round(med_total,          4),
            median_z_short          = round(med_short,          4),
            median_z_long           = round(med_long,           4),
            median_z_short_fraction = round(med_short_fraction, 4),
            call_total              = call_total,
            call_short              = call_short,
            call_long               = call_long,
            call_short_fraction     = call_short_fraction,
            ff_corrected            = round(ff_z, 4),
            call_ff_corrected       = call_ff,
        ))
        prev = bkp

    return segs


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def run_region_cnv(
    norm_tsv_path:          str,
    output_path:            str,
    zscore_gain_threshold:  float         = 3.0,
    zscore_loss_threshold:  float         = -3.0,
    min_size_mb:            float         = 1.0,
    penalty:                float         = 3.0,
    ff_percent:             Optional[float] = None,
    target_chroms:          Optional[list[str]] = None,
    report_normal:          bool          = False,
) -> list[Segment]:
    """
    normalized TSV → PELT segmentation → region CNV call TSV

    Parameters
    ----------
    norm_tsv_path         : run_normalize 출력 TSV
    output_path           : region CNV TSV 저장 경로
    zscore_gain_threshold : GAIN Z-score threshold (FF 보정 후 적용)
    zscore_loss_threshold : LOSS Z-score threshold (FF 보정 후 적용)
    min_size_mb           : GAIN/LOSS report 최소 region 크기 (Mb)
    penalty               : PELT l2 penalty (1.0=민감, 3.0=기본, 8.0=보수적)
    ff_percent            : fetal fraction (%)
                            있으면: ff_corrected_z = raw_z / (FF/200) 기준 판정
                            없으면: raw Z-score 기준 판정
    target_chroms         : 분석할 염색체 (None이면 standard 전체)
    report_normal         : True면 NORMAL segment도 출력

    Returns
    -------
    전체 Segment 리스트 (NORMAL 포함, plot 용도)
    """
    rows = read_tsv(norm_tsv_path)
    if not rows:
        raise ValueError(f"빈 TSV: {norm_tsv_path}")

    if ff_percent:
        log.info("FF 보정 Z-score 적용: FF=%.2f%% → expected_signal=%.4f",
                 ff_percent, ff_percent / 200.0)
    else:
        log.info("FF 없음 → raw Z-score 기준 판정")

    target = target_chroms or STANDARD_CHROMS
    chrom_bins: dict[str, list[dict]] = {c: [] for c in target}
    for r in rows:
        c = r.get("chrom", "")
        if c in chrom_bins:
            chrom_bins[c].append(r)

    all_segs: list[Segment] = []

    for chrom in target:
        bins = chrom_bins[chrom]
        if not bins:
            continue

        z_total = np.array([
            float(r["zscore_total"]) if r.get("zscore_total", "NA") != "NA" else 0.0
            for r in bins
        ])
        z_frac = np.array([
            float(r.get("zscore_short_fraction", "NA"))
            if r.get("zscore_short_fraction", "NA") != "NA" else 0.0
            for r in bins
        ])

        # total Z-score 기반 segmentation
        bkps_total = _pelt_ruptures(z_total, penalty=penalty, min_size=3)
        segs_total = _make_segments(
            chrom, bins, bkps_total,
            zscore_gain_threshold, zscore_loss_threshold,
            min_size_mb, ff_percent,
        )

        # short_fraction Z-score 기반 segmentation (보완적)
        bkps_frac = _pelt_ruptures(z_frac, penalty=penalty, min_size=3)
        segs_frac = _make_segments(
            chrom, bins, bkps_frac,
            zscore_gain_threshold, zscore_loss_threshold,
            min_size_mb, ff_percent,
        )

        # 병합: frac에서만 잡힌 non-NORMAL region 추가
        def _overlaps(s: Segment, existing: list[Segment]) -> bool:
            return any(s.start < e.end and s.end > e.start for e in existing)

        segs = list(segs_total)
        for s in segs_frac:
            if s.call_ff_corrected != CALL_NORMAL and not _overlaps(s, segs_total):
                segs.append(s)
                log.info("    [frac-only] %s %s:%d-%d (%.2fMb) frac_z=%.3f",
                         s.call_ff_corrected, chrom, s.start, s.end,
                         s.size_mb, s.median_z_short_fraction)

        # 로그
        n_ab = sum(1 for s in segs if s.call_ff_corrected != CALL_NORMAL)
        log.info("  %s: %d segments, %d GAIN/LOSS (≥%.1fMb, FF-corrected)",
                 chrom, len(segs), n_ab, min_size_mb)
        for s in segs:
            if s.call_ff_corrected != CALL_NORMAL:
                log.info(
                    "    %s %s:%d-%d (%.2fMb) "
                    "raw_z=%.3f ff_corr_z=%.3f",
                    s.call_ff_corrected, chrom, s.start, s.end, s.size_mb,
                    s.median_z_total, s.ff_corrected,
                )

        all_segs.extend(segs)

    # 출력: report_normal=False이면 FF 보정 기준 GAIN/LOSS만
    out_segs = all_segs if report_normal else [
        s for s in all_segs if s.call_ff_corrected != CALL_NORMAL
    ]
    write_tsv([s.to_dict() for s in out_segs], output_path, FIELDNAMES)
    log.info("region CNV 완료: %d GAIN/LOSS regions → %s", len(out_segs), output_path)
    return all_segs