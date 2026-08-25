"""
cnv/within_sample_cnv.py
========================
Within-sample reference bin 기반 CNV detection.

원리
----
각 target bin에 대해 다른 염색체에서 GC content와 mappability이 유사한
reference bin을 선택하고 within-sample Z-score를 계산한다.

reference selection:
  1. 같은 chromosome 제외
  2. chr13, chr15, chr18, chr21, chrX, chrY reference 제외
  3. GC + mappability 기반 KDTree nearest-neighbor
  4. reference index는 1회 계산 후 total/short/long/short_fraction에서 재사용

출력
----
1. region CNV 결과
2. bin-level within-sample normalized 결과

normalized bin columns:
  chrom, bin_start, bin_end, gc, mappability,
  total_gc, short_gc, long_gc, short_fraction_gc,
  wss_z_total, wss_z_short, wss_z_long, wss_z_short_fraction
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import ttest_ind

from .utils import read_tsv, write_tsv, STANDARD_CHROMS

log = logging.getLogger(__name__)

CALL_GAIN = "GAIN"
CALL_LOSS = "LOSS"
CALL_NORMAL = "NORMAL"

_EXCLUDED_FROM_REF = {"chr13", "chr15", "chr18", "chr21", "chrX", "chrY"}

FIELDNAMES = [
    "chrom", "start", "end", "size_mb", "n_bins",
    "median_zscore_total", "median_zscore_short", "median_zscore_long",
    "median_zscore_short_fraction",
    "call_total", "call_short", "call_long", "call_short_fraction",
    "ff_corrected_zscore_total", "call_ff_corrected",
]

NORMALIZED_BIN_FIELDNAMES = [
    "chrom", "bin_start", "bin_end", "gc", "mappability",
    "total_gc", "short_gc", "long_gc", "short_fraction_gc",
    "wss_z_total", "wss_z_short", "wss_z_long", "wss_z_short_fraction",
]


@dataclass
class BinData:
    chrom: str
    start: int
    end: int
    gc: float
    mappability: float
    total_gc: float
    short_gc: float
    long_gc: float
    short_fraction_gc: float


@dataclass
class WSSSegment:
    chrom: str
    start: int
    end: int
    n_bins: int
    median_z_total: float
    median_z_short: float
    median_z_long: float
    median_z_short_fraction: float
    call_total: str
    call_short: str
    call_long: str
    call_short_fraction: str
    ff_corrected: float
    call_ff_corrected: str

    @property
    def size_mb(self) -> float:
        return round((self.end - self.start) / 1e6, 4)

    def to_dict(self) -> dict:
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "size_mb": self.size_mb,
            "n_bins": self.n_bins,
            "median_zscore_total": f"{self.median_z_total:.4f}",
            "median_zscore_short": f"{self.median_z_short:.4f}",
            "median_zscore_long": f"{self.median_z_long:.4f}",
            "median_zscore_short_fraction": f"{self.median_z_short_fraction:.4f}",
            "call_total": self.call_total,
            "call_short": self.call_short,
            "call_long": self.call_long,
            "call_short_fraction": self.call_short_fraction,
            "ff_corrected_zscore_total": f"{self.ff_corrected:.4f}",
            "call_ff_corrected": self.call_ff_corrected,
        }


def _robust_center_scale(values: np.ndarray) -> tuple[float, float]:
    """Median + MAD 기반 robust center/scale."""
    if values.size == 0:
        return 0.0, 0.0
    center = float(np.median(values))
    mad = float(np.median(np.abs(values - center)))
    return center, 1.4826 * mad


def _build_reference_indices(
    bins: list[BinData],
    k: int = 100,
    gc_weight: float = 0.7,
    map_weight: float = 0.3,
    query_multiplier: int = 5,
    min_query_extra: int = 500,
) -> list[np.ndarray]:
    """
    모든 bin에 대한 reference index를 한 번만 계산.
    GC + mappability weighted Euclidean distance 기반 KDTree 사용.
    """
    n = len(bins)
    if n == 0:
        return []

    candidate_idx = np.asarray([
        i for i, b in enumerate(bins)
        if b.chrom not in _EXCLUDED_FROM_REF
        and b.total_gc > 0
        and np.isfinite(b.gc)
        and np.isfinite(b.mappability)
    ], dtype=np.int64)

    if candidate_idx.size == 0:
        log.warning("reference candidate가 없습니다.")
        return [np.empty(0, dtype=np.int64) for _ in bins]

    sqrt_gc = np.sqrt(gc_weight)
    sqrt_map = np.sqrt(map_weight)

    coords = np.column_stack([
        np.asarray([bins[i].gc for i in candidate_idx], dtype=np.float64) * sqrt_gc,
        np.asarray([bins[i].mappability for i in candidate_idx], dtype=np.float64) * sqrt_map,
    ])

    tree = cKDTree(coords)
    query_k = min(candidate_idx.size, max(k * query_multiplier, k + min_query_extra))

    log.info(
        "reference KDTree: bins=%d candidate=%d k=%d query_k=%d",
        n, candidate_idx.size, k, query_k,
    )

    references: list[np.ndarray] = []

    for i, b in enumerate(bins):
        if not np.isfinite(b.gc) or not np.isfinite(b.mappability):
            references.append(np.empty(0, dtype=np.int64))
            continue

        point = np.asarray([b.gc * sqrt_gc, b.mappability * sqrt_map], dtype=np.float64)
        _, positions = tree.query(point, k=query_k, workers=-1)
        positions = np.atleast_1d(positions)

        selected = []
        for pos in positions:
            pos = int(pos)
            if pos < 0 or pos >= candidate_idx.size:
                continue
            ref_idx = int(candidate_idx[pos])
            if ref_idx == i:
                continue
            if bins[ref_idx].chrom == b.chrom:
                continue
            selected.append(ref_idx)
            if len(selected) >= k:
                break

        references.append(np.asarray(selected, dtype=np.int64))

        if (i + 1) % 10000 == 0 or i + 1 == n:
            log.info("reference 탐색: %d/%d bins", i + 1, n)

    counts = np.asarray([len(x) for x in references])
    if counts.size:
        log.info(
            "reference 완료: median=%d min=%d max=%d",
            int(np.median(counts)), int(np.min(counts)), int(np.max(counts)),
        )

    return references


def _compute_within_sample_z(
    bins: list[BinData],
    ref_indices: list[np.ndarray],
    metric: str,
    robust: bool = True,
    min_reference_bins: int = 10,
    exclude_zero_reference: bool = True,
) -> np.ndarray:
    """bin별 within-sample Z-score 계산."""
    cov = np.asarray([getattr(b, metric) for b in bins], dtype=np.float64)
    z_scores = np.zeros(len(bins), dtype=np.float32)

    insufficient = 0
    zero_scale = 0

    for i, refs in enumerate(ref_indices):
        if refs.size < min_reference_bins:
            insufficient += 1
            continue

        ref_cov = cov[refs]
        valid = np.isfinite(ref_cov)

        if exclude_zero_reference:
            valid &= ref_cov > 0

        ref_cov = ref_cov[valid]

        if ref_cov.size < min_reference_bins:
            insufficient += 1
            continue

        target = cov[i]

        if not np.isfinite(target):
            continue

        if robust:
            center, scale = _robust_center_scale(ref_cov)
            if scale <= 1e-12:
                center = float(np.mean(ref_cov))
                scale = float(np.std(ref_cov, ddof=1))
        else:
            center = float(np.mean(ref_cov))
            scale = float(np.std(ref_cov, ddof=1))

        if not np.isfinite(scale) or scale <= 1e-12:
            zero_scale += 1
            continue

        z = (target - center) / scale

        if np.isfinite(z):
            z_scores[i] = z

    log.info(
        "%s Z-score 완료: insufficient=%d zero_scale=%d",
        metric, insufficient, zero_scale,
    )

    return z_scores


def _save_normalized_bins(
    bins: list[BinData],
    z_total: np.ndarray,
    z_short: np.ndarray,
    z_long: np.ndarray,
    z_sfrac: np.ndarray,
    output_path: str,
) -> None:
    """bin-level within-sample normalized 값을 TSV로 저장."""
    rows = []
    for i, b in enumerate(bins):
        rows.append({
            "chrom": b.chrom,
            "bin_start": b.start,
            "bin_end": b.end,
            "gc": f"{b.gc:.6f}",
            "mappability": f"{b.mappability:.6f}",
            "total_gc": f"{b.total_gc:.6f}",
            "short_gc": f"{b.short_gc:.6f}",
            "long_gc": f"{b.long_gc:.6f}",
            "short_fraction_gc": f"{b.short_fraction_gc:.6f}",
            "wss_z_total": f"{float(z_total[i]):.6f}",
            "wss_z_short": f"{float(z_short[i]):.6f}",
            "wss_z_long": f"{float(z_long[i]):.6f}",
            "wss_z_short_fraction": f"{float(z_sfrac[i]):.6f}",
        })

    write_tsv(rows, output_path, NORMALIZED_BIN_FIELDNAMES)
    log.info("normalized bins 저장: %d bins -> %s", len(rows), output_path)


def _cbs_segment(
    z_scores: np.ndarray,
    alpha: float = 0.01,
    min_width: int = 3,
) -> list[tuple[int, int]]:
    """간소화된 CBS segmentation."""
    n = len(z_scores)

    if n == 0:
        return []
    if n < min_width * 2:
        return [(0, n - 1)]

    segments: list[tuple[int, int]] = []

    def _find_cp(left: int, right: int) -> Optional[int]:
        seg = z_scores[left:right + 1]
        if len(seg) < min_width * 2:
            return None
        mean = seg.mean()
        cumsum = np.cumsum(seg - mean)
        inner = cumsum[min_width - 1:len(seg) - min_width]
        if inner.size == 0:
            return None
        cp_local = int(np.argmax(np.abs(inner))) + min_width
        return left + cp_local

    def _recurse(left: int, right: int) -> None:
        if right - left + 1 < min_width * 2:
            segments.append((left, right))
            return
        cp = _find_cp(left, right)
        if cp is None:
            segments.append((left, right))
            return
        seg_l = z_scores[left:cp]
        seg_r = z_scores[cp:right + 1]
        if len(seg_l) < 2 or len(seg_r) < 2:
            segments.append((left, right))
            return
        _, pval = ttest_ind(seg_l, seg_r, equal_var=False, nan_policy="omit")
        if np.isfinite(pval) and pval < alpha:
            _recurse(left, cp - 1)
            _recurse(cp, right)
        else:
            segments.append((left, right))

    _recurse(0, n - 1)
    segments.sort()
    return segments


def _call(z: float, gain_thr: float, loss_thr: float) -> str:
    if z >= gain_thr:
        return CALL_GAIN
    if z <= loss_thr:
        return CALL_LOSS
    return CALL_NORMAL


def _ff_correct(z: float, ff_percent: Optional[float]) -> tuple[float, bool]:
    if ff_percent is not None and ff_percent > 0:
        expected = ff_percent / 200.0
        if expected > 0:
            return z / expected, True
    return z, False


def run_within_sample_cnv(
    norm_tsv_path: str,
    output_path: str,
    zscore_gain_threshold: float = 3.0,
    zscore_loss_threshold: float = -3.0,
    min_size_mb: float = 1.0,
    k_ref: int = 100,
    cbs_alpha: float = 0.01,
    ff_percent: Optional[float] = None,
    target_chroms: Optional[list[str]] = None,
    report_normal: bool = False,
    normalized_bin_output_path: Optional[str] = None,
    robust_zscore: bool = True,
    min_reference_bins: int = 10,
    gc_weight: float = 0.7,
    map_weight: float = 0.3,
) -> list[WSSSegment]:
    """
    Within-sample normalization + CBS segmentation + CNV call.

    Parameters
    ----------
    norm_tsv_path
        normalize.py 출력 TSV (GC correction 결과)
    output_path
        segment-level CNV 결과 TSV
    normalized_bin_output_path
        bin-level WSS Z-score 저장 경로. None이면 자동 생성.
    robust_zscore
        True이면 median/MAD 기반, False이면 mean/std 기반.
    """
    rows = read_tsv(norm_tsv_path)

    if not rows:
        raise ValueError(f"빈 TSV: {norm_tsv_path}")

    log.info("within-sample CNV 시작: %s", norm_tsv_path)
    log.info("입력 bins=%d", len(rows))

    all_bins: list[BinData] = []
    for r in rows:
        total_gc = float(r.get("total_gc", 0) or 0)
        short_gc = float(r.get("short_gc", 0) or 0)
        long_gc  = float(r.get("long_gc", 0) or 0)
        sf_raw   = r.get("short_fraction_gc", "NA")
        sf_gc    = float(sf_raw) if sf_raw not in ("NA", "", None) else (
            short_gc / total_gc if total_gc > 0 else 0.0
        )
        all_bins.append(BinData(
            chrom             = r["chrom"],
            start             = int(r["bin_start"]),
            end               = int(r["bin_end"]),
            gc                = float(r.get("gc", 0) or 0),
            mappability       = float(r.get("mappability", 1) or 1),
            total_gc          = total_gc,
            short_gc          = short_gc,
            long_gc           = long_gc,
            short_fraction_gc = sf_gc,
        ))

    log.info("reference bin 탐색 시작")
    ref_indices = _build_reference_indices(
        all_bins,
        k              = k_ref,
        gc_weight      = gc_weight,
        map_weight     = map_weight,
    )

    log.info("within-sample Z-score 계산 시작")
    z_total = _compute_within_sample_z(
        all_bins, ref_indices, "total_gc",
        robust=robust_zscore, min_reference_bins=min_reference_bins,
        exclude_zero_reference=True,
    )
    z_short = _compute_within_sample_z(
        all_bins, ref_indices, "short_gc",
        robust=robust_zscore, min_reference_bins=min_reference_bins,
        exclude_zero_reference=True,
    )
    z_long = _compute_within_sample_z(
        all_bins, ref_indices, "long_gc",
        robust=robust_zscore, min_reference_bins=min_reference_bins,
        exclude_zero_reference=True,
    )
    z_sfrac = _compute_within_sample_z(
        all_bins, ref_indices, "short_fraction_gc",
        robust=robust_zscore, min_reference_bins=min_reference_bins,
        exclude_zero_reference=False,
    )

    # bin-level 저장
    if normalized_bin_output_path is None:
        base = output_path[:-4] if output_path.endswith(".tsv") else output_path
        normalized_bin_output_path = base + ".normalized_bins.tsv"

    _save_normalized_bins(all_bins, z_total, z_short, z_long, z_sfrac,
                          normalized_bin_output_path)

    log.info("CBS segmentation 시작")
    target  = target_chroms or STANDARD_CHROMS
    bin_idx = {chrom: [] for chrom in target}

    for i, b in enumerate(all_bins):
        if b.chrom in bin_idx:
            bin_idx[b.chrom].append(i)

    all_segs: list[WSSSegment] = []

    for chrom in target:
        idxs = bin_idx[chrom]
        if not idxs:
            continue

        z_chr = z_total[idxs]

        try:
            raw_segs = _cbs_segment(z_chr, alpha=cbs_alpha)
        except Exception as e:
            log.warning("%s CBS 실패: %s", chrom, e)
            raw_segs = [(0, len(idxs) - 1)]

        chrom_abnormal = 0

        for seg_start, seg_end in raw_segs:
            seg_idxs = idxs[seg_start:seg_end + 1]
            if not seg_idxs:
                continue

            med_total = float(np.median(z_total[seg_idxs]))
            med_short = float(np.median(z_short[seg_idxs]))
            med_long  = float(np.median(z_long [seg_idxs]))
            med_sfrac = float(np.median(z_sfrac[seg_idxs]))

            ff_z, ff_used = _ff_correct(med_total, ff_percent)

            call_total  = _call(med_total, zscore_gain_threshold, zscore_loss_threshold)
            call_short  = _call(med_short, zscore_gain_threshold, zscore_loss_threshold)
            call_long   = _call(med_long,  zscore_gain_threshold, zscore_loss_threshold)
            call_sfrac  = _call(med_sfrac, zscore_gain_threshold, zscore_loss_threshold)
            call_ff     = _call(ff_z,      zscore_gain_threshold, zscore_loss_threshold)

            primary_call = call_ff if ff_used else call_total

            s_start = all_bins[seg_idxs[0]].start
            s_end   = all_bins[seg_idxs[-1]].end
            size_mb = (s_end - s_start) / 1e6

            if primary_call != CALL_NORMAL and size_mb < min_size_mb:
                continue

            seg = WSSSegment(
                chrom                   = chrom,
                start                   = s_start,
                end                     = s_end,
                n_bins                  = len(seg_idxs),
                median_z_total          = round(med_total, 4),
                median_z_short          = round(med_short, 4),
                median_z_long           = round(med_long,  4),
                median_z_short_fraction = round(med_sfrac, 4),
                call_total              = call_total,
                call_short              = call_short,
                call_long               = call_long,
                call_short_fraction     = call_sfrac,
                ff_corrected            = round(ff_z, 4),
                call_ff_corrected       = call_ff,
            )
            all_segs.append(seg)

            if primary_call != CALL_NORMAL:
                chrom_abnormal += 1
                log.info(
                    "  %s %s:%d-%d %.2fMb n=%d z=%.3f ff_z=%.3f",
                    primary_call, chrom, s_start, s_end,
                    size_mb, len(seg_idxs), med_total, ff_z,
                )

        log.info("%s: %d abnormal segments", chrom, chrom_abnormal)

    # 출력 필터
    if report_normal:
        out_segs = all_segs
    elif ff_percent is not None and ff_percent > 0:
        out_segs = [s for s in all_segs if s.call_ff_corrected != CALL_NORMAL]
    else:
        out_segs = [s for s in all_segs if s.call_total != CALL_NORMAL]

    write_tsv([s.to_dict() for s in out_segs], output_path, FIELDNAMES)
    log.info(
        "within-sample CNV 완료: segments=%d output=%d -> %s",
        len(all_segs), len(out_segs), output_path,
    )
    log.info("bin normalized -> %s", normalized_bin_output_path)

    return all_segs