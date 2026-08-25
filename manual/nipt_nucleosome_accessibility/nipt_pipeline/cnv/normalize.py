"""
cnv/normalize.py
================
bin count TSV -> GC correction TSV

변경 이력
---------
이전:
  GC correction + genome-wide Z-score

현재:
  GC correction만 수행
  -> Z-score는 within_sample_cnv.py (WSS 방식)에서 계산

개선사항
--------
1. raw count 대신 log(count + pseudocount) 공간에서 robust LOWESS 수행
   - high-GC 영역의 heteroscedastic noise 영향 감소
   - chr19 같은 GC-rich chromosome 과보정 완화

2. 기본 LOESS span을 0.45로 증가
   - 0.3보다 local fluctuation을 덜 따라가도록 보수적으로 smoothing

3. short_fraction은 count와 다른 방식으로 correction
   - multiplicative correction 대신 additive residual correction
   - 0~1 범위 clip

4. chromosome별 GC correction diagnostic logging 추가
   - raw median
   - corrected median
   - correction factor
   - median GC

출력 컬럼
---------
chrom, bin_start, bin_end, gc, mappability,
total_raw, short_raw, long_raw,
total_gc, short_gc, long_gc,
short_fraction_gc
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

from .utils import read_tsv, write_tsv


logger = logging.getLogger(__name__)


FIELDNAMES = [
    "chrom",
    "bin_start",
    "bin_end",
    "gc",
    "mappability",
    "total_raw",
    "short_raw",
    "long_raw",
    "total_gc",
    "short_gc",
    "long_gc",
    "short_fraction_gc",
]


# =====================================================================
# 공통 helper
# =====================================================================

def _build_gc_fit_mask(
    values: np.ndarray,
    gc: np.ndarray,
    gc_trim_pct: float,
    require_positive: bool = True,
) -> tuple[np.ndarray, float, float]:
    """
    LOWESS fitting에 사용할 GC 범위와 mask 생성.

    Parameters
    ----------
    values
        fitting 대상 값.
    gc
        GC fraction.
    gc_trim_pct
        GC 분포 양 끝에서 제거할 비율.
        예: 0.01 -> lower/upper 1% 제거.
    require_positive
        True이면 values > 0인 값만 사용.
        count에 사용.
        fraction에는 False 사용.

    Returns
    -------
    fit_mask, gc_lo, gc_hi
    """

    values = np.asarray(values, dtype=np.float64)
    gc = np.asarray(gc, dtype=np.float64)

    valid = np.isfinite(values) & np.isfinite(gc)

    if require_positive:
        valid &= values > 0

    if valid.sum() < 10:
        return valid, np.nan, np.nan

    gc_valid = gc[valid]

    gc_lo = float(
        np.percentile(
            gc_valid,
            gc_trim_pct * 100.0,
        )
    )

    gc_hi = float(
        np.percentile(
            gc_valid,
            (1.0 - gc_trim_pct) * 100.0,
        )
    )

    fit_mask = (
        valid
        & (gc >= gc_lo)
        & (gc <= gc_hi)
    )

    return fit_mask, gc_lo, gc_hi


# =====================================================================
# Count용 LOESS GC correction
# =====================================================================

def _loess_gc_correct(
    counts: np.ndarray,
    gc: np.ndarray,
    frac: float = 0.45,
    max_correction_ratio: float = 3.0,
    gc_trim_pct: float = 0.01,
    pseudocount: float = 0.5,
) -> np.ndarray:
    """
    Robust LOESS GC bias correction.

    log(count + pseudocount) ~ GC 를 LOWESS fitting한 뒤
    GC-dependent trend만 제거한다.

    corrected_log
        = observed_log
        - predicted_log
        + median(predicted_log)

    이후 count scale로 복원한다.

    안전장치
    --------
    1. GC 극단값 상/하 gc_trim_pct 제거
    2. robust LOWESS it=3
    3. log-space fitting으로 high-count/high-GC 분산 영향 감소
    4. np.interp endpoint 사용으로 극단 GC extrapolation 방지
    5. correction ratio가 max_correction_ratio 초과하면 raw fallback
    """

    counts = np.asarray(
        counts,
        dtype=np.float64,
    )

    gc = np.asarray(
        gc,
        dtype=np.float64,
    )

    valid = (
        np.isfinite(counts)
        & np.isfinite(gc)
        & (counts > 0)
    )

    if valid.sum() < 10:
        logger.warning(
            "유효 bin이 너무 적어 GC correction 스킵"
        )
        return counts.copy()

    # ---------------------------------------------------------
    # fitting GC 범위
    # ---------------------------------------------------------
    fit_mask, gc_lo, gc_hi = _build_gc_fit_mask(
        counts,
        gc,
        gc_trim_pct=gc_trim_pct,
        require_positive=True,
    )

    if fit_mask.sum() < 10:
        logger.warning(
            "GC trim 후 적합 bin 부족 - GC correction 스킵"
        )
        return counts.copy()

    logger.debug(
        "GC LOWESS fit range: %.4f - %.4f (%d bins)",
        gc_lo,
        gc_hi,
        int(fit_mask.sum()),
    )

    # ---------------------------------------------------------
    # log count
    # ---------------------------------------------------------
    log_counts = np.log(
        counts + pseudocount
    )

    x = gc[fit_mask]
    y = log_counts[fit_mask]

    order = np.argsort(x)

    x_sorted = x[order]
    y_sorted = y[order]

    # ---------------------------------------------------------
    # robust LOWESS
    # ---------------------------------------------------------
    smoothed = lowess(
        y_sorted,
        x_sorted,
        frac=frac,
        it=3,
        return_sorted=False,
    )

    smoothed = np.asarray(
        smoothed,
        dtype=np.float64,
    )

    if not np.all(np.isfinite(smoothed)):
        logger.warning(
            "LOWESS 결과에 non-finite 값 존재 - GC correction 스킵"
        )
        return counts.copy()

    # ---------------------------------------------------------
    # 모든 GC 위치에서 predicted trend
    #
    # np.interp:
    #   GC fit 범위 밖에서는 endpoint 값을 사용.
    #   즉 극단 GC에서 무리한 extrapolation을 하지 않는다.
    # ---------------------------------------------------------
    predicted_log = np.interp(
        gc,
        x_sorted,
        smoothed,
    )

    center = float(
        np.median(
            predicted_log[fit_mask]
        )
    )

    if not np.isfinite(center):
        logger.warning(
            "GC predicted center 계산 실패 - GC correction 스킵"
        )
        return counts.copy()

    # ---------------------------------------------------------
    # GC trend 제거
    # ---------------------------------------------------------
    corrected_log = (
        log_counts
        - predicted_log
        + center
    )

    corrected = (
        np.exp(corrected_log)
        - pseudocount
    )

    corrected = np.maximum(
        corrected,
        0.0,
    )

    # ---------------------------------------------------------
    # correction ratio safety
    # ---------------------------------------------------------
    ratio = np.ones(
        len(counts),
        dtype=np.float64,
    )

    ratio[valid] = (
        corrected[valid]
        / counts[valid]
    )

    outlier_mask = (
        valid
        & (
            (ratio > max_correction_ratio)
            | (
                ratio
                < 1.0 / max_correction_ratio
            )
        )
    )

    n_outlier = int(
        outlier_mask.sum()
    )

    if n_outlier > 0:
        logger.warning(
            "GC correction 이상치 %d bin -> raw fallback",
            n_outlier,
        )

        corrected[outlier_mask] = (
            counts[outlier_mask]
        )

    # 0 / invalid count는 원래 값 유지
    corrected[~valid] = counts[~valid]

    return corrected


# =====================================================================
# Fraction용 LOESS correction
# =====================================================================

def _loess_fraction_correct(
    values: np.ndarray,
    gc: np.ndarray,
    frac: float = 0.45,
    gc_trim_pct: float = 0.01,
) -> np.ndarray:
    """
    0~1 fraction 값용 GC correction.

    count와 달리 multiplicative correction을 하지 않고:

        corrected
            = observed
            - predicted_by_gc
            + median(predicted)

    additive residual correction 수행.

    최종적으로 0~1 범위로 clip.
    """

    values = np.asarray(
        values,
        dtype=np.float64,
    )

    gc = np.asarray(
        gc,
        dtype=np.float64,
    )

    valid = (
        np.isfinite(values)
        & np.isfinite(gc)
    )

    if valid.sum() < 10:
        logger.warning(
            "fraction 유효 bin 부족 - GC correction 스킵"
        )
        return values.copy()

    fit_mask, gc_lo, gc_hi = _build_gc_fit_mask(
        values,
        gc,
        gc_trim_pct=gc_trim_pct,
        require_positive=False,
    )

    if fit_mask.sum() < 10:
        logger.warning(
            "fraction GC trim 후 적합 bin 부족 - raw 사용"
        )
        return values.copy()

    x = gc[fit_mask]
    y = values[fit_mask]

    order = np.argsort(x)

    x_sorted = x[order]
    y_sorted = y[order]

    smoothed = lowess(
        y_sorted,
        x_sorted,
        frac=frac,
        it=3,
        return_sorted=False,
    )

    smoothed = np.asarray(
        smoothed,
        dtype=np.float64,
    )

    if not np.all(np.isfinite(smoothed)):
        logger.warning(
            "fraction LOWESS 결과 non-finite - raw 사용"
        )
        return values.copy()

    predicted = np.interp(
        gc,
        x_sorted,
        smoothed,
    )

    center = float(
        np.median(
            predicted[fit_mask]
        )
    )

    corrected = values.copy()

    corrected[valid] = (
        values[valid]
        - predicted[valid]
        + center
    )

    corrected[valid] = np.clip(
        corrected[valid],
        0.0,
        1.0,
    )

    return corrected


# =====================================================================
# Diagnostic logging
# =====================================================================

def _log_gc_correction_by_chrom(
    chroms: np.ndarray,
    gc: np.ndarray,
    raw: np.ndarray,
    corrected: np.ndarray,
    label: str,
) -> None:
    """
    chromosome별 GC correction 상태 출력.

    chr19 같은 특정 chromosome에서 systematic하게
    correction factor가 낮아지는지 확인하기 위한 diagnostic.
    """

    logger.info(
        "%s GC correction chromosome diagnostics:",
        label,
    )

    chrom_order = [
        f"chr{i}" for i in range(1, 23)
    ] + [
        "chrX",
        "chrY",
    ]

    present = set(
        str(c)
        for c in chroms
    )

    ordered = [
        c for c in chrom_order
        if c in present
    ]

    # chr prefix 없는 데이터 대응
    if not ordered:
        chrom_order = [
            str(i) for i in range(1, 23)
        ] + [
            "X",
            "Y",
        ]

        ordered = [
            c for c in chrom_order
            if c in present
        ]

    # 예상 외 contig도 마지막에 추가
    extras = sorted(
        present - set(ordered)
    )

    ordered.extend(extras)

    for chrom in ordered:

        mask = (
            chroms.astype(str)
            == str(chrom)
        )

        valid = (
            mask
            & np.isfinite(raw)
            & np.isfinite(corrected)
            & np.isfinite(gc)
            & (raw > 0)
        )

        if valid.sum() == 0:
            continue

        raw_med = float(
            np.median(raw[valid])
        )

        corrected_med = float(
            np.median(corrected[valid])
        )

        gc_med = float(
            np.median(gc[valid])
        )

        factor = (
            corrected_med / raw_med
            if raw_med > 0
            else np.nan
        )

        logger.info(
            "  %-5s n=%5d GC=%.4f raw_med=%10.3f "
            "gc_med=%10.3f factor=%.4f",
            chrom,
            int(valid.sum()),
            gc_med,
            raw_med,
            corrected_med,
            factor,
        )


# =====================================================================
# 메인
# =====================================================================

def run_normalize(
    bin_tsv_path: str,
    output_path: str,
    excluded_chroms: set[str] | None = None,
    gc_correction: bool = True,
    gc_loess_span: float = 0.45,
    gc_trim_pct: float = 0.01,
    gc_pseudocount: float = 0.5,
) -> None:
    """
    GC correction만 수행.
    Z-score는 within_sample_cnv.py에서 WSS 방식으로 계산.

    Parameters
    ----------
    bin_tsv_path
        run_binning 출력 TSV.

    output_path
        GC correction 결과 TSV.

    excluded_chroms
        하위 호환용. 현재 GC fitting에서 별도 제외하지 않음.

    gc_correction
        LOESS GC correction 수행 여부.

    gc_loess_span
        LOWESS frac.
        기본값 0.45.

        low-pass NIPT + 후단 WSS를 고려해
        너무 local한 correction보다 smooth한 GC trend 제거를 목표로 한다.

    gc_trim_pct
        GC 양 끝 trim 비율.
        기본 0.01 = lower/upper 1%.

    gc_pseudocount
        log(count + pseudocount)에 사용할 값.
        기본 0.5.
    """

    rows = read_tsv(
        bin_tsv_path
    )

    if not rows:
        raise ValueError(
            f"빈 TSV: {bin_tsv_path}"
        )

    # ---------------------------------------------------------
    # input
    # ---------------------------------------------------------
    chroms = np.array(
        [
            r["chrom"]
            for r in rows
        ],
        dtype=object,
    )

    bin_starts = np.array(
        [
            int(r["bin_start"])
            for r in rows
        ],
        dtype=np.int64,
    )

    bin_ends = np.array(
        [
            int(r["bin_end"])
            for r in rows
        ],
        dtype=np.int64,
    )

    gc = np.array(
        [
            float(r["gc"])
            for r in rows
        ],
        dtype=np.float64,
    )

    mappability = np.array(
        [
            float(r["mappability"])
            for r in rows
        ],
        dtype=np.float64,
    )

    total_raw = np.array(
        [
            int(r["total"])
            for r in rows
        ],
        dtype=np.float64,
    )

    short_raw = np.array(
        [
            int(r["short"])
            for r in rows
        ],
        dtype=np.float64,
    )

    long_raw = np.array(
        [
            int(r["long"])
            for r in rows
        ],
        dtype=np.float64,
    )

    # ---------------------------------------------------------
    # GC correction
    # ---------------------------------------------------------
    if gc_correction:

        logger.info(
            "LOESS GC correction "
            "(span=%.3f, trim=%.3f, pseudocount=%.3f)",
            gc_loess_span,
            gc_trim_pct,
            gc_pseudocount,
        )

        total_gc = _loess_gc_correct(
            total_raw,
            gc,
            frac=gc_loess_span,
            gc_trim_pct=gc_trim_pct,
            pseudocount=gc_pseudocount,
        )

        short_gc = _loess_gc_correct(
            short_raw,
            gc,
            frac=gc_loess_span,
            gc_trim_pct=gc_trim_pct,
            pseudocount=gc_pseudocount,
        )

        long_gc = _loess_gc_correct(
            long_raw,
            gc,
            frac=gc_loess_span,
            gc_trim_pct=gc_trim_pct,
            pseudocount=gc_pseudocount,
        )

        # -----------------------------------------------------
        # chromosome diagnostic
        # -----------------------------------------------------
        _log_gc_correction_by_chrom(
            chroms,
            gc,
            total_raw,
            total_gc,
            "TOTAL",
        )

        _log_gc_correction_by_chrom(
            chroms,
            gc,
            short_raw,
            short_gc,
            "SHORT",
        )

        _log_gc_correction_by_chrom(
            chroms,
            gc,
            long_raw,
            long_gc,
            "LONG",
        )

    else:

        logger.info(
            "GC correction disabled"
        )

        total_gc = total_raw.copy()
        short_gc = short_raw.copy()
        long_gc = long_raw.copy()

    # ---------------------------------------------------------
    # short fraction
    #
    # 각각 short/total을 독립적으로 multiplicative correction 후
    # 나누지 않는다.
    #
    # raw fraction 자체의 GC-dependent trend를 additive 방식으로 제거.
    # ---------------------------------------------------------
    short_fraction_raw = np.full(
        len(total_raw),
        np.nan,
        dtype=np.float64,
    )

    total_positive = (
        total_raw > 0
    )

    short_fraction_raw[total_positive] = (
        short_raw[total_positive]
        / total_raw[total_positive]
    )

    if (
        gc_correction
        and np.isfinite(
            short_fraction_raw
        ).sum() >= 10
    ):

        short_fraction_gc = (
            _loess_fraction_correct(
                short_fraction_raw,
                gc,
                frac=gc_loess_span,
                gc_trim_pct=gc_trim_pct,
            )
        )

    else:

        if gc_correction:
            logger.warning(
                "short_fraction GC correction 불가 - raw 사용"
            )

        short_fraction_gc = (
            short_fraction_raw.copy()
        )

    # ---------------------------------------------------------
    # 결과 저장
    # ---------------------------------------------------------
    out_rows = []

    for i in range(
        len(rows)
    ):

        out_rows.append(
            {
                "chrom": chroms[i],

                "bin_start": int(
                    bin_starts[i]
                ),

                "bin_end": int(
                    bin_ends[i]
                ),

                "gc": (
                    f"{gc[i]:.4f}"
                ),

                "mappability": (
                    f"{mappability[i]:.4f}"
                ),

                "total_raw": int(
                    total_raw[i]
                ),

                "short_raw": int(
                    short_raw[i]
                ),

                "long_raw": int(
                    long_raw[i]
                ),

                "total_gc": (
                    f"{total_gc[i]:.4f}"
                ),

                "short_gc": (
                    f"{short_gc[i]:.4f}"
                ),

                "long_gc": (
                    f"{long_gc[i]:.4f}"
                ),

                "short_fraction_gc": (
                    f"{short_fraction_gc[i]:.6f}"
                    if np.isfinite(
                        short_fraction_gc[i]
                    )
                    else "NA"
                ),
            }
        )

    write_tsv(
        out_rows,
        output_path,
        FIELDNAMES,
    )

    logger.info(
        "GC correction 완료 -> %s (%d bins)",
        output_path,
        len(out_rows),
    )