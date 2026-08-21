"""
fetal_fraction.py
=================

Fetal Fraction 추정

1. short / long fragment ratio 기반 FF
2. chrY fraction 기반 male-fetus FF

중요:
- Y_MALE_THRESHOLD는 성별 판정 threshold
- chrY fraction -> FF 변환은 별도 calibration 필요
- ff = 2 * y_fraction 같은 단순식은 사용하지 않음
"""

from __future__ import annotations

import json
import logging
import os

import numpy as np
import pandas as pd

from nipt_fragmentomics.core.constants import (
    SEQFF_ALPHA,
    SEQFF_BETA,
    Y_MALE_THRESHOLD,
)

log = logging.getLogger(__name__)

SEX_CHROMS = {"chrX", "chrY"}

FMY_EXCLUDE_REGIONS = [
    (0, 2_781_479),
    (13_000_000, 15_000_000),
    (21_000_000, 21_699_574),
    (26_637_971, 30_000_000),
    (56_887_903, 57_217_415),
]

# 반드시 자체 male-reference cohort로 보정하는 값.
# 아래 값은 예시 기본값이며 기존 pipeline calibration 값이 있으면 그것을 사용해야 함.
Y_FF_SCALE_DEFAULT = 1000.0


def _is_fmy_region(start: int, end: int) -> bool:
    for rs, re in FMY_EXCLUDE_REGIONS:
        if start < re and end > rs:
            return True
    return False


def estimate_seqff(
    df: pd.DataFrame,
    alpha: float = SEQFF_ALPHA,
    beta: float = SEQFF_BETA,
) -> tuple[float, float]:
    SEQFF_EXCL = {"chr13", "chr18", "chr21", "chrX", "chrY"}

    mask = ~df["chrom"].isin(SEQFF_EXCL)

    if "mappability_pass" in df.columns:
        mask &= df["mappability_pass"].fillna(False).astype(bool)

    auto = df.loc[mask]

    if auto.empty:
        return float("nan"), float("nan")

    short_sum = float(auto["short_count"].fillna(0).sum())
    long_sum = float(auto["long_count"].fillna(0).sum())

    if long_sum <= 0:
        return float("nan"), float("nan")

    sl_ratio = short_sum / long_sum
    ff = alpha + beta * sl_ratio

    if not np.isfinite(ff):
        return float("nan"), float(sl_ratio)

    return float(np.clip(ff, 0.0, 1.0)), float(sl_ratio)


def estimate_y_chr(
    df: pd.DataFrame,
    threshold: float = Y_MALE_THRESHOLD,
    y_ff_scale: float = Y_FF_SCALE_DEFAULT,
) -> tuple[float, float, bool]:
    y_df = df[df["chrom"] == "chrY"].copy()

    if y_df.empty:
        log.warning("chrY bin 없음")
        return float("nan"), float("nan"), False

    fmy_mask = y_df.apply(
        lambda r: _is_fmy_region(int(r["start"]), int(r["end"])),
        axis=1,
    )

    y_clean = y_df.loc[~fmy_mask]

    if y_clean.empty:
        log.warning("chrY FMY 필터 후 유효 bin 없음")
        return float("nan"), float("nan"), False

    # denominator는 동일한 count definition으로 계산
    total_reads = float(
        df["short_count"].fillna(0).sum()
        + df["long_count"].fillna(0).sum()
    )

    if total_reads <= 0:
        return float("nan"), float("nan"), False

    y_reads = float(
        y_clean["short_count"].fillna(0).sum()
        + y_clean["long_count"].fillna(0).sum()
    )

    y_fraction = y_reads / total_reads
    is_male = bool(y_fraction > threshold)

    log.info(
        "ChrY fraction=%.8e threshold=%.8e Y_reads=%.0f total_reads=%.0f male=%s",
        y_fraction,
        threshold,
        y_reads,
        total_reads,
        is_male,
    )

    if not is_male:
        return float("nan"), float(y_fraction), False

    # calibration:
    # 예) y_fraction=1e-4, scale=1000 -> FF=0.10
    y_ff = float(np.clip(y_fraction * y_ff_scale, 0.0, 1.0))

    return y_ff, float(y_fraction), True


def run(
    corrected_path: str,
    out_path: str,
    alpha: float = SEQFF_ALPHA,
    beta: float = SEQFF_BETA,
    y_ff_scale: float = Y_FF_SCALE_DEFAULT,
) -> dict:
    df = pd.read_parquet(corrected_path)

    required = {
        "chrom",
        "start",
        "end",
        "short_count",
        "long_count",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"fetal_fraction 입력 컬럼 없음: {sorted(missing)}"
        )

    seqff_ff, sl_ratio = estimate_seqff(
        df,
        alpha=alpha,
        beta=beta,
    )

    y_ff, y_fraction, is_male = estimate_y_chr(
        df,
        threshold=Y_MALE_THRESHOLD,
        y_ff_scale=y_ff_scale,
    )

    log.info(
        "SeqFF: FF=%s short/long=%.6f",
        f"{seqff_ff:.4f}" if np.isfinite(seqff_ff) else "NA",
        sl_ratio if np.isfinite(sl_ratio) else float("nan"),
    )

    if is_male:
        log.info(
            "Y-based: FF=%.4f Y_fraction=%.8e scale=%.4f",
            y_ff,
            y_fraction,
            y_ff_scale,
        )
    else:
        log.info(
            "Y-based: female/undetermined Y_fraction=%.8e",
            y_fraction if np.isfinite(y_fraction) else float("nan"),
        )

    # 기존 결과 안정성을 위해:
    # male에서도 SeqFF를 버리지 않고 둘 다 저장.
    # consensus는 calibration이 완료된 Y FF가 있을 때 Y 우선.
    if is_male and np.isfinite(y_ff):
        consensus_ff = y_ff
        methods = ["y_chr"]

        if np.isfinite(seqff_ff):
            methods.append("seqff")

    elif np.isfinite(seqff_ff):
        consensus_ff = seqff_ff
        methods = ["seqff"]

    else:
        consensus_ff = float("nan")
        methods = []

    def _safe(v, digits=6):
        return round(float(v), digits) if np.isfinite(v) else None

    result = {
        "seqff_ff": _safe(seqff_ff, 4),
        "seqff_sl_ratio": _safe(sl_ratio, 6),
        "y_ff": _safe(y_ff, 4),
        "y_fraction": _safe(y_fraction, 8),
        "y_chr_pct": _safe(y_fraction, 8),
        "y_ff_scale": float(y_ff_scale),
        "is_male_fetus": bool(is_male),
        "consensus_ff": _safe(consensus_ff, 4),
        "method_used": methods,
    }

    os.makedirs(
        os.path.dirname(os.path.abspath(out_path)),
        exist_ok=True,
    )

    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    log.info(
        "fetal_fraction 저장: %s consensus=%s male=%s",
        out_path,
        result["consensus_ff"],
        is_male,
    )

    return result