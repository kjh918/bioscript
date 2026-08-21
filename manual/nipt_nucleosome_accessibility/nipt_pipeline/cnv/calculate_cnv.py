"""
cnv/calculate_cnv.py

normalized TSV → 염색체 단위 CNV call TSV

- 염색체 내 bin Z-score의 median으로 집계
- total / short / long 각각 독립 call
- 최소 bin 수 미달 시 INSUFFICIENT_BINS

출력 컬럼:
  chrom,
  n_bins,
  median_zscore_total, median_zscore_short, median_zscore_long,
  call_total, call_short, call_long
"""
from __future__ import annotations

import logging
from collections import defaultdict

import numpy as np

from .utils import read_tsv, write_tsv, STANDARD_CHROMS

logger = logging.getLogger(__name__)

CALL_GAIN        = "GAIN"
CALL_LOSS        = "LOSS"
CALL_NORMAL      = "NORMAL"
CALL_INSUFFICIENT = "INSUFFICIENT_BINS"

FIELDNAMES = [
    "chrom",
    "n_bins",
    "median_zscore_total", "median_zscore_short", "median_zscore_long",
    "call_total", "call_short", "call_long",
]


def _call(zscore: float | None, gain_thr: float, loss_thr: float) -> str:
    if zscore is None:
        return CALL_INSUFFICIENT
    if zscore >= gain_thr:
        return CALL_GAIN
    if zscore <= loss_thr:
        return CALL_LOSS
    return CALL_NORMAL


def run_calculate_cnv(
    norm_tsv_path: str,
    output_path: str,
    zscore_gain_threshold: float = 3.0,
    zscore_loss_threshold: float = -3.0,
    min_bins_per_chrom: int = 5,
    target_chroms: list[str] | None = None,
) -> None:
    """
    Parameters
    ----------
    norm_tsv_path         : run_normalize 출력 TSV
    output_path           : CNV call 결과 TSV
    zscore_gain_threshold : GAIN 판정 Z-score 상한
    zscore_loss_threshold : LOSS 판정 Z-score 하한
    min_bins_per_chrom    : call 산출을 위한 최소 유효 bin 수
    target_chroms         : call을 산출할 염색체 목록 (None이면 standard 전체)
    """
    rows = read_tsv(norm_tsv_path)
    if not rows:
        raise ValueError(f"빈 TSV: {norm_tsv_path}")

    # chrom별 Z-score 수집
    chrom_zscores: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {"total": [], "short": [], "long": []}
    )

    for r in rows:
        chrom = r["chrom"]
        for key in ("total", "short", "long"):
            val = r.get(f"zscore_{key}", "NA")
            if val != "NA":
                chrom_zscores[chrom][key].append(float(val))

    target = target_chroms or STANDARD_CHROMS
    out_rows = []

    for chrom in target:
        data = chrom_zscores.get(chrom, {"total": [], "short": [], "long": []})
        n_bins = len(data["total"])  # total 기준 (short/long도 동일 bin 수)

        if n_bins < min_bins_per_chrom:
            logger.warning("%s: 유효 bin %d개 < 최소 %d → INSUFFICIENT_BINS", chrom, n_bins, min_bins_per_chrom)
            out_rows.append({
                "chrom":                chrom,
                "n_bins":               n_bins,
                "median_zscore_total":  "NA",
                "median_zscore_short":  "NA",
                "median_zscore_long":   "NA",
                "call_total":           CALL_INSUFFICIENT,
                "call_short":           CALL_INSUFFICIENT,
                "call_long":            CALL_INSUFFICIENT,
            })
            continue

        med_total = float(np.median(data["total"]))
        med_short = float(np.median(data["short"])) if data["short"] else None
        med_long  = float(np.median(data["long"]))  if data["long"]  else None

        out_rows.append({
            "chrom":               chrom,
            "n_bins":              n_bins,
            "median_zscore_total": f"{med_total:.4f}",
            "median_zscore_short": f"{med_short:.4f}" if med_short is not None else "NA",
            "median_zscore_long":  f"{med_long:.4f}"  if med_long  is not None else "NA",
            "call_total":          _call(med_total, zscore_gain_threshold, zscore_loss_threshold),
            "call_short":          _call(med_short, zscore_gain_threshold, zscore_loss_threshold),
            "call_long":           _call(med_long,  zscore_gain_threshold, zscore_loss_threshold),
        })

    write_tsv(out_rows, output_path, FIELDNAMES)
    logger.info("CNV call 완료 → %s", output_path)
