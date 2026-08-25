"""
cnv/calculate_cnv.py
====================
WSS normalized bin TSV → 염색체 단위 CNV call TSV

입력: wss_cnv.normalized_bins.tsv
  컬럼: wss_z_total, wss_z_short, wss_z_long, wss_z_short_fraction

출력 컬럼:
  chrom, n_bins,
  median_zscore_total, median_zscore_short, median_zscore_long,
  median_zscore_short_fraction,
  call_total, call_short, call_long, call_short_fraction
"""
from __future__ import annotations

import logging
from collections import defaultdict

import numpy as np

from .utils import read_tsv, write_tsv, STANDARD_CHROMS

logger = logging.getLogger(__name__)

CALL_GAIN         = "GAIN"
CALL_LOSS         = "LOSS"
CALL_NORMAL       = "NORMAL"
CALL_INSUFFICIENT = "INSUFFICIENT_BINS"

FIELDNAMES = [
    "chrom", "n_bins",
    "median_zscore_total", "median_zscore_short", "median_zscore_long",
    "median_zscore_short_fraction",
    "call_total", "call_short", "call_long", "call_short_fraction",
]

# normalized_bins.tsv 컬럼 → 내부 키 매핑
_COL_MAP = {
    "total":          "wss_z_total",
    "short":          "wss_z_short",
    "long":           "wss_z_long",
    "short_fraction": "wss_z_short_fraction",
}


def _call(z: float | None, gain_thr: float, loss_thr: float) -> str:
    if z is None:
        return CALL_INSUFFICIENT
    if z >= gain_thr:
        return CALL_GAIN
    if z <= loss_thr:
        return CALL_LOSS
    return CALL_NORMAL


def run_calculate_cnv(
    norm_tsv_path:         str,
    output_path:           str,
    zscore_gain_threshold: float = 3.0,
    zscore_loss_threshold: float = -3.0,
    min_bins_per_chrom:    int   = 5,
    target_chroms:         list[str] | None = None,
) -> None:
    """
    Parameters
    ----------
    norm_tsv_path         : wss_cnv.normalized_bins.tsv (wss_z_* 컬럼)
    output_path           : 염색체 단위 CNV call TSV
    zscore_gain_threshold : GAIN threshold
    zscore_loss_threshold : LOSS threshold
    min_bins_per_chrom    : 최소 유효 bin 수
    target_chroms         : 분석 염색체 (None이면 standard 전체)
    """
    rows = read_tsv(norm_tsv_path)
    if not rows:
        raise ValueError(f"빈 TSV: {norm_tsv_path}")

    chrom_zscores: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {"total": [], "short": [], "long": [], "short_fraction": []}
    )

    for r in rows:
        chrom = r.get("chrom", "")
        for key, col in _COL_MAP.items():
            val = r.get(col, "NA")
            if val not in ("NA", "", None):
                try:
                    chrom_zscores[chrom][key].append(float(val))
                except ValueError:
                    pass

    target   = target_chroms or STANDARD_CHROMS
    out_rows = []

    for chrom in target:
        data   = chrom_zscores.get(chrom, {
            "total": [], "short": [], "long": [], "short_fraction": []
        })
        n_bins = len(data["total"])

        if n_bins < min_bins_per_chrom:
            logger.warning("%s: bin %d개 < 최소 %d → INSUFFICIENT_BINS",
                           chrom, n_bins, min_bins_per_chrom)
            out_rows.append({
                "chrom":                        chrom,
                "n_bins":                       n_bins,
                "median_zscore_total":          "NA",
                "median_zscore_short":          "NA",
                "median_zscore_long":           "NA",
                "median_zscore_short_fraction": "NA",
                "call_total":                   CALL_INSUFFICIENT,
                "call_short":                   CALL_INSUFFICIENT,
                "call_long":                    CALL_INSUFFICIENT,
                "call_short_fraction":          CALL_INSUFFICIENT,
            })
            continue

        med_total = float(np.median(data["total"]))
        med_short = float(np.median(data["short"]))          if data["short"]          else None
        med_long  = float(np.median(data["long"]))           if data["long"]           else None
        med_sfrac = float(np.median(data["short_fraction"])) if data["short_fraction"] else None

        out_rows.append({
            "chrom":                        chrom,
            "n_bins":                       n_bins,
            "median_zscore_total":          f"{med_total:.4f}",
            "median_zscore_short":          f"{med_short:.4f}" if med_short is not None else "NA",
            "median_zscore_long":           f"{med_long:.4f}"  if med_long  is not None else "NA",
            "median_zscore_short_fraction": f"{med_sfrac:.4f}" if med_sfrac is not None else "NA",
            "call_total":                   _call(med_total, zscore_gain_threshold, zscore_loss_threshold),
            "call_short":                   _call(med_short, zscore_gain_threshold, zscore_loss_threshold),
            "call_long":                    _call(med_long,  zscore_gain_threshold, zscore_loss_threshold),
            "call_short_fraction":          _call(med_sfrac, zscore_gain_threshold, zscore_loss_threshold),
        })

    write_tsv(out_rows, output_path, FIELDNAMES)
    logger.info("염색체 단위 CNV call 완료 → %s (%d chroms)", output_path, len(out_rows))