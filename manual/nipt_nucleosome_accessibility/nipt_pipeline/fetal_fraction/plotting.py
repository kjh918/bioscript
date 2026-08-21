"""
fetal_fraction/plotting.py
==========================
bin-level short_gc vs long_gc scatter + SeqFF 선형회귀선

- X축: long_gc  (GC-corrected long fragment count, ≥150bp)
- Y축: short_gc (GC-corrected short fragment count, <150bp)
- 선형회귀선: numpy polyfit (degree=1)
- SeqFF FF 추정값 / R² / 기울기 타이틀/범례 표기
- reference bin (aneuploidy 염색체 제외) 만 사용
"""
from __future__ import annotations

import logging
import os
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from cnv.utils import read_tsv

log = logging.getLogger(__name__)

_EXCLUDED = {"chr15", "chr18", "chr21", "chrX", "chrY"}


def plot_ff_correlation(
    normalized_tsv: str,
    output_path:    str,
    sample_id:      str           = "",
    ff_seqff:       Optional[float] = None,
    ff_chry:        Optional[float] = None,
) -> None:
    """
    bin-level short_gc vs long_gc scatter + 선형회귀선

    Parameters
    ----------
    normalized_tsv : run_normalize 출력 TSV
    output_path    : PNG 저장 경로
    sample_id      : 타이틀 표기용
    ff_seqff       : SeqFF FF 추정값 (타이틀 표기용, 없으면 생략)
    ff_chry        : chrY FF 추정값 (타이틀 표기용, 없으면 생략)
    """
    rows = read_tsv(normalized_tsv)
    if not rows:
        log.warning("normalized TSV 비어있음: %s", normalized_tsv)
        return

    long_vals:  list[float] = []
    short_vals: list[float] = []
    chrom_labels: list[str] = []

    for r in rows:
        chrom    = r["chrom"]
        total_gc = float(r.get("total_gc", 0) or 0)
        short_gc = float(r.get("short_gc",  0) or 0)
        long_gc  = float(r.get("long_gc",   0) or 0)

        # reference bin만 사용 (aneuploidy 염색체 제외, total_gc > 0)
        if chrom in _EXCLUDED or total_gc <= 0:
            continue

        long_vals.append(long_gc)
        short_vals.append(short_gc)
        chrom_labels.append(chrom)

    if not long_vals:
        log.warning("유효한 reference bin 없음 — FF correlation plot 스킵")
        return

    x = np.array(long_vals,  dtype=np.float32)
    y = np.array(short_vals, dtype=np.float32)

    # ── 선형회귀 ──────────────────────────────────────────────────────
    # short = slope * long + intercept
    coeffs = np.polyfit(x, y, deg=1)
    slope, intercept = float(coeffs[0]), float(coeffs[1])

    x_line = np.linspace(float(x.min()), float(x.max()), 500)
    y_line = slope * x_line + intercept

    # R² 계산
    y_pred = slope * x + intercept
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    r2     = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # ── 플롯 ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    ax.set_facecolor("white")

    # scatter (점 밀도가 높으므로 작게 + alpha)
    ax.scatter(
        x, y,
        s       = 3,
        alpha   = 0.35,
        color   = "#1E88E5",
        linewidths = 0,
        rasterized = True,
        label   = f"Bins (n={len(x):,})",
        zorder  = 2,
    )

    # 선형회귀선
    ax.plot(
        x_line, y_line,
        color     = "#E53935",
        linewidth = 1.8,
        label     = f"Linear fit  slope={slope:.4f}  intercept={intercept:.2f}",
        zorder    = 3,
    )

    # 축 레이블
    ax.set_xlabel("Long fragment count (GC-corrected, ≥150bp)", fontsize=10)
    ax.set_ylabel("Short fragment count (GC-corrected, <150bp)", fontsize=10)

    # 타이틀: sample_id + FF 값들
    ff_parts = []
    if ff_seqff is not None:
        ff_parts.append(f"SeqFF={ff_seqff:.2f}%")
    if ff_chry is not None:
        ff_parts.append(f"chrY FF={ff_chry:.2f}%")
    ff_str   = "  |  ".join(ff_parts) if ff_parts else ""
    title    = f"{sample_id} — Short vs Long bin count correlation"
    subtitle = f"R²={r2:.4f}  {ff_str}"
    ax.set_title(f"{title}\n{subtitle}", fontsize=10, fontweight="bold")

    ax.legend(fontsize=8, framealpha=0.8, loc="upper left")
    ax.tick_params(labelsize=8)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.spines[["top", "right"]].set_visible(False)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("FF correlation plot 저장: %s (R²=%.4f)", output_path, r2)