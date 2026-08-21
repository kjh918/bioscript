"""
cnv/plotting.py

normalized TSV + cnv call TSV → CNV Z-score scatter plot

기능
----
1. plot_cnv_scatter()     : genome-wide scatter (short / long 각각 1 PNG)
2. plot_cnv_scatter_by_chrom(): chromosome별 개별 PNG

공통 사항
---------
- bin_start 기준 오름차순 정렬 보장
- X축: bin position (Mbp)
- Y축: Z-score (IQR 기반 robust ylim)
- 색상: GAIN(빨강) / LOSS(파랑) / NORMAL(회색) — 염색체 단위 call
- median Z-score 수평선 오버레이
- gain/loss threshold 점선
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

from .utils import read_tsv, STANDARD_CHROMS

logger = logging.getLogger(__name__)

# region CNV call 색상 (bin-level과 구분되는 진한 색)
REGION_COLOR: dict[str, str] = {
    "GAIN": "#B71C1C",   # 진한 빨강
    "LOSS": "#0D47A1",   # 진한 파랑
}

CHROM_ORDER: list[str] = STANDARD_CHROMS

CALL_COLOR: dict[str, str] = {
    "GAIN":              "#D32F2F",
    "LOSS":              "#1565C0",
    "NORMAL":            "#BDBDBD",
    "INSUFFICIENT_BINS": "#9E9E9E",
}
BG_COLORS = ["#F5F5F5", "#EEEEEE"]


# ── 데이터 로더 ────────────────────────────────────────────────────────

def _load_region_cnv(
    region_cnv_path: Optional[str],
    call_col: str = "call_total",
) -> dict[str, list[tuple[int, int, str]]]:
    """
    region_cnv.tsv → chrom: [(start, end, call), ...] 딕셔너리.
    GAIN / LOSS region만 반환.
    """
    result: dict[str, list[tuple[int, int, str]]] = {}
    if not region_cnv_path:
        return result
    from pathlib import Path as _Path
    if not _Path(region_cnv_path).exists():
        return result
    rows = read_tsv(region_cnv_path)
    for r in rows:
        call = r.get(call_col, "NORMAL")
        if call not in ("GAIN", "LOSS"):
            continue
        chrom = r.get("chrom", "")
        try:
            start = int(r["start"])
            end   = int(r["end"])
        except (KeyError, ValueError):
            continue
        result.setdefault(chrom, []).append((start, end, call))
    return result
def _load_bin_data(
    norm_tsv_path: str,
    zscore_col: str,
) -> dict[str, list[tuple[int, float]]]:
    """
    chrom → [(bin_start, zscore), ...] 정렬된 리스트.
    bin_start 오름차순 정렬 보장.
    """
    rows = read_tsv(norm_tsv_path)
    result: dict[str, list[tuple[int, float]]] = {c: [] for c in CHROM_ORDER}

    for r in rows:
        chrom = r.get("chrom", "")
        if chrom not in result:
            continue
        val = r.get(zscore_col, "NA")
        if val == "NA":
            continue
        try:
            result[chrom].append((int(r["bin_start"]), float(val)))
        except (ValueError, KeyError):
            continue

    # bin_start 기준 정렬
    for chrom in result:
        result[chrom].sort(key=lambda x: x[0])

    return result


def _load_cnv_calls(cnv_tsv_path: str, call_col: str) -> dict[str, str]:
    rows = read_tsv(cnv_tsv_path)
    return {r["chrom"]: r.get(call_col, "NORMAL") for r in rows}


def _load_median_zscores(cnv_tsv_path: str, median_col: str) -> dict[str, Optional[float]]:
    rows = read_tsv(cnv_tsv_path)
    result: dict[str, Optional[float]] = {}
    for r in rows:
        val = r.get(median_col, "NA")
        result[r["chrom"]] = float(val) if val != "NA" else None
    return result


def _robust_ylim(
    y_vals: list[float],
    gain_thr: float,
    loss_thr: float,
) -> tuple[float, float]:
    """IQR 기반 robust ylim — 이상치 bin이 스케일 압축하는 문제 방지"""
    if not y_vals:
        return (loss_thr - 1.0, gain_thr + 1.0)
    arr = np.array(y_vals)
    q25, q75 = np.percentile(arr, [25, 75])
    iqr = q75 - q25 or 1.0
    med = float(np.median(arr))
    ymin = min(med - 3.5 * iqr, loss_thr - 0.5)
    ymax = max(med + 3.5 * iqr, gain_thr + 0.5)
    return ymin, ymax


def _make_legend(gain_thr: float, loss_thr: float) -> list:
    return [
        mpatches.Patch(color=CALL_COLOR["GAIN"],   label=f"GAIN  (Z ≥ {gain_thr})"),
        mpatches.Patch(color=CALL_COLOR["LOSS"],   label=f"LOSS  (Z ≤ {loss_thr})"),
        mpatches.Patch(color=CALL_COLOR["NORMAL"], label="NORMAL"),
        Line2D([0], [0], color="black", linewidth=1.5, label="Chrom median"),
    ]


# ── genome-wide scatter (1 PNG) ───────────────────────────────────────
def _draw_genome_scatter(
    ax:             plt.Axes,
    bin_data:       dict[str, list[tuple[int, float]]],
    cnv_calls:      dict[str, str],
    median_zscores: dict[str, Optional[float]],
    gain_thr:       float,
    loss_thr:       float,
    fragment_type:  str,
    sample_id:      str,
) -> None:
    ax.set_facecolor("white")

    x_pos:   list[float] = []
    y_vals:  list[float] = []
    colors:  list[str]   = []
    chrom_centers:    dict[str, float] = {}
    chrom_boundaries: list[float]      = [0.0]
    x_cursor = 0.0

    for idx, chrom in enumerate(CHROM_ORDER):
        bins  = bin_data.get(chrom, [])
        call  = cnv_calls.get(chrom, "NORMAL")
        color = CALL_COLOR.get(call, CALL_COLOR["NORMAL"])
        n     = len(bins)

        ax.axvspan(
            x_cursor, x_cursor + max(n, 1),
            color=BG_COLORS[idx % 2], alpha=0.5, zorder=0,
        )

        if n > 0:
            # X: 순서 index (bin_start 정렬 후)
            xs = np.arange(x_cursor, x_cursor + n)
            ys = [v for _, v in bins]
            x_pos.extend(xs.tolist())
            y_vals.extend(ys)
            colors.extend([color] * n)

            med = median_zscores.get(chrom)
            if med is not None:
                ax.hlines(med, x_cursor, x_cursor + n,
                          colors=color, linewidths=1.5, zorder=3)

        chrom_centers[chrom] = x_cursor + max(n, 1) / 2
        x_cursor += max(n, 1)
        chrom_boundaries.append(x_cursor)

    ax.scatter(x_pos, y_vals, c=colors, s=4, alpha=0.7, linewidths=0, zorder=2)
    ax.axhline(gain_thr, color="#D32F2F", linestyle="--", linewidth=0.9, alpha=0.8, zorder=4)
    ax.axhline(loss_thr, color="#1565C0", linestyle="--", linewidth=0.9, alpha=0.8, zorder=4)
    ax.axhline(0,        color="#757575", linestyle="-",  linewidth=0.6, alpha=0.5, zorder=1)

    for xb in chrom_boundaries[1:-1]:
        ax.axvline(xb, color="#E0E0E0", linewidth=0.5, zorder=1)

    ax.set_xticks([chrom_centers[c] for c in CHROM_ORDER if c in chrom_centers])
    ax.set_xticklabels(
        [c.replace("chr", "") for c in CHROM_ORDER if c in chrom_centers],
        fontsize=7, rotation=0,
    )
    ax.set_xlim(0, x_cursor)
    ax.set_xlabel("Chromosome", fontsize=9)
    ax.set_ylabel("Z-score", fontsize=9)
    ax.set_title(f"{sample_id} — {fragment_type} fragment Z-score",
                 fontsize=10, fontweight="bold")

    ymin, ymax = _robust_ylim(y_vals, gain_thr, loss_thr)
    ax.set_ylim(ymin, ymax)
    ax.legend(handles=_make_legend(gain_thr, loss_thr),
              fontsize=7, loc="upper right", framealpha=0.8)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)


# ── chromosome-level scatter (1 PNG per chrom) ────────────────────────
def _draw_chrom_scatter(
    ax:             plt.Axes,
    chrom:          str,
    bins:           list[tuple[int, float]],  # (bin_start, zscore) 정렬됨
    call:           str,
    median_z:       Optional[float],
    gain_thr:       float,
    loss_thr:       float,
    fragment_type:  str,
    sample_id:      str,
    regions:        Optional[list[tuple[int, int, str]]] = None,  # region CNV
) -> None:
    """
    regions: [(start, end, call), ...] — region_cnv.tsv GAIN/LOSS region
    """
    ax.set_facecolor("white")

    if not bins:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=10, color="#9E9E9E")
        return

    xs = np.array([b[0] for b in bins], dtype=np.float64)
    ys = np.array([b[1] for b in bins], dtype=np.float64)

    # ── region CNV 배경 하이라이트 ───────────────────────────────────
    if regions:
        for reg_start, reg_end, reg_call in regions:
            rc = REGION_COLOR.get(reg_call, "#E0E0E0")
            ax.axvspan(
                reg_start / 1e6, reg_end / 1e6,
                color=rc, alpha=0.18, zorder=0,
                label=f"Region {reg_call} ({reg_start//1e6:.1f}-{reg_end//1e6:.1f}Mb)",
            )
            # region 경계 수직선
            for xb in (reg_start / 1e6, reg_end / 1e6):
                ax.axvline(xb, color=rc, linewidth=1.0,
                           linestyle=":", alpha=0.8, zorder=2)

    # ── bin scatter: region 안에 있는 bin은 진한 색으로 ──────────────
    bin_colors = []
    for bstart, bz in bins:
        in_region = False
        if regions:
            for reg_start, reg_end, reg_call in regions:
                if reg_start <= bstart < reg_end:
                    bin_colors.append(REGION_COLOR.get(reg_call, CALL_COLOR["NORMAL"]))
                    in_region = True
                    break
        if not in_region:
            bin_colors.append(CALL_COLOR.get(call, CALL_COLOR["NORMAL"]))

    ax.scatter(xs / 1e6, ys, c=bin_colors, s=5, alpha=0.8, linewidths=0, zorder=3)

    if median_z is not None:
        chrom_color = CALL_COLOR.get(call, "#555555")
        ax.axhline(median_z, color=chrom_color, linewidth=1.5,
                   linestyle="-", zorder=4, label=f"Chrom median Z={median_z:.2f}")

    ax.axhline(gain_thr, color="#D32F2F", linestyle="--",
               linewidth=0.9, alpha=0.8, zorder=5)
    ax.axhline(loss_thr, color="#1565C0", linestyle="--",
               linewidth=0.9, alpha=0.8, zorder=5)
    ax.axhline(0, color="#757575", linestyle="-",
               linewidth=0.6, alpha=0.5, zorder=1)

    ax.set_xlabel(f"{chrom} position (Mbp)", fontsize=9)
    ax.set_ylabel("Z-score", fontsize=9)
    ax.set_title(
        f"{sample_id} — {chrom}  {fragment_type} fragment Z-score  [{call}]",
        fontsize=10, fontweight="bold",
        color=CALL_COLOR.get(call, "#333333"),
    )
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{v:.0f}")
    )
    ax.set_xlim(xs[0] / 1e6, xs[-1] / 1e6)

    ymin, ymax = _robust_ylim(list(ys), gain_thr, loss_thr)
    ax.set_ylim(ymin, ymax)

    # 범례 (region 있으면 추가)
    legend_handles = _make_legend(gain_thr, loss_thr)
    if regions:
        legend_handles += [
            mpatches.Patch(color=REGION_COLOR["GAIN"], alpha=0.5, label="Region GAIN"),
            mpatches.Patch(color=REGION_COLOR["LOSS"], alpha=0.5, label="Region LOSS"),
        ]
    ax.legend(handles=legend_handles, fontsize=7, loc="upper right", framealpha=0.8)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def plot_cnv_scatter(
    norm_tsv_path:         str,
    cnv_tsv_path:          str,
    output_short:          str,
    output_long:           str,
    output_combined:       str,
    sample_id:             str,
    zscore_gain_threshold: float = 3.0,
    zscore_loss_threshold: float = -3.0,
) -> None:
    """
    genome-wide scatter plot.
    - output_short   : short fragment 단독 PNG
    - output_long    : long fragment 단독 PNG
    - output_combined: total / short / long 3-panel 합친 PNG

    bin_start 기준 정렬 후 그림.
    """
    Path(output_short).parent.mkdir(parents=True, exist_ok=True)

    # ── 단독 PNG (short / long) ───────────────────────────────────────
    for output_path, frag_type, zscore_col, call_col, median_col in [
        (output_short, "Short", "zscore_short", "call_short", "median_zscore_short"),
        (output_long,  "Long",  "zscore_long",  "call_long",  "median_zscore_long"),
    ]:
        logger.info("%s genome-wide scatter → %s", frag_type, output_path)
        bin_data       = _load_bin_data(norm_tsv_path, zscore_col)
        cnv_calls      = _load_cnv_calls(cnv_tsv_path, call_col)
        median_zscores = _load_median_zscores(cnv_tsv_path, median_col)

        fig, ax = plt.subplots(figsize=(18, 4))
        _draw_genome_scatter(
            ax, bin_data, cnv_calls, median_zscores,
            zscore_gain_threshold, zscore_loss_threshold,
            frag_type, sample_id,
        )
        fig.tight_layout()
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        logger.info("저장 완료 → %s", output_path)

    # ── 3-panel 합친 PNG (total / short / long) ───────────────────────
    logger.info("combined 3-panel scatter → %s", output_combined)

    panels = [
        ("Total",  "zscore_total", "call_total", "median_zscore_total"),
        ("Short",  "zscore_short", "call_short", "median_zscore_short"),
        ("Long",   "zscore_long",  "call_long",  "median_zscore_long"),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(18, 10), sharex=True)
    fig.subplots_adjust(hspace=0.08)
    fig.suptitle(f"{sample_id} — CNV Z-score  (Total / Short / Long)",
                 fontsize=11, fontweight="bold", y=0.98)

    for ax, (frag_type, zscore_col, call_col, median_col) in zip(axes, panels):
        bin_data       = _load_bin_data(norm_tsv_path, zscore_col)
        cnv_calls      = _load_cnv_calls(cnv_tsv_path, call_col)
        median_zscores = _load_median_zscores(cnv_tsv_path, median_col)
        _draw_genome_scatter(
            ax, bin_data, cnv_calls, median_zscores,
            zscore_gain_threshold, zscore_loss_threshold,
            frag_type, sample_id,
        )
        # 합친 그림에서는 타이틀 대신 Y축 레이블로 구분
        ax.set_title("")
        ax.set_ylabel(f"{frag_type}
Z-score", fontsize=8)
        # X축 레이블은 마지막 패널만
        if ax is not axes[-1]:
            ax.set_xlabel("")

    fig.savefig(output_combined, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("combined scatter 저장 완료 → %s", output_combined)


def plot_cnv_scatter_by_chrom(
    norm_tsv_path:         str,
    cnv_tsv_path:          str,
    output_dir:            str,
    sample_id:             str,
    zscore_gain_threshold: float = 3.0,
    zscore_loss_threshold: float = -3.0,
    chroms:                Optional[list[str]] = None,
    region_cnv_path:       Optional[str] = None,
) -> list[str]:
    """
    chromosome별 total/short/long 3-panel scatter PNG 저장.
    bin_start 기준 정렬, X축 = position (Mbp).
    region_cnv_path 제공 시 GAIN/LOSS region 배경 하이라이트 + bin 색상 표기.
    """
    target = chroms or CHROM_ORDER
    out    = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved: list[str] = []

    # region CNV 로드
    region_map = _load_region_cnv(region_cnv_path, call_col="call_total")

    panels = [
        ("Total", "zscore_total", "call_total", "median_zscore_total"),
        ("Short", "zscore_short", "call_short", "median_zscore_short"),
        ("Long",  "zscore_long",  "call_long",  "median_zscore_long"),
    ]

    # 데이터 사전 로드 (chrom 루프마다 TSV 재읽기 방지)
    all_bin_data:       dict[str, dict] = {}
    all_cnv_calls:      dict[str, dict] = {}
    all_median_zscores: dict[str, dict] = {}
    for frag_type, zscore_col, call_col, median_col in panels:
        all_bin_data[frag_type]       = _load_bin_data(norm_tsv_path, zscore_col)
        all_cnv_calls[frag_type]      = _load_cnv_calls(cnv_tsv_path, call_col)
        all_median_zscores[frag_type] = _load_median_zscores(cnv_tsv_path, median_col)

    for chrom in target:
        regions    = region_map.get(chrom)
        call_total = all_cnv_calls["Total"].get(chrom, "NORMAL")

        fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
        fig.subplots_adjust(hspace=0.08)

        for ax, (frag_type, _, _, _) in zip(axes, panels):
            bins  = all_bin_data[frag_type].get(chrom, [])
            call  = all_cnv_calls[frag_type].get(chrom, "NORMAL")
            med_z = all_median_zscores[frag_type].get(chrom)

            _draw_chrom_scatter(
                ax, chrom, bins, call, med_z,
                zscore_gain_threshold, zscore_loss_threshold,
                frag_type, sample_id,
                regions=regions,
            )
            ax.set_title("")
            ax.set_ylabel(f"{frag_type}\nZ-score", fontsize=8)
            if ax is not axes[-1]:
                ax.set_xlabel("")

        # region 요약 타이틀 문자열
        region_str = ""
        if regions:
            region_str = "  |  " + "  ".join(
                f"{c}:{s//1_000_000:.1f}-{e//1_000_000:.1f}Mb"
                for s, e, c in regions
            )

        fig.suptitle(
            f"{sample_id} — {chrom}  CNV  [Total: {call_total}]{region_str}",
            fontsize=9, fontweight="bold",
            color=CALL_COLOR.get(call_total, "#333333"),
            y=0.99,
        )

        out_png = str(out / f"{sample_id}.{chrom}.cnv.png")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_png)
        logger.info("  ✓ %s (regions=%d) → %s", chrom, len(regions or []), out_png)

    logger.info("chromosome별 CNV scatter 완료: %d PNG → %s", len(saved), output_dir)
    return saved