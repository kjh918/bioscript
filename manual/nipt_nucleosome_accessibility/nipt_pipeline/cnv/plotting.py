"""
cnv/plotting.py

within-sample normalized bin TSV (wss_cnv.normalized_bins.tsv) 기반
CNV Z-score scatter plot

입력 컬럼 (wss_cnv.normalized_bins.tsv):
  chrom, bin_start, bin_end, gc, mappability,
  total_gc, short_gc, long_gc, short_fraction_gc,
  wss_z_total, wss_z_short, wss_z_long, wss_z_short_fraction

기능
----
1. plot_cnv_scatter()          : genome-wide 4-panel PNG
2. plot_cnv_scatter_by_chrom() : chromosome별 4-panel PNG
"""
from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
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

CHROM_ORDER: list[str] = STANDARD_CHROMS

CALL_GAIN   = "GAIN"
CALL_LOSS   = "LOSS"
CALL_NORMAL = "NORMAL"

CALL_COLOR: dict[str, str] = {
    "GAIN":              "#D32F2F",
    "LOSS":              "#1565C0",
    "NORMAL":            "#BDBDBD",
    "INSUFFICIENT_BINS": "#9E9E9E",
}
REGION_COLOR: dict[str, str] = {
    "GAIN": "#B71C1C",
    "LOSS": "#0D47A1",
}
BG_COLORS = ["#F5F5F5", "#EEEEEE"]

# wss bin TSV 컬럼 → panel 정의
# (label, wss_z_col, call_col_in_wss_tsv)
# call은 wss_cnv.tsv에서 별도 로드하므로 여기선 Z-score 컬럼만 정의
PANELS = [
    ("Total",          "wss_z_total"),
    ("Short (count)",  "wss_z_short"),
    ("Long (count)",   "wss_z_long"),
    ("Short fraction", "wss_z_short_fraction"),
]


# ── 데이터 로더 ────────────────────────────────────────────────────────

def _load_wss_bin_data(
    wss_bin_tsv: str,
) -> dict[str, dict[str, list[tuple[int, float]]]]:
    """
    wss_cnv.normalized_bins.tsv → {metric: {chrom: [(bin_start, z), ...]}}
    bin_start 오름차순 정렬 보장.
    """
    rows = read_tsv(wss_bin_tsv)
    metrics = ["wss_z_total", "wss_z_short", "wss_z_long", "wss_z_short_fraction"]

    result: dict[str, dict[str, list[tuple[int, float]]]] = {
        m: {c: [] for c in CHROM_ORDER} for m in metrics
    }

    for r in rows:
        chrom = r.get("chrom", "")
        if chrom not in CHROM_ORDER:
            continue
        try:
            start = int(r["bin_start"])
        except (KeyError, ValueError):
            continue

        for m in metrics:
            val = r.get(m, "NA")
            if val in ("NA", "", None):
                continue
            try:
                result[m][chrom].append((start, float(val)))
            except ValueError:
                continue

    for m in metrics:
        for chrom in result[m]:
            result[m][chrom].sort(key=lambda x: x[0])

    return result


def _load_wss_calls(
    wss_cnv_tsv: str,
) -> dict[str, dict[str, str]]:
    """
    wss_cnv.tsv → {call_col: {chrom: call}}
    segment 단위이므로 chrom별 대표 call = call_ff_corrected 우선
    """
    if not Path(wss_cnv_tsv).exists():
        return {}

    rows = read_tsv(wss_cnv_tsv)
    calls: dict[str, dict[str, str]] = {
        "call_total": {}, "call_short": {}, "call_long": {},
        "call_short_fraction": {}, "call_ff_corrected": {},
    }

    for r in rows:
        chrom = r.get("chrom", "")
        if not chrom:
            continue
        for col in calls:
            v = r.get(col, "NORMAL")
            # 이미 GAIN/LOSS가 있으면 유지, NORMAL이면 덮어쓰지 않음
            if calls[col].get(chrom, "NORMAL") == "NORMAL":
                calls[col][chrom] = v

    return calls


def _load_segment_regions(
    wss_cnv_tsv: str,
) -> dict[str, list]:
    """
    wss_cnv.tsv → {chrom: [SimpleNamespace(segment), ...]}
    GAIN/LOSS segment만 반환.
    """
    if not Path(wss_cnv_tsv).exists():
        return {}

    result: dict[str, list] = {}
    for r in read_tsv(wss_cnv_tsv):
        chrom    = r.get("chrom", "")
        call     = r.get("call_ff_corrected") or r.get("call_total", "NORMAL")
        if call not in ("GAIN", "LOSS"):
            continue
        try:
            seg = SimpleNamespace(
                chrom          = chrom,
                start          = int(r["start"]),
                end            = int(r["end"]),
                size_mb        = float(r.get("size_mb", 0)),
                median_z_total = float(r.get("median_zscore_total", 0)),
                call_total     = call,
            )
            result.setdefault(chrom, []).append(seg)
        except (KeyError, ValueError):
            continue
    return result


def _robust_ylim(
    y_vals: list[float],
    gain_thr: float,
    loss_thr: float,
) -> tuple[float, float]:
    if not y_vals:
        return (loss_thr - 1.0, gain_thr + 1.0)
    arr = np.array(y_vals)
    q25, q75 = np.percentile(arr, [25, 75])
    iqr  = q75 - q25 or 1.0
    med  = float(np.median(arr))
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


# ── genome-wide 단일 panel ─────────────────────────────────────────────
def _draw_genome_panel(
    ax:        plt.Axes,
    bins_dict: dict[str, list[tuple[int, float]]],  # chrom → [(start, z)]
    chrom_calls: dict[str, str],
    gain_thr:  float,
    loss_thr:  float,
    label:     str,
    sample_id: str,
    segments:  dict[str, list],  # chrom → [seg]
) -> None:
    ax.set_facecolor("white")

    x_pos:   list[float] = []
    y_vals:  list[float] = []
    colors:  list[str]   = []
    chrom_centers:    dict[str, float] = {}
    chrom_boundaries: list[float]      = [0.0]
    x_cursor = 0.0

    for idx, chrom in enumerate(CHROM_ORDER):
        bins  = bins_dict.get(chrom, [])
        call  = chrom_calls.get(chrom, "NORMAL")
        color = CALL_COLOR.get(call, CALL_COLOR["NORMAL"])
        n     = len(bins)

        ax.axvspan(x_cursor, x_cursor + max(n, 1),
                   color=BG_COLORS[idx % 2], alpha=0.5, zorder=0)

        if n > 0:
            xs = np.arange(x_cursor, x_cursor + n)
            ys = [v for _, v in bins]

            # segment별 bin 색상
            seg_list = segments.get(chrom, [])
            bin_colors = []
            for bstart, _ in bins:
                matched = False
                for seg in seg_list:
                    if seg.start <= bstart < seg.end:
                        bin_colors.append(REGION_COLOR.get(seg.call_total, color))
                        matched = True
                        break
                if not matched:
                    bin_colors.append(color)

            x_pos.extend(xs.tolist())
            y_vals.extend(ys)
            colors.extend(bin_colors)

            # segment median 수평선
            for seg in seg_list:
                seg_xs = [x for (s, _), x in zip(bins, xs) if seg.start <= s < seg.end]
                if seg_xs:
                    rc = REGION_COLOR.get(seg.call_total, color)
                    ax.hlines(seg.median_z_total, min(seg_xs), max(seg_xs),
                              colors=rc, linewidths=1.5, zorder=3)

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
        fontsize=7,
    )
    ax.set_xlim(0, x_cursor)
    ax.set_ylabel(f"{label}\nWSS Z-score", fontsize=8)

    ymin, ymax = _robust_ylim(y_vals, gain_thr, loss_thr)
    ax.set_ylim(ymin, ymax)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)


# ── chromosome-level 단일 panel ────────────────────────────────────────
def _draw_chrom_panel(
    ax:       plt.Axes,
    chrom:    str,
    bins:     list[tuple[int, float]],
    call:     str,
    gain_thr: float,
    loss_thr: float,
    label:    str,
    segments: list,
) -> None:
    ax.set_facecolor("white")

    if not bins:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="#9E9E9E")
        return

    xs = np.array([b[0] for b in bins], dtype=np.float64)
    ys = np.array([b[1] for b in bins], dtype=np.float64)

    # segment 띠 + median 수평선 + 라벨
    for seg in segments:
        if seg.call_total == CALL_NORMAL:
            continue
        rc   = REGION_COLOR.get(seg.call_total, "#E0E0E0")
        s_mb = seg.start / 1e6
        e_mb = seg.end   / 1e6
        med  = seg.median_z_total

        ax.axvspan(s_mb, e_mb, color=rc, alpha=0.18, zorder=0)
        for xb in (s_mb, e_mb):
            ax.axvline(xb, color=rc, linewidth=1.0, linestyle=":", alpha=0.8, zorder=2)
        ax.hlines(med, s_mb, e_mb, colors=rc, linewidth=2.0, alpha=0.85, zorder=4)
        mid = (s_mb + e_mb) / 2
        ax.text(mid, med + 0.3,
                f"{seg.call_total}\n(Z={med:.1f}, {seg.size_mb:.1f}Mb)",
                ha="center", va="bottom", fontsize=7,
                color=rc, fontweight="bold", zorder=5)

    # bin scatter
    bin_colors: list[str] = []
    for bstart, _ in bins:
        matched = False
        for seg in segments:
            if seg.call_total != CALL_NORMAL and seg.start <= bstart < seg.end:
                bin_colors.append(REGION_COLOR.get(seg.call_total, CALL_COLOR["NORMAL"]))
                matched = True
                break
        if not matched:
            bin_colors.append(CALL_COLOR.get(call, CALL_COLOR["NORMAL"]))

    ax.scatter(xs / 1e6, ys, c=bin_colors, s=5, alpha=0.75, linewidths=0, zorder=3)

    ax.axhline(gain_thr, color="#D32F2F", linestyle="--", linewidth=0.9, alpha=0.8, zorder=5)
    ax.axhline(loss_thr, color="#1565C0", linestyle="--", linewidth=0.9, alpha=0.8, zorder=5)
    ax.axhline(0, color="#757575", linestyle="-", linewidth=0.6, alpha=0.5, zorder=1)

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.0f}"))
    if len(xs):
        ax.set_xlim(xs[0] / 1e6, xs[-1] / 1e6)

    ymin, ymax = _robust_ylim(list(ys), gain_thr, loss_thr)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel(f"{label}\nWSS Z-score", fontsize=8)

    # 범례
    handles = _make_legend(gain_thr, loss_thr)
    n_loss = sum(1 for s in segments if s.call_total == CALL_LOSS)
    n_gain = sum(1 for s in segments if s.call_total == CALL_GAIN)
    if n_loss:
        handles.append(mpatches.Patch(color=REGION_COLOR["LOSS"], alpha=0.5, label="DEL segment"))
    if n_gain:
        handles.append(mpatches.Patch(color=REGION_COLOR["GAIN"], alpha=0.5, label="GAIN segment"))
    title = f"Segments: {len(segments)} | DEL={n_loss} GAIN={n_gain}" if segments else None
    ax.legend(handles=handles, fontsize=7, loc="upper right",
              framealpha=0.8, title=title, title_fontsize=7)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def plot_cnv_scatter(
    norm_tsv_path:         str,   # wss_cnv.normalized_bins.tsv
    cnv_tsv_path:          str,   # wss_cnv.tsv
    output_short:          str,   # 미사용 (하위 호환)
    output_long:           str,   # 미사용 (하위 호환)
    output_combined:       str,   # 4-panel genome-wide PNG
    sample_id:             str,
    zscore_gain_threshold: float = 3.0,
    zscore_loss_threshold: float = -3.0,
) -> None:
    """
    genome-wide 4-panel scatter (Total / Short / Long / Short fraction).
    output_combined에 저장.
    """
    Path(output_combined).parent.mkdir(parents=True, exist_ok=True)

    if not Path(norm_tsv_path).exists():
        logger.warning("WSS bin TSV 없음 → scatter plot 스킵: %s", norm_tsv_path)
        return

    logger.info("genome-wide scatter 로딩: %s", norm_tsv_path)
    all_data   = _load_wss_bin_data(norm_tsv_path)
    wss_calls  = _load_wss_calls(cnv_tsv_path)
    segments   = _load_segment_regions(cnv_tsv_path)

    fig, axes = plt.subplots(4, 1, figsize=(20, 13), sharex=True)
    fig.subplots_adjust(hspace=0.06)
    fig.suptitle(
        f"{sample_id} — Within-sample CNV Z-score"
        f"  (Total / Short / Long / Short fraction)",
        fontsize=11, fontweight="bold", y=0.995,
    )

    call_map = {
        "Total":          wss_calls.get("call_ff_corrected", {}),
        "Short (count)":  wss_calls.get("call_short", {}),
        "Long (count)":   wss_calls.get("call_long", {}),
        "Short fraction": wss_calls.get("call_short_fraction", {}),
    }

    for ax, (label, z_col) in zip(axes, PANELS):
        _draw_genome_panel(
            ax,
            bins_dict   = all_data.get(z_col, {}),
            chrom_calls = call_map.get(label, {}),
            gain_thr    = zscore_gain_threshold,
            loss_thr    = zscore_loss_threshold,
            label       = label,
            sample_id   = sample_id,
            segments    = segments,
        )
        if ax is not axes[-1]:
            ax.set_xlabel("")
        else:
            ax.set_xlabel("Chromosome", fontsize=9)

    fig.savefig(output_combined, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("genome-wide scatter 저장: %s", output_combined)

    # output_short / output_long 도 combined PNG 복사 (하위 호환)
    import shutil
    for p in (output_short, output_long):
        if p and p != output_combined:
            try:
                shutil.copy2(output_combined, p)
            except Exception:
                pass


def plot_cnv_scatter_by_chrom(
    norm_tsv_path:         str,   # wss_cnv.normalized_bins.tsv
    cnv_tsv_path:          str,   # wss_cnv.tsv
    output_dir:            str,
    sample_id:             str,
    zscore_gain_threshold: float = 3.0,
    zscore_loss_threshold: float = -3.0,
    chroms:                Optional[list[str]] = None,
    region_cnv_path:       Optional[str] = None,  # 하위 호환 (미사용)
) -> list[str]:
    """
    chromosome별 4-panel scatter PNG.
    각 panel: Total / Short / Long / Short fraction WSS Z-score.
    """
    target = chroms or CHROM_ORDER
    out    = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved: list[str] = []

    if not Path(norm_tsv_path).exists():
        logger.warning("WSS bin TSV 없음 → chrom scatter 스킵: %s", norm_tsv_path)
        return saved

    logger.info("chrom별 scatter 로딩: %s", norm_tsv_path)
    all_data  = _load_wss_bin_data(norm_tsv_path)
    wss_calls = _load_wss_calls(cnv_tsv_path)
    segments  = _load_segment_regions(cnv_tsv_path)

    call_map = {
        "Total":          wss_calls.get("call_ff_corrected", {}),
        "Short (count)":  wss_calls.get("call_short", {}),
        "Long (count)":   wss_calls.get("call_long", {}),
        "Short fraction": wss_calls.get("call_short_fraction", {}),
    }

    for chrom in target:
        call_total = call_map["Total"].get(chrom, "NORMAL")
        seg_list   = segments.get(chrom, [])

        fig, axes = plt.subplots(4, 1, figsize=(14, 12), sharex=True)
        fig.subplots_adjust(hspace=0.06)

        for ax, (label, z_col) in zip(axes, PANELS):
            bins = all_data.get(z_col, {}).get(chrom, [])
            call = call_map.get(label, {}).get(chrom, "NORMAL")

            _draw_chrom_panel(
                ax, chrom, bins, call,
                zscore_gain_threshold, zscore_loss_threshold,
                label, seg_list,
            )
            if ax is not axes[-1]:
                ax.set_xlabel("")
            else:
                ax.set_xlabel(f"{chrom} position (Mbp)", fontsize=9)

        # 타이틀
        region_str = ""
        if seg_list:
            region_str = "  |  " + "  ".join(
                f"{s.call_total}:{s.start//1_000_000:.1f}-{s.end//1_000_000:.1f}Mb"
                for s in seg_list
            )

        fig.suptitle(
            f"{sample_id} — {chrom}  WSS CNV  [FF-corr: {call_total}]{region_str}",
            fontsize=9, fontweight="bold",
            color=CALL_COLOR.get(call_total, "#333333"),
            y=0.995,
        )

        out_png = str(out / f"{sample_id}.{chrom}.cnv.png")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_png)
        logger.info("  ✓ %s [%s] segs=%d → %s",
                    chrom, call_total, len(seg_list), out_png)

    logger.info("chrom별 scatter 완료: %d PNG → %s", len(saved), output_dir)
    return saved