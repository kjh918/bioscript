"""
Plotly figure builders.
Pure data → Figure; no Dash dependency.
"""

from __future__ import annotations
from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..core.nipt_markers import NiptSyndrome, CALL_COLORS


# ── colour palette ──────────────────────────────────────────────────────────
_CN_COLOR   = "#4299E1"   # blue  – CN dots
_BAF_COLOR  = "#805AD5"   # purple – BAF dots
_CN_MED     = "#2B6CB0"
_BAF_MED    = "#553C9A"
_GRID       = "#EDF2F7"
_AXIS_COLOR = "#718096"
_FONT       = "Inter, Arial, sans-serif"

# CN reference lines
_CN_LINES = {1: "#FC8181", 2: "#68D391", 3: "#F6AD55", 4: "#FC8181"}

FEAT_ALPHA = {
    "TargetChromosome":      0.05,
    "PartialChromosome":     0.07,
    "PrimaryTargetRegion":   0.12,
    "CoreRegion":            0.20,
    "CoreGene":              0.32,
}

def _marker_overlays(
    fig: go.Figure,
    syndromes: list[NiptSyndrome],
    chrom: str,
    start_bp: int,
    end_bp: int,
    cn_row: int = 1,
    baf_row: int = 2,
) -> None:
    """Add shaded region + label for each syndrome feature on this chrom."""
    chrom_key = chrom.replace("chr", "")
    seen_labels: set[str] = set()

    def _hex_to_rgba(hex_color: str, alpha: float) -> str:
        h = hex_color.lstrip("#")
        r, g, b = int(h[0:2],16), int(h[2:4],16), int(h[4:6],16)
        return f"rgba({r},{g},{b},{alpha})"

    for syn in syndromes:
        for feat in syn.features:
            fc = feat.chromosome.replace("chr", "")
            if fc != chrom_key:
                continue
            s = max(feat.start, start_bp)
            e = min(feat.end,   end_bp)
            if e <= s:
                continue

            color      = syn.call_color
            fill_alpha = FEAT_ALPHA.get(feat.feature_type, 0.12)
            label_key  = f"{syn.nipt_id}:{feat.feature_name}"

            for row in (cn_row, baf_row):
                fig.add_vrect(
                    x0=s / 1e6, x1=e / 1e6,
                    fillcolor=_hex_to_rgba(color, fill_alpha),
                    line_color=_hex_to_rgba(color, 0.7),
                    line_width = 0 if feat.feature_type == "TargetChromosome" else (
                                0.6 if feat.feature_type == "PartialChromosome" else
                                1.0 if feat.feature_type == "PrimaryTargetRegion" else
                                1.2),
                    line_dash = "dot" if feat.feature_type in ("CoreGene",) else "solid",
                    opacity=1,
                    row=row, col=1,
                )

            # label once per feature (top of CN panel)
            if label_key not in seen_labels:
                seen_labels.add(label_key)
                short = feat.feature_name if len(feat.feature_name) <= 12 else feat.feature_name[:11] + "…"
                fig.add_annotation(
                    x=(s + e) / 2 / 1e6,
                    y=4.6,
                    text=f"<b>{short}</b>",
                    showarrow=False,
                    font={"size": 8, "color": color, "family": _FONT},
                    bgcolor="white",
                    bordercolor=color,
                    borderwidth=0.8,
                    borderpad=2,
                    opacity=0.92,
                    row=cn_row, col=1,
                )


def region_fig(
    df: pd.DataFrame,
    chrom: str,
    start_bp: int,
    end_bp: int,
    syndromes: Optional[list[NiptSyndrome]] = None,
) -> go.Figure:
    """
    CN + BAF scatter 서브플롯.

    Parameters
    ----------
    df        : DataFrame (pos, cn[, baf][, chrom])
    chrom     : chromosome without "chr"
    start_bp  : region start bp
    end_bp    : region end bp
    syndromes : NiptSyndrome 리스트 (오버레이용, 없으면 생략)
    """
    chrom = str(chrom).replace("chr", "")

    work = df.copy()
    if "chrom" in work.columns:
        work["chrom"] = work["chrom"].astype(str).str.replace("chr", "", regex=False)
        work = work[work["chrom"] == chrom]

    sub = work[(work["pos"] >= start_bp) & (work["pos"] <= end_bp)].copy()

    span_mb   = max((end_bp - start_bp) / 1e6, 0.001)
    marker_sz = 5 if span_mb <= 2 else (4 if span_mb <= 30 else 3)

    baf_col = next((c for c in ("baf", "vaf") if c in sub.columns), None)
    n_rows  = 2 if baf_col else 1
    titles  = [f"Copy Number  chr{chrom}", "B-Allele Frequency (BAF)"] if baf_col else [f"Copy Number  chr{chrom}"]

    fig = make_subplots(
        rows=n_rows, cols=1,
        shared_xaxes=True,
        row_heights=[0.55, 0.45] if baf_col else [1],
        vertical_spacing=0.06,
        subplot_titles=titles,
    )

    # ── CN reference lines ──────────────────────────────────────────────
    for cn_val, c in _CN_LINES.items():
        fig.add_hline(
            y=cn_val,
            line={"color": c, "width": 1 if cn_val == 2 else 0.6, "dash": "dot"},
            row=1, col=1,
        )

    # ── CN scatter ──────────────────────────────────────────────────────
    if not sub.empty:
        fig.add_trace(
            go.Scatter(
                x=sub["pos"] / 1e6,
                y=sub["cn"],
                mode="markers",
                name="CN",
                marker={
                    "size":    marker_sz,
                    "opacity": 0.75,
                    "color":   sub["cn"],
                    "colorscale": [
                        [0.0,  "#90CDF4"],
                        [0.33, "#4299E1"],
                        [0.5,  "#68D391"],
                        [0.66, "#F6AD55"],
                        [1.0,  "#FC8181"],
                    ],
                    "cmin": 0, "cmax": 4,
                    "showscale": False,
                },
                hovertemplate="%{x:.3f} Mb<br>CN = %{y:.3f}<extra></extra>",
            ),
            row=1, col=1,
        )
        cn_med = float(sub["cn"].median())
        fig.add_hline(
            y=cn_med,
            line={"color": _CN_MED, "width": 1.8, "dash": "solid"},
            row=1, col=1,
        )
        fig.add_annotation(
            x=end_bp / 1e6, y=cn_med,
            text=f"  median {cn_med:.2f}",
            xanchor="left", yanchor="middle", showarrow=False,
            font={"size": 9, "color": _CN_MED, "family": _FONT},
            row=1, col=1,
        )

    # ── BAF scatter ─────────────────────────────────────────────────────
    if baf_col and not sub.empty:
        fig.add_trace(
            go.Scatter(
                x=sub["pos"] / 1e6,
                y=sub[baf_col],
                mode="markers",
                name="BAF",
                marker={
                    "size":    marker_sz,
                    "opacity": 0.60,
                    "color":   _BAF_COLOR,
                },
                hovertemplate="%{x:.3f} Mb<br>BAF = %{y:.3f}<extra></extra>",
            ),
            row=2, col=1,
        )
        fig.add_hline(y=0.5, line={"color": "#A0AEC0", "width": 0.8, "dash": "dot"}, row=2, col=1)

    if sub.empty:
        fig.add_annotation(
            text="선택 구간에 데이터가 없습니다",
            x=0.5, y=0.5, xref="paper", yref="paper",
            showarrow=False,
            font={"size": 12, "color": "#A0AEC0", "family": _FONT},
        )

    # ── Syndrome overlays ────────────────────────────────────────────────
    if syndromes:
        _marker_overlays(
            fig, syndromes, chrom, start_bp, end_bp,
            cn_row=1, baf_row=2 if baf_col else 1,
        )

    # ── Axes styling ─────────────────────────────────────────────────────
    axis_common = {
        "showgrid":  True,
        "gridcolor": _GRID,
        "gridwidth": 1,
        "zeroline":  False,
        "tickfont":  {"size": 9, "color": _AXIS_COLOR, "family": _FONT},
        "linecolor": "#CBD5E0",
        "linewidth": 1,
        "showline":  True,
        "ticks":     "outside",
        "ticklen":   3,
    }
    _tf = {"size": 10, "color": _AXIS_COLOR, "family": _FONT}
    # axis_common에 title 없음 – 아래서 개별 지정
    fig.update_xaxes(**axis_common, range=[start_bp / 1e6, end_bp / 1e6])
    fig.update_yaxes(**axis_common,
                     title={"text": "Copy Number", "font": _tf},
                     range=[-0.1, 5.2], row=1, col=1)
    if baf_col:
        fig.update_yaxes(**axis_common,
                         title={"text": "BAF", "font": _tf},
                         range=[-0.05, 1.05], row=2, col=1)
    fig.update_xaxes(title={"text": "Genomic position (Mb)", "font": _tf},
                     row=n_rows, col=1)

    # subplot title styling
    for ann in fig.layout.annotations:
        ann.font = {"size": 10, "color": "#4A5568", "family": _FONT}

    fig.update_layout(
        height=400 if baf_col else 240,
        margin={"l": 54, "r": 80, "t": 36, "b": 36},
        showlegend=False,
        paper_bgcolor="white",
        plot_bgcolor="white",
        font={"family": _FONT, "size": 10, "color": "#2D3748"},
        hovermode="closest",
    )
    return fig


def syndrome_summary_fig(
    syndromes: dict[str, NiptSyndrome],
) -> go.Figure:
    """
    전체 syndrome 판정 결과를 가로 bar 차트로 시각화.
    각 bar = syndrome, 색상 = call 결과.
    """
    if not syndromes:
        fig = go.Figure()
        fig.add_annotation(text="마커 데이터 없음", x=0.5, y=0.5,
                           xref="paper", yref="paper", showarrow=False)
        return fig

    rows = sorted(syndromes.values(), key=lambda s: (s.group, s.syndrome))
    labels  = [s.syndrome for s in rows]
    calls   = [s.call for s in rows]
    colors  = [s.call_color for s in rows]
    cn_vals = [s.cn_value if s.cn_value is not None else 2.0 for s in rows]
    groups  = [s.group for s in rows]

    # call 텍스트 표시
    call_labels = [
        f"<b>{c}</b>" if c != "NORMAL" else c
        for c in calls
    ]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=cn_vals,
        y=labels,
        orientation="h",
        marker={
            "color":       colors,
            "opacity":     0.85,
            "line":        {"color": "white", "width": 0.5},
        },
        text=call_labels,
        textposition="inside",
        textfont={"size": 10, "color": "white", "family": _FONT},
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Call: %{text}<br>"
            "CN: %{x:.2f}<extra></extra>"
        ),
        customdata=groups,
    ))

    # CN=2 reference
    fig.add_vline(x=2.0, line={"color": "#68D391", "width": 1.5, "dash": "dash"})

    _tf = {"size": 10, "color": _AXIS_COLOR, "family": _FONT}
    axis_common = {
        "showgrid":  True, "gridcolor": _GRID, "gridwidth": 1,
        "tickfont":  {"size": 9, "color": _AXIS_COLOR, "family": _FONT},
        "linecolor": "#CBD5E0", "linewidth": 1, "showline": True,
    }
    fig.update_xaxes(**axis_common,
                     title={"text": "Copy Number Signal", "font": _tf}, range=[0, 4.5])
    fig.update_yaxes(**axis_common, autorange="reversed")

    fig.update_layout(
        height=max(260, len(rows) * 28 + 60),
        margin={"l": 10, "r": 20, "t": 30, "b": 40},
        paper_bgcolor="white",
        plot_bgcolor="white",
        font={"family": _FONT, "size": 10},
        showlegend=False,
        title={
            "text": "NIPT Syndrome Call Summary",
            "font": {"size": 12, "color": "#2D3748", "family": _FONT},
            "x": 0.01,
        },
    )
    return fig
