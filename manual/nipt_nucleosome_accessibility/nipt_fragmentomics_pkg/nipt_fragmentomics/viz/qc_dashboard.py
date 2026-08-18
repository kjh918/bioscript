"""
qc_dashboard.py
===============
샘플 QC 대시보드 (Plotly).

패널 구성
---------
  Row 1 col 1 : GC bias curve (GC vs log2 count, LOWESS fit 오버레이)
  Row 1 col 2 : Mappability 분포 (histogram)
  Row 2 col 1 : Breadth ratio 분포 (short / long)
  Row 2 col 2 : short_count vs long_count scatter (FF proxy)
  Row 3 col 1 : VAF 분포 (short / long)
  Row 3 col 2 : Fetal Fraction 요약 bar
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def plot_qc_dashboard(
    corrected_path:  str,
    ff_json_path:    str | None = None,
    out_html:        str | None = None,
    title:           str = "QC Dashboard",
    sample_id:       str = "",
    height:          int = 1100,
) -> go.Figure:
    """
    bins_corrected.parquet + fetal_fraction.json → QC 대시보드 Figure.
    """
    df = pd.read_parquet(corrected_path)

    # autosome only (QC 노이즈 줄임)
    auto = df[~df["chrom"].isin({"chrX","chrY"}) & df["mappability_pass"]].copy()

    fig = make_subplots(
        rows=3, cols=2,
        vertical_spacing=0.1,
        horizontal_spacing=0.1,
        subplot_titles=[
            "GC bias (short count)",      "Mappability 분포",
            "Breadth ratio 분포",          "Short vs Long count (FF proxy)",
            "VAF 분포 (short / long)",     "Fetal Fraction 추정",
        ],
    )

    # ── Row 1 col 1: GC bias ────────────────────────────────────
    if "gc" in auto.columns and "log2_corrected_short" in auto.columns:
        gc_valid = auto.dropna(subset=["gc","log2_corrected_short"])
        fig.add_trace(go.Scatter(
            x=gc_valid["gc"], y=gc_valid["log2_corrected_short"],
            mode="markers",
            marker=dict(size=3, color="rgba(100,150,220,0.4)"),
            name="bins (short)",
        ), row=1, col=1)

        if "gc_fit_short" in auto.columns:
            fit = auto.dropna(subset=["gc","gc_fit_short"]).sort_values("gc")
            fig.add_trace(go.Scatter(
                x=fit["gc"], y=fit["gc_fit_short"],
                mode="lines",
                line=dict(color="rgba(220,60,60,0.9)", width=2),
                name="LOWESS fit",
            ), row=1, col=1)

    fig.update_xaxes(title_text="GC 비율", row=1, col=1)
    fig.update_yaxes(title_text="log₂ CPM", row=1, col=1)

    # ── Row 1 col 2: Mappability 분포 ───────────────────────────
    if "mappability" in df.columns:
        mp_vals = df["mappability"].dropna()
        fig.add_trace(go.Histogram(
            x=mp_vals, nbinsx=40,
            marker_color="rgba(80,160,100,0.7)",
            name="Mappability",
        ), row=1, col=2)
        fig.add_vline(x=0.75, line_dash="dash",
                      line_color="red", row=1, col=2,
                      annotation_text="threshold")

    fig.update_xaxes(title_text="Mappability", row=1, col=2)

    # ── Row 2 col 1: Breadth ratio 분포 ─────────────────────────
    for col, color, name in [
        ("short_breadth", "rgba(220,100,50,0.6)", "Short breadth"),
        ("long_breadth",  "rgba(50,130,220,0.6)", "Long breadth"),
    ]:
        if col in auto.columns:
            fig.add_trace(go.Histogram(
                x=auto[col], nbinsx=40,
                marker_color=color, opacity=0.7,
                name=name,
            ), row=2, col=1)

    fig.update_xaxes(title_text="Breadth ratio", row=2, col=1)
    fig.update_layout(barmode="overlay")

    # ── Row 2 col 2: Short vs Long count scatter ─────────────────
    if "short_count" in auto.columns and "long_count" in auto.columns:
        sl = auto.dropna(subset=["short_count","long_count"])
        fig.add_trace(go.Scatter(
            x=sl["long_count"], y=sl["short_count"],
            mode="markers",
            marker=dict(size=3, color="rgba(180,80,180,0.4)"),
            name="short vs long",
            hovertemplate="long=%{x}<br>short=%{y}",
        ), row=2, col=2)

        # y=x reference line
        max_val = max(float(sl["long_count"].max()),
                      float(sl["short_count"].max())) * 1.05
        fig.add_trace(go.Scatter(
            x=[0, max_val], y=[0, max_val],
            mode="lines",
            line=dict(color="gray", dash="dash", width=1),
            name="y=x (equal counts)",
            showlegend=False,
        ), row=2, col=2)

    fig.update_xaxes(title_text="Long count", row=2, col=2)
    fig.update_yaxes(title_text="Short count", row=2, col=2)

    # ── Row 3 col 1: VAF 분포 ────────────────────────────────────
    for col, color, name in [
        ("short_vaf", "rgba(220,120,50,0.65)", "Short VAF"),
        ("long_vaf",  "rgba(50,160,220,0.65)", "Long VAF"),
    ]:
        if col in auto.columns:
            vaf_vals = auto[col].dropna()
            vaf_vals = vaf_vals[vaf_vals > 0]
            fig.add_trace(go.Histogram(
                x=vaf_vals, nbinsx=50,
                marker_color=color, opacity=0.7,
                name=name,
            ), row=3, col=1)

    fig.update_xaxes(title_text="VAF", row=3, col=1)

    # ── Row 3 col 2: Fetal Fraction bar ──────────────────────────
    ff_labels, ff_values, ff_colors = [], [], []

    if ff_json_path and Path(ff_json_path).exists():
        with open(ff_json_path) as f:
            ff_data = json.load(f)

        for key, label, color in [
            ("seqff_ff",     "SeqFF",     "rgba(50,160,100,0.8)"),
            ("y_ff",         "Y-chr",     "rgba(220,80,80,0.8)"),
            ("consensus_ff", "Consensus", "rgba(80,80,220,0.9)"),
        ]:
            val = ff_data.get(key)
            if val is not None:
                ff_labels.append(label)
                ff_values.append(val)
                ff_colors.append(color)

    if ff_labels:
        fig.add_trace(go.Bar(
            x=ff_labels, y=ff_values,
            marker_color=ff_colors,
            text=[f"{v:.1%}" for v in ff_values],
            textposition="auto",
            name="FF estimate",
        ), row=3, col=2)
        fig.add_hline(y=0.04, line_dash="dash",
                      line_color="gray", row=3, col=2,
                      annotation_text="FF=4% (QC 하한)")
    else:
        fig.add_annotation(
            text="FF 데이터 없음", x=0.5, y=0.5,
            xref="x6", yref="y6", showarrow=False,
        )

    fig.update_yaxes(title_text="Fetal Fraction", row=3, col=2)

    fig.update_layout(
        title=f"{title}  {('— ' + sample_id) if sample_id else ''}",
        height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11),
        showlegend=True,
    )

    if out_html:
        fig.write_html(out_html)
    return fig
