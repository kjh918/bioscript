"""
wps_profile.py
==============
WPS 관련 시각화.

1. plot_fragment_track()        : 단일 frag WPS/endpoint/coverage 3-track
2. plot_fragment_track_compare(): short vs long 6-track 비교
3. plot_wps_genome_track()      : 전장 유전체 bin-level WPS
4. plot_marker_wps()            : marker 단위 WPS 요약
5. plot_wps_peaks()             : peak score 분포

모든 short/long 비교 그래프는 y축을 공유 고정합니다.
outlier 영향 최소화를 위해 상하위 1% trimmed 기준으로 y 범위를 결정합니다.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

_CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ─────────────────────────────────────────────────────────────────────
# 공통 유틸
# ─────────────────────────────────────────────────────────────────────
def _trimmed_yrange(
    *arrays: np.ndarray,
    q:   float = 0.01,
    pad: float = 0.12,
) -> tuple[float, float]:
    """
    여러 배열을 합쳐 상하위 q% 제외 후 y 범위 계산.
    short/long 공유 y축에 사용.
    """
    combined = np.concatenate([a[np.isfinite(a)] for a in arrays if len(a) > 0])
    if len(combined) == 0:
        return (-1.0, 1.0)
    lo = float(np.quantile(combined, q))
    hi = float(np.quantile(combined, 1.0 - q))
    margin = max((hi - lo) * pad, 0.05)
    return (lo - margin, hi + margin)


def _genomic_x(df: pd.DataFrame) -> tuple[np.ndarray, dict[str, int]]:
    chroms  = [c for c in _CHROM_ORDER if c in df["chrom"].values]
    offsets: dict[str, int] = {}
    cursor  = 0
    for chrom in chroms:
        offsets[chrom] = cursor
        cursor += int(df[df["chrom"] == chrom]["end"].max())
    x = np.array([offsets.get(r.chrom, 0) + (r.start + r.end) // 2
                  for r in df.itertuples()])
    return x, offsets


# ─────────────────────────────────────────────────────────────────────
# 1. 단일 region 3-track (WPS / Endpoint / Coverage)
# ─────────────────────────────────────────────────────────────────────
def plot_fragment_track(
    corrected_path: str,
    chrom:          str,
    region_start:   int,
    region_end:     int,
    frag:           str = "long",
    wps_type:       str = "L",
    out_html:       str | None = None,
    title:          str | None = None,
    height:         int = 700,
    trim_q:         float = 0.01,
) -> go.Figure:
    """
    단일 frag × region 에 대해 WPS / Endpoint / Coverage 3-track 시각화.
    각 row y축 독립 고정 (trimmed 기준).
    """
    df = pd.read_parquet(corrected_path)
    df = df[(df["chrom"] == chrom) &
            (df["start"] >= region_start) &
            (df["end"]   <= region_end)].sort_values("start").reset_index(drop=True)

    if df.empty:
        fig = go.Figure()
        fig.update_layout(title=f"{chrom}:{region_start}-{region_end} — 데이터 없음")
        return fig

    x = ((df["start"] + df["end"]) / 2).values

    def _col(name):
        corr = f"{frag}_{name}_corrected"
        return corr if corr in df.columns else f"{frag}_{name}"

    wps_col = _col(f"wps_{wps_type}")
    ep_col  = _col("endpoint_density")
    cov_col = _col("coverage_median")

    _title = title or f"{frag} fragment  {chrom}:{region_start:,}-{region_end:,}"

    fig = make_subplots(
        rows=3, cols=1, shared_xaxes=True,
        vertical_spacing=0.06,
        row_heights=[0.38, 0.31, 0.31],
        subplot_titles=[
            f"WPS ({frag}, {wps_type}-WPS)",
            "Endpoint Density",
            "Coverage (median depth)",
        ],
    )

    track_specs = [
        (1, wps_col,  "WPS",      "rgba(50,120,220,0.85)",  "rgba(50,120,220,0.12)"),
        (2, ep_col,   "Endpoint", "rgba(220,80,80,0.85)",   "rgba(220,80,80,0.12)"),
        (3, cov_col,  "Coverage", "rgba(80,170,80,0.85)",   "rgba(80,170,80,0.12)"),
    ]

    for row_idx, col, name, color, fill in track_specs:
        if col not in df.columns:
            continue
        y = df[col].values.astype(float)
        # 각 row 독립 y 범위
        y_min, y_max = _trimmed_yrange(y, q=trim_q)

        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines",
            line=dict(color=color, width=1.2),
            name=name, fill="tozeroy", fillcolor=fill,
            hovertemplate=f"pos=%{{x:,.0f}}<br>{name}=%{{y:.3f}}",
        ), row=row_idx, col=1)

        fig.update_yaxes(range=[y_min, y_max], row=row_idx, col=1)
        if name in ("WPS", "Endpoint"):
            fig.add_hline(y=0, line_dash="dash",
                          line_color="black", line_width=0.7,
                          row=row_idx, col=1)

    fig.update_xaxes(title_text=f"Genomic position ({chrom})",
                     tickformat=",", row=3, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=2, col=1)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.3)")
    fig.update_layout(
        title=_title, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11), showlegend=True,
    )

    if out_html:
        fig.write_html(out_html)
    return fig


# ─────────────────────────────────────────────────────────────────────
# 2. short vs long 6-track 비교 (row별 y축 공유 고정)
# ─────────────────────────────────────────────────────────────────────
def plot_fragment_track_compare(
    corrected_path: str,
    chrom:          str,
    region_start:   int,
    region_end:     int,
    wps_type:       str = "L",
    out_html:       str | None = None,
    title:          str | None = None,
    height:         int = 900,
    trim_q:         float = 0.01,
) -> go.Figure:
    """
    short vs long 6-track 비교.
    같은 지표(WPS / Endpoint / Coverage)끼리 y축 공유 고정.
    """
    df = pd.read_parquet(corrected_path)
    df = df[(df["chrom"] == chrom) &
            (df["start"] >= region_start) &
            (df["end"]   <= region_end)].sort_values("start").reset_index(drop=True)

    if df.empty:
        fig = go.Figure()
        fig.update_layout(title="데이터 없음")
        return fig

    x = ((df["start"] + df["end"]) / 2).values

    def _col(frag, name):
        corr = f"{frag}_{name}_corrected"
        return corr if corr in df.columns else f"{frag}_{name}"

    _title = title or f"Short vs Long  {chrom}:{region_start:,}-{region_end:,}  ({wps_type}-WPS)"

    # (row, frag, metric_key, color, dash, label)
    tracks = [
        (1, "short", f"wps_{wps_type}", "rgba(50,120,220,0.9)",  "solid", "Short WPS"),
        (2, "long",  f"wps_{wps_type}", "rgba(180,80,220,0.9)",  "solid", "Long WPS"),
        (3, "short", "endpoint_density","rgba(220,80,80,0.9)",   "solid", "Short Endpoint"),
        (4, "long",  "endpoint_density","rgba(220,140,40,0.9)",  "solid", "Long Endpoint"),
        (5, "short", "coverage_median", "rgba(50,170,80,0.9)",   "solid", "Short Coverage"),
        (6, "long",  "coverage_median", "rgba(20,130,130,0.9)",  "solid", "Long Coverage"),
    ]

    subtitles = [t[5] for t in tracks]
    fig = make_subplots(
        rows=6, cols=1, shared_xaxes=True,
        vertical_spacing=0.03,
        row_heights=[0.17] * 6,
        subplot_titles=subtitles,
    )

    # 같은 지표끼리 y 범위 공유 계산
    # rows 1-2: WPS, rows 3-4: Endpoint, rows 5-6: Coverage
    metric_rows = {
        f"wps_{wps_type}":     (1, 2),
        "endpoint_density":    (3, 4),
        "coverage_median":     (5, 6),
    }
    shared_ranges: dict[str, tuple[float, float]] = {}
    for metric, (r1, r2) in metric_rows.items():
        arrs = []
        for frag in ("short", "long"):
            col = _col(frag, metric)
            if col in df.columns:
                arrs.append(df[col].values.astype(float))
        if arrs:
            shared_ranges[metric] = _trimmed_yrange(*arrs, q=trim_q)

    for row, frag, metric, color, dash, label in tracks:
        col = _col(frag, metric)
        if col not in df.columns:
            continue
        y          = df[col].values.astype(float)
        fill_color = color.replace("0.9)", "0.10)")
        y_range    = shared_ranges.get(metric, None)

        fig.add_trace(go.Scatter(
            x=x, y=y, mode="lines",
            line=dict(color=color, width=1.0, dash=dash),
            name=label, fill="tozeroy", fillcolor=fill_color,
            hovertemplate=f"{label}<br>pos=%{{x:,.0f}}<br>val=%{{y:.3f}}",
        ), row=row, col=1)

        if y_range:
            fig.update_yaxes(range=list(y_range), row=row, col=1)
        if "wps" in metric or "endpoint" in metric:
            fig.add_hline(y=0, line_dash="dash",
                          line_color="black", line_width=0.5,
                          row=row, col=1)

    fig.update_xaxes(title_text=f"Genomic position ({chrom})",
                     tickformat=",", row=6, col=1)
    for r in range(1, 6):
        fig.update_xaxes(showticklabels=False, row=r, col=1)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.3)")
    fig.update_layout(
        title=_title, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=10), showlegend=True,
    )

    if out_html:
        fig.write_html(out_html)
    return fig


# ─────────────────────────────────────────────────────────────────────
# 3. 전장 유전체 WPS track (short/long y축 공유)
# ─────────────────────────────────────────────────────────────────────
def plot_wps_genome_track(
    corrected_path: str,
    wps_cols: list[str] | None = None,
    out_html: str | None = None,
    title:    str = "WPS Genome Track",
    height:   int = 1000,
    trim_q:   float = 0.01,
    shared_y_within_group: bool = True,
) -> go.Figure:
    """
    전장 유전체 WPS 4트랙: Short L / Long L / Short S / Long S

    y축 정책
    --------
    shared_y_within_group=True:
      L-WPS 2트랙 공유, S-WPS 2트랙 공유 → L vs S 스케일 차이 유지하면서 short/long 비교
    shared_y_within_group=False:
      각 트랙 독립 y축

    clipping 방지
    -------------
    trimmed range 에 pad=0.2 적용. 데이터 없는 컬럼은 건너뜀.
    """
    if wps_cols is None:
        wps_cols = [
            "short_wps_L_corrected", "long_wps_L_corrected",
            "short_wps_S_corrected", "long_wps_S_corrected",
        ]
    # 실제 존재하는 컬럼만 사용
    df = pd.read_parquet(corrected_path)
    wps_cols = [c for c in wps_cols if c in df.columns]
    if not wps_cols:
        fig = go.Figure()
        fig.update_layout(title=f"{title} — WPS 컬럼 없음")
        return fig

    df = df[df["chrom"].isin(_CHROM_ORDER)].copy()

    # 불량 bin 마스킹
    bad = pd.Series(False, index=df.index)
    if "mappability_pass" in df.columns:
        bad |= ~df["mappability_pass"].astype(bool)
    if "gc" in df.columns:
        bad |= ~df["gc"].between(0.01, 0.99, inclusive="both").fillna(True)
    for col in wps_cols:
        df.loc[bad, col] = np.nan

    df["chrom"] = pd.Categorical(df["chrom"], categories=_CHROM_ORDER, ordered=True)
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    x, offsets = _genomic_x(df)

    # ── y 범위: L 그룹 / S 그룹 분리 ─────────────────────────────
    L_cols = [c for c in wps_cols if "_wps_L" in c]
    S_cols = [c for c in wps_cols if "_wps_S" in c]

    def _yrange(cols):
        arrs = [df[c].values.astype(float) for c in cols if c in df.columns]
        return _trimmed_yrange(*arrs, q=trim_q, pad=0.20) if arrs else (-1.0, 1.0)

    if shared_y_within_group:
        L_range = _yrange(L_cols) if L_cols else None
        S_range = _yrange(S_cols) if S_cols else None
        col_ranges = {}
        for c in wps_cols:
            if c in L_cols:   col_ranges[c] = L_range
            elif c in S_cols: col_ranges[c] = S_range
            else:             col_ranges[c] = _yrange([c])
    else:
        col_ranges = {c: _yrange([c]) for c in wps_cols}

    # ── 서브플롯 구성 ─────────────────────────────────────────────
    subtitle_map = {
        "short_wps_L_corrected": "Short L-WPS  (frag 120-180bp, k=120bp)",
        "long_wps_L_corrected":  "Long  L-WPS  (frag 120-180bp, k=120bp)",
        "short_wps_S_corrected": "Short S-WPS  (frag 35-80bp,   k=16bp)",
        "long_wps_S_corrected":  "Long  S-WPS  (frag 35-80bp,   k=16bp)",
    }
    subtitles = [subtitle_map.get(c, c) for c in wps_cols]

    # L=빨강 계열 / S=파랑 계열, short=진함 / long=연함
    color_map = {
        "short_wps_L_corrected": "rgba(210,60,60,0.80)",
        "long_wps_L_corrected":  "rgba(100,160,230,0.80)",
        "short_wps_S_corrected": "rgba(180,40,40,0.75)",
        "long_wps_S_corrected":  "rgba(40,100,210,0.75)",
    }

    n_rows = len(wps_cols)
    fig = make_subplots(
        rows=n_rows, cols=1, shared_xaxes=True,
        vertical_spacing=0.025,
        subplot_titles=subtitles,
        row_heights=[1.0] * n_rows,
    )

    for i, col in enumerate(wps_cols):
        y_vals = df[col].values.astype(float)
        y_min, y_max = col_ranges[col]
        color = color_map.get(col, "rgba(100,100,100,0.7)")

        fig.add_trace(go.Scatter(
            x=x, y=y_vals,
            mode="lines",
            line=dict(color=color, width=0.6),
            name=subtitle_map.get(col, col),
            hovertemplate="chrom=%{customdata}<br>WPS=%{y:.3f}",
            customdata=df["chrom"].astype(str),
        ), row=i + 1, col=1)

        fig.update_yaxes(
            range=[y_min, y_max],
            row=i + 1, col=1,
            showgrid=True,
            gridcolor="rgba(200,200,200,0.3)",
        )
        fig.add_hline(
            y=0, line_dash="dash",
            line_color="rgba(0,0,0,0.35)", line_width=0.6,
            row=i + 1, col=1,
        )

    # 염색체 경계선 + 레이블
    shapes = [dict(
        type="line", x0=off, x1=off, y0=0, y1=1,
        yref="paper", line=dict(color="lightgray", width=0.6, dash="dot"),
    ) for off in offsets.values()]

    annotations = []
    for chrom, off in offsets.items():
        sub = df[df["chrom"] == chrom]
        if not sub.empty:
            mid = off + int(sub["end"].max()) // 2
            annotations.append(dict(
                x=mid, y=1.005, xref="x", yref="paper",
                text=chrom.replace("chr", ""),
                showarrow=False, font=dict(size=8),
            ))

    fig.update_layout(
        title=title, height=height,
        shapes=shapes, annotations=annotations,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=10),
        showlegend=True,
        legend=dict(orientation="h", y=1.01, x=0, font=dict(size=9)),
    )
    fig.update_xaxes(showticklabels=False, showgrid=False)

    if out_html:
        fig.write_html(out_html)
    return fig


# ─────────────────────────────────────────────────────────────────────
# 4. Marker WPS 요약 (short/long y축 공유)
# ─────────────────────────────────────────────────────────────────────
def plot_marker_wps(
    marker_wps_path: str,
    out_html:  str | None = None,
    title:     str = "Marker WPS Summary",
    height:    int = 700,
    trim_q:    float = 0.01,
) -> go.Figure:
    df = pd.read_parquet(marker_wps_path)
    if df.empty:
        fig = go.Figure()
        fig.update_layout(title=f"{title} — 데이터 없음")
        if out_html:
            fig.write_html(out_html)
        return fig

    df = df[df["chrom"].isin(_CHROM_ORDER)].copy()
    df["chrom"] = pd.Categorical(df["chrom"], categories=_CHROM_ORDER, ordered=True)
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)

    offsets: dict[str, int] = {}
    cursor = 0
    for chrom in _CHROM_ORDER:
        if chrom not in df["chrom"].values:
            continue
        offsets[chrom] = cursor
        cursor += int(df[df["chrom"] == chrom]["end"].max())

    df["x_pos"] = df.apply(
        lambda r: offsets.get(r["chrom"], 0) + (r["start"] + r["end"]) // 2,
        axis=1,
    )

    short_col = "short_wps_L_corrected_mean"
    long_col  = "long_wps_L_corrected_mean"

    # short/long 공유 y 범위
    arrs = [df[c].values.astype(float) for c in (short_col, long_col) if c in df.columns]
    y_min, y_max = _trimmed_yrange(*arrs, q=trim_q) if arrs else (-1.0, 1.0)

    types  = df["marker_type"].unique().tolist() if "marker_type" in df.columns else ["marker"]
    colors = ["rgba(220,80,80,0.7)", "rgba(50,120,220,0.7)",
              "rgba(50,180,100,0.7)", "rgba(180,100,220,0.7)"]
    type_color = {t: colors[i % len(colors)] for i, t in enumerate(types)}

    fig = make_subplots(
        rows=3, cols=1, shared_xaxes=False,
        row_heights=[0.35, 0.35, 0.30],
        vertical_spacing=0.08,
        subplot_titles=[
            "Short WPS_L (marker trimmed median)",
            "Long WPS_L (marker trimmed median)",
            "WPS 분포 (marker_type 별)",
        ],
    )

    for row_idx, col in enumerate([short_col, long_col], start=1):
        if col not in df.columns:
            continue
        for mtype in types:
            sub = df[df["marker_type"] == mtype] if "marker_type" in df.columns else df
            fig.add_trace(go.Scatter(
                x=sub["x_pos"], y=sub[col],
                mode="markers",
                marker=dict(size=4, color=type_color.get(mtype, "gray"), opacity=0.6),
                name=mtype,
                showlegend=(row_idx == 1),
                hovertemplate="marker=%{customdata[0]}<br>chrom=%{customdata[1]}<br>val=%{y:.3f}",
                customdata=sub[["marker_id", "chrom"]].values,
            ), row=row_idx, col=1)

        # 공유 y 범위 적용
        fig.update_yaxes(range=[y_min, y_max], row=row_idx, col=1)
        fig.add_hline(y=0, line_dash="dash",
                      line_color="black", line_width=0.5,
                      row=row_idx, col=1)

    # Box plot
    for mtype in types:
        sub = df[df["marker_type"] == mtype] if "marker_type" in df.columns else df
        if short_col in sub.columns:
            fig.add_trace(go.Box(
                y=sub[short_col], name=f"{mtype} (short)",
                marker_color=type_color.get(mtype, "gray"),
                boxmean=True, showlegend=False,
            ), row=3, col=1)

    fig.update_layout(
        title=title, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11),
    )
    fig.update_xaxes(showticklabels=False, showgrid=False, row=1, col=1)
    fig.update_xaxes(showticklabels=False, showgrid=False, row=2, col=1)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.3)")

    if out_html:
        fig.write_html(out_html)
    return fig


# ─────────────────────────────────────────────────────────────────────
# 5. WPS peak score 분포
# ─────────────────────────────────────────────────────────────────────
def plot_wps_peaks(
    peaks_path: str,
    out_html:   str | None = None,
    title:      str = "WPS Peak Distribution",
    height:     int = 500,
) -> go.Figure:
    peaks = pd.read_parquet(peaks_path)
    if peaks.empty:
        fig = go.Figure()
        fig.update_layout(title=f"{title} — peaks 없음")
        if out_html:
            fig.write_html(out_html)
        return fig

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Peak score 분포 (by wps_type)", "염색체별 peak 수"],
    )

    wps_types = peaks["wps_type"].unique() if "wps_type" in peaks.columns else ["all"]
    colors    = ["rgba(220,80,80,0.7)", "rgba(50,120,220,0.7)",
                 "rgba(50,180,100,0.7)", "rgba(180,100,220,0.7)"]

    # peak score 공유 x 범위
    all_scores = peaks["peak_score"].values.astype(float)
    finite_sc  = all_scores[np.isfinite(all_scores)]
    if len(finite_sc):
        sc_lo = float(np.quantile(finite_sc, 0.01))
        sc_hi = float(np.quantile(finite_sc, 0.99))
    else:
        sc_lo, sc_hi = 0.0, 1.0

    for i, wt in enumerate(wps_types):
        sub = peaks[peaks["wps_type"] == wt] if "wps_type" in peaks.columns else peaks
        fig.add_trace(go.Histogram(
            x=sub["peak_score"],
            name=str(wt),
            marker_color=colors[i % len(colors)],
            opacity=0.7, nbinsx=40,
        ), row=1, col=1)

    fig.update_xaxes(range=[sc_lo, sc_hi], row=1, col=1)

    chrom_counts = peaks.groupby("chrom").size().reset_index(name="count")
    chrom_counts["chrom"] = pd.Categorical(
        chrom_counts["chrom"], categories=_CHROM_ORDER, ordered=True
    )
    chrom_counts = chrom_counts.sort_values("chrom")

    fig.add_trace(go.Bar(
        x=chrom_counts["chrom"].astype(str),
        y=chrom_counts["count"],
        marker_color="rgba(100,150,220,0.7)",
        name="Peak count",
    ), row=1, col=2)

    fig.update_layout(
        title=title, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        barmode="overlay", font=dict(size=11),
    )

    if out_html:
        fig.write_html(out_html)
    return fig
