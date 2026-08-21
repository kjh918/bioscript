"""
cnv_track.py
============
cnv_calls.parquet → 전장 유전체 CNV track.

Row 1: short log2_norm (LOO 정규화) + gain/loss 색상
Row 2: long  log2_norm
Row 3: short/long VAF (있을 경우)

y축: short/long row 공유 고정 (상하위 1% trimmed 기준)
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import logging as _log

_CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

_CALL_COLOR = {
    "gain":    "rgba(220,50,50,0.85)",
    "loss":    "rgba(50,100,220,0.85)",
    "normal":  "rgba(150,150,150,0.4)",
    "unknown": "rgba(200,200,200,0.2)",
}


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


def _trimmed_range(
    vals: np.ndarray,
    q:    float = 0.01,
    pad:  float = 0.15,
) -> tuple[float, float]:
    """상하위 1% 제외 후 y 범위 계산."""
    finite = vals[np.isfinite(vals)]
    if len(finite) == 0:
        return (-1.0, 1.0)
    lo = float(np.quantile(finite, q))
    hi = float(np.quantile(finite, 1.0 - q))
    margin = (hi - lo) * pad
    return (lo - margin, hi + margin)


def plot_cnv_track(
    cnv_path: str,
    out_html: str | None = None,
    title:    str   = "CNV Track",
    height:   int   = 900,
    trim_q:   float = 0.01,
    y_fixed:  float | None = 1.5,
) -> go.Figure:
    df = pd.read_parquet(cnv_path)
    df = df[df["chrom"].isin(_CHROM_ORDER)].copy()

    # mappability_pass=False / gc=0 / count=0 bin → 핵심 값 NaN 으로 마스킹
    # (gc_corrector 에서 처리됐어도 엣지케이스 방어)
    bad_mask = pd.Series(False, index=df.index)
    if "mappability_pass" in df.columns:
        bad_mask |= ~df["mappability_pass"].astype(bool)
    if "gc" in df.columns:
        bad_mask |= ~df["gc"].between(0.01, 0.99, inclusive="both").fillna(True)
    if "short_count" in df.columns:
        bad_mask |= (df["short_count"] + df.get("long_count", 0)) <= 0

    nan_cols = [c for c in df.columns
                if any(k in c for k in ("log2", "zscore", "copy_number",
                                        "baf", "wps", "breadth"))]
    df.loc[bad_mask, nan_cols] = np.nan
    n_masked = int(bad_mask.sum())
    if n_masked:
        import logging as _log
        _log.getLogger(__name__).info(
            "cnv_track: %d bins masked (mappability/gc/count 불량)", n_masked
        )

    df["chrom"] = pd.Categorical(df["chrom"], categories=_CHROM_ORDER, ordered=True)
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)

    x, offsets = _genomic_x(df)

    # y축: copy_number 우선 (정상=2.0), fallback log2_norm → log2_corrected
    use_cn = ("short_copy_number" in df.columns and "long_copy_number" in df.columns)

    if use_cn:
        short_col  = "short_copy_number"
        long_col   = "long_copy_number"
        y_label    = "Copy Number"
        y_normal   = 2.0          # 정상 diploid
        y_gain_ref = 2.585        # trisomy ~2 + FF/2
        y_loss_ref = 1.415        # monosomy ~2 - FF/2
        y_fmt      = ".2f"
        # y축 범위: copy number 기준 ±대칭 (정상=2.0 중심)
        if y_fixed is not None:
            y_min = 2.0 - abs(float(y_fixed))
            y_max = 2.0 + abs(float(y_fixed))
        else:
            combined = np.concatenate([
                df[short_col].dropna().values,
                df[long_col].dropna().values,
            ])
            if len(combined):
                lo = float(np.quantile(combined[np.isfinite(combined)], trim_q))
                hi = float(np.quantile(combined[np.isfinite(combined)], 1.0 - trim_q))
                bound = max(abs(lo - 2.0), abs(hi - 2.0)) * 1.15
            else:
                bound = 1.5
            y_min, y_max = 2.0 - bound, 2.0 + bound
    else:
        short_col  = "short_log2_norm"   if "short_log2_norm"   in df.columns else "log2_corrected_short"
        long_col   = "long_log2_norm"    if "long_log2_norm"    in df.columns else "log2_corrected_long"
        y_label    = "log₂ norm"
        y_normal   = 0.0
        y_gain_ref =  0.585
        y_loss_ref = -0.585
        y_fmt      = ".3f"
        combined = np.concatenate([
            df[short_col].values if short_col in df.columns else np.array([]),
            df[long_col].values  if long_col  in df.columns else np.array([]),
        ])
        if y_fixed is not None:
            y_min, y_max = -abs(float(y_fixed)), abs(float(y_fixed))
        else:
            y_min, y_max = _trimmed_range(combined.astype(float), q=trim_q)
        bound = max(abs(y_min), abs(y_max))
        y_min, y_max = -bound, bound

    short_call = "short_cnv_call"
    long_call  = "long_cnv_call"

    n_rows = 3
    has_baf = any(c in df.columns for c in
                  ["combined_baf_median","short_baf_median","long_baf_median"])
    if has_baf:
        n_rows = 4
    row_heights = ([0.30, 0.30, 0.20, 0.20] if n_rows == 4
                   else [0.42, 0.42, 0.16])[:n_rows]

    cn_label = "Copy Number" if use_cn else "log₂ norm (LOO)"
    subtitles = ([f"Short fragment {cn_label}",
                  f"Long fragment {cn_label}",
                  "BAF (combined)",
                  "VAF"] if n_rows == 4
                 else [f"Short fragment {cn_label}",
                       f"Long fragment {cn_label}",
                       "VAF"])[:n_rows]

    fig = make_subplots(
        rows=n_rows, cols=1, shared_xaxes=True,
        row_heights=row_heights,
        vertical_spacing=0.04,
        subplot_titles=subtitles,
    )

    for row_idx, (y_col, call_col, label) in enumerate([
        (short_col, short_call, "Short"),
        (long_col,  long_call,  "Long"),
    ], start=1):
        if y_col not in df.columns:
            continue
        colors = (df[call_col].map(_CALL_COLOR).fillna(_CALL_COLOR["unknown"])
                  if call_col in df.columns
                  else [_CALL_COLOR["normal"]] * len(df))

        # CN 컬럼이면 hovertemplate에 CN 표기
        cn_or_log = "CN" if use_cn else "log₂"
        fig.add_trace(go.Scatter(
            x=x, y=df[y_col].values,
            mode="markers",
            marker=dict(size=2.5, color=colors, opacity=0.75),
            name=f"{label} {cn_or_log}",
            hovertemplate=(
                "chrom=%{customdata[0]}<br>"
                "pos=%{customdata[1]:,}<br>"
                f"{cn_or_log}=%{{y:{y_fmt}}}<br>"
                "call=%{customdata[2]}"
            ),
            customdata=np.column_stack([
                df["chrom"].astype(str),
                (df["start"] + df["end"]) // 2,
                df[call_col].fillna("unknown") if call_col in df.columns
                else ["unknown"] * len(df),
            ]),
        ), row=row_idx, col=1)

        fig.update_yaxes(range=[y_min, y_max], row=row_idx, col=1,
                         title_text=y_label)
        # 정상 기준선
        fig.add_hline(y=y_normal,   line_dash="dash",
                      line_color="black", line_width=0.8, row=row_idx, col=1)
        # gain/loss 참고선 (trisomy/monosomy 경계)
        fig.add_hline(y=y_gain_ref, line_dash="dot",
                      line_color="rgba(220,50,50,0.4)", line_width=0.6,
                      row=row_idx, col=1)
        fig.add_hline(y=y_loss_ref, line_dash="dot",
                      line_color="rgba(50,100,220,0.4)", line_width=0.6,
                      row=row_idx, col=1)

    # BAF row (combined → short → long 순서로 오버레이)
    baf_row = 3
    if has_baf:
        for baf_col, color, name in [
            ("combined_baf_median", "rgba(40,40,40,0.8)",    "BAF combined"),
            ("short_baf_median",    "rgba(220,80,80,0.6)",   "BAF short"),
            ("long_baf_median",     "rgba(50,120,220,0.6)",  "BAF long"),
        ]:
            if baf_col not in df.columns:
                continue
            fig.add_trace(go.Scatter(
                x=x, y=df[baf_col].values,
                mode="markers",
                marker=dict(size=2.5, color=color, opacity=0.7),
                name=name,
                hovertemplate=(
                    f"{name}<br>chrom=%{{customdata[0]}}<br>"
                    "pos=%{customdata[1]:,}<br>BAF=%{y:.3f}"
                ),
                customdata=np.column_stack([
                    df["chrom"].astype(str),
                    (df["start"] + df["end"]) // 2,
                ]),
            ), row=baf_row, col=1)

        # BAF y축: 0~1 고정, 0.5 기준선
        fig.update_yaxes(range=[0.3, 0.7], row=baf_row, col=1)
        fig.add_hline(y=0.5, line_dash="dash",
                      line_color="black", line_width=0.8,
                      row=baf_row, col=1)
        baf_row = 4   # VAF 는 다음 row

    # VAF row
    vaf_row = baf_row
    if any(c in df.columns for c in ["short_vaf", "long_vaf"]):
        for vaf_col, color, name in [
            ("short_vaf", "rgba(220,120,50,0.6)", "Short VAF"),
            ("long_vaf",  "rgba(50,160,220,0.6)", "Long VAF"),
        ]:
            if vaf_col not in df.columns:
                continue
            fig.add_trace(go.Scatter(
                x=x, y=df[vaf_col].values,
                mode="markers",
                marker=dict(size=1.5, color=color),
                name=name,
            ), row=vaf_row, col=1)
        fig.update_yaxes(range=[0, 1], row=vaf_row, col=1)

    # 염색체 경계선
    shapes, annotations = [], []
    for chrom, offset in offsets.items():
        shapes.append(dict(
            type="line", x0=offset, x1=offset, y0=0, y1=1,
            yref="paper", line=dict(color="lightgray", width=0.7, dash="dot"),
        ))
        sub = df[df["chrom"] == chrom]
        if not sub.empty:
            mid = offset + int(sub["end"].max()) // 2
            annotations.append(dict(
                x=mid, y=1.01, xref="x", yref="paper",
                text=chrom.replace("chr", ""),
                showarrow=False, font=dict(size=9),
            ))

    fig.update_layout(
        title=title, height=height,
        showlegend=True,
        shapes=shapes, annotations=annotations,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11),
    )
    fig.update_xaxes(showticklabels=False, showgrid=False)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.3)")

    if out_html:
        fig.write_html(out_html)
        # PNG / PDF 동시 저장
        base = out_html.rsplit(".", 1)[0]
        try:
            fig.write_image(f"{base}.png", width=1800, height=900, scale=2)
        except Exception as e:
            _log.warning("PNG 저장 실패 (kaleido 미설치 가능): %s", e)
        try:
            fig.write_image(f"{base}.pdf")
        except Exception as e:
            _log.warning("PDF 저장 실패: %s", e)

    return fig