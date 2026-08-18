"""
pipeline.py — nipt_fragmentomics 파이프라인 오케스트레이터.

Step 1  bin_extractor    : BAM → bins_raw.parquet       (100kb, CNV/FF용)
Step 2  gc_corrector     : bins_raw → bins_corrected.parquet
Step 3  marker_extractor : BAM → marker_stats.parquet   (marker BED 기준)
Step 4  fetal_fraction   : bins_corrected → fetal_fraction.json
Step 5  cnv_caller       : bins_corrected → cnv_calls.parquet
Step 6  baf_calculator   : BAM + VCF → bins_baf.parquet → cnv_baf.parquet
"""
from __future__ import annotations

import json, logging, os, time
from datetime import datetime
from typing import Optional

from nipt_fragmentomics.core.constants import DEFAULT_BIN_SIZE, FNAME
from nipt_fragmentomics.steps import (
    bin_extractor, gc_corrector, marker_extractor,
    fetal_fraction, cnv_caller, baf_calculator,
)
from nipt_fragmentomics.viz import cnv_track, qc_dashboard

log = logging.getLogger(__name__)


class Paths:
    def __init__(self, out_dir: str):
        self.out_dir        = out_dir
        self.bins_raw       = os.path.join(out_dir, FNAME["bins_raw"])
        self.bins_corrected = os.path.join(out_dir, FNAME["bins_corrected"])
        self.marker_stats   = os.path.join(out_dir, FNAME["marker_stats"])
        self.fetal_fraction = os.path.join(out_dir, FNAME["fetal_fraction"])
        self.cnv_calls      = os.path.join(out_dir, FNAME["cnv_calls"])
        self.bins_baf       = os.path.join(out_dir, FNAME["bins_baf"])
        self.cnv_baf        = os.path.join(out_dir, FNAME["cnv_baf"])
        self.manifest       = os.path.join(out_dir, FNAME["manifest"])
        self.viz_dir        = os.path.join(out_dir, "viz")


def _e(path: str) -> bool:
    return os.path.exists(path)


def _save_manifest(path: str, params: dict, timings: dict) -> None:
    with open(path, "w") as f:
        json.dump({
            "created_at":  datetime.now().isoformat(),
            "params":      {k: str(v) for k, v in params.items()},
            "timings_sec": timings,
            "total_sec":   round(sum(timings.values()), 2),
        }, f, indent=2)


def run(
    bam_path:    str,
    out_dir:     str,
    bin_bed:     Optional[str] = None,
    marker_bed:  Optional[str] = None,
    fasta_path:  Optional[str] = None,
    bw_path:     Optional[str] = None,
    vcf_path:    Optional[str] = None,
    bin_size:    int   = DEFAULT_BIN_SIZE,
    min_mapq:    int   = 20,
    min_baseq:   int   = 20,
    min_mappability: float = 0.75,
    zscore_gain: float = 3.0,
    zscore_loss: float = -3.0,
    baf_af_min:    float = 0.2,
    baf_af_max:    float = 0.8,
    baf_min_depth: int   = 5,
    n_jobs:    int  = 4,
    resume:    bool = False,
    make_viz:  bool = True,
    sample_id: str  = "",
) -> Paths:
    os.makedirs(out_dir, exist_ok=True)
    p = Paths(out_dir)
    os.makedirs(p.viz_dir, exist_ok=True)

    timings: dict[str, float] = {}
    params = dict(
        bam_path=bam_path, bin_bed=bin_bed, marker_bed=marker_bed,
        fasta_path=fasta_path, bw_path=bw_path, vcf_path=vcf_path,
        bin_size=bin_size, min_mapq=min_mapq, min_baseq=min_baseq,
        min_mappability=min_mappability,
        zscore_gain=zscore_gain, zscore_loss=zscore_loss,
        baf_af_min=baf_af_min, baf_af_max=baf_af_max,
        baf_min_depth=baf_min_depth, n_jobs=n_jobs,
    )

    log.info("=" * 60)
    log.info("nipt_fragmentomics  sample=%s", sample_id or "—")
    log.info("CNV bin  : %s", bin_bed or f"auto {bin_size:,} bp")
    log.info("marker   : %s", marker_bed or "없음")
    log.info("BAF      : %s", "VCF 있음" if vcf_path else "생략")
    log.info("=" * 60)

    # ── Step 1: CNV용 bin 스캔 ───────────────────────────────────
    if resume and _e(p.bins_raw):
        log.info("[Step 1] 건너뜀 (resume)")
    else:
        t0 = time.time()
        bin_extractor.run(
            bam_path=bam_path, out_path=p.bins_raw,
            bed_path=bin_bed, fasta_path=fasta_path,
            bw_path=bw_path, bin_size=bin_size,
            min_mapq=min_mapq, n_jobs=n_jobs,
        )
        timings["step1_bin_extractor"] = round(time.time() - t0, 2)
        log.info("[Step 1] 완료  %.1f s", timings["step1_bin_extractor"])

    # ── Step 2: GC 보정 ──────────────────────────────────────────
    if resume and _e(p.bins_corrected):
        log.info("[Step 2] 건너뜀 (resume)")
    else:
        t0 = time.time()
        gc_corrector.run(
            raw_path=p.bins_raw, out_path=p.bins_corrected,
            min_mappability=min_mappability,
        )
        timings["step2_gc_corrector"] = round(time.time() - t0, 2)
        log.info("[Step 2] 완료  %.1f s", timings["step2_gc_corrector"])

    # ── Step 3: marker 통계 (marker BED 있을 때만) ───────────────
    if marker_bed and os.path.exists(marker_bed):
        if resume and _e(p.marker_stats):
            log.info("[Step 3] 건너뜀 (resume)")
        else:
            t0 = time.time()
            marker_extractor.run(
                bam_path=bam_path,
                marker_bed=marker_bed,
                out_path=p.marker_stats,
                fasta_path=fasta_path,
                bw_path=bw_path,
                min_mapq=min_mapq,
                n_jobs=n_jobs,
            )
            timings["step3_marker_extractor"] = round(time.time() - t0, 2)
            log.info("[Step 3] 완료  %.1f s", timings["step3_marker_extractor"])
    else:
        log.info("[Step 3] marker BED 없음 — 생략")

    # ── Step 4: Fetal Fraction ───────────────────────────────────
    if resume and _e(p.fetal_fraction):
        log.info("[Step 4] 건너뜀 (resume)")
    else:
        t0 = time.time()
        ff_result = fetal_fraction.run(
            corrected_path=p.bins_corrected,
            out_path=p.fetal_fraction,
        )
        timings["step4_fetal_fraction"] = round(time.time() - t0, 2)
        log.info("[Step 4] 완료  consensus FF=%s  %.1f s",
                 ff_result.get("consensus_ff"), timings["step4_fetal_fraction"])

    # ── Step 5: CNV ───────────────────────────────────────────────
    if resume and _e(p.cnv_calls):
        log.info("[Step 5] 건너뜀 (resume)")
    else:
        t0 = time.time()
        cnv_caller.run(
            corrected_path=p.bins_corrected,
            out_path=p.cnv_calls,
            ff_json_path=p.fetal_fraction,
            zscore_gain=zscore_gain, zscore_loss=zscore_loss,
        )
        timings["step5_cnv_caller"] = round(time.time() - t0, 2)
        log.info("[Step 5] 완료  %.1f s", timings["step5_cnv_caller"])

    # ── Step 6: BAF ───────────────────────────────────────────────
    if vcf_path and os.path.exists(vcf_path):
        if resume and _e(p.cnv_baf):
            log.info("[Step 6] 건너뜀 (resume)")
        else:
            t0 = time.time()
            baf_calculator.run(
                bam_path=bam_path, vcf_path=vcf_path,
                bin_path=p.bins_raw, out_path=p.bins_baf,
                min_mapq=min_mapq, min_baseq=min_baseq,
                min_depth=baf_min_depth,
                af_min=baf_af_min, af_max=baf_af_max,
                n_jobs=n_jobs,
            )
            if _e(p.bins_baf) and _e(p.cnv_calls):
                baf_calculator.merge_into_cnv(
                    cnv_path=p.cnv_calls,
                    baf_path=p.bins_baf,
                    out_path=p.cnv_baf,
                )
            timings["step6_baf_calculator"] = round(time.time() - t0, 2)
            log.info("[Step 6] 완료  %.1f s", timings["step6_baf_calculator"])
    elif vcf_path:
        log.warning("[Step 6] VCF 없음: %s", vcf_path)
    else:
        log.info("[Step 6] --vcf 미지정 — 생략")

    # ── Viz ───────────────────────────────────────────────────────
    if make_viz:
        t0 = time.time()
        _make_viz(p, sample_id=sample_id)
        timings["viz"] = round(time.time() - t0, 2)
        log.info("[Viz] 완료  %.1f s", timings["viz"])

    _save_manifest(p.manifest, params, timings)
    log.info("완료  total=%.1f s  →  %s", sum(timings.values()), out_dir)
    return p


def _plot_marker_stats(
    marker_stats_path: str,
    out_html: str,
    sample_id: str = "",
) -> None:
    """
    marker_stats.parquet 에서 short/long WPS 및 count 비교 시각화.
    """
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    df = pd.read_parquet(marker_stats_path)
    if df.empty:
        return

    types   = df["marker_type"].unique().tolist() if "marker_type" in df.columns else ["marker"]
    colors  = ["rgba(220,80,80,0.7)", "rgba(50,120,220,0.7)",
               "rgba(50,180,100,0.7)", "rgba(180,80,220,0.7)"]

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            "Short vs Long count (marker별)",
            "Short ratio (fetal fraction proxy)",
            "L-WPS: Short vs Long",
            "S-WPS: Short vs Long",
        ],
        vertical_spacing=0.12, horizontal_spacing=0.1,
    )

    for i, mtype in enumerate(types):
        sub   = df[df["marker_type"] == mtype] if "marker_type" in df.columns else df
        color = colors[i % len(colors)]
        name  = str(mtype)

        # Row1 Col1: short/long count scatter
        fig.add_trace(go.Scatter(
            x=sub["long_count"], y=sub["short_count"],
            mode="markers", marker=dict(size=3, color=color, opacity=0.5),
            name=name, showlegend=True,
            hovertemplate="marker=%{customdata}<br>long=%{x}<br>short=%{y}",
            customdata=sub["marker_id"].astype(str),
        ), row=1, col=1)

        # Row1 Col2: short_ratio box
        fig.add_trace(go.Box(
            y=sub["short_ratio"], name=name,
            marker_color=color, boxmean=True, showlegend=False,
        ), row=1, col=2)

        # Row2 Col1: L-WPS short vs long scatter
        if "short_wps_L" in sub.columns and "long_wps_L" in sub.columns:
            fig.add_trace(go.Scatter(
                x=sub["long_wps_L"], y=sub["short_wps_L"],
                mode="markers", marker=dict(size=3, color=color, opacity=0.5),
                name=name, showlegend=False,
                hovertemplate="marker=%{customdata}<br>long_L=%{x:.3f}<br>short_L=%{y:.3f}",
                customdata=sub["marker_id"].astype(str),
            ), row=2, col=1)

        # Row2 Col2: S-WPS short vs long scatter
        if "short_wps_S" in sub.columns and "long_wps_S" in sub.columns:
            fig.add_trace(go.Scatter(
                x=sub["long_wps_S"], y=sub["short_wps_S"],
                mode="markers", marker=dict(size=3, color=color, opacity=0.5),
                name=name, showlegend=False,
                hovertemplate="marker=%{customdata}<br>long_S=%{x:.3f}<br>short_S=%{y:.3f}",
                customdata=sub["marker_id"].astype(str),
            ), row=2, col=2)

    # y=x 기준선 (WPS 패널)
    for row, col, xcol, ycol in [(2,1,"long_wps_L","short_wps_L"),
                                   (2,2,"long_wps_S","short_wps_S")]:
        if xcol in df.columns and ycol in df.columns:
            mn = float(min(df[xcol].min(), df[ycol].min()))
            mx = float(max(df[xcol].max(), df[ycol].max()))
            fig.add_trace(go.Scatter(
                x=[mn, mx], y=[mn, mx], mode="lines",
                line=dict(color="gray", dash="dash", width=1),
                showlegend=False,
            ), row=row, col=col)

    fig.update_xaxes(title_text="Long count",   row=1, col=1)
    fig.update_yaxes(title_text="Short count",  row=1, col=1)
    fig.update_yaxes(title_text="Short ratio",  row=1, col=2)
    fig.update_xaxes(title_text="Long L-WPS",   row=2, col=1)
    fig.update_yaxes(title_text="Short L-WPS",  row=2, col=1)
    fig.update_xaxes(title_text="Long S-WPS",   row=2, col=2)
    fig.update_yaxes(title_text="Short S-WPS",  row=2, col=2)

    fig.update_layout(
        title=f"Marker WPS & Count  {sample_id}",
        height=800,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11),
    )
    fig.write_html(out_html)


def _make_viz(p: Paths, sample_id: str = "") -> None:
    cnv_viz = p.cnv_baf if _e(p.cnv_baf) else p.cnv_calls
    viz_tasks = []

    if _e(p.bins_corrected):
        viz_tasks.append(("qc_dashboard", lambda: qc_dashboard.plot_qc_dashboard(
            corrected_path=p.bins_corrected,
            ff_json_path=p.fetal_fraction if _e(p.fetal_fraction) else None,
            out_html=os.path.join(p.viz_dir, "qc_dashboard.html"),
            sample_id=sample_id,
        )))

    if _e(cnv_viz):
        viz_tasks.append(("cnv_track", lambda: cnv_track.plot_cnv_track(
            cnv_path=cnv_viz,
            out_html=os.path.join(p.viz_dir, "cnv_track.html"),
            title=f"CNV Track  {sample_id}",
        )))

    if _e(p.marker_stats):
        viz_tasks.append(("marker_wps", lambda: _plot_marker_stats(
            marker_stats_path=p.marker_stats,
            out_html=os.path.join(p.viz_dir, "marker_wps.html"),
            sample_id=sample_id,
        )))

    for name, fn in viz_tasks:
        try:
            fn()
            log.info("[Viz] %s 완료", name)
        except Exception as exc:
            log.warning("[Viz] %s 실패 — %s: %s", name, type(exc).__name__, exc)