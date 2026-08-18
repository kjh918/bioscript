"""
pipeline.py — nipt_fragmentomics 파이프라인 오케스트레이터.

Step 1   bin_extractor (CNV)  : BAM → bins_raw.parquet        (100kb)
Step 2   gc_corrector         : bins_raw → bins_corrected.parquet
Step 3a  bin_extractor (WPS)  : BAM → bins_wps_raw.parquet    (1kb, --wps-bin-size)
Step 3b  wps_normalizer       : bins_wps_raw → bins_wps_norm.parquet
Step 3c  marker_extractor     : BAM + marker BED → marker_stats.parquet + profiles.npy
Step 4   fetal_fraction       : bins_corrected → fetal_fraction.json
Step 5   cnv_caller           : bins_corrected → cnv_calls.parquet
Step 6   baf_calculator       : BAM + VCF → cnv_baf.parquet

WPS 두 가지 모드
-----------------
  Genome-wide (Step 3a+3b):
    1kb bin 단위 WPS median → local median 차감 / depth 보정
    → bins_wps_norm.parquet (전장 유전체 WPS track)

  Marker-wise (Step 3c):
    marker 영역별 bp 단위 WPS → (WPS - region_median) / coverage * 100
    → marker_stats.parquet (summary)
    → marker_stats_profiles.npy (bp 단위 profile, 시각화용)
"""
from __future__ import annotations

import json, logging, os, time
from datetime import datetime
from typing import Optional

from nipt_fragmentomics.core.constants import DEFAULT_BIN_SIZE, FNAME
from nipt_fragmentomics.steps import (
    bin_extractor, gc_corrector, wps_normalizer,
    marker_extractor, fetal_fraction, cnv_caller, baf_calculator,
)
from nipt_fragmentomics.viz import cnv_track, qc_dashboard

log = logging.getLogger(__name__)

WPS_BIN_SIZE_DEFAULT = 1_000   # 1kb


class Paths:
    def __init__(self, out_dir: str):
        self.out_dir          = out_dir
        self.bins_raw         = os.path.join(out_dir, FNAME["bins_raw"])
        self.bins_corrected   = os.path.join(out_dir, FNAME["bins_corrected"])
        self.bins_wps_raw     = os.path.join(out_dir, FNAME["bins_wps_raw"])
        self.bins_wps_norm    = os.path.join(out_dir, FNAME["bins_wps_norm"])
        self.marker_stats     = os.path.join(out_dir, FNAME["marker_stats"])
        self.marker_profiles  = os.path.join(out_dir, FNAME["marker_profiles"])
        self.fetal_fraction   = os.path.join(out_dir, FNAME["fetal_fraction"])
        self.cnv_calls        = os.path.join(out_dir, FNAME["cnv_calls"])
        self.bins_baf         = os.path.join(out_dir, FNAME["bins_baf"])
        self.cnv_baf          = os.path.join(out_dir, FNAME["cnv_baf"])
        self.manifest         = os.path.join(out_dir, FNAME["manifest"])
        self.viz_dir          = os.path.join(out_dir, "viz")


def _e(path: str) -> bool:
    return os.path.exists(path)


def _save_manifest(path, params, timings):
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
    wps_bin_size: int  = WPS_BIN_SIZE_DEFAULT,
    min_mapq:    int   = 20,
    min_baseq:   int   = 20,
    min_mappability: float = 0.75,
    zscore_gain: float = 3.0,
    zscore_loss: float = -3.0,
    baf_af_min:    float = 0.2,
    baf_af_max:    float = 0.8,
    baf_min_depth: int   = 5,
    n_jobs:   int  = 4,
    resume:   bool = False,
    make_viz: bool = True,
    sample_id: str = "",
) -> Paths:
    os.makedirs(out_dir, exist_ok=True)
    p = Paths(out_dir)
    os.makedirs(p.viz_dir, exist_ok=True)

    timings: dict[str, float] = {}
    params = dict(
        bam_path=bam_path, bin_bed=bin_bed, marker_bed=marker_bed,
        fasta_path=fasta_path, bw_path=bw_path, vcf_path=vcf_path,
        bin_size=bin_size, wps_bin_size=wps_bin_size,
        min_mapq=min_mapq, min_baseq=min_baseq,
        min_mappability=min_mappability,
        zscore_gain=zscore_gain, zscore_loss=zscore_loss,
        baf_af_min=baf_af_min, baf_af_max=baf_af_max,
        baf_min_depth=baf_min_depth, n_jobs=n_jobs,
    )

    log.info("=" * 60)
    log.info("nipt_fragmentomics  sample=%s", sample_id or "—")
    log.info("CNV bin  : %s", bin_bed or f"auto {bin_size:,} bp")
    log.info("WPS bin  : %d bp (genome-wide)", wps_bin_size)
    log.info("marker   : %s", marker_bed or "없음")
    log.info("BAF      : %s", "VCF 있음" if vcf_path else "생략")
    log.info("=" * 60)

    # ── Step 1: CNV bin (100kb) ──────────────────────────────────
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

    # ── Step 3a: Genome-wide WPS bin (1kb) ───────────────────────
    if resume and _e(p.bins_wps_raw):
        log.info("[Step 3a] 건너뜀 (resume)")
    else:
        t0 = time.time()
        log.info("[Step 3a] WPS 1kb bin 스캔 시작")
        bin_extractor.run(
            bam_path=bam_path, out_path=p.bins_wps_raw,
            bed_path=None,           # 자동 grid
            fasta_path=fasta_path,
            bw_path=bw_path,
            bin_size=wps_bin_size,
            min_mapq=min_mapq, n_jobs=n_jobs,
        )
        timings["step3a_wps_bin"] = round(time.time() - t0, 2)
        log.info("[Step 3a] 완료  %.1f s", timings["step3a_wps_bin"])

    # ── Step 3b: Genome-wide WPS normalization ────────────────────
    if resume and _e(p.bins_wps_norm):
        log.info("[Step 3b] 건너뜀 (resume)")
    else:
        t0 = time.time()
        wps_normalizer.run(
            raw_path=p.bins_wps_raw,
            out_path=p.bins_wps_norm,
            bin_len=wps_bin_size,
            min_mappability=min_mappability,
        )
        timings["step3b_wps_norm"] = round(time.time() - t0, 2)
        log.info("[Step 3b] 완료  %.1f s", timings["step3b_wps_norm"])

    # ── Step 3c: Marker WPS (영역별 bp 단위) ─────────────────────
    if marker_bed and os.path.exists(marker_bed):
        if resume and _e(p.marker_stats):
            log.info("[Step 3c] 건너뜀 (resume)")
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
                save_profiles=True,
            )
            timings["step3c_marker"] = round(time.time() - t0, 2)
            log.info("[Step 3c] 완료  %.1f s", timings["step3c_marker"])
    else:
        log.info("[Step 3c] marker BED 없음 — 생략")

    # ── Step 4: Fetal Fraction ────────────────────────────────────
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
            timings["step6_baf"] = round(time.time() - t0, 2)
            log.info("[Step 6] 완료  %.1f s", timings["step6_baf"])
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

    if _e(p.bins_wps_norm):
        viz_tasks.append(("wps_genome_track", lambda: _plot_wps_genome_track(
            wps_norm_path=p.bins_wps_norm,
            out_html=os.path.join(p.viz_dir, "wps_genome_track.html"),
            sample_id=sample_id,
        )))

    if _e(p.marker_profiles):
        viz_tasks.append(("marker_wps_profile", lambda: _plot_marker_profiles(
            profiles_npy=p.marker_profiles,
            marker_stats_path=p.marker_stats if _e(p.marker_stats) else None,
            out_html=os.path.join(p.viz_dir, "marker_wps_profiles.html"),
            sample_id=sample_id,
        )))

    for name, fn in viz_tasks:
        try:
            fn()
            log.info("[Viz] %s 완료", name)
        except Exception as exc:
            log.warning("[Viz] %s 실패 — %s: %s", name, type(exc).__name__, exc)


def _plot_wps_genome_track(
    wps_norm_path: str,
    out_html: str,
    sample_id: str = "",
) -> None:
    """bins_wps_norm.parquet → genome-wide WPS track (short/long × L/S)."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pandas as pd
    import numpy as np

    CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    df = pd.read_parquet(wps_norm_path)
    df = df[df["chrom"].isin(CHROM_ORDER)].copy()
    df["chrom"] = pd.Categorical(df["chrom"], categories=CHROM_ORDER, ordered=True)
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)

    # 유전체 x 좌표
    offsets: dict[str, int] = {}
    cursor = 0
    for chrom in CHROM_ORDER:
        if chrom not in df["chrom"].values: continue
        offsets[chrom] = cursor
        cursor += int(df[df["chrom"] == chrom]["end"].max())
    x = np.array([offsets.get(r.chrom, 0) + (r.start + r.end) // 2
                  for r in df.itertuples()])

    cols = [c for c in [
        "short_wps_L_norm", "long_wps_L_norm",
        "short_wps_S_norm", "long_wps_S_norm",
    ] if c in df.columns]

    if not cols:
        log.warning("WPS norm 컬럼 없음 — genome track 생략")
        return

    colors = ["rgba(220,80,80,0.7)", "rgba(50,120,220,0.7)",
              "rgba(180,60,60,0.7)", "rgba(40,100,200,0.7)"]
    labels = {"short_wps_L_norm": "Short L-WPS",
              "long_wps_L_norm":  "Long  L-WPS",
              "short_wps_S_norm": "Short S-WPS",
              "long_wps_S_norm":  "Long  S-WPS"}

    # shared y 범위 (상하위 1%)
    combined = np.concatenate([df[c].dropna().values for c in cols])
    if len(combined):
        lo = float(np.quantile(combined, 0.01))
        hi = float(np.quantile(combined, 0.99))
        pad = (hi - lo) * 0.15
        y_min, y_max = lo - pad, hi + pad
    else:
        y_min, y_max = -1.0, 1.0

    n = len(cols)
    fig = make_subplots(rows=n, cols=1, shared_xaxes=True,
                        vertical_spacing=0.02,
                        subplot_titles=[labels.get(c, c) for c in cols])

    for i, col in enumerate(cols):
        fig.add_trace(go.Scatter(
            x=x, y=df[col].values,
            mode="lines",
            line=dict(color=colors[i % len(colors)], width=0.6),
            name=labels.get(col, col),
        ), row=i+1, col=1)
        fig.update_yaxes(range=[y_min, y_max], row=i+1, col=1,
                         showgrid=True, gridcolor="rgba(200,200,200,0.3)")
        fig.add_hline(y=0, line_dash="dash",
                      line_color="rgba(0,0,0,0.3)", line_width=0.6,
                      row=i+1, col=1)

    shapes = [dict(type="line", x0=off, x1=off, y0=0, y1=1,
                   yref="paper",
                   line=dict(color="lightgray", width=0.5, dash="dot"))
              for off in offsets.values()]

    fig.update_layout(
        title=f"WPS Genome Track (1kb)  {sample_id}",
        height=900, shapes=shapes,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=10), showlegend=True,
    )
    fig.update_xaxes(showticklabels=False, showgrid=False)
    fig.write_html(out_html)


def _plot_marker_profiles(
    profiles_npy:      str,
    marker_stats_path: Optional[str],
    out_html:          str,
    sample_id:         str = "",
    max_markers:       int = 5000,
) -> None:
    """
    marker_stats_profiles.npy → marker별 WPS 프로파일.

    marker_type 별로 평균 프로파일 오버레이.
    x축 = marker center 기준 상대 위치 (bp).
    y축 = adjusted WPS = (WPS - median) / coverage * 100.
    short/long 각각 별도 row.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pandas as pd
    import numpy as np
    from scipy.signal import savgol_filter

    profiles = np.load(profiles_npy, allow_pickle=True).item()
    if not profiles:
        log.warning("profiles 없음")
        return

    # marker_type 매핑
    type_map: dict[str, str] = {}
    if marker_stats_path and os.path.exists(marker_stats_path):
        ms = pd.read_parquet(marker_stats_path, columns=["marker_id", "marker_type"])
        type_map = dict(zip(ms["marker_id"].astype(str),
                            ms["marker_type"].astype(str)))

    # marker_type별 프로파일 수집
    from collections import defaultdict
    type_profiles: dict[str, dict[str, list]] = defaultdict(lambda: {
        "short_wps_L": [], "long_wps_L": [],
        "short_wps_S": [], "long_wps_S": [],
        "pos": None,
    })

    sampled = list(profiles.items())[:max_markers]
    for mid, pdata in sampled:
        mtype = type_map.get(str(mid), "marker")
        pos   = pdata.get("pos")
        if pos is None:
            continue
        if type_profiles[mtype]["pos"] is None:
            type_profiles[mtype]["pos"] = pos

        for key in ("short_wps_L", "long_wps_L", "short_wps_S", "long_wps_S"):
            arr = pdata.get(key)
            if arr is not None and len(arr) == len(pos):
                type_profiles[mtype][key].append(arr)

    if not type_profiles:
        log.warning("type_profiles 없음")
        return

    mtypes = list(type_profiles.keys())
    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#eda100", "#4a3aa7"]

    track_pairs = [
        ("short_wps_L", "long_wps_L", "L-WPS (뉴클레오솜)"),
        ("short_wps_S", "long_wps_S", "S-WPS (TF footprint)"),
    ]

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=[t[2] for t in track_pairs],
    )

    # 공유 y 범위 계산
    all_vals = []
    for mtype, tdata in type_profiles.items():
        for key in ("short_wps_L", "long_wps_L", "short_wps_S", "long_wps_S"):
            arrs = tdata[key]
            if arrs:
                med = np.nanmedian(np.array(arrs, dtype=np.float32), axis=0)
                all_vals.extend(med[np.isfinite(med)].tolist())
    if all_vals:
        lo = float(np.quantile(all_vals, 0.01))
        hi = float(np.quantile(all_vals, 0.99))
        pad = (hi - lo) * 0.15
        y_min, y_max = lo - pad, hi + pad
    else:
        y_min, y_max = -5.0, 5.0

    def _smooth(arr):
        finite = np.where(np.isfinite(arr), arr, 0.0)
        n = len(finite)
        wl = min(21, n)
        if wl % 2 == 0: wl -= 1
        if wl >= 3:
            return savgol_filter(finite, window_length=wl, polyorder=2)
        return finite

    for row_idx, (short_key, long_key, title) in enumerate(track_pairs, start=1):
        for ci, mtype in enumerate(mtypes):
            tdata  = type_profiles[mtype]
            pos    = tdata["pos"]
            color  = colors[ci % len(colors)]
            if pos is None:
                continue

            for key, dash, label_suffix in [
                (short_key, "solid", "short"),
                (long_key,  "dot",   "long"),
            ]:
                arrs = tdata[key]
                if not arrs:
                    continue
                mat  = np.array(arrs, dtype=np.float32)
                # 상하위 1% 제외 trimmed median
                med  = np.nanmedian(mat, axis=0)
                smooth = _smooth(med)

                fig.add_trace(go.Scatter(
                    x=pos.tolist(),
                    y=smooth.tolist(),
                    mode="lines",
                    line=dict(color=color, width=1.5 if dash=="solid" else 1.0,
                              dash=dash),
                    name=f"{mtype} ({label_suffix}, n={len(arrs)})",
                    showlegend=(row_idx == 1),
                    hovertemplate=f"{mtype}/{label_suffix}<br>pos=%{{x}}<br>adjWPS=%{{y:.2f}}",
                ), row=row_idx, col=1)

        fig.update_yaxes(range=[y_min, y_max], row=row_idx, col=1,
                         showgrid=True, gridcolor="rgba(200,200,200,0.3)",
                         title_text="adjusted WPS")
        fig.add_vline(x=0, line_dash="dash",
                      line_color="black", line_width=1.0,
                      annotation_text="center", row=row_idx, col=1)
        fig.add_hline(y=0, line_dash="dot",
                      line_color="gray", line_width=0.6,
                      row=row_idx, col=1)

    fig.update_xaxes(title_text="Distance from marker center (bp)", row=2, col=1)
    fig.update_layout(
        title=f"Marker WPS Profile  {sample_id}",
        height=700,
        plot_bgcolor="white", paper_bgcolor="white",
        font=dict(size=11), showlegend=True,
        legend=dict(orientation="h", y=1.04, x=0),
    )
    fig.write_html(out_html)
