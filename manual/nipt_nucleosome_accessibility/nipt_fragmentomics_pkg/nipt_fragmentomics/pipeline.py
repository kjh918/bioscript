"""
pipeline.py — nipt_fragmentomics 파이프라인 오케스트레이터.

Step 1   bin_extractor        : BAM → bins_raw.parquet
Step 2   gc_corrector         : bins_raw → bins_corrected.parquet
Step 3   wps_compute          : BAM → genome-wide normalized WPS (L/S)
Step 3c  marker_extractor     : BAM + marker BED → marker_stats.parquet + profiles.npy
Step 4   fetal_fraction       : bins_corrected → fetal_fraction.json
Step 5   cnv_caller           : bins_corrected → cnv_calls.parquet
Step 6   baf_calculator       : BAM + VCF → cnv_baf.parquet

Genome-wide WPS 결과
--------------------
  out_dir/wps/L/chr1.npy ... chrY.npy   (adjusted WPS)
  out_dir/wps/S/chr1.npy ... chrY.npy
  out_dir/wps/wps_results.json
  out_dir/viz/wps_genome_track.png
"""
from __future__ import annotations

import json, logging, os, time
from datetime import datetime
from typing import Optional

from nipt_fragmentomics.core.constants import DEFAULT_BIN_SIZE, FNAME
from nipt_fragmentomics.steps import (
    bin_extractor, gc_corrector,
    marker_extractor, fetal_fraction, cnv_caller, baf_calculator,
    wps_compute,
)          # ← scripts 패키지
from nipt_fragmentomics.viz import cnv_track, qc_dashboard

log = logging.getLogger(__name__)

WPS_WIN_DEFAULT = 1_000   # adjusted WPS windowMedian 크기 (bp)


# ─────────────────────────────────────────────────────────────────────
# 경로 관리
# ─────────────────────────────────────────────────────────────────────
class Paths:
    def __init__(self, out_dir: str):
        self.out_dir          = out_dir
        self.bins_raw         = os.path.join(out_dir, FNAME["bins_raw"])
        self.bins_corrected   = os.path.join(out_dir, FNAME["bins_corrected"])
        self.marker_stats     = os.path.join(out_dir, FNAME["marker_stats"])
        self.marker_profiles  = os.path.join(out_dir, FNAME["marker_profiles"])
        self.fetal_fraction   = os.path.join(out_dir, FNAME["fetal_fraction"])
        self.cnv_calls        = os.path.join(out_dir, FNAME["cnv_calls"])
        self.bins_baf         = os.path.join(out_dir, FNAME["bins_baf"])
        self.cnv_baf          = os.path.join(out_dir, FNAME["cnv_baf"])
        self.manifest         = os.path.join(out_dir, FNAME["manifest"])
        self.viz_dir          = os.path.join(out_dir, "viz")
        self.wps_dir          = os.path.join(out_dir, "wps")
        # WPS 결과 manifest (Step 3에서 생성)
        self.wps_results      = os.path.join(self.wps_dir, "wps_results.json")
        self.wps_genome_png   = os.path.join(self.viz_dir, "wps_genome_track.png")


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


# ─────────────────────────────────────────────────────────────────────
# 파이프라인
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:    str,
    out_dir:     str,
    bin_bed:     Optional[str] = None,
    marker_bed:  Optional[str] = None,
    fasta_path:  Optional[str] = None,
    bw_path:     Optional[str] = None,
    vcf_path:    Optional[str] = None,
    bin_size:    int   = DEFAULT_BIN_SIZE,
    wps_win:     int   = WPS_WIN_DEFAULT,   # adjusted WPS windowMedian (bp)
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
    os.makedirs(p.wps_dir, exist_ok=True)

    timings: dict[str, float] = {}
    params = dict(
        bam_path=bam_path, bin_bed=bin_bed, marker_bed=marker_bed,
        fasta_path=fasta_path, bw_path=bw_path, vcf_path=vcf_path,
        bin_size=bin_size, wps_win=wps_win,
        min_mapq=min_mapq, min_baseq=min_baseq,
        min_mappability=min_mappability,
        zscore_gain=zscore_gain, zscore_loss=zscore_loss,
        baf_af_min=baf_af_min, baf_af_max=baf_af_max,
        baf_min_depth=baf_min_depth, n_jobs=n_jobs,
    )

    log.info("=" * 60)
    log.info("nipt_fragmentomics  sample=%s", sample_id or "—")
    log.info("CNV bin  : %s", bin_bed or f"auto {bin_size:,} bp")
    log.info("WPS win  : %d bp", wps_win)
    log.info("marker   : %s", marker_bed or "없음")
    log.info("BAF      : %s", "VCF 있음" if vcf_path else "생략")
    log.info("=" * 60)

    # ── Step 1 ────────────────────────────────────────────────────
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

    # ── Step 2 ────────────────────────────────────────────────────
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

    # ── Step 3: Genome-wide WPS (L / S) ──────────────────────────
    # L 과 S 를 ProcessPoolExecutor 로 동시 실행하여 속도 최적화
    # 각 mode 내부도 염색체 병렬 → 총 n_jobs 프로세스로 L/S 동시 처리
    if resume and _e(p.wps_results):
        log.info("[Step 3] 건너뜀 (resume)")
        with open(p.wps_results) as f:
            wps_results = json.load(f)
    else:
        from concurrent.futures import ProcessPoolExecutor as _PPE, as_completed as _ac
        t0 = time.time()
        wps_results: dict[str, dict] = {}

        # L/S 각각에 n_jobs // 2 프로세스 배분 (최소 1)
        jobs_per_mode = max(1, n_jobs // 2)
        modes_to_run  = ["L", "S"]

        log.info("[Step 3] WPS 시작: L/S 동시 병렬  jobs_per_mode=%d", jobs_per_mode)
        for mode in modes_to_run:
            prm = wps_compute.WPS_PARAMS[mode]
            log.info("  mode=%s  frag=%d-%dbp  k=%d",
                     mode, prm["frag_min"], prm["frag_max"], prm["window"])

        # ThreadPoolExecutor 로 L/S 동시 submit (각 run() 내부는 ProcessPool 사용)
        import concurrent.futures as _cf
        with _cf.ThreadPoolExecutor(max_workers=2) as tex:
            mode_futures = {
                tex.submit(
                    wps_compute.run,
                    bam_path = bam_path,
                    out_dir  = p.wps_dir,
                    mode     = mode,
                    chroms   = None,
                    min_mapq = min_mapq,
                    win_size = wps_win,
                    n_jobs   = jobs_per_mode,
                    resume   = resume,
                    prefix   = sample_id,   # 파일명: {sample_id}.{chrom}.raw/cov/norm.npy
                ): mode
                for mode in modes_to_run
            }
            for fut in _cf.as_completed(mode_futures):
                mode = mode_futures[fut]
                try:
                    result = fut.result()
                    if result:
                        wps_results[mode] = result
                    log.info("[Step 3] mode=%s 완료", mode)
                except Exception as exc:
                    log.error("[Step 3] mode=%s 실패: %s", mode, exc)

        # wps_results.json 저장
        with open(p.wps_results, "w") as f:
            json.dump(wps_results, f, indent=2)

        timings["step3_wps"] = round(time.time() - t0, 2)
        log.info("[Step 3] 완료  modes=%s  %.1f s",
                 ",".join(wps_results.keys()), timings["step3_wps"])

    # ── Step 3c: Marker WPS ───────────────────────────────────────
    if marker_bed and os.path.exists(marker_bed):
        # marker_stats.parquet 있으면 resume
        if resume and _e(p.marker_stats):
            log.info("[Step 3c] 건너뜀 (resume)")
        else:
            t0 = time.time()
            marker_extractor.run(
                marker_bed    = marker_bed,
                out_path      = p.marker_stats,
                wps_dir       = p.wps_dir,
                bam_path      = bam_path,
                fasta_path    = fasta_path,
                bw_path       = bw_path,
                min_mapq      = min_mapq,
                n_jobs        = n_jobs,
                save_profiles = True,
                make_plots    = True,
                plot_dir      = os.path.join(p.viz_dir, "marker_plots"),
                profile_flank = 1000,
                prefix        = sample_id,
            )
            timings["step3c_marker"] = round(time.time() - t0, 2)
            log.info("[Step 3c] 완료  %.1f s", timings["step3c_marker"])
    else:
        log.info("[Step 3c] marker BED 없음 — 생략")

    # ── Step 4 ────────────────────────────────────────────────────
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

    # ── Step 5 ────────────────────────────────────────────────────
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

    # ── Step 6 ────────────────────────────────────────────────────
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
        _make_viz(p, wps_results=wps_results if _e(p.wps_results) else {},
                  sample_id=sample_id)
        timings["viz"] = round(time.time() - t0, 2)
        log.info("[Viz] 완료  %.1f s", timings["viz"])

    _save_manifest(p.manifest, params, timings)
    log.info("완료  total=%.1f s  →  %s", sum(timings.values()), out_dir)
    return p


# ─────────────────────────────────────────────────────────────────────
# 시각화
# ─────────────────────────────────────────────────────────────────────
def _make_viz(
    p: Paths,
    wps_results: dict,
    sample_id: str = "",
) -> None:
    cnv_viz = p.cnv_baf if _e(p.cnv_baf) else p.cnv_calls
    tasks   = []

    if _e(p.bins_corrected):
        _bc = p.bins_corrected; _ff = p.fetal_fraction
        _oh = os.path.join(p.viz_dir, "qc_dashboard.html"); _s = sample_id
        tasks.append(("qc_dashboard",
            lambda b=_bc, f=_ff, o=_oh, s=_s:
                qc_dashboard.plot_qc_dashboard(
                    corrected_path=b,
                    ff_json_path=f if _e(f) else None,
                    out_html=o, sample_id=s)))

    if _e(cnv_viz):
        _cv = cnv_viz; _oh = os.path.join(p.viz_dir, "cnv_track.html"); _s = sample_id
        tasks.append(("cnv_track",
            lambda c=_cv, o=_oh, s=_s:
                cnv_track.plot_cnv_track(
                    cnv_path=c, out_html=o, title=f"CNV Track  {s}")))

    if _e(p.wps_results):
        _wj = p.wps_results; _wp = p.wps_genome_png; _s = sample_id
        tasks.append(("wps_genome_track",
            lambda j=_wj, o=_wp, s=_s:
                _plot_wps_genome_track(wps_results_json=j, out_png=o, sample_id=s)))

    if _e(p.marker_profiles):
        _npy = p.marker_profiles
        _ms  = p.marker_stats if _e(p.marker_stats) else None
        _oh  = os.path.join(p.viz_dir, "marker_wps_profiles.png")
        _s   = sample_id
        tasks.append(("marker_wps_profile",
            lambda n=_npy, m=_ms, o=_oh, s=_s:
                _plot_marker_profiles(profiles_npy=n, marker_stats_path=m,
                                      out_png=o, sample_id=s)))

    for name, fn in tasks:
        try:
            fn()
            log.info("[Viz] %s 완료", name)
        except Exception as exc:
            log.warning("[Viz] %s 실패 — %s: %s", name, type(exc).__name__, exc)


# ─────────────────────────────────────────────────────────────────────
# WPS genome track 시각화
# ─────────────────────────────────────────────────────────────────────
def _robust_ylim(
    values: "np.ndarray",
    q_low: float = 0.01,
    q_high: float = 0.99,
    default_half_range: float = 1.0,
    positive_only: bool = False,
) -> tuple:
    import numpy as np
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (0.0, default_half_range) if positive_only else (-default_half_range, default_half_range)
    lo = float(np.quantile(arr, q_low))
    hi = float(np.quantile(arr, q_high))
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return (0.0, default_half_range) if positive_only else (-default_half_range, default_half_range)
    if np.isclose(lo, hi):
        center = float(np.nanmedian(arr))
        pad = max(abs(center) * 0.2, default_half_range)
        if positive_only:
            return 0.0, center + pad
        return center - pad, center + pad
    pad = max((hi - lo) * 0.15, 1e-6)
    if positive_only:
        return max(0.0, lo - pad), hi + pad
    return lo - pad, hi + pad


def _plot_wps_genome_track(
    wps_results_json: str,
    out_png:  str,
    sample_id: str = "",
    window_size: int = 100_000,
) -> None:
    """
    wps_results.json → genome-wide WPS PNG.
    mode별 (raw WPS / coverage / adjusted WPS) 3트랙씩 그립니다.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def _aggregate(arr, window):
        xs, ys = [], []
        for s in range(0, int(arr.size), window):
            e   = min(s + window, int(arr.size))
            blk = np.asarray(arr[s:e], dtype=np.float32)
            fin = blk[np.isfinite(blk)]
            xs.append((s + e) // 2)
            ys.append(float(np.mean(fin)) if fin.size else np.nan)
        return np.asarray(xs, np.int64), np.asarray(ys, np.float32)

    with open(wps_results_json) as f:
        results = json.load(f)

    CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    modes = [m for m in ("L", "S") if m in results]
    if not modes:
        log.warning("[Viz] wps_results에 데이터 없음")
        return

    # metric별 path 수집 + 염색체 크기
    chrom_sizes: dict[str, int] = {}
    # path_map[mode][metric][chrom] = path
    path_map: dict[str, dict[str, dict[str, str]]] = {}

    for mode in modes:
        path_map[mode] = {"raw": {}, "cov": {}, "norm": {}, "frag_cov": {}}
        for chrom, entry in results[mode].get("chromosomes", {}).items():
            if isinstance(entry, str):
                path_map[mode]["norm"][chrom] = entry
            elif isinstance(entry, dict):
                for metric in ("raw", "cov", "norm", "frag_cov"):
                    p = entry.get(metric, "")
                    if p and os.path.isfile(p):
                        path_map[mode][metric][chrom] = p
            if chrom not in chrom_sizes:
                for metric in ("norm", "raw", "cov", "frag_cov"):
                    p = path_map[mode][metric].get(chrom, "")
                    if p and os.path.isfile(p):
                        arr = np.load(p, mmap_mode="r", allow_pickle=False)
                        chrom_sizes[chrom] = int(arr.size)
                        break

    present = [c for c in CHROM_ORDER if c in chrom_sizes]
    if not present:
        log.warning("[Viz] WPS plot용 염색체 없음")
        return

    offsets: dict[str, int] = {}
    cursor = 0
    for chrom in present:
        offsets[chrom] = cursor
        cursor += chrom_sizes[chrom]

    # 패널 정의: (mode, metric, label, positive_only, color)
    panel_defs = []
    colors     = {"L": "#3278d6", "S": "#dc5050"}
    cov_color  = "#888888"   # 전체 coverage는 회색으로 구분
    metric_labels = {
        "cov":      "Coverage (all)",
        "frag_cov": "Coverage (frag)",
        "raw":      "Raw WPS",
        "norm":     "Adjusted WPS",
    }
    for mode in modes:
        frag = "120-180bp" if mode == "L" else "35-80bp"
        # (mode, metric, label, positive_only, color)
        panel_defs.append((mode, "cov",      f"{mode} Coverage all  (total)",            True,  cov_color))
        panel_defs.append((mode, "frag_cov", f"{mode} Coverage frag ({frag})",            True,  colors[mode]))
        panel_defs.append((mode, "raw",      f"{mode}-WPS Raw WPS  ({frag})",             False, colors[mode]))
        panel_defs.append((mode, "norm",     f"{mode}-WPS Adjusted WPS  ({frag})",        False, colors[mode]))

    # 데이터 수집
    panel_data: dict = {}
    for mode, metric, *_ in panel_defs:
        xs_all, ys_all = [], []
        for chrom in present:
            p = path_map.get(mode, {}).get(metric, {}).get(chrom, "")
            if not p or not os.path.isfile(p):
                continue
            arr = np.load(p, mmap_mode="r", allow_pickle=False)
            xl, yl = _aggregate(arr, window_size)
            if xl.size == 0:
                continue
            xs_all.append(xl + offsets[chrom])
            ys_all.append(yl)
            xs_all.append(np.array([offsets[chrom] + chrom_sizes[chrom]], np.int64))
            ys_all.append(np.array([np.nan], np.float32))
        if xs_all:
            panel_data[(mode, metric)] = (
                np.concatenate(xs_all),
                np.concatenate(ys_all),
            )

    active = [(m, mt, lb, po, co) for (m, mt, lb, po, co) in panel_defs
              if (m, mt) in panel_data]
    if not active:
        log.warning("[Viz] WPS plot 데이터 없음")
        return

    n = len(active)
    fig, axes = plt.subplots(n, 1, figsize=(22, 2.8 * n),
                             sharex=True, facecolor="white")
    if n == 1:
        axes = [axes]

    for ax, (mode, metric, label, positive_only, color) in zip(axes, active):
        x, y = panel_data[(mode, metric)]
        ax.fill_between(x, 0, y,
                        where=np.isfinite(y) & (y >= 0), color=color, alpha=0.7, lw=0)
        ax.fill_between(x, 0, y,
                        where=np.isfinite(y) & (y < 0),  color=color, alpha=0.4, lw=0)
        ax.axhline(0, linewidth=0.5, linestyle="--", color="black", alpha=0.4)

        fin = y[np.isfinite(y)]
        if fin.size:
            lo = float(np.quantile(fin, 0.005))
            hi = float(np.quantile(fin, 0.995))
            pad = (hi - lo) * 0.1
            if positive_only:
                ax.set_ylim(0, max(hi + pad, 1.0))
            else:
                bound = max(abs(lo - pad), abs(hi + pad))
                ax.set_ylim(-bound, bound)

        ax.set_ylabel(label, fontsize=8)
        ax.set_facecolor("white")
        ax.grid(axis="y", linewidth=0.3, alpha=0.3)
        for sp in ax.spines.values(): sp.set_linewidth(0.4)

        for chrom in present:
            ax.axvline(offsets[chrom], linewidth=0.3, linestyle=":", alpha=0.4, color="gray")

    centers = [offsets[c] + chrom_sizes[c] // 2 for c in present]
    axes[-1].set_xticks(centers)
    axes[-1].set_xticklabels([c.replace("chr", "") for c in present], fontsize=8)
    axes[-1].set_xlabel("Chromosome", fontsize=10)

    fig.suptitle("Genome-wide WPS  " + sample_id, fontsize=12)
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info("[Viz] genome-wide WPS 저장: %s", out_png)



# ─────────────────────────────────────────────────────────────────────
# Marker WPS 프로파일 시각화
# ─────────────────────────────────────────────────────────────────────
def _plot_marker_profiles(
    profiles_npy:      str,
    marker_stats_path: Optional[str],
    out_png:           str,
    sample_id:         str = "",
) -> None:
    """marker_stats_profiles.npy -> aggregate WPS profile PNG."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        data = np.load(profiles_npy, allow_pickle=True).item()
    except Exception as e:
        log.error("profiles npy load failed: %s", e)
        return

    pos           = data.get("pos")
    marker_count  = data.get("marker_count", 0)
    if pos is None:
        log.warning("profiles npy: pos 없음")
        return
    pos = np.asarray(pos, dtype=np.int32)

    panel_order = [
        ("L_cov",  "L Coverage (120-180bp)",  True,  "#2a78d6"),
        ("L_raw",  "L Raw WPS",               False, "#2a78d6"),
        ("L_norm", "L Adjusted WPS",          False, "#2a78d6"),
        ("S_cov",  "S Coverage (35-80bp)",    True,  "#eb6834"),
        ("S_raw",  "S Raw WPS",               False, "#eb6834"),
        ("S_norm", "S Adjusted WPS",          False, "#eb6834"),
    ]

    fig, axes = plt.subplots(6, 1, figsize=(14, 12), sharex=True, facecolor="white")
    for ax, (key, ylabel, positive_only, color) in zip(axes, panel_order):
        arr = data.get(key)
        if arr is None:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
            ax.set_ylabel(ylabel, fontsize=9)
            continue
        arr    = np.asarray(arr, dtype=np.float32)
        finite = np.isfinite(arr)
        if not finite.any():
            ax.text(0.5, 0.5, "All NaN", transform=ax.transAxes, ha="center", va="center")
            ax.set_ylabel(ylabel, fontsize=9)
            continue
        ax.plot(pos, arr, linewidth=1.0, color=color)
        ax.axvline(0, linewidth=0.8, linestyle="--", alpha=0.6, color="black")
        vals = arr[finite]
        vmin, vmax = float(vals.min()), float(vals.max())
        if positive_only:
            ax.set_ylim(0, max(vmax * 1.1, 1.0))
        else:
            pad = max((vmax - vmin) * 0.1, 0.01) if not np.isclose(vmin, vmax) else 0.01
            ax.set_ylim(vmin - pad, vmax + pad)
            ax.axhline(0, linewidth=0.5, linestyle="--", alpha=0.4, color="gray")
        n_arr = data.get(key + "_n")
        n_max = int(n_arr.max()) if n_arr is not None else marker_count
        label = ylabel + " (n<=" + str(n_max) + ")"
        ax.set_ylabel(label, fontsize=9)
        ax.grid(axis="y", linewidth=0.3, alpha=0.25)
        ax.set_facecolor("white")
        for sp in ax.spines.values():
            sp.set_linewidth(0.4)

    axes[-1].set_xlabel("Distance from marker center (bp)", fontsize=10)
    title = "Aggregate Marker WPS  " + sample_id + "  (markers=" + str(marker_count) + ")"
    fig.suptitle(title, fontsize=13)
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info("marker WPS profile saved: %s", out_png)