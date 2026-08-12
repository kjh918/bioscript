"""
pipeline.py — nipt_fragmentomics 파이프라인 오케스트레이터.

Step 1  bin_extractor  : BAM → bins_raw.parquet       (count/breadth/gc/mappability)
Step 2  gc_corrector   : bins_raw → bins_corrected.parquet
Step 3  wps_compute    : BAM → bp 단위 WPS bigWig + npy  (논문 방식, --wps-mode 지정 시)
Step 3m wps_processor  : npy + marker BED → marker_wps_summary.parquet
Step 4  fetal_fraction : bins_corrected → fetal_fraction.json
Step 5  cnv_caller     : bins_corrected → cnv_calls.parquet
Step 6  baf_calculator : BAM + VCF → bins_baf.parquet → cnv_baf.parquet
"""
from __future__ import annotations

import json, logging, os, time
from datetime import datetime
from typing import Optional

from nipt_fragmentomics.core.constants import DEFAULT_BIN_SIZE, FNAME
from nipt_fragmentomics.steps import (
    bin_extractor, gc_corrector, wps_processor,
    fetal_fraction, cnv_caller, baf_calculator,
)
from nipt_fragmentomics.scripts import wps_compute
from nipt_fragmentomics.viz import cnv_track, wps_profile, qc_dashboard

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────
# 경로 관리
# ─────────────────────────────────────────────────────────────────────
class Paths:
    def __init__(self, out_dir: str):
        self.out_dir        = out_dir
        self.bins_raw       = os.path.join(out_dir, FNAME["bins_raw"])
        self.bins_corrected = os.path.join(out_dir, FNAME["bins_corrected"])
        self.marker_wps     = os.path.join(out_dir, FNAME["marker_wps"])
        self.fetal_fraction = os.path.join(out_dir, FNAME["fetal_fraction"])
        self.cnv_calls      = os.path.join(out_dir, FNAME["cnv_calls"])
        self.bins_baf       = os.path.join(out_dir, FNAME["bins_baf"])
        self.cnv_baf        = os.path.join(out_dir, FNAME["cnv_baf"])
        self.manifest       = os.path.join(out_dir, FNAME["manifest"])
        self.viz_dir        = os.path.join(out_dir, "viz")
        self.wps_dir        = os.path.join(out_dir, "wps")

    def wps_npy(self, mode: str, frag: str) -> str:
        """wps_compute 출력 npy 경로."""
        prefix = os.path.join(self.wps_dir, "wps")
        return f"{prefix}.wps_{mode}_{frag}.npy"

    def wps_prefix(self) -> str:
        return os.path.join(self.wps_dir, "wps")


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


# ─────────────────────────────────────────────────────────────────────
# 파이프라인
# ─────────────────────────────────────────────────────────────────────
def run(
    # 필수
    bam_path:   str,
    out_dir:    str,
    # 참조 파일
    bin_bed:    Optional[str] = None,
    marker_bed: Optional[str] = None,
    fasta_path: Optional[str] = None,
    bw_path:    Optional[str] = None,
    vcf_path:   Optional[str] = None,
    # bin 파라미터
    bin_size:   int   = DEFAULT_BIN_SIZE,
    min_mapq:   int   = 20,
    min_baseq:  int   = 20,
    min_mappability: float = 0.75,
    # WPS 파라미터 (wps_compute)
    wps_modes:  list[str] = None,   # ["L", "S"] 등. None 이면 WPS 생략
    wps_frags:  list[str] = None,   # ["short", "long"]. None 이면 ["long"]
    wps_extend: int = 2000,         # marker profile ±extend bp
    # CNV 파라미터
    zscore_gain: float = 3.0,
    zscore_loss: float = -3.0,
    # BAF 파라미터
    baf_af_min:    float = 0.2,
    baf_af_max:    float = 0.8,
    baf_min_depth: int   = 5,
    # 실행 옵션
    n_jobs:    int  = 4,
    resume:    bool = False,
    make_viz:  bool = True,
    sample_id: str  = "",
) -> Paths:
    """
    전체 파이프라인 실행.

    WPS 파라미터
    ------------
    wps_modes : 계산할 WPS 모드 목록
                None → WPS 생략
                ["L"] → L-WPS만 (뉴클레오솜, k=120bp, frag 120-180bp)
                ["L", "S"] → L+S WPS 모두
    wps_frags : fragment 필터
                None → ["long"] 사용
                ["short", "long"] → 둘 다 계산
    """
    os.makedirs(out_dir, exist_ok=True)
    p = Paths(out_dir)
    os.makedirs(p.viz_dir,  exist_ok=True)
    os.makedirs(p.wps_dir,  exist_ok=True)

    if wps_modes is None:
        wps_modes = []
    if wps_frags is None:
        wps_frags = ["long"]

    timings: dict[str, float] = {}
    params = dict(
        bam_path=bam_path, bin_bed=bin_bed, marker_bed=marker_bed,
        fasta_path=fasta_path, bw_path=bw_path, vcf_path=vcf_path,
        bin_size=bin_size, min_mapq=min_mapq, min_baseq=min_baseq,
        min_mappability=min_mappability,
        wps_modes=wps_modes, wps_frags=wps_frags, wps_extend=wps_extend,
        zscore_gain=zscore_gain, zscore_loss=zscore_loss,
        baf_af_min=baf_af_min, baf_af_max=baf_af_max,
        baf_min_depth=baf_min_depth, n_jobs=n_jobs,
    )

    log.info("=" * 60)
    log.info("nipt_fragmentomics  sample=%s", sample_id or "—")
    log.info("out_dir  : %s", out_dir)
    log.info("CNV bin  : %s", bin_bed or f"auto {bin_size:,} bp")
    log.info("WPS      : modes=%s  frags=%s", wps_modes or "생략", wps_frags)
    log.info("BAF      : %s", "VCF 있음" if vcf_path else "생략")
    log.info("=" * 60)

    # ── Step 1: BAM 스캔 (count/breadth/gc/mappability) ──────────
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

    # ── Step 3: WPS (bp 단위, 논문 방식) ─────────────────────────
    last_npy = None   # marker WPS 에 넘길 npy 경로 (마지막으로 계산된 것)

    if not wps_modes:
        log.info("[Step 3] WPS 생략 (--wps-mode 미지정)")
    else:
        for mode in wps_modes:
            for frag in wps_frags:
                npy_path = p.wps_npy(mode, frag)
                step_key = f"step3_wps_{mode}_{frag}"

                if resume and _e(npy_path):
                    log.info("[Step 3] 건너뜀 (resume): WPS %s-%s", mode, frag)
                    last_npy = npy_path
                    continue

                t0 = time.time()
                log.info("[Step 3] WPS 계산: mode=%s  frag=%s", mode, frag)
                wps_compute.run(
                    bam_path    = bam_path,
                    out_prefix  = p.wps_prefix(),
                    mode        = mode,
                    frag_filter = frag,
                    min_mapq    = min_mapq,
                    n_jobs      = n_jobs,
                )
                last_npy = npy_path
                timings[step_key] = round(time.time() - t0, 2)
                log.info("[Step 3] WPS %s-%s 완료  %.1f s",
                         mode, frag, timings[step_key])

    # ── Step 3m: marker WPS (npy + marker BED) ───────────────────
    if marker_bed and os.path.exists(marker_bed):
        marker_wps_done = _e(p.marker_wps)
        if resume and marker_wps_done:
            log.info("[Step 3m] 건너뜀 (resume)")
        else:
            t0 = time.time()
            wps_processor.run(
                wps_npy_path        = last_npy,
                out_marker_wps_path = p.marker_wps,
                marker_bed_path     = marker_bed,
                extend              = wps_extend,
            )
            timings["step3m_marker_wps"] = round(time.time() - t0, 2)
            log.info("[Step 3m] 완료  %.1f s", timings["step3m_marker_wps"])
    else:
        log.info("[Step 3m] marker BED 없음 — 생략")

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

    # ── Step 6: BAF (VCF 있을 때만) ──────────────────────────────
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
        log.warning("[Step 6] VCF 파일 없음: %s", vcf_path)
    else:
        log.info("[Step 6] --vcf 미지정 — BAF 생략")

    # ── Viz ───────────────────────────────────────────────────────
    if make_viz:
        t0 = time.time()
        _make_viz(p, last_npy=last_npy, sample_id=sample_id)
        timings["viz"] = round(time.time() - t0, 2)
        log.info("[Viz] 완료  %.1f s", timings["viz"])

    _save_manifest(p.manifest, params, timings)
    log.info("완료  total=%.1f s  →  %s", sum(timings.values()), out_dir)
    return p


# ─────────────────────────────────────────────────────────────────────
# 시각화
# ─────────────────────────────────────────────────────────────────────
def _make_viz(p: Paths, last_npy: Optional[str] = None, sample_id: str = "") -> None:
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

    if _e(p.marker_wps):
        viz_tasks.append(("marker_wps", lambda: wps_profile.plot_marker_wps(
            marker_wps_path=p.marker_wps,
            out_html=os.path.join(p.viz_dir, "marker_wps.html"),
            title=f"Marker WPS  {sample_id}",
        )))

    for name, fn in viz_tasks:
        try:
            fn()
            log.info("[Viz] %s 완료", name)
        except Exception as exc:
            log.warning("[Viz] %s 실패 — %s: %s", name, type(exc).__name__, exc)