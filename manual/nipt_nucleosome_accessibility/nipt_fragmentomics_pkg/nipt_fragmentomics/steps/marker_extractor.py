"""
marker_extractor.py
===================
marker BED 영역별로 fragment 통계와 WPS 를 계산합니다.

설계 원칙
---------
- count (short/long) : fragment midpoint 가 marker 영역 내에 있을 때 귀속
- coverage (breadth) : fragment 가 marker 영역과 overlap 될 때 누적
- WPS               : fragment 가 marker 영역과 overlap 될 때 spanning/endpoints 누적
                      (midpoint 귀속과 다름 — WPS 는 영역 전체에서 계산)

같은 염색체의 marker 를 1회 BAM fetch 로 일괄 처리합니다.
각 read 를:
  - midpoint 기반으로 단일 marker 에 count 귀속
  - overlap 기반으로 모든 겹치는 marker 에 WPS/coverage 누적

출력 (marker_stats.parquet)
---------------------------
  chrom, start, end, marker_id, marker_type,
  short_count, long_count, total_count, short_ratio,
  short_wps_L, long_wps_L, short_wps_S, long_wps_S,
  short_breadth, long_breadth,
  gc, mappability
"""

from __future__ import annotations

import bisect
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pandas as pd
import pysam

from nipt_fragmentomics.core.constants import MIN_MAPQ, WPS_PARAMS, SHORT_MAX
from nipt_fragmentomics.core.schema import FragmentScore

log = logging.getLogger(__name__)

try:
    import pyBigWig as _pbw
    _HAS_BW = True
except ImportError:
    _HAS_BW = False


# ─────────────────────────────────────────────────────────────────────
# marker BED 로더
# ─────────────────────────────────────────────────────────────────────
def load_marker_bed(bed_path: str) -> pd.DataFrame:
    df = pd.read_csv(
        bed_path, sep="\t", header=None, comment="#",
        names=["chrom", "start", "end", "marker_id", "marker_type"],
    )
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    if "marker_id" not in df.columns or df["marker_id"].isna().all():
        df["marker_id"] = (df["chrom"].astype(str) + ":" +
                           df["start"].astype(str) + "-" + df["end"].astype(str))
    if "marker_type" not in df.columns or df["marker_type"].isna().all():
        df["marker_type"] = "marker"
    log.info("marker BED: %d markers [%s]", len(df), bed_path)
    return df


# ─────────────────────────────────────────────────────────────────────
# marker 단위 누적기
# ─────────────────────────────────────────────────────────────────────
class _MarkerAccumulator:
    """
    단일 marker region 의 지표를 누적합니다.

    count  : midpoint 기반 (marker 에 딱 하나만 귀속)
    WPS    : overlap 기반 (fragment 가 marker 영역에 걸치면 누적)
    breadth: overlap 기반 (coverage mask)
    """
    __slots__ = (
        "marker_id", "marker_type", "chrom", "start", "end", "marker_len",
        "short_count", "long_count",
        "short_cov", "long_cov",
        "s_L_span", "s_L_ep",
        "s_S_span", "s_S_ep",
        "l_L_span", "l_L_ep",
        "l_S_span", "l_S_ep",
    )

    def __init__(self, marker_id: str, marker_type: str,
                 chrom: str, start: int, end: int):
        self.marker_id   = marker_id
        self.marker_type = marker_type
        self.chrom       = chrom
        self.start       = start
        self.end         = end
        self.marker_len  = end - start

        self.short_count = 0
        self.long_count  = 0
        self.short_cov   = np.zeros(self.marker_len, dtype=bool)
        self.long_cov    = np.zeros(self.marker_len, dtype=bool)

        for attr in ("s_L_span", "s_L_ep", "s_S_span", "s_S_ep",
                     "l_L_span", "l_L_ep", "l_S_span", "l_S_ep"):
            setattr(self, attr, np.zeros(self.marker_len, dtype=np.int32))

    def add_count(self, is_short: bool) -> None:
        """midpoint 기반 count 귀속."""
        if is_short:
            self.short_count += 1
        else:
            self.long_count  += 1

    def add_overlap(self, fs: FragmentScore) -> None:
        """
        overlap 기반 WPS / coverage 누적.
        fragment 가 marker 영역과 겹칠 때 호출됩니다.
        """
        bs, be = self.start, self.end

        # coverage
        lo = max(fs.frag_start, bs) - bs
        hi = min(fs.frag_end,   be) - bs
        if lo < hi:
            if fs.is_short:
                self.short_cov[lo:hi] = True
            else:
                self.long_cov[lo:hi]  = True

        # WPS: L / S 모드 각각
        for mode, prm in WPS_PARAMS.items():
            if not (prm["frag_min"] <= fs.frag_len <= prm["frag_max"]):
                continue
            hk = prm["window"] // 2

            if fs.is_short:
                span_arr = getattr(self, f"s_{mode}_span")
                ep_arr   = getattr(self, f"s_{mode}_ep")
            else:
                span_arr = getattr(self, f"l_{mode}_span")
                ep_arr   = getattr(self, f"l_{mode}_ep")

            # spanning
            sp_lo = max(fs.frag_start + hk + 1, bs) - bs
            sp_hi = min(fs.frag_end   - hk,     be) - bs
            if sp_lo < sp_hi:
                span_arr[sp_lo:sp_hi] += 1

            # endpoints ±half_k 윈도우
            for ep in (fs.frag_start, fs.frag_end - 1):
                w_lo = max(ep - hk,     bs) - bs
                w_hi = min(ep + hk + 1, be) - bs
                if w_lo < w_hi:
                    ep_arr[w_lo:w_hi] += 1

    def to_row(self, gc: float, mappability: float) -> dict:
        ml    = float(self.marker_len)
        total = self.short_count + self.long_count

        def _wps_median(span, ep) -> float:
            arr = (span - ep).astype(np.float32)
            return float(np.median(arr))

        return {
            "chrom":         self.chrom,
            "start":         self.start,
            "end":           self.end,
            "marker_id":     self.marker_id,
            "marker_type":   self.marker_type,
            "short_count":   self.short_count,
            "long_count":    self.long_count,
            "total_count":   total,
            "short_ratio":   float(self.short_count / total) if total > 0 else float("nan"),
            "short_wps_L":   _wps_median(self.s_L_span, self.s_L_ep),
            "long_wps_L":    _wps_median(self.l_L_span, self.l_L_ep),
            "short_wps_S":   _wps_median(self.s_S_span, self.s_S_ep),
            "long_wps_S":    _wps_median(self.l_S_span, self.l_S_ep),
            "short_breadth": float(self.short_cov.sum() / ml),
            "long_breadth":  float(self.long_cov.sum()  / ml),
            "gc":            gc,
            "mappability":   mappability,
        }

    def to_profile(self) -> dict[str, np.ndarray]:
        """
        bp 단위 adjusted WPS + raw WPS + coverage 배열 반환.
        x축 = marker center 기준 상대 위치.

        adjusted_WPS = (WPS(bp) - median(WPS_region)) / coverage(bp) * 100
          Snyder et al. normalizeWPSwigs.py 방식
          coverage = 0 → NaN
        """
        center = self.marker_len // 2
        pos    = np.arange(self.marker_len) - center

        short_cov_int = self.short_cov.astype(np.int32)
        long_cov_int  = self.long_cov.astype(np.int32)

        def _adj(span_arr, ep_arr, cov_arr):
            wps      = (span_arr - ep_arr).astype(np.float64)
            med      = float(np.median(wps))
            adjusted = wps - med
            out      = np.full(len(wps), np.nan, dtype=np.float32)
            mask     = cov_arr > 0
            out[mask] = (adjusted[mask] / cov_arr[mask] * 100.0).astype(np.float32)
            return out

        return {
            "pos":             pos,
            # adjusted WPS (normalized)
            "short_wps_L":    _adj(self.s_L_span, self.s_L_ep, short_cov_int),
            "long_wps_L":     _adj(self.l_L_span, self.l_L_ep, long_cov_int),
            "short_wps_S":    _adj(self.s_S_span, self.s_S_ep, short_cov_int),
            "long_wps_S":     _adj(self.l_S_span, self.l_S_ep, long_cov_int),
            # raw WPS (비교용)
            "short_wps_L_raw":(self.s_L_span - self.s_L_ep).astype(np.float32),
            "long_wps_L_raw": (self.l_L_span - self.l_L_ep).astype(np.float32),
            "short_cov":      short_cov_int,
            "long_cov":       long_cov_int,
        }


# ─────────────────────────────────────────────────────────────────────
# 어노테이션 헬퍼
# ─────────────────────────────────────────────────────────────────────
def _gc_of_region(fasta, chrom, start, end) -> float:
    if fasta is None:
        return float("nan")
    try:
        seq = fasta.fetch(chrom, start, end).upper()
        return (seq.count("G") + seq.count("C")) / len(seq) if seq else float("nan")
    except Exception:
        return float("nan")


def _map_of_region(bw, chrom, start, end) -> float:
    if bw is None:
        return float("nan")
    try:
        val = bw.stats(chrom, start, end, type="mean")[0]
        return float(val) if val is not None else float("nan")
    except Exception:
        return float("nan")


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 스캔 (1회 fetch)
# ─────────────────────────────────────────────────────────────────────
def _scan_chrom_markers(
    chrom:      str,
    markers_df: pd.DataFrame,
    bam_path:   str,
    fasta_path: Optional[str],
    bw_path:    Optional[str],
    min_mapq:   int,
) -> tuple[list[dict], dict[str, dict]]:
    """
    같은 염색체의 모든 marker 를 1회 BAM fetch 로 처리합니다.

    각 read 처리:
      1. count  : midpoint 가 속한 단일 marker 에만 귀속
      2. WPS    : fragment 가 겹치는 모든 marker 에 overlap 누적
      3. breadth: WPS 와 동일 (overlap 기반)
    """
    if markers_df.empty:
        return [], {}

    # Accumulator 초기화
    accumulators: dict[str, _MarkerAccumulator] = {}
    for row in markers_df.itertuples(index=False):
        mid_id = str(row.marker_id)
        accumulators[mid_id] = _MarkerAccumulator(
            marker_id   = mid_id,
            marker_type = str(row.marker_type),
            chrom       = chrom,
            start       = int(row.start),
            end         = int(row.end),
        )

    # GC / Mappability 사전 계산
    fasta_h = pysam.FastaFile(fasta_path) if fasta_path else None
    bw_h    = (_pbw.open(bw_path) if (_HAS_BW and bw_path) else None)
    gc_vals  = {str(r.marker_id): _gc_of_region(fasta_h, chrom, int(r.start), int(r.end))
                for r in markers_df.itertuples(index=False)}
    map_vals = {str(r.marker_id): _map_of_region(bw_h, chrom, int(r.start), int(r.end))
                for r in markers_df.itertuples(index=False)}
    if fasta_h: fasta_h.close()
    if bw_h:    bw_h.close()

    # bisect 구조 (정렬된 marker 목록)
    sdf = markers_df.sort_values("start").reset_index(drop=True)
    mk_starts  = sdf["start"].tolist()
    mk_ends    = sdf["end"].tolist()
    mk_ids     = sdf["marker_id"].astype(str).tolist()

    # fetch 범위: 전체 marker 범위 + WPS window padding
    max_wps_window = max(p["window"] for p in WPS_PARAMS.values())
    fetch_start = max(0, int(sdf["start"].min()) - max_wps_window)
    fetch_end   = int(sdf["end"].max()) + max_wps_window

    seen: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
        for read in bam.fetch(chrom, fetch_start, fetch_end):
            if (read.is_unmapped or read.is_duplicate
                    or read.is_secondary or read.is_supplementary):
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.is_paired and not read.is_read1:
                continue

            qname = read.query_name
            if qname in seen:
                continue
            seen.add(qname)

            fs = FragmentScore.from_read(read)
            if fs is None:
                continue

            # ── 1. count: midpoint 기반 단일 marker 귀속 ────────────
            mid = fs.midpoint
            cnt_idx = bisect.bisect_right(mk_starts, mid) - 1
            if 0 <= cnt_idx < len(mk_ids):
                if mk_starts[cnt_idx] <= mid < mk_ends[cnt_idx]:
                    accumulators[mk_ids[cnt_idx]].add_count(fs.is_short)

            # ── 2. WPS / breadth: overlap 기반 모든 marker 누적 ─────
            # fragment [frag_start, frag_end) 와 겹치는 marker 탐색
            # 겹침 조건: mk_start < frag_end AND mk_end > frag_start
            lo_idx = bisect.bisect_left(mk_ends,    fs.frag_start + 1)
            hi_idx = bisect.bisect_right(mk_starts, fs.frag_end   - 1)

            for idx in range(lo_idx, hi_idx):
                accumulators[mk_ids[idx]].add_overlap(fs)

    rows     = [acc.to_row(gc_vals[mid], map_vals[mid])
                for mid, acc in accumulators.items()]
    profiles = {mid: acc.to_profile()
                for mid, acc in accumulators.items()}

    return rows, profiles


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker (ProcessPoolExecutor pickle 호환)
# ─────────────────────────────────────────────────────────────────────
def _chrom_worker(args: tuple) -> tuple[list[dict], dict]:
    chrom, markers_df, bam_path, fasta_path, bw_path, min_mapq = args
    return _scan_chrom_markers(
        chrom, markers_df, bam_path, fasta_path, bw_path, min_mapq
    )


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:     str,
    marker_bed:   str,
    out_path:     str,
    fasta_path:   Optional[str] = None,
    bw_path:      Optional[str] = None,
    min_mapq:     int = MIN_MAPQ,
    n_jobs:       int = 4,
    save_profiles: bool = True,   # profile npy 저장 여부
) -> pd.DataFrame:
    """
    marker BED → marker_stats.parquet (+ profiles.npy)

    marker_stats.parquet:
      marker 별 count/ratio/WPS median/breadth/gc/mappability

    profiles.npy (save_profiles=True):
      {marker_id: {"pos", "short_wps_L", "long_wps_L", "short_wps_S",
                   "long_wps_S", "short_cov", "long_cov"}}
      x축 = marker center 기준 상대 위치
    """
    marker_df = load_marker_bed(marker_bed)
    if marker_df.empty:
        log.warning("marker 없음")
        return pd.DataFrame()

    chrom_groups = {c: g.reset_index(drop=True)
                    for c, g in marker_df.groupby("chrom")}
    log.info("marker 스캔: %d markers, %d chroms",
             len(marker_df), len(chrom_groups))

    task_args = [
        (chrom, cdf, bam_path, fasta_path, bw_path, min_mapq)
        for chrom, cdf in chrom_groups.items()
    ]

    all_rows: list[dict]        = []
    all_profiles: dict[str, dict] = {}

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_chrom_worker, a): a[0] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                rows, profiles = fut.result()
                all_rows.extend(rows)
                all_profiles.update(profiles)
                log.info("  ✓ %s (%d markers)", chrom, len(rows))
            except Exception as e:
                log.error("  ✗ %s: %s", chrom, e)

    df = pd.DataFrame(all_rows)
    if df.empty:
        log.warning("결과 없음")
        return df

    for c in ["short_count", "long_count", "total_count"]:
        if c in df: df[c] = df[c].astype("int32")
    for c in ["short_ratio", "short_wps_L", "long_wps_L",
              "short_wps_S", "long_wps_S",
              "short_breadth", "long_breadth", "gc", "mappability"]:
        if c in df: df[c] = df[c].astype("float32")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("marker_stats 저장: %s (%d rows)", out_path, len(df))

    # profile npy 저장
    if save_profiles and all_profiles:
        npy_path = out_path.replace(".parquet", "_profiles.npy")
        np.save(npy_path, all_profiles, allow_pickle=True)
        log.info("profiles npy 저장: %s (%.1f MB)",
                 npy_path,
                 os.path.getsize(npy_path) / 1e6 if os.path.exists(npy_path) else 0)

    return df
