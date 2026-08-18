"""
marker_extractor.py
===================
marker BED 영역별로 fragment 통계와 WPS 를 계산합니다.

bin_extractor 와 동일한 로직을 사용하지만
입력이 100kb grid 가 아닌 marker BED (크기 무관) 입니다.

출력 (marker_stats.parquet)
---------------------------
  chrom, start, end, marker_id, marker_type,
  short_count, long_count, total_count, short_ratio,
  short_wps_L, long_wps_L,
  short_wps_S, long_wps_S,
  short_breadth, long_breadth,
  gc, mappability
"""

from __future__ import annotations

import bisect
import logging
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pandas as pd
import pysam

from nipt_fragmentomics.core.constants import MIN_MAPQ, WPS_PARAMS
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
    """
    marker BED 로드.
    포맷: chrom  start  end  marker_id  marker_type  (header 없음)
    """
    df = pd.read_csv(
        bed_path, sep="\t", header=None, comment="#",
        names=["chrom", "start", "end", "marker_id", "marker_type"],
    )
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)

    if "marker_id" not in df.columns or df["marker_id"].isna().all():
        df["marker_id"] = (df["chrom"].astype(str) + ":" +
                           df["start"].astype(str) + "-" +
                           df["end"].astype(str))
    if "marker_type" not in df.columns or df["marker_type"].isna().all():
        df["marker_type"] = "marker"

    log.info("marker BED 로드: %d markers  [%s]", len(df), bed_path)
    return df


# ─────────────────────────────────────────────────────────────────────
# marker 단위 누적기 (bin_extractor._BinAccumulator 와 동일 구조)
# ─────────────────────────────────────────────────────────────────────
class _MarkerAccumulator:
    """
    단일 marker region 의 지표를 누적합니다.
    WPS 는 marker 크기가 작을 수 있으므로
    spanning/endpoints 배열을 marker_len 으로 초기화합니다.
    """
    __slots__ = (
        "marker_id", "marker_type", "chrom", "start", "end", "marker_len",
        "short_count", "long_count",
        "short_cov", "long_cov",
        "s_L_span", "s_L_ep", "s_S_span", "s_S_ep",
        "l_L_span", "l_L_ep", "l_S_span", "l_S_ep",
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

    def add_fragment(self, fs: FragmentScore) -> None:
        if fs.is_short:
            self.short_count += 1
        else:
            self.long_count  += 1

        lo = max(fs.frag_start, self.start) - self.start
        hi = min(fs.frag_end,   self.end)   - self.start
        if lo < hi:
            if fs.is_short:
                self.short_cov[lo:hi] = True
            else:
                self.long_cov[lo:hi]  = True

        self._accum_wps(fs)

    def _accum_wps(self, fs: FragmentScore) -> None:
        bs, be = self.start, self.end
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

            sp_lo = max(fs.frag_start + hk + 1, bs) - bs
            sp_hi = min(fs.frag_end   - hk,     be) - bs
            if sp_lo < sp_hi:
                span_arr[sp_lo:sp_hi] += 1

            for ep in (fs.frag_start, fs.frag_end - 1):
                lo = max(ep - hk,     bs) - bs
                hi = min(ep + hk + 1, be) - bs
                if lo < hi:
                    ep_arr[lo:hi] += 1

    def to_row(self, gc: float, mappability: float) -> dict:
        ml = float(self.marker_len)

        def _wps_median(span, ep) -> float:
            return float(np.median((span - ep).astype(np.float32)))

        total = self.short_count + self.long_count
        return {
            "chrom":        self.chrom,
            "start":        self.start,
            "end":          self.end,
            "marker_id":    self.marker_id,
            "marker_type":  self.marker_type,
            "short_count":  self.short_count,
            "long_count":   self.long_count,
            "total_count":  total,
            "short_ratio":  float(self.short_count / total) if total > 0 else float("nan"),
            "short_wps_L":  _wps_median(self.s_L_span, self.s_L_ep),
            "long_wps_L":   _wps_median(self.l_L_span, self.l_L_ep),
            "short_wps_S":  _wps_median(self.s_S_span, self.s_S_ep),
            "long_wps_S":   _wps_median(self.l_S_span, self.l_S_ep),
            "short_breadth": float(self.short_cov.sum() / ml),
            "long_breadth":  float(self.long_cov.sum()  / ml),
            "gc":            gc,
            "mappability":   mappability,
        }


# ─────────────────────────────────────────────────────────────────────
# 어노테이션 헬퍼
# ─────────────────────────────────────────────────────────────────────
def _gc_of_region(fasta: Optional[pysam.FastaFile],
                  chrom: str, start: int, end: int) -> float:
    if fasta is None:
        return float("nan")
    try:
        seq = fasta.fetch(chrom, start, end).upper()
        return (seq.count("G") + seq.count("C")) / len(seq) if seq else float("nan")
    except Exception:
        return float("nan")


def _mappability_of_region(bw, chrom: str, start: int, end: int) -> float:
    if bw is None:
        return float("nan")
    try:
        val = bw.stats(chrom, start, end, type="mean")[0]
        return float(val) if val is not None else float("nan")
    except Exception:
        return float("nan")


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 스캔 (같은 chrom marker 를 1회 fetch 로 처리)
# ─────────────────────────────────────────────────────────────────────
def _scan_chrom_markers(
    chrom:      str,
    markers_df: pd.DataFrame,   # 해당 chrom 의 marker 행들
    bam_path:   str,
    fasta_path: Optional[str],
    bw_path:    Optional[str],
    min_mapq:   int,
) -> list[dict]:
    """
    같은 염색체의 marker 를 한 번의 BAM fetch 로 처리합니다.

    전략
    ----
    - 모든 marker region 을 커버하는 범위 한 번 fetch
    - 각 read 의 midpoint 가 속한 marker 에 귀속
      (bin_extractor 와 동일: midpoint 기반 단일 귀속)
    - WPS 는 marker 경계 기준으로 spanning/endpoints 누적
    """
    if markers_df.empty:
        return []

    # Accumulator 초기화
    accumulators: dict[str, _MarkerAccumulator] = {}
    for row in markers_df.itertuples(index=False):
        acc = _MarkerAccumulator(
            marker_id   = str(row.marker_id),
            marker_type = str(row.marker_type),
            chrom       = chrom,
            start       = int(row.start),
            end         = int(row.end),
        )
        accumulators[str(row.marker_id)] = acc

    # GC / Mappability 사전 계산
    fasta_h = pysam.FastaFile(fasta_path) if fasta_path else None
    bw_h    = (_pbw.open(bw_path) if (_HAS_BW and bw_path) else None)

    gc_vals  = {str(row.marker_id): _gc_of_region(fasta_h, chrom, int(row.start), int(row.end))
                for row in markers_df.itertuples(index=False)}
    map_vals = {str(row.marker_id): _mappability_of_region(bw_h, chrom, int(row.start), int(row.end))
                for row in markers_df.itertuples(index=False)}

    if fasta_h: fasta_h.close()
    if bw_h:    bw_h.close()

    # bisect 용 정렬된 구조 (midpoint 기반 귀속)
    sorted_rows = markers_df.sort_values("start").reset_index(drop=True)
    marker_starts = sorted_rows["start"].tolist()
    marker_ends   = sorted_rows["end"].tolist()
    marker_ids    = sorted_rows["marker_id"].astype(str).tolist()

    # fetch 범위
    fetch_start = max(0, int(sorted_rows["start"].min()) - 200)
    fetch_end   = int(sorted_rows["end"].max()) + 200

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

            # midpoint 기반 단일 marker 귀속
            mid = fs.midpoint
            idx = bisect.bisect_right(marker_starts, mid) - 1
            if idx < 0 or idx >= len(marker_starts):
                continue
            if not (marker_starts[idx] <= mid < marker_ends[idx]):
                continue

            marker_id = marker_ids[idx]
            if marker_id not in accumulators:
                continue

            accumulators[marker_id].add_fragment(fs)

    return [acc.to_row(gc_vals[mid], map_vals[mid])
            for mid, acc in accumulators.items()]


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker
# ─────────────────────────────────────────────────────────────────────
def _chrom_worker(args: tuple) -> list[dict]:
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
) -> pd.DataFrame:
    """
    marker BED → marker_stats.parquet

    Parameters
    ----------
    bam_path   : 인덱스된 BAM
    marker_bed : marker BED (chrom start end marker_id marker_type)
    out_path   : 저장 경로
    fasta_path : GC 계산용 FASTA (선택)
    bw_path    : Mappability bigWig (선택)
    """
    marker_df = load_marker_bed(marker_bed)
    if marker_df.empty:
        log.warning("marker 없음")
        return pd.DataFrame()

    # 염색체별 그룹화
    chrom_groups = {c: g.reset_index(drop=True)
                    for c, g in marker_df.groupby("chrom")}
    log.info("marker 스캔: %d markers, %d chroms",
             len(marker_df), len(chrom_groups))

    all_rows: list[dict] = []

    task_args = [
        (chrom, cdf, bam_path, fasta_path, bw_path, min_mapq)
        for chrom, cdf in chrom_groups.items()
    ]

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_chrom_worker, a): a[0] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                rows = fut.result()
                all_rows.extend(rows)
                log.info("  ✓ %s (%d markers)", chrom, len(rows))
            except Exception as e:
                log.error("  ✗ %s: %s", chrom, e)

    df = pd.DataFrame(all_rows)
    if df.empty:
        log.warning("결과 없음")
        return df

    # dtype 정리
    for c in ["short_count", "long_count", "total_count"]:
        if c in df: df[c] = df[c].astype("int32")
    for c in ["short_ratio", "short_wps_L", "long_wps_L",
              "short_wps_S", "long_wps_S",
              "short_breadth", "long_breadth", "gc", "mappability"]:
        if c in df: df[c] = df[c].astype("float32")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("marker_stats 저장: %s (%d rows)", out_path, len(df))
    return df