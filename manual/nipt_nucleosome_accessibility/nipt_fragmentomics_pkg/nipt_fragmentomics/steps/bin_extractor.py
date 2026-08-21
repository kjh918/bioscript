"""
bin_extractor.py
================
BAM 파일을 염색체 단위로 1회 순차 스캔하여
bin 단위 원시 지표를 수집합니다.

수집 지표 (bin 당)
------------------
  short / long count     — midpoint 귀속
  short / long breadth   — coverage 비율
  short / long wps_L/S   — WPS median (position별 spanning-endpoints 배열의 median)
  gc                     — FASTA 참조 GC 비율
  mappability            — bigWig 평균

WPS 계산 원리
-------------
  WPS(i) = N_spanning(i) - N_endpoints(i)

  N_spanning(i) : position i 를 중심으로 한 window [i-k/2, i+k/2] 를
                  완전히 가로지르는 fragment 수
                  → fragment_start < i-k/2  AND  fragment_end > i+k/2

  N_endpoints(i): 해당 window 내에 5' 또는 3' 말단이 있는 fragment 수

  bin 내 모든 position 에 대해 WPS(i) 를 계산하고
  bin 대표값 = median(WPS 배열)

  중요: fragment 가 bin 과 overlap 되면 해당 bin 의 WPS 배열에 기여함
        (midpoint 귀속이 아님 — bin 경계를 넘는 fragment 도 포함)
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

from nipt_fragmentomics.core.constants import DEFAULT_BIN_SIZE, MIN_MAPQ, WPS_PARAMS, SHORT_MAX

log = logging.getLogger(__name__)

try:
    import pyBigWig as _pbw
    _HAS_BW = True
except ImportError:
    _HAS_BW = False


# ─────────────────────────────────────────────────────────────────────
# bin 단위 누적기
# ─────────────────────────────────────────────────────────────────────
class _BinAccumulator:
    """
    단일 bin 의 count / breadth / WPS 지표를 누적합니다.

    count  : midpoint 기반 (단일 bin 귀속)
    breadth: overlap 기반
    WPS    : position별 배열 누적 (overlap 기반)
              → bin 내 각 bp position 마다 spanning/endpoints 카운트
              → 최종: median(spanning - endpoints)
    """
    __slots__ = (
        "chrom", "start", "end", "bin_len",
        "short_count", "long_count",
        "short_cov", "long_cov",
        # position별 WPS 배열 (bin_len 크기)
        "s_L_span", "s_L_ep",
        "s_S_span", "s_S_ep",
        "l_L_span", "l_L_ep",
        "l_S_span", "l_S_ep",
    )

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom   = chrom
        self.start   = start
        self.end     = end
        self.bin_len = end - start

        self.short_count = 0
        self.long_count  = 0
        self.short_cov   = np.zeros(self.bin_len, dtype=bool)
        self.long_cov    = np.zeros(self.bin_len, dtype=bool)

        for attr in ("s_L_span", "s_L_ep", "s_S_span", "s_S_ep",
                     "l_L_span", "l_L_ep", "l_S_span", "l_S_ep"):
            setattr(self, attr, np.zeros(self.bin_len, dtype=np.int32))

    def add_count_and_breadth(self, frag_start: int, frag_end: int,
                               is_short: bool, midpoint: int) -> None:
        """count(midpoint 기반) + breadth(overlap 기반) 누적."""
        # count: midpoint 가 이 bin 에 있을 때만
        bs, be = self.start, self.end
        if bs <= midpoint < be:
            if is_short: self.short_count += 1
            else:        self.long_count  += 1

        # breadth: overlap
        lo = max(frag_start, bs) - bs
        hi = min(frag_end,   be) - bs
        if lo < hi:
            if is_short: self.short_cov[lo:hi] = True
            else:        self.long_cov[lo:hi]  = True

    def add_wps(self, frag_start: int, frag_end: int, frag_len: int,
                is_short: bool) -> None:
        """
        WPS position 배열에 누적.
        fragment 가 bin 과 overlap 될 때 호출됩니다.

        spanning(i): fragment_start < i - k/2  AND  fragment_end > i + k/2
        endpoints(i): 5' 또는 3' 말단이 window [i-k/2, i+k/2] 내에 있음
        """
        bs, be = self.start, self.end

        for mode, prm in WPS_PARAMS.items():
            if not (prm["frag_min"] <= frag_len <= prm["frag_max"]):
                continue
            hk = prm["window"] // 2

            if is_short:
                span_arr = getattr(self, f"s_{mode}_span")
                ep_arr   = getattr(self, f"s_{mode}_ep")
            else:
                span_arr = getattr(self, f"l_{mode}_span")
                ep_arr   = getattr(self, f"l_{mode}_ep")

            # spanning: position i 에서 fragment 가 window 전체를 덮음
            # 조건: frag_start < i - hk  AND  frag_end > i + hk
            # → i 범위: (frag_start + hk, frag_end - hk) exclusive
            sp_lo = max(frag_start + hk + 1, bs) - bs
            sp_hi = min(frag_end   - hk,     be) - bs
            if sp_lo < sp_hi:
                span_arr[sp_lo:sp_hi] += 1

            # endpoints: 5' 말단(frag_start), 3' 말단(frag_end-1) 각각
            # position i 에서 말단이 window [i-hk, i+hk] 안에 있으면
            # → 말단 ep 에 대해, i 범위: [ep-hk, ep+hk]
            for ep in (frag_start, frag_end - 1):
                w_lo = max(ep - hk,     bs) - bs
                w_hi = min(ep + hk + 1, be) - bs
                if w_lo < w_hi:
                    ep_arr[w_lo:w_hi] += 1

    def to_row(self, gc: float, mappability: float) -> dict:
        bl = float(self.bin_len)

        def _wps_med(span, ep) -> float:
            return float(np.median((span - ep).astype(np.float32)))

        return {
            "chrom":         self.chrom,
            "start":         self.start,
            "end":           self.end,
            "short_count":   self.short_count,
            "long_count":    self.long_count,
            "short_breadth": float(self.short_cov.sum() / bl),
            "long_breadth":  float(self.long_cov.sum()  / bl),
            "short_wps_L":   _wps_med(self.s_L_span, self.s_L_ep),
            "long_wps_L":    _wps_med(self.l_L_span, self.l_L_ep),
            "short_wps_S":   _wps_med(self.s_S_span, self.s_S_ep),
            "long_wps_S":    _wps_med(self.l_S_span, self.l_S_ep),
            "gc":            gc,
            "mappability":   mappability,
        }


# ─────────────────────────────────────────────────────────────────────
# 어노테이션 헬퍼
# ─────────────────────────────────────────────────────────────────────
def _gc_of_region(fasta, chrom, start, end) -> float:
    if fasta is None: return float("nan")
    try:
        seq = fasta.fetch(chrom, start, end).upper()
        return (seq.count("G") + seq.count("C")) / len(seq) if seq else float("nan")
    except Exception: return float("nan")


def _map_of_region(bw, chrom, start, end) -> float:
    if bw is None: return float("nan")
    try:
        val = bw.stats(chrom, start, end, type="mean")[0]
        return float(val) if val is not None else float("nan")
    except Exception: return float("nan")


# ─────────────────────────────────────────────────────────────────────
# BED 로더
# ─────────────────────────────────────────────────────────────────────
def _read_bed(bed_path: str) -> pd.DataFrame:
    import gzip
    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt") as f:
        first = f.readline().rstrip("\n")
    try:
        int(first.split("\t")[1])
        has_header = False
    except (ValueError, IndexError):
        has_header = True

    bed = pd.read_csv(
        bed_path, sep="\t",
        header=0 if has_header else None,
        comment="#",
        compression="gzip" if bed_path.endswith(".gz") else None,
    )
    if not has_header:
        cols = list(bed.columns)
        cols[0], cols[1], cols[2] = "chrom", "start", "end"
        bed.columns = cols
    else:
        bed = bed.rename(columns={
            bed.columns[0]: "chrom",
            bed.columns[1]: "start",
            bed.columns[2]: "end",
        })

    bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
    bed["end"]   = pd.to_numeric(bed["end"],   errors="coerce")
    bed = bed.dropna(subset=["start", "end"]).copy()
    bed["start"] = bed["start"].astype(int)
    bed["end"]   = bed["end"].astype(int)
    log.info("BED 로드: %d bins [%s]", len(bed), bed_path)
    return bed


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 스캔
# ─────────────────────────────────────────────────────────────────────
def _scan_chrom(
    chrom:       str,
    bins:        list[tuple[int, int]],
    bam_path:    str,
    fasta_path:  Optional[str],
    bw_path:     Optional[str],
    min_mapq:    int,
    gc_precomp:  Optional[list[float]] = None,
    map_precomp: Optional[list[float]] = None,
) -> list[dict]:
    if not bins: return []

    bin_starts   = [s for s, _ in bins]
    bin_ends     = [e for _, e in bins]
    accumulators = [_BinAccumulator(chrom, s, e) for s, e in bins]

    # GC / Mappability
    if gc_precomp is not None:
        gc_vals = gc_precomp
    else:
        fasta_h = pysam.FastaFile(fasta_path) if fasta_path else None
        gc_vals = [_gc_of_region(fasta_h, chrom, s, e) for s, e in bins]
        if fasta_h: fasta_h.close()

    if map_precomp is not None:
        map_vals = map_precomp
    else:
        bw_h = (_pbw.open(bw_path) if (_HAS_BW and bw_path) else None)
        map_vals = [_map_of_region(bw_h, chrom, s, e) for s, e in bins]
        if bw_h: bw_h.close()

    # WPS 최대 window 패딩 (fetch 범위 확장용)
    max_hk = max(p["window"] // 2 for p in WPS_PARAMS.values())
    max_frag = max(p["frag_max"] for p in WPS_PARAMS.values())
    fetch_pad = max_frag + max_hk

    chrom_start = max(0, bin_starts[0] - fetch_pad)
    chrom_end   = bin_ends[-1] + fetch_pad

    seen: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
        for read in bam.fetch(chrom, chrom_start, chrom_end):
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

            if read.is_paired and read.template_length != 0:
                tlen = abs(read.template_length)
                # 비정상 template_length 필터 (수백만 bp 등)
                if tlen > 1000 or tlen < (read.query_length or 0):
                    frag_start = read.reference_start
                    frag_end   = read.reference_end
                else:
                    frag_start = read.reference_start
                    frag_end   = frag_start + tlen
            else:
                frag_start = read.reference_start
                frag_end   = read.reference_end

            frag_len = frag_end - frag_start
            if frag_len <= 0 or frag_len > 1000:
                continue

            is_short = (frag_len <= SHORT_MAX)
            midpoint = (frag_start + frag_end) // 2

            # ── count + breadth: midpoint 기반 단일 bin ──────────
            cnt_idx = bisect.bisect_right(bin_starts, midpoint) - 1
            if 0 <= cnt_idx < len(bins):
                s, e = bins[cnt_idx]
                if s <= midpoint < e:
                    accumulators[cnt_idx].add_count_and_breadth(
                        frag_start, frag_end, is_short, midpoint)

            # ── WPS: fragment 와 overlap 되는 모든 bin ───────────
            # WPS 계산에 필요한 fetch 범위: fragment ±max_hk
            # overlap 조건: bin_start < frag_end+hk AND bin_end > frag_start-hk
            max_hk_wps = max((WPS_PARAMS[m]["window"] // 2
                              for m in WPS_PARAMS
                              if WPS_PARAMS[m]["frag_min"] <= frag_len <= WPS_PARAMS[m]["frag_max"]),
                             default=0)
            if max_hk_wps == 0:
                continue   # 이 fragment 는 어떤 WPS 모드에도 해당 안 됨

            eff_start = frag_start - max_hk_wps
            eff_end   = frag_end   + max_hk_wps

            lo_idx = bisect.bisect_left(bin_ends,    eff_start + 1)
            hi_idx = bisect.bisect_right(bin_starts, eff_end   - 1)

            for idx in range(lo_idx, hi_idx):
                accumulators[idx].add_wps(frag_start, frag_end, frag_len, is_short)

    return [acc.to_row(gc_vals[i], map_vals[i])
            for i, acc in enumerate(accumulators)]


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker
# ─────────────────────────────────────────────────────────────────────
def _scan_chrom_worker(args: tuple) -> list[dict]:
    chrom, bins, bam_path, fasta_path, bw_path, min_mapq, gc_pre, map_pre = args
    return _scan_chrom(chrom, bins, bam_path, fasta_path, bw_path,
                       min_mapq, gc_pre, map_pre)


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:   str,
    out_path:   str,
    bed_path:   Optional[str] = None,
    fasta_path: Optional[str] = None,
    bw_path:    Optional[str] = None,
    bin_size:   int           = DEFAULT_BIN_SIZE,
    min_mapq:   int           = MIN_MAPQ,
    n_jobs:     int           = 4,
    vcf_path:   Optional[str] = None,
    min_baseq:  int           = 20,
) -> pd.DataFrame:
    """BAM → bins_raw.parquet"""

    if bed_path and os.path.exists(bed_path):
        bed = _read_bed(bed_path)
    else:
        log.info("BED 없음 — 자동 grid (bin_size=%d)", bin_size)
        rows = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for sq in bam.header["SQ"]:
                chrom, length = sq["SN"], sq["LN"]
                for s in range(0, length, bin_size):
                    rows.append({"chrom": chrom, "start": s,
                                 "end": min(s + bin_size, length)})
        bed = pd.DataFrame(rows)

    bed = bed.sort_values(["chrom", "start"]).reset_index(drop=True)

    # --fasta / --bw 있으면 BED 컬럼 무시하고 재계산
    has_gc  = ("gc"          in bed.columns) and (fasta_path is None)
    has_map = ("mappability" in bed.columns) and (bw_path    is None)
    if not has_gc  and fasta_path: log.info("FASTA로 gc 재계산")
    if not has_map and bw_path:    log.info("bigWig로 mappability 재계산")

    chrom_data: dict[str, dict] = {}
    for chrom, g in bed.groupby("chrom"):
        chrom_data[chrom] = {
            "bins": list(zip(g["start"], g["end"])),
            "gc":   g["gc"].tolist()          if has_gc  else None,
            "map":  g["mappability"].tolist() if has_map else None,
        }

    log.info("총 %d bins, %d chroms", len(bed), len(chrom_data))

    task_args = [
        (chrom, chrom_data[chrom]["bins"], bam_path,
         None if has_gc  else fasta_path,
         None if has_map else bw_path,
         min_mapq,
         chrom_data[chrom]["gc"],
         chrom_data[chrom]["map"])
        for chrom in chrom_data
    ]

    all_rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {ex.submit(_scan_chrom_worker, a): a[0] for a in task_args}
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                rows = fut.result()
                all_rows.extend(rows)
                log.info("  ✓ %s (%d bins)", chrom, len(rows))
            except Exception as exc:
                log.error("  ✗ %s: %s", chrom, exc)

    df = pd.DataFrame(all_rows)

    # BED 필터 컬럼 보존
    for col in ("is_filtered", "is_low_mappability",
                "is_blacklisted", "is_abnormal_gc"):
        if col in bed.columns:
            bed_key = bed.set_index(["chrom", "start", "end"])[col]
            df[col] = df.set_index(["chrom", "start", "end"]).index.map(
                lambda k: bed_key.get(k, False)
            )
            df[col] = df[col].fillna(False).astype(bool)

    # dtype
    for c in ["short_count", "long_count"]:
        if c in df: df[c] = df[c].astype("int32")
    for c in ["short_breadth", "long_breadth",
              "short_wps_L", "long_wps_L",
              "short_wps_S", "long_wps_S",
              "gc", "mappability"]:
        if c in df: df[c] = df[c].astype("float32")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("bins_raw 저장: %s (%d rows)", out_path, len(df))
    return df