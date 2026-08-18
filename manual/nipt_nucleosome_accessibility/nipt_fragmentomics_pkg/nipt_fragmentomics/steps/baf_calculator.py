"""
baf_calculator.py
=================
BAM + population SNP VCF 에서 bin 단위 BAF (B-Allele Frequency) 를 계산합니다.

처리 순서
---------
1. VCF 에서 population heterozygous site 추출
   (AF 기준: hetero_min ≤ pop_af ≤ hetero_max)
2. BAM pileup 으로 각 site 의 ref/alt depth 수집
   - short (≤150 bp) / long (≥151 bp) fragment 분리
3. bin 단위 집계
   - short/long/combined 각각 BAF median, std, n_sites

BAF 해석
--------
  정상 heterozygous site: BAF ≈ 0.5
  이수성 (gain): alt 과잉 → BAF > 0.5 또는 < 0.5 편향 (imbalance)
  태아 fraction 이 낮으면 편향 신호도 약함

출력 (bins_baf.parquet)
-----------------------
  chrom, start, end,
  short_baf_median, short_baf_std, short_n_sites,
  long_baf_median,  long_baf_std,  long_n_sites,
  combined_baf_median, combined_baf_std, combined_n_sites,
  short_baf_mad,    long_baf_mad,   combined_baf_mad
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

log = logging.getLogger(__name__)

SHORT_MAX = 150

# Population AF 기준 — heterozygous site 선별
POP_AF_MIN = 0.2
POP_AF_MAX = 0.8

# pileup QC
MIN_MAPQ  = 20
MIN_BASEQ = 20
MIN_SITE_DEPTH = 5   # site 당 최소 depth (미만이면 해당 site 제외)


# ─────────────────────────────────────────────────────────────────────
# VCF 로더: chrom 전체 heterozygous site 추출
# ─────────────────────────────────────────────────────────────────────
def _load_vcf_sites(
    vcf_path: str,
    chrom:    str,
    start:    int,
    end:      int,
    af_min:   float = POP_AF_MIN,
    af_max:   float = POP_AF_MAX,
) -> dict[int, tuple[str, str]]:
    """
    VCF 에서 pos → (ref, alt) 를 반환합니다.
    population AF 가 [af_min, af_max] 범위인 SNP 만 포함합니다.

    Returns
    -------
    {pos(1-based): (ref, alt)}
    """
    sites: dict[int, tuple[str, str]] = {}
    try:
        with pysam.VariantFile(vcf_path) as vcf:
            for rec in vcf.fetch(chrom, start, end):
                # SNP 만 (ref/alt 모두 단일 염기)
                if len(rec.ref) != 1:
                    continue
                if not rec.alts or len(rec.alts[0]) != 1:
                    continue

                kova   = rec.info.get("KOVA_AF")
                gnomad = rec.info.get("GNOMAD_AF")
                kaf = (kova[0]   if isinstance(kova,   (list,tuple)) else kova)   or 0.0
                gaf = (gnomad[0] if isinstance(gnomad, (list,tuple)) else gnomad) or 0.0
                pop_af = kaf if kaf > 0 else gaf

                if af_min <= pop_af <= af_max:
                    sites[rec.pos] = (rec.ref, rec.alts[0])
    except Exception as e:
        log.warning("VCF fetch 실패 [%s:%d-%d]: %s", chrom, start, end, e)
    return sites


# ─────────────────────────────────────────────────────────────────────
# 염색체 단위 BAF 계산
# ─────────────────────────────────────────────────────────────────────
def _calc_chrom_baf(
    chrom:    str,
    bins:     list[tuple[int, int]],
    bam_path: str,
    vcf_path: str,
    min_mapq:  int = MIN_MAPQ,
    min_baseq: int = MIN_BASEQ,
    min_depth: int = MIN_SITE_DEPTH,
    af_min:    float = POP_AF_MIN,
    af_max:    float = POP_AF_MAX,
) -> list[dict]:
    """
    chrom 의 bins 에 대해 BAF 를 계산합니다.
    BAM 1회 스캔으로 모든 bin 처리.
    """
    if not bins:
        return []

    bin_starts = [s for s, _ in bins]
    chrom_start = bins[0][0]
    chrom_end   = bins[-1][1]

    # VCF site 로드 (chrom 전체)
    site_lookup = _load_vcf_sites(vcf_path, chrom,
                                  chrom_start, chrom_end,
                                  af_min, af_max)
    if not site_lookup:
        log.debug("[BAF] %s: VCF site 없음", chrom)
        return [_empty_row(chrom, s, e) for s, e in bins]

    # site → bin index 매핑
    pos_to_bin: dict[int, int] = {}
    for pos in site_lookup:
        idx = bisect.bisect_right(bin_starts, pos - 1) - 1   # 0-based pos 보정
        if 0 <= idx < len(bins):
            s, e = bins[idx]
            if s < pos <= e:   # VCF pos는 1-based, bin은 0-based half-open
                pos_to_bin[pos] = idx

    # bin 별 site 누적기: {idx: {pos: {"sr","sa","lr","la"}}}
    # sr=short_ref, sa=short_alt, lr=long_ref, la=long_alt
    bin_sites: dict[int, dict[int, dict]] = {
        i: {} for i in range(len(bins))
    }
    for pos, bin_idx in pos_to_bin.items():
        bin_sites[bin_idx][pos] = {"sr": 0, "sa": 0, "lr": 0, "la": 0}

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

            # fragment size → short/long
            if read.is_paired and read.template_length != 0:
                flen = abs(read.template_length)
            else:
                flen = read.query_length or 0
            is_short = (flen <= SHORT_MAX)

            # aligned pair 순회
            if read.query_qualities is None or read.query_sequence is None:
                continue

            for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
                g_pos = r_pos + 1   # 1-based
                if g_pos not in pos_to_bin:
                    continue
                if read.query_qualities[q_pos] < min_baseq:
                    continue

                base    = read.query_sequence[q_pos].upper()
                ref, alt = site_lookup[g_pos]
                bin_idx  = pos_to_bin[g_pos]

                if bin_idx not in bin_sites or g_pos not in bin_sites[bin_idx]:
                    continue

                stat = bin_sites[bin_idx][g_pos]
                if is_short:
                    if   base == alt.upper(): stat["sa"] += 1
                    elif base == ref.upper(): stat["sr"] += 1
                else:
                    if   base == alt.upper(): stat["la"] += 1
                    elif base == ref.upper(): stat["lr"] += 1

    # bin 단위 집계
    rows = []
    for i, (s, e) in enumerate(bins):
        site_data = bin_sites[i]

        short_bafs, long_bafs, comb_bafs = [], [], []

        for pos, stat in site_data.items():
            s_total = stat["sa"] + stat["sr"]
            l_total = stat["la"] + stat["lr"]
            c_total = s_total + l_total

            if s_total >= min_depth:
                short_bafs.append(stat["sa"] / s_total)
            if l_total >= min_depth:
                long_bafs.append(stat["la"] / l_total)
            if c_total >= min_depth:
                comb_bafs.append((stat["sa"] + stat["la"]) / c_total)

        rows.append(_summarize_bafs(chrom, s, e, short_bafs, long_bafs, comb_bafs))

    return rows


def _summarize_bafs(
    chrom: str, start: int, end: int,
    short_bafs: list[float],
    long_bafs:  list[float],
    comb_bafs:  list[float],
) -> dict:
    """BAF 리스트 → bin row dict."""
    def _stats(bafs: list[float], prefix: str) -> dict:
        if not bafs:
            return {
                f"{prefix}_baf_median": np.nan,
                f"{prefix}_baf_std":    np.nan,
                f"{prefix}_baf_mad":    np.nan,
                f"{prefix}_n_sites":    0,
            }
        arr = np.array(bafs, dtype=np.float32)
        med = float(np.median(arr))
        return {
            f"{prefix}_baf_median": med,
            f"{prefix}_baf_std":    float(arr.std()),
            f"{prefix}_baf_mad":    float(np.median(np.abs(arr - med))),
            f"{prefix}_n_sites":    len(arr),
        }

    row = {"chrom": chrom, "start": start, "end": end}
    row.update(_stats(short_bafs, "short"))
    row.update(_stats(long_bafs,  "long"))
    row.update(_stats(comb_bafs,  "combined"))
    return row


def _empty_row(chrom: str, start: int, end: int) -> dict:
    row = {"chrom": chrom, "start": start, "end": end}
    for prefix in ("short", "long", "combined"):
        row.update({
            f"{prefix}_baf_median": np.nan,
            f"{prefix}_baf_std":    np.nan,
            f"{prefix}_baf_mad":    np.nan,
            f"{prefix}_n_sites":    0,
        })
    return row


# ─────────────────────────────────────────────────────────────────────
# 모듈 레벨 worker (ProcessPoolExecutor pickle 호환)
# ─────────────────────────────────────────────────────────────────────
def _baf_chrom_worker(chrom, bins, bam_path, vcf_path,
                      min_mapq, min_baseq, min_depth, af_min, af_max):
    return _calc_chrom_baf(chrom, bins, bam_path, vcf_path,
                           min_mapq, min_baseq, min_depth, af_min, af_max)


# ─────────────────────────────────────────────────────────────────────
# 공개 인터페이스
# ─────────────────────────────────────────────────────────────────────
def run(
    bam_path:   str,
    vcf_path:   str,
    bin_path:   str,   # bins_raw.parquet 또는 bins_corrected.parquet (bin 좌표 참조)
    out_path:   str,
    min_mapq:   int   = MIN_MAPQ,
    min_baseq:  int   = MIN_BASEQ,
    min_depth:  int   = MIN_SITE_DEPTH,
    af_min:     float = POP_AF_MIN,
    af_max:     float = POP_AF_MAX,
    n_jobs:     int   = 4,
) -> pd.DataFrame:
    """
    BAM + VCF + bin 좌표 → bins_baf.parquet

    Parameters
    ----------
    bin_path : bins_raw.parquet 또는 bins_corrected.parquet
               chrom/start/end 컬럼만 사용 (bin 좌표 참조용)
    """
    if not os.path.exists(vcf_path):
        log.error("VCF 파일 없음: %s", vcf_path)
        return pd.DataFrame()

    bin_df = pd.read_parquet(bin_path, columns=["chrom", "start", "end"])
    chrom_groups = {
        c: list(zip(g["start"], g["end"]))
        for c, g in bin_df.groupby("chrom")
    }
    log.info("BAF 계산: %d bins, %d chroms", len(bin_df), len(chrom_groups))

    all_rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futures = {
            ex.submit(_baf_chrom_worker, c, b,
                      bam_path, vcf_path,
                      min_mapq, min_baseq, min_depth, af_min, af_max): c
            for c, b in chrom_groups.items()
        }
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                rows = fut.result()
                all_rows.extend(rows)
                n_valid = sum(1 for r in rows if r["combined_n_sites"] > 0)
                log.info("  ✓ %s  (%d bins, %d with sites)", chrom, len(rows), n_valid)
            except Exception as exc:
                log.error("  ✗ %s: %s", chrom, exc)

    df = pd.DataFrame(all_rows)
    if df.empty:
        log.warning("BAF 결과 없음")
        return df

    # dtype 최적화
    for col in df.columns:
        if "_n_sites" in col:
            df[col] = df[col].astype("int32")
        elif col not in ("chrom",):
            try:
                df[col] = df[col].astype("float32")
            except Exception:
                pass

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    df.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("bins_baf 저장: %s (%d rows)", out_path, len(df))
    return df


def merge_into_cnv(
    cnv_path: str,
    baf_path: str,
    out_path: str,
) -> pd.DataFrame:
    """
    cnv_calls.parquet 에 bins_baf.parquet 를 join 하여
    BAF 컬럼이 추가된 최종 결과를 저장합니다.
    """
    cnv_df = pd.read_parquet(cnv_path)
    baf_df = pd.read_parquet(baf_path)

    baf_cols = ["chrom", "start", "end"] + [
        c for c in baf_df.columns if c not in ("chrom", "start", "end")
    ]
    merged = cnv_df.merge(
        baf_df[baf_cols],
        on=["chrom", "start", "end"],
        how="left",
    )
    merged.to_parquet(out_path, engine="pyarrow", index=False)
    log.info("CNV+BAF 병합 저장: %s (%d rows)", out_path, len(merged))
    return merged
