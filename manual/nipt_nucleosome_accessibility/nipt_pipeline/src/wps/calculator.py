"""
calculator.py - WPS 핵심 계산 모듈

Short-range WPS (120 bp window)와 Long-range WPS (1000 bp window)를
BAM 파일에서 직접 계산한다.

알고리즘 (Snyder et al. 2016):
    WPS(t, w) = (윈도우 내 완전 포함 fragment 수) - (윈도우 경계 걸침 fragment 수)
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pysam

from preprocess import PreprocessResult
from utils import get_logger, log_progress

# ── 상수 ──────────────────────────────────────────────────────────────────────

SHORT_WINDOW: int = 120    # short-range WPS 윈도우 크기 (bp)
LONG_WINDOW: int = 1000   # long-range WPS 윈도우 크기 (bp)
HALF_SHORT: int = SHORT_WINDOW // 2
HALF_LONG: int = LONG_WINDOW // 2


# ── 데이터 클래스 ─────────────────────────────────────────────────────────────

@dataclass
class ChromWPS:
    """단일 염색체의 WPS 결과.

    Attributes:
        chrom: 염색체 이름
        length: 염색체 길이
        short_wps: short-range WPS 배열 (length 크기)
        long_wps: long-range WPS 배열 (length 크기)
        n_fragments: 처리된 fragment 수
    """
    chrom: str
    length: int
    short_wps: np.ndarray
    long_wps: np.ndarray
    n_fragments: int


@dataclass
class WPSResult:
    """전체 게놈 WPS 계산 결과.

    Attributes:
        chroms: 염색체별 WPS 결과 딕셔너리
        total_fragments: 전체 처리된 fragment 수
    """
    chroms: dict[str, ChromWPS] = field(default_factory=dict)
    total_fragments: int = 0


# ── 공개 함수 ─────────────────────────────────────────────────────────────────

def run_calculate(
    prep: PreprocessResult,
    logger: Optional[logging.Logger] = None,
) -> WPSResult:
    """전체 분석 대상 염색체에 대해 WPS를 계산한다.

    Args:
        prep: 전처리 결과 (PreprocessResult)
        logger: Logger 인스턴스

    Returns:
        WPSResult: 염색체별 WPS 배열을 담은 결과 객체
    """
    log = logger or get_logger(__name__)
    result = WPSResult()

    with pysam.AlignmentFile(str(prep.bam_path), "rb") as bam:
        for chrom in prep.chroms:
            chrom_len = prep.chrom_lengths[chrom]
            log.info(f"WPS 계산 시작: {chrom} ({chrom_len:,} bp)")

            chrom_wps = _calculate_chrom_wps(
                bam=bam,
                chrom=chrom,
                chrom_len=chrom_len,
                min_frag_len=prep.min_frag_len,
                max_frag_len=prep.max_frag_len,
                min_mapq=prep.min_mapq,
                logger=log,
            )
            result.chroms[chrom] = chrom_wps
            result.total_fragments += chrom_wps.n_fragments
            log.info(
                f"{chrom} 완료: {chrom_wps.n_fragments:,} fragments 처리"
            )

    log.info(f"전체 처리 fragment 수: {result.total_fragments:,}")
    return result


# ── 내부 함수 ─────────────────────────────────────────────────────────────────

def _calculate_chrom_wps(
    bam: pysam.AlignmentFile,
    chrom: str,
    chrom_len: int,
    min_frag_len: int,
    max_frag_len: int,
    min_mapq: int,
    logger: logging.Logger,
) -> ChromWPS:
    """단일 염색체에 대한 short/long WPS를 계산한다.

    각 염색체 포지션 t에서:
    - complete[t]: fragment가 [t-half, t+half] 윈도우에 완전히 포함되면 +1
    - spanning[t]: fragment가 윈도우 경계를 걸치면 +1

    WPS = complete - spanning

    Args:
        bam: 열린 pysam AlignmentFile
        chrom: 염색체 이름
        chrom_len: 염색체 길이
        min_frag_len: fragment 최소 길이
        max_frag_len: fragment 최대 길이
        min_mapq: 최소 MapQ

    Returns:
        ChromWPS: 해당 염색체의 WPS 결과
    """
    # complete/spanning 카운터 (short/long 각각)
    s_complete = np.zeros(chrom_len, dtype=np.int32)
    s_spanning = np.zeros(chrom_len, dtype=np.int32)
    l_complete = np.zeros(chrom_len, dtype=np.int32)
    l_spanning = np.zeros(chrom_len, dtype=np.int32)

    n_fragments = 0

    for read in bam.fetch(chrom):
        if not _is_valid_read(read, min_mapq):
            continue

        frag_start, frag_end = _get_fragment_coords(read)
        frag_len = frag_end - frag_start

        if not (min_frag_len <= frag_len <= max_frag_len):
            continue

        n_fragments += 1

        # Short-range WPS 기여
        _accumulate_wps(
            frag_start, frag_end, chrom_len,
            HALF_SHORT, s_complete, s_spanning,
        )
        # Long-range WPS 기여
        _accumulate_wps(
            frag_start, frag_end, chrom_len,
            HALF_LONG, l_complete, l_spanning,
        )

        log_progress(n_fragments, 0, chrom, logger, interval=500_000)

    short_wps = (s_complete - s_spanning).astype(np.float32)
    long_wps = (l_complete - l_spanning).astype(np.float32)

    return ChromWPS(
        chrom=chrom,
        length=chrom_len,
        short_wps=short_wps,
        long_wps=long_wps,
        n_fragments=n_fragments,
    )


def _accumulate_wps(
    frag_start: int,
    frag_end: int,
    chrom_len: int,
    half_window: int,
    complete: np.ndarray,
    spanning: np.ndarray,
) -> None:
    """하나의 fragment가 각 포지션의 complete/spanning 카운터에 기여하도록 누적한다.

    Fragment [frag_start, frag_end)가 윈도우 [t-half, t+half]에 대해:
    - 완전 포함 조건: frag_start >= t-half AND frag_end <= t+half
    - 경계 걸침 조건: 포함도 외부도 아닌 경우

    벡터화를 위해 각 조건을 구간 범위로 변환하여 누적 덧셈.

    Args:
        frag_start: fragment 시작 포지션 (0-based)
        frag_end: fragment 끝 포지션 (exclusive)
        chrom_len: 염색체 길이
        half_window: 윈도우 절반 크기
        complete: 완전 포함 카운터 배열 (in-place 수정)
        spanning: 경계 걸침 카운터 배열 (in-place 수정)
    """
    # 완전 포함: frag_start >= t-half AND frag_end <= t+half
    # → t >= frag_start - half AND t <= frag_end - 1 + half
    # → t ∈ [frag_start - half, frag_end - 1 + half]
    # 단, t ∈ [half, chrom_len - half) 범위 내
    c_lo = max(frag_start - half_window, half_window)
    c_hi = min(frag_end - 1 + half_window, chrom_len - half_window)
    if c_lo <= c_hi:
        complete[c_lo: c_hi + 1] += 1

    # 경계 걸침 (왼쪽): frag_start이 윈도우 내에 있고 frag_end는 밖
    # → t-half <= frag_start < t+half AND frag_end > t+half
    # → t ∈ (frag_end - 1 - half, frag_start + half]
    sp_lo = max(frag_start - half_window + 1, half_window)
    sp_hi = min(frag_end - 1 + half_window - 1, chrom_len - half_window)
    if sp_lo <= sp_hi:
        # 위 구간에서 complete 구간을 제외한 부분이 spanning
        # (complete 이미 처리했으므로 전체 overlap 구간에서 complete 제외)
        # 간단화: complete 카운터와 달리 spanning은 fragment 끝/시작이
        # 윈도우 경계와 교차하는 t를 카운트
        # left-spanning: frag_end > t + half → t < frag_end - half
        ls_hi = min(frag_end - half_window - 1, chrom_len - half_window)
        ls_lo = max(frag_start - half_window + 1, half_window)
        if ls_lo <= ls_hi:
            spanning[ls_lo: ls_hi + 1] += 1

        # right-spanning: frag_start < t - half → t > frag_start + half
        rs_lo = max(frag_start + half_window + 1, half_window)
        rs_hi = min(frag_end - 1 + half_window - 1, chrom_len - half_window)
        if rs_lo <= rs_hi:
            spanning[rs_lo: rs_hi + 1] += 1


def _is_valid_read(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    """Read가 WPS 계산에 유효한지 확인한다.

    Args:
        read: pysam AlignedSegment
        min_mapq: 최소 매핑 품질

    Returns:
        유효하면 True
    """
    return (
        read.is_paired
        and read.is_read1           # pair 중 read1만 처리 (fragment 중복 방지)
        and not read.is_unmapped
        and not read.mate_is_unmapped
        and not read.is_duplicate
        and not read.is_secondary
        and not read.is_supplementary
        and read.mapping_quality >= min_mapq
        and read.is_proper_pair
    )


def _get_fragment_coords(read: pysam.AlignedSegment) -> tuple[int, int]:
    """Read에서 fragment 시작/끝 좌표를 추출한다.

    template_length의 부호로 방향을 판단하여 항상 (start < end) 반환.

    Args:
        read: 유효한 pysam AlignedSegment (read1)

    Returns:
        (frag_start, frag_end): 0-based fragment 좌표
    """
    tlen = read.template_length
    if tlen > 0:
        return read.reference_start, read.reference_start + tlen
    else:
        return read.reference_start + tlen, read.reference_start