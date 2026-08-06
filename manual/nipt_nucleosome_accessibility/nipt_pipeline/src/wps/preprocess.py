"""
preprocess.py - WPS 분석 전처리 모듈

BAM 파일 유효성 검증, 분석 대상 염색체 목록 추출,
fragment 필터 기준 적용 후 정제된 입력 정보를 반환한다.
분석 계산 로직은 포함하지 않는다.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pysam

from utils import check_bam_index, check_file_exists, get_logger

# ── 상수 ──────────────────────────────────────────────────────────────────────

DEFAULT_MIN_FRAG_LEN: int = 120   # fragment 최소 길이 (bp)
DEFAULT_MAX_FRAG_LEN: int = 180   # fragment 최대 길이 (bp) — 모노뉴클레오솜
DEFAULT_MIN_MAPQ: int = 20        # 최소 매핑 품질

AUTOSOMAL_CHROMS_HG38 = [f"chr{i}" for i in range(1, 23)]
AUTOSOMAL_CHROMS_HG19 = [str(i) for i in range(1, 23)]


# ── 데이터 클래스 ─────────────────────────────────────────────────────────────

@dataclass
class PreprocessResult:
    """전처리 결과를 담는 데이터 클래스.

    Attributes:
        bam_path: 검증된 BAM 파일 경로
        chroms: 분석 대상 염색체 목록
        min_frag_len: fragment 최소 길이
        max_frag_len: fragment 최대 길이
        min_mapq: 최소 MapQ
        chrom_lengths: 염색체별 길이 딕셔너리
    """
    bam_path: Path
    chroms: list[str]
    min_frag_len: int
    max_frag_len: int
    min_mapq: int
    chrom_lengths: dict[str, int] = field(default_factory=dict)


# ── 공개 함수 ─────────────────────────────────────────────────────────────────

def run_preprocess(
    bam_path: str,
    min_frag_len: int = DEFAULT_MIN_FRAG_LEN,
    max_frag_len: int = DEFAULT_MAX_FRAG_LEN,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    chroms: Optional[list[str]] = None,
    logger: Optional[logging.Logger] = None,
) -> PreprocessResult:
    """BAM 파일을 검증하고 분석 준비 정보를 반환한다.

    Args:
        bam_path: 입력 BAM 파일 경로
        min_frag_len: fragment 최소 길이 (bp)
        max_frag_len: fragment 최대 길이 (bp)
        min_mapq: 최소 매핑 품질 스코어
        chroms: 분석 대상 염색체 목록 (None이면 자동 감지)
        logger: Logger 인스턴스

    Returns:
        PreprocessResult: 검증된 분석 입력 정보

    Raises:
        FileNotFoundError: BAM 파일 또는 인덱스가 없을 경우
        ValueError: fragment 길이 파라미터가 잘못된 경우
    """
    log = logger or get_logger(__name__)

    # 파일/인덱스 검증
    bam_p = check_file_exists(bam_path, "BAM 파일")
    check_bam_index(bam_p)
    log.info(f"BAM 파일 확인: {bam_p}")

    # 파라미터 검증
    _validate_frag_params(min_frag_len, max_frag_len)
    log.info(f"Fragment 필터: {min_frag_len}–{max_frag_len} bp, MapQ ≥ {min_mapq}")

    # BAM 헤더에서 염색체 정보 추출
    with pysam.AlignmentFile(str(bam_p), "rb") as bam:
        chrom_lengths = _extract_chrom_lengths(bam)
        available = set(chrom_lengths.keys())

    # 분석 대상 염색체 결정
    target_chroms = _resolve_chroms(chroms, available, log)
    log.info(f"분석 염색체 수: {len(target_chroms)}")

    return PreprocessResult(
        bam_path=bam_p,
        chroms=target_chroms,
        min_frag_len=min_frag_len,
        max_frag_len=max_frag_len,
        min_mapq=min_mapq,
        chrom_lengths={c: chrom_lengths[c] for c in target_chroms},
    )


# ── 내부 함수 ─────────────────────────────────────────────────────────────────

def _validate_frag_params(min_len: int, max_len: int) -> None:
    """Fragment 길이 파라미터 유효성을 검사한다.

    Args:
        min_len: 최소 fragment 길이
        max_len: 최대 fragment 길이

    Raises:
        ValueError: 파라미터가 유효하지 않은 경우
    """
    if min_len <= 0 or max_len <= 0:
        raise ValueError("Fragment 길이는 양수여야 합니다.")
    if min_len >= max_len:
        raise ValueError(
            f"min_frag_len({min_len}) < max_frag_len({max_len}) 이어야 합니다."
        )


def _extract_chrom_lengths(bam: pysam.AlignmentFile) -> dict[str, int]:
    """BAM 헤더에서 염색체 이름과 길이를 추출한다.

    Args:
        bam: 열린 pysam AlignmentFile 객체

    Returns:
        {chrom_name: length} 딕셔너리
    """
    return {
        sq["SN"]: sq["LN"]
        for sq in bam.header.to_dict().get("SQ", [])
    }


def _resolve_chroms(
    requested: Optional[list[str]],
    available: set[str],
    logger: logging.Logger,
) -> list[str]:
    """분석 대상 염색체 목록을 결정한다.

    명시적으로 요청된 염색체가 없으면 BAM에서 감지된 autosomal 염색체를
    자동으로 선택한다 (hg38 chr1-22 또는 hg19 1-22 형식).

    Args:
        requested: 사용자가 요청한 염색체 목록 (None이면 자동)
        available: BAM 헤더에 존재하는 염색체 집합
        logger: Logger 인스턴스

    Returns:
        정렬된 분석 대상 염색체 목록

    Raises:
        ValueError: 요청된 염색체가 BAM에 없는 경우
    """
    if requested:
        missing = set(requested) - available
        if missing:
            raise ValueError(f"BAM에 없는 염색체: {sorted(missing)}")
        return requested

    # 자동 감지: hg38 → hg19 순서로 시도
    for candidates in (AUTOSOMAL_CHROMS_HG38, AUTOSOMAL_CHROMS_HG19):
        found = [c for c in candidates if c in available]
        if found:
            logger.info(f"염색체 형식 감지: {'hg38' if 'chr1' in found[0] else 'hg19'}")
            return found

    # fallback: BAM에 있는 전체 중 숫자 포함 항목
    fallback = sorted(
        [c for c in available if any(ch.isdigit() for ch in c)],
        key=lambda x: int("".join(filter(str.isdigit, x)) or 0),
    )
    if not fallback:
        raise ValueError("분석 가능한 autosomal 염색체를 찾을 수 없습니다.")
    logger.warning(f"알 수 없는 염색체 형식 — fallback 사용: {fallback[:3]}...")
    return fallback