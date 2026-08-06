"""
utils.py - WPS 분석 공통 유틸리티

모든 모듈에서 공유하는 로거, 타이머, 파일 헬퍼를 제공한다.
분석 도메인 로직은 포함하지 않는다.
"""

import logging
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Generator, Optional


# ── 상수 ──────────────────────────────────────────────────────────────────────

LOG_FORMAT = "[%(asctime)s] %(levelname)-8s %(name)s - %(message)s"
LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"


# ── 로거 ──────────────────────────────────────────────────────────────────────

def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """표준 포맷 로거를 생성하여 반환한다.

    Args:
        name: 로거 이름 (보통 __name__ 전달)
        level: 로깅 레벨 (기본 INFO)

    Returns:
        설정된 Logger 인스턴스
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt=LOG_DATEFMT))
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger


# ── 타이머 ────────────────────────────────────────────────────────────────────

@contextmanager
def timer(label: str, logger: Optional[logging.Logger] = None) -> Generator:
    """코드 블록 실행 시간을 측정하는 컨텍스트 매니저.

    Args:
        label: 출력할 레이블 문자열
        logger: 사용할 Logger (None이면 stdout print)

    Yields:
        None

    Example:
        with timer("BAM 파싱", logger):
            fragments = parse_bam(bam_path)
    """
    start = time.perf_counter()
    yield
    elapsed = time.perf_counter() - start
    msg = f"{label} 완료: {elapsed:.2f}s"
    if logger:
        logger.info(msg)
    else:
        print(msg)


# ── 파일 유틸리티 ─────────────────────────────────────────────────────────────

def ensure_dir(path: str | Path) -> Path:
    """디렉토리가 없으면 생성하고 Path 객체를 반환한다.

    Args:
        path: 생성할 디렉토리 경로

    Returns:
        생성된 Path 객체
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def check_file_exists(path: str | Path, label: str = "파일") -> Path:
    """파일 존재 여부를 확인하고 Path 객체를 반환한다.

    Args:
        path: 확인할 파일 경로
        label: 에러 메시지에 표시할 파일 설명

    Returns:
        존재하는 파일의 Path 객체

    Raises:
        FileNotFoundError: 파일이 없을 경우
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{label}을 찾을 수 없습니다: {p}")
    return p


def check_bam_index(bam_path: str | Path) -> None:
    """BAM/CRAM 인덱스 파일 존재 여부를 확인한다.

    Args:
        bam_path: BAM 또는 CRAM 파일 경로

    Raises:
        FileNotFoundError: 인덱스 파일이 없을 경우
    """
    p = Path(bam_path)
    candidates = [
        Path(str(p) + ".bai"),
        Path(str(p) + ".crai"),
        p.with_suffix(".bai"),
    ]
    if not any(c.exists() for c in candidates):
        raise FileNotFoundError(
            f"BAM 인덱스가 없습니다. 먼저 실행하세요: samtools index {p}"
        )


# ── 진행 표시 ─────────────────────────────────────────────────────────────────

def log_progress(
    current: int,
    total: int,
    label: str,
    logger: logging.Logger,
    interval: int = 500_000,
) -> None:
    """일정 간격으로 진행상황을 로그에 출력한다.

    Args:
        current: 현재 처리 건수
        total: 전체 건수 (0이면 퍼센트 미출력)
        label: 출력할 레이블
        logger: Logger 인스턴스
        interval: 출력 간격 (기본 500,000 reads)
    """
    if current % interval == 0 and current > 0:
        if total > 0:
            pct = current / total * 100
            logger.info(f"{label}: {current:,}/{total:,} ({pct:.1f}%)")
        else:
            logger.info(f"{label}: {current:,} 처리 완료")