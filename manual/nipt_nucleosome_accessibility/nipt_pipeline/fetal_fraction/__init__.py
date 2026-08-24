"""
fetal_fraction/__init__.py
==========================
chrY 기반 FF 추정 (SeqFF 제거)

SeqFF 제거 근거:
- Kim 2015 SeqFF는 single-end sequencing 기반 multivariate regression
- 현재 데이터는 paired-end → chrY coverage 기반이 더 적합
- female fetus는 chrY FF=None으로 반환 (한계 명시)
"""
from __future__ import annotations

import logging
from typing import Optional

from .chry  import estimate as _chry_estimate
from .utils import FetalFractionResult, write_json

log = logging.getLogger(__name__)


def run(
    bincount_tsv:   str,
    normalized_tsv: str,        # 시그니처 호환 유지 (미사용)
    output_json:    str,
    sample_id:      str          = "",
    plot_path:      Optional[str] = None,
    **kwargs,                   # 구버전 호환 (seqff_intercept 등 무시)
) -> FetalFractionResult:
    """
    chrY 기반 FF 추정 → JSON 저장

    Parameters
    ----------
    bincount_tsv   : run_binning 출력 TSV (chrY coverage 기반)
    normalized_tsv : 미사용 (인터페이스 호환용)
    output_json    : 결과 JSON 저장 경로
    sample_id      : 샘플 ID
    plot_path      : FF correlation plot (normalized_tsv 필요 시 사용)

    Returns
    -------
    FetalFractionResult
      - male fetus:   ff_chry = 추정값 (%)
      - female fetus: ff_chry = None (chrY coverage 없음)
    """
    log.info("=== Fetal Fraction 추정 (chrY 기반): %s ===", sample_id)

    result = _chry_estimate(bincount_tsv, sample_id)
    write_json(result, output_json)

    if result.ff_chry is not None:
        log.info("FF (chrY): %.2f%%", result.ff_chry)
    else:
        log.info("FF (chrY): None (female fetus 또는 chrY coverage 없음)")

    # FF correlation plot (short/long count 기반, normalized_tsv 필요)
    if plot_path and normalized_tsv:
        try:
            from .plotting import plot_ff_correlation
            plot_ff_correlation(
                normalized_tsv = normalized_tsv,
                output_path    = plot_path,
                sample_id      = sample_id,
                ff_seqff       = None,
                ff_chry        = result.ff_chry,
            )
        except Exception as e:
            log.warning("FF correlation plot 실패: %s", e)

    return result