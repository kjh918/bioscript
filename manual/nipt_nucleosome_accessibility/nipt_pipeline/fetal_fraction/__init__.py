"""
fetal_fraction/__init__.py
==========================
두 FF 방법을 병합하여 단일 FetalFractionResult로 반환하는
run() 함수 제공
"""
from __future__ import annotations

import logging
from typing import Optional

from .chry     import estimate as _chry_estimate
from .seqff    import estimate as _seqff_estimate
from .utils    import FetalFractionResult, write_json
from .plotting import plot_ff_correlation

log = logging.getLogger(__name__)


def run(
    bincount_tsv:   str,
    normalized_tsv: str,
    output_json:    str,
    sample_id:      str            = "",
    seqff_intercept: Optional[float] = None,
    seqff_slope:     Optional[float] = None,
    plot_path:       Optional[str]   = None,   # FF correlation plot PNG 경로
) -> FetalFractionResult:
    """
    chrY + SeqFF 두 방법 모두 실행 후 결과 병합 → JSON 저장

    Parameters
    ----------
    bincount_tsv    : run_binning 출력 (chrY FF 입력)
    normalized_tsv  : run_normalize 출력 (SeqFF 입력)
    output_json     : 결과 JSON 저장 경로
    sample_id       : 샘플 ID
    seqff_intercept : SeqFF 절편 오버라이드 (None이면 기본값)
    seqff_slope     : SeqFF 기울기 오버라이드 (None이면 기본값)
    plot_path       : FF correlation plot PNG 경로 (None이면 plot 스킵)

    Returns
    -------
    두 방법 필드가 합쳐진 FetalFractionResult
    """
    from .seqff import _DEFAULT_INTERCEPT, _DEFAULT_SLOPE

    log.info("=== Fetal Fraction 추정 시작: %s ===", sample_id)

    # chrY 기반
    chry_result = _chry_estimate(bincount_tsv, sample_id)

    # SeqFF 기반
    intercept    = seqff_intercept if seqff_intercept is not None else _DEFAULT_INTERCEPT
    slope        = seqff_slope     if seqff_slope     is not None else _DEFAULT_SLOPE
    seqff_result = _seqff_estimate(normalized_tsv, sample_id, intercept, slope)

    # 두 결과 병합
    merged = FetalFractionResult(
        sample_id           = sample_id,
        ff_chry             = chry_result.ff_chry,
        chry_mean_cov       = chry_result.chry_mean_cov,
        autosome_mean_cov   = chry_result.autosome_mean_cov,
        chry_bin_count      = chry_result.chry_bin_count,
        ff_seqff            = seqff_result.ff_seqff,
        short_fraction_mean = seqff_result.short_fraction_mean,
        seqff_intercept     = seqff_result.seqff_intercept,
        seqff_slope         = seqff_result.seqff_slope,
        method_note         = (
            f"chrY={chry_result.method_note} | "
            f"SeqFF={seqff_result.method_note}"
        ),
    )

    write_json(merged, output_json)
    log.info(
        "FF 결과: ff_chry=%.2f%% ff_seqff=%.2f%% → %s",
        merged.ff_chry  or 0.0,
        merged.ff_seqff or 0.0,
        output_json,
    )

    # FF correlation plot
    if plot_path:
        plot_ff_correlation(
            normalized_tsv = normalized_tsv,
            output_path    = plot_path,
            sample_id      = sample_id,
            ff_seqff       = merged.ff_seqff,
            ff_chry        = merged.ff_chry,
        )

    return merged