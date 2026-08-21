"""
fetal_fraction/utils.py
=======================
FF 결과 dataclass + JSON 저장
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional


@dataclass
class FetalFractionResult:
    sample_id:    str

    # chrY 기반
    ff_chry:           Optional[float]  # % 단위
    chry_mean_cov:     Optional[float]  # chrY 평균 bin coverage
    autosome_mean_cov: Optional[float]  # autosome 평균 bin coverage (chrY 제외)
    chry_bin_count:    Optional[int]    # 사용된 chrY bin 수

    # SeqFF 기반
    ff_seqff:             Optional[float]  # % 단위
    short_fraction_mean:  Optional[float]  # reference bin short/total 평균
    seqff_intercept:      float = 0.0
    seqff_slope:          float = 0.0

    # 메타
    method_note: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        # float 소수점 정리
        for k, v in d.items():
            if isinstance(v, float):
                d[k] = round(v, 4) if v == v else None  # NaN → None
        return d


def write_json(result: FetalFractionResult, output_path: str) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(result.to_dict(), fh, indent=2)