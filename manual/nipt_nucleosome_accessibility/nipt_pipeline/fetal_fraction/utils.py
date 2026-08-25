"""
fetal_fraction/utils.py
=======================
FF 결과 dataclass + JSON 저장

SeqFF 제거 이유:
- SeqFF (Kim 2015)는 25,312명 훈련 데이터 기반 고차원 multivariate regression
- 단순 intercept/slope 선형 모델이 아님
- single-end sequencing 기반 → paired-end 데이터에 적합하지 않음
- chrY 기반 FF가 paired-end에서 더 신뢰도 높음
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional


@dataclass
class FetalFractionResult:
    sample_id:    str

    # chrY 기반 (paired-end 표준)
    ff_chry:           Optional[float]  # % 단위 (male fetus만 유효)
    chry_mean_cov:     Optional[float]  # chrY 평균 bin coverage
    autosome_mean_cov: Optional[float]  # autosome 평균 bin coverage
    chry_bin_count:    Optional[int]    # 사용된 chrY bin 수
    ff_short_fragments_ratio: Optional[float]

    # 메타
    method_note: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        for k, v in d.items():
            if isinstance(v, float):
                d[k] = round(v, 4) if v == v else None
        return d


def write_json(result: FetalFractionResult, output_path: str) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(result.to_dict(), fh, indent=2)