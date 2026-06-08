from dataclasses import dataclass, field
from typing import Dict, Any
from core.exceptions import MissingParameterError

@dataclass
class DiseaseDocument:
    disease_id: str
    content: str
    metadata: Dict[str, Any]

    def __post_init__(self):
        # [Rule 0] Default 값 의존 배제 및 엄격한 파라미터 검증
        if not self.disease_id:
            raise MissingParameterError("disease_id는 필수 값입니다. (예: C0001080)")
        if not self.content:
            raise MissingParameterError("문서의 content(텍스트 데이터)가 비어있습니다.")
        if not isinstance(self.metadata, dict):
            raise TypeError("metadata는 반드시 dictionary 형태여야 합니다.")

@dataclass
class SearchQuery:
    query_text: str
    target_disease_id: str
    top_k: int

    def __post_init__(self):
        if not self.query_text:
            raise MissingParameterError("검색할 query_text가 누락되었습니다.")
        if not self.target_disease_id:
            raise MissingParameterError("필터링을 위한 target_disease_id가 누락되었습니다.")
        if not isinstance(self.top_k, int) or self.top_k <= 0:
            raise ValueError("top_k는 1 이상의 정수여야 합니다.")
