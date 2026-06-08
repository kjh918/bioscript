import os

def create_file(path: str, content: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content.strip() + "\n")

# 1. Core Exceptions (Rule 0 강제)
CORE_EXCEPTIONS = """
class MissingParameterError(ValueError):
    pass

class DatabaseConnectionError(RuntimeError):
    pass
"""

# 2. Components Schemas (데이터 유효성 검증)
COMPONENTS_SCHEMAS = """
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
"""

# 3. Core VectorDB Interface (명확한 규약 설정)
CORE_VECTOR_DB = """
from typing import List
from components.schemas import DiseaseDocument, SearchQuery
from core.exceptions import MissingParameterError

class VectorStoreInterface:
    def insert_documents(self, documents: List[DiseaseDocument]) -> None:
        if not documents:
            raise MissingParameterError("insert할 documents 리스트가 비어있습니다.")
        raise NotImplementedError("하위 클래스에서 구현해야 합니다.")

    def search_by_disease(self, search_query: SearchQuery) -> List[DiseaseDocument]:
        if not search_query:
            raise MissingParameterError("search_query 객체가 누락되었습니다.")
        raise NotImplementedError("하위 클래스에서 구현해야 합니다.")
"""

if __name__ == "__main__":
    base_dir = "disease_rag_system"
    
    # Init files
    create_file(f"{base_dir}/core/__init__.py", "")
    create_file(f"{base_dir}/components/__init__.py", "")
    create_file(f"{base_dir}/tasks/__init__.py", "")
    create_file(f"{base_dir}/pipelines/__init__.py", "")
    create_file(f"{base_dir}/executors/__init__.py", "")
    
    # Core logic files
    create_file(f"{base_dir}/core/exceptions.py", CORE_EXCEPTIONS)
    create_file(f"{base_dir}/core/vector_db.py", CORE_VECTOR_DB)
    create_file(f"{base_dir}/components/schemas.py", COMPONENTS_SCHEMAS)
    
    # Placeholder for other modules
    create_file(f"{base_dir}/tasks/json_textifier.py", "# JSON을 자연어 청크로 변환하는 단일 작업 로직")
    create_file(f"{base_dir}/tasks/document_chunker.py", "# 논문 텍스트를 문맥 단위로 분할하는 로직")
    create_file(f"{base_dir}/pipelines/ingestion_pipeline.py", "# 데이터 전처리 및 DB 적재 파이프라인")
    create_file(f"{base_dir}/pipelines/retrieval_pipeline.py", "# 메타데이터 필터링 기반 RAG 검색 파이프라인")
    create_file(f"{base_dir}/executors/web_api.py", "# FastAPI 기반 Web API 엔드포인트")
    
    print(f"✅ RAG 시스템 기본 디렉토리 및 Core 구조({base_dir})가 성공적으로 생성되었습니다.")