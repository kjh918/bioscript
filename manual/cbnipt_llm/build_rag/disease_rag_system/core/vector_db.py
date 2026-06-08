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
