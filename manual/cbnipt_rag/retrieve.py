"""
retrieve.py
쿼리 → Qdrant 검색 → 결과 출력 (RAG retrieval 검증용)

Usage:
    python retrieve.py "DiGeorge syndrome cardiac defects"
    python retrieve.py "TBX1 haploinsufficiency" --top 5 --syndrome "DiGeorge syndrome"
    python retrieve.py "immune deficiency thymus" --type abstract
"""

import argparse
import json
from qdrant_client import QdrantClient
from qdrant_client.models import Filter, FieldCondition, MatchValue
from sentence_transformers import SentenceTransformer

QDRANT_PATH = "./qdrant_local"
COLLECTION  = "nipt_rag"
EMBED_MODEL = "pritamdeka/PubMedBERT-mnli-snli-scinli-scitail-mednli-stsb"


def build_filter(syndrome: str | None, chunk_type: str | None) -> Filter | None:
    conditions = []
    if syndrome:
        conditions.append(FieldCondition(key="syndrome", match=MatchValue(value=syndrome)))
    if chunk_type:
        conditions.append(FieldCondition(key="chunk_type", match=MatchValue(value=chunk_type)))
    if not conditions:
        return None
    return Filter(must=conditions)


def retrieve(query: str, top_k: int = 5, syndrome: str = None, chunk_type: str = None):
    model  = SentenceTransformer(EMBED_MODEL)
    client = QdrantClient(path=QDRANT_PATH)

    q_vec = model.encode([query], normalize_embeddings=True)[0].tolist()
    f     = build_filter(syndrome, chunk_type)

    results = client.query_points(
        collection_name=COLLECTION,
        query=q_vec,
        limit=top_k,
        query_filter=f,
        with_payload=True,
    ).points

    print(f"\n🔍 Query: {query!r}")
    if syndrome:   print(f"   Filter syndrome   : {syndrome}")
    if chunk_type: print(f"   Filter chunk_type : {chunk_type}")
    print(f"   Top-{top_k} results\n{'─'*70}")

    for i, hit in enumerate(results, 1):
        p = hit.payload
        print(f"\n[{i}] score={hit.score:.4f}  type={p.get('chunk_type')}  "
              f"gene={p.get('gene_symbol')}  pmid={p.get('pmid','')}")
        print(f"    syndrome: {p.get('syndrome')}")
        if p.get("title"):
            print(f"    title   : {p['title']}")
        # text 앞 300자만 미리보기
        preview = p.get("text", "")[:300].replace("\n", " ")
        print(f"    text    : {preview}…")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="검색 쿼리")
    parser.add_argument("--top",      type=int, default=5,  help="반환 개수")
    parser.add_argument("--syndrome", type=str, default=None)
    parser.add_argument("--type",     type=str, default=None,
                        choices=["gene_meta", "gene_summary", "abstract"])
    parser.add_argument("--model", default=EMBED_MODEL,
                        help="HuggingFace 모델명 또는 로컬 경로")
    args = parser.parse_args()
    EMBED_MODEL = args.model

    retrieve(args.query, top_k=args.top, syndrome=args.syndrome, chunk_type=args.type)