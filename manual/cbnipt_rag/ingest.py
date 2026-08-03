"""
ingest.py
JSON 파일(들) → 청킹 → 임베딩 → Qdrant persist DB 적재

Usage:
    python ingest.py data/digeorge_syndrome.json [data/trisomy21.json ...]
    python ingest.py data/          # 디렉토리 내 전체 JSON 일괄 처리
"""

import sys
import uuid
from pathlib import Path

from qdrant_client import QdrantClient
from qdrant_client.models import (
    Distance, VectorParams, PointStruct, PayloadSchemaType
)
from sentence_transformers import SentenceTransformer

from parse_json import parse_disease_json

# ── 설정 ─────────────────────────────────────────────────────────────
QDRANT_PATH    = "./qdrant_local"
COLLECTION     = "nipt_rag"
EMBED_MODEL    = "pritamdeka/PubMedBERT-mnli-snli-scinli-scitail-mednli-stsb"
# 로컬 모델 경로 사용 시 위 값 대신 아래처럼 지정
# EMBED_MODEL = "/path/to/local/pubmedbert"
VECTOR_DIM     = 768
BATCH_SIZE     = 32


def id_to_int(hex_id: str) -> int:
    """MD5 hex → Qdrant uint64 (상위 8바이트)"""
    return int(hex_id[:16], 16)


def get_json_files(paths: list[str]) -> list[Path]:
    files = []
    for p in paths:
        path = Path(p)
        if path.is_dir():
            files.extend(sorted(path.glob("*.json")))
        elif path.suffix == ".json":
            files.append(path)
    return files


def init_collection(client: QdrantClient):
    existing = [c.name for c in client.get_collections().collections]
    if COLLECTION not in existing:
        client.create_collection(
            collection_name=COLLECTION,
            vectors_config=VectorParams(size=VECTOR_DIM, distance=Distance.COSINE),
        )
        print(f"✅ 컬렉션 생성: {COLLECTION}")
    else:
        print(f"ℹ️  기존 컬렉션 재사용: {COLLECTION}")


def ingest(json_paths: list[str]):
    # 1) 파일 수집
    files = get_json_files(json_paths)
    if not files:
        print("❌ JSON 파일을 찾을 수 없습니다.")
        sys.exit(1)
    print(f"\n📂 처리 대상: {[f.name for f in files]}\n")

    # 2) 전체 chunk 수집
    all_chunks = []
    for f in files:
        chunks = parse_disease_json(str(f))
        print(f"  {f.name}: {len(chunks)}개 chunk")
        all_chunks.extend(chunks)
    print(f"\n📝 총 chunk: {len(all_chunks)}개\n")

    # 3) 임베딩 모델 로드
    print(f"🔄 임베딩 모델 로드: {EMBED_MODEL}")
    model = SentenceTransformer(EMBED_MODEL)

    # 4) 배치 임베딩
    texts = [c["text"] for c in all_chunks]
    print(f"🔄 임베딩 중... (배치 크기={BATCH_SIZE})")
    vectors = model.encode(
        texts,
        batch_size=BATCH_SIZE,
        show_progress_bar=True,
        normalize_embeddings=True,   # Cosine 거리용
    )
    print(f"✅ 임베딩 완료: shape={vectors.shape}\n")

    # 5) Qdrant 적재
    client = QdrantClient(path=QDRANT_PATH)
    init_collection(client)

    points = []
    for chunk, vec in zip(all_chunks, vectors):
        points.append(PointStruct(
            id=id_to_int(chunk["id"]),
            vector=vec.tolist(),
            payload={**chunk["metadata"], "text": chunk["text"]},
        ))

    # 배치 upsert
    for i in range(0, len(points), BATCH_SIZE):
        batch = points[i: i + BATCH_SIZE]
        client.upsert(collection_name=COLLECTION, points=batch)

    info = client.get_collection(COLLECTION)
    print(f"✅ Qdrant 적재 완료")
    print(f"   저장 경로 : {QDRANT_PATH}")
    print(f"   컬렉션   : {COLLECTION}")
    print(f"   총 포인트 : {info.points_count}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", default=["data/"])
    parser.add_argument("--model", default=EMBED_MODEL,
                        help="HuggingFace 모델명 또는 로컬 경로")
    args = parser.parse_args()
    EMBED_MODEL = args.model
    ingest(args.paths)