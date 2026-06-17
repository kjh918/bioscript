import os
from pathlib import Path

def setup_fastapi_rag(base_dir: str):
    if not base_dir: raise ValueError("base_dir이 누락되었습니다.")

    base_path = Path(base_dir)
    for folder in ["core", "api", "scripts", "data"]:
        (base_path / folder).mkdir(parents=True, exist_ok=True)
        if folder != "data": (base_path / folder / "__init__.py").touch()

    # 1. 오프라인 인덱스 빌더 (스크립트)
    # [MODIFIED] 매번 JSON을 파싱하지 않고 FAISS 바이너리 인덱스로 저장
    build_index_code = '''import argparse
from pathlib import Path
from langchain_community.vectorstores import FAISS
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_core.documents import Document
import json

def build_offline_index(json_path: str, output_dir: str, embed_model: str):
    if not json_path: raise ValueError("json_path가 필요합니다.")
    if not output_dir: raise ValueError("output_dir이 필요합니다.")
    
    print(f"1. Loading Embedding Model: {embed_model}")
    embeddings = HuggingFaceEmbeddings(model_name=embed_model, model_kwargs={'device': 'cpu'})
    
    print(f"2. Parsing JSON: {json_path}")
    raw_data = []
    with open(json_path, "r", encoding="utf-8") as f:
        if json_path.endswith(".jsonl"):
            for line in f:
                if line.strip(): raw_data.append(json.loads(line))
        else:
            raw_data = json.load(f)

    print("3. Chunking Data")
    text_splitter = RecursiveCharacterTextSplitter(chunk_size=500, chunk_overlap=50)
    chunked_docs = []
    
    for doc in raw_data:
        parts = []
        if "title" in doc: parts.append(f"Title: {doc['title']}")
        if "abstract" in doc: parts.append(f"Abstract: {doc['abstract']}")
        if "text" in doc: parts.append(f"Text: {doc['text']}")
        page_content = "\\n".join(parts)
        
        metadata = {"pmid": doc.get("pmid", ""), "accession": doc.get("accession", "")}
        for split_text in text_splitter.split_text(page_content):
            chunked_docs.append(Document(page_content=split_text, metadata=metadata))

    print("4. Building FAISS Index & Saving to Disk...")
    vectorstore = FAISS.from_documents(chunked_docs, embeddings)
    vectorstore.save_local(output_dir)
    print(f"[SUCCESS] Index saved to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_path", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--embed_model", type=str, default="BAAI/bge-m3")
    args = parser.parse_args()
    build_offline_index(args.json_path, args.output_dir, args.embed_model)
'''

    # 2. 시스템 프롬프트 및 로직 코어
    core_logic_code = '''SYSTEM_PROMPT = """You are an AI assistant specialized in fetal and embryonic genetic disorders.
You must NOT make final medical diagnoses or treatment decisions.
You must:
- Use only the provided database evidence and retrieved references.
- Prioritize ACMG/AMP guideline-based interpretation principles.
- Always provide responses in Korean.
- Mention evidence source categories: ClinVar, PubMed, ACMG guideline.

Output format must follow:
[질환 요약]
[유전적 특징]
[임상적 영향]
[병원성/위험도]
[근거 문헌]
[ACMG/가이드라인 해석]
"""

def format_context(retrieved_docs) -> str:
    if not retrieved_docs:
        return "No reliable evidence was retrieved."
    
    blocks = []
    for i, (doc, score) in enumerate(retrieved_docs, start=1):
        pmid = doc.metadata.get("pmid", "")
        blocks.append(f"[Evidence {i}]\\npmid: {pmid}\\nretrieval_score: {score}\\n\\n{doc.page_content}")
    return "\\n\\n".join(blocks)
'''

    # 3. FastAPI 서버 (메모리 상주)
    # [MODIFIED] 애플리케이션 시작 시 Model과 FAISS DB를 한 번만 로드하여 메모리에 유지
    server_code = '''from contextlib import asynccontextmanager
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM
from langchain_community.vectorstores import FAISS
from langchain_huggingface import HuggingFaceEmbeddings
from core.logic import SYSTEM_PROMPT, format_context

# 전역 상태를 관리할 딕셔너리
ml_models = {}

@asynccontextmanager
async def lifespan(app: FastAPI):
    # [설정] 환경변수나 Config에서 주입받는 것이 좋으나 직관성을 위해 하드코딩
    MODEL_ID = "Qwen/Qwen2.5-0.5B-Instruct"
    EMBED_MODEL = "BAAI/bge-m3"
    INDEX_DIR = "./data/faiss_index"
    
    print("🚀 [1/3] Loading Embedding Model & Offline DB...")
    try:
        embeddings = HuggingFaceEmbeddings(model_name=EMBED_MODEL, model_kwargs={'device': 'cpu'})
        ml_models["db"] = FAISS.load_local(INDEX_DIR, embeddings, allow_dangerous_deserialization=True)
    except Exception as e:
        print(f"DB 로드 실패 (미리 인덱스를 빌드했는지 확인하세요): {e}")
        ml_models["db"] = None

    print(f"🚀 [2/3] Loading LLM Tokenizer ({MODEL_ID})...")
    ml_models["tokenizer"] = AutoTokenizer.from_pretrained(MODEL_ID)
    
    print(f"🚀 [3/3] Loading LLM Model ({MODEL_ID}) into Memory...")
    ml_models["llm"] = AutoModelForCausalLM.from_pretrained(
        MODEL_ID, torch_dtype=torch.float32, device_map="cpu"
    )
    
    print("✅ Server is Ready!")
    yield
    
    print("Shutdown: Clearing Memory...")
    ml_models.clear()

app = FastAPI(lifespan=lifespan, title="Bioinformatics RAG API")

class QueryRequest(BaseModel):
    prompt: str
    top_k: int = 5
    max_new_tokens: int = 512
    temperature: float = 0.3
    top_p: float = 0.9

@app.post("/api/v1/generate")
async def generate_answer(req: QueryRequest):
    if not ml_models.get("db") or not ml_models.get("llm"):
        raise HTTPException(status_code=503, detail="모델이나 DB가 아직 로드되지 않았습니다.")
        
    if not req.prompt:
        raise HTTPException(status_code=400, detail="prompt가 누락되었습니다.")

    # 1. 고속 벡터 검색
    results = ml_models["db"].similarity_search_with_score(req.prompt, k=req.top_k)
    context = format_context(results)
    
    # 2. 프롬프트 구성 (사용자 코드 반영)
    full_prompt = f"<|im_start|>system\\n{SYSTEM_PROMPT}<|im_end|>\\n<|im_start|>user\\nContext: {context}\\nQuestion: {req.prompt}<|im_end|>\\n<|im_start|>assistant\\n"
    
    # 3. 인퍼런스 (이미 메모리에 올라간 모델 사용)
    tokenizer = ml_models["tokenizer"]
    model = ml_models["llm"]
    
    inputs = tokenizer(full_prompt, return_tensors="pt").to(model.device)
    
    with torch.no_grad():
        outputs = model.generate(
            **inputs, 
            max_new_tokens=req.max_new_tokens,
            temperature=req.temperature,
            top_p=req.top_p
        )
        
    input_length = inputs.input_ids.shape[1]
    generated_tokens = outputs[0][input_length:]
    response_text = tokenizer.decode(generated_tokens, skip_special_tokens=True)
    
    return {
        "status": "success",
        "retrieved_context": context,
        "response": response_text
    }
'''

    with open(base_path / "scripts" / "build_index.py", "w", encoding="utf-8") as f: f.write(build_index_code)
    with open(base_path / "core" / "logic.py", "w", encoding="utf-8") as f: f.write(core_logic_code)
    with open(base_path / "api" / "main.py", "w", encoding="utf-8") as f: f.write(server_code)

    print(f"[SUCCESS] {base_dir} 경로에 백엔드 시스템 구축 완료.")

if __name__ == "__main__":
    setup_fastapi_rag("bio_api_server")