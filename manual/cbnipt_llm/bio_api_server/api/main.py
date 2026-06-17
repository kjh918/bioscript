from contextlib import asynccontextmanager
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
    full_prompt = f"<|im_start|>system\n{SYSTEM_PROMPT}<|im_end|>\n<|im_start|>user\nContext: {context}\nQuestion: {req.prompt}<|im_end|>\n<|im_start|>assistant\n"
    
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
