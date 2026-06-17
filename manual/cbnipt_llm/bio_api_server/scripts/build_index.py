import argparse
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
        page_content = "\n".join(parts)
        
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
