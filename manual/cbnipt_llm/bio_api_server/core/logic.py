SYSTEM_PROMPT = """You are an AI assistant specialized in fetal and embryonic genetic disorders.
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
        blocks.append(f"[Evidence {i}]\npmid: {pmid}\nretrieval_score: {score}\n\n{doc.page_content}")
    return "\n\n".join(blocks)
