"""
rag_query.py
NIPT 분석 결과 TSV + Qdrant RAG → llama-cpp-python → 임상 해석

Usage:
    python rag_query.py --tsv /path/to/result.tsv
    python rag_query.py --tsv /path/to/result.tsv --syndrome "Down syndrome"
    python rag_query.py --tsv /path/to/result.tsv --all
"""

import argparse
import csv
from pathlib import Path

from qdrant_client import QdrantClient
from qdrant_client.models import Filter, FieldCondition, MatchValue
from sentence_transformers import SentenceTransformer
from llama_cpp import Llama

# ── 고정 경로 ─────────────────────────────────────────────────────────
MODEL_PATH  = "/storage/home/jhkim/models/Qwen2.5-32B-Instruct-Q4_K_M.gguf"
EMBED_MODEL = "/storage/home/jhkim/models/pubmedbert-base-embeddings"
QDRANT_PATH = "./qdrant_local"
COLLECTION  = "nipt_rag"

# ── 추론 설정 ─────────────────────────────────────────────────────────
TOP_K      = 5
MAX_TOKENS = 1024
N_THREADS  = 16
CTX_SIZE   = 4096

# ── 모델 캐시 (한 번만 로드) ─────────────────────────────────────────
_llm         = None
_embed_model = None


def get_llm(model_path: str) -> Llama:
    global _llm
    if _llm is None:
        print(f"🔄 LLM 로드 중: {Path(model_path).name}")
        _llm = Llama(
            model_path=model_path,
            n_ctx=CTX_SIZE,
            n_threads=N_THREADS,
            n_gpu_layers=0,   # CPU only
            verbose=False,
        )
        print("✅ LLM 로드 완료")
    return _llm


def get_embed_model(path: str) -> SentenceTransformer:
    global _embed_model
    if _embed_model is None:
        _embed_model = SentenceTransformer(path)
    return _embed_model


# ── TSV 파싱 ─────────────────────────────────────────────────────────
def load_tsv(tsv_path: str) -> dict[str, list[dict]]:
    result = {}
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            syndrome = row["SYNDROME"]
            result.setdefault(syndrome, []).append(dict(row))
    return result


def summarize_syndrome(rows: list[dict]) -> dict:
    positive  = [r for r in rows if r["DIAGNOSIS"] == "POSITIVE"]
    suspicious = [r for r in rows if r["DIAGNOSIS"] == "SUSPICIOUS"]

    key_findings = []
    for r in positive + suspicious:
        key_findings.append({
            "feature":     r["FEATURE_NAME"],
            "type":        r["FEATURE_TYPE"],
            "diagnosis":   r["DIAGNOSIS"],
            "evidence":    r["EVIDENCE"],
            "observed_cn": r.get("OBSERVED_CN", ""),
            "chromosome":  r["CHROMOSOME"],
        })

    overall = "NEGATIVE"
    if positive:
        overall = "POSITIVE"
    elif suspicious:
        overall = "SUSPICIOUS"

    return {
        "syndrome":          rows[0]["SYNDROME"],
        "nipt_group":        rows[0]["NIPT_GROUP"],
        "overall_diagnosis": overall,
        "positive_count":    len(positive),
        "suspicious_count":  len(suspicious),
        "negative_count":    len([r for r in rows if r["DIAGNOSIS"] == "NEGATIVE"]),
        "key_findings":      key_findings,
    }


# ── RAG Retrieval ─────────────────────────────────────────────────────
def retrieve(syndrome: str, summary: dict, embed_model_path: str) -> list[dict]:
    model  = get_embed_model(embed_model_path)
    client = QdrantClient(path=QDRANT_PATH)

    findings_str = ", ".join(
        f["feature"] for f in summary["key_findings"]
    ) if summary["key_findings"] else syndrome

    query = f"{syndrome} {findings_str} clinical diagnosis genetic evidence"
    q_vec = model.encode([query], normalize_embeddings=True)[0].tolist()

    # syndrome 필터 우선, 없으면 전체 검색
    try:
        f = Filter(must=[
            FieldCondition(key="syndrome", match=MatchValue(value=syndrome))
        ])
        results = client.query_points(
            collection_name=COLLECTION,
            query=q_vec,
            limit=TOP_K,
            query_filter=f,
            with_payload=True,
        ).points
        if not results:
            raise ValueError("empty")
    except Exception:
        results = client.query_points(
            collection_name=COLLECTION,
            query=q_vec,
            limit=TOP_K,
            with_payload=True,
        ).points

    return [hit.payload for hit in results]


def build_prompt(summary: dict[str, Any], contexts: list[dict[str, Any]]) -> str:
    """
    cbNIPT 결과와 RAG 검색 문서를 기반으로 임상 해석 prompt를 생성합니다.

    권장 context chunk_type:
        - disease_description
        - sop
        - sop_classification
        - sop_limitation
        - guideline
        - abstract
        - literature
        - database

    각 context에는 가능한 경우 다음 metadata를 포함합니다.
        - chunk_type
        - reference_id
        - title
        - text
        - source
        - pmid
        - version
    """

    syndrome = summary.get("syndrome", "Unknown syndrome")
    overall_diagnosis = str(
        summary.get("overall_diagnosis", "UNKNOWN")
    ).upper()

    positive_count = summary.get("positive_count", 0)
    suspicious_count = summary.get("suspicious_count", 0)
    negative_count = summary.get("negative_count", 0)

    # ------------------------------------------------------------------
    # 1. RAG context를 목적별로 분류
    # ------------------------------------------------------------------
    disease_blocks: list[str] = []
    sop_blocks: list[str] = []
    limitation_blocks: list[str] = []
    literature_blocks: list[str] = []

    for i, ctx in enumerate(contexts, start=1):
        chunk_type = str(
            ctx.get("chunk_type")
            or ctx.get("document_type")
            or "literature"
        ).lower()

        reference_id = (
            ctx.get("reference_id")
            or ctx.get("ref_id")
            or f"REF-{i:03d}"
        )

        title = str(ctx.get("title", "")).strip()
        text = str(ctx.get("text", "")).strip()
        source = str(ctx.get("source", "")).strip()
        pmid = str(ctx.get("pmid", "")).strip()
        version = str(ctx.get("version", "")).strip()

        if not text:
            continue

        metadata_parts = [
            f"Reference ID: {reference_id}",
            f"Document Type: {chunk_type}",
        ]

        if title:
            metadata_parts.append(f"Title: {title}")

        if source:
            metadata_parts.append(f"Source: {source}")

        if pmid:
            metadata_parts.append(f"PMID: {pmid}")

        if version:
            metadata_parts.append(f"Version: {version}")

        metadata_parts.append(f"Content: {text}")

        block = "\n".join(metadata_parts)

        if chunk_type in {
            "disease_description",
            "disease",
            "clinical_description",
            "syndrome_description",
        }:
            disease_blocks.append(block)

        elif chunk_type in {
            "sop",
            "sop_classification",
            "classification_criteria",
            "interpretation_rule",
            "decision_rule",
        }:
            sop_blocks.append(block)

        elif chunk_type in {
            "sop_limitation",
            "test_limitation",
            "limitation",
            "limitations",
        }:
            limitation_blocks.append(block)

        else:
            literature_blocks.append(block)

    disease_context_str = (
        "\n\n---\n\n".join(disease_blocks)
        if disease_blocks
        else "No disease-description evidence was retrieved."
    )

    sop_context_str = (
        "\n\n---\n\n".join(sop_blocks)
        if sop_blocks
        else "No syndrome-specific SOP criteria were retrieved."
    )

    limitation_context_str = (
        "\n\n---\n\n".join(limitation_blocks)
        if limitation_blocks
        else "No syndrome-specific test limitations were retrieved."
    )

    literature_context_str = (
        "\n\n---\n\n".join(literature_blocks)
        if literature_blocks
        else "No additional supporting literature was retrieved."
    )

    # ------------------------------------------------------------------
    # 2. Marker-level finding 구성
    # ------------------------------------------------------------------
    key_findings = summary.get("key_findings", [])

    if key_findings:
        finding_blocks = []

        for index, finding in enumerate(key_findings, start=1):
            feature = finding.get("feature", "-")
            feature_type = finding.get("type", "-")
            diagnosis = finding.get("diagnosis", "-")
            evidence = finding.get("evidence", "-")
            observed_cn = finding.get("observed_cn", "-")

            feature_rank = finding.get(
                "feature_rank",
                finding.get("feat_rank", "-"),
            )
            genomic_region = finding.get(
                "genomic_region",
                finding.get("region", "-"),
            )
            z_score = finding.get("z_score", "-")
            coverage = finding.get("coverage", "-")
            qc_status = finding.get("qc_status", "-")

            finding_blocks.append(
                "\n".join(
                    [
                        f"Finding {index}:",
                        f"- Feature: {feature}",
                        f"- Feature type: {feature_type}",
                        f"- Feature rank: {feature_rank}",
                        f"- Genomic region: {genomic_region}",
                        f"- Marker diagnosis: {diagnosis}",
                        f"- Observed copy number: {observed_cn}",
                        f"- Z-score: {z_score}",
                        f"- Coverage: {coverage}",
                        f"- QC status: {qc_status}",
                        f"- Evidence: {evidence}",
                    ]
                )
            )

        findings_str = "\n\n".join(finding_blocks)

    else:
        findings_str = "No positive or suspicious marker-level findings."

    # ------------------------------------------------------------------
    # 3. Prompt 생성
    # ------------------------------------------------------------------
    return (
        "<|im_start|>system\n"
        "You are a clinical genomics expert specializing in cbNIPT, "
        "copy-number variation, and prenatal diagnosis.\n\n"

        "Interpret the provided screening result using only the marker-level "
        "findings and retrieved RAG evidence.\n\n"

        "Mandatory interpretation rules:\n"
        "1. cbNIPT is a screening test and must not be described as a "
        "definitive fetal diagnosis.\n"
        "2. Do not invent clinical facts, thresholds, SOP criteria, or "
        "references.\n"
        "3. Use the retrieved SOP evidence as the primary basis for explaining "
        "the final classification.\n"
        "4. Use disease-description and literature evidence to explain the "
        "clinical relevance of the syndrome or genomic region.\n"
        "5. When the result is SUSPICIOUS, AMBIGUOUS, INDETERMINATE, or "
        "INCONCLUSIVE, explain specifically why it is uncertain.\n"
        "6. For an uncertain result, distinguish criteria that were met from "
        "criteria that were not met.\n"
        "7. Do not state only that the result is ambiguous. Connect the "
        "observed marker findings to the relevant SOP rule.\n"
        "8. Cite evidence using the exact provided Reference ID in square "
        "brackets, for example [SOP-PWS-001].\n"
        "9. Never cite a Reference ID that is absent from the retrieved "
        "context.\n"
        "10. If the evidence required for a conclusion is missing, explicitly "
        "state '근거 자료가 충분하지 않음'.\n"
        "11. Respond in Korean.\n"
        "12. Keep the response concise, clear, and suitable for a clinical "
        "screening report.\n"
        "<|im_end|>\n"

        "<|im_start|>user\n"
        "## cbNIPT Analysis Result\n"
        f"Syndrome: {syndrome}\n"
        f"Overall Diagnosis: {overall_diagnosis}\n"
        f"Positive Markers: {positive_count}\n"
        f"Suspicious Markers: {suspicious_count}\n"
        f"Negative Markers: {negative_count}\n\n"

        "## Marker-Level Findings\n"
        f"{findings_str}\n\n"

        "## Retrieved Disease Description\n"
        f"{disease_context_str}\n\n"

        "## Retrieved SOP Classification Criteria\n"
        f"{sop_context_str}\n\n"

        "## Retrieved Test Limitations\n"
        f"{limitation_context_str}\n\n"

        "## Retrieved Supporting Literature\n"
        f"{literature_context_str}\n\n"

        "## Task\n"
        "Provide a brief clinical interpretation using the following format.\n\n"

        "### 1. 결과 요약\n"
        "- 대상 질환과 최종 screening classification을 기술한다.\n"
        "- 핵심 marker 또는 genomic finding을 간략히 요약한다.\n\n"

        "### 2. 질환 설명\n"
        "- 검색된 질환 설명을 기준으로 질환의 원인 영역, 관련 유전자 또는 "
        "주요 임상적 특징을 1~3문장으로 설명한다.\n"
        "- 사용한 문장 끝에 Reference ID를 표시한다.\n\n"

        "### 3. 판정 근거\n"
        "- 관찰된 marker 결과와 SOP 판정 조건을 연결한다.\n"
        "- 어떤 feature가 판정에 가장 큰 영향을 주었는지 설명한다.\n"
        "- FEAT_RANK 또는 marker 우선순위가 제공된 경우 이를 반영한다.\n"
        "- 사용한 SOP Reference ID를 표시한다.\n\n"

        "### 4. 모호성 설명\n"
        "- 결과가 SUSPICIOUS, AMBIGUOUS, INDETERMINATE 또는 "
        "INCONCLUSIVE인 경우 작성한다.\n"
        "- 먼저 충족된 기준을 설명한다.\n"
        "- 다음으로 충족되지 않은 양성 기준을 설명한다.\n"
        "- marker 간 불일치, 핵심 marker 미검출, 보조 marker만 검출, "
        "경계성 copy-number signal, 낮은 coverage, QC 저하, mosaic 가능성, "
        "비특이적 신호 또는 SOP 기준 미충족 중 해당 원인을 명시한다.\n"
        "- 실제 marker 결과와 SOP 문장을 직접 연결해 설명한다.\n"
        "- 모호하지 않은 결과이면 '해당 없음'으로 작성한다.\n\n"

        "### 5. 참고 근거\n"
        "- 실제 해석에 사용한 Reference ID, title, source만 나열한다.\n"
        "- 검색되었지만 사용하지 않은 문헌은 포함하지 않는다.\n"
        "<|im_end|>\n"

        "<|im_start|>assistant\n"
    )

# ── LLM 추론 ─────────────────────────────────────────────────────────
def run_llm(prompt: str, model_path: str) -> str:
    llm = get_llm(model_path)
    response = llm(
        prompt,
        max_tokens=MAX_TOKENS,
        temperature=0.1,
        repeat_penalty=1.1,
        echo=False,
    )
    return response["choices"][0]["text"].strip()


# ── 메인 ─────────────────────────────────────────────────────────────
def process_syndrome(syndrome: str, rows: list[dict], args) -> str:
    summary = summarize_syndrome(rows)

    print(f"\n{'─'*60}")
    print(f"[{syndrome}] overall={summary['overall_diagnosis']} "
          f"(+{summary['positive_count']} ?{summary['suspicious_count']})")

    if args.skip_negative and summary["overall_diagnosis"] == "NEGATIVE":
        print("  → NEGATIVE, skip")
        return ""

    contexts = retrieve(syndrome, summary, args.embed_model)
    print(f"  → {len(contexts)}개 RAG chunk 검색됨")

    prompt = build_prompt(summary, contexts)
    print(f"  → LLM 추론 중...")
    output = run_llm(prompt, args.model_path)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv",           required=True)
    parser.add_argument("--syndrome",      default=None)
    parser.add_argument("--all",           action="store_true")
    parser.add_argument("--skip-negative", action="store_true", default=True)
    parser.add_argument("--embed-model",   default=EMBED_MODEL)
    parser.add_argument("--model-path",    default=MODEL_PATH)
    parser.add_argument("--output-dir",    default="./reports")
    parser.add_argument("--threads",       type=int, default=32)
    args = parser.parse_args()

    N_THREADS = args.threads
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    data = load_tsv(args.tsv)
    print(f"✅ TSV 로드: {len(data)}개 syndrome")

    if args.syndrome:
        if args.syndrome not in data:
            print(f"❌ '{args.syndrome}' not found in TSV")
            return
        targets = {args.syndrome: data[args.syndrome]}
    elif args.all:
        targets = data
    else:
        targets = {
            s: rows for s, rows in data.items()
            if any(r["DIAGNOSIS"] in ("POSITIVE", "SUSPICIOUS") for r in rows)
        }
        print(f"   POSITIVE/SUSPICIOUS: {list(targets.keys())}")

    for syndrome, rows in targets.items():
        output = process_syndrome(syndrome, rows, args)
        if output:
            safe_name = syndrome.replace(" ", "_").replace("/", "-")
            out_path  = Path(args.output_dir) / f"{safe_name}_report.txt"
            out_path.write_text(output, encoding="utf-8")
            print(f"  💾 저장: {out_path}")
            print(f"\n{output[:300]}...\n")


if __name__ == "__main__":
    main()