"""
medgen_llm_extractor.py
────────────────────────
Syndrome 이름 → MedGen HTML 파싱 → Local LLM 분석
→ 유전자/염색체 영역 구조화 JSON 출력

파이프라인:
  [1] MedGen ESearch → 최상위 UID 확보
  [2] MedGen HTML 파싱
        - definition, synonyms, cross-refs (MONDO/OMIM/Orphanet)
        - Gene(location) 직접 추출
        - GeneReview 발췌, OMIM 발췌
        - 카테고리별 논문 제목 목록 (etiology/diagnosis/therapy/prognosis)
  [3] 수집 텍스트 → Local LLM 프롬프트
        - 유전자 심볼 + 근거 추출
        - 염색체 영역 + 근거 추출
        - causative vs associated 분류
        - confidence (high/medium/low)
  [4] JSON 파싱 + 후처리 → 저장

Usage:
  python medgen_llm_extractor.py \\
      --syndrome "DiGeorge syndrome" \\
      --model Qwen/Qwen2.5-3B-Instruct \\
      --output output/digeorge.json

  python medgen_llm_extractor.py \\
      --syndrome "Down syndrome" \\
      --model microsoft/Phi-3-mini-4k-instruct \\
      --load-in-4bit \\
      --output output/down.json

모델 권장:
  GPU  8GB+  : Qwen/Qwen2.5-7B-Instruct / BioMistral/BioMistral-7B
  GPU  4-8GB : Qwen/Qwen2.5-3B-Instruct / microsoft/Phi-3-mini-4k-instruct
  CPU / RAM  : Qwen/Qwen2.5-3B-Instruct (--load-in-4bit)
               google/gemma-2-2b-it
"""

import re
import sys
import json
import time
import argparse
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from bs4 import BeautifulSoup

# ── 상수 ──────────────────────────────────────────────────────
EUTILS    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MEDGEN_BASE = "https://www.ncbi.nlm.nih.gov/medgen"

SESSION: requests.Session | None = None


def _make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=1.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://",  HTTPAdapter(max_retries=retry))
    s.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                      "AppleWebKit/537.36 Chrome/124.0.0.0 Safari/537.36"
    })
    return s


def _get(url: str, params: dict | None = None,
         delay: float = 0.35) -> requests.Response:
    time.sleep(delay)
    r = SESSION.get(url, params=params, timeout=20)
    r.raise_for_status()
    return r


# ════════════════════════════════════════════════════════════════
# Step 1+2. MedGen 파서
# ════════════════════════════════════════════════════════════════

def fetch_medgen_full(syndrome_name: str, email: str) -> dict[str, Any]:
    """
    MedGen에서 syndrome의 모든 텍스트 정보를 수집.
    LLM 입력에 최대한 많은 컨텍스트를 제공하기 위해
    모든 섹션을 raw text로 보존.
    """
    print(f"[MedGen] '{syndrome_name}' 검색 중...")

    # ESearch
    r = _get(f"{EUTILS}/esearch.fcgi",
             params={"db": "medgen", "term": f'"{syndrome_name}"',
                     "retmax": 1, "sort": "relevance",
                     "retmode": "json", "email": email})
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise ValueError(f"MedGen에서 '{syndrome_name}' 검색 결과 없음")

    uid = ids[0]
    url = f"{MEDGEN_BASE}/{uid}"
    print(f"  UID={uid}  URL={url}")

    time.sleep(1.5)
    resp = SESSION.get(url, timeout=20)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "lxml")

    result: dict[str, Any] = {
        "medgen_uid": uid,
        "syndrome": syndrome_name,
        "medgen_url": url,
    }

    # ── medgenTable 전체 파싱 ─────────────────────────────────
    table = soup.find("table", class_="medgenTable")
    xrefs: dict[str, str] = {}
    gene_locations: list[dict] = []
    synonyms: list[str] = []

    if table:
        for tr in table.find_all("tr"):
            tds = tr.find_all("td")
            if len(tds) < 2:
                continue
            label = tds[0].get_text(strip=True).rstrip(":").lower()
            a_tags = tds[1].find_all("a")
            first_a = a_tags[0] if a_tags else None
            val = (first_a.get_text(strip=True) if first_a
                   else tds[1].get_text(strip=True))

            if "synonym" in label:
                raw = tds[1].get_text(separator=";", strip=True)
                synonyms = [s.strip() for s in raw.split(";") if s.strip()]

            elif "gene" in label and "location" in label:
                for a in a_tags:
                    sym = a.get_text(strip=True)
                    if not sym:
                        continue
                    href = a.get("href", "")
                    gid_m = re.search(r"/gene/(\d+)", href)
                    cytoband = ""
                    for sib in a.next_siblings:
                        cb_m = re.search(
                            r"\(([0-9XY]+[pq][\d.]+(?:[- ][pq]?[\d.]+)?)\)",
                            str(sib))
                        if cb_m:
                            cytoband = cb_m.group(1).strip()
                            break
                    gene_locations.append({
                        "gene_symbol":  sym,
                        "ncbi_gene_id": gid_m.group(1) if gid_m else "",
                        "cytoband":     cytoband,
                    })

            elif "monarch" in label or "mondo" in label:
                xrefs["MONDO"] = val
            elif "omim" in label:
                xrefs["OMIM"] = val
            elif "orphanet" in label:
                xrefs["Orphanet"] = re.sub(r"^ORPHA", "", val)
            elif "snomed" in label:
                m = re.search(r"\((\d+)\)", val)
                xrefs["SNOMED-CT"] = m.group(1) if m else val

    result["synonyms"] = synonyms
    result["cross_references"] = xrefs
    result["gene_locations"] = gene_locations
    print(f"  Gene(location): {[g['gene_symbol'] for g in gene_locations]}")

    # ── Definition (portlet#ID_100) ───────────────────────────
    definition = ""
    portlet = soup.find("div", id="ID_100")
    if portlet:
        content = portlet.find("div", class_="portlet_content")
        if content:
            definition = re.sub(r"\s+", " ",
                                content.get_text(separator=" ", strip=True))
    if not definition:
        meta = soup.find("meta", attrs={"name": "description"})
        if meta and meta.get("content"):
            definition = meta["content"].strip()
    result["definition"] = definition

    # ── GeneReview 발췌 ──────────────────────────────────────
    genereview_text = ""
    genereview_url  = ""
    # Disease characteristics 섹션에서 발췌 텍스트 추출
    for portlet_id in ["ID_101", "ID_102"]:
        p = soup.find("div", id=portlet_id)
        if p:
            content = p.find("div", class_="portlet_content")
            if content:
                txt = re.sub(r"\s+", " ",
                             content.get_text(separator=" ", strip=True))
                if txt and len(txt) > 50:
                    genereview_text += txt + " "
            # GeneReview 링크
            for a in p.find_all("a", href=True):
                if "books.ncbi" in a["href"] or "genereviews" in a["href"].lower():
                    if not genereview_url:
                        genereview_url  = a["href"]

    # Additional descriptions (OMIM 발췌) 섹션
    for h in soup.find_all(["h3", "h4", "strong", "b", "p"]):
        txt_h = h.get_text(strip=True).lower()
        if "from omim" in txt_h or "omim description" in txt_h:
            parts = []
            node = h.find_next_sibling()
            while node and node.name not in ("h2","h3","h4"):
                t = re.sub(r"\s+", " ",
                           node.get_text(separator=" ", strip=True))
                if t:
                    parts.append(t)
                node = node.find_next_sibling()
            if parts:
                genereview_text += " ".join(parts)
            break

    result["genereview_text"] = genereview_text.strip()
    result["genereview_url"]  = genereview_url

    # ── 카테고리별 논문 제목 목록 ─────────────────────────────
    SUBHEAD_MAP = {
        "etiology": "etiology",
        "diagnosis": "diagnosis",
        "therapy": "therapy",
        "prognosis": "prognosis",
        "clinical prediction": "clinical_prediction",
    }
    literature: dict[str, list[str]] = {}

    # Professional guidelines (ID_105)
    for portlet_id, cat in [("ID_105", "professional_guidelines"),
                              ("ID_103", "clinical_studies"),
                              ("ID_104", "systematic_reviews")]:
        p = soup.find("div", id=portlet_id)
        if not p:
            continue
        content = p.find("div", class_="portlet_content")
        if not content:
            continue

        if portlet_id == "ID_103":
            # 서브헤딩별 분류
            current = "general"
            for child in content.children:
                if not isinstance(child, BeautifulSoup.__class__) and \
                   not hasattr(child, 'name'):
                    continue
                if not hasattr(child, 'name'):
                    continue
                if child.name == "h3" and "subhead" in child.get("class", []):
                    h3t = child.get_text(strip=True).lower()
                    current = next(
                        (v for k, v in SUBHEAD_MAP.items() if k in h3t),
                        h3t.replace(" ", "_")
                    )
                elif child.name == "div" and "nl" in child.get("class", []):
                    a = child.find("a")
                    if a:
                        title = re.sub(r"\s+", " ", a.get_text(strip=True))
                        see_m = re.match(r"See all \((\d+)\)", title)
                        if see_m:
                            literature.setdefault(
                                f"{current}_total", []
                            ).append(f"Total: {see_m.group(1)}")
                        elif title:
                            literature.setdefault(current, []).append(title)
        else:
            titles = []
            for a in content.find_all("a"):
                title = re.sub(r"\s+", " ", a.get_text(strip=True))
                see_m = re.match(r"See all \((\d+)\)", title)
                if see_m:
                    titles.append(f"[Total: {see_m.group(1)}]")
                elif title and len(title) > 10:
                    titles.append(title)
            if titles:
                literature[cat] = titles

    result["literature"] = literature

    # ── LLM용 컨텍스트 텍스트 조립 ──────────────────────────
    result["llm_context"] = _build_llm_context(result)

    n_papers = sum(len(v) for v in literature.values())
    print(f"  definition={'O' if definition else 'X'} "
          f"genereview={'O' if genereview_text else 'X'} "
          f"papers={n_papers}편")
    return result


def _build_llm_context(data: dict) -> str:
    """
    LLM 프롬프트에 넣을 컨텍스트 문자열.
    핵심 정보를 구조화된 섹션으로 구성.
    """
    lines = [f"=== {data['syndrome']} ===", ""]

    # 기본 정보
    if data.get("synonyms"):
        lines.append(f"Synonyms: {', '.join(data['synonyms'][:5])}")
    xr = data.get("cross_references", {})
    if xr:
        lines.append("External IDs: " +
                     ", ".join(f"{k}:{v}" for k, v in xr.items()))
    lines.append("")

    # Gene(location) - MedGen 직접 명시
    if data.get("gene_locations"):
        lines.append("Genes directly listed in MedGen:")
        for g in data["gene_locations"]:
            lines.append(f"  - {g['gene_symbol']} ({g['cytoband']})")
    lines.append("")

    # Definition
    if data.get("definition"):
        lines.append("Definition:")
        lines.append(data["definition"][:1500])
    lines.append("")

    # GeneReview/OMIM 발췌
    if data.get("genereview_text"):
        lines.append("Disease Description (from GeneReview/OMIM):")
        lines.append(data["genereview_text"][:2000])
    lines.append("")

    # 논문 제목 목록 (유전자/영역 단서 제공)
    lit = data.get("literature", {})
    if lit:
        lines.append("Related Literature Titles:")
        for cat, titles in lit.items():
            if not titles or "total" in cat:
                continue
            lines.append(f"  [{cat}]")
            for t in titles[:5]:
                lines.append(f"    - {t}")
    lines.append("")

    return "\n".join(lines)


# ════════════════════════════════════════════════════════════════
# Step 3. Local LLM 분석
# ════════════════════════════════════════════════════════════════

EXTRACTION_PROMPT = """\
You are a clinical genomics expert. Analyze the following disease information and extract ALL genes and chromosomal regions associated with this disease.

DISEASE INFORMATION:
{context}

TASK:
Extract all genes and chromosomal regions mentioned or implied in the text above.
For each item, provide:
1. gene_symbol OR region_name (e.g., "TBX1" or "22q11.21")
2. type: "gene" or "region"  
3. role: "causative" | "associated" | "modifier" | "candidate"
4. confidence: "high" | "medium" | "low"
   - high: directly stated as causative, haploinsufficiency proven
   - medium: mentioned in context of the disease, evidence exists
   - low: mentioned peripherally, candidate only
5. cytoband: chromosomal location (e.g., "22q11.21", "chr22")
6. evidence: brief text excerpt from the source that supports this (1-2 sentences max)

IMPORTANT RULES:
- Extract ONLY human genes (official HGNC symbols, e.g., TBX1, not Tbx1)
- Include chromosomal regions/bands (e.g., 22q11.2, 15q11-q13)
- Causative genes get confidence=high
- If a gene is mentioned as "critical", "key", "primary" → confidence=high
- Do NOT hallucinate genes not mentioned in the text

Respond ONLY with valid JSON, no explanation, no markdown:
{{
  "disease": "<disease name>",
  "genes": [
    {{
      "gene_symbol": "<SYMBOL>",
      "type": "gene",
      "role": "<role>",
      "confidence": "<confidence>",
      "cytoband": "<location>",
      "evidence": "<supporting text>"
    }}
  ],
  "regions": [
    {{
      "region_name": "<band notation>",
      "type": "region",
      "role": "<role>",
      "confidence": "<confidence>",
      "chromosome": "<chrN>",
      "evidence": "<supporting text>"
    }}
  ],
  "summary": "<1 sentence summary of key genomic findings>"
}}"""


class LLMExtractor:
    """
    HuggingFace transformers 기반 local LLM 래퍼.
    모델명만 바꾸면 교체 가능.
    """

    def __init__(self, model_name: str, load_in_4bit: bool = False,
                 device: str = "auto", max_new_tokens: int = 2048):
        self.model_name     = model_name
        self.max_new_tokens = max_new_tokens
        self._load(model_name, load_in_4bit, device)

    def _load(self, model_name: str, load_in_4bit: bool, device: str):
        print(f"[LLM] 모델 로드: {model_name}")
        from transformers import AutoTokenizer, AutoModelForCausalLM
        import torch

        kwargs: dict[str, Any] = {
            "trust_remote_code": True,
            "torch_dtype": torch.float16 if torch.cuda.is_available() else torch.float32,
        }

        if load_in_4bit:
            try:
                from transformers import BitsAndBytesConfig
                kwargs["quantization_config"] = BitsAndBytesConfig(
                    load_in_4bit=True,
                    bnb_4bit_compute_dtype=torch.float16,
                )
                print("  4bit 양자화 적용")
            except ImportError:
                print("  bitsandbytes 없음 → 일반 로드")

        if device == "auto":
            kwargs["device_map"] = "auto"

        self.tokenizer = AutoTokenizer.from_pretrained(
            model_name, trust_remote_code=True
        )
        self.model = AutoModelForCausalLM.from_pretrained(
            model_name, **kwargs
        )
        self.model.eval()
        print(f"  로드 완료")

    def _build_chat_prompt(self, user_message: str) -> str:
        """
        모델별 chat template 자동 적용.
        tokenizer.chat_template 있으면 사용, 없으면 직접 포맷.
        """
        messages = [
            {"role": "system",
             "content": "You are a clinical genomics expert. "
                        "Always respond with valid JSON only."},
            {"role": "user", "content": user_message},
        ]
        if hasattr(self.tokenizer, "apply_chat_template") and \
           self.tokenizer.chat_template:
            return self.tokenizer.apply_chat_template(
                messages,
                tokenize=False,
                add_generation_prompt=True,
            )
        # fallback
        return (f"<|system|>\nYou are a clinical genomics expert. "
                f"Always respond with valid JSON only.<|end|>\n"
                f"<|user|>\n{user_message}<|end|>\n<|assistant|>\n")

    def extract(self, context: str, syndrome: str) -> dict:
        """
        컨텍스트 텍스트 → LLM 분석 → 구조화 JSON 반환.
        """
        import torch

        prompt_text = EXTRACTION_PROMPT.format(context=context)
        full_prompt = self._build_chat_prompt(prompt_text)

        inputs = self.tokenizer(
            full_prompt,
            return_tensors="pt",
            truncation=True,
            max_length=4096,
        )
        input_len = inputs["input_ids"].shape[1]
        print(f"  [LLM] 입력 토큰: {input_len}")

        if torch.cuda.is_available():
            inputs = {k: v.cuda() for k, v in inputs.items()}

        with torch.no_grad():
            outputs = self.model.generate(
                **inputs,
                max_new_tokens=self.max_new_tokens,
                do_sample=False,          # greedy → JSON 안정성
                temperature=1.0,
                pad_token_id=self.tokenizer.eos_token_id,
                eos_token_id=self.tokenizer.eos_token_id,
            )

        # 새로 생성된 부분만 디코딩
        new_tokens = outputs[0][input_len:]
        raw_output = self.tokenizer.decode(new_tokens, skip_special_tokens=True)
        print(f"  [LLM] 출력 토큰: {len(new_tokens)}")
        print(f"  [LLM] raw output 앞 200자: {raw_output[:200]}")

        return self._parse_json_output(raw_output, syndrome)

    def _parse_json_output(self, raw: str, syndrome: str) -> dict:
        """
        LLM 출력에서 JSON 추출. 여러 방어 전략 적용.
        """
        # 1. 마크다운 코드블록 제거
        raw = re.sub(r"```json\s*", "", raw)
        raw = re.sub(r"```\s*", "", raw)
        raw = raw.strip()

        # 2. 첫 번째 { ... } 블록 추출
        start = raw.find("{")
        end   = raw.rfind("}")
        if start != -1 and end != -1 and end > start:
            json_str = raw[start:end+1]
            try:
                parsed = json.loads(json_str)
                print(f"  [LLM] JSON 파싱 성공: "
                      f"genes={len(parsed.get('genes',[]))} "
                      f"regions={len(parsed.get('regions',[]))}")
                return parsed
            except json.JSONDecodeError as e:
                print(f"  [LLM] JSON 파싱 실패: {e}")
                print(f"  [LLM] 시도한 JSON: {json_str[:300]}")

        # 3. 파싱 실패 → 빈 구조 + raw 보존
        print(f"  [LLM] JSON 추출 실패 → 빈 구조 반환")
        return {
            "disease": syndrome,
            "genes": [],
            "regions": [],
            "summary": "",
            "raw_output": raw[:2000],
            "parse_error": True,
        }


# ════════════════════════════════════════════════════════════════
# Step 4. 후처리 + 병합
# ════════════════════════════════════════════════════════════════

def postprocess(medgen_data: dict, llm_result: dict) -> dict:
    """
    MedGen 구조 + LLM 결과를 병합하고 정제.

    1. MedGen Gene(location) → 최우선 (source=medgen_direct)
    2. LLM 추출 결과 → source=llm_extracted
    3. 중복 제거 (gene_symbol 기준)
    4. MedGen에서 직접 확인된 항목은 confidence 보정
    """
    # MedGen 직접 명시 유전자 집합
    medgen_direct_genes = {
        g["gene_symbol"]: g
        for g in medgen_data.get("gene_locations", [])
    }

    genes: list[dict] = []
    seen_genes: set[str] = set()

    # MedGen direct 먼저
    for sym, info in medgen_direct_genes.items():
        genes.append({
            "gene_symbol":  sym,
            "type":         "gene",
            "role":         "causative",
            "confidence":   "high",
            "cytoband":     info.get("cytoband", ""),
            "ncbi_gene_id": info.get("ncbi_gene_id", ""),
            "evidence":     f"Directly listed in MedGen Gene(location) field "
                            f"for {medgen_data['syndrome']}",
            "source":       "medgen_direct",
        })
        seen_genes.add(sym)

    # LLM 추출 유전자 병합
    for g in llm_result.get("genes", []):
        sym = g.get("gene_symbol", "").strip()
        if not sym or not re.match(r"^[A-Z][A-Z0-9\-]{1,15}$", sym):
            continue  # 유효하지 않은 심볼 제거
        if sym in seen_genes:
            # 이미 있으면 evidence만 보완
            existing = next(x for x in genes if x["gene_symbol"] == sym)
            if g.get("evidence") and not existing.get("llm_evidence"):
                existing["llm_evidence"] = g["evidence"]
            # MedGen direct가 있으면 cytoband 보완
            if g.get("cytoband") and not existing.get("cytoband"):
                existing["cytoband"] = g["cytoband"]
            continue
        seen_genes.add(sym)
        genes.append({
            "gene_symbol":  sym,
            "type":         "gene",
            "role":         g.get("role", "associated"),
            "confidence":   g.get("confidence", "medium"),
            "cytoband":     g.get("cytoband", ""),
            "ncbi_gene_id": "",
            "evidence":     g.get("evidence", ""),
            "source":       "llm_extracted",
        })

    # 염색체 영역 (LLM만)
    regions: list[dict] = []
    seen_regions: set[str] = set()
    for r in llm_result.get("regions", []):
        name = r.get("region_name", "").strip()
        if not name or name in seen_regions:
            continue
        # 유효한 cytogenetic 표기법 검증
        if not re.match(r"^\d{1,2}[pq][\d.]+", name) and \
           not re.match(r"^[XY][pq][\d.]+", name):
            continue
        seen_regions.add(name)
        regions.append({
            "region_name": name,
            "type":        "region",
            "role":        r.get("role", "associated"),
            "confidence":  r.get("confidence", "medium"),
            "chromosome":  r.get("chromosome", ""),
            "evidence":    r.get("evidence", ""),
            "source":      "llm_extracted",
        })

    # confidence 정렬: high > medium > low
    conf_order = {"high": 0, "medium": 1, "low": 2}
    genes.sort(key=lambda x: conf_order.get(x["confidence"], 3))
    regions.sort(key=lambda x: conf_order.get(x["confidence"], 3))

    return {
        "syndrome":          medgen_data["syndrome"],
        "medgen_uid":        medgen_data["medgen_uid"],
        "medgen_url":        medgen_data["medgen_url"],
        "cross_references":  medgen_data.get("cross_references", {}),
        "synonyms":          medgen_data.get("synonyms", []),
        "definition":        medgen_data.get("definition", ""),
        "genereview_url":    medgen_data.get("genereview_url", ""),
        "llm_summary":       llm_result.get("summary", ""),
        "genes":             genes,
        "regions":           regions,
        "literature":        medgen_data.get("literature", {}),
        "stats": {
            "total_genes":         len(genes),
            "medgen_direct_genes": len(medgen_direct_genes),
            "llm_extracted_genes": len([g for g in genes if g["source"] == "llm_extracted"]),
            "high_confidence":     len([g for g in genes + regions if g["confidence"] == "high"]),
            "total_regions":       len(regions),
        },
        "collected_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
    }


# ════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="MedGen → Local LLM → 유전자/영역 구조화 추출",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
모델 권장:
  GPU  8GB+ : Qwen/Qwen2.5-7B-Instruct
  GPU  4-8GB: Qwen/Qwen2.5-3B-Instruct  또는  microsoft/Phi-3-mini-4k-instruct
  CPU       : Qwen/Qwen2.5-3B-Instruct --load-in-4bit  또는  google/gemma-2-2b-it

예시:
  python medgen_llm_extractor.py \\
      --syndrome "DiGeorge syndrome" \\
      --model Qwen/Qwen2.5-3B-Instruct \\
      --output output/digeorge.json

  python medgen_llm_extractor.py \\
      --syndrome "Down syndrome" \\
      --model Qwen/Qwen2.5-3B-Instruct \\
      --load-in-4bit \\
      --output output/down.json

  # 모든 19개 질환 일괄 처리
  python medgen_llm_extractor.py \\
      --syndrome-list syndromes.txt \\
      --model Qwen/Qwen2.5-3B-Instruct \\
      --output-dir output/
        """
    )
    parser.add_argument("--syndrome",     type=str, help="단일 syndrome 이름")
    parser.add_argument("--syndrome-list",type=str, help="syndrome 목록 파일 (줄당 1개)")
    parser.add_argument("--model",        required=True,
                        help="HuggingFace 모델 ID")
    parser.add_argument("--load-in-4bit", action="store_true",
                        help="4bit 양자화 (메모리 절약)")
    parser.add_argument("--device",       default="auto")
    parser.add_argument("--max-new-tokens", type=int, default=2048)
    parser.add_argument("--output",       type=str, default=None,
                        help="단일 출력 파일 (--syndrome 사용 시)")
    parser.add_argument("--output-dir",   type=str, default="output",
                        help="출력 디렉토리 (--syndrome-list 사용 시)")
    parser.add_argument("--email",        default="your@email.com")
    parser.add_argument("--skip-llm",     action="store_true",
                        help="LLM 없이 MedGen 파싱만 (디버그용)")
    args = parser.parse_args()

    if not args.syndrome and not args.syndrome_list:
        parser.error("--syndrome 또는 --syndrome-list 필수")

    # syndrome 목록 수집
    syndromes: list[str] = []
    if args.syndrome:
        syndromes = [args.syndrome]
    else:
        with open(args.syndrome_list, encoding="utf-8") as f:
            syndromes = [l.strip() for l in f if l.strip()]

    # 세션 초기화
    global SESSION
    SESSION = _make_session()

    # LLM 로드 (한 번만)
    extractor = None
    if not args.skip_llm:
        extractor = LLMExtractor(
            model_name=args.model,
            load_in_4bit=args.load_in_4bit,
            device=args.device,
            max_new_tokens=args.max_new_tokens,
        )

    # 출력 디렉토리
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for i, syndrome in enumerate(syndromes, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(syndromes)}] {syndrome}")
        print(f"{'='*60}")

        try:
            # Step 1+2: MedGen 파싱
            medgen_data = fetch_medgen_full(syndrome, args.email)

            # Step 3: LLM 분석
            if extractor:
                print(f"\n[LLM] 분석 중...")
                llm_result = extractor.extract(
                    medgen_data["llm_context"], syndrome
                )
            else:
                print("[LLM] skip (--skip-llm)")
                llm_result = {"disease": syndrome, "genes": [], "regions": [],
                              "summary": "LLM skipped"}

            # Step 4: 후처리 + 병합
            final = postprocess(medgen_data, llm_result)
            results.append(final)

            # 출력 파일 결정
            if args.output and len(syndromes) == 1:
                out_path = Path(args.output)
            else:
                safe_name = re.sub(r"[^a-z0-9]+", "_", syndrome.lower()).strip("_")
                out_path = out_dir / f"{safe_name}.json"

            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(final, f, ensure_ascii=False, indent=2)

            # 요약 출력
            stats = final["stats"]
            print(f"\n  ✓ {syndrome}")
            print(f"    genes  : {stats['total_genes']}개 "
                  f"(MedGen직접={stats['medgen_direct_genes']}, "
                  f"LLM추출={stats['llm_extracted_genes']})")
            print(f"    regions: {stats['total_regions']}개")
            print(f"    high_confidence: {stats['high_confidence']}개")
            print(f"    LLM summary: {final['llm_summary'][:80]}")
            print(f"    저장: {out_path}")

        except Exception as e:
            print(f"  [Error] {syndrome}: {e}")
            import traceback
            traceback.print_exc()

        time.sleep(1.0)

    # 전체 요약 (복수 질환)
    if len(syndromes) > 1:
        summary_path = out_dir / "summary.json"
        summary = [
            {
                "syndrome":     r["syndrome"],
                "medgen_uid":   r["medgen_uid"],
                "n_genes":      r["stats"]["total_genes"],
                "n_regions":    r["stats"]["total_regions"],
                "high_conf":    r["stats"]["high_confidence"],
                "top_genes":    [g["gene_symbol"] for g in r["genes"]
                                 if g["confidence"] == "high"][:5],
                "top_regions":  [rg["region_name"] for rg in r["regions"]
                                 if rg["confidence"] == "high"][:3],
            }
            for r in results
        ]
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        print(f"\n▶ 전체 요약: {summary_path}")


if __name__ == "__main__":
    main()