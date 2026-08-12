"""
parsers/uniprot_parser.py
──────────────────────────
UniProt 파싱 - 내부망 환경 대응 3가지 방법:

  Method 1 (권장): 로컬 텍스트 파일
    - https://www.uniprot.org/uniprotkb?query=gene_exact:{SYMBOL}&format=txt 직접 다운로드
    - parse_uniprot_txt(path) 로 파싱
    - 배치: uniprot_sprot_human.dat (Swiss-Prot 전체 human, ~200MB, 연 1-2회 갱신)

  Method 2: Biopython ExPASy
    - Bio.ExPASy.get_sprot_raw(accession) 사용
    - 내부망이라도 ExPASy 포트(443)가 열려있으면 작동

  Method 3: REST API (외부망 있을 때)
    - rest.uniprot.org REST API 직접 호출

사용 우선순위: 로컬 파일 → Biopython → REST API
"""

import re
import os
import gzip
import time
from io import StringIO
from pathlib import Path
from typing import Any

# ── UniProt 텍스트 포맷 파서 (Method 1/2 공통) ───────────────

def _parse_txt_entry(txt: str) -> dict[str, Any]:
    """
    UniProt .txt (Swiss-Prot flat file) 단일 entry 파싱.
    Bio.SwissProt.read()가 있으면 사용, 없으면 직접 파싱.
    """
    try:
        from Bio import SwissProt
        entry = SwissProt.read(StringIO(txt))
        return _from_biopython(entry)
    except ImportError:
        return _parse_txt_manual(txt)
    except Exception as e:
        print(f"  [UniProt] Biopython 파싱 실패, 수동 파싱으로: {e}")
        return _parse_txt_manual(txt)


def _from_biopython(entry) -> dict[str, Any]:
    """Bio.SwissProt entry → dict"""
    # 기능 설명
    function = ""
    for comment in entry.comments:
        if comment.startswith("FUNCTION:"):
            function = comment[len("FUNCTION:"):].strip()
            break

    # 도메인
    domains = []
    for feature in entry.features:
        if feature.type == "DOMAIN":
            domains.append({
                "name":  str(feature.qualifiers.get("note", [""])[0])
                         if hasattr(feature, "qualifiers") else "",
                "start": int(feature.location.start) + 1,
                "end":   int(feature.location.end),
            })

    # 질환 연관
    diseases = []
    for comment in entry.comments:
        if comment.startswith("DISEASE:"):
            raw = comment[len("DISEASE:"):].strip()
            mim = re.search(r"MIM:(\d+)", raw)
            diseases.append({
                "name":  re.match(r"([^[]+)", raw).group(1).strip()
                         if re.match(r"([^[]+)", raw) else raw[:80],
                "omim":  mim.group(1) if mim else "",
                "description": raw,
            })

    # Subcellular location
    subcell = []
    for comment in entry.comments:
        if comment.startswith("SUBCELLULAR LOCATION:"):
            raw = comment[len("SUBCELLULAR LOCATION:"):].strip()
            subcell = [s.strip() for s in raw.split(".") if s.strip()][:5]

    return {
        "uniprot_accession":    entry.accessions[0] if entry.accessions else "",
        "uniprot_id":           entry.entry_name,
        "protein_name":         entry.description,
        "gene_names":           entry.gene_name if hasattr(entry, "gene_name") else "",
        "organism":             entry.organism,
        "sequence_length":      entry.sequence_length,
        "protein_function":     function,
        "protein_domains":      domains,
        "disease_associations": diseases,
        "subcellular_location": subcell,
        "keywords":             list(entry.keywords)[:15],
        "uniprot_url":          "",  # accession으로 구성
    }


def _parse_txt_manual(txt: str) -> dict[str, Any]:
    """
    SwissProt flat file 수동 파싱.
    Biopython 없을 때 fallback.
    """
    lines = txt.splitlines()
    acc, entry_name = "", ""
    description, organism = "", ""
    function_lines: list[str] = []
    disease_lines:  list[str] = []
    subcell_lines:  list[str] = []
    keywords:       list[str] = []
    domains:        list[dict] = []
    seq_len = 0

    in_cc_function = False
    in_cc_disease  = False
    in_cc_subcell  = False

    for line in lines:
        tag = line[:2]
        val = line[5:].rstrip() if len(line) > 5 else ""

        if tag == "AC" and not acc:
            acc = val.split(";")[0].strip()
        elif tag == "ID" and not entry_name:
            entry_name = val.split()[0]
        elif tag == "DE":
            if "RecName: Full=" in val and not description:
                description = val.replace("RecName: Full=", "").strip(";").strip()
        elif tag == "OS" and not organism:
            organism = val.strip(".")
        elif tag == "KW":
            keywords.extend([k.strip().strip(";") for k in val.split(";")])
        elif tag == "SQ":
            m = re.search(r"(\d+) AA", val)
            if m:
                seq_len = int(m.group(1))
        elif tag == "FT":
            ft_match = re.match(r"(\S+)\s+(\d+)\.\.(\d+)", val)
            if ft_match and ft_match.group(1) == "DOMAIN":
                domains.append({
                    "name":  "",
                    "start": int(ft_match.group(2)),
                    "end":   int(ft_match.group(3)),
                })
            if "/note=" in val and domains:
                note = re.search(r'/note="([^"]+)"', val)
                if note:
                    domains[-1]["name"] = note.group(1)
        elif tag == "CC":
            if "-!- FUNCTION:" in val:
                in_cc_function = True
                in_cc_disease  = False
                in_cc_subcell  = False
                function_lines.append(val.replace("-!- FUNCTION:", "").strip())
            elif "-!- DISEASE:" in val:
                in_cc_disease  = True
                in_cc_function = False
                in_cc_subcell  = False
                disease_lines.append(val.replace("-!- DISEASE:", "").strip())
            elif "-!- SUBCELLULAR LOCATION:" in val:
                in_cc_subcell  = True
                in_cc_function = False
                in_cc_disease  = False
                subcell_lines.append(val.replace("-!- SUBCELLULAR LOCATION:", "").strip())
            elif val.startswith("-!-"):
                in_cc_function = in_cc_disease = in_cc_subcell = False
            elif in_cc_function:
                function_lines.append(val.strip())
            elif in_cc_disease:
                disease_lines.append(val.strip())
            elif in_cc_subcell:
                subcell_lines.append(val.strip())

    function_text = " ".join(function_lines).strip()
    disease_text  = " ".join(disease_lines).strip()

    # 질환 파싱 (DISEASE: Name [...])
    diseases = []
    for raw in re.split(r"(?=-!- DISEASE:)", disease_text):
        name_m = re.match(r"([A-Z][^[{]+)", raw.strip())
        mim_m  = re.search(r"MIM:(\d+)", raw)
        if name_m:
            diseases.append({
                "name":        name_m.group(1).strip(),
                "omim":        mim_m.group(1) if mim_m else "",
                "description": raw.strip()[:300],
            })

    return {
        "uniprot_accession":    acc,
        "uniprot_id":           entry_name,
        "protein_name":         description,
        "gene_names":           "",
        "organism":             organism,
        "sequence_length":      seq_len,
        "protein_function":     function_text,
        "protein_domains":      domains,
        "disease_associations": diseases,
        "subcellular_location": [s.strip() for s in
                                  " ".join(subcell_lines).split(".")
                                  if s.strip()][:5],
        "keywords":             [k for k in keywords if k][:15],
        "uniprot_url":          f"https://www.uniprot.org/uniprotkb/{acc}",
    }


# ── Method 1: 로컬 파일 ───────────────────────────────────────

class UniProtLocalDB:
    """
    Swiss-Prot human flat file (uniprot_sprot_human.dat[.gz]) 로드.
    gene_symbol → entry 매핑을 메모리에 캐싱.

    다운로드:
      wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/
           knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz
    """

    def __init__(self, dat_path: str):
        self.dat_path = dat_path
        self._gene_index: dict[str, str] | None = None  # symbol → accession
        self._entries:    dict[str, dict]         = {}   # accession → parsed

    def _load_index(self):
        """GN (Gene Name) 라인에서 gene_symbol → accession 인덱스 생성."""
        if self._gene_index is not None:
            return
        print(f"  [UniProt] 인덱스 빌드 중: {self.dat_path}")
        index: dict[str, str] = {}
        current_acc: list[str] = []
        open_fn = gzip.open if self.dat_path.endswith(".gz") else open
        with open_fn(self.dat_path, "rt", encoding="utf-8",
                     errors="ignore") as f:
            for line in f:
                tag = line[:2]
                val = line[5:].rstrip()
                if tag == "AC":
                    if not current_acc:
                        current_acc = [v.strip() for v in
                                       val.split(";") if v.strip()]
                elif tag == "GN":
                    # GN   Name=TBX1; Synonyms=...
                    for m in re.finditer(r"Name=([A-Z0-9\-]+)", val):
                        sym = m.group(1)
                        if current_acc:
                            index[sym] = current_acc[0]
                    for m in re.finditer(r"Synonyms=([^;]+)", val):
                        for s in m.group(1).split(","):
                            s = s.strip()
                            if s and current_acc:
                                index.setdefault(s, current_acc[0])
                elif tag == "//":
                    current_acc = []
        self._gene_index = index
        print(f"  [UniProt] 인덱스 완료: {len(index)}개 심볼")

    def _load_entry_by_acc(self, accession: str) -> dict:
        """accession으로 entry 텍스트 추출 후 파싱."""
        if accession in self._entries:
            return self._entries[accession]
        open_fn = gzip.open if self.dat_path.endswith(".gz") else open
        entry_lines: list[str] = []
        in_entry = False
        with open_fn(self.dat_path, "rt", encoding="utf-8",
                     errors="ignore") as f:
            for line in f:
                tag = line[:2]
                val = line[5:].rstrip()
                if tag == "AC" and accession in val and not in_entry:
                    in_entry = True
                if in_entry:
                    entry_lines.append(line.rstrip())
                    if tag == "//":
                        break
        if not entry_lines:
            return {}
        result = _parse_txt_entry("\n".join(entry_lines))
        result["uniprot_url"] = (
            f"https://www.uniprot.org/uniprotkb/{accession}")
        self._entries[accession] = result
        return result

    def get_by_gene(self, gene_symbol: str) -> dict:
        """gene_symbol → UniProt entry dict"""
        self._load_index()
        acc = self._gene_index.get(gene_symbol, "")
        if not acc:
            return {}
        return self._load_entry_by_acc(acc)


# ── Method 2: Biopython ExPASy ───────────────────────────────

def fetch_uniprot_expasy(gene_symbol: str) -> dict:
    """
    Biopython Bio.ExPASy + Bio.SwissProt 사용.
    ExPASy 서버(https://www.expasy.org)에서 직접 조회.
    내부망에서도 443 포트가 열려있으면 작동.
    """
    try:
        from Bio import ExPASy, SwissProt
        # gene name으로 accession 검색 (ExPASy search)
        # 직접 accession을 모르면 먼저 esearch로 accession 확인 필요
        # 여기서는 gene_symbol이 accession처럼 동작하는 경우를 처리
        handle = ExPASy.get_sprot_raw(gene_symbol)
        record = SwissProt.read(handle)
        result = _from_biopython(record)
        result["uniprot_url"] = (
            f"https://www.uniprot.org/uniprotkb/{result.get('uniprot_accession','')}")
        return result
    except Exception as e:
        print(f"  [UniProt/ExPASy] {gene_symbol} 실패: {e}")
        return {}


def search_uniprot_accession_via_ncbi(gene_symbol: str,
                                       email: str) -> str:
    """
    NCBI Gene → UniProt accession 매핑 (elink 또는 ESummary dblinks).
    외부 UniProt API 없이도 accession을 얻는 우회 방법.
    """
    try:
        from nipt_pipeline.utils import http_get, EUTILS
        # ESearch: gene_symbol → NCBI Gene ID
        r = http_get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "gene",
                             "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
                             "retmax": 1, "retmode": "json", "email": email})
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return ""
        gene_id = ids[0]

        # ESummary: gene_id → UniProt accession (dblinks 필드)
        time.sleep(0.35)
        r2 = http_get(f"{EUTILS}/esummary.fcgi",
                      params={"db": "gene", "id": gene_id,
                              "retmode": "json", "email": email})
        data = r2.json().get("result", {}).get(gene_id, {})
        # UniProt accession은 nomenclatureauthorityid 또는 별도 필드
        for xref in data.get("gene_other_alias", []):
            if re.match(r"[OPQ][0-9][A-Z0-9]{3}[0-9]", str(xref)):
                return xref
    except Exception as e:
        print(f"  [UniProt/NCBI] {gene_symbol} accession 조회 실패: {e}")
    return ""


# ── Method 3: REST API ────────────────────────────────────────

def fetch_uniprot_rest(gene_symbol: str) -> dict:
    """
    UniProt REST API (외부망 있을 때).
    rest.uniprot.org
    """
    try:
        from nipt_pipeline.utils import http_get, UNIPROT
        query = (f"gene_exact:{gene_symbol} AND "
                 f"organism_id:9606 AND reviewed:true")
        r = http_get(UNIPROT,
                     params={"query": query, "format": "json", "size": 1})
        results = r.json().get("results", [])
        if not results:
            return {}
        entry = results[0]
        acc = entry.get("primaryAccession", "")

        function = next(
            (c["texts"][0]["value"]
             for c in entry.get("comments", [])
             if c.get("commentType") == "FUNCTION" and c.get("texts")),
            "")

        domains = [
            {"name":  f.get("description", ""),
             "start": f.get("location", {}).get("start", {}).get("value", ""),
             "end":   f.get("location", {}).get("end",   {}).get("value", "")}
            for f in entry.get("features", [])
            if f.get("type") == "Domain"
        ]

        diseases = [
            {"name":        c.get("disease", {}).get("diseaseId", ""),
             "omim":        c.get("disease", {}).get(
                                "diseaseCrossReference", {}).get("id", ""),
             "description": (c.get("texts") or [{}])[0].get("value", "")}
            for c in entry.get("comments", [])
            if c.get("commentType") == "DISEASE"
        ]

        subcell = [
            loc.get("location", {}).get("value", "")
            for c in entry.get("comments", [])
            if c.get("commentType") == "SUBCELLULAR LOCATION"
            for loc in c.get("subcellularLocations", [])
            if loc.get("location", {}).get("value")
        ]

        return {
            "uniprot_accession":    acc,
            "uniprot_id":           entry.get("uniProtkbId", ""),
            "protein_name":         (entry.get("proteinDescription", {})
                                     .get("recommendedName", {})
                                     .get("fullName", {}).get("value", "")),
            "gene_names":           gene_symbol,
            "organism":             "Homo sapiens",
            "sequence_length":      0,
            "protein_function":     function,
            "protein_domains":      domains,
            "disease_associations": diseases,
            "subcellular_location": subcell,
            "keywords":             [kw.get("name", "")
                                     for kw in entry.get("keywords", [])][:15],
            "uniprot_url":          f"https://www.uniprot.org/uniprotkb/{acc}",
        }
    except Exception as e:
        print(f"  [UniProt/REST] {gene_symbol} 실패: {e}")
        return {}


# ── 통합 진입점 ───────────────────────────────────────────────

def fetch_uniprot(gene_symbol: str,
                  local_db: "UniProtLocalDB | None" = None,
                  email: str = "") -> dict:
    """
    우선순위: 로컬 파일 DB → Biopython ExPASy → REST API

    Parameters
    ----------
    gene_symbol : str
    local_db    : UniProtLocalDB 인스턴스 (로컬 파일 방식)
    email       : NCBI API 이메일 (accession 조회 우회에 사용)
    """
    # 1. 로컬 파일 DB
    if local_db is not None:
        result = local_db.get_by_gene(gene_symbol)
        if result:
            print(f"    [UniProt] {gene_symbol} → 로컬 파일")
            return result

    # 2. Biopython ExPASy
    result = fetch_uniprot_expasy(gene_symbol)
    if result:
        print(f"    [UniProt] {gene_symbol} → ExPASy")
        return result

    # 3. REST API
    result = fetch_uniprot_rest(gene_symbol)
    if result:
        print(f"    [UniProt] {gene_symbol} → REST API")
        return result

    print(f"    [UniProt] {gene_symbol} → 모든 방법 실패, 빈 dict 반환")
    return {}