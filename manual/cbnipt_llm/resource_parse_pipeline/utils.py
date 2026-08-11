"""
utils.py — 공통 유틸리티

공유 상수, HTTP 세션, GRCh38 염색체 좌표
"""

import time
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from resource_parse_pipeline.reference import hg38


# ── API 엔드포인트 ─────────────────────────────────────────────
EUTILS   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUBMED   = "https://pubmed.ncbi.nlm.nih.gov"
PMC_BASE = "https://www.ncbi.nlm.nih.gov/pmc/articles"
UNIPROT  = "https://rest.uniprot.org/uniprotkb"

# ── GRCh38 전체 염색체 좌표 ────────────────────────────────────
# 출처: NCBI GRCh38.p14
GRCH38: dict[str, dict] = {
    "chr1":  {"start": 1, "end": 248956422},
    "chr2":  {"start": 1, "end": 242193529},
    "chr3":  {"start": 1, "end": 198295559},
    "chr4":  {"start": 1, "end": 190214555},
    "chr5":  {"start": 1, "end": 181538259},
    "chr6":  {"start": 1, "end": 170805979},
    "chr7":  {"start": 1, "end": 159345973},
    "chr8":  {"start": 1, "end": 145138636},
    "chr9":  {"start": 1, "end": 138394717},
    "chr10": {"start": 1, "end": 133797422},
    "chr11": {"start": 1, "end": 135086622},
    "chr12": {"start": 1, "end": 133275309},
    "chr13": {"start": 1, "end": 114364328},
    "chr14": {"start": 1, "end": 107043718},
    "chr15": {"start": 1, "end": 101991189},
    "chr16": {"start": 1, "end": 90338345},
    "chr17": {"start": 1, "end": 83257441},
    "chr18": {"start": 1, "end": 80373285},
    "chr19": {"start": 1, "end": 58617616},
    "chr20": {"start": 1, "end": 64444167},
    "chr21": {"start": 1, "end": 46709983},
    "chr22": {"start": 1, "end": 50818468},
    "chrX":  {"start": 1, "end": 156040895},
    "chrY":  {"start": 1, "end": 57227415},
}

# ── HTTP 세션 ─────────────────────────────────────────────────
_session: requests.Session | None = None
_email_global: str = "your@email.com"


def init_session(email: str) -> None:
    global _session, _email_global
    _email_global = email
    s = requests.Session()
    retry = Retry(
        total=5, backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://",  HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": f"nipt-pipeline/2.0 (mailto:{email})"})
    _session = s


def get_session() -> requests.Session:
    if _session is None:
        init_session(_email_global)
    return _session


def http_get(url: str, params: dict | None = None,
             delay: float = 0.35) -> requests.Response:
    time.sleep(delay)
    r = get_session().get(url, params=params, timeout=20)
    r.raise_for_status()
    return r


def chrom_coords(chrom: str) -> dict | None:
    """염색체 이름 → GRCh38 좌표 (없으면 None)"""
    key = chrom if chrom.startswith("chr") else f"chr{chrom.upper()}"
    coords = GRCH38.get(key)
    if coords:
        return {"chromosome": key, **coords, "assembly": "GRCh38"}
    return None