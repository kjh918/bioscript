"""
utils.py
─────────
공통 유틸리티: config 로딩, HTTP 세션, 염색체 좌표 헬퍼

references/
  hg38.json        → 염색체 좌표 (p_arm/centromere/q_arm)
  api_config.json  → API 엔드포인트, HTTP 설정,
                     로컬 파일 경로 (local_files),
                     원격 다운로드 URL (remote_urls)

local_files 우선 사용: 경로가 있으면 로컬 파일, 없으면 remote_urls 또는 API 직접 호출.
"""

import json
import time
import pathlib
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ── references 경로 ────────────────────────────────────────────
_REFS_DIR = pathlib.Path(__file__).parent / "references"

def _load_json(path: pathlib.Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)

_api  = _load_json(_REFS_DIR / "api_config.json")
_hg38 = _load_json(_REFS_DIR / "hg38.json")

# ── API 엔드포인트 ─────────────────────────────────────────────
EUTILS      = _api["ncbi"]["eutils"]
PUBMED      = _api["ncbi"]["pubmed"]
PMC_BASE    = _api["ncbi"]["pmc_base"]
UNIPROT     = _api["uniprot"]["rest"]
MEDGEN_BASE = _api["ncbi"]["medgen"]
HTTP_CFG    = _api["http"]

# ── 로컬 파일 / 원격 URL ──────────────────────────────────────
LOCAL_FILES:  dict[str, str] = _api.get("local_files",  {})
REMOTE_URLS:  dict[str, str] = _api.get("remote_urls",  {})


def get_file_path(key: str) -> str | None:
    """
    key에 해당하는 로컬 파일 경로 반환.
    local_files[key]가 비어있거나 파일이 없으면 None.

    Examples
    --------
    get_file_path("clingen_gene")   → "/data/ClinGen_gene_curation_list_GRCh38.tsv"
    get_file_path("gencc")          → None  (경로 미설정)
    """
    path_str = LOCAL_FILES.get(key, "").strip()
    if not path_str:
        return None
    p = pathlib.Path(path_str)
    if p.exists():
        return str(p)
    print(f"  [config] local_files[{key!r}]={path_str!r} 파일 없음 → None")
    return None


def get_remote_url(key: str) -> str | None:
    """key에 해당하는 원격 다운로드 URL 반환."""
    return REMOTE_URLS.get(key) or None


def resolve_path(key: str) -> tuple[str | None, str]:
    """
    로컬 파일 우선, 없으면 remote_url 반환.

    Returns
    -------
    (path_or_url, source)
      source = "local" | "remote" | "none"

    Examples
    --------
    path, src = resolve_path("clingen_gene")
    if src == "local":
        df = parse_clingen_gene(path)
    elif src == "remote":
        # 다운로드 필요
        print(f"wget {path}")
    """
    local = get_file_path(key)
    if local:
        return local, "local"
    remote = get_remote_url(key)
    if remote:
        return remote, "remote"
    return None, "none"


# ── hg38 좌표 ─────────────────────────────────────────────────
ASSEMBLY       = _hg38["assembly"]
CHROMS: dict[str, dict] = _hg38["chromosomes"]


def chrom_info(chrom: str) -> dict | None:
    """
    염색체 이름 → {p_arm, centromere, q_arm, start, end, assembly}

    자동 변환: "21" → "chr21", "X" → "chrX"
    없는 염색체 → None
    """
    key = chrom if chrom.startswith("chr") else f"chr{chrom.upper()}"
    c = CHROMS.get(key)
    if not c:
        return None
    return {
        **c,
        "start":    c["p_arm"]["start"],
        "end":      c["q_arm"]["end"],
        "assembly": ASSEMBLY,
    }


def arm_coords(chrom: str, arm: str) -> dict | None:
    """
    arm = 'p' | 'q' | 'centromere' → {start, end}
    """
    info = chrom_info(chrom)
    if not info:
        return None
    key = {"p": "p_arm", "q": "q_arm", "centromere": "centromere"}.get(arm)
    return info.get(key)


def full_chrom_pos(chrom: str) -> dict | None:
    """
    TargetChromosome 마커용: 전체 염색체 좌표.
    → {"chromosome": "chr21", "start": 1, "end": 46709983, "assembly": "GRCh38"}
    """
    info = chrom_info(chrom)
    if not info:
        return None
    key = chrom if chrom.startswith("chr") else f"chr{chrom.upper()}"
    return {
        "chromosome": key,
        "start":      info["start"],
        "end":        info["end"],
        "assembly":   ASSEMBLY,
    }


# ── HTTP 세션 ─────────────────────────────────────────────────
_session: requests.Session | None = None
_email_global: str = "your@email.com"


def init_session(email: str) -> None:
    global _session, _email_global
    _email_global = email
    s = requests.Session()
    retry = Retry(
        total            = HTTP_CFG["retry_total"],
        backoff_factor   = HTTP_CFG["retry_backoff"],
        status_forcelist = HTTP_CFG["retry_status"],
        allowed_methods  = ["GET"],
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
             delay: float | None = None) -> requests.Response:
    d = delay if delay is not None else HTTP_CFG["delay_default"]
    time.sleep(d)
    r = get_session().get(url, params=params,
                          timeout=HTTP_CFG["timeout"])
    r.raise_for_status()
    return r

# ── hg38 cytoband 좌표 ──────────────────────────────────────────

def cytoband_pos(chrom: str, band: str) -> dict | None:
    """
    cytoband → GRCh38 genomic position.

    hg38.json 구조 예:
        "chr22": {
            "p_arm": {...},
            "centromere": {...},
            "q_arm": {...},
            "cytobands": {
                "q11.21": {
                    "start": 18100001,
                    "end": 21400000,
                    "stain": "gneg"
                },
                ...
            }
        }

    지원:
        cytoband_pos("chr22", "q11.21")
        cytoband_pos("22", "22q11.21")
        cytoband_pos("chr22", "22q11.2")

    정확한 band가 없고 상위 band가 들어오면 하위 band 범위를 병합:
        q11.2
          → q11.21 + q11.22 + q11.23
          → min(start) ~ max(end)

    Returns
    -------
    {
        "chromosome": "chr22",
        "cytoband": "22q11.2",
        "start": ...,
        "end": ...,
        "assembly": "GRCh38"
    }

    찾지 못하면 None.
    """

    # chromosome 이름 정규화
    chrom_key = chrom if chrom.startswith("chr") else f"chr{chrom.upper()}"

    chrom_data = CHROMS.get(chrom_key)
    if not chrom_data:
        return None

    bands = chrom_data.get("cytobands", {})
    if not bands:
        return None

    # ---------------------------------------------------------
    # band 정규화
    #
    # "22q11.21" → "q11.21"
    # "chr22q11.21" → "q11.21"
    # "q11.21" → "q11.21"
    # ---------------------------------------------------------
    normalized_band = str(band).strip()

    normalized_band = normalized_band.removeprefix("chr")

    chrom_name = chrom_key.removeprefix("chr")

    if normalized_band.upper().startswith(chrom_name.upper()):
        normalized_band = normalized_band[len(chrom_name):]

    normalized_band = normalized_band.strip()

    if not normalized_band:
        return None

    # ---------------------------------------------------------
    # 1. 정확한 cytoband
    # ---------------------------------------------------------
    exact = bands.get(normalized_band)

    if exact:
        return {
            "chromosome": chrom_key,
            "cytoband": f"{chrom_name}{normalized_band}",
            "start": exact["start"],
            "end": exact["end"],
            "assembly": ASSEMBLY,
        }

    # ---------------------------------------------------------
    # 2. 상위 cytoband → 하위 band 병합
    #
    # q11.2
    #   q11.21
    #   q11.22
    #   q11.23
    # ---------------------------------------------------------
    matched = [
        info
        for band_name, info in bands.items()
        if band_name.startswith(normalized_band)
    ]

    if not matched:
        return None

    start = min(x["start"] for x in matched)
    end = max(x["end"] for x in matched)

    return {
        "chromosome": chrom_key,
        "cytoband": f"{chrom_name}{normalized_band}",
        "start": start,
        "end": end,
        "assembly": ASSEMBLY,
    }


def cytoband_coords(cytoband: str) -> dict | None:
    """
    chromosome이 포함된 cytoband 문자열 하나만 받아 좌표 반환.

    Examples
    --------
    cytoband_coords("22q11.21")
    cytoband_coords("22q11.2")
    cytoband_coords("1p36.33")
    cytoband_coords("Xp22.31")

    Returns
    -------
    {
        "chromosome": "chr22",
        "cytoband": "22q11.21",
        "start": ...,
        "end": ...,
        "assembly": "GRCh38"
    }
    """

    import re

    value = str(cytoband).strip()

    # chr22q11.21 / 22q11.21 / Xp22.31
    m = re.match(
        r"^(?:chr)?(\d+|X|Y)([pq].+)$",
        value,
        re.IGNORECASE,
    )

    if not m:
        return None

    chrom = m.group(1).upper()
    band = m.group(2)

    return cytoband_pos(
        chrom=f"chr{chrom}",
        band=band,
    )
