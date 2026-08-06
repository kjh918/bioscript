"""
parsers/medgen_parser.py
─────────────────────────
MedGen HTML 파싱:
  - 질환 정의, synonyms, cross-refs (MONDO/OMIM/Orphanet/SNOMED)
  - Gene(location) 직접 추출
  - 카테고리별 논문 (etiology/diagnosis/therapy/prognosis/systematic_review)
  - GeneReview 발췌, OMIM 발췌
"""

import re
import time
from typing import Any

from bs4 import BeautifulSoup, Tag

from resource_parse_pipeline.utils import EUTILS, PUBMED, http_get

MEDGEN_BASE = "https://www.ncbi.nlm.nih.gov/medgen"

SUBHEAD_MAP = {
    "etiology":           "etiology",
    "diagnosis":          "diagnosis",
    "therapy":            "therapy",
    "prognosis":          "prognosis",
    "clinical prediction":"clinical_prediction",
}

CHROMOSOME_DISEASE_MAP = {
    # autosomal trisomy
    "down syndrome": "chr21",
    "patau syndrome": "chr13",
    "edwards syndrome": "chr18",

    # sex chromosome abnormalities
    "turner syndrome": "chrX",
    "klinefelter syndrome": "chrX",
    "jacobs syndrome": "chrY",
    "triple x syndrome": "chrX",
}


def _normalize_chromosome(value: str) -> str:
    """21, chr21, X 등을 chr21, chrX 형식으로 통일."""
    value = str(value).strip()
    value = re.sub(r"^chr", "", value, flags=re.I)

    if value.upper() in {"X", "Y"}:
        return f"chr{value.upper()}"

    if value.isdigit() and 1 <= int(value) <= 22:
        return f"chr{int(value)}"

    return ""


def _extract_chromosome_from_text(text: str) -> list[str]:
    """
    질환명 또는 synonym에서 chromosome 추출.

    지원 예:
      Trisomy 21
      chromosome 21
      chr21
      T21
      trisomy X
      monosomy X
      1p36 deletion
      22q11.2 deletion
    """
    if not text:
        return []

    found: list[str] = []
    lowered = text.lower().strip()

    # 대표 질환명 직접 매핑
    for disease_name, chromosome in CHROMOSOME_DISEASE_MAP.items():
        if disease_name in lowered:
            found.append(chromosome)

    patterns = [
        # Trisomy 21, Monosomy X, Chromosome 21
        r"\b(?:trisomy|monosomy|chromosome)\s*([1-9]|1\d|2[0-2]|x|y)\b",

        # T21, T18, T13
        r"\bT([1-9]|1\d|2[0-2])\b",

        # chr21
        r"\bchr([1-9]|1\d|2[0-2]|x|y)\b",

        # 1p36, 22q11.2
        r"(?<![A-Za-z0-9])([1-9]|1\d|2[0-2]|x|y)[pq]\d",
    ]

    for pattern in patterns:
        for match in re.finditer(pattern, text, flags=re.I):
            chromosome = _normalize_chromosome(match.group(1))
            if chromosome:
                found.append(chromosome)

    return list(dict.fromkeys(found))


def _extract_regions_from_text(text: str) -> list[dict[str, str]]:
    """
    1p36, 22q11.2, 8q23-q24.1 등의 cytoband 추출.
    """
    if not text:
        return []

    pattern = re.compile(
        r"(?<![A-Za-z0-9])"
        r"([1-9]|1\d|2[0-2]|X|Y)"
        r"([pq]\d+(?:\.\d+)?"
        r"(?:[-–](?:[1-9]|1\d|2[0-2]|X|Y)?[pq]?\d+(?:\.\d+)?)?)",
        flags=re.I,
    )

    regions: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for match in pattern.finditer(text):
        chromosome_number = match.group(1).upper()
        region_suffix = match.group(2).replace("–", "-")
        region_name = f"{chromosome_number}{region_suffix}"
        chromosome = _normalize_chromosome(chromosome_number)

        key = (chromosome, region_name.lower())
        if chromosome and key not in seen:
            seen.add(key)
            regions.append({
                "chromosome": chromosome,
                "region": region_name,
            })

    return regions

# ── ESearch ───────────────────────────────────────────────────

def search_medgen_uid(syndrome_name: str, email: str,
                      max_results: int = 1) -> list[str]:
    """질환명 → MedGen UID 목록 (relevance 순)"""
    r = http_get(f"{EUTILS}/esearch.fcgi",
                 params={"db": "medgen", "term": f'"{syndrome_name}"',
                         "retmax": max_results, "sort": "relevance",
                         "retmode": "json", "email": email})
    return r.json().get("esearchresult", {}).get("idlist", [])


# ── HTML 파싱 ─────────────────────────────────────────────────

def parse_medgen_page(html: str) -> dict[str, Any]:
    """
    MedGen HTML → 구조화된 dict.
    반환 필드:
      medgen_uid, definition, synonyms, cross_references,
      gene_locations, genereview_text, genereview_url,
      omim_description, literature
    """
    soup = BeautifulSoup(html, "lxml")
    result: dict[str, Any] = {}

    result["definition"]       = _parse_definition(soup)
    result["synonyms"]         = _parse_synonyms(soup)
    result["cross_references"] = _parse_xrefs(soup)
    result["gene_locations"]   = _parse_gene_locations(soup)
    gr_text, gr_url            = _parse_genereview(soup)
    result["genereview_text"]  = gr_text
    result["genereview_url"]   = gr_url
    result["omim_description"] = _parse_omim_description(soup)
    result["literature"]       = _parse_literature(soup)

    return result


def _parse_definition(soup: BeautifulSoup) -> str:
    portlet = soup.find("div", id="ID_100")
    if portlet:
        content = portlet.find("div", class_="portlet_content")
        if content:
            return re.sub(r"\s+", " ",
                          content.get_text(separator=" ", strip=True))
    meta = soup.find("meta", attrs={"name": "description"})
    if meta and meta.get("content"):
        return meta["content"].strip()
    return ""


def _parse_synonyms(soup: BeautifulSoup) -> list[str]:
    table = soup.find("table", class_="medgenTable")
    if not table:
        return []
    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) >= 2 and "synonym" in tds[0].get_text().lower():
            raw = tds[1].get_text(separator=";", strip=True)
            return [s.strip() for s in raw.split(";") if s.strip()]
    return []


def _parse_xrefs(soup: BeautifulSoup) -> dict[str, str]:
    xrefs: dict[str, str] = {}
    table = soup.find("table", class_="medgenTable")
    if not table:
        return xrefs
    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 2:
            continue
        label = tds[0].get_text().strip().rstrip(":").lower()
        a_tags = tds[1].find_all("a")
        first_a = a_tags[0] if a_tags else None
        val = (first_a.get_text(strip=True) if first_a
               else tds[1].get_text(strip=True))

        if "monarch" in label or "mondo" in label:
            xrefs["MONDO"] = val
        elif "omim" in label:
            xrefs["OMIM"] = val
        elif "orphanet" in label:
            xrefs["Orphanet"] = re.sub(r"^ORPHA", "", val)
        elif "snomed" in label:
            m = re.search(r"\((\d+)\)", val)
            xrefs["SNOMED-CT"] = m.group(1) if m else val
    return xrefs


def _parse_gene_locations(soup: BeautifulSoup) -> list[dict]:
    """medgenTable의 Gene(location) 행 파싱."""
    gene_locations: list[dict] = []
    table = soup.find("table", class_="medgenTable")
    if not table:
        return gene_locations
    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 2:
            continue
        label = tds[0].get_text().strip().rstrip(":").lower()
        if not ("gene" in label and "location" in label):
            continue
        for a in tds[1].find_all("a"):
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
                "source":       "medgen_direct",
            })
    return gene_locations


def _parse_genereview(soup: BeautifulSoup) -> tuple[str, str]:
    text, url = "", ""
    for portlet_id in ["ID_101", "ID_102"]:
        p = soup.find("div", id=portlet_id)
        if not p:
            continue
        content = p.find("div", class_="portlet_content")
        if content:
            t = re.sub(r"\s+", " ",
                       content.get_text(separator=" ", strip=True))
            if t and len(t) > 50:
                text += t + " "
        for a in p.find_all("a", href=True):
            if ("books.ncbi" in a["href"] or
                    "genereviews" in a["href"].lower()):
                if not url:
                    url = a["href"]
    return text.strip(), url


def _parse_omim_description(soup: BeautifulSoup) -> str:
    for h in soup.find_all(["h3", "h4", "strong", "b", "p"]):
        if "from omim" in h.get_text(strip=True).lower():
            parts = []
            node = h.find_next_sibling()
            while node and node.name not in ("h2", "h3", "h4"):
                if isinstance(node, Tag):
                    t = re.sub(r"\s+", " ",
                               node.get_text(separator=" ", strip=True))
                    if t:
                        parts.append(t)
                node = node.find_next_sibling()
            return " ".join(parts).strip()
    return ""


def _parse_literature(soup: BeautifulSoup) -> dict[str, Any]:
    """
    카테고리별 논문 파싱.
    반환: {cat: {"papers": [...], "see_all": {...}}}
    papers 항목: {pmid, title, authors, journal, year, pmc_free, url}
    """
    literature: dict[str, Any] = {}

    def _paper_from_blocks(nl_div: Tag,
                           detail_div: Tag | None) -> dict | None:
        a = nl_div.find("a")
        if not a:
            return None
        title = re.sub(r"\s+", " ", a.get_text(strip=True))
        if re.match(r"See all \(\d+\)", title):
            return None
        pmid = ""
        m = re.search(r"/pubmed/(\d+)", a.get("href", ""))
        if m:
            pmid = m.group(1)
        authors, journal, year = "", "", ""
        if detail_div:
            asp = detail_div.find("span", class_="medgenPMauthor")
            if asp:
                authors = re.sub(r"\s+", " ", asp.get_text(strip=True))
            jsp = detail_div.find("span", class_="medgenPMjournal")
            if jsp:
                journal = jsp.get_text(strip=True)
            ym = re.search(r"\b(19|20)\d{2}\b", detail_div.get_text())
            if ym:
                year = ym.group(0)
            if not pmid:
                pa = detail_div.find("a", attrs={"data-ga-action": "PMID"})
                if pa:
                    pmid = pa.get_text(strip=True)
        pmc_free = bool(detail_div and
                        detail_div.find("a", class_="PubMedFree"))
        return {
            "pmid":     pmid,
            "title":    title,
            "authors":  authors,
            "journal":  journal,
            "year":     year,
            "pmc_free": pmc_free,
            "url":      f"{PUBMED}/{pmid}/" if pmid else "",
        }

    def _see_all(container: Tag) -> dict:
        for a in container.find_all(
                "a", attrs={"data-ga-label": re.compile(r"See all", re.I)}):
            label = a.get("data-ga-label", "")
            m = re.search(r"See all \((\d+)\)", label)
            if m:
                href = a.get("href", "")
                if href.startswith("/"):
                    href = "https://www.ncbi.nlm.nih.gov" + href
                return {"count": int(m.group(1)), "url": href}
        return {}

    def _parse_portlet(portlet_id: str, default_cat: str):
        portlet = soup.find("div", id=portlet_id)
        if not portlet:
            return
        content = portlet.find("div", class_="portlet_content")
        if not content:
            return
        current = default_cat
        children = list(content.children)
        i = 0
        while i < len(children):
            node = children[i]
            if not isinstance(node, Tag):
                i += 1
                continue
            if node.name == "h3" and "subhead" in node.get("class", []):
                h3t = node.get_text(strip=True).lower()
                current = next(
                    (v for k, v in SUBHEAD_MAP.items() if k in h3t),
                    h3t.replace(" ", "_")
                )
                if current not in literature:
                    literature[current] = {"papers": [], "see_all": {}}
                i += 1
                continue
            if node.name == "div" and "nl" in node.get("class", []):
                sa = _see_all(node)
                if sa:
                    literature.setdefault(current, {"papers": [], "see_all": {}})
                    literature[current]["see_all"] = sa
                    i += 1
                    continue
                detail = None
                j = i + 1
                while j < len(children):
                    nxt = children[j]
                    if not isinstance(nxt, Tag):
                        j += 1
                        continue
                    if ("portlet_content" in nxt.get("class", []) and
                            "ln" in nxt.get("class", [])):
                        detail = nxt
                        i = j
                    break
                paper = _paper_from_blocks(node, detail)
                if paper:
                    literature.setdefault(
                        current, {"papers": [], "see_all": {}})
                    literature[current]["papers"].append(paper)
            i += 1

    # Professional guidelines (ID_105)
    portlet_105 = soup.find("div", id="ID_105")
    if portlet_105:
        content = portlet_105.find("div", class_="portlet_content")
        if content:
            papers = []
            children = list(content.children)
            i = 0
            while i < len(children):
                node = children[i]
                if not isinstance(node, Tag):
                    i += 1
                    continue
                if node.name == "div" and "nl" in node.get("class", []):
                    sa = _see_all(node)
                    if sa:
                        literature.setdefault(
                            "professional_guidelines",
                            {"papers": [], "see_all": {}})
                        literature["professional_guidelines"]["see_all"] = sa
                        i += 1
                        continue
                    detail = None
                    j = i + 1
                    while j < len(children):
                        nxt = children[j]
                        if not isinstance(nxt, Tag):
                            j += 1
                            continue
                        if ("portlet_content" in nxt.get("class", []) and
                                "ln" in nxt.get("class", [])):
                            detail = nxt
                            i = j
                        break
                    paper = _paper_from_blocks(node, detail)
                    if paper:
                        papers.append(paper)
                i += 1
            if papers:
                literature["professional_guidelines"] = {
                    "papers": papers, "see_all":
                    literature.get(
                        "professional_guidelines", {}).get("see_all", {})
                }

    # Recent clinical studies (ID_103)
    _parse_portlet("ID_103", "general")

    # Systematic reviews (ID_104)
    portlet_104 = soup.find("div", id="ID_104")
    if portlet_104:
        content = portlet_104.find("div", class_="portlet_content")
        if content:
            papers, sa = [], {}
            children = list(content.children)
            i = 0
            while i < len(children):
                node = children[i]
                if not isinstance(node, Tag):
                    i += 1
                    continue
                if node.name == "div" and "nl" in node.get("class", []):
                    s = _see_all(node)
                    if s:
                        sa = s
                        i += 1
                        continue
                    detail = None
                    j = i + 1
                    while j < len(children):
                        nxt = children[j]
                        if not isinstance(nxt, Tag):
                            j += 1
                            continue
                        if ("portlet_content" in nxt.get("class", []) and
                                "ln" in nxt.get("class", [])):
                            detail = nxt
                            i = j
                        break
                    paper = _paper_from_blocks(node, detail)
                    if paper:
                        papers.append(paper)
                i += 1
            if papers or sa:
                literature["systematic_reviews"] = {
                    "papers": papers, "see_all": sa}

    return literature


# ── 질환명 기반 feature 생성 ────────────────────────────────────

def _chromosome_from_cytoband(cytoband: str) -> str:
    """22q11.2, Xp21.1 등의 cytoband에서 chr명을 반환한다."""
    match = re.match(
        r"^(?:chr)?(\d{1,2}|X|Y|MT)",
        str(cytoband or "").strip(),
        re.I,
    )
    return f"chr{match.group(1).upper()}" if match else ""


def _extract_primary_targets(syndrome_name: str) -> dict[str, list[str]]:
    """
    질환명에 직접 명시된 chromosome/cytoband만 표적 locus로 추출한다.

    MedGen의 Gene(location)은 연관 유전자 정보이므로 TargetChromosome이나
    CoreRegion을 생성하는 근거로 사용하지 않는다.
    """
    text = str(syndrome_name or "").strip()
    chromosomes: list[str] = []
    regions: list[str] = []

    # 1p36, 22q11.2, 8q23-q24.1, Xp21.1
    region_pattern = re.compile(
        r"(?<![A-Za-z0-9])(?:chr)?(\d{1,2}|X|Y)"
        r"([pq]\d+(?:\.\d+)?(?:[-–][pq]?\d+(?:\.\d+)?)?)",
        re.I,
    )
    for match in region_pattern.finditer(text):
        chrom_token = match.group(1).upper()
        region = f"{chrom_token}{match.group(2).replace('–', '-')}"
        chromosome = f"chr{chrom_token}"
        if chromosome not in chromosomes:
            chromosomes.append(chromosome)
        if region not in regions:
            regions.append(region)

    # Trisomy 21, Monosomy X, Chromosome 18 syndrome 등
    if not chromosomes:
        whole_patterns = (
            r"\b(?:trisomy|monosomy)\s+(\d{1,2}|X|Y)\b",
            r"\bchromosome\s+(\d{1,2}|X|Y)\b",
        )
        for pattern in whole_patterns:
            match = re.search(pattern, text, re.I)
            if match:
                chromosomes.append(f"chr{match.group(1).upper()}")
                break

    return {"chromosomes": chromosomes, "regions": regions}
def build_discovered_features(
    medgen: dict[str, Any],
    syndrome_name: str,
) -> dict[str, list[dict]]:
    """
    MedGen 결과에서 annotation feature 생성.

    중요:
      - TargetChromosome과 CoreRegion은 질환명/synonym에서만 생성
      - gene_locations의 chromosome은 TargetChromosome으로 승격하지 않음
      - gene_locations는 CoreGene 생성에만 사용
    """
    discovered = {
        "target_chromosomes": [],
        "partial_chromosomes": [],
        "core_regions": [],
        "core_genes": [],
    }

    synonyms = medgen.get("synonyms", []) or []

    # syndrome과 synonym을 합쳐서 target 정보 탐색
    target_texts = [
        syndrome_name,
        medgen.get("syndrome", ""),
        *synonyms,
    ]

    target_chromosomes: list[str] = []
    region_candidates: list[dict[str, str]] = []

    for text in target_texts:
        target_chromosomes.extend(
            _extract_chromosome_from_text(str(text))
        )
        region_candidates.extend(
            _extract_regions_from_text(str(text))
        )

    target_chromosomes = list(dict.fromkeys(target_chromosomes))

    # CoreRegion에서 발견된 chromosome도 target에 포함
    for region in region_candidates:
        chromosome = region["chromosome"]
        if chromosome not in target_chromosomes:
            target_chromosomes.append(chromosome)

    for chromosome in target_chromosomes:
        discovered["target_chromosomes"].append({
            "feature_name": chromosome,
            "chromosome": chromosome,
            "start": 0,
            "end": 0,
            "size_mb": 0.0,
            "feature_type": "TargetChromosome",
            "source": "MedGen syndrome/synonym",
        })

    seen_regions: set[tuple[str, str]] = set()

    for region in region_candidates:
        key = (
            region["chromosome"],
            region["region"].lower(),
        )
        if key in seen_regions:
            continue

        seen_regions.add(key)
        discovered["core_regions"].append({
            "feature_name": region["region"],
            "chromosome": region["chromosome"],
            "start": 0,
            "end": 0,
            "size_mb": 0.0,
            "feature_type": "CoreRegion",
            "source": "MedGen syndrome/synonym",
        })

    # MedGen Gene(location)은 CoreGene에만 사용
    allowed_chromosomes = set(target_chromosomes)
    seen_genes: set[str] = set()

    for gene in medgen.get("gene_locations", []) or []:
        symbol = str(gene.get("gene_symbol", "")).strip()
        if not symbol:
            continue

        symbol_key = symbol.upper()
        if symbol_key in seen_genes:
            continue

        cytoband = str(gene.get("cytoband", "")).strip()
        gene_chromosome = ""

        chromosome_match = re.match(
            r"^(?:chr)?([1-9]|1\d|2[0-2]|X|Y)",
            cytoband,
            flags=re.I,
        )
        if chromosome_match:
            gene_chromosome = _normalize_chromosome(
                chromosome_match.group(1)
            )

        # 질환 표적 chromosome이 명확하면 다른 chromosome gene 제외
        if (
            allowed_chromosomes
            and gene_chromosome
            and gene_chromosome not in allowed_chromosomes
        ):
            continue

        seen_genes.add(symbol_key)

        discovered["core_genes"].append({
            "feature_name": symbol,
            "chromosome": gene_chromosome,
            "cytoband": cytoband,
            "start": 0,
            "end": 0,
            "size_mb": 0.0,
            "feature_type": "CoreGene",
            "ncbi_gene_id": gene.get("ncbi_gene_id", ""),
            "source": gene.get("source", "MedGen"),
        })

    return discovered
# ── 메인 수집 함수 ─────────────────────────────────────────────
def fetch_medgen(syndrome_name: str, email: str) -> dict[str, Any]:
    print(f"  [MedGen] '{syndrome_name}' 검색")

    ids = search_medgen_uid(syndrome_name, email)
    if not ids:
        raise ValueError(
            f"MedGen에서 '{syndrome_name}' 결과 없음"
        )

    uid = ids[0]
    url = f"{MEDGEN_BASE}/{uid}"
    print(f"  [MedGen] UID={uid}")

    time.sleep(1.5)

    from resource_parse_pipeline.utils import get_session

    resp = get_session().get(
        url,
        timeout=20,
        headers={"User-Agent": "Mozilla/5.0 Chrome/124"},
    )
    resp.raise_for_status()

    parsed = parse_medgen_page(resp.text)
    parsed["medgen_uid"] = uid
    parsed["syndrome"] = syndrome_name
    parsed["medgen_url"] = url

    # syndrome + synonyms 기반 feature 생성
    parsed["discovered_features"] = build_discovered_features(
        medgen=parsed,
        syndrome_name=syndrome_name,
    )

    targets = parsed["discovered_features"]["target_chromosomes"]
    regions = parsed["discovered_features"]["core_regions"]
    genes = parsed["discovered_features"]["core_genes"]

    print(
        "  [MedGen] "
        f"targets={[x['chromosome'] for x in targets]} "
        f"regions={[x['feature_name'] for x in regions]} "
        f"genes={[x['feature_name'] for x in genes]}"
    )

    return parsed