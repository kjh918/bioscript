"""
parsers/medgen_parser.py
─────────────────────────
MedGen HTML 파싱 → disease_info JSON

INPUT : syndrome_name (str)
OUTPUT: disease_info dict
  {
    medgen_uid, syndrome, medgen_url,
    definition, synonyms,
    db_ids: {MONDO, OMIM, Orphanet, SNOMED-CT},
    gene_locations: [{gene_symbol, ncbi_gene_id, cytoband}],
    genomic_targets: [{type, chromosome, cytoband, event, mechanism}],
    genereview_url,
    literature: {
      etiology:   {papers: [{pmid, title, authors, journal, year, url}], see_all: {}},
      diagnosis:  {...},
      therapy:    {...},
      prognosis:  {...},
      professional_guidelines: {...},
      systematic_reviews: {...},
    }
  }

단독 실행:
  python -m nipt_pipeline.parsers.medgen_parser \\
      --syndrome "Down syndrome" \\
      --output disease_info.json
"""

import re
import sys
import json
import time
import argparse
from typing import Any

from bs4 import BeautifulSoup, Tag

from nipt_pipeline.utils import EUTILS, PUBMED, http_get, get_session

MEDGEN_BASE = "https://www.ncbi.nlm.nih.gov/medgen"

SUBHEAD_MAP = {
    "etiology":           "etiology",
    "diagnosis":          "diagnosis",
    "therapy":            "therapy",
    "prognosis":          "prognosis",
    "clinical prediction":"clinical_prediction",
}

# 염색체 수 이상 패턴
_GAIN_PAT = re.compile(r"\b(trisomy|duplication|duplicated|gain|tetrasomy)\b", re.I)
_LOSS_PAT = re.compile(r"\b(monosomy|deletion|deleted|loss)\b", re.I)
_TRISOMY_PAT = re.compile(
    r"\b(?:trisomy|monosomy|tetrasomy)\s*(?:of\s+chromosome\s*)?(\d{1,2}|x|y)\b", re.I)
_BAND_PAT = re.compile(
    r"(?<![a-z0-9])"
    r"((?:[1-9]|1\d|2[0-2]|x|y)[pq](?:\d+(?:\.\d+)?)?(?:\s*[-–]\s*(?:(?:[1-9]|1\d|2[0-2]|x|y)?[pq]?\d+(?:\.\d+)?))?)"
    r"(?![a-z0-9])", re.I)


def _event(text: str) -> str:
    if _GAIN_PAT.search(text):
        return "gain"
    if _LOSS_PAT.search(text):
        return "loss"
    return "unknown"


def extract_genomic_targets(syndrome_name: str,
                             synonyms: list[str],
                             definition: str) -> list[dict]:
    """
    질환명·동의어·정의 텍스트에서 염색체/band target 추출.
    position 정보는 여기서 채우지 않음 (chromosome_annotator에서 처리).
    """
    combined = " ; ".join(filter(None, [syndrome_name, *synonyms, definition]))
    low = combined.lower()
    ev  = _event(low)
    targets: dict[tuple, dict] = {}

    # trisomy/monosomy + 숫자
    for m in _TRISOMY_PAT.finditer(low):
        chrom = m.group(1).upper()
        key = (f"chr{chrom}", chrom, ev)
        targets[key] = {
            "type":       "TargetChromosome",
            "chromosome": f"chr{chrom}",
            "cytoband":   chrom,
            "event":      "gain" if "trisomy" in m.group(0).lower() else "loss",
            "mechanism":  m.group(0).split()[0].lower(),
            "source":     "name_or_synonym",
        }

    # sex chromosome 패턴: 45,X / 47,XXY / 47,XYY
    for pat, chroms, mech, ev2 in [
        (r"45\s*,?\s*x(?!\s*[xy])\b", ["chrX"], "monosomy", "loss"),
        (r"47\s*,?\s*xxy\b",           ["chrX", "chrY"], "trisomy", "gain"),
        (r"47\s*,?\s*xyy\b",           ["chrY"], "trisomy", "gain"),
        (r"47\s*,?\s*xxx\b",           ["chrX"], "trisomy", "gain"),
    ]:
        if re.search(pat, low):
            for chrom in chroms:
                key = (chrom, chrom.replace("chr",""), ev2)
                targets[key] = {
                    "type":       "TargetChromosome",
                    "chromosome": chrom,
                    "cytoband":   chrom.replace("chr",""),
                    "event":      ev2,
                    "mechanism":  mech,
                    "source":     "karyotype_notation",
                }

    # cytogenetic band (22q11.2, 15q11-q13 등)
    for m in _BAND_PAT.finditer(combined, re.I):
        band = re.sub(r"\s+", "", m.group(1)).replace("–", "-")
        cm = re.match(r"(\d{1,2}|X|Y)", band, re.I)
        if not cm:
            continue
        chrom_num = cm.group(1).upper()
        key = (f"chr{chrom_num}", band.upper(), ev)
        # 이미 TargetChromosome으로 잡은 건 더 구체적인 CoreRegion으로 업데이트
        targets[key] = {
            "type":       "CoreRegion",
            "chromosome": f"chr{chrom_num}",
            "cytoband":   band,
            "event":      ev,
            "mechanism":  "deletion" if ev == "loss" else "duplication",
            "source":     "name_synonym_or_definition",
        }

    # 더 구체적인 band가 있으면 같은 염색체의 whole-chrom은 제거
    specific_chroms = {v["chromosome"] for v in targets.values()
                       if v["type"] == "CoreRegion"}
    result = [v for v in targets.values()
              if not (v["type"] == "TargetChromosome"
                      and v["chromosome"] in specific_chroms)]
    return result


def _safe_text(el: Tag | None) -> str:
    return re.sub(r"\s+", " ", el.get_text(separator=" ", strip=True)) if el else ""


def _parse_paper(nl_div: Tag, detail_div: Tag | None) -> dict | None:
    """div.nl(제목) + div.portlet_content.ln(메타) → paper dict"""
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

    return {
        "pmid":    pmid,
        "title":   title,
        "authors": authors,
        "journal": journal,
        "year":    year,
        "url":     f"{PUBMED}/{pmid}/" if pmid else "",
    }


def _walk_portlet(content: Tag, default_cat: str) -> dict[str, Any]:
    """portlet_content → {cat: {papers, see_all}} 파싱"""
    lit: dict[str, Any] = {}
    current = default_cat
    children = list(content.children)
    i = 0
    while i < len(children):
        node = children[i]
        if not isinstance(node, Tag):
            i += 1
            continue

        # 서브헤딩
        if node.name == "h3" and "subhead" in node.get("class", []):
            h3t = node.get_text(strip=True).lower()
            current = next((v for k, v in SUBHEAD_MAP.items() if k in h3t),
                           h3t.replace(" ", "_"))
            lit.setdefault(current, {"papers": [], "see_all": {}})
            i += 1
            continue

        if node.name == "div" and "nl" in node.get("class", []):
            # "See all" 링크 처리
            sa_a = node.find("a", attrs={"data-ga-label":
                                          re.compile(r"See all", re.I)})
            if sa_a:
                lbl = sa_a.get("data-ga-label", "")
                m_n = re.search(r"\((\d+)\)", lbl)
                if m_n:
                    href = sa_a.get("href", "")
                    if href.startswith("/"):
                        href = "https://www.ncbi.nlm.nih.gov" + href
                    lit.setdefault(current, {"papers": [], "see_all": {}})
                    lit[current]["see_all"] = {
                        "count": int(m_n.group(1)), "url": href}
                i += 1
                continue

            # 논문 + detail
            detail = None
            j = i + 1
            while j < len(children):
                nxt = children[j]
                if not isinstance(nxt, Tag):
                    j += 1
                    continue
                if ("portlet_content" in nxt.get("class", [])
                        and "ln" in nxt.get("class", [])):
                    detail = nxt
                    i = j
                break

            paper = _parse_paper(node, detail)
            if paper:
                lit.setdefault(current, {"papers": [], "see_all": {}})
                lit[current]["papers"].append(paper)
        i += 1
    return lit


def parse_medgen_html(html: str, syndrome_name: str) -> dict[str, Any]:
    """
    MedGen HTML → disease_info dict.
    이 함수가 medgen_parser의 핵심. fetch_medgen()이 호출.
    """
    soup = BeautifulSoup(html, "lxml")

    # ── definition ───────────────────────────────────────────
    definition = ""
    portlet_100 = soup.find("div", id="ID_100")
    if portlet_100:
        c = portlet_100.find("div", class_="portlet_content")
        if c:
            definition = re.sub(r"\s+", " ", c.get_text(separator=" ", strip=True))
    if not definition:
        meta = soup.find("meta", attrs={"name": "description"})
        if meta and meta.get("content"):
            definition = meta["content"].strip()

    # ── medgenTable 순회 ──────────────────────────────────────
    synonyms:       list[str] = []
    db_ids:         dict[str, str] = {}
    gene_locations: list[dict] = []

    table = soup.find("table", class_="medgenTable")
    if table:
        for tr in table.find_all("tr"):
            tds = tr.find_all("td")
            if len(tds) < 2:
                continue
            label    = tds[0].get_text(strip=True).rstrip(":").lower()
            a_tags   = tds[1].find_all("a")
            first_a  = a_tags[0] if a_tags else None
            val      = (first_a.get_text(strip=True) if first_a
                        else tds[1].get_text(strip=True))

            if "synonym" in label:
                raw = tds[1].get_text(separator=";", strip=True)
                synonyms = [s.strip() for s in raw.split(";") if s.strip()]

            elif "gene" in label and "location" in label:
                for a in a_tags:
                    sym = a.get_text(strip=True)
                    if not sym:
                        continue
                    href  = a.get("href", "")
                    gid_m = re.search(r"/gene/(\d+)", href)
                    cytoband = ""
                    for sib in a.next_siblings:
                        cb = re.search(
                            r"\(([0-9XY]+[pq][\d.]+(?:[- ][pq]?[\d.]+)?)\)",
                            str(sib))
                        if cb:
                            cytoband = cb.group(1).strip()
                            break
                    gene_locations.append({
                        "gene_symbol":  sym,
                        "ncbi_gene_id": gid_m.group(1) if gid_m else None,
                        "cytoband":     cytoband or None,
                        "source":       "medgen_direct",
                    })

            elif "monarch" in label or "mondo" in label:
                db_ids["MONDO"] = val
            elif "omim" in label:
                db_ids["OMIM"] = val
            elif "orphanet" in label:
                db_ids["Orphanet"] = re.sub(r"^ORPHA", "", val)
            elif "snomed" in label:
                m2 = re.search(r"\((\d+)\)", val)
                db_ids["SNOMED-CT"] = m2.group(1) if m2 else val

    # ── GeneReview URL ────────────────────────────────────────
    genereview_url = None
    for a in soup.find_all("a", href=True):
        if "books.ncbi" in a["href"] or "genereviews" in a["href"].lower():
            title_txt = a.get_text(strip=True)
            if title_txt and len(title_txt) > 5:
                genereview_url = a["href"]
                break

    # ── literature ────────────────────────────────────────────
    literature: dict[str, Any] = {}

    # ID_105: Professional guidelines
    p105 = soup.find("div", id="ID_105")
    if p105:
        c = p105.find("div", class_="portlet_content")
        if c:
            sub = _walk_portlet(c, "professional_guidelines")
            literature.update(sub)

    # ID_103: Recent clinical studies (카테고리별)
    p103 = soup.find("div", id="ID_103")
    if p103:
        c = p103.find("div", class_="portlet_content")
        if c:
            sub = _walk_portlet(c, "general")
            literature.update(sub)

    # ID_104: Systematic reviews
    p104 = soup.find("div", id="ID_104")
    if p104:
        c = p104.find("div", class_="portlet_content")
        if c:
            sub = _walk_portlet(c, "systematic_reviews")
            literature.update(sub)

    # ── genomic_targets ───────────────────────────────────────
    genomic_targets = extract_genomic_targets(
        syndrome_name, synonyms, definition)

    return {
        "syndrome":        syndrome_name,
        "definition":      definition,
        "synonyms":        synonyms,
        "db_ids":          db_ids,
        "gene_locations":  gene_locations,
        "genomic_targets": genomic_targets,
        "genereview_url":  genereview_url,
        "literature":      literature,
    }


def fetch_medgen(syndrome_name: str, email: str) -> dict[str, Any]:
    """
    syndrome_name → MedGen 검색 + HTML 파싱 → disease_info 반환.
    """
    print(f"  [MedGen] '{syndrome_name}' 검색")

    r = http_get(f"{EUTILS}/esearch.fcgi",
                 params={"db": "medgen", "term": f'"{syndrome_name}"',
                         "retmax": 1, "sort": "relevance",
                         "retmode": "json", "email": email})
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        print(f"  [MedGen] '{syndrome_name}' 결과 없음")
        return {"syndrome": syndrome_name, "definition": None,
                "synonyms": [], "db_ids": {}, "gene_locations": [],
                "genomic_targets": [], "genereview_url": None,
                "literature": {}}

    uid = ids[0]
    url = f"{MEDGEN_BASE}/{uid}"
    print(f"  [MedGen] UID={uid}")

    time.sleep(1.2)
    resp = get_session().get(url, timeout=20,
                             headers={"User-Agent": "Mozilla/5.0 Chrome/124"})
    resp.raise_for_status()

    result = parse_medgen_html(resp.text, syndrome_name)
    result["medgen_uid"] = uid
    result["medgen_url"] = url

    print(f"  [MedGen] gene_loc={[g['gene_symbol'] for g in result['gene_locations']]} "
          f"targets={[t['cytoband'] for t in result['genomic_targets']]} "
          f"lit_cats={list(result['literature'].keys())}")
    return result


# ── 단독 실행 ─────────────────────────────────────────────────
def _main():
    ap = argparse.ArgumentParser(
        description="MedGen 파싱 → disease_info JSON",
        epilog="python -m nipt_pipeline.parsers.medgen_parser "
               "--syndrome 'Down syndrome' --output disease_info.json")
    ap.add_argument("--syndrome", required=True)
    ap.add_argument("--email",    default="your@email.com")
    ap.add_argument("--output",   default="-",
                    help="출력 JSON 파일 경로 (- = stdout)")
    args = ap.parse_args()

    from nipt_pipeline.utils import init_session
    init_session(args.email)
    result = fetch_medgen(args.syndrome, args.email)

    out = json.dumps(result, ensure_ascii=False, indent=2)
    if args.output == "-":
        print(out)
    else:
        import pathlib
        pathlib.Path(args.output).write_text(out, encoding="utf-8")
        print(f"저장: {args.output}", file=sys.stderr)


import json  # noqa: E402 (필요한 위치)

if __name__ == "__main__":
    _main()