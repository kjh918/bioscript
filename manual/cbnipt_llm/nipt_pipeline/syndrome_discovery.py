"""
syndrome_discovery.py
──────────────────────
질환명 하나 → 최종 JSON 하나

각 parser/annotator 단계가 완료될 때마다 중간 파일로 저장.
이미 중간 파일이 있으면 해당 step을 건너뜀 (resume 지원).

출력 디렉토리 구조:
  {output_dir}/
  └─ {safe_syndrome_name}/
       ├─ 01_disease_info.json      ← Step 1: MedGen
       ├─ 02_articles.json          ← Step 2: PubMed fetch/discover
       ├─ 03_markers.json           ← Step 3+4: DB 로드 + 마커 구성
       └─ discovery.json            ← 최종 통합 결과

단독 실행:
  python -m nipt_pipeline.syndrome_discovery \\
      --syndrome "Down syndrome" \\
      --email you@example.com \\
      --output-dir ./output

  # resume: 이미 완료된 step은 건너뜀
  python -m nipt_pipeline.syndrome_discovery \\
      --syndrome "Down syndrome" \\
      --email you@example.com \\
      --output-dir ./output \\
      --resume

  # 특정 step부터 재실행 (중간 파일 무시)
  python -m nipt_pipeline.syndrome_discovery \\
      --syndrome "Down syndrome" \\
      --email you@example.com \\
      --output-dir ./output \\
      --from-step 2
"""

import re
import sys
import json
import argparse
import pathlib
import datetime
from typing import Any


# ── 중간 파일 이름 ────────────────────────────────────────────
STEP_FILES = {
    1: "01_disease_info.json",
    2: "02_articles.json",
    3: "03_markers.json",
}
FINAL_FILE = "discovery.json"


def _safe_name(syndrome_name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", syndrome_name.lower()).strip("_")


def _save(path: pathlib.Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2),
                    encoding="utf-8")
    size_kb = path.stat().st_size / 1024
    print(f"  → 저장: {path.name} ({size_kb:.1f} KB)")


def _load(path: pathlib.Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


# ── Step 함수 ─────────────────────────────────────────────────

def step1_medgen(syndrome_name: str, email: str) -> dict:
    """MedGen 파싱 → disease_info"""
    from nipt_pipeline.parsers.medgen_parser import fetch_medgen
    return fetch_medgen(syndrome_name, email)


def step2_pubmed(
    syndrome_name: str,
    email: str,
    disease_info: dict,
    pubmed_limits: dict | None = None,
    verify_pmc: bool = False,
) -> dict:
    """PubMed discover + fetch → {discovery, pmids, articles}"""
    from nipt_pipeline.parsers.pubmed_parser import discover_and_fetch
    return discover_and_fetch(
        disease_name=syndrome_name,
        synonyms=disease_info.get("synonyms"),
        email=email,
        medgen_literature=disease_info.get("literature"),
        limits=pubmed_limits,
        verify_pmc=verify_pmc,
    )


def step3_markers(
    disease_info: dict,
    email: str,
    clingen_gene_path: str = "",
    gencc_path: str = "",
    hpo_path: str = "",
    fetch_ncbi: bool = True,
) -> dict:
    """외부 DB 로드 + 마커 구성 → {cg_db, gc_db, hp_db, markers}"""
    from nipt_pipeline.utils import resolve_path
    from nipt_pipeline.annotators.marker_builder import build_markers

    cg_db, gc_db, hp_db = {}, {}, {}

    # ClinGen
    cg_path = clingen_gene_path or resolve_path("clingen_gene")[0]
    if cg_path:
        from nipt_pipeline.parsers.clingen_parser import parse_clingen_gene
        cg_db = parse_clingen_gene(cg_path)

    # GenCC
    gc_path = gencc_path or resolve_path("gencc")[0]
    if gc_path:
        from nipt_pipeline.parsers.gencc_parser import load_gencc
        gc_db = load_gencc(gc_path)

    # HPO
    hp_path = hpo_path or resolve_path("hpo")[0]
    if hp_path:
        from nipt_pipeline.parsers.hpo_parser import load_hpo
        hp_db = load_hpo(hp_path)

    markers = build_markers(
        disease_info=disease_info,
        clingen_gene_db=cg_db,
        gencc_db=gc_db,
        hpo_db=hp_db,
    )

    return {
        "db_sources": {
            "clingen": cg_path or None,
            "gencc":   gc_path or None,
            "hpo":     hp_path or None,
        },
        "markers": markers,
    }


# ── Orchestrator ──────────────────────────────────────────────

def discover(
    syndrome_name:      str,
    email:              str,
    output_dir:         str  = ".",
    clingen_gene_path:  str  = "",
    gencc_path:         str  = "",
    hpo_path:           str  = "",
    uniprot_dat:        str  = "",
    pubmed_limits:      dict | None = None,
    verify_pmc:         bool = False,
    fetch_ncbi:         bool = True,
    resume:             bool = True,
    from_step:          int  = 1,
) -> dict[str, Any]:
    """
    syndrome_name → 최종 discovery dict.

    Parameters
    ----------
    resume    : True면 중간 파일이 있을 때 해당 step 건너뜀
    from_step : 이 step부터 재실행 (이전 중간 파일은 재사용)
    """
    from nipt_pipeline.utils import init_session
    init_session(email)

    # 작업 디렉토리: output_dir/{safe_syndrome_name}/
    work_dir = pathlib.Path(output_dir) / _safe_name(syndrome_name)
    work_dir.mkdir(parents=True, exist_ok=True)

    def _should_run(step: int) -> bool:
        if step < from_step:
            return False                    # from_step 이전은 무조건 스킵
        if not resume:
            return True                     # resume=False면 항상 실행
        path = work_dir / STEP_FILES[step]
        if path.exists():
            print(f"  [Step {step}] 중간 파일 존재 → 건너뜀: {path.name}")
            return False
        return True

    def _get(step: int) -> Any:
        """중간 파일 로드"""
        return _load(work_dir / STEP_FILES[step])

    # ── Step 1: MedGen ────────────────────────────────────────
    print(f"\n[1/3] MedGen 파싱")
    if _should_run(1):
        try:
            disease_info = step1_medgen(syndrome_name, email)
            _save(work_dir / STEP_FILES[1], disease_info)
        except Exception as e:
            print(f"  [MedGen] 오류: {e}")
            disease_info = {
                "syndrome": syndrome_name, "synonyms": [],
                "db_ids": {}, "gene_locations": [],
                "genomic_targets": [], "literature": {},
            }
            _save(work_dir / STEP_FILES[1], disease_info)
    else:
        disease_info = _get(1)

    # ── Step 2: PubMed ────────────────────────────────────────
    print(f"\n[2/3] PubMed 문헌 수집")
    if _should_run(2):
        try:
            pub_result = step2_pubmed(
                syndrome_name, email, disease_info,
                pubmed_limits=pubmed_limits,
                verify_pmc=verify_pmc,
            )
            _save(work_dir / STEP_FILES[2], pub_result)
        except Exception as e:
            print(f"  [PubMed] 오류: {e}")
            pub_result = {"discovery": {}, "pmids": [], "articles": {}}
            _save(work_dir / STEP_FILES[2], pub_result)
    else:
        pub_result = _get(2)

    # ── Step 3: 외부 DB + 마커 구성 ──────────────────────────
    print(f"\n[3/3] 마커 후보 구성")
    if _should_run(3):
        try:
            marker_result = step3_markers(
                disease_info=disease_info,
                email=email,
                clingen_gene_path=clingen_gene_path,
                gencc_path=gencc_path,
                hpo_path=hpo_path,
                fetch_ncbi=fetch_ncbi,
            )
            _save(work_dir / STEP_FILES[3], marker_result)
        except Exception as e:
            print(f"  [marker] 오류: {e}")
            marker_result = {"db_sources": {}, "markers": []}
            _save(work_dir / STEP_FILES[3], marker_result)
    else:
        marker_result = _get(3)

    # ── 최종 통합 ─────────────────────────────────────────────
    markers    = marker_result.get("markers", [])
    articles   = pub_result.get("articles", {})

    result: dict[str, Any] = {
        "disease_info":        disease_info,
        "db_ids":              disease_info.get("db_ids", {}),
        "markers":             markers,
        "evidence_literature": list(articles.values()),
        "pubmed_discovery":    pub_result.get("discovery", {}),
        "meta": {
            "syndrome":     syndrome_name,
            "collected_at": datetime.datetime.utcnow().isoformat() + "Z",
            "work_dir":     str(work_dir),
            "steps": {
                "01_disease_info": str(work_dir / STEP_FILES[1]),
                "02_articles":     str(work_dir / STEP_FILES[2]),
                "03_markers":      str(work_dir / STEP_FILES[3]),
            },
            "annotation_sources": {
                "medgen":    "NCBI MedGen HTML (공개)",
                "pubmed":    "NCBI PubMed ESearch+EFetch (공개)",
                "ncbi_gene": "NCBI Gene ESearch+ESummary" if fetch_ncbi else "skip",
                **marker_result.get("db_sources", {}),
            },
        },
    }

    final_path = work_dir / FINAL_FILE
    _save(final_path, result)

    # ── 요약 ─────────────────────────────────────────────────
    n_gene   = sum(1 for m in markers if m["feature_type"] == "CoreGene")
    n_region = sum(1 for m in markers if m["feature_type"] == "CoreRegion")
    n_chrom  = sum(1 for m in markers if m["feature_type"] == "TargetChromosome")
    print(f"\n▶ 완료: {syndrome_name}")
    print(f"  markers     : CoreGene={n_gene} CoreRegion={n_region} TargetChromosome={n_chrom}")
    print(f"  literature  : {len(articles)}편")
    print(f"  db_ids      : {result['db_ids']}")
    print(f"  work_dir    : {work_dir}/")
    print(f"    ├─ {STEP_FILES[1]}")
    print(f"    ├─ {STEP_FILES[2]}")
    print(f"    ├─ {STEP_FILES[3]}")
    print(f"    └─ {FINAL_FILE}")

    return result


# ── 단독 실행 ─────────────────────────────────────────────────
def _main():
    ap = argparse.ArgumentParser(
        description="질환명 → discovery JSON (중간 파일 저장 지원)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
출력 구조:
  output-dir/
  └─ {syndrome}/
       ├─ 01_disease_info.json   ← MedGen
       ├─ 02_articles.json       ← PubMed
       ├─ 03_markers.json        ← 마커 후보
       └─ discovery.json         ← 최종 통합

예시:
  # 처음 실행
  python -m nipt_pipeline.syndrome_discovery \\
      --syndrome "DiGeorge syndrome" --email you@example.com \\
      --output-dir ./output

  # 재실행 (완료된 step 건너뜀)
  python -m nipt_pipeline.syndrome_discovery \\
      --syndrome "DiGeorge syndrome" --email you@example.com \\
      --output-dir ./output --resume

  # Step 2부터 재실행 (01은 재사용, 02~부터 덮어씀)
  python -m nipt_pipeline.syndrome_discovery \\
      --syndrome "DiGeorge syndrome" --email you@example.com \\
      --output-dir ./output --from-step 2
        """)
    ap.add_argument("--syndrome",     required=True)
    ap.add_argument("--email",        default="your@email.com")
    ap.add_argument("--clingen-gene", default="")
    ap.add_argument("--gencc",        default="")
    ap.add_argument("--hpo",          default="")
    ap.add_argument("--uniprot-dat",  default="")
    ap.add_argument("--no-ncbi",      action="store_true")
    ap.add_argument("--verify-pmc",   action="store_true")
    ap.add_argument("--output-dir",   default=".")
    ap.add_argument("--resume",       action="store_true",
                    help="중간 파일 있으면 해당 step 건너뜀")
    ap.add_argument("--from-step",    type=int, default=1, choices=[1, 2, 3],
                    help="이 step부터 재실행 (기본: 1)")
    args = ap.parse_args()

    discover(
        syndrome_name=args.syndrome,
        email=args.email,
        output_dir=args.output_dir,
        clingen_gene_path=args.clingen_gene,
        gencc_path=args.gencc,
        hpo_path=args.hpo,
        uniprot_dat=args.uniprot_dat,
        verify_pmc=args.verify_pmc,
        fetch_ncbi=not args.no_ncbi,
        resume=args.resume,
        from_step=args.from_step,
    )


if __name__ == "__main__":
    _main()