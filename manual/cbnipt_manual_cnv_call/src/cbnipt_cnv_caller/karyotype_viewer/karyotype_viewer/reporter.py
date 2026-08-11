"""
reporter.py
===========
HTML 템플릿 조합 + asset 복사 + manifest 생성.
실제 CSS/HTML/JS 내용은 templates/ 패키지에 분리.

Public API
----------
build_html(mode, inline_manifest, extra_inline_js, inline_assets, asset_base)
save_report_dir(output_dir, sample, syndromes, cnv_data, ...)
save_report(path, ...)   ← legacy shim
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Optional

import pandas as pd

from .core.models      import SampleData
from .core.annot       import build_iscn, build_raw_annots
from .core.nipt_markers import NiptSyndrome, GROUP_COLORS
from .core.reference   import CHROM_SIZES
from .core.data        import demo_dataframe

from .templates        import CSS, html_head, html_body, JS_RENDER, JS_IDEOGRAM, JS_CNV

_GROUP_ORDER = [
    "Autosome Abnormality",
    "Sex Chromosome Abnormality",
    "Micro Deletion",
]

# CHROM_SIZES JS 문자열 (ideogram brush 범위용)
_CHROM_SIZES_JS = (
    "var CHROM_SIZES = " +
    json.dumps({k: v for k, v in CHROM_SIZES.items()}, ensure_ascii=False) + ";"
)


# ---------------------------------------------------------------------------
# Manifest builder
# ---------------------------------------------------------------------------
def _overall_call(syndromes: dict) -> str:
    calls = [s.call for s in syndromes.values()]
    if "HIGH_RISK"   in calls: return "HIGH_RISK"
    if "SUSPECTED" in calls: return "SUSPECTED"
    return "LOW_RISK"


def _build_manifest(
    sample: SampleData,
    syndromes: dict[str, NiptSyndrome],
    affected_chroms: list[str],
    report_id: str,
    report_date: str,
    mode: str = "fetch",
) -> dict:
    iscn       = build_iscn(sample)
    raw_annots = build_raw_annots(sample)

    syn_list = []
    for grp in _GROUP_ORDER:
        for s in sorted(syndromes.values(), key=lambda x: x.syndrome):
            if s.group != grp:
                continue
            features = [
                {
                    "name":  f.feature_name,
                    "type":  f.feature_type,
                    "chrom": f.chromosome.replace("chr", ""),
                    "start": f.start,
                    "end":   f.end,
                }
                for f in s.features
            ]
            syn_list.append({
                "nipt_id":       s.nipt_id,
                "syndrome":      s.syndrome,
                "group":         s.group,
                "call":          s.call,
                "cn_value":      round(s.cn_value, 4) if s.cn_value is not None else None,
                "primary_chrom": s.primary_chrom,
                "cnv_file":      (
                    f"/api/cnv/{s.primary_chrom}" if mode == "inline"
                    else f"data/cnv/chr{s.primary_chrom}.tsv"
                ) if s.primary_chrom in affected_chroms else None,
                "features": features,
            })

    return {
        "report_id":       report_id,
        "report_date":     report_date,
        "test_date":       report_date,
        "sample": {
            "id":            sample.id,
            "sex":           sample.sex,
            "iscn":          iscn,
            "display_chroms": sample.display_chroms,
            "pipeline":      "cbNIPT v2",
        },
        "maternal": getattr(sample, "maternal", {}),
        "qc":       getattr(sample, "qc", {}),
        "institution": getattr(sample, "institution", {}),
        "signatures":  getattr(sample, "signatures", {}),
        "overall_call":    _overall_call(syndromes),
        "affected_chroms": affected_chroms,
        "syndromes":       syn_list,
        "events": [
            {
                "chr":   ev.chr,  "type":  ev.type,
                "cn":    ev.cn,   "iscn":  ev.iscn,
                "start": ev.start,"stop":  ev.stop,
                "color": ev.color,
            }
            for ev in sample.events
        ],
        "raw_annots":  raw_annots,
        "chrom_sizes": {k: v for k, v in CHROM_SIZES.items()},
    }


# ---------------------------------------------------------------------------
# build_html — templates 조합
# ---------------------------------------------------------------------------
def build_html(
    mode:            str  = "static",
    inline_manifest: dict = None,
    extra_inline_js: str  = "",
    inline_assets:   dict = None,
    asset_base:      str  = "assets",
) -> str:
    """
    mode='static' → <script src> 로 ideogram 로드 (file:// 지원)
    mode='inline' → manifest + CNV 전부 JS inline, ideogram은 /assets/ 엔드포인트
    """
    ia = inline_assets or {}

    # ── plotly tag ──────────────────────────────────────────────────────
    if "plotly" in ia and ia["plotly"]:
        plotly_tag = f"<script>{ia['plotly']}</script>"
    else:
        plotly_tag = f'<script src="{asset_base}/plotly.min.js"></script>'

    # ── ideogram tag — 항상 src (f-string 중괄호 충돌 방지) ─────────────
    _IDEO_PH  = "___IDEOGRAM_SCRIPT___"
    if mode == "static":
        ideogram_real = '<script src="data/assets/ideogram.min.js"></script>'
    else:
        ideogram_real = f'<script src="{asset_base}/ideogram.min.js"></script>'

    # ── manifest + cnv vars ──────────────────────────────────────────────
    if mode == "inline":
        manifest_js    = (
            f"var MANIFEST = {json.dumps(inline_manifest, ensure_ascii=False)};\n"
            f"{_CHROM_SIZES_JS}\n"
            f"{extra_inline_js}"
        )
        extra_scripts  = ""
        init_call      = "render();"
    else:
        manifest_js   = _CHROM_SIZES_JS
        extra_scripts = '<script src="data/manifest.js"></script>'
        init_call     = "render();"

    # ── 조합 ────────────────────────────────────────────────────────────
    head = html_head(
        plotly_tag   = plotly_tag,
        ideogram_tag = _IDEO_PH,        # placeholder
        extra_scripts= extra_scripts,
        css          = CSS,
    )
    body = html_body(
        manifest_js = manifest_js,
        init_call   = init_call,
        js_blocks   = [JS_RENDER, JS_IDEOGRAM, JS_CNV],
    )
    return (head + body).replace(_IDEO_PH, ideogram_real)


# ---------------------------------------------------------------------------
# save_report_dir
# ---------------------------------------------------------------------------
def save_report_dir(
    output_dir:   "str | Path",
    sample:       SampleData,
    syndromes:    Optional[dict[str, NiptSyndrome]] = None,
    cnv_data:     Optional[dict[str, pd.DataFrame]] = None,
    report_id:    str  = "GCX-REPORT",
    report_date:  str  = "",
    affected_only: bool = True,
) -> Path:
    import datetime
    syndromes = syndromes or {}
    cnv_data  = cnv_data  or {}
    if not report_date:
        report_date = datetime.date.today().isoformat()

    out = Path(output_dir)
    (out / "assets").mkdir(parents=True, exist_ok=True)
    (out / "data" / "assets").mkdir(parents=True, exist_ok=True)
    (out / "data" / "cnv").mkdir(parents=True, exist_ok=True)

    # 1. plotly.min.js → assets/
    import plotly as _plotly
    plotly_src = Path(_plotly.__file__).parent / "package_data" / "plotly.min.js"
    if plotly_src.exists():
        shutil.copy2(plotly_src, out / "assets" / "plotly.min.js")

    # 2. ideogram.min.js → data/assets/ (<script src> 상대경로)
    ideo_src = Path(__file__).parent / "assets" / "ideogram.min.js"
    if ideo_src.exists():
        shutil.copy2(ideo_src, out / "data" / "assets" / "ideogram.min.js")
        print(f"[reporter] ideogram.min.js → data/assets/")
    else:
        print(f"[reporter] WARNING: ideogram.min.js not found at {ideo_src}")

    # 3. affected chroms
    affected = sorted(
        {s.primary_chrom for s in syndromes.values() if s.primary_chrom in CHROM_SIZES},
        key=lambda c: (int(c) if c.isdigit() else 100 + ord(c[0])),
    )
    chroms_to_export = affected if affected_only else list(CHROM_SIZES.keys())

    # 4. CNV → .tsv + .js (var 전역변수)
    for chrom in chroms_to_export:
        df = cnv_data[chrom] if chrom in cnv_data else demo_dataframe(sample, chrom)
        df.to_csv(out / "data" / "cnv" / f"chr{chrom}.tsv", sep="\t", index=False)
        rows = df.to_dict(orient="records")
        js   = f"var CNV_CHR{chrom.upper()} = {json.dumps(rows, ensure_ascii=False)};"
        (out / "data" / "cnv" / f"chr{chrom}.js").write_text(js, encoding="utf-8")
        print(f"[reporter] cnv/chr{chrom}.js  ({len(df):,} rows)")

    # 5. manifest.js (var)
    manifest = _build_manifest(sample, syndromes, affected, report_id, report_date,
                                mode="static")
    manifest_js = f"var MANIFEST = {json.dumps(manifest, ensure_ascii=False, indent=2)};"
    (out / "data" / "manifest.js").write_text(manifest_js, encoding="utf-8")
    (out / "data" / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # 6. plotly inline 읽기
    plotly_js = plotly_src.read_text(encoding="utf-8", errors="replace") if plotly_src.exists() else ""

    # 7. index.html
    html = build_html(
        mode          = "static",
        asset_base    = "assets",
        inline_assets = {"plotly": plotly_js},
    )
    (out / "index.html").write_text(html, encoding="utf-8")

    print(f"[reporter] → {out}/index.html")
    return out / "index.html"


# ---------------------------------------------------------------------------
# Legacy shim
# ---------------------------------------------------------------------------
def save_report(path, sample, syndromes=None, cnv_data=None,
                report_id="GCX-REPORT", report_date="") -> Path:
    out_dir = Path(path).parent / Path(path).stem
    return save_report_dir(out_dir, sample, syndromes, cnv_data, report_id, report_date)
