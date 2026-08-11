"""
FastAPI 단일 페이지 서버.
plotly + manifest + CNV → inline embed
ideogram → /assets/ideogram.min.js 엔드포인트 (패키지 내부 파일)
"""

from __future__ import annotations
from typing import Optional
import json
from pathlib import Path

import pandas as pd
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, PlainTextResponse

from .core.models       import SampleData
from .core.nipt_markers import NiptSyndrome
from .core.reference    import CHROM_SIZES
from .reporter          import build_html, _build_manifest

_ASSETS_DIR = Path(__file__).parent / "assets"


def create_app(
    sample:        SampleData,
    initial_chrom: str = "21",
    title:         str = "Karyotype Viewer",
    cnv_data:      Optional[dict[str, pd.DataFrame]] = None,
    syndromes:     Optional[dict[str, NiptSyndrome]]  = None,
) -> FastAPI:
    import datetime
    from .core.data import demo_dataframe

    cnv_data  = cnv_data  or {}
    syndromes = syndromes or {}

    if initial_chrom not in sample.display_chroms:
        initial_chrom = sample.display_chroms[0]

    report_date = datetime.date.today().isoformat()
    affected = sorted(
        {s.primary_chrom for s in syndromes.values() if s.primary_chrom in CHROM_SIZES},
        key=lambda c: (int(c) if c.isdigit() else 100 + ord(c[0])),
    )

    manifest = _build_manifest(sample, syndromes, affected, title, report_date,
                                mode="inline")

    # CNV → var 전역변수 (window에 바인딩)
    cnv_vars = "\n".join(
        f"var CNV_CHR{chrom.upper()} = "
        f"{json.dumps((cnv_data[chrom] if chrom in cnv_data else demo_dataframe(sample, chrom)).to_dict(orient='records'), ensure_ascii=False)};"
        for chrom in affected
    )

    # plotly inline
    import plotly as _plotly
    plotly_js = (
        Path(_plotly.__file__).parent / "package_data" / "plotly.min.js"
    ).read_text(encoding="utf-8", errors="replace")

    html = build_html(
        mode            = "inline",
        inline_manifest = manifest,
        extra_inline_js = cnv_vars,
        inline_assets   = {"plotly": plotly_js},
        asset_base      = "/assets",
    )

    fastapi_app = FastAPI(title=title)

    @fastapi_app.get("/", response_class=HTMLResponse)
    async def _root():
        return HTMLResponse(html)

    @fastapi_app.get("/assets/{filename}")
    async def _assets(filename: str):
        p = _ASSETS_DIR / filename
        if not p.exists():
            raise HTTPException(404, f"{filename} not found")
        return PlainTextResponse(
            p.read_text(encoding="utf-8", errors="replace"),
            media_type="application/javascript",
        )

    return fastapi_app
