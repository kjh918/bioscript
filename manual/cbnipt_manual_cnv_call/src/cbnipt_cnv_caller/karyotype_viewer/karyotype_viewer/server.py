"""
FastAPI host server factory.
"""

from __future__ import annotations
from typing import Optional

import pandas as pd

from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from starlette.middleware.wsgi import WSGIMiddleware

from .core.models import SampleData
from .core.nipt_markers import NiptSyndrome
from .apps.overview import create_overview_app
from .apps.detail import create_detail_app


_HOST_HTML_TEMPLATE = r"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  * {{ box-sizing: border-box; }}
  html, body {{ margin: 0; padding: 0; background: #edf2f7; font-family: Arial, sans-serif; }}
  iframe {{ display: block; width: 100%; border: 0; margin: 0; padding: 0; background: #edf2f7; }}
  #detail-frame {{ height: 1100px; }}
</style>
</head>
<body>
<div>
  <iframe id="overview-frame" src="/overview/" title="Karyotype overview" scrolling="no"></iframe>
  <iframe id="detail-frame"   src="/detail/"   title="Chromosome detail"  scrolling="no"></iframe>
</div>
<script>
(function () {{
  const detailFrame   = document.getElementById('detail-frame');
  const overviewFrame = document.getElementById('overview-frame');
  let lastChrom = '{initial_chrom}';

  function resizeOverview() {{
    try {{
      const h = overviewFrame.contentDocument.body.scrollHeight;
      if (h > 100) overviewFrame.style.height = h + 'px';
    }} catch(e) {{}}
  }}
  overviewFrame.addEventListener('load', function () {{
    resizeOverview();
    setInterval(resizeOverview, 800);
  }});

  function forward(chrom) {{
    lastChrom = String(chrom || lastChrom);
    if (detailFrame && detailFrame.contentWindow) {{
      detailFrame.contentWindow.postMessage(
        {{type: 'karyotype-chromosome-selected', chrom: lastChrom}},
        window.location.origin
      );
    }}
  }}

  window.addEventListener('message', function (event) {{
    if (event.origin !== window.location.origin) return;
    const data = event.data || {{}};
    if (data.type !== 'karyotype-chromosome-selected') return;
    forward(data.chrom);
  }});

  detailFrame.addEventListener('load', function () {{
    setTimeout(function () {{ forward(lastChrom); }}, 150);
  }});
}})();
</script>
</body>
</html>
"""


def create_app(
    sample: SampleData,
    initial_chrom: str = "21",
    title: str = "Karyotype Viewer",
    cnv_data: Optional[dict[str, pd.DataFrame]] = None,
    syndromes: Optional[dict[str, NiptSyndrome]] = None,
) -> FastAPI:
    if initial_chrom not in sample.display_chroms:
        initial_chrom = sample.display_chroms[0]

    overview_dash = create_overview_app(
        sample, initial_chrom=initial_chrom, requests_prefix="/overview/"
    )
    detail_dash = create_detail_app(
        sample,
        initial_chrom=initial_chrom,
        requests_prefix="/detail/",
        cnv_data=cnv_data,
        syndromes=syndromes,
    )

    host_html = _HOST_HTML_TEMPLATE.format(
        title=title, initial_chrom=initial_chrom,
    )

    fastapi_app = FastAPI(title=title)

    @fastapi_app.get("/", response_class=HTMLResponse)
    async def _root():
        return HTMLResponse(host_html)

    fastapi_app.mount("/overview", WSGIMiddleware(overview_dash.server))
    fastapi_app.mount("/detail",   WSGIMiddleware(detail_dash.server))
    return fastapi_app
