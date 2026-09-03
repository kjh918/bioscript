"""
HTML 골격 — sections/ 를 순서대로 조합.

html_head(plotly_tag, ideogram_tag, extra_scripts, css) → str
html_body(manifest_js, init_call, js_blocks, sections)  → str
"""

from .sections import ALL_SECTIONS


def html_head(
    plotly_tag:    str,
    ideogram_tag:  str,
    extra_scripts: str,
    css:           str,
) -> str:
    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Karyotype Report</title>
{plotly_tag}
{ideogram_tag}
{extra_scripts}
<style>
{css}
</style>
</head>"""


def html_body(
    manifest_js: str,
    init_call:   str,
    js_blocks:   list[str],
    sections:    list[str] | None = None,
) -> str:
    """
    sections=None → ALL_SECTIONS 순서대로 전체 렌더
    sections=[...] → 지정한 section 문자열만 렌더 (커스텀 순서/선택 가능)
    """
    body_sections = "\n".join(sections if sections is not None else ALL_SECTIONS)
    js_all        = "\n\n".join(js_blocks)

    return f"""
<body>
<div class="page">

{body_sections}

</div><!-- /page -->

<script>
// ── Data ─────────────────────────────────────────────────────────────────
{manifest_js}

// ── Application ──────────────────────────────────────────────────────────
{js_all}

// ── Init ─────────────────────────────────────────────────────────────────
{init_call}
</script>
</body>
</html>"""
