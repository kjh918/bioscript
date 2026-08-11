"""
S02 — 종합 판정 배너
ABNORMAL / SUSPICIOUS / NORMAL 색상 배너 + 요약 chips
JS: renderBanner(), renderChips()
"""

SECTION_BANNER = """
<!-- ═══════════════════════════════════════════════════════ S02 BANNER -->
<div id="banner" class="banner NORMAL">
  <div class="bicon">…</div>
  <div>
    <div class="rlbl">종합 판정</div>
    <div class="bres" id="banner-result">로드 중…</div>
    <div class="bdesc" id="banner-desc"></div>
  </div>
</div>

<!-- Summary chips -->
<div class="chips" id="chips"></div>
"""
