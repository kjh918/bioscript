"""
S06 — Karyotype Overview
Ideogram.js rows=2 전체 핵형 + proband annotation
JS: renderIdeogram(), onClickChromosome()
"""

SECTION_KARYOTYPE = """
<!-- ═══════════════════════════════════════════════════════ S06 KARYOTYPE -->
<div class="card">
  <div class="ch">
    <span class="ct">Karyotype Overview</span>
    <div style="display:flex;align-items:center;gap:10px;">
      <!-- Overview annotation toggle -->
      <span id="ideo-selected"
            style="font-family:var(--mono);color:var(--navy);font-weight:600;font-size:11px;"></span>
      <span class="no-print" style="font-size:10px;color:var(--text-muted);">
        염색체를 클릭하면 상세 보기로 이동합니다
      </span>
    </div>
  </div>

  <!-- Ideogram 렌더 타깃 — overflow-x:auto로 좁은 화면 대응 -->
  <div id="ideogram-wrap"
       style="width:100%;overflow-x:auto;overflow-y:visible;
              padding:10px 8px 20px;background:white;min-height:300px;"></div>

  <!-- Legend row -->
  <div style="display:flex;gap:16px;flex-wrap:wrap;padding:6px 14px 10px;
              border-top:1px solid var(--border);background:var(--surface2);font-size:11px;">
    <span style="color:var(--text-muted);font-weight:600;">범례:</span>
    <span><span style="display:inline-block;width:12px;height:12px;background:#E53E3E;
                        border-radius:2px;margin-right:4px;vertical-align:middle;"></span>ABNORMAL</span>
    <span><span style="display:inline-block;width:12px;height:12px;background:#DD6B20;
                        border-radius:2px;margin-right:4px;vertical-align:middle;"></span>SUSPICIOUS</span>
    <span><span style="display:inline-block;width:12px;height:12px;background:#38A169;
                        border-radius:2px;margin-right:4px;vertical-align:middle;"></span>NORMAL</span>
    <span><span style="display:inline-block;width:10px;height:10px;background:#6B46C1;
                        border-radius:50%;margin-right:4px;vertical-align:middle;"></span>Gene</span>
  </div>
</div>
"""
