"""
S04 — 검사 결과 요약
성별, 종합 판정, ISCN — 선명하게 강조
JS: renderSampleInfo()
"""

SECTION_SAMPLE = """
<!-- ═══════════════════════════════════════════════════════ S04 RESULT SUMMARY -->
<div class="card" id="result-summary-card" style="margin-bottom:0.75rem;">
  <div class="ch" id="result-summary-ch">
    <span class="ct">검사 결과 요약</span>
    <span style="font-size:10px;color:var(--text-muted);">Analysis Result Summary</span>
  </div>

  <div style="padding:18px 22px;display:grid;
              grid-template-columns:130px 1px 1fr 1px 240px;
              gap:24px;align-items:center;">

    <!-- 태아 성별 -->
    <div style="text-align:center;">
      <div style="font-size:10px;color:var(--text-muted);margin-bottom:6px;
                  letter-spacing:.08em;text-transform:uppercase;">태아 성별 / Fetal Sex</div>
      <div id="s-sex-icon" style="font-size:52px;line-height:1;">—</div>
      <div id="s-fetal-sex" style="font-size:13px;font-weight:700;margin-top:8px;">—</div>
    </div>

    <div style="background:var(--border);width:1px;height:80px;align-self:center;"></div>

    <!-- 종합 판정 -->
    <div style="text-align:center;">
      <div style="font-size:10px;color:var(--text-muted);margin-bottom:10px;
                  letter-spacing:.08em;text-transform:uppercase;">종합 판정 / Overall Call</div>
      <div id="s-call-badge"
           style="display:inline-flex;align-items:center;gap:8px;
                  padding:8px 28px;border-radius:99px;margin-bottom:10px;
                  border:2px solid var(--border);">
        <span id="s-call-icon" style="font-size:18px;"></span>
        <span id="s-call-text" style="font-size:22px;font-weight:800;">—</span>
      </div>
      <div id="s-events-count" style="font-size:12px;color:var(--text-sub);"></div>
    </div>

    <div style="background:var(--border);width:1px;height:80px;align-self:center;"></div>

    <!-- ISCN -->
    <div>
      <div style="font-size:10px;color:var(--text-muted);margin-bottom:6px;
                  letter-spacing:.08em;text-transform:uppercase;">핵형 / ISCN</div>
      <div id="s-iscn"
           style="font-family:var(--mono);font-size:12px;font-weight:700;
                  color:var(--navy);background:var(--navy-l);
                  padding:6px 10px;border-radius:6px;margin-bottom:6px;
                  word-break:break-all;display:inline-block;max-width:100%;">—</div>
    </div>

  </div>
</div>
"""