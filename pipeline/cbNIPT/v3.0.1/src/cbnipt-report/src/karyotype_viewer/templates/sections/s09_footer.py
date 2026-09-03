"""
S09 — Footer / 서명란
기관명, 담당자, 서명란, 면책조항
JS: renderFooter()
"""

SECTION_FOOTER = """
<!-- ═══════════════════════════════════════════════════════ S09 FOOTER -->

<!-- 면책조항 -->
<div style="background:var(--surface2);border:1px solid var(--border);
            border-radius:var(--r);padding:12px 16px;margin-bottom:1rem;font-size:11px;
            color:var(--text-muted);line-height:1.7;">
  <strong style="color:var(--text-sub);">※ 본 검사 결과의 해석에 관한 주의사항</strong><br>
  본 검사는 임산부 혈액에서 분리한 순환태아세포를 이용한 선별검사로, 확진 검사가 아닙니다.
  양성 결과 시 반드시 양수 검사 또는 융모막 검사 등 확진 검사를 시행하여야 합니다.
  본 검사 결과는 임상 증상 및 기타 검사 결과와 종합하여 해석하여야 하며,
  검사 결과 해석에 대한 최종 책임은 의뢰 의사에게 있습니다.
</div>

<!-- 서명란 -->
<div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:1rem;
            margin-bottom:1rem;">
  <div style="border:1px solid var(--border);border-radius:var(--r);
              padding:12px 14px;background:white;">
    <div style="font-size:10px;color:var(--text-muted);margin-bottom:4px;">검사자 (Analyst)</div>
    <div id="sig-analyst" style="font-size:13px;font-weight:500;min-height:24px;">—</div>
    <div style="border-top:1px solid var(--border);margin-top:12px;padding-top:4px;
                font-size:9px;color:var(--text-muted);">서명 / Signature</div>
  </div>
  <div style="border:1px solid var(--border);border-radius:var(--r);
              padding:12px 14px;background:white;">
    <div style="font-size:10px;color:var(--text-muted);margin-bottom:4px;">검토자 (Reviewer)</div>
    <div id="sig-reviewer" style="font-size:13px;font-weight:500;min-height:24px;">—</div>
    <div style="border-top:1px solid var(--border);margin-top:12px;padding-top:4px;
                font-size:9px;color:var(--text-muted);">서명 / Signature</div>
  </div>
  <div style="border:1px solid var(--border);border-radius:var(--r);
              padding:12px 14px;background:white;">
    <div style="font-size:10px;color:var(--text-muted);margin-bottom:4px;">책임자 (Director)</div>
    <div id="sig-director" style="font-size:13px;font-weight:500;min-height:24px;">—</div>
    <div style="border-top:1px solid var(--border);margin-top:12px;padding-top:4px;
                font-size:9px;color:var(--text-muted);">서명 / Signature</div>
  </div>
</div>

<!-- Footer bar -->
<div class="footer">
  <div>
    <span id="ft-institution">Genecurix Inc.</span>
  </div>
  <div id="ft-right" style="text-align:right;"></div>
</div>
"""
