"""
S01 — Report Header
병원명(기관명), Report ID, 검사일, 보고일
JS에서 MANIFEST.report_id / report_date / institution 으로 채움.
"""

SECTION_HEADER = """
<!-- ═══════════════════════════════════════════════════════ S01 HEADER -->
<div class="hdr">
  <div>
    <div class="hdr-brand" id="hdr-brand">GENECURIX</div>
    <h1 id="hdr-title">순환태아세포 기반 염색체 이상 선별검사 결과보고서</h1>
    <div class="hdr-sub" id="hdr-sub">cbNIPT · Cell-based Non-Invasive Prenatal Testing</div>
  </div>
  <div class="hdr-right">
    <div class="hdr-lbl">Report No.</div>
    <div class="hdr-val" id="hdr-id">—</div>
    <div style="display:grid;grid-template-columns:1fr 1fr;gap:0 16px;margin-top:8px;">
      <div>
        <div class="hdr-lbl">검사일</div>
        <div class="hdr-val" style="font-size:11px;" id="hdr-test-date">—</div>
      </div>
      <div>
        <div class="hdr-lbl">보고일</div>
        <div class="hdr-val" style="font-size:11px;" id="hdr-report-date">—</div>
      </div>
    </div>
  </div>
</div>
"""
