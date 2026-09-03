"""
S04 — 검체 / 분석 정보 + QC
Sample ID, 성별, ISCN, Fetal Sex
QC: 회수된 세포 수(mapped reads), Fetal Fraction, 분석 파이프라인
JS: renderSampleInfo()  ← MANIFEST.sample + MANIFEST.qc
"""

SECTION_SAMPLE = """
<!-- ═══════════════════════════════════════════════════════ S04 SAMPLE & QC -->
<div style="display:grid;grid-template-columns:1fr 1fr;gap:1rem;margin-bottom:1rem;">

  <!-- 검체 정보 -->
  <div class="card" style="margin-bottom:0;">
    <div class="ch">
      <span class="ct">검체 / 분석 정보</span>
      <span style="font-size:10px;color:var(--text-muted);">Sample Information</span>
    </div>
    <div class="igrid" id="sample-grid">
      <div class="iitem">
        <div class="il">Sample ID</div>
        <div class="iv" id="s-id" style="font-family:var(--mono);font-size:12px;">—</div>
      </div>
      <div class="iitem">
        <div class="il">태아 성별 (Fetal Sex)</div>
        <div class="iv" id="s-fetal-sex">—</div>
      </div>
      <div class="iitem" style="grid-column:1/-1;">
        <div class="il">핵형 (ISCN)</div>
        <div class="iv" id="s-iscn"
             style="font-family:var(--mono);font-size:11px;word-break:break-all;">—</div>
      </div>
      <div class="iitem">
        <div class="il">이상 소견 수</div>
        <div class="iv" id="s-events">—</div>
      </div>
      <div class="iitem">
        <div class="il">분석 파이프라인</div>
        <div class="iv" id="s-pipeline">cbNIPT v2</div>
      </div>
    </div>
  </div>

  <!-- QC 정보 -->
  <div class="card" style="margin-bottom:0;">
    <div class="ch">
      <span class="ct">QC 정보</span>
      <span style="font-size:10px;color:var(--text-muted);">Quality Control</span>
    </div>
    <div class="igrid" id="qc-grid">
      <div class="iitem">
        <div class="il">회수된 세포 수 (Mapped Reads)</div>
        <div class="iv" id="qc-mapped-reads">—</div>
      </div>
      <div class="iitem">
        <div class="il">태아 분율 (Fetal Fraction)</div>
        <div class="iv" id="qc-fetal-fraction">—</div>
      </div>
      <div class="iitem">
        <div class="il">GC Bias</div>
        <div class="iv" id="qc-gc-bias">—</div>
      </div>
      <div class="iitem">
        <div class="il">QC 결과</div>
        <div class="iv" id="qc-pass">—</div>
      </div>
      <div class="iitem">
        <div class="il">중앙값 NCV</div>
        <div class="iv" id="qc-ncv">—</div>
      </div>
      <div class="iitem">
        <div class="il">분석일</div>
        <div class="iv" id="qc-analysis-date">—</div>
      </div>
    </div>
  </div>

</div>
"""
