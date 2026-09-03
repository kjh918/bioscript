"""
S08 — 검사 소견 (2컬럼)
LEFT:  Syndrome 판정 테이블 (전체 그룹별)
RIGHT: CN / BAF 탭 + Plotly 그래프
JS: renderSynTable(), renderCnvTabs(), loadCnv(), drawPlot()
"""

SECTION_FINDINGS = """
<!-- ═══════════════════════════════════════════════════════ S08 FINDINGS -->
<div class="cols">

  <!-- LEFT: Syndrome 판정 테이블 -->
  <div>
    <div class="card" style="height:fit-content;">
      <div class="ch">
        <span class="ct">Syndrome 판정 결과</span>
        <span style="font-size:10px;color:var(--text-muted);">클릭 시 상세 보기</span>
      </div>
      <div id="syn-summary-bar"
           style="display:flex;gap:8px;padding:8px 12px;border-bottom:1px solid var(--border);
                  background:var(--surface2);flex-wrap:wrap;">
        <span id="syn-count-ab" style="font-size:11px;font-weight:700;color:var(--red);"></span>
        <span id="syn-count-su" style="font-size:11px;font-weight:700;color:var(--amber);"></span>
        <span id="syn-count-ok" style="font-size:11px;font-weight:700;color:var(--teal);"></span>
      </div>
      <table class="stbl">
        <thead>
          <tr>
            <th>Syndrome</th><th>판정</th><th>CN Value</th><th>Target Region</th>
          </tr>
        </thead>
        <tbody id="syn-tbody"></tbody>
      </table>
    </div>
  </div>

  <!-- RIGHT: Chromosome detail + CN/BAF -->
  <div>

    <!-- S07 Chromosome Detail (chromosome 클릭 시 표시) -->
    <div class="card" id="chrom-detail-card" style="display:none;margin-bottom:1rem;">
      <div class="ch">
        <div style="display:flex;align-items:center;gap:8px;">
          <span class="ct">Chromosome Detail</span>
          <span id="chrom-detail-badge"
                style="font-family:var(--mono);font-weight:700;font-size:15px;
                       color:var(--text);margin-left:4px;"></span>
        </div>
        <div style="display:flex;align-items:center;gap:6px;" class="no-print">
          <span style="font-size:10px;color:var(--text-muted);">Annotation:</span>
          <button id="annot-btn-band"    onclick="setAnnotMode('band')"
                  style="padding:2px 8px;font-size:10px;font-weight:600;border-radius:99px;
                         cursor:pointer;border:1px solid var(--navy);
                         background:var(--navy);color:white;">Band</button>
          <button id="annot-btn-proband" onclick="setAnnotMode('proband')"
                  style="padding:2px 8px;font-size:10px;font-weight:600;border-radius:99px;
                         cursor:pointer;border:1px solid var(--border);
                         background:var(--surface);color:var(--text-muted);">Proband</button>
          <button id="annot-btn-both"    onclick="setAnnotMode('both')"
                  style="padding:2px 8px;font-size:10px;font-weight:600;border-radius:99px;
                         cursor:pointer;border:1px solid var(--border);
                         background:var(--surface);color:var(--text-muted);">Both</button>
        </div>
      </div>
      <!-- Syndrome chips -->
      <div id="chrom-detail-syn-chips"
           style="display:flex;gap:6px;flex-wrap:wrap;padding:6px 14px;
                  border-bottom:1px solid var(--border);background:var(--surface2);min-height:32px;">
      </div>
      <!-- Horizontal ideogram + brush -->
      <div id="chrom-detail-wrap"
           style="width:100%;overflow-x:auto;overflow-y:visible;
                  padding:10px 10px 6px;background:white;min-height:100px;"></div>
      <!-- Marker region 클릭 버튼 목록 -->
      <div id="marker-legend-placeholder"
           style="display:flex;gap:6px;flex-wrap:wrap;padding:6px 14px;
                  border-top:1px solid var(--border);background:white;min-height:0;">
      </div>
      <!-- Brush range -->
      <div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap;
                  padding:5px 14px 8px;border-top:1px solid var(--border);
                  background:var(--surface2);">
        <span style="font-size:10px;color:var(--text-muted);">선택 구간:</span>
        <span id="brush-range"
              style="font-family:var(--mono);font-size:11px;
                     color:var(--navy);font-weight:600;">전체 염색체</span>
        <span style="font-size:10px;color:var(--text-muted);">
          — brush를 드래그하면 CN/BAF 그래프가 해당 구간으로 확대됩니다
        </span>
        <button onclick="resetBrush()" class="no-print"
                style="margin-left:auto;padding:2px 8px;font-size:10px;cursor:pointer;
                       border:1px solid var(--border);border-radius:4px;
                       background:var(--surface);color:var(--text-muted);">초기화</button>
      </div>
    </div>

    <!-- CN / BAF Analysis -->
    <div class="card">
      <div class="ch">
        <span class="ct">CN / BAF Analysis</span>
        <span id="cnv-status"
              style="font-size:10px;color:var(--text-muted);font-family:var(--mono);"></span>
      </div>
      <div class="cb">
        <div class="chrom-tabs" id="chrom-tabs"></div>
        <div id="cnv-panels"></div>
      </div>
    </div>

  </div>
</div>
"""
