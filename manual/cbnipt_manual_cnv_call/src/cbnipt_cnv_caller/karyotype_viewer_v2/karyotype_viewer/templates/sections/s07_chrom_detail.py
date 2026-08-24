"""
S07 — Single Chromosome Detail Card
chromosome 클릭 시 표시:
  - annotation toggle (기본 band / proband annotation)
  - horizontal ideogram + brush (CN 구간 파싱)
  - syndrome chips
  - brush 범위 표시
JS: renderChromDetail(), onBrushMove(), toggleAnnotation()
"""

SECTION_CHROM_DETAIL = """
<!-- ═══════════════════════════════════════════════════════ S07 CHROM DETAIL -->
<div class="card" id="chrom-detail-card" style="display:none;margin-bottom:1rem;">

  <!-- Card header -->
  <div class="ch">
    <div style="display:flex;align-items:center;gap:8px;">
      <span class="ct">Chromosome Detail</span>
      <span id="chrom-detail-badge"
            style="font-family:var(--mono);font-weight:700;font-size:15px;
                   color:var(--text);margin-left:4px;"></span>
    </div>
    <!-- Annotation toggle -->
    <div style="display:flex;align-items:center;gap:6px;" class="no-print">
      <span style="font-size:10px;color:var(--text-muted);">Annotation:</span>
      <button id="annot-btn-band"
              onclick="setAnnotMode('band')"
              style="padding:3px 10px;font-size:10px;font-weight:600;border-radius:99px;
                     cursor:pointer;border:1px solid var(--navy);
                     background:var(--navy);color:white;">
        Band
      </button>
      <button id="annot-btn-proband"
              onclick="setAnnotMode('proband')"
              style="padding:3px 10px;font-size:10px;font-weight:600;border-radius:99px;
                     cursor:pointer;border:1px solid var(--border);
                     background:var(--surface);color:var(--text-muted);">
        Proband
      </button>
      <button id="annot-btn-both"
              onclick="setAnnotMode('both')"
              style="padding:3px 10px;font-size:10px;font-weight:600;border-radius:99px;
                     cursor:pointer;border:1px solid var(--border);
                     background:var(--surface);color:var(--text-muted);">
        Both
      </button>
    </div>
  </div>

  <!-- Syndrome chips for this chromosome -->
  <div id="chrom-detail-syn-chips"
       style="display:flex;gap:6px;flex-wrap:wrap;
              padding:6px 14px;border-bottom:1px solid var(--border);
              background:var(--surface2);min-height:32px;">
  </div>

  <!-- Horizontal ideogram + brush -->
  <div id="chrom-detail-wrap"
       style="width:100%;overflow-x:auto;overflow-y:visible;
              padding:10px 10px 6px;background:white;min-height:100px;">
  </div>

  <!-- Brush range + hint -->
  <div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap;
              padding:5px 14px 8px;border-top:1px solid var(--border);
              background:var(--surface2);">
    <span style="font-size:10px;color:var(--text-muted);">선택 구간:</span>
    <span id="brush-range"
          style="font-family:var(--mono);font-size:11px;
                 color:var(--navy);font-weight:600;">전체 염색체</span>
    <span style="font-size:10px;color:var(--text-muted);">
      — brush를 드래그하면 아래 CN/BAF 그래프가 해당 구간으로 확대됩니다
    </span>
    <button onclick="resetBrush()"
            class="no-print"
            style="margin-left:auto;padding:2px 8px;font-size:10px;cursor:pointer;
                   border:1px solid var(--border);border-radius:4px;
                   background:var(--surface);color:var(--text-muted);">
      초기화
    </button>
  </div>

</div>
"""
