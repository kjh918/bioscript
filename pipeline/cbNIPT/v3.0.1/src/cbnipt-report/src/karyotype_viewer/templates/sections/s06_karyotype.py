"""
S06 — Abnormal Karyotype Overview

전체 karyotype는 표시하지 않음.

HIGH_RISK / SUSPECTED가 존재하는 chromosome만 표시:
  - chromosome horizontal ideogram
  - syndrome / call / CN
  - PrimaryTargetRegion
  - CoreRegion
  - CoreGene

chromosome 또는 syndrome chip 클릭:
  → 기존 CN/BAF chromosome plot으로 이동

S07 Chromosome Detail은 사용하지 않음.
"""


SECTION_KARYOTYPE = """
<!-- ═══════════════════════════════════════════════════════ S06 KARYOTYPE -->
<div class="card" id="karyotype-card">

  <!-- Header -->
  <div class="ch">

    <div style="display:flex;align-items:center;gap:8px;">
      <span class="ct">Karyotype Overview</span>

      <span id="karyotype-count"
            style="
              font-family:var(--mono);
              font-size:10px;
              font-weight:600;
              color:var(--text-muted);
            ">
      </span>
    </div>

    <span class="no-print"
          style="
            font-size:10px;
            color:var(--text-muted);
          ">
      이상 의심 염색체만 표시합니다
    </span>

  </div>


  <!-- HIGH_RISK / SUSPECTED chromosome ideograms -->
  <div id="abnormal-ideogram-list"
       style="
         display:flex;
         flex-direction:column;
         gap:10px;
         padding:10px 12px 12px;
         background:white;
         min-height:80px;
       ">
  </div>


  <!-- No abnormal finding -->
  <div id="karyotype-normal-message"
       style="
         display:none;
         padding:24px 16px;
         text-align:center;
         color:var(--text-muted);
         font-size:12px;
         background:white;
       ">

    검사 대상 질환에서 HIGH RISK 또는 SUSPECTED 염색체가 없습니다.

  </div>


  <!-- Legend -->
  <div style="
    display:flex;
    gap:16px;
    flex-wrap:wrap;
    align-items:center;
    padding:7px 14px 10px;
    border-top:1px solid var(--border);
    background:var(--surface2);
    font-size:10px;
  ">

    <span style="
      color:var(--text-muted);
      font-weight:600;
    ">
      범례:
    </span>


    <!-- HIGH RISK -->
    <span style="display:inline-flex;align-items:center;gap:4px;">

      <span style="
        display:inline-block;
        width:12px;
        height:8px;
        background:#E53E3E;
        border-radius:2px;
      "></span>

      HIGH RISK

    </span>


    <!-- SUSPECTED -->
    <span style="display:inline-flex;align-items:center;gap:4px;">

      <span style="
        display:inline-block;
        width:12px;
        height:8px;
        background:#DD6B20;
        border-radius:2px;
      "></span>

      SUSPECTED

    </span>


    <!-- Core Region -->
    <span style="display:inline-flex;align-items:center;gap:4px;">

      <span style="
        display:inline-block;
        width:12px;
        height:8px;
        background:#3182CE;
        border-radius:2px;
      "></span>

      Core Region

    </span>


    <!-- Core Gene -->
    <span style="display:inline-flex;align-items:center;gap:4px;">

      <span style="
        display:inline-block;
        width:10px;
        height:10px;
        background:#6B46C1;
        border-radius:50%;
      "></span>

      Core Gene

    </span>

  </div>

</div>
"""