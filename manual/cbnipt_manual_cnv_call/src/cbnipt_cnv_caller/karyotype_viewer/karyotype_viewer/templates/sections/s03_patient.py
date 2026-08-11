"""
S03 — 의뢰기관 정보 + 검사대상자 정보 + 정도관리
레이아웃 참조: cbNIPT 결과보고서 양식
"""

SECTION_PATIENT = """
<!-- ═══════════════════════════════════════════════════════ S03 PATIENT -->

<!-- 의뢰기관 / 환자 정보 (3컬럼 × 5행 단일 그리드) -->
<div class="card" style="margin-bottom:0.5rem;">
  <div style="display:grid;grid-template-columns:1fr 1fr 1fr;grid-template-rows:repeat(5,auto);">

    <!-- Row 1 -->
    <div class="iitem">
      <div class="il">의뢰기관</div>
      <div class="iv" id="pt-institution">—</div>
    </div>
    <div class="iitem">
      <div class="il">기관기호</div>
      <div class="iv" id="pt-institution-code">—</div>
    </div>
    <div class="iitem" style="border-right:none;">
      <div class="il">접수번호</div>
      <div class="iv" id="pt-accession">—</div>
    </div>

    <!-- Row 2 -->
    <div class="iitem">
      <div class="il">성명</div>
      <div class="iv" id="pt-name" style="font-size:18px;font-weight:700;">—</div>
    </div>
    <div class="iitem">
      <div class="il">주민번호</div>
      <div class="iv" id="pt-dob" style="font-size:18px;font-weight:700;">—</div>
    </div>
    <div class="iitem" style="border-right:none;">
      <div class="il">검체채취일</div>
      <div class="iv" id="pt-collection-date">—</div>
    </div>

    <!-- Row 3 -->
    <div class="iitem">
      <div class="il">등록번호</div>
      <div class="iv" id="pt-reg-no">—</div>
    </div>
    <div class="iitem">
      <div class="il">진료과</div>
      <div class="iv" id="pt-department">—</div>
    </div>
    <div class="iitem" style="border-right:none;">
      <div class="il">검사의뢰일</div>
      <div class="iv" id="pt-request-date">—</div>
    </div>

    <!-- Row 4 -->
    <div class="iitem">
      <div class="il">검사종류</div>
      <div class="iv" id="pt-panel">—</div>
    </div>
    <div class="iitem">
      <div class="il">병동/주치의</div>
      <div class="iv" id="pt-physician">—</div>
    </div>
    <div class="iitem" style="border-right:none;">
      <div class="il">검사일</div>
      <div class="iv" id="pt-test-date">—</div>
    </div>

    <!-- Row 5 -->
    <div class="iitem" style="border-bottom:none;">
      <div class="il">임상정보/기타</div>
      <div class="iv" id="pt-clinical-info">—</div>
    </div>
    <div class="iitem" style="border-bottom:none;">
      <div class="il">임신주수/체중</div>
      <div class="iv" id="pt-ga-weight">—</div>
    </div>
    <div class="iitem" style="border-right:none;border-bottom:none;">
      <div class="il">결과보고일</div>
      <div class="iv" id="pt-report-date">—</div>
    </div>

  </div>
</div>

<!-- 검사대상자 정보 + 정도관리 (2컬럼) -->
<div style="display:grid;grid-template-columns:1fr 1fr;gap:0.5rem;margin-bottom:1rem;">

  <!-- 검사대상자 정보 -->
  <div class="card" style="margin-bottom:0;">
    <div class="ch" style="justify-content:center;">
      <span class="ct" style="font-size:11px;letter-spacing:.05em;">검사대상자 정보</span>
    </div>
    <div style="display:grid;grid-template-columns:repeat(5,1fr);">
      <div class="iitem" style="text-align:center;border-bottom:none;">
        <div class="il" style="text-align:center;">초음파<br>특이소견</div>
        <div class="iv" id="pt-ultrasound" style="text-align:center;font-size:12px;">—</div>
      </div>
      <div class="iitem" style="text-align:center;border-bottom:none;">
        <div class="il" style="text-align:center;">목덜미투명대<br>(NT)</div>
        <div class="iv" id="pt-nt" style="text-align:center;font-size:12px;">—</div>
      </div>
      <div class="iitem" style="text-align:center;border-bottom:none;">
        <div class="il" style="text-align:center;">모체혈청<br>선별검사 소견</div>
        <div class="iv" id="pt-maternal-serum" style="text-align:center;font-size:12px;">—</div>
      </div>
      <div class="iitem" style="text-align:center;border-bottom:none;">
        <div class="il" style="text-align:center;">시험관아기<br>시술 여부</div>
        <div class="iv" id="pt-ivf" style="text-align:center;font-size:12px;">—</div>
      </div>
      <div class="iitem" style="text-align:center;border-right:none;border-bottom:none;">
        <div class="il" style="text-align:center;">태아 수</div>
        <div class="iv" id="pt-fetus" style="text-align:center;font-size:12px;">—</div>
      </div>
    </div>
  </div>

  <!-- 정도관리 -->
  <div class="card" style="margin-bottom:0;">
    <div class="ch" style="justify-content:center;">
      <span class="ct" style="font-size:11px;letter-spacing:.05em;">정도관리</span>
    </div>
    <div style="display:grid;grid-template-columns:repeat(3,1fr);">
      <div class="iitem" style="text-align:center;border-bottom:none;">
        <div class="il" style="text-align:center;">CFC 회수</div>
        <div class="iv" id="qc-cfdna" style="text-align:center;font-size:12px;">—</div>
      </div>
      <div class="iitem" style="text-align:center;border-bottom:none;">
        <div class="il" style="text-align:center;">NGS 검사<br>데이터 품질</div>
        <div class="iv" id="qc-ngs" style="text-align:center;font-size:12px;">—</div>
      </div>
      <div class="iitem" style="text-align:center;border-right:none;border-bottom:none;">
        <div class="il" style="text-align:center;">표준물질<br>검사결과</div>
        <div class="iv" id="qc-ref-material" style="text-align:center;font-size:12px;">—</div>
      </div>
    </div>
  </div>

</div>
"""