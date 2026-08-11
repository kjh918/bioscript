"""
JS_RENDER: 데이터 → DOM 렌더링 함수들.
  - render()           페이지 전체 초기 렌더
  - renderBanner()     종합 판정 배너
  - renderSampleGrid() Sample Information 그리드
  - renderEvCards()    Chromosomal Events 카드
  - renderSynTable()   Syndrome 판정 테이블
"""

JS_RENDER = r"""
// ── Constants ─────────────────────────────────────────────────────────────
const CALL_ICON  = {HIGH_RISK:'⚠', SUSPECTED:'〜', LOW_RISK:'✓', UNKNOWN:'?'};
const CALL_TEXT  = {
  HIGH_RISK: 'HIGH RISK — 고위험',
  SUSPECTED: 'SUSPECTED — 추가 확인 필요',
  LOW_RISK:  'LOW RISK — 저위험',
};
const PILL_CLS  = {HIGH_RISK:'phr', SUSPECTED:'pmo', LOW_RISK:'plr', UNKNOWN:'pnc'};
const TAB_CLS   = {HIGH_RISK:'abn', SUSPECTED:'sus', LOW_RISK:'nml'};
const GROUP_ORDER = [
  'Autosome Abnormality',
  'Sex Chromosome Abnormality',
  'Micro Deletion',
];
const GROUP_COLOR = {
  'Autosome Abnormality':       '#FC8181',
  'Sex Chromosome Abnormality': '#90CDF4',
  'Micro Deletion':             '#F6AD55',
};
const ALL_CHROMS = [
  '1','2','3','4','5','6','7','8','9','10','11',
  '12','13','14','15','16','17','18','19','20','21','22','X','Y'
];

// ── render ────────────────────────────────────────────────────────────────
function render() {
  var M = MANIFEST;
  var S = M.sample;

  document.title = 'Karyotype Report · ' + (S.id || '');

  renderBanner();          // S02
  renderChips();           // S02
  renderPatientInfo();     // S03
  renderSampleInfo();      // S04
  renderEvCards();         // S05
  renderIdeogram();        // S06
  renderSynSummaryBar();   // S08 summary bar
  renderSynTable();        // S08 left
  renderCnvTabs();         // S08 right
  renderFooter();          // S09
}

// ── renderBanner ──────────────────────────────────────────────────────────
function renderBanner() {
  var M    = MANIFEST;
  var call = M.overall_call;
  var abn  = M.syndromes.filter(function(s) { return s.call === 'HIGH_RISK'; }).map(function(s) { return s.syndrome; });
  var sus  = M.syndromes.filter(function(s) { return s.call === 'SUSPECTED'; }).map(function(s) { return s.syndrome; });
  var desc = '';
  if (abn.length) desc += '<strong>' + abn.join(', ') + '</strong> — 이상 소견 확인. ';
  if (sus.length) desc += '추가 확인 필요: ' + sus.join(', ') + '. ';
  if (!desc)      desc  = '검사한 전 항목에서 이상 소견이 관찰되지 않았습니다.';

  var el = document.getElementById('banner');
  if (!el) return;
  el.className = 'banner ' + call;
  el.innerHTML =
    '<div class="bicon">' + (CALL_ICON[call] || '?') + '</div>' +
    '<div><div class="bres">' + (CALL_TEXT[call] || call) + '</div>' +
    '<div class="bdesc">' + desc + '</div></div>';
}

// ── renderChips ───────────────────────────────────────────────────────────
function renderChips() {
  var M   = MANIFEST;
  var S   = M.sample;
  var abn = M.syndromes.filter(function(s) { return s.call === 'HIGH_RISK'; }).map(function(s) { return s.syndrome; });
  var sus = M.syndromes.filter(function(s) { return s.call === 'SUSPECTED'; }).map(function(s) { return s.syndrome; });
  var el  = document.getElementById('chips');
  if (!el) return;
  el.innerHTML = '';
  abn.forEach(function(n) { el.innerHTML += '<span class="chip chip-ab">🔴 ' + n + ': HIGH RISK</span>'; });
  sus.forEach(function(n) { el.innerHTML += '<span class="chip chip-su">🟠 ' + n + ': SUSPECTED</span>'; });
  if (!abn.length && !sus.length)
    el.innerHTML += '<span class="chip chip-ok">✅ 전 항목 LOW_RISK</span>';
  el.innerHTML +=
    '<span class="chip chip-ev">Sex: ' + (S.sex === 'female' ? '♀ Female' : '♂ Male') + '</span>' +
    '<span class="chip chip-ev">ISCN: ' + S.iscn + '</span>';
}

// ── renderSampleGrid ──────────────────────────────────────────────────────
function renderSampleGrid() {
  const S  = MANIFEST.sample;
  const el = document.getElementById('sample-grid');
  [
    ['Sample ID', S.id],
    ['Sex',       S.sex === 'female' ? '♀ Female' : '♂ Male'],
    ['ISCN',      S.iscn],
    ['CNV Events',MANIFEST.events.length],
  ].forEach(function([l, v]) {
    el.innerHTML +=
      '<div class="iitem"><div class="il">' + l + '</div>' +
      '<div class="iv" style="' +
        (l === 'ISCN' ? 'font-family:var(--mono);font-size:11px;' : '') +
      '">' + v + '</div></div>';
  });
}

// ── renderEvCards ─────────────────────────────────────────────────────────
function renderEvCards() {
  var el = document.getElementById('ev-cards');
  if (!el) return;
  var TYPE_LABEL = {
    trisomy:      '삼염색체 (Trisomy)',
    monosomy:     '단염색체 (Monosomy)',
    partial_gain: '부분 증가 (Partial Gain)',
    partial_loss: '부분 결실 (Partial Loss)',
  };
  el.style.cssText = 'display:grid;grid-template-columns:repeat(auto-fill,minmax(190px,1fr));gap:10px;';
  el.innerHTML = '';
  (MANIFEST.events || []).forEach(function(ev) {
    var c = ev.color || '#718096';
    var label = TYPE_LABEL[ev.type] || ev.type.replace(/_/g,' ');
    el.innerHTML +=
      '<div style="border:1px solid ' + c + '44;border-left:5px solid ' + c + ';' +
      'border-radius:8px;padding:13px 14px;background:' + c + '0d;">' +
        '<div style="font-size:26px;font-weight:800;font-family:var(--mono);' +
             'color:' + c + ';margin-bottom:5px;">' + ev.iscn + '</div>' +
        '<div style="font-size:11px;color:var(--text-sub);margin-bottom:9px;">' + label + '</div>' +
        '<div style="display:flex;gap:7px;align-items:center;">' +
          '<span style="font-size:11px;font-weight:700;font-family:var(--mono);' +
               'background:' + c + '22;padding:2px 7px;border-radius:4px;color:' + c + ';">chr' + ev.chr + '</span>' +
          '<span style="font-size:11px;color:var(--text-sub);font-family:var(--mono);">CN = ' + ev.cn + '</span>' +
        '</div>' +
      '</div>';
  });
}

// ── renderSynTable ────────────────────────────────────────────────────────
function renderSynTable() {
  const tbody = document.getElementById('syn-tbody');
  GROUP_ORDER.forEach(function(grp) {
    const items = MANIFEST.syndromes.filter(s => s.group === grp);
    if (!items.length) return;
    const gc = GROUP_COLOR[grp] || '#CBD5E0';
    tbody.innerHTML +=
      '<tr class="catrow"><td colspan="4" style="border-left:3px solid ' + gc + ';">' +
      grp + '</td></tr>';
    items.sort((a, b) => a.syndrome.localeCompare(b.syndrome)).forEach(function(s) {
      const pc   = PILL_CLS[s.call] || 'pnc';
      const cn   = s.cn_value != null ? s.cn_value.toFixed(3) : '—';
      const prim = (s.features || []).find(function(f) {
        return f.type === 'TargetChromosome' || f.type === 'PrimaryTargetRegion';
      }) || (s.features || [])[0];
      const reg  = prim
        ? 'chr' + prim.chrom + ' ' +
          (prim.start / 1e6).toFixed(1) + '–' + (prim.end / 1e6).toFixed(1) + ' Mb'
        : '—';

      const row = document.createElement('tr');
      row.innerHTML =
        '<td><div class="sname">' + s.syndrome + '</div>' +
            '<div class="ssub">' + s.nipt_id + '</div></td>' +
        '<td><span class="pill ' + pc + '">' + s.call + '</span></td>' +
        '<td style="font-family:var(--mono);font-size:11px;">' + cn + '</td>' +
        '<td style="font-size:10px;color:var(--text-muted);">' + reg + '</td>';

      // 행 클릭 → 해당 염색체 탭 + detail ideogram 전환
      row.addEventListener('click', function() {
        if (s.primary_chrom) {
          switchToChrom(s.primary_chrom);
          renderChromDetail(s.primary_chrom);
        }
      });
      tbody.appendChild(row);
    });
  });
}
"""

# ── S03 renderPatientInfo ──────────────────────────────────────────────────
JS_RENDER += r"""
// ── renderPatientInfo ─────────────────────────────────────────────────────
// ── renderPatientInfo ─────────────────────────────────────────────────────
function renderPatientInfo() {
  var M = MANIFEST.maternal || {};
  var Q = MANIFEST.qc       || {};
  var R = MANIFEST.report_date || '—';

  var fields = {
    // 컬럼1
    'pt-institution':     M.institution      || '—',
    'pt-name':            M.name             || '—',
    'pt-reg-no':          M.reg_no           || '—',
    'pt-panel':           M.panel            || 'cbNIPT',
    'pt-clinical-info':   M.clinical_info    || '—',
    // 컬럼2
    'pt-institution-code':M.institution_code || '—',
    'pt-dob':             M.dob              || '—',
    'pt-department':      M.department       || '—',
    'pt-physician':       M.physician        || '—',
    'pt-ga-weight':       (M.ga_weeks != null
                            ? M.ga_weeks + 'w' + (M.ga_days != null ? '+' + M.ga_days + 'd' : '')
                            : '—') +
                          ' / ' + (M.weight != null ? M.weight + 'kg' : '—'),
    // 컬럼3
    'pt-accession':       M.accession        || '—',
    'pt-collection-date': M.collection_date  || '—',
    'pt-request-date':    M.request_date     || '—',
    'pt-test-date':       M.test_date        || MANIFEST.report_date || '—',
    'pt-report-date':     R,
    // 검사대상자 정보
    'pt-ultrasound':      M.ultrasound       || '정보없음',
    'pt-nt':              M.nt               || '—',
    'pt-maternal-serum':  M.maternal_serum   || '없음',
    'pt-ivf':             M.ivf              || '없음',
    'pt-fetus':           M.fetus_count != null ? M.fetus_count + '태아' : '단태아',
    // 정도관리
    'qc-cfdna':           Q.cfdna_quality    || (Q.qc_pass === true ? '적합' : Q.qc_pass === false ? '부적합' : '—'),
    'qc-ngs':             Q.ngs_quality      || (Q.qc_pass === true ? '적합' : '—'),
    'qc-ref-material':    Q.ref_material     || '적합',
  };

  Object.entries(fields).forEach(function([id, val]) {
    var el = document.getElementById(id);
    if (el) el.textContent = val;
  });
}

// ── renderSampleInfo ──────────────────────────────────────────────────────
function renderSampleInfo() {
  var S = MANIFEST.sample || {};
  var M = MANIFEST.maternal || {};
  var call = MANIFEST.overall_call || 'LOW_RISK';

  var CALL_COLOR = {
    HIGH_RISK: {bg:'#FDF0F0', fg:'#9B1B1B', border:'#EFA5A5', icon:'⚠'},
    SUSPECTED: {bg:'#FEF8EC', fg:'#7A4800', border:'#F0C060', icon:'〜'},
    LOW_RISK:     {bg:'#E8F7F4', fg:'#0A6E5C', border:'#6EC4B5', icon:'✓'},
  };
  var cc = CALL_COLOR[call] || CALL_COLOR['LOW_RISK'];

  // 결과 요약 카드 배경 색상 반영
  var ch = document.getElementById('result-summary-ch');
  if (ch) { ch.style.background = cc.bg; ch.style.borderBottomColor = cc.border; }
  var card = document.getElementById('result-summary-card');
  if (card) { card.style.border = '1.5px solid ' + cc.border; }

  // 성별
  var isFemale = S.sex === 'female';
  var sexIcon  = document.getElementById('s-sex-icon');
  var sexText  = document.getElementById('s-fetal-sex');
  if (sexIcon) { sexIcon.textContent = isFemale ? '♀' : '♂';
                 sexIcon.style.color = isFemale ? '#E53E3E' : '#3182CE'; }
  if (sexText) { sexText.textContent = isFemale ? '여아 (Female)' : '남아 (Male)';
                 sexText.style.color = isFemale ? '#E53E3E' : '#3182CE'; }

  // 판정 배지
  var badge = document.getElementById('s-call-badge');
  var icon  = document.getElementById('s-call-icon');
  var text  = document.getElementById('s-call-text');
  if (badge) { badge.style.background = cc.bg; badge.style.borderColor = cc.border; }
  if (icon)  icon.textContent  = cc.icon;
  if (text)  { text.textContent = call; text.style.color = cc.fg; }

  // 이벤트 수
  var cnt = document.getElementById('s-events-count');
  if (cnt) {
    var n = (MANIFEST.events || []).length;
    cnt.textContent = n > 0 ? n + '개 이상 소견 검출' : '이상 소견 없음';
    cnt.style.color = n > 0 ? cc.fg : 'var(--text-muted)';
  }

  // ISCN / Sample ID
  var iscnEl = document.getElementById('s-iscn');
  var idEl   = document.getElementById('s-id');
  if (iscnEl) iscnEl.textContent = S.iscn || '—';
  if (idEl)   idEl.textContent   = S.id   || '—';
}

// ── renderSynSummaryBar ───────────────────────────────────────────────────
function renderSynSummaryBar() {
  var syns = MANIFEST.syndromes || [];
  var ab   = syns.filter(function(s) { return s.call === 'HIGH_RISK'; }).length;
  var su   = syns.filter(function(s) { return s.call === 'SUSPECTED'; }).length;
  var ok   = syns.filter(function(s) { return s.call === 'LOW_RISK'; }).length;

  var abEl = document.getElementById('syn-count-ab');
  var suEl = document.getElementById('syn-count-su');
  var okEl = document.getElementById('syn-count-ok');
  if (abEl) abEl.textContent = ab ? '🔴 HIGH RISK: ' + ab + '개' : '';
  if (suEl) suEl.textContent = su ? '🟠 SUSPECTED: ' + su + '개' : '';
  if (okEl) okEl.textContent = '✅ LOW RISK: ' + ok + '개';
}

// ── renderFooter ──────────────────────────────────────────────────────────
function renderFooter() {
  var M  = MANIFEST;
  var In = M.institution || {};
  var Sg = M.signatures  || {};

  var instEl = document.getElementById('ft-institution');
  if (instEl && In.name) instEl.textContent = In.name;

  var ctEl = document.getElementById('ft-contact');
  if (ctEl && In.contact) ctEl.textContent = In.contact;

  var frEl = document.getElementById('ft-right');
  if (frEl) frEl.textContent =
    'Pipeline: karyotype_viewer · 생성: ' + (M.report_date || '');

  var hdrBrand = document.getElementById('hdr-brand');
  if (hdrBrand && In.name) hdrBrand.textContent = In.name + '';

  // 서명란
  [['sig-analyst',  Sg.analyst ],
   ['sig-reviewer', Sg.reviewer],
   ['sig-director', Sg.director]].forEach(function([id, val]) {
    var el = document.getElementById(id);
    if (el && val) el.textContent = val;
  });

  // Header dates
  var testDateEl   = document.getElementById('hdr-test-date');
  var reportDateEl = document.getElementById('hdr-report-date');
  var hdrIdEl      = document.getElementById('hdr-id');
  if (testDateEl)   testDateEl.textContent   = M.test_date    || M.report_date || '—';
  if (reportDateEl) reportDateEl.textContent = M.report_date  || '—';
  if (hdrIdEl)      hdrIdEl.textContent      = M.report_id    || '—';
}
"""
