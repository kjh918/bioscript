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
  HIGH_RISK:  'HIGH RISK — 고위험',
  SUSPECTED:'SUSPECTED — 추가 확인 필요',
  LOW_RISK:    'LOW RISK — 저위험',
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
  const M   = MANIFEST;
  const call = M.overall_call;
  const abn  = M.syndromes.filter(s => s.call === 'HIGH_RISK').map(s => s.syndrome);
  const sus  = M.syndromes.filter(s => s.call === 'SUSPECTED').map(s => s.syndrome);
  let desc   = '';
  if (abn.length) desc += '<strong>' + abn.join(', ') + '</strong> — 이상 소견 확인. ';
  if (sus.length) desc += '추가 확인 필요: ' + sus.join(', ') + '. ';
  if (!desc)      desc  = '검사한 전 항목에서 이상 소견이 관찰되지 않았습니다.';

  const el = document.getElementById('banner');
  if (!el) return;
  el.className = 'banner ' + call;
  el.innerHTML =
    '<div class="bicon">' + (CALL_ICON[call] || '?') + '</div>' +
    '<div><div class="bres">' + (CALL_TEXT[call] || call) + '</div>' +
    '<div class="bdesc">' + desc + '</div></div>';
}

// ── renderChips ───────────────────────────────────────────────────────────
function renderChips() {
  const M   = MANIFEST;
  const S   = M.sample;
  const abn = M.syndromes.filter(s => s.call === 'HIGH_RISK').map(s => s.syndrome);
  const sus = M.syndromes.filter(s => s.call === 'SUSPECTED').map(s => s.syndrome);
  const el  = document.getElementById('chips');
  if (!el) return;
  el.innerHTML = '';
  abn.forEach(n => el.innerHTML += '<span class="chip chip-ab">🔴 ' + n + ': ABNORMAL</span>');
  sus.forEach(n => el.innerHTML += '<span class="chip chip-su">🟠 ' + n + ': SUSPICIOUS</span>');
  if (!abn.length && !sus.length)
    el.innerHTML += '<span class="chip chip-ok">✅ 전 항목 LOW RISK</span>';
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
  const el = document.getElementById('ev-cards');
  MANIFEST.events.forEach(function(ev) {
    const c = ev.color || '#718096';
    el.innerHTML +=
      '<div class="ev-card" style="border:1px solid ' + c + '88;' +
      'border-left:4px solid ' + c + ';background:' + c + '0d;">' +
      '<div style="font-weight:700;font-family:var(--mono);font-size:13px;color:' + c + ';">' +
        ev.iscn + '</div>' +
      '<div style="font-size:11px;color:var(--text-sub);margin-top:2px;">' +
        'chr' + ev.chr + ' · ' + ev.type.replace(/_/g, ' ') + '</div>' +
      '<div style="font-size:11px;color:var(--text-muted);">CN = ' + ev.cn + '</div>' +
      '</div>';
  });
}

// ── renderSynTable ────────────────────────────────────────────────────────
function renderSynTable() {
  var tbody = document.getElementById('syn-tbody');
  GROUP_ORDER.forEach(function(grp) {
    var items = MANIFEST.syndromes.filter(function(s) { return s.group === grp; });
    if (!items.length) return;
    var gc = GROUP_COLOR[grp] || '#CBD5E0';

    // catrow — createElement 사용 (innerHTML+= 하면 기존 이벤트 리스너 사라짐)
    var catrow = document.createElement('tr');
    catrow.className = 'catrow';
    catrow.innerHTML = '<td colspan="4" style="border-left:3px solid ' + gc + ';">' + grp + '</td>';
    tbody.appendChild(catrow);

    items.sort(function(a, b) { return a.syndrome.localeCompare(b.syndrome); })
    .forEach(function(s) {
      var pc  = PILL_CLS[s.call] || 'pnc';
      var cn  = s.cn_value != null ? s.cn_value.toFixed(3) : '—';
      var prim = (s.features || []).find(function(f) {
        return f.type === 'TargetChromosome' || f.type === 'PrimaryTargetRegion';
      }) || (s.features || [])[0];
      var reg = prim
        ? 'chr' + prim.chrom + ' ' +
          (prim.start / 1e6).toFixed(1) + '–' + (prim.end / 1e6).toFixed(1) + ' Mb'
        : '—';

      var row = document.createElement('tr');
      row.innerHTML =
        '<td><div class="sname">' + s.syndrome + '</div>' +
            '<div class="ssub">' + s.nipt_id + '</div></td>' +
        '<td><span class="pill ' + pc + '">' + s.call + '</span></td>' +
        '<td style="font-family:var(--mono);font-size:11px;">' + cn + '</td>' +
        '<td style="font-size:10px;color:var(--text-muted);">' + reg + '</td>';

      var chrom = String(s.primary_chrom || '').replace(/^chr/i, '');
      row.style.cursor = 'pointer';
      row.addEventListener('click', function() {
        if (chrom) {
          switchToChrom(chrom);
          renderChromDetail(chrom);
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
  var Q = MANIFEST.qc     || {};

  // Sample fields
  var sFields = {
    's-id':         S.id   || '—',
    's-fetal-sex':  S.sex === 'female' ? '♀ Female' : S.sex === 'male' ? '♂ Male' : '—',
    's-iscn':       S.iscn || '—',
    's-events':     (MANIFEST.events || []).length + '개',
    's-pipeline':   S.pipeline || 'cbNIPT v2',
  };
  Object.entries(sFields).forEach(function([id, val]) {
    var el = document.getElementById(id);
    if (el) el.textContent = val;
  });

  // QC fields
  var qFields = {
    'qc-mapped-reads':  Q.mapped_reads != null
                          ? Number(Q.mapped_reads).toLocaleString() + ' reads'
                          : '—',
    'qc-fetal-fraction':Q.fetal_fraction != null
                          ? (Q.fetal_fraction * 100).toFixed(1) + ' %'
                          : '—',
    'qc-gc-bias':       Q.gc_bias != null ? Q.gc_bias.toFixed(4) : '—',
    'qc-pass':          Q.qc_pass === true  ? '✅ PASS'
                      : Q.qc_pass === false ? '❌ FAIL' : '—',
    'qc-ncv':           Q.median_ncv != null ? Q.median_ncv.toFixed(3) : '—',
    'qc-analysis-date': Q.analysis_date || (MANIFEST.report_date || '—'),
  };
  Object.entries(qFields).forEach(function([id, val]) {
    var el = document.getElementById(id);
    if (el) el.textContent = val;
  });
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
  if (abEl) abEl.textContent = ab ? '🔴 HIGH_RISK: ' + ab + '개' : '';
  if (suEl) suEl.textContent = su ? '🟠 SUSPECTED: ' + su + '개' : '';
  if (okEl) okEl.textContent = '✅ LOW_RISK: ' + ok + '개';
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
  if (hdrBrand && In.name) hdrBrand.textContent = In.name + ' · 임상유전체연구소';

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
