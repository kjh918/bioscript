"""
JS_CNV: CN/BAF 탭 렌더링 + Plotly 그래프 함수들.
  - renderCnvTabs()   전체 염색체 탭 생성 (ALL_CHROMS 기준)
  - switchTab()       탭 전환
  - switchToChrom()   외부에서 염색체 지정 전환
  - loadCnv()         CNV 데이터 로드 (var 전역변수)
  - drawPlot()        Plotly CN + BAF scatter
"""

JS_CNV = r"""
// ── CNV state ─────────────────────────────────────────────────────────────
var cnvCache  = {};
var plotCache = {};

// ── renderCnvTabs ─────────────────────────────────────────────────────────
// 전체 염색체 탭 생성 — affected는 색상 강조, 나머지는 NORMAL 표시
function renderCnvTabs() {
  var M        = MANIFEST;
  var affected = M.affected_chroms || [];

  // syndrome 맵: chrom → worst call
  var chromCall = {};
  (M.syndromes || []).forEach(function(s) {
    var c    = s.primary_chrom;
    var rank = {HIGH_RISK:3, SUSPECTED:2, LOW_RISK:1, UNKNOWN:0};
    if ((rank[s.call] || 0) > (rank[chromCall[c]] || 0)) chromCall[c] = s.call;
  });

  // 표시 순서: affected 먼저, 나머지 순서대로
  var displayChroms = M.sample.display_chroms || ALL_CHROMS;
  var tabChroms = affected.concat(
    displayChroms.filter(function(c) { return affected.indexOf(c) === -1; })
  );

  var tabsEl   = document.getElementById('chrom-tabs');
  var panelsEl = document.getElementById('cnv-panels');

  tabChroms.forEach(function(chrom, i) {
    var call = chromCall[chrom] || 'LOW_RISK';
    var tc   = TAB_CLS[call] || 'nml';
    var isAffected = affected.indexOf(chrom) !== -1;

    // Tab button
    var btn      = document.createElement('button');
    btn.className = 'ctab ' + tc + (i === 0 ? ' active' : '');
    btn.textContent  = 'chr' + chrom;
    btn.dataset.chrom = chrom;
    btn.title    = call;
    btn.addEventListener('click', function() { switchTab(chrom); });
    tabsEl.appendChild(btn);

    // Panel
    var panel       = document.createElement('div');
    panel.className = 'cnv-panel' + (i === 0 ? ' active' : '');
    panel.id        = 'panel-' + chrom;

    // syndrome chips (affected만)
    var syns   = (M.syndromes || []).filter(function(s) { return s.primary_chrom === chrom; });
    var schips = '';
    syns.forEach(function(s) {
      var pc = PILL_CLS[s.call] || 'pnc';
      var cn = s.cn_value != null ? ' · CN ' + s.cn_value.toFixed(3) : '';
      var cl = s.call === 'HIGH_RISK'   ? 'var(--red)' :
               s.call === 'SUSPECTED' ? 'var(--amber)' : 'var(--teal)';
      schips +=
        '<span style="display:inline-flex;align-items:center;gap:5px;padding:4px 10px;' +
        'border-radius:99px;font-size:11px;font-weight:600;margin-right:4px;' +
        'background:var(--surface2);border:1px solid ' + cl + '66;color:' + cl + ';">' +
        s.syndrome + ' · <span class="pill ' + pc + '">' + s.call + '</span>' + cn + '</span>';
    });

    // primary region 요약
    var allFeats = syns.flatMap(function(s) { return s.features || []; })
      .filter(function(f) {
        return f.type === 'PrimaryTargetRegion' || f.type === 'TargetChromosome';
      });
    var regionStr = allFeats
      .map(function(f) {
        return f.name + ': ' + (f.start/1e6).toFixed(1) + '–' + (f.end/1e6).toFixed(1) + ' Mb';
      }).join('  |  ');

    panel.innerHTML =
      '<div class="syn-chips">' +
        (schips || '<span style="color:var(--text-muted);font-size:11px;">' +
          (isAffected ? '이상 소견' : 'LOW_RISK') + '</span>') +
      '</div>' +
      (regionStr ? '<div class="region-list">' + regionStr + '</div>' : '') +
      '<div id="plot-' + chrom + '" style="min-height:160px;">' +
        '<div class="loading"><div class="spinner"></div>데이터 로드 중…</div>' +
      '</div>';

    panelsEl.appendChild(panel);
  });

  // 첫 탭 로드
  if (tabChroms.length) loadCnv(tabChroms[0]);
}

// ── switchTab ─────────────────────────────────────────────────────────────
function switchTab(chrom) {
  document.querySelectorAll('.ctab').forEach(function(b) {
    b.classList.toggle('active', b.dataset.chrom === chrom);
  });
  document.querySelectorAll('.cnv-panel').forEach(function(p) {
    p.classList.toggle('active', p.id === 'panel-' + chrom);
  });
  loadCnv(chrom);
}

// ── switchToChrom ─────────────────────────────────────────────────────────
// 외부(ideogram 클릭, syndrome 행 클릭)에서 호출
function switchToChrom(chrom) {
  chrom = String(chrom).replace(/^chr/i, '');
  // data-chrom 대소문자 모두 시도
  var btn = document.querySelector('.ctab[data-chrom="' + chrom + '"]')
         || document.querySelector('.ctab[data-chrom="' + chrom.toUpperCase() + '"]')
         || document.querySelector('.ctab[data-chrom="' + chrom.toLowerCase() + '"]');
  if (!btn) return;
  switchTab(btn.dataset.chrom);
  btn.scrollIntoView({behavior:'smooth', block:'nearest', inline:'nearest'});
}

// ── loadCnv ───────────────────────────────────────────────────────────────
// var CNV_CHR{N} 전역변수에서 바로 읽기 (inline 모드)
// 없으면 data/cnv/chr{N}.js <script> 동적 로드 (static 모드)
function loadCnv(chrom) {
  if (plotCache[chrom]) return;

  var syns    = (MANIFEST.syndromes || []).filter(function(s) {
    return s.primary_chrom === chrom;
  });
  var plotDiv = document.getElementById('plot-' + chrom);
  if (!plotDiv) return;

  var varName = 'CNV_CHR' + chrom.toUpperCase();

  // 이미 전역변수로 있으면 바로 사용
  if (window[varName] && window[varName].length) {
    plotCache[chrom] = true;
    var rows = window[varName];
    cnvCache[chrom] = rows;
    var statusEl = document.getElementById('cnv-status');
    if (statusEl) statusEl.textContent = 'chr' + chrom + ' (' + rows.length.toLocaleString() + ' rows)';
    drawPlot(chrom, rows, syns, plotDiv);
    return;
  }

  // 전역변수 없음 → data/cnv/chrN.js 동적 로드
  var jsFile = 'data/cnv/chr' + chrom + '.js';
  var statusEl = document.getElementById('cnv-status');
  if (statusEl) statusEl.textContent = '로드 중: ' + jsFile;

  var s    = document.createElement('script');
  s.src    = jsFile;
  s.onload = function() {
    plotCache[chrom] = true;
    var rows = window[varName] || [];
    cnvCache[chrom] = rows;
    if (statusEl) {
      statusEl.textContent = rows.length
        ? 'chr' + chrom + ' (' + rows.length.toLocaleString() + ' rows)'
        : 'chr' + chrom + ' 데이터 없음';
    }
    if (rows.length) {
      drawPlot(chrom, rows, syns, plotDiv);
    } else {
      plotDiv.innerHTML =
        '<p style="color:var(--text-muted);padding:16px;font-size:12px;">데이터 없음</p>';
    }
  };
  s.onerror = function() {
    plotCache[chrom] = true;
    plotDiv.innerHTML =
      '<p style="color:var(--text-muted);padding:16px;font-size:12px;">' +
      'chr' + chrom + ' — CNV 데이터 없음 (LOW RISK 추정)</p>';
    if (statusEl) statusEl.textContent = '';
  };
  document.head.appendChild(s);
}

// ── drawPlot ──────────────────────────────────────────────────────────────
function drawPlot(chrom, rows, syns, plotDiv) {
  var posCols = ['pos','start','position','chromstart'];
  var cnCols  = ['cn','copy_number_signal','copy_number','copynumber','ratio','log2','observed_cn'];
  var bafCols = ['baf','bin_baf','b_allele_freq','allele_freq','vaf'];
  var posCol  = posCols.find(function(c) { return rows[0] && c in rows[0]; });
  var cnCol   = cnCols.find(function(c)  { return rows[0] && c in rows[0]; });
  var bafCol  = bafCols.find(function(c) { return rows[0] && c in rows[0]; });

  if (!posCol || !cnCol) {
    plotDiv.innerHTML = '<p style="padding:16px;color:var(--text-muted);">pos/cn 컬럼 없음</p>';
    return;
  }

  var pos     = rows.map(function(r) { return r[posCol] / 1e6; });
  var cn      = rows.map(function(r) { return r[cnCol]; });
  var baf     = bafCol ? rows.map(function(r) { return r[bafCol]; }) : null;
  var cnColor = cn.map(function(v) {
    return v >= 2.8 ? '#FC8181' : v <= 1.3 ? '#90CDF4' : '#4A5568';
  });

  // ── 마커 regions 수집 ─────────────────────────────────────────────────
  var shapes = [], markerList = [];
  (syns || []).forEach(function(s) {
    var callC = s.call === 'HIGH_RISK'   ? 'rgba(159,27,27,'  :
                s.call === 'SUSPECTED' ? 'rgba(122,72,0,'   :
                                          'rgba(10,110,92,';
    var hexC  = s.call === 'HIGH_RISK'   ? '#E53E3E' :
                s.call === 'SUSPECTED' ? '#DD6B20' : '#38A169';
    (s.features || []).forEach(function(f) {
      if (String(f.chrom) !== String(chrom)) return;
      var alpha = f.type === 'CoreGene'           ? 0.15 :
                  f.type === 'CoreRegion'          ? 0.12 :
                  f.type === 'PrimaryTargetRegion' ? 0.09 : 0.05;
      var x0 = f.start / 1e6, x1 = f.end / 1e6;
      shapes.push({
        type:'rect', xref:'x', yref:'paper',
        x0:x0, x1:x1, y0:0, y1:1,
        fillcolor: callC + alpha + ')',
        line: {
          color: callC + '0.6)',
          width: f.type === 'CoreGene' ? 1.0 : 0.7,
          dash:  f.type === 'CoreGene' ? 'dot' : 'solid',
        },
        layer:'below',
      });
      if (f.type !== 'TargetChromosome') {
        markerList.push({
          name:  f.name,
          type:  f.type,
          start: f.start, end: f.end,
          color: hexC,
          syn:   s.syndrome,
          call:  s.call,
        });
      }
    });
  });

  // CN 참조선
  [{v:1,c:'#90CDF4'},{v:2,c:'#68D391',w:1.5},{v:3,c:'#F6AD55'}].forEach(function(ref) {
    shapes.push({
      type:'line', xref:'x', yref:'y',
      x0:pos[0]||0, x1:pos[pos.length-1]||1,
      y0:ref.v, y1:ref.v,
      line:{color:ref.c, width:ref.w||0.8, dash:'dot'},
    });
  });

  // ── Plotly traces ──────────────────────────────────────────────────────
  var traces = [{
    x:pos, y:cn, mode:'markers', name:'CN',
    marker:{size:3, opacity:0.75, color:cnColor},
    hovertemplate:'%{x:.3f} Mb<br>CN=%{y:.3f}<extra></extra>',
    xaxis:'x', yaxis:'y',
  }];
  if (baf) {
    traces.push({
      x:pos, y:baf, mode:'markers', name:'BAF',
      marker:{size:3, opacity:0.6, color:'#805AD5'},
      hovertemplate:'%{x:.3f} Mb<br>BAF=%{y:.3f}<extra></extra>',
      xaxis:'x', yaxis:'y2',
    });
  }

  var axCommon = {
    showgrid:true, gridcolor:'#EDF2F7',
    tickfont:{size:9}, linecolor:'#CBD5E0', showline:true,
    fixedrange: true,   // y축 드래그/줌 고정
  };
  var layout = {
    xaxis: Object.assign({}, axCommon, {
      title:{text:'Genomic position (Mb)', font:{size:10}},
      fixedrange: false,  // x축은 자유롭게
    }),
    yaxis: Object.assign({}, axCommon, {
      title:{text:'Copy Number', font:{size:10}},
      range:[-0.1,5.5],
      domain: baf ? [0.45,1] : [0,1],
    }),
    shapes,
    margin:{l:54, r:16, t:20, b:42},
    height: baf ? 400 : 240,
    paper_bgcolor:'white', plot_bgcolor:'white',
    showlegend:false, hovermode:'closest',
    font:{family:'Inter,Arial,sans-serif', size:10},
  };
  if (baf) {
    layout.yaxis2 = Object.assign({}, axCommon, {
      title:{text:'BAF', font:{size:10}},
      range:[-0.05,1.05], domain:[0,0.4],
    });
  }

  // ── 그래프 + 라벨 패널을 flex 컨테이너로 감쌈 ────────────────────────
  plotDiv.innerHTML = '';
  var wrapper = document.createElement('div');
  wrapper.style.cssText = 'display:flex;align-items:flex-start;gap:0;';

  var graphDiv = document.createElement('div');
  graphDiv.style.cssText = 'flex:1;min-width:0;';
  wrapper.appendChild(graphDiv);

  // 마커 라벨 패널 (우측)
  if (markerList.length) {
    var labelPanel = document.createElement('div');
    labelPanel.style.cssText =
      'width:140px;flex-shrink:0;padding:' + (baf ? '8px' : '4px') + ' 6px;' +
      'font-size:10px;display:flex;flex-direction:column;gap:4px;';

    markerList.forEach(function(m) {
      var btn = document.createElement('div');
      btn.title = m.syn + ' · ' + m.call +
                  '\n' + (m.start/1e6).toFixed(2) + '–' + (m.end/1e6).toFixed(2) + ' Mb';
      btn.style.cssText =
        'display:flex;align-items:center;gap:4px;cursor:pointer;' +
        'padding:3px 6px;border-radius:4px;border:1px solid ' + m.color + '44;' +
        'background:' + m.color + '0d;' +
        'transition:background .1s;';
      btn.innerHTML =
        '<span style="width:7px;height:7px;border-radius:1px;flex-shrink:0;' +
        'background:' + m.color + ';display:inline-block;"></span>' +
        '<span style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap;' +
        'color:' + m.color + ';font-weight:600;">' + m.name + '</span>';
      btn.addEventListener('mouseenter', function() {
        btn.style.background = m.color + '22';
      });
      btn.addEventListener('mouseleave', function() {
        btn.style.background = m.color + '0d';
      });
      btn.addEventListener('click', function() {
        // graphDiv는 plotDiv 안의 Plotly div
        var gd = plotDiv.querySelector('.js-plotly-plot') || graphDiv;
        if (window.Plotly) {
          Plotly.relayout(gd, {
            'xaxis.range':      [m.start/1e6, m.end/1e6],
            'xaxis.autorange':  false,
            'yaxis.range':      [-0.1, 5.5],
            'yaxis.autorange':  false,
            'yaxis2.range':     [-0.05, 1.05],
            'yaxis2.autorange': false,
          });
        }
        if (typeof zoomToRegion === 'function') zoomToRegion(m.start, m.end);
      });
      labelPanel.appendChild(btn);
    });
    wrapper.appendChild(labelPanel);
  }

  plotDiv.appendChild(wrapper);
  Plotly.newPlot(graphDiv, traces, layout, {
    responsive:true, displayModeBar:true,
    modeBarButtonsToRemove:['select2d','lasso2d'],
    toImageButtonOptions:{format:'png', scale:2},
  });
}
"""
