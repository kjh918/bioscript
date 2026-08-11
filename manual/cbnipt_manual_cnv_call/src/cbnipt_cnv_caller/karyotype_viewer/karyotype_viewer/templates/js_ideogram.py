"""
JS_IDEOGRAM: Ideogram.js 관련 함수들.

brush API (ideogram.js 1.47):
  - onBrushMove / onBrushEnd 는 인수 없이 호출됨
  - 결과는 ideogram 인스턴스의 selectedRegion.from / .to 에 저장
  - 따라서 detailIdeogram.selectedRegion 을 직접 읽어야 함
"""

JS_IDEOGRAM = r"""
// ── Ideogram state ────────────────────────────────────────────────────────
var currentChrom      = null;
var detailIdeogram    = null;
var currentAnnotMode  = 'proband';   // 기본값: proband annotation 표시
var overviewAnnotMode = 'band';

// ── buildAnnotations ──────────────────────────────────────────────────────
// MANIFEST.raw_annots + syndrome call 색상 오버레이
function buildAnnotations() {
  var ra = MANIFEST.raw_annots;
  if (!ra || !ra.annots) return null;
  var callColor = {};
  (MANIFEST.syndromes || []).forEach(function(s) {
    var c = s.call === 'HIGH_RISK'   ? '#E53E3E' :
            s.call === 'SUSPECTED' ? '#DD6B20' : '#38A169';
    (s.features || []).forEach(function(f) {
      callColor[f.chrom + '_' + f.start] = c;
    });
  });
  return {
    keys: ra.keys,
    annots: (ra.annots || []).map(function(chrObj) {
      return {
        chr: chrObj.chr,
        annots: (chrObj.annots || []).map(function(a) {
          var ov = callColor[chrObj.chr + '_' + a[1]];
          return ov ? [a[0], a[1], a[2], a[3], ov] : a;
        }),
      };
    }),
  };
}

// ── buildChromMarkers ─────────────────────────────────────────────────────
// 특정 염색체의 syndrome marker feature 목록 반환
function buildChromMarkers(chrom) {
  var markers = [];
  (MANIFEST.syndromes || []).forEach(function(s) {
    (s.features || []).forEach(function(f) {
      if (String(f.chrom) !== String(chrom)) return;
      markers.push({
        syndrome:  s.syndrome,
        nipt_id:   s.nipt_id,
        call:      s.call,
        cn_value:  s.cn_value,
        feat_name: f.name,
        feat_type: f.type,
        start:     f.start,
        end:       f.end,
        color:     s.call === 'HIGH_RISK'   ? '#E53E3E' :
                   s.call === 'SUSPECTED' ? '#DD6B20' : '#38A169',
      });
    });
  });
  // PrimaryTargetRegion / TargetChromosome 우선, CoreGene 나중
  markers.sort(function(a, b) {
    var rank = {TargetChromosome:0, PrimaryTargetRegion:1, CoreRegion:2, CoreGene:3};
    return (rank[a.feat_type] || 9) - (rank[b.feat_type] || 9);
  });
  return markers;
}

// ── _setToggleBtnStyle ────────────────────────────────────────────────────
function _setToggleBtnStyle(prefix, active) {
  ['band','proband','both'].forEach(function(m) {
    var btn = document.getElementById(prefix + '-btn-' + m);
    if (!btn) return;
    var on = (m === active);
    btn.style.background  = on ? 'var(--navy)' : 'var(--surface)';
    btn.style.color       = on ? 'white'       : 'var(--text-muted)';
    btn.style.borderColor = on ? 'var(--navy)' : 'var(--border)';
  });
}

// ── renderIdeogram (overview) ─────────────────────────────────────────────
function renderIdeogram() {
  if (typeof Ideogram === 'undefined') {
    document.getElementById('ideogram-wrap').innerHTML =
      '<p style="color:var(--text-muted);padding:20px;font-size:12px;">ideogram.js 로드 실패</p>';
    return;
  }
  _renderOverviewWithMode(overviewAnnotMode);
}

function _renderOverviewWithMode(mode) {
  var wrap = document.getElementById('ideogram-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';
  wrap.dataset.clickInstalled = '';

  var S = MANIFEST.sample;
  // overview: annotation 없이 band only (proband는 detail에서만)
  var config = {
    organism:             'human',
    assembly:             'GRCh38',
    container:            '#ideogram-wrap',
    orientation:          'vertical',
    chromosomes:          S.display_chroms,
    rows:                 1,
    chrHeight:            220,
    chrWidth:             12,
    chrMargin:            8,
    rotatable:            false,
    showBandLabels:       false,
    showChromosomeLabels: true,
    resolution:           550,
    showAnnotTooltip:     true,
  };
  try {
    new Ideogram(config);
    _installChromClickHandler();
  } catch(e) { console.error('Overview ideogram error:', e); }
}

function setOverviewAnnotMode(mode) {
  overviewAnnotMode = mode;
  _setToggleBtnStyle('ov-annot', mode);
  _renderOverviewWithMode(mode);
}

// ── _installChromClickHandler ─────────────────────────────────────────────
function _installChromClickHandler() {
  var ALLOWED = new Set([
    '1','2','3','4','5','6','7','8','9','10','11','12',
    '13','14','15','16','17','18','19','20','21','22','X','Y'
  ]);
  function parseChrom(val) {
    if (!val) return null;
    var m = String(val).match(/(?:chr)?([0-9]{1,2}|X|Y)(?:[^A-Za-z0-9]|$)/i);
    return (m && ALLOWED.has(m[1].toUpperCase())) ? m[1].toUpperCase() : null;
  }
  function chromFromNode(target) {
    var root = document.getElementById('ideogram-wrap');
    var node = target;
    while (node && node !== root) {
      var cands = [
        node.id,
        node.getAttribute && node.getAttribute('class'),
        node.getAttribute && node.getAttribute('data-chr'),
      ];
      for (var i = 0; i < cands.length; i++) {
        var c = parseChrom(cands[i]); if (c) return c;
      }
      if (node.tagName && node.tagName.toLowerCase() === 'text') {
        var ct = parseChrom(node.textContent); if (ct) return ct;
      }
      node = node.parentElement;
    }
    var g = target && target.closest ? target.closest('g') : null;
    if (g) {
      var texts = g.querySelectorAll('text');
      for (var j = 0; j < texts.length; j++) {
        var cv = parseChrom(texts[j].textContent); if (cv) return cv;
      }
    }
    return null;
  }
  var tries = 0;
  var timer = setInterval(function() {
    tries++;
    var wrap = document.getElementById('ideogram-wrap');
    if (!wrap || wrap.dataset.clickInstalled === '1') {
      if (wrap && wrap.dataset.clickInstalled === '1') clearInterval(timer);
      if (tries > 100) clearInterval(timer);
      return;
    }
    if (!wrap.querySelector('svg')) { if (tries > 100) clearInterval(timer); return; }
    wrap.dataset.clickInstalled = '1';
    wrap.style.cursor = 'pointer';
    wrap.addEventListener('click', function(ev) {
      var chrom = chromFromNode(ev.target);
      if (chrom) onClickChromosome(chrom);
    }, true);
    clearInterval(timer);
  }, 150);
}

// ── onClickChromosome ─────────────────────────────────────────────────────
function onClickChromosome(chrom) {
  var c = String(chrom).replace('chr','').toUpperCase();
  if (!c) return;
  var el = document.getElementById('ideo-selected');
  if (el) el.textContent = 'chr' + c + ' selected';
  switchToChrom(c);
  renderChromDetail(c);
}

// ── renderChromDetail ─────────────────────────────────────────────────────
function renderChromDetail(chrom) {
  chrom = String(chrom).replace('chr','').toUpperCase();
  if (!chrom) return;
  currentChrom = chrom;

  var card  = document.getElementById('chrom-detail-card');
  var badge = document.getElementById('chrom-detail-badge');
  if (card)  card.style.display = '';
  if (badge) badge.textContent  = 'chr' + chrom;

  // Syndrome chips
  var chipsEl = document.getElementById('chrom-detail-syn-chips');
  if (chipsEl) {
    var syns = (MANIFEST.syndromes || []).filter(function(s) {
      return s.primary_chrom === chrom;
    });
    if (!syns.length) {
      chipsEl.innerHTML =
        '<span style="font-size:11px;color:var(--text-muted);">이 염색체에 마커 없음 (LOW_RISK)</span>';
    } else {
      chipsEl.innerHTML = syns.map(function(s) {
        var pc = {HIGH_RISK:'phr',SUSPECTED:'pmo',LOW_RISK:'plr',UNKNOWN:'pnc'}[s.call]||'pnc';
        var cn = s.cn_value != null ? ' CN ' + s.cn_value.toFixed(3) : '';
        return '<span style="font-size:11px;font-weight:600;cursor:pointer;" ' +
               'onclick="zoomToSyndrome(\'' + s.nipt_id + '\')">' +
               s.syndrome + ' <span class="pill ' + pc + '">' + s.call + '</span>' + cn + '</span>';
      }).join('');
    }
  }

  _renderDetailIdeogram(chrom, currentAnnotMode);
  _renderMarkerLegend(chrom);
}

// ── _renderMarkerLegend ───────────────────────────────────────────────────
function _renderMarkerLegend(chrom) {
  var placeholder = document.getElementById('marker-legend-placeholder');
  if (!placeholder) return;

  var markers = buildChromMarkers(chrom).filter(function(m) {
    return m.feat_type === 'PrimaryTargetRegion' || m.feat_type === 'CoreGene';
  });
  if (!markers.length) { placeholder.innerHTML = ''; return; }

  placeholder.innerHTML =
    '<span style="font-size:10px;color:var(--text-muted);align-self:center;">마커 클릭:</span>' +
    markers.map(function(m) {
      var c = m.color;
      return '<span onclick="zoomToRegion(' + m.start + ',' + m.end + ')" ' +
             'title="' + m.syndrome + '\n' + (m.start/1e6).toFixed(2) + '–' + (m.end/1e6).toFixed(2) + ' Mb" ' +
             'style="display:inline-flex;align-items:center;gap:4px;cursor:pointer;' +
             'padding:2px 8px;border-radius:99px;font-size:10px;font-weight:600;' +
             'background:' + c + '18;border:1px solid ' + c + '88;color:' + c + ';">' +
             '<span style="width:7px;height:7px;border-radius:2px;background:' + c + ';display:inline-block;"></span>' +
             m.feat_name + '</span>';
    }).join('');
}

// ── _renderDetailIdeogram ─────────────────────────────────────────────────
function _renderDetailIdeogram(chrom, mode) {
  var wrap = document.getElementById('chrom-detail-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  if (typeof Ideogram === 'undefined') {
    wrap.innerHTML = '<p style="padding:12px;color:var(--text-muted);">ideogram.js 없음</p>';
    return;
  }

  var chromSize = (MANIFEST.chrom_sizes || {})[chrom] || 250000000;

  // Syndrome marker annotation (PrimaryTargetRegion + CoreGene)
  var markers  = buildChromMarkers(chrom);
  var annotObj = null;
  if (markers.length && mode !== 'band') {
    var annotList = markers.map(function(m) {
      // keys: name, start, length, trackIndex, color
      var trackIdx = (m.feat_type === 'CoreGene' || m.feat_type === 'CoreRegion') ? 1 : 0;
      return [m.feat_name, m.start, Math.max(1, m.end - m.start), trackIdx, m.color];
    });
    annotObj = {
      keys:   ['name','start','length','trackIndex','color'],
      annots: [{chr: chrom, annots: annotList}],
    };
  }

  var cfg = {
    organism:             'human',
    assembly:             'GRCh38',
    container:            '#chrom-detail-wrap',
    orientation:          'horizontal',
    chromosomes:          [chrom],
    rows:                 1,
    chrHeight:            600,
    chrWidth:             20,
    chrMargin:            20,
    rotatable:            false,
    showBandLabels:       true,
    showChromosomeLabels: true,
    resolution:           850,
    showAnnotTooltip:     true,
    brush:                'chr' + chrom + ':1-' + chromSize,
    onBrushMove:          _onBrushCallback,
    onBrushEnd:           _onBrushCallback,
  };

  if (annotObj) {
    cfg.annotations       = annotObj;
    cfg.annotationsLayout = 'tracks';
    cfg.annotationHeight  = 10;
    cfg.annotationTracks  = [
      {id:'region', displayName:'Target Region', color:'#E53E3E', shape:'rectangle'},
      {id:'gene',   displayName:'Core Gene',      color:'#6B46C1', shape:'circle'},
    ];
  }

  try {
    detailIdeogram = new Ideogram(cfg);
  } catch(e) {
    console.error('Detail ideogram error:', e);
  }
}

// ── _onBrushCallback ──────────────────────────────────────────────────────
// ideogram.js 1.47: brush 콜백은 인수 없이 호출됨
// selectedRegion 은 ideogram 인스턴스에 저장됨
function _onBrushCallback() {
  if (!detailIdeogram || !currentChrom) return;
  var region = detailIdeogram.selectedRegion;
  if (!region) return;
  var s = region.from || 0;
  var e = region.to   || 0;
  if (e <= s) return;
  _applyBrushRange(s, e);
}

function _applyBrushRange(startBp, endBp) {
  var rangeEl = document.getElementById('brush-range');
  if (rangeEl) {
    rangeEl.textContent =
      'chr' + currentChrom + ': ' +
      (startBp / 1e6).toFixed(3) + ' – ' + (endBp / 1e6).toFixed(3) + ' Mb' +
      '  (' + ((endBp - startBp) / 1e6).toFixed(3) + ' Mb)';
  }

  // plotDiv는 flex wrapper 안의 첫 번째 자식 div (graphDiv)
  var outerDiv = document.getElementById('plot-' + currentChrom);
  if (!outerDiv || !window.Plotly) return;
  // drawPlot에서 plotDiv 안에 wrapper > graphDiv 구조로 만들어짐
  var graphDiv = outerDiv.querySelector('.js-plotly-plot') || outerDiv;

  Plotly.relayout(graphDiv, {
    'xaxis.range':        [startBp / 1e6, endBp / 1e6],
    'xaxis.autorange':    false,
    // y축 고정
    'yaxis.range':        [-0.1, 5.5],
    'yaxis.autorange':    false,
    'yaxis2.range':       [-0.05, 1.05],
    'yaxis2.autorange':   false,
  });
}

// ── zoomToSyndrome ────────────────────────────────────────────────────────
// syndrome chip 클릭 → primary region 확대
function zoomToSyndrome(nipt_id) {
  var syn = (MANIFEST.syndromes || []).find(function(s) { return s.nipt_id === nipt_id; });
  if (!syn) return;
  var prim = (syn.features || []).find(function(f) {
    return f.type === 'PrimaryTargetRegion' || f.type === 'TargetChromosome';
  }) || syn.features[0];
  if (!prim) return;
  zoomToRegion(prim.start, prim.end);
}

// ── zoomToRegion ──────────────────────────────────────────────────────────
// 마커 라벨 클릭 → CN/BAF 그래프 해당 구간 확대
function zoomToRegion(startBp, endBp) {
  _applyBrushRange(startBp, endBp);
}

// ── setAnnotMode ──────────────────────────────────────────────────────────
function setAnnotMode(mode) {
  currentAnnotMode = mode;
  _setToggleBtnStyle('annot', mode);
  if (currentChrom) _renderDetailIdeogram(currentChrom, mode);
}

// ── resetBrush ────────────────────────────────────────────────────────────
function resetBrush() {
  if (!currentChrom) return;
  var chromSize = (MANIFEST.chrom_sizes || {})[currentChrom] || 250000000;
  var rangeEl   = document.getElementById('brush-range');
  if (rangeEl) rangeEl.textContent = '전체 염색체';

  var outerDiv = document.getElementById('plot-' + currentChrom);
  if (!outerDiv || !window.Plotly) return;
  var graphDiv = outerDiv.querySelector('.js-plotly-plot') || outerDiv;
  Plotly.relayout(graphDiv, {
    'xaxis.range':      [0, chromSize / 1e6],
    'xaxis.autorange':  false,
    'yaxis.range':      [-0.1, 5.5],
    'yaxis.autorange':  false,
    'yaxis2.range':     [-0.05, 1.05],
    'yaxis2.autorange': false,
  });
}
"""
