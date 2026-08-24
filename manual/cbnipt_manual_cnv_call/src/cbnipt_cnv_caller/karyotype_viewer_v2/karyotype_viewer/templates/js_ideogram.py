"""
JS_IDEOGRAM

전체 핵형 제거 버전.

HIGH_RISK / SUSPECTED chromosome만 horizontal ideogram으로 표시.

구조
────────────────────────────────────────────────────────

renderIdeogram()
    │
    ├─ getAbnormalChromosomes()
    │
    └─ renderAbnormalChromosomes()
           │
           └─ createAbnormalChromCard()
                   │
                   ├─ renderSingleAbnormalIdeogram()
                   │
                   └─ renderAbnormalMarkerLegend()

"""


JS_IDEOGRAM = r"""

// ═════════════════════════════════════════════════════════════════════════════
// Constants
// ═════════════════════════════════════════════════════════════════════════════

var IDEO_CHROM_ORDER = [
  '1','2','3','4','5','6','7','8','9','10','11','12',
  '13','14','15','16','17','18','19','20','21','22',
  'X','Y'
];


// ═════════════════════════════════════════════════════════════════════════════
// Utility
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// normalizeChrom
//
// chr21 → 21
// 21    → 21
// chrX  → X
// ─────────────────────────────────────────────────────────────────────────────

function normalizeChrom(chrom) {

  if (chrom === null || chrom === undefined) {
    return '';
  }

  return String(chrom)
    .replace(/^chr/i, '')
    .trim()
    .toUpperCase();
}


// ─────────────────────────────────────────────────────────────────────────────
// chromosome order
// ─────────────────────────────────────────────────────────────────────────────

function chromosomeOrder(chrom) {

  chrom = normalizeChrom(chrom);

  var idx = IDEO_CHROM_ORDER.indexOf(chrom);

  return idx >= 0 ? idx : 999;
}


// ─────────────────────────────────────────────────────────────────────────────
// getCallColor
// ─────────────────────────────────────────────────────────────────────────────

function getCallColor(call) {

  if (call === 'HIGH_RISK') {
    return '#E53E3E';
  }

  if (call === 'SUSPECTED') {
    return '#DD6B20';
  }

  if (call === 'LOW_RISK') {
    return '#38A169';
  }

  return '#718096';
}


// ─────────────────────────────────────────────────────────────────────────────
// getCallPillClass
// ─────────────────────────────────────────────────────────────────────────────

function getCallPillClass(call) {

  if (call === 'HIGH_RISK') {
    return 'phr';
  }

  if (call === 'SUSPECTED') {
    return 'pmo';
  }

  if (call === 'LOW_RISK') {
    return 'plr';
  }

  return 'pnc';
}


// ─────────────────────────────────────────────────────────────────────────────
// escapeHtml
// ─────────────────────────────────────────────────────────────────────────────

function escapeHtml(value) {

  return String(value == null ? '' : value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#039;');
}


// ═════════════════════════════════════════════════════════════════════════════
// Syndrome / chromosome mapping
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// syndromeHasChromosome
//
// syndrome의 primary_chrom 또는 feature.chrom에 해당 chromosome이 있는지 확인
// ─────────────────────────────────────────────────────────────────────────────

function syndromeHasChromosome(syndrome, chrom) {

  chrom = normalizeChrom(chrom);

  var primary = normalizeChrom(
    syndrome ? syndrome.primary_chrom : ''
  );

  if (primary && primary === chrom) {
    return true;
  }

  var features = syndrome && syndrome.features
    ? syndrome.features
    : [];

  return features.some(function(f) {

    return normalizeChrom(f.chrom) === chrom;

  });
}


// ─────────────────────────────────────────────────────────────────────────────
// getAbnormalSyndromesForChrom
//
// 특정 chromosome의 HIGH_RISK / SUSPECTED syndrome
// ─────────────────────────────────────────────────────────────────────────────

function getAbnormalSyndromesForChrom(chrom) {

  chrom = normalizeChrom(chrom);

  return (MANIFEST.syndromes || [])
    .filter(function(s) {

      if (
        s.call !== 'HIGH_RISK' &&
        s.call !== 'SUSPECTED'
      ) {
        return false;
      }

      return syndromeHasChromosome(s, chrom);

    })
    .sort(function(a, b) {

      // HIGH_RISK 우선
      if (
        a.call === 'HIGH_RISK' &&
        b.call !== 'HIGH_RISK'
      ) {
        return -1;
      }

      if (
        b.call === 'HIGH_RISK' &&
        a.call !== 'HIGH_RISK'
      ) {
        return 1;
      }

      return String(a.syndrome || '')
        .localeCompare(String(b.syndrome || ''));

    });
}


// ─────────────────────────────────────────────────────────────────────────────
// getAbnormalChromosomes
//
// HIGH_RISK / SUSPECTED가 존재하는 chromosome만 반환
//
// 예:
// [
//   {chrom:'21', syndromes:[...]},
//   {chrom:'22', syndromes:[...]}
// ]
// ─────────────────────────────────────────────────────────────────────────────

function getAbnormalChromosomes() {

  var chromSet = {};

  (MANIFEST.syndromes || []).forEach(function(s) {

    if (
      s.call !== 'HIGH_RISK' &&
      s.call !== 'SUSPECTED'
    ) {
      return;
    }


    // primary chromosome
    var primary = normalizeChrom(
      s.primary_chrom
    );

    if (primary) {
      chromSet[primary] = true;
    }


    // features chromosome
    (s.features || []).forEach(function(f) {

      var c = normalizeChrom(
        f.chrom
      );

      if (c) {
        chromSet[c] = true;
      }

    });

  });


  return Object.keys(chromSet)
    .sort(function(a, b) {

      return chromosomeOrder(a) -
             chromosomeOrder(b);

    })
    .map(function(chrom) {

      return {
        chrom: chrom,
        syndromes: getAbnormalSyndromesForChrom(chrom)
      };

    });
}


// ═════════════════════════════════════════════════════════════════════════════
// Marker builder
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// buildChromMarkers
//
// HIGH_RISK / SUSPECTED syndrome의 feature만 annotation으로 사용.
//
// Feature type:
//   TargetChromosome
//   PrimaryTargetRegion
//   CoreRegion
//   CoreGene
//
// 색상:
//   TargetChromosome / PrimaryTargetRegion
//       HIGH_RISK  = red
//       SUSPECTED  = orange
//
//   CoreRegion = blue
//   CoreGene   = purple
// ─────────────────────────────────────────────────────────────────────────────

function buildChromMarkers(chrom) {

  chrom = normalizeChrom(chrom);

  var markers = [];


  (MANIFEST.syndromes || []).forEach(function(s) {

    if (
      s.call !== 'HIGH_RISK' &&
      s.call !== 'SUSPECTED'
    ) {
      return;
    }


    var callColor =
      getCallColor(s.call);


    (s.features || []).forEach(function(f) {

      var featureChrom =
        normalizeChrom(f.chrom);

      if (featureChrom !== chrom) {
        return;
      }


      var start = Number(f.start);
      var end   = Number(f.end);


      if (!Number.isFinite(start)) {
        return;
      }

      if (!Number.isFinite(end)) {
        end = start + 1;
      }

      if (end < start) {

        var tmp = start;
        start = end;
        end   = tmp;

      }


      var featureColor =
        callColor;


      if (f.type === 'CoreRegion') {

        featureColor = '#3182CE';

      } else if (f.type === 'CoreGene') {

        featureColor = '#6B46C1';

      }


      markers.push({

        chromosome: chrom,

        syndrome:
          s.syndrome || '',

        nipt_id:
          s.nipt_id || '',

        call:
          s.call || '',

        cn_value:
          s.cn_value,

        feat_name:
          f.name || f.type || 'Marker',

        feat_type:
          f.type || '',

        start:
          start,

        end:
          Math.max(start + 1, end),

        color:
          featureColor,

        call_color:
          callColor

      });

    });

  });


  // feature priority
  var rank = {
    TargetChromosome:    0,
    PrimaryTargetRegion: 1,
    CoreRegion:          2,
    CoreGene:            3
  };


  markers.sort(function(a, b) {

    var ra =
      rank[a.feat_type] !== undefined
        ? rank[a.feat_type]
        : 9;

    var rb =
      rank[b.feat_type] !== undefined
        ? rank[b.feat_type]
        : 9;


    if (ra !== rb) {
      return ra - rb;
    }


    if (a.start !== b.start) {
      return a.start - b.start;
    }


    return String(a.feat_name)
      .localeCompare(String(b.feat_name));

  });


  return markers;
}


// ═════════════════════════════════════════════════════════════════════════════
// Main render
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// renderIdeogram
//
// JS_RENDER.render()에서 호출
// ─────────────────────────────────────────────────────────────────────────────

function renderIdeogram() {

  var list =
    document.getElementById(
      'abnormal-ideogram-list'
    );

  var normalMessage =
    document.getElementById(
      'karyotype-normal-message'
    );


  if (!list) {
    return;
  }


  list.innerHTML = '';


  if (normalMessage) {
    normalMessage.style.display = 'none';
  }


  // Ideogram.js loading check
  if (typeof Ideogram === 'undefined') {

    var tries = 0;


    var timer = setInterval(function() {

      tries++;


      if (typeof Ideogram !== 'undefined') {

        clearInterval(timer);

        renderAbnormalChromosomes();

      } else if (tries >= 50) {

        clearInterval(timer);

        list.innerHTML =
          '<div style="' +
            'padding:20px;' +
            'font-size:11px;' +
            'text-align:center;' +
            'color:var(--text-muted);' +
          '">' +
            'Ideogram.js 로드 실패' +
          '</div>';

      }

    }, 100);


    return;
  }


  renderAbnormalChromosomes();
}


// ═════════════════════════════════════════════════════════════════════════════
// Abnormal chromosome renderer
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// renderAbnormalChromosomes
// ─────────────────────────────────────────────────────────────────────────────

function renderAbnormalChromosomes() {

  var list =
    document.getElementById(
      'abnormal-ideogram-list'
    );

  var normalMessage =
    document.getElementById(
      'karyotype-normal-message'
    );

  var countEl =
    document.getElementById(
      'karyotype-count'
    );


  if (!list) {
    return;
  }


  list.innerHTML = '';


  var abnormalChroms =
    getAbnormalChromosomes();


  // ───────────────────────────────────────────────────────────────────────────
  // No abnormal chromosome
  // ───────────────────────────────────────────────────────────────────────────

  if (!abnormalChroms.length) {

    list.style.display = 'none';


    if (normalMessage) {
      normalMessage.style.display = '';
    }


    if (countEl) {
      countEl.textContent = '0 abnormal chromosomes';
    }


    return;
  }


  list.style.display = 'flex';


  if (normalMessage) {
    normalMessage.style.display = 'none';
  }


  if (countEl) {

    countEl.textContent =
      abnormalChroms.length +
      ' abnormal chromosome' +
      (abnormalChroms.length > 1 ? 's' : '');

  }


  abnormalChroms.forEach(
    function(item, idx) {

      createAbnormalChromCard(
        item.chrom,
        item.syndromes,
        idx,
        list
      );

    }
  );
}


// ═════════════════════════════════════════════════════════════════════════════
// Chromosome card
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// createAbnormalChromCard
// ─────────────────────────────────────────────────────────────────────────────

function createAbnormalChromCard(
  chrom,
  syndromes,
  idx,
  parent
) {

  chrom =
    normalizeChrom(chrom);


  // ───────────────────────────────────────────────────────────────────────────
  // Card
  // ───────────────────────────────────────────────────────────────────────────

  var card =
    document.createElement('div');


  card.className =
    'abnormal-chrom-card';


  card.style.cssText =
    'border:1px solid var(--border);' +
    'border-radius:7px;' +
    'background:white;' +
    'overflow:hidden;' +
    'page-break-inside:avoid;';


  // ───────────────────────────────────────────────────────────────────────────
  // Header
  // ───────────────────────────────────────────────────────────────────────────

  var header =
    document.createElement('div');


  header.style.cssText =
    'display:flex;' +
    'align-items:center;' +
    'gap:8px;' +
    'flex-wrap:wrap;' +
    'padding:7px 10px;' +
    'border-bottom:1px solid var(--border);' +
    'background:var(--surface2);';


  // chromosome name
  var chromLabel =
    document.createElement('span');


  chromLabel.textContent =
    'chr' + chrom;


  chromLabel.title =
    'chr' + chrom + ' CN/BAF 보기';


  chromLabel.style.cssText =
    'font-family:var(--mono);' +
    'font-size:14px;' +
    'font-weight:700;' +
    'color:var(--navy);' +
    'cursor:pointer;' +
    'min-width:42px;';


  chromLabel.addEventListener(
    'click',
    function() {

      onClickChromosome(chrom);

    }
  );


  header.appendChild(
    chromLabel
  );


  // vertical separator
  var separator =
    document.createElement('span');


  separator.style.cssText =
    'width:1px;' +
    'height:18px;' +
    'background:var(--border);' +
    'display:inline-block;';


  header.appendChild(
    separator
  );


  // ───────────────────────────────────────────────────────────────────────────
  // Syndrome chips
  // ───────────────────────────────────────────────────────────────────────────

  (syndromes || []).forEach(
    function(s) {

      var color =
        getCallColor(s.call);


      var pillClass =
        getCallPillClass(s.call);


      var cnText = '';

      if (
        s.cn_value !== null &&
        s.cn_value !== undefined &&
        Number.isFinite(Number(s.cn_value))
      ) {

        cnText =
          ' · CN ' +
          Number(s.cn_value).toFixed(3);

      }


      var chip =
        document.createElement('span');


      chip.style.cssText =
        'display:inline-flex;' +
        'align-items:center;' +
        'gap:5px;' +
        'padding:2px 7px;' +
        'border-radius:99px;' +
        'border:1px solid ' +
          color +
          '66;' +
        'background:' +
          color +
          '0D;' +
        'font-size:10px;' +
        'cursor:pointer;';


      chip.innerHTML =
        '<span style="' +
          'font-weight:700;' +
          'color:var(--text);' +
        '">' +
          escapeHtml(s.syndrome || '') +
        '</span>' +

        '<span class="pill ' +
          pillClass +
          '" style="' +
            'font-size:8px;' +
            'padding:1px 5px;' +
          '">' +
          escapeHtml(s.call || '') +
        '</span>' +

        (
          cnText
            ? '<span style="' +
                'font-family:var(--mono);' +
                'font-size:9px;' +
                'color:' +
                color +
                ';font-weight:600;' +
              '">' +
                escapeHtml(cnText) +
              '</span>'
            : ''
        );


      chip.title =
        (s.syndrome || '') +
        ' · ' +
        (s.call || '') +
        cnText;


      chip.addEventListener(
        'click',
        function() {

          onClickChromosome(chrom);

        }
      );


      header.appendChild(
        chip
      );

    }
  );


  card.appendChild(
    header
  );


  // ───────────────────────────────────────────────────────────────────────────
  // Ideogram container
  // ───────────────────────────────────────────────────────────────────────────

  var containerId =
    'abnormal-ideo-' +
    chrom +
    '-' +
    idx;


  var ideoWrap =
    document.createElement('div');


  ideoWrap.id =
    containerId;


  ideoWrap.style.cssText =
    'width:100%;' +
    'overflow-x:auto;' +
    'overflow-y:visible;' +
    'padding:7px 12px 4px;' +
    'background:white;' +
    'min-height:105px;';


  card.appendChild(
    ideoWrap
  );


  // ───────────────────────────────────────────────────────────────────────────
  // Annotation marker legend
  // ───────────────────────────────────────────────────────────────────────────

  var markerLegend =
    document.createElement('div');


  markerLegend.style.cssText =
    'display:flex;' +
    'align-items:center;' +
    'gap:5px;' +
    'flex-wrap:wrap;' +
    'padding:5px 10px 7px;' +
    'border-top:1px dashed var(--border);' +
    'background:#fff;';


  card.appendChild(
    markerLegend
  );


  parent.appendChild(
    card
  );


  // DOM이 생성된 후 Ideogram 렌더
  setTimeout(
    function() {

      renderSingleAbnormalIdeogram(
        chrom,
        containerId
      );


      renderAbnormalMarkerLegend(
        chrom,
        markerLegend
      );

    },
    0
  );
}


// ═════════════════════════════════════════════════════════════════════════════
// Ideogram renderer
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// renderSingleAbnormalIdeogram
// ─────────────────────────────────────────────────────────────────────────────

function renderSingleAbnormalIdeogram(
  chrom,
  containerId
) {

  chrom =
    normalizeChrom(chrom);


  if (
    typeof Ideogram === 'undefined'
  ) {
    return;
  }


  var markers =
    buildChromMarkers(chrom);


  var annotList = [];


  markers.forEach(
    function(m) {

      var trackIdx = 0;


      if (
        m.feat_type === 'CoreRegion'
      ) {

        trackIdx = 1;

      } else if (
        m.feat_type === 'CoreGene'
      ) {

        trackIdx = 2;

      } else {

        // TargetChromosome /
        // PrimaryTargetRegion
        trackIdx = 0;

      }


      annotList.push([

        m.feat_name,

        Number(m.start),

        Math.max(
          1,
          Number(m.end) -
          Number(m.start)
        ),

        trackIdx,

        m.color

      ]);

    }
  );


  // ───────────────────────────────────────────────────────────────────────────
  // Config
  // ───────────────────────────────────────────────────────────────────────────

  var cfg = {

    organism:
      'human',

    assembly:
      'GRCh38',

    container:
      '#' + containerId,

    orientation:
      'horizontal',

    chromosomes:
      [chrom],

    rows:
      1,

    // horizontal chromosome 길이
    chrHeight:
      650,

    chrWidth:
      18,

    chrMargin:
      18,

    rotatable:
      false,

    showBandLabels:
      true,

    showChromosomeLabels:
      true,

    resolution:
      850,

    showAnnotTooltip:
      true

  };


  // ───────────────────────────────────────────────────────────────────────────
  // Annotation
  // ───────────────────────────────────────────────────────────────────────────

  if (annotList.length) {

    cfg.annotations = {

      keys: [
        'name',
        'start',
        'length',
        'trackIndex',
        'color'
      ],

      annots: [
        {
          chr: chrom,
          annots: annotList
        }
      ]

    };


    cfg.annotationsLayout =
      'tracks';


    cfg.annotationHeight =
      9;


    cfg.annotationTracks = [

      {
        id:
          'target',

        displayName:
          'Target Region',

        color:
          '#E53E3E',

        shape:
          'rectangle'
      },

      {
        id:
          'core-region',

        displayName:
          'Core Region',

        color:
          '#3182CE',

        shape:
          'rectangle'
      },

      {
        id:
          'gene',

        displayName:
          'Core Gene',

        color:
          '#6B46C1',

        shape:
          'circle'
      }

    ];

  }


  try {

    new Ideogram(
      cfg
    );

  } catch(e) {

    console.error(
      'Abnormal chromosome ideogram error:',
      chrom,
      e
    );


    var wrap =
      document.getElementById(
        containerId
      );


    if (wrap) {

      wrap.innerHTML =
        '<div style="' +
          'padding:12px;' +
          'font-size:10px;' +
          'color:var(--text-muted);' +
        '">' +
          'chr' +
          escapeHtml(chrom) +
          ' ideogram 렌더링 실패' +
        '</div>';

    }

  }
}


// ═════════════════════════════════════════════════════════════════════════════
// Marker legend
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// renderAbnormalMarkerLegend
// ─────────────────────────────────────────────────────────────────────────────

function renderAbnormalMarkerLegend(
  chrom,
  element
) {

  chrom =
    normalizeChrom(chrom);


  if (!element) {
    return;
  }


  element.innerHTML = '';


  var markers =
    buildChromMarkers(chrom);


  // TargetChromosome은 annotation에서 표현하되
  // marker chip에서는 너무 큰 의미 없는 label이면 제외 가능.
  //
  // 현재는 모든 marker 표시.
  var displayMarkers =
    markers;


  if (!displayMarkers.length) {

    element.style.display =
      'none';

    return;
  }


  element.style.display =
    'flex';


  var label =
    document.createElement('span');


  label.textContent =
    'Annotations:';


  label.style.cssText =
    'font-size:9px;' +
    'font-weight:600;' +
    'color:var(--text-muted);' +
    'margin-right:2px;';


  element.appendChild(
    label
  );


  displayMarkers.forEach(
    function(m) {

      var startMb =
        Number(m.start) / 1000000;


      var endMb =
        Number(m.end) / 1000000;


      var chip =
        document.createElement('span');


      var markerShape = '';

      if (
        m.feat_type === 'CoreGene'
      ) {

        markerShape =
          '<span style="' +
            'width:7px;' +
            'height:7px;' +
            'border-radius:50%;' +
            'background:' +
            m.color +
            ';display:inline-block;' +
            'flex-shrink:0;' +
          '"></span>';

      } else {

        markerShape =
          '<span style="' +
            'width:9px;' +
            'height:6px;' +
            'border-radius:2px;' +
            'background:' +
            m.color +
            ';display:inline-block;' +
            'flex-shrink:0;' +
          '"></span>';

      }


      chip.innerHTML =

        markerShape +

        '<strong>' +
          escapeHtml(m.feat_name) +
        '</strong>' +

        '<span style="' +
          'color:var(--text-muted);' +
          'font-weight:400;' +
        '">' +
          escapeHtml(m.feat_type) +
        '</span>';


      chip.title =

        (m.syndrome || '') +

        '\nchr' +
        chrom +
        ':' +
        Number(m.start).toLocaleString() +
        '-' +
        Number(m.end).toLocaleString() +

        '\n' +
        startMb.toFixed(3) +
        '-' +
        endMb.toFixed(3) +
        ' Mb';


      chip.style.cssText =
        'display:inline-flex;' +
        'align-items:center;' +
        'gap:4px;' +
        'padding:2px 6px;' +
        'border-radius:99px;' +
        'font-size:9px;' +
        'font-weight:600;' +
        'border:1px solid ' +
          m.color +
          '44;' +
        'background:' +
          m.color +
          '0D;' +
        'color:' +
          m.color +
        ';';


      element.appendChild(
        chip
      );

    }
  );
}


// ═════════════════════════════════════════════════════════════════════════════
// Chromosome click
// ═════════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// onClickChromosome
//
// S07 renderChromDetail() 호출하지 않음.
//
// 기존 S08 CN/BAF chromosome tab으로 이동.
// ─────────────────────────────────────────────────────────────────────────────

function onClickChromosome(chrom) {

  var c =
    normalizeChrom(chrom);


  if (!c) {
    return;
  }


  if (
    typeof switchToChrom === 'function'
  ) {

    switchToChrom(c);

  }

}

"""