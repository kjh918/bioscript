"""
Karyotype Viewer - FastAPI + 2 isolated Dash apps (v5)

Architecture
------------
FastAPI host :8050
  /overview/  -> Dash App 1 (whole karyotype, one Dash-Bio Ideogram)
  /detail/    -> Dash App 2 (single chromosome, one horizontal Dash-Bio Ideogram + brush)

Important changes
-----------------
1) No Flask session / no polling callback.
2) No iframe auto-height loop (prevents cumulative vertical drift).
3) Click a chromosome directly in Overview -> postMessage -> Detail Store.
4) Overview and Detail remain separate documents, so Ideogram.js instances do not share DOM.
5) Detail uses one horizontal Ideogram for chromosome + annotations + brush.

Install
-------
pip install fastapi uvicorn dash dash-bio plotly pandas numpy

Run
---
python karyotype_fastapi_two_dash_click_v5.py
Open http://localhost:8050
"""

import io
import json
import base64

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from starlette.middleware.wsgi import WSGIMiddleware

from dash import Dash, html, dcc, Input, Output, State, no_update, ctx
import dash_bio as dashbio


# =============================================================================
# Reference data
# =============================================================================
CHROM_SIZES = {
    '1':248956422,'2':242193529,'3':198295559,'4':190214555,
    '5':181538259,'6':170805979,'7':159345973,'8':145138636,
    '9':138394717,'10':133797422,'11':135086622,'12':133275309,
    '13':114364328,'14':107043718,'15':101991189,'16':90338345,
    '17':83257441,'18':80373285,'19':58617616,'20':64444167,
    '21':46709983,'22':50818468,'X':156040895,'Y':57227415,
}
ALL_CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']
FEMALE_CHROMS = [str(i) for i in range(1, 23)] + ['X']
MALE_CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']

SAMPLE = {
    'id': 'SAMPLE_001',
    'sex': 'female',
    'events': [
        {'chr':'21','type':'trisomy','cn':3,'iscn':'+21',
         'start':1,'stop':CHROM_SIZES['21'],'color':'#FC8181'},
        {'chr':'X','type':'monosomy','cn':1,'iscn':'-X',
         'start':1,'stop':CHROM_SIZES['X'],'color':'#90CDF4'},
        {'chr':'5','type':'partial_loss','cn':1,'iscn':'del(5p)',
         'start':1,'stop':45_368_512,'color':'#90CDF4'},
        {'chr':'17','type':'partial_gain','cn':3,'iscn':'dup(17q)',
         'start':43_044_000,'stop':CHROM_SIZES['17'],'color':'#FBD38D'},
    ],
    'genes': [
        {'chr':'21','name':'DYRK1A','start':37_700_000,'stop':37_865_335,'color':'#B794F4'},
        {'chr':'21','name':'RUNX1','start':34_787_801,'stop':36_004_954,'color':'#B794F4'},
        {'chr':'17','name':'TP53','start':7_661_779,'stop':7_687_538,'color':'#68D391'},
        {'chr':'17','name':'BRCA1','start':43_044_295,'stop':43_125_483,'color':'#68D391'},
        {'chr':'17','name':'ERBB2','start':39_687_914,'stop':39_730_426,'color':'#68D391'},
        {'chr':'5','name':'CTNND2','start':11_080_000,'stop':11_690_000,'color':'#68D391'},
        {'chr':'X','name':'MECP2','start':154_021_573,'stop':154_137_103,'color':'#68D391'},
        {'chr':'X','name':'DMD','start':31_119_221,'stop':33_339_609,'color':'#68D391'},
        {'chr':'13','name':'BRCA2','start':32_315_086,'stop':32_400_266,'color':'#68D391'},
        {'chr':'7','name':'CFTR','start':117_480_025,'stop':117_668_665,'color':'#68D391'},
    ],
}
DISPLAY_CHROMS = FEMALE_CHROMS if SAMPLE['sex'] == 'female' else MALE_CHROMS


# =============================================================================
# Helpers
# =============================================================================
def parse_bp(value, default=0):
    if value is None or value == '':
        return default
    if isinstance(value, (int, float, np.integer, np.floating)):
        return int(value)
    try:
        return int(float(str(value).replace(',', '').replace('bp', '').strip()))
    except (TypeError, ValueError):
        return default


def build_iscn(sample):
    total = 46
    parts = []
    for ev in sample['events']:
        if ev['type'] == 'trisomy':
            total += 1
        elif ev['type'] == 'monosomy':
            total -= 1
        parts.append(ev['iscn'])
    sex = 'XX' if sample['sex'] == 'female' else 'XY'
    return f"{total},{sex}" + (',' + ','.join(parts) if parts else '')


ISCN = build_iscn(SAMPLE)


def events_for_chrom(chrom):
    return [x for x in SAMPLE['events'] if x['chr'] == str(chrom)]


def genes_for_chrom(chrom):
    return [x for x in SAMPLE['genes'] if x['chr'] == str(chrom)]


def make_raw_annots():
    """Stable rawAnnots object shared by both isolated Dash apps."""
    by_chr = {c: [] for c in ALL_CHROMS}

    for ev in SAMPLE['events']:
        by_chr[ev['chr']].append([
            ev['iscn'],
            ev['start'],
            max(1, ev['stop'] - ev['start']),
            0,
            ev['color'],
        ])

    for gene in SAMPLE['genes']:
        by_chr[gene['chr']].append([
            gene['name'],
            gene['start'],
            max(50_000, gene['stop'] - gene['start']),
            1,
            gene['color'],
        ])

    return {
        'keys': ['name', 'start', 'length', 'trackIndex', 'color'],
        'annots': [
            {'chr': c, 'annots': by_chr[c]}
            for c in ALL_CHROMS
            if by_chr[c]
        ],
    }


RAW_ANNOTS = make_raw_annots()

ANNOTATION_TRACKS = [
    {'id':'cnv',  'displayName':'CNV',  'color':'#E53E3E', 'shape':'rectangle'},
    {'id':'gene', 'displayName':'Gene', 'color':'#6B46C1', 'shape':'circle'},
]

LEGEND = [{'name':'Legend','rows':[
    {'name':'Trisomy/Gain','color':'#FC8181','shape':'rectangle'},
    {'name':'Monosomy/Loss','color':'#90CDF4','shape':'rectangle'},
    {'name':'Partial Gain','color':'#FBD38D','shape':'rectangle'},
    {'name':'Gene','color':'#B794F4','shape':'circle'},
]}]


def annotation_options(chrom):
    chrom = str(chrom)
    options = []

    for ev in events_for_chrom(chrom):
        options.append({
            'label': (
                f"CNV | {ev['iscn']} | "
                f"{ev['start']/1e6:.3f}-{ev['stop']/1e6:.3f} Mb"
            ),
            'value': json.dumps({'start':ev['start'], 'end':ev['stop']}),
        })

    for gene in genes_for_chrom(chrom):
        options.append({
            'label': (
                f"Gene | {gene['name']} | "
                f"{gene['start']/1e6:.3f}-{gene['stop']/1e6:.3f} Mb"
            ),
            'value': json.dumps({'start':gene['start'], 'end':gene['stop']}),
        })

    return options


def demo_data(chrom, start_bp=None, end_bp=None):
    chrom = str(chrom)
    size = CHROM_SIZES[chrom]
    rng = np.random.default_rng(sum(ord(x) for x in chrom) + 111)

    # Full chromosome: moderate density. Narrow brush: regenerate enough demo
    # points inside the selected range so CN/BAF never visually disappear.
    region_mode = start_bp is not None and end_bp is not None
    if region_mode:
        start_bp = max(1, int(start_bp))
        end_bp = min(size, int(end_bp))
        if end_bp <= start_bp:
            end_bp = min(size, start_bp + 1_000_000)
        n = 420
        pos = np.linspace(start_bp, end_bp, n).astype(int)
    else:
        n = 1600
        pos = np.linspace(1, size, n).astype(int)
    cn = np.full(n, 2.0)
    baf = np.full(n, 0.5)

    for ev in events_for_chrom(chrom):
        mask = (pos >= ev['start']) & (pos <= ev['stop'])
        if ev['type'] in ('trisomy', 'partial_gain'):
            cn[mask] = 3.0
            baf[mask] = rng.choice([0.33, 0.67], mask.sum())
        elif ev['type'] in ('monosomy', 'partial_loss'):
            cn[mask] = 1.0
            baf[mask] = rng.choice([0.03, 0.97], mask.sum())

    return pd.DataFrame({
        'chrom': chrom,
        'pos': pos,
        'cn': np.clip(cn + rng.normal(0, 0.18, n), 0, 6),
        'baf': np.clip(baf + rng.normal(0, 0.035, n), 0, 1),
    })


def region_fig(df, chrom, start_bp, end_bp):
    chrom = str(chrom)
    size = CHROM_SIZES[chrom]
    work = df.copy()

    if 'chrom' in work.columns:
        work['chrom'] = work['chrom'].astype(str).str.replace('chr', '', regex=False)
        work = work[work['chrom'] == chrom]

    sub = work[(work['pos'] >= start_bp) & (work['pos'] <= end_bp)].copy()

    span_mb = max((end_bp - start_bp) / 1e6, 0.001)
    marker_size = 6 if span_mb <= 2 else (5 if span_mb <= 20 else 3.5)

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.5, 0.5],
        vertical_spacing=0.10,
        subplot_titles=[f'Copy Number - chr{chrom}', f'BAF - chr{chrom}'],
    )

    fig.add_trace(
        go.Scatter(
            x=sub['pos']/1e6,
            y=sub['cn'],
            mode='markers',
            marker={'size':marker_size, 'opacity':0.82, 'color':'#4A5568'},
            hovertemplate='%{x:.3f} Mb<br>CN=%{y:.3f}<extra></extra>',
        ),
        row=1,
        col=1,
    )

    baf_col = 'baf' if 'baf' in sub.columns else ('vaf' if 'vaf' in sub.columns else None)
    if baf_col:
        fig.add_trace(
            go.Scatter(
                x=sub['pos']/1e6,
                y=sub[baf_col],
                mode='markers',
                marker={'size':marker_size, 'opacity':0.82, 'color':'#718096'},
                hovertemplate='%{x:.3f} Mb<br>BAF=%{y:.3f}<extra></extra>',
            ),
            row=2,
            col=1,
        )

    # Keep expected CN/BAF levels visible even when the brush becomes narrow.
    if not sub.empty:
        cn_med = float(sub['cn'].median())
        fig.add_hline(y=cn_med, line={'color':'#3182CE','width':1.4,'dash':'dot'}, row=1, col=1)
        fig.add_annotation(
            x=end_bp/1e6, y=cn_med, text=f'CN~{cn_med:.2f}',
            xanchor='right', yanchor='bottom', showarrow=False,
            font={'size':10,'color':'#2B6CB0'}, row=1, col=1,
        )
        if baf_col:
            baf_med = float(sub[baf_col].median())
            fig.add_hline(y=baf_med, line={'color':'#718096','width':1,'dash':'dot'}, row=2, col=1)
    else:
        fig.add_annotation(
            text='선택 구간에 CN/BAF 데이터가 없습니다',
            x=0.5, y=0.5, xref='paper', yref='paper', showarrow=False,
            font={'size':11,'color':'#A0AEC0'},
        )

    fig.update_yaxes(title_text='CN', range=[0, 5], row=1, col=1)
    fig.update_yaxes(title_text='BAF', range=[0, 1], row=2, col=1)
    fig.update_xaxes(range=[start_bp/1e6, end_bp/1e6], row=1, col=1)
    fig.update_xaxes(title_text='Position (Mb)', range=[start_bp/1e6, end_bp/1e6], row=2, col=1)
    fig.update_layout(
        height=340,
        margin={'l':52, 'r':18, 't':42, 'b':34},
        showlegend=False,
        paper_bgcolor='white',
        plot_bgcolor='white',
        font={'family':'Inter,Arial,sans-serif', 'size':10, 'color':'#2D3748'},
    )
    return fig


# =============================================================================
# Shared UI styles
# =============================================================================
PAGE = {
    'fontFamily':'Inter,Arial,sans-serif',
    'background':'#EDF2F7',
    'color':'#2D3748',
    'padding':'8px',
    'margin':'0',
    'fontSize':'12px',
    'boxSizing':'border-box',
}

CARD = {
    'background':'white',
    'border':'1px solid #D9E2EC',
    'borderRadius':'7px',
    'marginBottom':'8px',
    'overflow':'visible',
}

DASHBAR = {
    'height':'30px',
    'boxSizing':'border-box',
    'padding':'5px 10px',
    'display':'flex',
    'alignItems':'center',
    'gap':'8px',
    'borderBottom':'1px solid #E2E8F0',
    'background':'#F7FAFC',
    'fontSize':'10px',
    'fontWeight':'700',
    'letterSpacing':'.08em',
    'color':'#718096',
}


def dashbar(title, right=None):
    children = [
        html.Span('● ● ●', style={
            'fontSize':'7px', 'color':'#CBD5E0', 'letterSpacing':'2px'
        }),
        html.Span(title),
    ]
    if right is not None:
        children.append(html.Div(
            right,
            style={
                'marginLeft':'auto',
                'display':'flex',
                'alignItems':'center',
                'gap':'8px',
            },
        ))
    return html.Div(children, style=DASHBAR)


def card(title, body, right=None, body_style=None):
    style = {'padding':'10px 12px', 'boxSizing':'border-box'}
    if body_style:
        style.update(body_style)
    return html.Div(
        [dashbar(title, right), html.Div(body, style=style)],
        style=CARD,
    )


# =============================================================================
# Dash App 1 - Overview
# =============================================================================
overview_dash = Dash(
    __name__ + '_overview',
    requests_pathname_prefix='/overview/',
    routes_pathname_prefix='/',
    suppress_callback_exceptions=True,
    title='Karyotype Overview',
)

INITIAL_CHROM = '5'

# The Store is set from injected JavaScript when an SVG chromosome is clicked.
overview_dash.layout = html.Div([
    dcc.Store(id='overview-selected-chrom', data=INITIAL_CHROM),

    card(
        'SAMPLE INFO',
        html.Div([
            html.Div([
                html.Div('Sample', style={'fontSize':'9px','color':'#A0AEC0'}),
                html.Div(SAMPLE['id'], style={'fontWeight':'700'}),
            ]),
            html.Div([
                html.Div('Sex', style={'fontSize':'9px','color':'#A0AEC0'}),
                html.Div('Female' if SAMPLE['sex']=='female' else 'Male', style={'fontWeight':'700'}),
            ]),
            html.Div([
                html.Div('ISCN', style={'fontSize':'9px','color':'#A0AEC0'}),
                html.Div(ISCN, style={'fontWeight':'700','fontFamily':'monospace'}),
            ]),
        ], style={
            'display':'grid',
            'gridTemplateColumns':'repeat(3,minmax(0,1fr))',
            'gap':'12px',
        }),
        body_style={'padding':'7px 10px'},
    ),

    card(
        'KARYOTYPE OVERVIEW',
        [
            html.Div(
                dashbio.Ideogram(
                    id='overview-ideo',
                    organism='human',
                    assembly='GRCh38',
                    orientation='vertical',
                    chromosomes=DISPLAY_CHROMS,
                    rows=1,
                    chrHeight=300,
                    chrWidth=15,
                    chrMargin=15,
                    rotatable=False,
                    showBandLabels=True,
                    showChromosomeLabels=True,
                    resolution=550,
                    annotations=RAW_ANNOTS,
                    annotationsLayout='tracks',
                    annotationTracks=ANNOTATION_TRACKS,
                    annotationHeight=8,
                    barWidth=3,
                    histogramScaling='relative',
                    showAnnotTooltip=True,
                    legend=LEGEND,
                ),
                id='overview-ideo-wrap',
                style={
                    'width':'100%',
                    'height':'330px',
                    'overflowX':'auto',
                    'overflowY':'visible',
                    'paddingTop':'2px',
                    'paddingBottom':'10px',
                    'boxSizing':'border-box',
                },
            ),
            html.Div([
                html.Span('Chromosome을 직접 클릭하세요. ', style={'color':'#718096'}),
                html.Span(id='overview-status', style={
                    'fontFamily':'monospace',
                    'color':'#3182CE',
                    'fontWeight':'600',
                }),
            ], style={'fontSize':'10px','paddingTop':'2px'}),
        ],
        right=html.Span('click chromosome → detail', style={
            'fontSize':'10px', 'fontWeight':'400', 'letterSpacing':'0', 'color':'#A0AEC0'
        }),
        body_style={'height':'365px','paddingTop':'5px','paddingBottom':'4px'},
    ),
], style=PAGE)


@overview_dash.callback(
    Output('overview-status','children'),
    Input('overview-selected-chrom','data'),
)
def overview_status(chrom):
    return f'chr{chrom} selected'


# Inject chromosome click bridge into Overview document.
# Dash-Bio itself exposes `rotated` only as a boolean, not the chromosome name.
OVERVIEW_CLICK_JS = r"""
<script>
(function () {
  const allowed = new Set([
    '1','2','3','4','5','6','7','8','9','10','11','12',
    '13','14','15','16','17','18','19','20','21','22','X','Y'
  ]);

  function normaliseChrom(value) {
    if (!value) return null;
    const s = String(value).trim();

    const patterns = [
      /(?:^|[^A-Za-z0-9])chr(?:omosome)?[-_: ]*([0-9]{1,2}|X|Y)(?:[^A-Za-z0-9]|$)/i,
      /(?:^|[^A-Za-z0-9])([0-9]{1,2}|X|Y)[-_]9606(?:[^A-Za-z0-9]|$)/i,
      /^([0-9]{1,2}|X|Y)$/i
    ];

    for (const re of patterns) {
      const m = s.match(re);
      if (m) {
        const c = m[1].toUpperCase();
        if (allowed.has(c)) return c;
      }
    }
    return null;
  }

  function chromFromNode(target) {
    const root = document.getElementById('overview-ideo');
    let node = target;

    while (node && node !== root && node !== document.body) {
      const candidates = [];
      if (node.id) candidates.push(node.id);
      if (node.getAttribute) {
        candidates.push(node.getAttribute('class'));
        candidates.push(node.getAttribute('data-chr'));
        candidates.push(node.getAttribute('data-chromosome'));
        candidates.push(node.getAttribute('aria-label'));
      }

      for (const value of candidates) {
        const c = normaliseChrom(value);
        if (c) return c;
      }
      node = node.parentElement;
    }

    // Ideogram chromosome groups also contain the chromosome label text.
    const g = target && target.closest ? target.closest('g') : null;
    if (g) {
      const texts = Array.from(g.querySelectorAll('text'));
      for (const t of texts) {
        const c = normaliseChrom(t.textContent);
        if (c) return c;
      }
    }

    // Direct click on chromosome label text.
    if (target && target.tagName && target.tagName.toLowerCase() === 'text') {
      return normaliseChrom(target.textContent);
    }

    return null;
  }

  function publish(chrom) {
    if (!chrom) return;

    // Update Overview Dash state without a server round trip.
    if (window.dash_clientside && window.dash_clientside.set_props) {
      window.dash_clientside.set_props('overview-selected-chrom', {data: chrom});
    }

    // Send to FastAPI host page; host forwards it to Detail iframe.
    window.parent.postMessage({
      type: 'karyotype-chromosome-selected',
      chrom: chrom
    }, window.location.origin);
  }

  function install() {
    const root = document.getElementById('overview-ideo');
    if (!root || root.dataset.chromClickInstalled === '1') return false;

    root.dataset.chromClickInstalled = '1';
    root.style.cursor = 'pointer';

    root.addEventListener('click', function (event) {
      const chrom = chromFromNode(event.target);
      if (chrom) publish(chrom);
    }, true);

    return true;
  }

  let tries = 0;
  const timer = setInterval(function () {
    tries += 1;
    if (install() || tries > 80) clearInterval(timer);
  }, 100);
})();
</script>
"""

overview_dash.index_string = overview_dash.index_string.replace(
    '</body>', OVERVIEW_CLICK_JS + '\n</body>'
)


# =============================================================================
# Dash App 2 - Detail
# =============================================================================
detail_dash = Dash(
    __name__ + '_detail',
    requests_pathname_prefix='/detail/',
    routes_pathname_prefix='/',
    suppress_callback_exceptions=True,
    title='Chromosome Detail',
)

INITIAL_SIZE = CHROM_SIZES[INITIAL_CHROM]


detail_dash.layout = html.Div([
    # Updated directly by postMessage bridge; no Interval / no session polling.
    dcc.Store(id='detail-selected-chrom', data=INITIAL_CHROM),
    dcc.Store(
        id='brush-region',
        data={'chrom':INITIAL_CHROM, 'start':1, 'end':INITIAL_SIZE},
    ),
    dcc.Store(id='uploaded-data', data=None),

    card(
        'CHROMOSOME DETAIL',
        [
            html.Div([
                html.Div([
                    html.Div('Annotation / Brush selector', style={
                        'fontSize':'10px',
                        'fontWeight':'600',
                        'color':'#4A5568',
                        'marginBottom':'3px',
                    }),
                    dcc.Dropdown(
                        id='annotation-select',
                        options=annotation_options(INITIAL_CHROM),
                        placeholder='CNV / Gene annotation 선택',
                        clearable=True,
                        style={'fontSize':'11px'},
                    ),
                ], style={'flex':'1', 'minWidth':'280px'}),

                html.Div('Brush는 하나만 표시됩니다', style={'fontSize':'10px','color':'#A0AEC0','paddingBottom':'7px'}),
            ], style={
                'display':'flex',
                'gap':'16px',
                'alignItems':'flex-end',
                'flexWrap':'wrap',
                'marginBottom':'6px',
            }),

            # Exactly ONE Ideogram in Detail document.
            html.Div(
                dashbio.Ideogram(
                    id='detail-ideo',
                    organism='human',
                    assembly='GRCh38',
                    orientation='horizontal',
                    chromosomes=[INITIAL_CHROM],
                    rows=1,
                    chrHeight=300,
                    chrWidth=16,
                    chrMargin=18,
                    rotatable=False,
                    showBandLabels=True,
                    showChromosomeLabels=True,
                    resolution=850,
                    annotations=None,
                    showAnnotTooltip=False,
                    brush=f'chr{INITIAL_CHROM}:1-{INITIAL_SIZE}',
                ),
                style={
                    'width':'100%',
                    'height':'145px',
                    'overflowX':'auto',
                    'overflowY':'visible',
                    'paddingTop':'4px',
                    'paddingBottom':'8px',
                    'boxSizing':'border-box',
                },
            ),

            html.Div([
                html.Span(id='brush-chip', style={
                    'fontFamily':'monospace',
                    'fontSize':'10px',
                    'color':'#3182CE',
                }),
                html.Span('  |  Brush handle을 드래그하거나 annotation을 선택하세요.', style={
                    'fontSize':'10px',
                    'color':'#A0AEC0',
                }),
            ]),
        ],
        right=html.Span(id='detail-chrom-badge', children='chr5', style={
            'fontFamily':'monospace',
            'fontWeight':'600',
            'letterSpacing':'0',
        }),
        body_style={'height':'230px', 'paddingBottom':'6px'},
    ),

    card(
        'CN / BAF REGION DETAIL',
        [
            html.Div([
                dcc.Upload(
                    id='upload-cnv',
                    children=html.Button('TSV/CSV 업로드', style={
                        'fontSize':'10px', 'padding':'4px 8px', 'cursor':'pointer'
                    }),
                    multiple=False,
                ),
                html.Button('Demo data', id='demo-btn', n_clicks=0, style={
                    'fontSize':'10px', 'padding':'4px 8px', 'cursor':'pointer'
                }),
                html.Span(id='upload-status', style={
                    'fontSize':'10px', 'color':'#718096'
                }),
            ], style={
                'display':'flex',
                'gap':'6px',
                'alignItems':'center',
                'flexWrap':'wrap',
                'marginBottom':'2px',
            }),

            dcc.Graph(
                id='region-graph',
                figure=region_fig(demo_data(INITIAL_CHROM), INITIAL_CHROM, 1, INITIAL_SIZE),
                config={
                    'displayModeBar':True,
                    'responsive':True,
                    'modeBarButtonsToRemove':['select2d','lasso2d'],
                },
                style={'height':'340px'},
            ),
        ],
        right=html.Span(id='region-title-chip', children='chr5: full', style={
            'fontFamily':'monospace',
            'fontWeight':'400',
            'letterSpacing':'0',
        }),
        body_style={'height':'375px', 'paddingBottom':'3px'},
    ),
], style=PAGE)


# Detail iframe receives chromosome messages and updates the Dash Store directly.
DETAIL_MESSAGE_JS = r"""
<script>
(function () {
  const allowed = new Set([
    '1','2','3','4','5','6','7','8','9','10','11','12',
    '13','14','15','16','17','18','19','20','21','22','X','Y'
  ]);

  function setChrom(chrom, attempt) {
    chrom = String(chrom || '').toUpperCase();
    if (!allowed.has(chrom)) return;

    if (window.dash_clientside && window.dash_clientside.set_props) {
      window.dash_clientside.set_props('detail-selected-chrom', {data: chrom});
      return;
    }

    if ((attempt || 0) < 50) {
      setTimeout(function () { setChrom(chrom, (attempt || 0) + 1); }, 100);
    }
  }

  window.addEventListener('message', function (event) {
    if (event.origin !== window.location.origin) return;
    const data = event.data || {};
    if (data.type !== 'karyotype-chromosome-selected') return;
    setChrom(data.chrom, 0);
  });
})();
</script>
"""

detail_dash.index_string = detail_dash.index_string.replace(
    '</body>', DETAIL_MESSAGE_JS + '\n</body>'
)


# =============================================================================
# Detail callbacks
# =============================================================================
@detail_dash.callback(
    Output('detail-ideo','chromosomes'),
    Output('detail-ideo','brush'),
    Output('annotation-select','options'),
    Output('annotation-select','value'),
    Output('detail-chrom-badge','children'),
    Output('brush-region','data'),
    Output('brush-chip','children'),
    Input('detail-selected-chrom','data'),
)
def detail_change_chrom(chrom):
    chrom = str(chrom or INITIAL_CHROM)
    if chrom not in CHROM_SIZES:
        chrom = INITIAL_CHROM

    size = CHROM_SIZES[chrom]
    brush = f'chr{chrom}:1-{size}'
    region = {'chrom':chrom, 'start':1, 'end':size}
    chip = f'chr{chrom}: 0.000 - {size/1e6:.3f} Mb'

    # Keep annotations / orientation / sizing fixed; change only chromosome + brush.
    return (
        [chrom],
        brush,
        annotation_options(chrom),
        None,
        f'chr{chrom}',
        region,
        chip,
    )


@detail_dash.callback(
    Output('detail-ideo','brush', allow_duplicate=True),
    Output('brush-region','data', allow_duplicate=True),
    Output('brush-chip','children', allow_duplicate=True),
    Input('annotation-select','value'),
    State('detail-selected-chrom','data'),
    prevent_initial_call=True,
)
def annotation_to_brush(value, chrom):
    if not value or not chrom:
        return no_update, no_update, no_update

    try:
        item = json.loads(value)
        chrom = str(chrom)
        start_bp = max(1, parse_bp(item.get('start'), 1))
        end_bp = min(
            CHROM_SIZES[chrom],
            parse_bp(item.get('end'), CHROM_SIZES[chrom]),
        )
    except Exception:
        return no_update, no_update, no_update

    brush = f'chr{chrom}:{start_bp}-{end_bp}'
    region = {'chrom':chrom, 'start':start_bp, 'end':end_bp}
    chip = f'chr{chrom}: {start_bp/1e6:.3f} - {end_bp/1e6:.3f} Mb'
    return brush, region, chip


@detail_dash.callback(
    Output('brush-region','data', allow_duplicate=True),
    Output('brush-chip','children', allow_duplicate=True),
    Input('detail-ideo','brushData'),
    State('detail-selected-chrom','data'),
    prevent_initial_call=True,
)
def brush_to_region(data, chrom):
    if not data or not chrom:
        return no_update, no_update

    chrom = str(chrom)
    start_bp = parse_bp(data.get('start'), 1)
    end_bp = parse_bp(data.get('end'), CHROM_SIZES[chrom])

    if end_bp < start_bp:
        start_bp, end_bp = end_bp, start_bp

    start_bp = max(1, start_bp)
    end_bp = min(CHROM_SIZES[chrom], end_bp)

    if end_bp <= start_bp:
        return no_update, no_update

    region = {'chrom':chrom, 'start':start_bp, 'end':end_bp}
    chip = f'chr{chrom}: {start_bp/1e6:.3f} - {end_bp/1e6:.3f} Mb'
    return region, chip


@detail_dash.callback(
    Output('uploaded-data','data'),
    Output('upload-status','children'),
    Input('upload-cnv','contents'),
    State('upload-cnv','filename'),
    prevent_initial_call=True,
)
def upload_data(contents, filename):
    if not contents:
        return no_update, no_update

    try:
        _, payload = contents.split(',', 1)
        decoded = base64.b64decode(payload).decode('utf-8')
        sep = '\t' if str(filename or '').lower().endswith('.tsv') else ','
        df = pd.read_csv(io.StringIO(decoded), sep=sep)
        df.columns = [str(c).strip().lower() for c in df.columns]

        rename = {}
        for c in ['position','start','chromstart','chrom_start']:
            if c in df.columns and 'pos' not in df.columns:
                rename[c] = 'pos'
                break
        for c in ['copy_number','copynumber','ratio','log2']:
            if c in df.columns and 'cn' not in df.columns:
                rename[c] = 'cn'
                break
        for c in ['vaf','b_allele_freq','allele_freq']:
            if c in df.columns and 'baf' not in df.columns:
                rename[c] = 'baf'
                break

        df = df.rename(columns=rename)
        required = [c for c in ['pos','cn'] if c not in df.columns]
        if required:
            return None, f'필수 컬럼 없음: {required}'

        keep = [c for c in ['chrom','pos','cn','baf'] if c in df.columns]
        df = df[keep].dropna(subset=['pos','cn'])

        return df.to_json(orient='split'), f'{filename} / {len(df):,} rows'
    except Exception as exc:
        return None, f'업로드 오류: {exc}'


@detail_dash.callback(
    Output('region-graph','figure'),
    Output('region-title-chip','children'),
    Input('brush-region','data'),
    Input('demo-btn','n_clicks'),
    Input('uploaded-data','data'),
)
def update_region_graph(region, _demo_clicks, uploaded_json):
    if not region:
        return no_update, no_update

    chrom = str(region['chrom'])
    start_bp = parse_bp(region['start'], 1)
    end_bp = parse_bp(region['end'], CHROM_SIZES[chrom])

    use_demo = ctx.triggered_id == 'demo-btn' or not uploaded_json
    if use_demo:
        df = demo_data(chrom, start_bp, end_bp)
    else:
        try:
            df = pd.read_json(io.StringIO(uploaded_json), orient='split')
        except Exception:
            df = demo_data(chrom, start_bp, end_bp)

    title = f'chr{chrom}: {start_bp/1e6:.3f}-{end_bp/1e6:.3f} Mb'
    return region_fig(df, chrom, start_bp, end_bp), title


# =============================================================================
# FastAPI host
# =============================================================================
fastapi_app = FastAPI(title='Karyotype Viewer')

HOST_HTML = r"""
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Karyotype Viewer</title>
<style>
  * { box-sizing: border-box; }
  html, body {
    margin: 0;
    padding: 0;
    background: #edf2f7;
    font-family: Arial, sans-serif;
  }
  .shell { padding: 0; }
  iframe {
    display: block;
    width: 100%;
    border: 0;
    margin: 0;
    padding: 0;
    background: #edf2f7;
  }
  /* Fixed sizes: no scrollHeight polling, so no cumulative line-by-line drift. */
  #overview-frame { height: 470px; }
  #detail-frame   { height: 675px; }
</style>
</head>
<body>
<div class="shell">
  <iframe
    id="overview-frame"
    src="/overview/"
    title="Karyotype overview"
    scrolling="no"></iframe>

  <iframe
    id="detail-frame"
    src="/detail/"
    title="Chromosome detail"
    scrolling="no"></iframe>
</div>

<script>
(function () {
  const detailFrame = document.getElementById('detail-frame');
  let lastChrom = '5';

  function forward(chrom) {
    lastChrom = String(chrom || lastChrom);
    if (detailFrame && detailFrame.contentWindow) {
      detailFrame.contentWindow.postMessage({
        type: 'karyotype-chromosome-selected',
        chrom: lastChrom
      }, window.location.origin);
    }
  }

  window.addEventListener('message', function (event) {
    if (event.origin !== window.location.origin) return;
    const data = event.data || {};
    if (data.type !== 'karyotype-chromosome-selected') return;
    forward(data.chrom);
  });

  // Re-send the current selection only when the detail iframe initially loads.
  // This does NOT reload the iframe and does not resize anything.
  detailFrame.addEventListener('load', function () {
    setTimeout(function () { forward(lastChrom); }, 150);
  });
})();
</script>
</body>
</html>
"""


@fastapi_app.get('/', response_class=HTMLResponse)
async def root():
    return HTMLResponse(HOST_HTML)


# Dash is WSGI; FastAPI/Starlette mounts each isolated Dash server under its path.
fastapi_app.mount('/overview', WSGIMiddleware(overview_dash.server))
fastapi_app.mount('/detail', WSGIMiddleware(detail_dash.server))

# Conventional ASGI variable for: uvicorn module:app
app = fastapi_app


if __name__ == '__main__':
    import uvicorn
    uvicorn.run(app, host='0.0.0.0', port=8050, reload=False)