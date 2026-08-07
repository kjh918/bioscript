"""
Karyotype Viewer — FastAPI + 2 isolated Dash apps

Architecture
------------
FastAPI host :8050
  /           → host HTML (두 iframe 포함)
  /overview/  → Dash App 1  전체 karyogram (vertical, click → postMessage)
  /detail/    → Dash App 2  단일 염색체 (horizontal + brush + CN/BAF)

Key design
----------
- 두 Dash 앱이 완전히 분리된 DOM → Ideogram.js 인스턴스 충돌 없음
- display:none 없이 visibility:hidden 사용 → brush 크기 0 문제 없음
- postMessage 로 Overview → Detail 염색체 전달 (서버 round-trip 없음)
- brush / annotation-select 모두 Detail 앱에서 독립 처리

Install
-------
pip install fastapi uvicorn dash dash-bio plotly pandas numpy

Run
---
python karyotype_app.py
http://localhost:8050
"""

import io, json, base64, re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from starlette.middleware.wsgi import WSGIMiddleware

from dash import Dash, html, dcc, Input, Output, State, no_update, ctx
import dash_bio as dashbio


# ══════════════════════════════════════════════════════════════════════
# 1. 참조 데이터
# ══════════════════════════════════════════════════════════════════════
CHROM_SIZES = {
    '1':248956422,'2':242193529,'3':198295559,'4':190214555,
    '5':181538259,'6':170805979,'7':159345973,'8':145138636,
    '9':138394717,'10':133797422,'11':135086622,'12':133275309,
    '13':114364328,'14':107043718,'15':101991189,'16':90338345,
    '17':83257441,'18':80373285,'19':58617616,'20':64444167,
    '21':46709983,'22':50818468,'X':156040895,'Y':57227415,
}
ALL_CHROMS    = [str(i) for i in range(1, 23)] + ['X', 'Y']
FEMALE_CHROMS = [str(i) for i in range(1, 23)] + ['X']
MALE_CHROMS   = [str(i) for i in range(1, 23)] + ['X', 'Y']

CENT_MB = {
    '1':123.4,'2':93.9,'3':90.9,'4':50.2,'5':48.8,
    '6':61.0,'7':59.9,'8':45.2,'9':43.0,'10':39.8,
    '11':53.4,'12':35.5,'13':17.7,'14':17.2,'15':19.0,
    '16':36.8,'17':25.1,'18':18.5,'19':26.2,'20':28.1,
    '21':12.0,'22':15.0,'X':61.0,'Y':10.4,
}

# ══════════════════════════════════════════════════════════════════════
# 2. 샘플 데이터  ← 파이프라인 결과로 교체
# ══════════════════════════════════════════════════════════════════════
SAMPLE = {
    'id':  'SAMPLE_001',
    'sex': 'female',
    'events': [
        {'chr':'21','type':'trisomy',      'cn':3,'iscn':'+21',
         'start':1,          'stop':46709983,  'color':'#FC8181'},
        {'chr':'X', 'type':'monosomy',     'cn':1,'iscn':'-X',
         'start':1,          'stop':156040895, 'color':'#90CDF4'},
        {'chr':'5', 'type':'partial_loss', 'cn':1,'iscn':'del(5p)',
         'start':1,          'stop':45368512,  'color':'#90CDF4'},
        {'chr':'17','type':'partial_gain', 'cn':3,'iscn':'dup(17q)',
         'start':43044000,   'stop':83257441,  'color':'#FBD38D'},
    ],
    'genes': [
        {'chr':'21','name':'DYRK1A','start':37700000, 'stop':37865335, 'color':'#B794F4'},
        {'chr':'21','name':'RUNX1', 'start':34787801, 'stop':36004954, 'color':'#B794F4'},
        {'chr':'17','name':'TP53',  'start':7661779,  'stop':7687538,  'color':'#68D391'},
        {'chr':'17','name':'BRCA1', 'start':43044295, 'stop':43125483, 'color':'#68D391'},
        {'chr':'17','name':'ERBB2', 'start':39687914, 'stop':39730426, 'color':'#68D391'},
        {'chr':'5', 'name':'CTNND2','start':11080000, 'stop':11690000, 'color':'#68D391'},
        {'chr':'X', 'name':'MECP2', 'start':154021573,'stop':154137103,'color':'#68D391'},
        {'chr':'X', 'name':'DMD',   'start':31119221, 'stop':33339609, 'color':'#68D391'},
        {'chr':'13','name':'BRCA2', 'start':32315086, 'stop':32400266, 'color':'#68D391'},
        {'chr':'7', 'name':'CFTR',  'start':117480025,'stop':117668665,'color':'#68D391'},
    ],
}

DISPLAY_CHROMS = FEMALE_CHROMS if SAMPLE['sex'] == 'female' else MALE_CHROMS
INITIAL_CHROM  = '21'


# ══════════════════════════════════════════════════════════════════════
# 3. 공통 유틸
# ══════════════════════════════════════════════════════════════════════
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
        if ev['type'] == 'trisomy':    total += 1
        elif ev['type'] == 'monosomy': total -= 1
        parts.append(ev['iscn'])
    sex = 'XX' if sample['sex'] == 'female' else 'XY'
    return f"{total},{sex}" + (',' + ','.join(parts) if parts else '')


ISCN = build_iscn(SAMPLE)


def events_for(chrom):
    return [e for e in SAMPLE['events'] if e['chr'] == str(chrom)]


def genes_for(chrom):
    return [g for g in SAMPLE['genes'] if g['chr'] == str(chrom)]


def _rgba(h, a):
    c = h.lstrip('#')
    r, g, b = int(c[0:2],16), int(c[2:4],16), int(c[4:6],16)
    return f'rgba({r},{g},{b},{a})'


# ══════════════════════════════════════════════════════════════════════
# 4. Ideogram rawAnnots
# ══════════════════════════════════════════════════════════════════════
def make_raw_annots():
    by_chr = {c: [] for c in ALL_CHROMS}
    for ev in SAMPLE['events']:
        by_chr[ev['chr']].append([
            ev['iscn'], ev['start'],
            max(1, ev['stop'] - ev['start']),
            0, ev['color'],
        ])
    for g in SAMPLE['genes']:
        by_chr[g['chr']].append([
            g['name'], g['start'],
            max(50_000, g['stop'] - g['start']),
            1, g['color'],
        ])
    return {
        'keys': ['name', 'start', 'length', 'trackIndex', 'color'],
        'annots': [{'chr': c, 'annots': by_chr[c]} for c in ALL_CHROMS if by_chr[c]],
    }


RAW_ANNOTS = make_raw_annots()

LEGEND = [{'name':'Legend','rows':[
    {'name':'Trisomy/Gain',  'color':'#FC8181','shape':'rectangle'},
    {'name':'Monosomy/Loss', 'color':'#90CDF4','shape':'rectangle'},
    {'name':'Partial Gain',  'color':'#FBD38D','shape':'rectangle'},
    {'name':'Gene',          'color':'#B794F4','shape':'circle'},
]}]


def annotation_options(chrom):
    opts = []
    for ev in events_for(chrom):
        opts.append({
            'label': f"CNV | {ev['iscn']} | {ev['start']/1e6:.3f}–{ev['stop']/1e6:.3f} Mb",
            'value': json.dumps({'start': ev['start'], 'end': ev['stop']}),
        })
    for g in genes_for(chrom):
        opts.append({
            'label': f"Gene | {g['name']} | {g['start']/1e6:.3f}–{g['stop']/1e6:.3f} Mb",
            'value': json.dumps({'start': g['start'], 'end': g['stop']}),
        })
    return opts


# ══════════════════════════════════════════════════════════════════════
# 5. Demo 데이터 / Region figure
# ══════════════════════════════════════════════════════════════════════
def demo_data(chrom, start_bp=None, end_bp=None):
    chrom = str(chrom)
    size  = CHROM_SIZES[chrom]
    rng   = np.random.default_rng(sum(ord(x) for x in chrom) + 111)

    if start_bp is not None and end_bp is not None:
        start_bp = max(1, int(start_bp))
        end_bp   = min(size, int(end_bp))
        if end_bp <= start_bp:
            end_bp = min(size, start_bp + 1_000_000)
        n   = 420
        pos = np.linspace(start_bp, end_bp, n).astype(int)
    else:
        n   = 1600
        pos = np.linspace(1, size, n).astype(int)

    cn  = np.full(n, 2.0)
    baf = np.full(n, 0.5)
    for ev in events_for(chrom):
        mask = (pos >= ev['start']) & (pos <= ev['stop'])
        if ev['type'] in ('trisomy', 'partial_gain'):
            cn[mask]  = 3.0
            baf[mask] = rng.choice([0.33, 0.67], mask.sum())
        elif ev['type'] in ('monosomy', 'partial_loss'):
            cn[mask]  = 1.0
            baf[mask] = rng.choice([0.03, 0.97], mask.sum())

    return pd.DataFrame({
        'chrom': chrom,
        'pos':   pos,
        'cn':    np.clip(cn  + rng.normal(0, 0.18, n), 0, 6),
        'baf':   np.clip(baf + rng.normal(0, 0.035, n), 0, 1),
    })


def region_fig(df, chrom, start_bp, end_bp):
    chrom = str(chrom)
    size  = CHROM_SIZES[chrom]

    work = df.copy()
    if 'chrom' in work.columns:
        work['chrom'] = work['chrom'].astype(str).str.replace('chr', '', regex=False)
        work = work[work['chrom'] == chrom]

    sub = work[(work['pos'] >= start_bp) & (work['pos'] <= end_bp)].copy()

    span_mb = max((end_bp - start_bp) / 1e6, 0.001)
    msize   = 5.5 if span_mb <= 5 else (4.5 if span_mb <= 30 else 3.2)
    baf_col = next((c for c in ['baf', 'vaf'] if c in sub.columns), None)
    has_baf = baf_col is not None
    nrows   = 2 if has_baf else 1

    titles = [f'Copy Number — chr{chrom}']
    if has_baf:
        titles.append(f'BAF — chr{chrom}')

    fig = make_subplots(
        rows=nrows, cols=1, shared_xaxes=True,
        row_heights=[0.52, 0.48] if has_baf else [1.0],
        vertical_spacing=0.10,
        subplot_titles=titles,
    )

    def pt_color(p):
        for ev in events_for(chrom):
            if ev['start'] <= p <= ev['stop']:
                return ev['color']
        return '#94A3B8'

    colors = sub['pos'].apply(pt_color).tolist() if not sub.empty else []

    if sub.empty:
        fig.add_annotation(
            text='선택 구간에 데이터가 없습니다',
            x=0.5, y=0.5, xref='paper', yref='paper',
            showarrow=False, font={'size':12,'color':'#A0AEC0'},
        )
    else:
        fig.add_trace(go.Scattergl(
            x=sub['pos']/1e6, y=sub['cn'], mode='markers',
            marker=dict(size=msize, color=colors, opacity=0.85),
            hovertemplate='%{x:.3f} Mb<br>CN: %{y:.3f}<extra></extra>',
            name='CN',
        ), row=1, col=1)
        if has_baf:
            fig.add_trace(go.Scattergl(
                x=sub['pos']/1e6, y=sub[baf_col], mode='markers',
                marker=dict(size=msize, color=colors, opacity=0.75),
                hovertemplate='%{x:.3f} Mb<br>BAF: %{y:.3f}<extra></extra>',
                name='BAF',
            ), row=2, col=1)

    # 이벤트 구간 shade
    for ev in events_for(chrom):
        x0 = max(start_bp, ev['start']) / 1e6
        x1 = min(end_bp,   ev['stop'])  / 1e6
        if x1 > x0:
            h = ev['color'].lstrip('#')
            fill = f"rgba({int(h[0:2],16)},{int(h[2:4],16)},{int(h[4:6],16)},0.10)"
            for row in range(1, nrows + 1):
                fig.add_vrect(x0=x0, x1=x1, fillcolor=fill,
                              line_color=ev['color'], line_width=1.2,
                              row=row, col=1)

    # 기준선
    for v, op in [(1,0.4),(2,0.7),(3,0.4)]:
        fig.add_hline(y=v, line=dict(color='#64748B',width=1,dash='dot'), opacity=op, row=1, col=1)
    if has_baf:
        for v in [0.33, 0.5, 0.67]:
            fig.add_hline(y=v, line=dict(color='#64748B',width=1,dash='dot'), opacity=0.4, row=2, col=1)

    ax = dict(gridcolor='#E2E8F0', zeroline=False, linecolor='#E2E8F0',
              tickfont=dict(size=10, color='#64748B'))
    fig.update_xaxes(**ax, range=[start_bp/1e6, end_bp/1e6])
    fig.update_yaxes(**ax)
    fig.update_yaxes(title_text='CN', range=[0, 5.5], row=1, col=1)
    if has_baf:
        fig.update_yaxes(
            title_text='BAF', range=[-0.05, 1.05],
            tickvals=[0,0.25,0.33,0.5,0.67,0.75,1],
            ticktext=['0','.25','⅓','.5','⅔','.75','1'],
            row=2, col=1,
        )
    fig.update_xaxes(title_text='Position (Mb)', row=nrows, col=1)
    fig.update_annotations(font=dict(size=11, color='#334155'))
    fig.update_layout(
        height=380 if has_baf else 260,
        margin=dict(l=56, r=28, t=46, b=38),
        paper_bgcolor='white', plot_bgcolor='white',
        showlegend=False, hovermode='x unified',
        font=dict(family='Inter,system-ui,sans-serif', size=11),
        hoverlabel=dict(bgcolor='white', bordercolor='#E2E8F0',
                        font=dict(size=11, color='#1E293B')),
    )
    return fig


EMPTY_FIG = go.Figure().update_layout(
    paper_bgcolor='white', plot_bgcolor='white', height=300,
    margin=dict(t=20,b=20,l=20,r=20),
    annotations=[dict(
        text='← 염색체를 클릭하면 CN / BAF가 표시됩니다',
        x=0.5, y=0.5, xref='paper', yref='paper',
        showarrow=False, font=dict(size=12, color='#A0AEC0'),
    )],
)


# ══════════════════════════════════════════════════════════════════════
# 6. 공통 UI 헬퍼
# ══════════════════════════════════════════════════════════════════════
_PAGE = {
    'fontFamily': 'Inter,system-ui,sans-serif',
    'background':  '#EDF2F7',
    'color':       '#1E293B',
    'padding':     '10px 12px',
    'margin':      '0',
    'fontSize':    '13px',
    'boxSizing':   'border-box',
}
_CARD = {
    'background':   'white',
    'border':       '1px solid #E2E8F0',
    'borderRadius': '8px',
    'marginBottom': '10px',
    'boxShadow':    '0 1px 3px rgba(0,0,0,.06)',
    'overflow':     'visible',
}
_HEAD_S = {
    'padding':'8px 14px','borderBottom':'1px solid #E2E8F0',
    'background':'#F8FAFC','display':'flex','alignItems':'center','gap':'8px',
    'fontSize':'10px','fontWeight':'700','color':'#718096',
    'letterSpacing':'.1em','textTransform':'uppercase',
}


def dashbar(title, right=None):
    ch = [
        html.Span('● ● ●', style={'fontSize':'7px','color':'#CBD5E0','letterSpacing':'2px'}),
        html.Span(title),
    ]
    if right is not None:
        ch.append(html.Div(right, style={'marginLeft':'auto','display':'flex',
                                         'alignItems':'center','gap':'8px'}))
    return html.Div(ch, style=_HEAD_S)


def card(title, body, right=None, body_style=None):
    bs = {'padding':'12px 14px','boxSizing':'border-box'}
    if body_style:
        bs.update(body_style)
    return html.Div([dashbar(title, right), html.Div(body, style=bs)], style=_CARD)


def badge(text, color='#3182CE'):
    return html.Span(text, style={
        'background':_rgba(color,0.1),'color':color,
        'border':f'1px solid {_rgba(color,0.3)}',
        'borderRadius':'4px','padding':'2px 9px',
        'fontSize':'12px','fontWeight':'700','fontFamily':'monospace',
    })


def ev_card(ev):
    c = ev['color']
    return html.Div(style={
        'border':f'1px solid {c}88','borderLeft':f'4px solid {c}',
        'borderRadius':'6px','padding':'8px 12px',
        'background':_rgba(c,0.06),'minWidth':'140px',
    }, children=[
        html.Div(ev['iscn'], style={
            'fontWeight':'800','fontFamily':'monospace','fontSize':'13px','color':c,
        }),
        html.Div(f"chr{ev['chr']}  {ev['type'].replace('_',' ')}",
                 style={'fontSize':'11px','color':'#718096','marginTop':'2px'}),
        html.Div(f"CN = {ev['cn']}",
                 style={'fontSize':'11px','color':'#4A5568','marginTop':'1px'}),
    ])


# ══════════════════════════════════════════════════════════════════════
# 7. Dash App 1 — Overview
# ══════════════════════════════════════════════════════════════════════
overview_app = Dash(
    __name__ + '_overview',
    requests_pathname_prefix='/overview/',
    routes_pathname_prefix='/',
    suppress_callback_exceptions=True,
    title='Karyotype Overview',
)

_sex_lbl = '♀ Female' if SAMPLE['sex'] == 'female' else '♂ Male'
_total   = int(ISCN.split(',')[0])

overview_app.layout = html.Div([
    dcc.Store(id='ov-chrom', data=INITIAL_CHROM),

    # Sample Info
    card('SAMPLE INFO', html.Div(style={
        'display':'grid',
        'gridTemplateColumns':'repeat(auto-fill,minmax(160px,1fr))',
        'gap':'12px',
    }, children=[
        html.Div([
            html.Div('Sample ID', style={'fontSize':'9px','color':'#A0AEC0','fontWeight':'700',
                                         'textTransform':'uppercase','letterSpacing':'.06em'}),
            html.Div(SAMPLE['id'], style={'fontWeight':'700','fontSize':'14px','marginTop':'2px',
                                          'fontFamily':'monospace'}),
        ]),
        html.Div([
            html.Div('Sex', style={'fontSize':'9px','color':'#A0AEC0','fontWeight':'700',
                                    'textTransform':'uppercase','letterSpacing':'.06em'}),
            html.Div(_sex_lbl, style={
                'fontWeight':'700','fontSize':'14px','marginTop':'2px',
                'color':'#E53E3E' if SAMPLE['sex']=='female' else '#3182CE',
            }),
        ]),
        html.Div([
            html.Div('Karyotype (ISCN)', style={'fontSize':'9px','color':'#A0AEC0','fontWeight':'700',
                                                 'textTransform':'uppercase','letterSpacing':'.06em'}),
            html.Div(badge(ISCN), style={'marginTop':'5px'}),
        ]),
        html.Div([
            html.Div('Total Chromosomes', style={'fontSize':'9px','color':'#A0AEC0','fontWeight':'700',
                                                  'textTransform':'uppercase','letterSpacing':'.06em'}),
            html.Div(str(_total), style={
                'fontWeight':'800','fontSize':'22px','marginTop':'2px',
                'color':'#E53E3E' if _total != 46 else '#1E293B',
            }),
        ]),
        html.Div(style={'gridColumn':'1/-1'}, children=[
            html.Div('Chromosomal Events', style={
                'fontSize':'9px','color':'#A0AEC0','fontWeight':'700',
                'textTransform':'uppercase','letterSpacing':'.06em','marginBottom':'8px',
            }),
            html.Div([ev_card(ev) for ev in SAMPLE['events']],
                     style={'display':'flex','gap':'8px','flexWrap':'wrap'}),
        ]),
    ]), body_style={'padding':'12px 14px'}),

    # Karyogram
    card('KARYOGRAM', [
        html.Div(
            dashbio.Ideogram(
                id='karyogram',
                organism='human', assembly='GRCh38',
                orientation='vertical',
                chromosomes=DISPLAY_CHROMS,
                rows=1,
                chrHeight=300, chrWidth=15, chrMargin=15,
                rotatable=False,
                showBandLabels=True, showChromosomeLabels=True,
                resolution=550,
                annotations=RAW_ANNOTS,
                annotationsLayout='overlay',
                annotationHeight=5,
                showAnnotTooltip=True,
                legend=LEGEND,
            ),
            style={'width':'100%','overflowX':'auto','overflowY':'visible',
                   'paddingTop':'4px','paddingBottom':'8px'},
        ),
        html.Div([
            html.Span('염색체를 직접 클릭하세요 → ',
                      style={'color':'#A0AEC0'}),
            html.Span(id='ov-status',
                      style={'fontFamily':'monospace','color':'#3182CE','fontWeight':'600'}),
        ], style={'fontSize':'10px','marginTop':'4px'}),
    ],
    right=html.Span('click → detail',
                    style={'fontSize':'10px','fontWeight':'400','letterSpacing':'0','color':'#A0AEC0'}),
    body_style={'padding':'10px 14px 8px'}),

], style=_PAGE)


@overview_app.callback(
    Output('ov-status', 'children'),
    Input('ov-chrom', 'data'),
)
def ov_status(chrom):
    return f'chr{chrom} selected'


_OVERVIEW_JS = r"""
<script>
(function(){
  const ALLOWED=new Set(['1','2','3','4','5','6','7','8','9','10','11','12',
    '13','14','15','16','17','18','19','20','21','22','X','Y']);

  function norm(v){
    if(!v) return null;
    const s=String(v).trim();
    const pats=[
      /(?:^|[^A-Za-z0-9])chr(?:omosome)?[-_: ]*([0-9]{1,2}|X|Y)(?:[^A-Za-z0-9]|$)/i,
      /(?:^|[^A-Za-z0-9])([0-9]{1,2}|X|Y)[-_]9606(?:[^A-Za-z0-9]|$)/i,
      /^([0-9]{1,2}|X|Y)$/i,
    ];
    for(const re of pats){
      const m=s.match(re);
      if(m){const c=m[1].toUpperCase();if(ALLOWED.has(c))return c;}
    }
    return null;
  }

  function fromNode(target){
    const root=document.getElementById('karyogram');
    let node=target;
    while(node&&node!==root&&node!==document.body){
      for(const v of[
        node.id,
        node.getAttribute&&node.getAttribute('class'),
        node.getAttribute&&node.getAttribute('data-chr'),
        node.getAttribute&&node.getAttribute('aria-label'),
      ]){const c=norm(v);if(c)return c;}
      node=node.parentElement;
    }
    const g=target&&target.closest&&target.closest('g');
    if(g){for(const t of g.querySelectorAll('text')){const c=norm(t.textContent);if(c)return c;}}
    if(target&&target.tagName&&target.tagName.toLowerCase()==='text')return norm(target.textContent);
    return null;
  }

  function publish(chrom){
    if(!chrom)return;
    if(window.dash_clientside&&window.dash_clientside.set_props)
      window.dash_clientside.set_props('ov-chrom',{data:chrom});
    window.parent.postMessage({type:'karyotype-chrom-selected',chrom},window.location.origin);
  }

  function install(){
    const root=document.getElementById('karyogram');
    if(!root||root.dataset.installed)return false;
    root.dataset.installed='1';
    root.style.cursor='pointer';
    root.addEventListener('click',e=>{const c=fromNode(e.target);if(c)publish(c);},true);
    return true;
  }

  let tries=0;
  const t=setInterval(()=>{if(install()||++tries>80)clearInterval(t);},100);
})();
</script>"""

overview_app.index_string = overview_app.index_string.replace(
    '</body>', _OVERVIEW_JS + '\n</body>'
)


# ══════════════════════════════════════════════════════════════════════
# 8. Dash App 2 — Detail
# ══════════════════════════════════════════════════════════════════════
detail_app = Dash(
    __name__ + '_detail',
    requests_pathname_prefix='/detail/',
    routes_pathname_prefix='/',
    suppress_callback_exceptions=True,
    title='Chromosome Detail',
)

_INIT_SIZE = CHROM_SIZES[INITIAL_CHROM]

detail_app.layout = html.Div([
    dcc.Store(id='det-chrom',  data=INITIAL_CHROM),
    dcc.Store(id='det-region', data={'chrom':INITIAL_CHROM,'start':1,'end':_INIT_SIZE}),
    dcc.Store(id='det-upload', data=None),

    # Chromosome + Brush
    card('CHROMOSOME DETAIL', [

        html.Div(style={'display':'flex','gap':'12px','alignItems':'flex-end',
                        'flexWrap':'wrap','marginBottom':'8px'}, children=[
            html.Div([
                html.Div('Annotation / Brush 구간 선택',
                         style={'fontSize':'10px','fontWeight':'600','color':'#4A5568','marginBottom':'3px'}),
                dcc.Dropdown(
                    id='annot-select',
                    options=annotation_options(INITIAL_CHROM),
                    placeholder='CNV / Gene annotation 선택 → brush 자동 이동',
                    clearable=True,
                    style={'fontSize':'11px','minWidth':'340px'},
                ),
            ], style={'flex':'1'}),
            html.Div(id='det-badge', style={
                'fontFamily':'monospace','fontWeight':'700',
                'fontSize':'14px','color':'#3182CE','paddingBottom':'6px',
            }, children=f'chr{INITIAL_CHROM}'),
        ]),

        html.Div(
            dashbio.Ideogram(
                id='det-ideo',
                organism='human', assembly='GRCh38',
                orientation='horizontal',
                chromosomes=[INITIAL_CHROM],
                rows=1,
                chrHeight=110, chrWidth=1140,
                rotatable=False,
                showBandLabels=True, showChromosomeLabels=True,
                resolution=850,
                annotations=None,
                showAnnotTooltip=False,
                brush=f'chr{INITIAL_CHROM}:1-{_INIT_SIZE}',
            ),
            style={'width':'100%','overflowX':'auto','overflowY':'visible',
                   'paddingTop':'4px','paddingBottom':'6px'},
        ),

        html.Div(id='brush-chip', style={
            'fontSize':'10px','fontFamily':'monospace','color':'#3182CE','marginTop':'2px',
        }),
        html.Div('← → 핸들 드래그로 구간 선택  ·  annotation 드롭다운으로 구간 이동',
                 style={'fontSize':'10px','color':'#A0AEC0','marginTop':'2px'}),

    ], body_style={'padding':'10px 14px 10px'}),

    # CN / BAF
    card('CN / BAF REGION DETAIL', [

        html.Div(style={'display':'flex','gap':'8px','alignItems':'center',
                        'flexWrap':'wrap','marginBottom':'6px'}, children=[
            dcc.Upload(
                id='det-upload-ctrl',
                children=html.Button('📂 TSV/CSV 업로드', style={
                    'fontSize':'11px','padding':'4px 10px','cursor':'pointer',
                    'background':'#F8FAFC','border':'1px solid #CBD5E0','borderRadius':'4px',
                }),
                multiple=False,
            ),
            html.Button('Demo', id='det-demo-btn', n_clicks=0, style={
                'fontSize':'11px','padding':'4px 12px','cursor':'pointer',
                'background':'#3182CE','color':'white','border':'none','borderRadius':'4px',
            }),
            html.Span(id='det-upload-status', style={'fontSize':'11px','color':'#718096'}),
            html.Span(id='det-region-chip', style={
                'marginLeft':'auto','fontFamily':'monospace',
                'fontSize':'11px','color':'#3182CE','fontWeight':'600',
            }),
        ]),

        dcc.Graph(
            id='det-graph',
            figure=region_fig(demo_data(INITIAL_CHROM), INITIAL_CHROM, 1, _INIT_SIZE),
            config={'displayModeBar':True,'responsive':True,
                    'modeBarButtonsToRemove':['select2d','lasso2d']},
        ),

    ], body_style={'padding':'10px 14px 10px'}),

], style=_PAGE)


_DETAIL_JS = r"""
<script>
(function(){
  const ALLOWED=new Set(['1','2','3','4','5','6','7','8','9','10','11','12',
    '13','14','15','16','17','18','19','20','21','22','X','Y']);

  function set(chrom,attempt){
    chrom=String(chrom||'').toUpperCase();
    if(!ALLOWED.has(chrom))return;
    if(window.dash_clientside&&window.dash_clientside.set_props){
      window.dash_clientside.set_props('det-chrom',{data:chrom});
      return;
    }
    if((attempt||0)<50)setTimeout(()=>set(chrom,(attempt||0)+1),100);
  }

  window.addEventListener('message',e=>{
    if(e.origin!==window.location.origin)return;
    const d=e.data||{};
    if(d.type!=='karyotype-chrom-selected')return;
    set(d.chrom,0);
  });
})();
</script>"""

detail_app.index_string = detail_app.index_string.replace(
    '</body>', _DETAIL_JS + '\n</body>'
)


# ── Detail callbacks ─────────────────────────────────────────────────

@detail_app.callback(
    Output('det-ideo',    'chromosomes'),
    Output('det-ideo',    'brush'),
    Output('annot-select','options'),
    Output('annot-select','value'),
    Output('det-badge',   'children'),
    Output('det-region',  'data'),
    Output('brush-chip',  'children'),
    Input('det-chrom',    'data'),
)
def chrom_changed(chrom):
    chrom = str(chrom or INITIAL_CHROM)
    if chrom not in CHROM_SIZES:
        chrom = INITIAL_CHROM
    size   = CHROM_SIZES[chrom]
    brush  = f'chr{chrom}:1-{size}'
    region = {'chrom':chrom,'start':1,'end':size}
    chip   = f'chr{chrom}: 0.000 – {size/1e6:.3f} Mb  (전체)'
    return [chrom], brush, annotation_options(chrom), None, f'chr{chrom}', region, chip


@detail_app.callback(
    Output('det-ideo',   'brush',           allow_duplicate=True),
    Output('det-region', 'data',            allow_duplicate=True),
    Output('brush-chip', 'children',        allow_duplicate=True),
    Input('annot-select','value'),
    State('det-chrom',   'data'),
    prevent_initial_call=True,
)
def annot_to_brush(value, chrom):
    if not value or not chrom:
        return no_update, no_update, no_update
    try:
        item     = json.loads(value)
        chrom    = str(chrom)
        start_bp = max(1, parse_bp(item.get('start'), 1))
        end_bp   = min(CHROM_SIZES[chrom], parse_bp(item.get('end'), CHROM_SIZES[chrom]))
    except Exception:
        return no_update, no_update, no_update
    brush  = f'chr{chrom}:{start_bp}-{end_bp}'
    region = {'chrom':chrom,'start':start_bp,'end':end_bp}
    chip   = f'chr{chrom}: {start_bp/1e6:.3f} – {end_bp/1e6:.3f} Mb'
    return brush, region, chip


@detail_app.callback(
    Output('det-region', 'data',     allow_duplicate=True),
    Output('brush-chip', 'children', allow_duplicate=True),
    Input('det-ideo',    'brushData'),
    State('det-chrom',   'data'),
    prevent_initial_call=True,
)
def brush_moved(data, chrom):
    if not data or not chrom:
        return no_update, no_update
    chrom    = str(chrom)
    start_bp = parse_bp(data.get('start'), 1)
    end_bp   = parse_bp(data.get('end'),   CHROM_SIZES[chrom])
    if end_bp < start_bp:
        start_bp, end_bp = end_bp, start_bp
    start_bp = max(1, start_bp)
    end_bp   = min(CHROM_SIZES[chrom], end_bp)
    if end_bp <= start_bp:
        return no_update, no_update
    region = {'chrom':chrom,'start':start_bp,'end':end_bp}
    chip   = f'chr{chrom}: {start_bp/1e6:.3f} – {end_bp/1e6:.3f} Mb'
    return region, chip


@detail_app.callback(
    Output('det-upload',        'data'),
    Output('det-upload-status', 'children'),
    Input('det-upload-ctrl',    'contents'),
    State('det-upload-ctrl',    'filename'),
    prevent_initial_call=True,
)
def upload_file(contents, filename):
    if not contents:
        return no_update, no_update
    try:
        _, payload = contents.split(',', 1)
        decoded    = base64.b64decode(payload).decode('utf-8')
        sep        = '\t' if str(filename or '').lower().endswith('.tsv') else ','
        df         = pd.read_csv(io.StringIO(decoded), sep=sep)
        df.columns = [str(c).strip().lower() for c in df.columns]
        for src, dst in [
            ('position','pos'),('start','pos'),('chromstart','pos'),('chrom_start','pos'),
            ('copy_number','cn'),('copynumber','cn'),('copy_number_signal','cn'),
            ('log2','cn'),('ratio','cn'),
            ('vaf','baf'),('b_allele_freq','baf'),('bin_baf','baf'),
        ]:
            if src in df.columns and dst not in df.columns:
                df = df.rename(columns={src: dst})
        missing = [c for c in ['pos','cn'] if c not in df.columns]
        if missing:
            return None, f'❌ 컬럼 없음: {missing}'
        keep = [c for c in ['chrom','pos','cn','baf'] if c in df.columns]
        df   = df[keep].dropna(subset=['pos','cn'])
        return df.to_json(orient='split'), f'✅ {filename} ({len(df):,} rows)'
    except Exception as ex:
        return None, f'❌ {ex}'


@detail_app.callback(
    Output('det-graph',       'figure'),
    Output('det-region-chip', 'children'),
    Input('det-region',    'data'),
    Input('det-demo-btn',  'n_clicks'),
    Input('det-upload',    'data'),
)
def render_graph(region, _demo, uploaded_json):
    if not region:
        return EMPTY_FIG, '—'
    chrom    = str(region['chrom'])
    start_bp = parse_bp(region['start'], 1)
    end_bp   = parse_bp(region['end'],   CHROM_SIZES[chrom])

    use_demo = (ctx.triggered_id == 'det-demo-btn') or not uploaded_json
    if use_demo:
        df = demo_data(chrom, start_bp, end_bp)
    else:
        try:
            df = pd.read_json(io.StringIO(uploaded_json), orient='split')
        except Exception:
            df = demo_data(chrom, start_bp, end_bp)

    chip = f'chr{chrom}: {start_bp/1e6:.3f}–{end_bp/1e6:.3f} Mb'
    return region_fig(df, chrom, start_bp, end_bp), chip


# ══════════════════════════════════════════════════════════════════════
# 9. FastAPI host
# ══════════════════════════════════════════════════════════════════════
fastapi_app = FastAPI(title='Karyotype Viewer')

_HOST_HTML = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Karyotype Viewer</title>
<style>
  * {{ box-sizing:border-box; }}
  html,body {{ margin:0;padding:0;background:#EDF2F7;font-family:Inter,Arial,sans-serif; }}
  iframe {{ display:block;width:100%;border:0;margin:0;padding:0;background:#EDF2F7; }}
  #overview-frame {{ height:660px; }}
  #detail-frame   {{ height:740px; }}
</style>
</head>
<body>
  <iframe id="overview-frame" src="/overview/" title="Karyotype overview" scrolling="no"></iframe>
  <iframe id="detail-frame"   src="/detail/"   title="Chromosome detail"  scrolling="no"></iframe>

<script>
(function(){{
  const det  = document.getElementById('detail-frame');
  let   last = '{initial}';

  function fwd(chrom){{
    last = String(chrom||last);
    if(det&&det.contentWindow)
      det.contentWindow.postMessage(
        {{type:'karyotype-chrom-selected',chrom:last}},
        window.location.origin
      );
  }}

  window.addEventListener('message',e=>{{
    if(e.origin!==window.location.origin)return;
    const d=e.data||{{}};
    if(d.type!=='karyotype-chrom-selected')return;
    fwd(d.chrom);
  }});

  det.addEventListener('load',()=>{{ setTimeout(()=>fwd(last),150); }});
}})();
</script>
</body>
</html>""".format(initial=INITIAL_CHROM)


@fastapi_app.get('/', response_class=HTMLResponse)
async def root():
    return HTMLResponse(_HOST_HTML)


fastapi_app.mount('/overview', WSGIMiddleware(overview_app.server))
fastapi_app.mount('/detail',   WSGIMiddleware(detail_app.server))

app = fastapi_app  # uvicorn karyotype_app:app


if __name__ == '__main__':
    import uvicorn
    uvicorn.run(app, host='0.0.0.0', port=8050, reload=False)