"""
Karyotype Viewer
────────────────
구조:
  1. Sample Info   — 성별, 총 염색체 수, ISCN 핵형
  2. Karyogram     — 전체 염색체 overview (vertical, 클릭으로 선택)
  3. Chromosome    — 선택 염색체 전체 horizontal + brush 드래그
  4. Region Detail — brush 구간 CN / VAF scatter + 업로드 파일

pip install dash dash-bio plotly pandas numpy
python karyotype_viewer.py → http://localhost:8050
"""

import re, io, base64
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash import Dash, html, dcc, Input, Output, State, callback, no_update, clientside_callback
import dash_bio as dashbio

# ══════════════════════════════════════════════════════════════════
# 참조 데이터
# ══════════════════════════════════════════════════════════════════
CHROM_SIZES = {
    '1':248956422,'2':242193529,'3':198295559,'4':190214555,
    '5':181538259,'6':170805979,'7':159345973,'8':145138636,
    '9':138394717,'10':133797422,'11':135086622,'12':133275309,
    '13':114364328,'14':107043718,'15':101991189,'16':90338345,
    '17':83257441,'18':80373285,'19':58617616,'20':64444167,
    '21':46709983,'22':50818468,'X':156040895,'Y':57227415,
}
CENT_MB = {
    '1':123.4,'2':93.9,'3':90.9,'4':50.2,'5':48.8,
    '6':61.0,'7':59.9,'8':45.2,'9':43.0,'10':39.8,
    '11':53.4,'12':35.5,'13':17.7,'14':17.2,'15':19.0,
    '16':36.8,'17':25.1,'18':18.5,'19':26.2,'20':28.1,
    '21':12.0,'22':15.0,'X':61.0,'Y':10.4,
}
ALL_CHROMS   = [str(i) for i in range(1, 23)] + ['X', 'Y']
FEMALE_CHROMS = [str(i) for i in range(1, 23)] + ['X']
MALE_CHROMS   = [str(i) for i in range(1, 23)] + ['X', 'Y']

# ══════════════════════════════════════════════════════════════════
# 샘플 정의  ← 파이프라인 결과로 교체
# ══════════════════════════════════════════════════════════════════
SAMPLE = {
    'id':   'SAMPLE_001',
    'sex':  'female',   # 'male' | 'female'
    'events': {         # chrom → {type, cn, iscn, region}
        '21': {'type':'trisomy',      'cn':3,    'iscn':'+21',      'region':None,        'color':'#FC8181'},
        'X':  {'type':'monosomy',     'cn':1,    'iscn':'-X',       'region':None,        'color':'#90CDF4'},
        '5':  {'type':'partial_loss', 'cn':None, 'iscn':'del(5p)',  'region':(0.0,0.25),  'color':'#90CDF4'},
        '17': {'type':'partial_gain', 'cn':None, 'iscn':'dup(17q)', 'region':(0.50,1.0),  'color':'#FBD38D'},
    },
    'genes': [
        {'chr':'21','name':'DYRK1A','start':37700000,'stop':37865335,'color':'#B794F4'},
        {'chr':'21','name':'RUNX1', 'start':34787801,'stop':36004954,'color':'#B794F4'},
        {'chr':'17','name':'TP53',  'start':7661779, 'stop':7687538, 'color':'#68D391'},
        {'chr':'17','name':'BRCA1', 'start':43044295,'stop':43125483,'color':'#68D391'},
        {'chr':'17','name':'ERBB2', 'start':39687914,'stop':39730426,'color':'#68D391'},
        {'chr':'5', 'name':'CTNND2','start':11080000,'stop':11690000,'color':'#68D391'},
        {'chr':'X', 'name':'MECP2', 'start':154021573,'stop':154137103,'color':'#68D391'},
        {'chr':'X', 'name':'DMD',   'start':31119221,'stop':33339609, 'color':'#68D391'},
        {'chr':'13','name':'BRCA2', 'start':32315086,'stop':32400266, 'color':'#68D391'},
        {'chr':'7', 'name':'CFTR',  'start':117480025,'stop':117668665,'color':'#68D391'},
    ],
}

def build_iscn(sample):
    sex_str = 'XX' if sample['sex'] == 'female' else 'XY'
    tri = sum(1 for e in sample['events'].values() if e['type'] == 'trisomy')
    mon = sum(1 for e in sample['events'].values() if e['type'] == 'monosomy')
    total = 46 + tri - mon
    parts = []
    for c, ev in sorted(sample['events'].items()):
        if   ev['type'] == 'trisomy':  parts.append(f"+{c}")
        elif ev['type'] == 'monosomy': parts.append(f"-{c}")
        elif ev['type'] in ('partial_gain','partial_loss'): parts.append(ev['iscn'])
    return f"{total},{sex_str}" + (',' + ','.join(parts) if parts else '')

ISCN = build_iscn(SAMPLE)

# ══════════════════════════════════════════════════════════════════
# Ideogram annotation (rawAnnots 형식)
# ══════════════════════════════════════════════════════════════════
def make_raw_annots(sample):
    cd = {c: [] for c in ALL_CHROMS}

    for chrom, ev in sample['events'].items():
        size = CHROM_SIZES.get(chrom, 1)
        if ev['region'] is None:
            s, length = 1, size - 1
        else:
            s      = int(size * ev['region'][0]) + 1
            length = int(size * (ev['region'][1] - ev['region'][0]))
        cd[chrom].append([ev['iscn'], s, max(1, length), 0, ev['color']])

    for g in sample['genes']:
        length = max(100_000, g['stop'] - g['start'])
        cd[g['chr']].append([g['name'], g['start'], length, 1, g['color']])

    # histogram density bins (500 kb)
    for chrom, ev in sample['events'].items():
        size = CHROM_SIZES.get(chrom, 1)
        s    = 1 if ev['region'] is None else int(size * ev['region'][0]) + 1
        e    = size if ev['region'] is None else int(size * ev['region'][1])
        pos  = s
        while pos < e:
            be = min(pos + 500_000, e)
            cd[chrom].append([f'__d_{chrom}', pos, be - pos, 2, ev['color']])
            pos += 500_000

    return {
        'keys': ['name','start','length','trackIndex','color'],
        'annots': [{'chr':c,'annots':cd[c]} for c in ALL_CHROMS if cd[c]],
    }

RAW_ANNOTS = make_raw_annots(SAMPLE)

LEGEND = [{'name':'Legend','rows':[
    {'name':'Trisomy/Gain',  'color':'#FC8181','shape':'rectangle'},
    {'name':'Monosomy/Loss', 'color':'#90CDF4','shape':'rectangle'},
    {'name':'Partial Gain',  'color':'#FBD38D','shape':'rectangle'},
    {'name':'Gene',          'color':'#B794F4','shape':'rectangle'},
]}]

# ══════════════════════════════════════════════════════════════════
# 유틸
# ══════════════════════════════════════════════════════════════════
def parse_bp(v, default=0):
    if v is None: return default
    if isinstance(v, (int, float)): return int(v)
    try: return int(float(str(v).replace(',','').replace('bp','').strip()))
    except: return default

def parse_annot_html(raw):
    if not raw or not isinstance(raw, str):
        return None, None, None, None
    m = re.search(r'chr([0-9XY]+):([0-9,]+)-([0-9,]+)', raw)
    if not m: return None, None, None, None
    chrom = m.group(1)
    start = int(m.group(2).replace(',',''))
    stop  = int(m.group(3).replace(',',''))
    name  = re.sub(r'<[^>]+>','', raw[:m.start()]).strip().rstrip('<br/')
    return chrom, start, stop, name or f'chr{chrom}'

def _rgba(h, a):
    c = h.lstrip('#'); r,g,b = int(c[0:2],16),int(c[2:4],16),int(c[4:6],16)
    return f'rgba({r},{g},{b},{a})'

def arm(chrom, bp):
    return 'p' if bp/1e6 < CENT_MB.get(str(chrom), 50) else 'q'

# ══════════════════════════════════════════════════════════════════
# CN / VAF demo 데이터 (파이프라인 연동 전)
# ══════════════════════════════════════════════════════════════════
def demo_data(chrom, start_bp, end_bp):
    rng = np.random.default_rng(abs(hash(chrom)) % 2**31)
    n   = 400
    pos = np.linspace(start_bp, end_bp, n).astype(int)
    ev  = SAMPLE['events'].get(str(chrom))
    size = CHROM_SIZES.get(str(chrom), 1)

    cn  = np.full(n, 2.0)
    vaf = np.full(n, 0.5)

    if ev:
        frac = pos / size
        if ev['region'] is None:
            mask = np.ones(n, bool)
        else:
            mask = (frac >= ev['region'][0]) & (frac <= ev['region'][1])
        t = ev['type']
        if t in ('trisomy','partial_gain'):
            cn[mask]  = 3.0
            vaf[mask] = rng.choice([0.33,0.67], mask.sum())
        elif t in ('monosomy','partial_loss'):
            cn[mask]  = 1.0
            vaf[mask] = rng.choice([0.02,0.98], mask.sum())

    return pd.DataFrame({
        'pos': pos,
        'cn':  np.clip(cn  + rng.normal(0,0.2,n), 0, 6),
        'vaf': np.clip(vaf + rng.normal(0,0.04,n), 0, 1),
    })

def region_fig(df, chrom, start_bp, end_bp):
    sub = df[(df['pos']>=start_bp)&(df['pos']<=end_bp)].copy()
    if sub.empty: sub = df.copy()
    xmb = sub['pos']/1e6
    ev  = SAMPLE['events'].get(str(chrom))

    # 점 색상: annotation 구간 내 → event 색, 외부 → 회색
    def pt_c(p):
        if ev:
            frac = p / CHROM_SIZES.get(str(chrom),1)
            inside = (ev['region'] is None or
                      ev['region'][0] <= frac <= ev['region'][1])
            if inside: return ev['color']
        return '#CBD5E0'
    colors = sub['pos'].apply(pt_c)

    has_vaf = 'vaf' in sub.columns
    nrows = 2 if has_vaf else 1
    titles = [f'Copy Number — chr{chrom}',
              f'VAF — chr{chrom}'] if has_vaf else [f'Copy Number — chr{chrom}']

    fig = make_subplots(rows=nrows, cols=1, shared_xaxes=True,
                        row_heights=[0.5,0.5] if has_vaf else [1.0],
                        vertical_spacing=0.08, subplot_titles=titles)

    # CN
    fig.add_trace(go.Scatter(
        x=xmb, y=sub['cn'], mode='markers',
        marker=dict(size=3.5, color=colors, opacity=0.85),
        hovertemplate='%{x:.3f} Mb<br>CN: %{y:.3f}<extra></extra>',
    ), row=1, col=1)
    for v,c in [(1,_rgba('#3182CE',0.4)),(2,_rgba('#4A5568',0.5)),(3,_rgba('#E53E3E',0.4))]:
        fig.add_hline(y=v, line=dict(color=c,width=1,dash='dot'), row=1, col=1)

    # VAF
    if has_vaf:
        fig.add_trace(go.Scatter(
            x=xmb, y=sub['vaf'], mode='markers',
            marker=dict(size=3.5, color=colors, opacity=0.85),
            hovertemplate='%{x:.3f} Mb<br>VAF: %{y:.3f}<extra></extra>',
        ), row=2, col=1)
        for v,c in [(0.33,_rgba('#718096',0.4)),(0.5,'#4A5568'),(0.67,_rgba('#718096',0.4))]:
            fig.add_hline(y=v, line=dict(color=c,width=1,dash='dot'), row=2, col=1)

    # 이상 구간 shade
    if ev and ev['region']:
        size = CHROM_SIZES.get(str(chrom), 1)
        x0 = max(start_bp, size*ev['region'][0]) / 1e6
        x1 = min(end_bp,   size*ev['region'][1]) / 1e6
        for row in range(1, nrows+1):
            fig.add_vrect(x0=x0, x1=x1,
                          fillcolor=_rgba(ev['color'],0.12),
                          line_color=ev['color'], line_width=1.5,
                          row=row, col=1)
    elif ev and ev['region'] is None:
        for row in range(1, nrows+1):
            fig.add_vrect(x0=start_bp/1e6, x1=end_bp/1e6,
                          fillcolor=_rgba(ev['color'],0.08),
                          line_color='rgba(0,0,0,0)', line_width=0,
                          row=row, col=1)

    ax = dict(gridcolor='#E2E8F0', zeroline=False,
              tickfont=dict(size=10,color='#4A5568'), linecolor='#E2E8F0')
    fig.update_layout(
        paper_bgcolor='#F7FAFC', plot_bgcolor='white',
        font=dict(family='Inter,sans-serif',color='#2D3748',size=11),
        margin=dict(t=48,b=30,l=58,r=80),
        showlegend=False, hovermode='x unified', height=440,
        hoverlabel=dict(bgcolor='white',bordercolor='#E2E8F0',
                        font=dict(color='#2D3748',size=11)),
    )
    fig.update_xaxes(**ax)
    fig.update_yaxes(**ax)
    fig.update_yaxes(title_text='CN', range=[0,5.5], row=1, col=1)
    if has_vaf:
        fig.update_yaxes(title_text='VAF', range=[-0.05,1.05],
                         tickvals=[0,0.25,0.33,0.5,0.67,0.75,1],
                         ticktext=['0','0.25','⅓','0.5','⅔','0.75','1'],
                         row=2, col=1)
    fig.update_xaxes(title_text='Position (Mb)', row=nrows, col=1)
    return fig

EMPTY_FIG = go.Figure().update_layout(
    paper_bgcolor='#F7FAFC', plot_bgcolor='white', height=440,
    margin=dict(t=20,b=20,l=20,r=20),
    annotations=[dict(text='← 염색체를 선택하고 Brush로 구간을 드래그하세요',
                      x=0.5,y=0.5,xref='paper',yref='paper',
                      showarrow=False,font=dict(size=13,color='#A0AEC0'))],
)


# ── 헬퍼 ─────────────────────────────────────────────────────────
def _label(text):
    return html.Div(text, style={'fontSize':'10px','fontWeight':'700',
                                 'color':'#A0AEC0','letterSpacing':'.06em',
                                 'textTransform':'uppercase'})

def _chip(text, color, bg):
    return html.Span(text, style={
        'background':bg,'color':color,'padding':'1px 7px',
        'borderRadius':'4px','fontSize':'11px','fontWeight':'600',
    })

def _kv(k, v):
    return html.Div([
        html.Span(k+': ', style={'color':'#718096','fontWeight':'600'}),
        html.Span(v,      style={'color':'#1A202C'}),
    ], style={'fontSize':'12px'})

# ══════════════════════════════════════════════════════════════════
# Layout
# ══════════════════════════════════════════════════════════════════
app = Dash(__name__, title='Karyotype Viewer')

_BD   = '1px solid #E2E8F0'
_CARD = dict(background='white', border=_BD, borderRadius='8px', overflow='hidden',
             boxShadow='0 1px 3px rgba(0,0,0,0.06)')
_HEAD = dict(padding='9px 16px', borderBottom=_BD, fontSize='10px',
             fontWeight='700', color='#718096', letterSpacing='.1em',
             textTransform='uppercase', background='#F7FAFC',
             display='flex', alignItems='center', gap='8px')

def card(hdr, body):
    return html.Div(style=_CARD, children=[
        html.Div(style=_HEAD, children=hdr),
        html.Div(style={'padding':'14px 16px'}, children=body),
    ])

def badge(text, color='#3182CE'):
    return html.Span(text, style={
        'background':_rgba(color,0.1),'color':color,
        'border':f'1px solid {_rgba(color,0.3)}',
        'borderRadius':'4px','padding':'2px 9px',
        'fontSize':'12px','fontWeight':'600','fontFamily':'monospace',
    })

sex_label = '♀ Female' if SAMPLE['sex']=='female' else '♂ Male'
display_chroms = FEMALE_CHROMS if SAMPLE['sex']=='female' else MALE_CHROMS
n_chroms = 46 if SAMPLE['sex']=='female' else 46  # 실제 총 수는 ISCN에서

app.layout = html.Div(style={
    'fontFamily':'Inter,system-ui,sans-serif',
    'background':'#EDF2F7','minHeight':'100vh','color':'#2D3748','fontSize':'13px',
}, children=[

    # ── 헤더 ─────────────────────────────────────────────────────
    html.Div(style={
        'background':'white','borderBottom':_BD,'padding':'13px 24px',
        'boxShadow':'0 1px 3px rgba(0,0,0,.05)',
        'display':'flex','alignItems':'center','gap':'10px',
    }, children=[
        html.Span('🧬', style={'fontSize':'18px'}),
        html.H1('Karyotype Viewer',
                style={'margin':0,'fontSize':'15px','fontWeight':'700'}),
    ]),

    html.Div(style={'padding':'14px 24px','display':'flex',
                    'flexDirection':'column','gap':'12px'}, children=[

        # ══ 1. SAMPLE INFO ══════════════════════════════════════
        card(
            hdr=['SAMPLE INFO'],
            body=[
                html.Div(style={
                    'display':'grid',
                    'gridTemplateColumns':'repeat(auto-fill,minmax(200px,1fr))',
                    'gap':'16px',
                }, children=[

                    # 기본 정보
                    html.Div([
                        _label('Sample ID'),
                        html.Div(SAMPLE['id'],
                                 style={'fontWeight':'700','fontSize':'14px','marginTop':'2px'}),
                    ]),
                    html.Div([
                        _label('Sex'),
                        html.Div(sex_label,
                                 style={'fontWeight':'700','fontSize':'14px','marginTop':'2px',
                                        'color':'#E53E3E' if SAMPLE['sex']=='female' else '#3182CE'}),
                    ]),
                    html.Div([
                        _label('Karyotype (ISCN)'),
                        html.Div(badge(ISCN), style={'marginTop':'4px'}),
                    ]),
                    html.Div([
                        _label('Total chromosomes'),
                        html.Div(
                            str(int(ISCN.split(',')[0])),
                            style={'fontWeight':'700','fontSize':'20px',
                                   'color': '#E53E3E' if int(ISCN.split(',')[0])!=46
                                           else '#2D3748','marginTop':'2px'},
                        ),
                    ]),

                    # CN 이상 이벤트 목록
                    html.Div(style={'gridColumn':'1/-1'}, children=[
                        _label('Chromosomal Events'),
                        html.Div(style={
                            'display':'flex','gap':'8px','flexWrap':'wrap','marginTop':'6px',
                        }, children=[
                            html.Div(style={
                                'border':f"1px solid {ev['color']}88",
                                'borderLeft':f"4px solid {ev['color']}",
                                'borderRadius':'5px','padding':'7px 12px',
                                'background':_rgba(ev['color'],0.06),
                                'minWidth':'140px',
                            }, children=[
                                html.Div(ev['iscn'],
                                         style={'fontWeight':'700','fontFamily':'monospace',
                                                'fontSize':'13px','color':ev['color']}),
                                html.Div(
                                    f"chr{chrom}  {ev['type'].replace('_',' ')}",
                                    style={'fontSize':'11px','color':'#718096','marginTop':'2px'},
                                ),
                                html.Div(
                                    f"CN = {ev['cn'] if ev['cn'] else 'partial'}",
                                    style={'fontSize':'11px','color':'#4A5568','marginTop':'1px'},
                                ),
                            ])
                            for chrom, ev in sorted(SAMPLE['events'].items())
                        ]),
                    ]),
                ]),
            ],
        ),

        # ══ 2. KARYOGRAM ════════════════════════════════════════
        card(
            hdr=[
                'KARYOGRAM',
                html.Span(id='karyogram-hint',
                          children='염색체를 클릭하여 선택하세요',
                          style={'marginLeft':'auto','fontWeight':'400',
                                 'fontSize':'11px','color':'#A0AEC0',
                                 'textTransform':'none','letterSpacing':'0'}),
                # 선택 해제 버튼 (선택 시에만 표시)
                html.Button(
                    '✕ 선택 해제',
                    id='clear-chrom-btn',
                    n_clicks=0,
                    style={
                        'display':'none',
                        'fontSize':'11px','padding':'3px 10px',
                        'background':'white','color':'#718096',
                        'border':'1px solid #CBD5E0','borderRadius':'4px',
                        'cursor':'pointer','marginLeft':'8px',
                    },
                ),
            ],
            body=[
                # 1줄 전체 표시: overflowX + minWidth 으로 가로 스크롤
                html.Div(style={'overflowX':'auto', 'overflowY':'hidden',
                                'minWidth':'0'}, children=[
                    dashbio.Ideogram(
                        id='karyogram',
                        organism='human', assembly='GRCh38',
                        orientation='vertical',
                        chromosomes=display_chroms,
                        chrHeight=350,            # 충분한 높이
                        chrWidth=15,
                        chrMargin=15,
                        rows=1,                   # ← 1줄로 전부 표기
                        rotatable=True,
                        showBandLabels=True,       # ← p/q band label 표시
                        showChromosomeLabels=True,
                        resolution=550,
                        annotations=RAW_ANNOTS,
                        annotationsLayout='overlay',
                        annotationHeight=4,
                        showAnnotTooltip=True,
                        legend=LEGEND,
                    ),
                ]),
            ],
        ),

        # ══ 3. CHROMOSOME VIEW (brush) ══════════════════════════
        html.Div(id='chrom-panel', style={'visibility':'hidden','height':'0','overflow':'hidden','marginBottom':'0'}, children=[
            card(
                hdr=[
                    html.Span(id='chrom-panel-title', children='CHROMOSOME VIEW'),
                    html.Span(id='brush-range-chip', style={
                        'marginLeft':'auto','fontFamily':'monospace',
                        'color':'#3182CE','fontSize':'12px',
                        'fontWeight':'500','textTransform':'none','letterSpacing':'0',
                    }),
                ],
                body=[
                    # 선택 염색체 전체를 가로로 표시 + brush
                    html.Div(style={'overflowX':'auto','overflowY':'hidden'}, children=[
                        dashbio.Ideogram(
                            id='chrom-ideo',
                            organism='human', assembly='GRCh38',
                            orientation='horizontal',  # ← 처음부터 고정
                            chrHeight=120,
                            chrWidth=1140,             # 컨테이너 가득 채움
                            rows=1,
                            rotatable=False,
                            showBandLabels=True,
                            showChromosomeLabels=True,
                            resolution=850,            # 고해상도 band
                            annotations=RAW_ANNOTS,
                            annotationsLayout='overlay',
                            annotationHeight=14,
                            showAnnotTooltip=True,
                            chromosomes=['1'],
                            brush='chr1:1-248956422',
                        ),
                    ]),
                    html.Div(
                        '← → 핸들을 드래그해서 구간 선택 · 드래그로 이동',
                        style={'marginTop':'6px','fontSize':'11px',
                               'color':'#A0AEC0','textAlign':'center'},
                    ),
                ],
            ),
        ]),

        # ══ 4. REGION DETAIL ═══════════════════════════════════
        html.Div(id='detail-panel', style={'visibility':'hidden','height':'0','overflow':'hidden','marginBottom':'0'}, children=[
            card(
                hdr=[
                    html.Span(id='detail-title', children='REGION DETAIL'),
                    # 파일 업로드
                    html.Div(style={'marginLeft':'auto','display':'flex',
                                    'gap':'8px','alignItems':'center'}, children=[
                        dcc.Upload(
                            id='upload',
                            children=html.Div([
                                html.Span('📂 TSV/CSV 업로드',
                                          style={'fontSize':'11px','fontWeight':'600'}),
                                html.Span(' (pos, cn, vaf)',
                                          style={'fontSize':'10px','color':'#A0AEC0'}),
                            ]),
                            style={
                                'border':'1px dashed #CBD5E0','borderRadius':'5px',
                                'padding':'4px 10px','cursor':'pointer',
                                'background':'#F7FAFC','textTransform':'none',
                                'letterSpacing':'0',
                            },
                            multiple=False,
                        ),
                        html.Button('Demo', id='demo-btn', n_clicks=0,
                                    style={
                                        'fontSize':'11px','padding':'4px 12px',
                                        'background':'#3182CE','color':'white',
                                        'border':'none','borderRadius':'4px',
                                        'cursor':'pointer','textTransform':'none',
                                    }),
                        html.Span(id='upload-status',
                                  style={'fontSize':'11px','color':'#718096',
                                         'textTransform':'none','letterSpacing':'0'}),
                    ]),
                ],
                body=[
                    dcc.Graph(id='detail-graph', figure=EMPTY_FIG,
                              config={'displayModeBar':True,
                                      'modeBarButtonsToRemove':['select2d','lasso2d'],
                                      'responsive':True}),
                ],
            ),
        ]),

    ]),

    # ── hidden stores ────────────────────────────────────────────
    dcc.Store(id='selected-chrom', data=None),  # 현재 선택 염색체
    dcc.Store(id='brush-region',   data=None),  # {chrom,start,end}
    dcc.Store(id='uploaded-df',    data=None),  # JSON
])



# ══════════════════════════════════════════════════════════════════
# Callbacks
# ══════════════════════════════════════════════════════════════════

# ── 2→3: Karyogram 클릭 / 선택 해제 ─────────────────────────────
@callback(
    Output('selected-chrom',    'data'),
    Output('chrom-panel',       'style'),
    Output('detail-panel',      'style',   allow_duplicate=True),
    Output('chrom-panel-title', 'children'),
    Output('chrom-ideo',        'chromosomes'),
    Output('chrom-ideo',        'brush'),
    Output('clear-chrom-btn',   'style'),   # 해제 버튼 표시/숨김
    Output('karyogram-hint',    'children'),
    Input('karyogram',          'annotationsData'),
    Input('karyogram',          'rotated'),
    Input('clear-chrom-btn',    'n_clicks'),
    State('selected-chrom',     'data'),
    prevent_initial_call=True,
)
def select_chrom(annot_html, rotated, clear_clicks, prev_chrom):
    from dash import ctx

    # ── 선택 해제 ──────────────────────────────────────────────
    if ctx.triggered_id == 'clear-chrom-btn':
        btn_hidden = {
            'display':'none','fontSize':'11px','padding':'3px 10px',
            'background':'white','color':'#718096',
            'border':'1px solid #CBD5E0','borderRadius':'4px',
            'cursor':'pointer','marginLeft':'8px',
        }
        return (None,                        # selected-chrom
                {'visibility':'hidden','height':'0','overflow':'hidden','marginBottom':'0'},  # chrom-panel 숨김
                {'visibility':'hidden','height':'0','overflow':'hidden','marginBottom':'0'},  # detail-panel 숨김
                'CHROMOSOME VIEW',           # title 초기화
                no_update,                   # chromosomes (그대로)
                no_update,                   # brush (그대로)
                btn_hidden,                  # 버튼 숨김
                '염색체를 클릭하여 선택하세요')  # 힌트 초기화

    # triggered prop으로 hover vs 클릭 정확히 구분
    triggered_prop = ctx.triggered[0]['prop_id'] if ctx.triggered else ''

    # ── hover: chrom만 store에 저장, 패널은 열지 않음 ──────────
    if 'annotationsData' in triggered_prop:
        chrom, _, _, _ = parse_annot_html(annot_html)
        if not chrom:
            return (no_update,)*8
        return (chrom,) + (no_update,)*7   # store만 업데이트

    # ── rotated(클릭): 패널 열기 ───────────────────────────────
    chrom = prev_chrom   # hover로 저장된 chrom 사용
    if not chrom:
        return (no_update,)*8

    # ── 클릭: 패널 열기 ───────────────────────────────────────
    size  = CHROM_SIZES.get(chrom, 100_000_000)
    brush = f'chr{chrom}:1-{size}'   # ← 전체 염색체 범위
    ev    = SAMPLE['events'].get(chrom)

    ev_span = html.Span(
        f"  [{ev['iscn']}]" if ev else '',
        style={'color': ev['color'] if ev else 'transparent',
               'fontFamily':'monospace','fontSize':'12px','fontWeight':'600'},
    )
    title = [f'CHROMOSOME VIEW — chr{chrom}', ev_span,
             html.Span('  ← → 드래그로 구간 선택',
                       style={'color':'#A0AEC0','fontWeight':'400',
                              'textTransform':'none','letterSpacing':'0'})]

    btn_visible = {
        'display':'inline-block','fontSize':'11px','padding':'3px 10px',
        'background':'white','color':'#718096',
        'border':'1px solid #CBD5E0','borderRadius':'4px',
        'cursor':'pointer','marginLeft':'8px',
    }
    hint = f'chr{chrom} 선택됨'

    return (chrom,
            {'visibility':'visible','height':'auto','overflow':'visible','marginBottom':'0'},   # chrom-panel
            no_update,             # detail-panel (brush 드래그 전까지 유지)
            title,
            [chrom],               # chrom-ideo.chromosomes
            brush,                 # chrom-ideo.brush = 전체 범위
            btn_visible,           # 해제 버튼 표시
            hint)


# ── 3→4: Brush 드래그 → region store + detail panel ──────────────
@callback(
    Output('brush-region',    'data'),
    Output('brush-range-chip','children'),
    Output('detail-panel',    'style',    allow_duplicate=True),
    Output('detail-title',    'children'),
    Input('chrom-ideo',       'brushData'),
    State('selected-chrom',   'data'),
    prevent_initial_call=True,
)
def on_brush(brush_data, chrom):
    if not brush_data or not chrom:
        return no_update, no_update, no_update, no_update

    start_bp = parse_bp(brush_data.get('start'), 0)
    end_bp   = parse_bp(brush_data.get('end'),   0)
    extent   = parse_bp(brush_data.get('extent'), max(0, end_bp - start_bp))

    if end_bp <= start_bp:
        return no_update, no_update, no_update, no_update

    chip = (f"chr{chrom}: {start_bp/1e6:.3f} – {end_bp/1e6:.3f} Mb"
            f"  ({extent/1e6:.3f} Mb)")
    ev   = SAMPLE['events'].get(chrom)

    title = [f'REGION DETAIL — chr{chrom}: {start_bp/1e6:.2f}–{end_bp/1e6:.2f} Mb']
    if ev:
        title.append(html.Span(
            f"  {ev['iscn']}",
            style={'color':ev['color'],'fontFamily':'monospace',
                   'fontSize':'12px','fontWeight':'600'},
        ))

    store = {'chrom': chrom, 'start': start_bp, 'end': end_bp}
    return store, chip, {'visibility':'visible','height':'auto','overflow':'visible','marginBottom':'0'}, title


# ── 4: 파일 업로드 ─────────────────────────────────────────────────
@callback(
    Output('uploaded-df',   'data'),
    Output('upload-status', 'children'),
    Input('upload',         'contents'),
    State('upload',         'filename'),
    prevent_initial_call=True,
)
def upload_file(contents, filename):
    if not contents: return no_update, no_update
    _, b64 = contents.split(',')
    decoded = base64.b64decode(b64)
    try:
        sep = '\t' if (filename or '').endswith('.tsv') else ','
        df  = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep=sep)
        df.columns = [c.lower().strip() for c in df.columns]
        # pos 컬럼 표준화
        for alt in ['position','start','chromstart']:
            if alt in df.columns and 'pos' not in df.columns:
                df = df.rename(columns={alt:'pos'})
        # cn 컬럼 표준화
        for alt in ['copy_number','log2','ratio']:
            if alt in df.columns and 'cn' not in df.columns:
                df = df.rename(columns={alt:'cn'})
        # vaf/baf
        for alt in ['baf','b_allele_freq','allele_freq']:
            if alt in df.columns and 'vaf' not in df.columns:
                df = df.rename(columns={alt:'vaf'})
        missing = [c for c in ['pos','cn'] if c not in df.columns]
        if missing: return None, f'❌ 컬럼 없음: {missing}'
        cols = ['pos','cn'] + (['vaf'] if 'vaf' in df.columns else [])
        df = df[cols].dropna()
        return df.to_json(), f'✅ {filename} ({len(df):,} rows)'
    except Exception as ex:
        return None, f'❌ {ex}'


# ── 4: Region Detail 그래프 렌더 ──────────────────────────────────
@callback(
    Output('detail-graph', 'figure'),
    Input('brush-region',  'data'),
    Input('demo-btn',      'n_clicks'),
    State('uploaded-df',   'data'),
    prevent_initial_call=True,
)
def render_detail(region, demo_clicks, df_json):
    from dash import ctx
    if not region:
        return EMPTY_FIG

    chrom    = region['chrom']
    start_bp = region['start']
    end_bp   = region['end']

    use_demo = (ctx.triggered_id == 'demo-btn') or not df_json
    if use_demo:
        df = demo_data(chrom, start_bp, end_bp)
    else:
        try:
            df = pd.read_json(df_json)
        except:
            df = demo_data(chrom, start_bp, end_bp)

    return region_fig(df, chrom, start_bp, end_bp)


# ── 패널이 보인 직후 resize 이벤트 발생 → ideogram brush 재계산 ──
clientside_callback(
    '''
    function(style) {
        if (style && style.visibility === 'visible') {
            setTimeout(function() {
                window.dispatchEvent(new Event('resize'));
            }, 150);
        }
        return window.dash_clientside.no_update;
    }
    ''',
    Output('chrom-ideo', 'id'),          # 더미 output (id는 변경 안 됨)
    Input('chrom-panel', 'style'),
)

if __name__ == '__main__':
    app.run(debug=True, port=8052)
