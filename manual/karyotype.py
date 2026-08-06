"""
pip install dash dash-bio plotly pandas numpy
python karyotype_dash.py  →  http://localhost:8050
"""

import json
import numpy as np
import pandas as pd
from dash import Dash, html, dcc, Input, Output, State, callback, no_update, ctx
import dash_bio as dashbio
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ══════════════════════════════════════════════════════════════════
# 1. 샘플 데이터  ← 파이프라인 결과로 교체
# ══════════════════════════════════════════════════════════════════
CHROM_SIZES = {
    '1':248956422,'2':242193529,'3':198295559,'4':190214555,
    '5':181538259,'6':170805979,'7':159345973,'8':145138636,
    '9':138394717,'10':133797422,'11':135086622,'12':133275309,
    '13':114364328,'14':107043718,'15':101991189,'16':90338345,
    '17':83257441,'18':80373285,'19':58617616,'20':64444167,
    '21':46709983,'22':50818468,'X':156040895,'Y':57227415,
}

# dash_bio annotation color = hex string, shape = 'circle'|'triangle'|'rectangle'
EVENTS = {
    '21': {'type':'trisomy',      'cn':3,    'iscn':'+21',      'region':None,        'color':'#E53E3E','shape':'triangle'},
    '5':  {'type':'partial_loss', 'cn':None, 'iscn':'del(5p)',  'region':(0.0, 0.25), 'color':'#3182CE','shape':'rectangle'},
    '17': {'type':'partial_gain', 'cn':None, 'iscn':'dup(17q)', 'region':(0.5, 1.0),  'color':'#DD6B20','shape':'rectangle'},
    'X':  {'type':'monosomy',     'cn':1,    'iscn':'-X',       'region':None,        'color':'#3182CE','shape':'triangle'},
}

SEX       = 'female'
SAMPLE_ID = 'SAMPLE_001'
ASSEMBLY  = 'GRCh38'

TYPE_LABEL = {
    'trisomy':      'Trisomy (CN=3)',
    'monosomy':     'Monosomy (CN=1)',
    'partial_gain': 'Partial Gain',
    'partial_loss': 'Partial Loss',
}

# ══════════════════════════════════════════════════════════════════
# 2. Annotation 목록
# ══════════════════════════════════════════════════════════════════
def build_annotations():
    anns = []
    for chrom, ev in EVENTS.items():
        size = CHROM_SIZES.get(chrom, 100_000_000)
        if ev['region'] is None:
            start, stop = 1, size
        else:
            start = int(size * ev['region'][0]) + 1
            stop  = int(size * ev['region'][1])
        anns.append({
            'chr':   chrom,
            'name':  ev['iscn'],
            'start': start,
            'stop':  stop,
            'color': ev['color'],
            'shape': ev['shape'],
        })
    return anns

ANNOTATIONS = build_annotations()

LEGEND = [{
    'name': 'CN Status',
    'rows': [
        {'name':'Trisomy/Gain',  'color':'#E53E3E','shape':'triangle'},
        {'name':'Monosomy/Loss', 'color':'#3182CE','shape':'triangle'},
        {'name':'Partial Gain',  'color':'#DD6B20','shape':'rectangle'},
        {'name':'Partial Loss',  'color':'#3182CE','shape':'rectangle'},
    ],
}]

# ══════════════════════════════════════════════════════════════════
# 3. 더미 CN/BAF  ← 실제 bin TSV로 교체
# ══════════════════════════════════════════════════════════════════
def generate_data(chrom: str) -> pd.DataFrame:
    rng  = np.random.default_rng(abs(hash(chrom)) % (2**31))
    size = CHROM_SIZES.get(chrom, 100_000_000)
    n    = 400
    pos  = np.linspace(1, size, n).astype(int)
    ev   = EVENTS.get(chrom)

    base_cn = np.full(n, 2.0)
    baf     = np.full(n, 0.5)
    labels  = ['normal'] * n

    if ev:
        frac = pos / size
        mask = (np.ones(n, bool) if ev['region'] is None
                else ((frac >= ev['region'][0]) & (frac <= ev['region'][1])))
        t = ev['type']
        if t == 'trisomy':
            base_cn[mask] = 3.0
            baf[mask] = rng.choice([0.33, 0.67], size=mask.sum())
        elif t == 'monosomy':
            base_cn[mask] = 1.0
            baf[mask] = rng.choice([0.03, 0.97], size=mask.sum())
        elif t == 'partial_gain':
            base_cn[mask] = 3.0
            baf[mask] = rng.choice([0.33, 0.67], size=mask.sum())
        elif t == 'partial_loss':
            base_cn[mask] = 1.0
            baf[mask] = rng.choice([0.04, 0.96], size=mask.sum())
        for i in range(n):
            if mask[i]:
                labels[i] = t

    return pd.DataFrame({
        'pos':     pos,
        'cn':      np.clip(base_cn + rng.normal(0, 0.22, n), 0, 6),
        'baf':     np.clip(baf     + rng.normal(0, 0.04, n), 0, 1),
        'anomaly': labels,
    })

# ══════════════════════════════════════════════════════════════════
# 4. Plotly 도우미
# ══════════════════════════════════════════════════════════════════
PT_COLOR = {
    'trisomy':      '#E53E3E',
    'monosomy':     '#3182CE',
    'partial_gain': '#DD6B20',
    'partial_loss': '#3182CE',
    'normal':       '#CBD5E0',
}

_PLOT_BG  = '#FFFFFF'
_PAPER_BG = '#F7FAFC'
_GRID     = '#E2E8F0'
_TICK     = '#4A5568'
_AXIS_BASE = dict(gridcolor=_GRID, zeroline=False,
                  linecolor='#CBD5E0', tickfont=dict(size=11, color=_TICK))

def _axis(**kw): return {**_AXIS_BASE, **kw}


def make_detail_fig(chrom: str) -> go.Figure:
    df   = generate_data(chrom)
    ev   = EVENTS.get(chrom)
    size = CHROM_SIZES.get(chrom, 100_000_000)
    xmb  = df['pos'] / 1e6
    colors = df['anomaly'].map(PT_COLOR).fillna('#CBD5E0')

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        row_heights=[0.5, 0.5],
        vertical_spacing=0.10,
        subplot_titles=[
            f'Copy Number — chr{chrom}',
            f'B-Allele Frequency — chr{chrom}',
        ],
    )

    # CN scatter
    fig.add_trace(go.Scatter(
        x=xmb, y=df['cn'], mode='markers',
        marker=dict(size=3.5, color=colors, opacity=0.85),
        hovertemplate='%{x:.2f} Mb<br>CN: %{y:.3f}<extra></extra>',
        name='CN',
    ), row=1, col=1)

    # CN 기준선
    for v, c, lbl in [(1,'#3182CE','CN=1'),(2,'#718096','CN=2'),(3,'#E53E3E','CN=3')]:
        fig.add_hline(y=v, line=dict(color=c, width=1, dash='dot'),
                      annotation_text=lbl,
                      annotation_font=dict(size=9, color=c),
                      annotation_position='right',
                      row=1, col=1)

    # BAF scatter
    fig.add_trace(go.Scatter(
        x=xmb, y=df['baf'], mode='markers',
        marker=dict(size=3.5, color=colors, opacity=0.85),
        hovertemplate='%{x:.2f} Mb<br>BAF: %{y:.3f}<extra></extra>',
        name='BAF',
    ), row=2, col=1)

    for v, c, lbl in [(0.33,'#718096','1/3'),(0.5,'#2D3748','0.5'),(0.67,'#718096','2/3')]:
        fig.add_hline(y=v, line=dict(color=c, width=1, dash='dot'),
                      annotation_text=lbl,
                      annotation_font=dict(size=9, color=c),
                      annotation_position='right',
                      row=2, col=1)

    # partial 구간 shade
    if ev and ev['region']:
        x0 = size * ev['region'][0] / 1e6
        x1 = size * ev['region'][1] / 1e6
        ec = ev['color']
        for row in [1, 2]:
            fig.add_vrect(x0=x0, x1=x1,
                          fillcolor=ec + '22', line_color=ec,
                          line_width=1.5, row=row, col=1)

    fig.update_layout(
        paper_bgcolor=_PAPER_BG,
        plot_bgcolor=_PLOT_BG,
        font=dict(family='Inter, system-ui, sans-serif', color='#2D3748', size=12),
        margin=dict(t=50, b=36, l=60, r=80),
        showlegend=False,
        hovermode='x unified',
        hoverlabel=dict(bgcolor='white', bordercolor='#CBD5E0',
                        font=dict(color='#2D3748', size=11)),
        height=480,
    )
    fig.update_xaxes(**_axis())
    fig.update_yaxes(**_axis())
    fig.update_yaxes(title_text='CN', range=[0, 5.5], row=1, col=1)
    fig.update_yaxes(
        title_text='BAF', range=[-0.05, 1.05],
        tickvals=[0, 0.25, 0.33, 0.5, 0.67, 0.75, 1],
        ticktext=['0','0.25','⅓','0.5','⅔','0.75','1'],
        row=2, col=1,
    )
    fig.update_xaxes(title_text='Position (Mb)', row=2, col=1)
    fig.update_annotations(font=dict(size=13, color='#2D3748'))  # subplot titles
    return fig


def make_region_fig(chrom: str) -> go.Figure:
    """두 번째 클릭 시: 이상 구간 줌인 CN+BAF."""
    df   = generate_data(chrom)
    ev   = EVENTS.get(chrom)
    size = CHROM_SIZES.get(chrom, 100_000_000)

    if not ev or not ev['region']:
        # 전체 염색체 이상이면 그냥 전체 보여줌
        return make_detail_fig(chrom)

    r0, r1 = ev['region']
    x0_mb = size * r0 / 1e6
    x1_mb = size * r1 / 1e6
    padding = (x1_mb - x0_mb) * 0.3
    mask = ((df['pos'] / 1e6) >= (x0_mb - padding)) & \
           ((df['pos'] / 1e6) <= (x1_mb + padding))
    sub = df[mask]
    xmb = sub['pos'] / 1e6
    colors = sub['anomaly'].map(PT_COLOR).fillna('#CBD5E0')

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True,
        row_heights=[0.5, 0.5], vertical_spacing=0.10,
        subplot_titles=[
            f'CN 줌인 — chr{chrom} {ev["iscn"]} 구간',
            f'BAF 줌인 — chr{chrom} {ev["iscn"]} 구간',
        ],
    )

    fig.add_trace(go.Scatter(
        x=xmb, y=sub['cn'], mode='markers',
        marker=dict(size=5, color=colors, opacity=0.9),
        hovertemplate='%{x:.2f} Mb<br>CN: %{y:.3f}<extra></extra>',
        name='CN',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=xmb, y=sub['baf'], mode='markers',
        marker=dict(size=5, color=colors, opacity=0.9),
        hovertemplate='%{x:.2f} Mb<br>BAF: %{y:.3f}<extra></extra>',
        name='BAF',
    ), row=2, col=1)

    ec = ev['color']
    for row in [1, 2]:
        fig.add_vrect(x0=x0_mb, x1=x1_mb,
                      fillcolor=ec + '18', line_color=ec,
                      line_width=2, row=row, col=1)

    for v, c, lbl in [(1,'#3182CE','CN=1'),(2,'#718096','CN=2'),(3,'#E53E3E','CN=3')]:
        fig.add_hline(y=v, line=dict(color=c, width=1, dash='dot'),
                      annotation_text=lbl,
                      annotation_font=dict(size=9, color=c),
                      annotation_position='right', row=1, col=1)
    for v, c, lbl in [(0.33,'#718096','1/3'),(0.5,'#2D3748','0.5'),(0.67,'#718096','2/3')]:
        fig.add_hline(y=v, line=dict(color=c, width=1, dash='dot'),
                      annotation_text=lbl,
                      annotation_font=dict(size=9, color=c),
                      annotation_position='right', row=2, col=1)

    # 구간 정보 annotation
    region_bp  = int(size * (r1 - r0))
    region_mb  = region_bp / 1e6
    fig.add_annotation(
        x=0.5, y=1.04, xref='paper', yref='paper',
        text=f"<b>{ev['iscn']}</b>  |  {ev['type'].replace('_',' ').title()}  |  "
             f"{x0_mb:.1f}–{x1_mb:.1f} Mb  ({region_mb:.1f} Mb)",
        showarrow=False, align='center',
        font=dict(size=12, color=ec),
        bgcolor='white', bordercolor=ec,
        borderwidth=1, borderpad=5,
    )

    fig.update_layout(
        paper_bgcolor=_PAPER_BG, plot_bgcolor=_PLOT_BG,
        font=dict(family='Inter, system-ui, sans-serif', color='#2D3748', size=12),
        margin=dict(t=70, b=36, l=60, r=80),
        showlegend=False,
        hovermode='x unified',
        hoverlabel=dict(bgcolor='white', bordercolor='#CBD5E0',
                        font=dict(color='#2D3748', size=11)),
        height=480,
    )
    fig.update_xaxes(**_axis())
    fig.update_yaxes(**_axis())
    fig.update_yaxes(title_text='CN',  range=[0, 5.5], row=1, col=1)
    fig.update_yaxes(
        title_text='BAF', range=[-0.05, 1.05],
        tickvals=[0, 0.25, 0.33, 0.5, 0.67, 0.75, 1],
        ticktext=['0','0.25','⅓','0.5','⅔','0.75','1'],
        row=2, col=1,
    )
    fig.update_xaxes(title_text='Position (Mb)', row=2, col=1)
    fig.update_annotations(font=dict(size=13, color='#2D3748'))
    return fig


EMPTY_FIG = go.Figure().update_layout(
    paper_bgcolor=_PAPER_BG, plot_bgcolor=_PLOT_BG, height=480,
    margin=dict(t=20, b=20, l=20, r=20),
    annotations=[dict(
        text='← 염색체를 클릭하면 CN / BAF 상세 정보가 표시됩니다',
        x=0.5, y=0.5, xref='paper', yref='paper',
        showarrow=False, font=dict(size=14, color='#A0AEC0'),
    )],
)

# ══════════════════════════════════════════════════════════════════
# 5. ISCN notation
# ══════════════════════════════════════════════════════════════════
def build_iscn():
    total, sx = 46, 'XX' if SEX == 'female' else 'XY'
    parts = []
    for c, ev in sorted(EVENTS.items()):
        if ev['type'] == 'trisomy':  total += 1; parts.append(f'+{c}')
        elif ev['type'] == 'monosomy': total -= 1; parts.append(f'-{c}')
    return f"{total},{sx}" + ((',' + ','.join(parts)) if parts else '')

ISCN = build_iscn()

# ══════════════════════════════════════════════════════════════════
# 6. Layout
# ══════════════════════════════════════════════════════════════════
DEFAULT_CHROMS = [str(i) for i in range(1, 23)] + ['X']
CHROM_OPTIONS  = [{'label': str(i), 'value': str(i)} for i in range(1, 23)] + \
                 [{'label': 'X', 'value': 'X'}, {'label': 'Y', 'value': 'Y'}]

# Dash external stylesheets — Google Fonts
app = Dash(
    __name__,
    title='Karyotype Viewer',
    external_stylesheets=['https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap'],
)

_BD = '1px solid #E2E8F0'

def card(children, style=None):
    base = {
        'background': 'white',
        'border': _BD,
        'borderRadius': '8px',
        'overflow': 'hidden',
        'boxShadow': '0 1px 3px rgba(0,0,0,0.08)',
    }
    if style:
        base.update(style)
    return html.Div(style=base, children=children)

def card_header(children, extra_style=None):
    base = {
        'padding': '10px 16px',
        'borderBottom': _BD,
        'background': '#F7FAFC',
        'fontSize': '11px',
        'fontWeight': '700',
        'color': '#4A5568',
        'letterSpacing': '.08em',
        'textTransform': 'uppercase',
        'display': 'flex',
        'alignItems': 'center',
        'gap': '10px',
    }
    if extra_style:
        base.update(extra_style)
    return html.Div(style=base, children=children)

def badge(text, color='#3182CE'):
    return html.Span(text, style={
        'background': color + '18',
        'color': color,
        'border': f'1px solid {color}44',
        'borderRadius': '4px',
        'padding': '1px 8px',
        'fontSize': '11px',
        'fontWeight': '600',
        'fontFamily': 'monospace',
    })

app.layout = html.Div(
    style={
        'background': '#EDF2F7',
        'minHeight': '100vh',
        'fontFamily': 'Inter, system-ui, sans-serif',
        'color': '#2D3748',
        'fontSize': '13px',
    },
    children=[

    # ── Header ──────────────────────────────────────────────────
    html.Div(style={
        'background': 'white',
        'borderBottom': _BD,
        'padding': '14px 24px',
        'display': 'flex',
        'alignItems': 'center',
        'gap': '12px',
        'flexWrap': 'wrap',
        'boxShadow': '0 1px 2px rgba(0,0,0,0.05)',
    }, children=[
        html.Div(style={'display':'flex','alignItems':'center','gap':'8px'}, children=[
            html.Div(style={
                'width':'8px','height':'8px','borderRadius':'50%',
                'background':'#48BB78',
            }),
            html.H1('Karyotype Viewer', style={
                'fontSize':'15px','fontWeight':'700','margin':0,'color':'#1A202C',
            }),
        ]),
        badge(ISCN, '#3182CE'),
        html.Div(style={'marginLeft':'auto','display':'flex','gap':'16px',
                        'alignItems':'center','flexWrap':'wrap'}, children=[
            *[html.Div(style={'display':'flex','alignItems':'center',
                              'gap':'5px','fontSize':'11px','color':'#4A5568'}, children=[
                html.Div(style={'width':'9px','height':'9px','borderRadius':'2px',
                                'background':c, 'flexShrink':'0'}),
                lbl,
            ]) for c, lbl in [
                ('#E53E3E','Trisomy/Gain'),
                ('#3182CE','Monosomy/Loss'),
                ('#DD6B20','Partial Gain'),
                ('#3182CE','Partial Loss'),
            ]],
        ]),
    ]),

    html.Div(style={'padding':'16px 24px','display':'flex',
                    'flexDirection':'column','gap':'14px'}, children=[

        # ── Ideogram Panel ──────────────────────────────────────
        card([
            card_header([
                'Ideogram',
                dcc.Dropdown(
                    id='chrom-dropdown',
                    options=CHROM_OPTIONS,
                    multi=True,
                    value=DEFAULT_CHROMS,
                    placeholder='염색체 선택…',
                    style={'width':'380px','fontSize':'12px'},
                ),
                html.Span(id='selected-chip', style={
                    'marginLeft':'auto','color':'#3182CE',
                    'fontFamily':'monospace','fontSize':'12px',
                    'fontWeight':'500',
                }, children='— 염색체 클릭'),
            ]),
            html.Div(style={'padding':'12px 16px 8px','overflowX':'auto'}, children=[
                dashbio.Ideogram(
                    id='ideogram',
                    organism='human',
                    assembly=ASSEMBLY,
                    orientation='vertical',
                    chromosomes=DEFAULT_CHROMS,
                    annotations=ANNOTATIONS,
                    annotationsLayout='overlay',
                    annotationHeight=8,
                    showBandLabels=False,
                    showAnnotTooltip=True,
                    rotatable=True,
                    chrHeight=280,
                    chrWidth=14,
                    chrMargin=12,
                    rows=1,
                    legend=LEGEND,
                ),
            ]),
        ]),

        # ── Detail Panel ────────────────────────────────────────
        card([
            card_header([], {'display':'flex'}, ),  # placeholder; overwritten by id
            html.Div(id='detail-card-header', style={
                'padding':'10px 16px',
                'borderBottom': _BD,
                'background':'#F7FAFC',
                'fontSize':'11px','fontWeight':'700','color':'#4A5568',
                'letterSpacing':'.08em','textTransform':'uppercase',
                'display':'flex','alignItems':'center','gap':'8px',
            }, children=['DETAIL VIEW']),
            dcc.Graph(
                id='detail-graph',
                figure=EMPTY_FIG,
                config={'displayModeBar': True,
                        'modeBarButtonsToRemove': ['select2d','lasso2d','autoScale2d'],
                        'responsive': True},
            ),
            # 클릭 힌트
            html.Div(id='click-hint', style={
                'padding':'6px 16px 10px',
                'fontSize':'11px','color':'#A0AEC0','textAlign':'center',
            }),
        ]),

    ]),

    # ── Footer ──────────────────────────────────────────────────
    html.Div(style={
        'background':'white','borderTop':_BD,
        'padding':'8px 24px','fontSize':'11px','color':'#718096',
        'display':'flex','gap':'24px','flexWrap':'wrap',
    }, children=[
        html.Div(['Sample  ', html.B(SAMPLE_ID, style={'color':'#2D3748'})]),
        html.Div(['Sex  ',    html.B('Female (XX)' if SEX=='female' else 'Male (XY)',
                                    style={'color':'#2D3748'})]),
        html.Div(['Assembly  ', html.B(ASSEMBLY, style={'color':'#2D3748'})]),
        html.Div(['Events  ',   html.B(str(len(EVENTS)), style={'color':'#2D3748'})]),
    ]),

    # ── Hidden stores ───────────────────────────────────────────
    dcc.Store(id='current-chrom',  data=None),   # 현재 선택 염색체
    dcc.Store(id='click-count',    data=0),       # 같은 염색체 클릭 횟수
    dcc.Store(id='last-rotated',   data=None),    # 직전 rotated 값
])

# ══════════════════════════════════════════════════════════════════
# 7. Callbacks
# ══════════════════════════════════════════════════════════════════

@callback(
    Output('ideogram', 'chromosomes'),
    Input('chrom-dropdown', 'value'),
)
def update_chromosomes(value):
    return value or DEFAULT_CHROMS


@callback(
    Output('current-chrom', 'data'),
    Output('selected-chip', 'children'),
    Input('ideogram', 'annotationsData'),
    prevent_initial_call=True,
)
def on_hover(ann_data):
    """Annotation hover → chrom store 업데이트."""
    if not ann_data:
        return no_update, no_update

    if isinstance(ann_data, str):
        try:
            ann_data = json.loads(ann_data)
        except Exception:
            return no_update, no_update

    chrom = str(ann_data.get('chromosome') or ann_data.get('chr', '')).replace('chr','')
    name  = ann_data.get('name', '')
    ev    = EVENTS.get(chrom, {})
    cn    = ev.get('cn', '—') if ev else '2'
    chip  = f"chr{chrom}  {name}  CN={cn}  ← 클릭"
    return chrom, chip


@callback(
    Output('detail-graph',       'figure'),
    Output('detail-card-header', 'children'),
    Output('click-hint',         'children'),
    Output('click-count',        'data'),
    Output('last-rotated',       'data'),
    Input('ideogram', 'rotated'),
    State('current-chrom',  'data'),
    State('click-count',    'data'),
    State('last-rotated',   'data'),
    prevent_initial_call=True,
)
def on_click(rotated, chrom, click_count, last_rotated):
    """
    rotated prop은 클릭마다 True/False 토글.
    같은 염색체를 두 번 클릭하면 → 구간 줌인.
    다른 염색체 클릭 → 초기화.
    """
    if chrom is None:
        return no_update, no_update, no_update, no_update, no_update

    # rotated 값이 바뀌었는지 확인 (실제 클릭 발생)
    if rotated == last_rotated:
        return no_update, no_update, no_update, no_update, no_update

    ev = EVENTS.get(chrom)

    # 클릭 횟수 관리
    new_count = (click_count or 0) + 1

    if new_count == 1:
        # ── 1st click: 전체 염색체 CN/BAF ──
        fig    = make_detail_fig(chrom)
        ev_lbl = f"  [{ev['iscn']}  {TYPE_LABEL.get(ev['type'],'')}]" if ev else '  [정상]'
        header = [
            f'DETAIL  —  chr{chrom}' + ev_lbl,
        ]
        if ev and ev['color']:
            header = [
                f'DETAIL  —  chr{chrom}',
                html.Span(f"  {ev['iscn']}", style={
                    'color': ev['color'], 'fontFamily':'monospace',
                    'fontSize':'12px', 'fontWeight':'600',
                }),
                html.Span(f"  {TYPE_LABEL.get(ev['type'],'')}",
                          style={'color':'#718096','fontWeight':'400',
                                 'textTransform':'none','letterSpacing':'0'}),
            ]
        hint = '같은 염색체를 한 번 더 클릭하면 이상 구간을 줌인합니다.' if ev else ''
        return fig, header, hint, new_count, rotated

    else:
        # ── 2nd+ click: 구간 줌인 ──
        fig    = make_region_fig(chrom)
        ev_lbl = f"  [{ev['iscn']}]  구간 상세" if ev else '  전체 보기'
        header = [
            f'REGION DETAIL  —  chr{chrom}',
            html.Span(f"  {ev['iscn']}" if ev else '', style={
                'color': ev['color'] if ev else '#718096',
                'fontFamily':'monospace','fontSize':'12px','fontWeight':'600',
            }),
            html.Span('  구간 상세',
                      style={'color':'#718096','fontWeight':'400',
                             'textTransform':'none','letterSpacing':'0'}),
        ]
        hint = '다른 염색체를 클릭하면 돌아옵니다.'
        # 3번째 클릭 이상이면 count 리셋 → 다시 1st click 동작
        next_count = 0 if new_count >= 2 else new_count
        return fig, header, hint, next_count, rotated


if __name__ == '__main__':
    app.run(debug=True, port=8050)