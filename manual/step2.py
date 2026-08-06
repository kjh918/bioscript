"""
STEP 2 — CN 이상 + Gene + p/q Band annotation 실제 동작 버전

핵심 발견:
  dash_bio Ideogram의 annotations prop에 단순 list를 주면 drawAnnots()에
  원시배열이 들어가 렌더가 안 됨.
  실제로는 rawAnnots 형식 {keys, annots} Object를 줘야
  afterRawAnnots() 경로로 처리되어 정상 렌더됨.

  rawAnnots 형식:
    {
      "keys": ["name", "start", "length", "trackIndex", "color", "shape"],
      "annots": [{"chr": "21", "annots": [[name, start, length, trackIndex, color], ...]}, ...]
    }

  annotationTracks로 트랙 레이블/색상 지정,
  trackIndex 0 = CN 이상, trackIndex 1 = Gene

pip install dash dash-bio
python step2_annotations.py → http://localhost:8051
"""

import json
from dash import Dash, html, dcc, Input, Output, callback
import dash_bio as dashbio


def parse_bp(value, default=0):
    """Ideogram callback의 '7,275,128' 같은 좌표 문자열을 안전하게 정수로 변환."""
    if value is None or value == '':
        return default
    if isinstance(value, (int, float)):
        return int(value)
    text = str(value).strip().replace(',', '').replace('bp', '').strip()
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return default

# ─── 염색체 사이즈 (GRCh38) ───────────────────────────────────────
CHROM_SIZES = {
    '1':248956422,'2':242193529,'3':198295559,'4':190214555,
    '5':181538259,'6':170805979,'7':159345973,'8':145138636,
    '9':138394717,'10':133797422,'11':135086622,'12':133275309,
    '13':114364328,'14':107043718,'15':101991189,'16':90338345,
    '17':83257441,'18':80373285,'19':58617616,'20':64444167,
    '21':46709983,'22':50818468,'X':156040895,'Y':57227415,
}

# GRCh38 centromere 위치 (Mb) — p/q arm 판별용
CENTROMERE_MB = {
    '1':123.4,'2':93.9,'3':90.9,'4':50.2,'5':48.8,
    '6':61.0,'7':59.9,'8':45.2,'9':43.0,'10':39.8,
    '11':53.4,'12':35.5,'13':17.7,'14':17.2,'15':19.0,
    '16':36.8,'17':25.1,'18':18.5,'19':26.2,'20':28.1,
    '21':12.0,'22':15.0,'X':61.0,'Y':10.4,
}

# ─── CN 이상 이벤트 (trackIndex=0) ───────────────────────────────
CN_EVENTS = [
    {'chr':'21','name':'+21 Trisomy', 'start':1,        'stop':46709983,  'color':'#E53E3E'},
    {'chr':'X', 'name':'-X Monosomy', 'start':1,        'stop':156040895, 'color':'#3182CE'},
    {'chr':'5', 'name':'del(5p)',     'start':1,        'stop':45368512,  'color':'#3182CE'},
    {'chr':'17','name':'dup(17q)',    'start':43044000, 'stop':83257441,  'color':'#DD6B20'},
]

# ─── Gene annotation (trackIndex=1) ──────────────────────────────
GENE_ANNOTS = [
    {'chr':'21','name':'DYRK1A', 'start':37700000, 'stop':37865335, 'color':'#6B46C1'},
    {'chr':'21','name':'RUNX1',  'start':34787801, 'stop':36004954, 'color':'#6B46C1'},
    {'chr':'17','name':'TP53',   'start':7661779,  'stop':7687538,  'color':'#276749'},
    {'chr':'17','name':'BRCA1',  'start':43044295, 'stop':43125483, 'color':'#276749'},
    {'chr':'17','name':'ERBB2',  'start':39687914, 'stop':39730426, 'color':'#276749'},
    {'chr':'5', 'name':'CTNND2', 'start':11080000, 'stop':11690000, 'color':'#276749'},
    {'chr':'X', 'name':'MECP2',  'start':154021573,'stop':154137103,'color':'#276749'},
    {'chr':'X', 'name':'DMD',    'start':31119221, 'stop':33339609, 'color':'#276749'},
    {'chr':'13','name':'BRCA2',  'start':32315086, 'stop':32400266, 'color':'#276749'},
    {'chr':'7', 'name':'CFTR',   'start':117480025,'stop':117668665,'color':'#276749'},
]

def make_raw_annots(cn_events, gene_annots):
    """
    dash_bio Ideogram rawAnnots 형식으로 변환.
    keys: [name, start, length, trackIndex, color, shape]
    annots: [{chr, annots: [[...], ...]}]
    trackIndex 0=CN이상, 1=Gene
    """
    all_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']
    chrom_dict = {c: [] for c in all_chroms}

    for ev in cn_events:
        c = ev['chr']
        chrom_dict[c].append([
            ev['name'],
            max(0, ev['start'] - 1),   # Ideogram 좌표: 0-based start
            max(1, ev['stop'] - ev['start'] + 1),  # 1-based inclusive -> length
            0,                          # trackIndex
            ev['color'],
            'rectangle',                 # amplification/deletion을 구간형으로 표시
        ])

    for g in gene_annots:
        c = g['chr']
        chrom_dict[c].append([
            g['name'],
            max(0, g['start'] - 1),    # Ideogram 좌표: 0-based start
            max(1, g['stop'] - g['start'] + 1),
            1,
            g['color'],
            'rectangle',                 # gene도 점이 아닌 genomic interval로 표시
        ])

    annots = [
        {'chr': c, 'annots': chrom_dict[c]}
        for c in all_chroms
        if chrom_dict[c]
    ]
    return {
        'keys': ['name', 'start', 'length', 'trackIndex', 'color', 'shape'],
        'annots': annots,
    }

RAW_ANNOTS = make_raw_annots(CN_EVENTS, GENE_ANNOTS)

ANNOTATION_TRACKS = [
    {'id': 'cn',   'displayName': 'Amplification / Deletion', 'color': '#E53E3E', 'shape': 'rectangle'},
    {'id': 'gene', 'displayName': 'Gene interval',             'color': '#6B46C1', 'shape': 'rectangle'},
]

LEGEND = [{
    'name': 'Track',
    'rows': [
        {'name': 'Trisomy/Gain',  'color': '#E53E3E', 'shape': 'rectangle'},
        {'name': 'Monosomy/Loss', 'color': '#3182CE', 'shape': 'rectangle'},
        {'name': 'Partial Gain',  'color': '#DD6B20', 'shape': 'rectangle'},
        {'name': 'Partial Loss',  'color': '#3182CE', 'shape': 'rectangle'},
        {'name': 'Gene interval', 'color': '#6B46C1', 'shape': 'rectangle'},
    ],
}]

DEFAULT_CHROMS = [str(i) for i in range(1, 23)] + ['X']

# ─────────────────────────────────────────────────────────────────
app = Dash(__name__, title='Karyotype Step 2')

_BD = '1px solid #E2E8F0'

app.layout = html.Div(style={
    'fontFamily': 'Inter, system-ui, sans-serif',
    'background': '#EDF2F7', 'minHeight': '100vh', 'color': '#2D3748',
}, children=[

    # ── 헤더 ──────────────────────────────────────────────────────
    html.Div(style={
        'background': 'white', 'borderBottom': _BD,
        'padding': '14px 24px', 'boxShadow': '0 1px 3px rgba(0,0,0,0.06)',
    }, children=[
        html.H2('Karyotype Viewer — Step 2', style={
            'margin': 0, 'fontSize': '16px', 'fontWeight': '700',
        }),
        html.Div('CN 이상 + Gene + Band label 표시', style={
            'fontSize': '12px', 'color': '#718096', 'marginTop': '2px',
        }),
    ]),

    # ── 컨트롤 ────────────────────────────────────────────────────
    html.Div(style={
        'background': 'white', 'borderBottom': _BD,
        'padding': '12px 24px', 'display': 'flex',
        'gap': '32px', 'alignItems': 'flex-start', 'flexWrap': 'wrap',
    }, children=[

        html.Div([
            html.Label('표시 염색체', style={'fontSize':'12px','fontWeight':'600','color':'#4A5568','display':'block','marginBottom':'4px'}),
            dcc.Dropdown(
                id='chrom-dd',
                options=[{'label': str(i), 'value': str(i)} for i in range(1,23)]
                      + [{'label':'X','value':'X'},{'label':'Y','value':'Y'}],
                multi=True, value=DEFAULT_CHROMS,
                style={'width':'380px','fontSize':'12px'},
            ),
        ]),

        html.Div([
            html.Label('Annotation 레이아웃', style={'fontSize':'12px','fontWeight':'600','color':'#4A5568','display':'block','marginBottom':'4px'}),
            dcc.RadioItems(
                id='layout-radio',
                options=[
                    {'label': ' interval tracks (Gene / CN 구간)', 'value': 'tracks'},
                    {'label': ' chromosome overlay (염색체 위 구간)', 'value': 'overlay'},
                    {'label': ' density histogram (위치별 annotation 수)', 'value': 'histogram'},
                ],
                value='tracks',
                labelStyle={'display':'block','marginBottom':'3px','fontSize':'13px'},
            ),
        ]),

        html.Div([
            html.Label('Histogram 설정', style={'fontSize':'12px','fontWeight':'600','color':'#4A5568','display':'block','marginBottom':'4px'}),
            dcc.RadioItems(
                id='hist-scaling',
                options=[
                    {'label': ' chromosome별 상대 높이', 'value': 'relative'},
                    {'label': ' 전체 chromosome 절대 높이', 'value': 'absolute'},
                ],
                value='relative',
                labelStyle={'display':'block','marginBottom':'3px','fontSize':'13px'},
            ),
            html.Label('Bar width', style={'fontSize':'11px','color':'#718096','display':'block','marginTop':'6px'}),
            dcc.Slider(id='hist-bar-width', min=2, max=14, step=1, value=6,
                       marks={2:'2',6:'6',10:'10',14:'14'}, tooltip={'placement':'bottom'}),
        ], style={'width':'240px'}),

        html.Div([
            html.Label('Brush 영역 선택', style={'fontSize':'12px','fontWeight':'600','color':'#4A5568','display':'block','marginBottom':'4px'}),
            dcc.Checklist(
                id='brush-enable',
                options=[{'label':' 단일 chromosome Brush 활성화', 'value':'on'}],
                value=[],
                labelStyle={'fontSize':'13px'},
            ),
            dcc.Dropdown(
                id='brush-chrom', clearable=False, value='17',
                options=[{'label':f'chr{i}','value':str(i)} for i in range(1,23)]
                      + [{'label':'chrX','value':'X'},{'label':'chrY','value':'Y'}],
                style={'width':'150px','fontSize':'12px','marginTop':'6px'},
            ),
            html.Div(style={'display':'flex','gap':'6px','marginTop':'6px'}, children=[
                dcc.Input(id='brush-start-mb', type='number', value=35, min=0, step=0.1,
                          placeholder='Start Mb', style={'width':'90px','fontSize':'12px'}),
                dcc.Input(id='brush-end-mb', type='number', value=50, min=0, step=0.1,
                          placeholder='End Mb', style={'width':'90px','fontSize':'12px'}),
            ]),
        ]),

        html.Div([
            html.Label('표시 옵션', style={'fontSize':'12px','fontWeight':'600','color':'#4A5568','display':'block','marginBottom':'4px'}),
            dcc.Checklist(
                id='options-check',
                options=[
                    {'label': ' Band label (p/q 위치)', 'value': 'band'},
                    {'label': ' Annotation tooltip',   'value': 'tooltip'},
                ],
                value=['band', 'tooltip'],
                labelStyle={'display':'block','marginBottom':'3px','fontSize':'13px'},
            ),
        ]),

    ]),

    # ── Ideogram ──────────────────────────────────────────────────
    html.Div(style={
        'background': 'white', 'margin': '16px 24px 0',
        'border': _BD, 'borderRadius': '8px', 'overflow': 'hidden',
    }, children=[
        html.Div('IDEOGRAM', style={
            'padding': '8px 16px', 'borderBottom': _BD,
            'fontSize': '10px', 'fontWeight': '700', 'color': '#718096',
            'letterSpacing': '.1em', 'background': '#F7FAFC',
        }),
        html.Div(style={'padding': '12px 16px', 'overflowX': 'auto'}, children=[
            dashbio.Ideogram(
                id='ideo',
                organism='human',
                assembly='GRCh38',
                orientation='vertical',
                chrHeight=300,
                chrWidth=12,
                chrMargin=14,
                rows=1,
                rotatable=True,
                showBandLabels=True,
                showChromosomeLabels=True,
                resolution=550,
                # ← rawAnnots 형식으로 전달
                annotations=RAW_ANNOTS,
                annotationsLayout='tracks',
                annotationTracks=ANNOTATION_TRACKS,
                annotationHeight=8,
                barWidth=6,
                histogramScaling='relative',
                brush=None,
                showAnnotTooltip=True,
                legend=LEGEND,
            ),
        ]),
    ]),

    # ── Hover 정보 카드 ───────────────────────────────────────────
    html.Div(style={
        'margin': '12px 24px 0',
        'background': 'white', 'border': _BD, 'borderRadius': '8px',
        'overflow': 'hidden',
    }, children=[
        html.Div('ANNOTATION INFO', style={
            'padding': '8px 16px', 'borderBottom': _BD,
            'fontSize': '10px', 'fontWeight': '700', 'color': '#718096',
            'letterSpacing': '.1em', 'background': '#F7FAFC',
        }),
        html.Div(id='annot-info', style={'padding': '14px 18px', 'minHeight': '52px'}),
    ]),

    # ── Brush 선택 정보 ─────────────────────────────────────────
    html.Div(style={
        'margin': '12px 24px 0',
        'background': 'white', 'border': _BD, 'borderRadius': '8px',
        'overflow': 'hidden',
    }, children=[
        html.Div('BRUSH SELECTION', style={
            'padding': '8px 16px', 'borderBottom': _BD,
            'fontSize': '10px', 'fontWeight': '700', 'color': '#718096',
            'letterSpacing': '.1em', 'background': '#F7FAFC',
        }),
        html.Div(id='brush-info', children='Brush를 활성화하면 선택 범위가 표시됩니다.',
                 style={'padding':'12px 16px','fontSize':'13px','color':'#718096'}),
    ]),

    # ── p/q Arm 테이블 ───────────────────────────────────────────
    html.Div(style={
        'margin': '12px 24px 20px',
        'background': 'white', 'border': _BD, 'borderRadius': '8px',
        'overflow': 'hidden',
    }, children=[
        html.Div('ANNOTATION SUMMARY — p/q Arm', style={
            'padding': '8px 16px', 'borderBottom': _BD,
            'fontSize': '10px', 'fontWeight': '700', 'color': '#718096',
            'letterSpacing': '.1em', 'background': '#F7FAFC',
        }),
        html.Div(id='arm-table', style={'padding': '12px 16px'}),
    ]),

])


# ─── Callbacks ────────────────────────────────────────────────────

@callback(
    Output('ideo','chromosomes'),
    Output('ideo','orientation'),
    Output('ideo','brush'),
    Output('ideo','chrHeight'),
    Input('chrom-dd','value'),
    Input('brush-enable','value'),
    Input('brush-chrom','value'),
    Input('brush-start-mb','value'),
    Input('brush-end-mb','value'),
)
def upd_view(chroms, brush_enable, brush_chrom, start_mb, end_mb):
    """Brush는 한 chromosome에서 가장 안정적으로 동작하므로 자동 단일 보기로 전환."""
    if 'on' not in (brush_enable or []):
        return chroms or DEFAULT_CHROMS, 'vertical', None, 300

    chrom = str(brush_chrom or '17')
    chrom_size = CHROM_SIZES[chrom]
    start = max(1, int(float(start_mb or 0) * 1_000_000))
    end = min(chrom_size, int(float(end_mb or chrom_size / 1e6) * 1_000_000))
    if end <= start:
        end = min(chrom_size, start + 1_000_000)
    return [chrom], 'horizontal', f'chr{chrom}:{start}-{end}', 850


@callback(
    Output('ideo','annotationsLayout'),
    Output('ideo','barWidth'),
    Output('ideo','histogramScaling'),
    Input('layout-radio','value'),
    Input('hist-bar-width','value'),
    Input('hist-scaling','value'),
)
def upd_layout(layout, bar_width, scaling):
    return layout, int(bar_width or 6), scaling or 'relative'

@callback(
    Output('ideo','showBandLabels'),
    Output('ideo','showAnnotTooltip'),
    Input('options-check','value'),
)
def upd_options(v):
    v = v or []
    return 'band' in v, 'tooltip' in v


@callback(Output('brush-info','children'), Input('ideo','brushData'))
def show_brush(data):
    if not data:
        return html.Span('Brush를 드래그하면 선택 구간이 표시됩니다.',
                         style={'color':'#A0AEC0'})
    if isinstance(data, str):
        try:
            data = json.loads(data)
        except Exception:
            return html.Pre(str(data), style={'fontSize':'11px'})

    start = parse_bp(data.get('start'), 0)
    end = parse_bp(data.get('end'), 0)
    extent = parse_bp(data.get('extent'), max(0, end - start))
    return html.Div(style={'display':'flex','gap':'28px','flexWrap':'wrap'}, children=[
        _kv('Start', f'{start:,} bp ({start/1e6:.3f} Mb)'),
        _kv('End', f'{end:,} bp ({end/1e6:.3f} Mb)'),
        _kv('Extent', f'{extent:,} bp ({extent/1e6:.3f} Mb)'),
    ])


@callback(Output('annot-info','children'), Input('ideo','annotationsData'))
def show_hover(data):
    if not data:
        return html.Span('← Annotation에 마우스를 올리면 정보가 표시됩니다',
                         style={'color':'#A0AEC0','fontSize':'13px'})

    if isinstance(data, str):
        try:
            data = json.loads(data)
        except Exception:
            return html.Pre(str(data), style={'fontSize':'11px'})

    name  = data.get('name','—')
    chrom = str(data.get('chromosome') or data.get('chr','')).replace('chr','')
    # annotationsData도 버전에 따라 쉼표 문자열을 반환할 수 있음
    start0 = parse_bp(data.get('start'), 0)
    stop0 = parse_bp(data.get('stop'), 0)
    length = parse_bp(data.get('length'), 0)
    if stop0 <= start0 and length > 0:
        stop0 = start0 + length

    # 표시 좌표는 사용자가 입력한 1-based 좌표로 되돌림
    start = start0 + 1 if start0 >= 0 else start0
    stop = stop0
    s_mb = f"{start/1e6:.3f} Mb"
    e_mb = f"{stop/1e6:.3f} Mb"
    len_mb = f"{max(0, stop-start+1)/1e6:.3f} Mb"

    # p/q arm 판별
    arm = '—'
    try:
        mid_mb = (start + stop) / 2 / 1e6
        cent   = CENTROMERE_MB.get(chrom, 50)
        arm    = 'p arm' if mid_mb < cent else 'q arm'
    except Exception:
        pass

    # track 구분
    track_idx = data.get('trackIndex', data.get('index', None))
    track_lbl = {0:'CN Anomaly', 1:'Gene'}.get(track_idx, '—')
    track_color = {'CN Anomaly':'#E53E3E','Gene':'#6B46C1'}.get(track_lbl,'#718096')

    return html.Div([
        html.Div(style={'display':'flex','alignItems':'center','gap':'8px','marginBottom':'8px'}, children=[
            html.Span(name, style={'fontWeight':'700','fontSize':'15px'}),
            html.Span(f'chr{chrom}', style={
                'background':'#EBF8FF','color':'#2B6CB0',
                'padding':'2px 8px','borderRadius':'4px','fontSize':'11px','fontWeight':'600',
            }),
            html.Span(track_lbl, style={
                'background': track_color+'18', 'color': track_color,
                'padding':'2px 8px','borderRadius':'4px','fontSize':'11px','fontWeight':'600',
            }),
        ]),
        html.Div(style={'display':'flex','gap':'28px','fontSize':'12px','color':'#4A5568'}, children=[
            _kv('Start',  s_mb),
            _kv('End',    e_mb),
            _kv('Length', len_mb),
            _kv('Arm',    arm),
            _kv('Chr',    f'chr{chrom}'),
        ]),
    ])


def _kv(k, v):
    return html.Div([
        html.Span(k+': ', style={'color':'#718096','fontWeight':'600'}),
        html.Span(v,      style={'color':'#1A202C'}),
    ])


@callback(Output('arm-table','children'), Input('ideo','chromosomes'))
def build_arm_table(_):
    """모든 annotation을 p/q arm 기준으로 정리한 테이블."""
    all_events = (
        [{'track':'CN', **e} for e in CN_EVENTS] +
        [{'track':'Gene', **g} for g in GENE_ANNOTS]
    )

    rows = []
    for ev in all_events:
        chrom = ev['chr']
        mid_mb = (ev['start'] + ev['stop']) / 2 / 1e6
        cent   = CENTROMERE_MB.get(chrom, 50)
        arm    = 'p' if mid_mb < cent else 'q'
        band_approx = f"{chrom}{arm}{mid_mb:.0f}"

        color = ev['color']
        rows.append(html.Tr([
            html.Td(f"chr{chrom}", style=td_s()),
            html.Td(html.Span(arm + ' arm', style={
                'background': '#EBF8FF' if arm=='p' else '#F0FFF4',
                'color': '#2B6CB0' if arm=='p' else '#276749',
                'padding':'1px 8px','borderRadius':'4px',
                'fontSize':'11px','fontWeight':'600',
            })),
            html.Td(ev['name'], style=td_s()),
            html.Td(html.Span(ev['track'], style={
                'background': color+'22', 'color': color,
                'padding':'1px 8px','borderRadius':'4px',
                'fontSize':'11px','fontWeight':'600',
            })),
            html.Td(f"{ev['start']/1e6:.2f} Mb", style=td_s()),
            html.Td(f"{ev['stop']/1e6:.2f} Mb",  style=td_s()),
            html.Td(f"{(ev['stop']-ev['start'])/1e6:.3f} Mb", style=td_s()),
        ], style={'borderBottom': _BD}))

    header = html.Tr([
        html.Th(h, style={
            'padding':'8px 12px','textAlign':'left',
            'fontSize':'10px','fontWeight':'700',
            'color':'#718096','letterSpacing':'.06em',
            'textTransform':'uppercase','background':'#F7FAFC',
            'borderBottom': _BD,
        })
        for h in ['Chr','Arm','Name','Track','Start','End','Length']
    ])

    return html.Table(
        [html.Thead(header), html.Tbody(rows)],
        style={'width':'100%','borderCollapse':'collapse','fontSize':'13px'},
    )


def td_s():
    return {'padding':'8px 12px','color':'#2D3748'}


if __name__ == '__main__':
    app.run(debug=True, port=8051)