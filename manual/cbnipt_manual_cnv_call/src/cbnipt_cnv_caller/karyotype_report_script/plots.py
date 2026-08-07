import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.io import to_html
from config import CHROM_SIZES, DEFAULT_COLORS

def _hex_to_rgba(c, a=0.09):
    h = c.lstrip('#')
    return f'rgba({int(h[0:2],16)},{int(h[2:4],16)},{int(h[4:6],16)},{a})'

_BG    = '#FFFFFF'
_PAPER = '#F8FAFC'
_GRID  = '#E2E8F0'
_TICK  = '#64748B'
_AXIS  = dict(gridcolor=_GRID, zeroline=False, linecolor=_GRID,
              tickfont=dict(size=10, color=_TICK),
              title_font=dict(size=11, color=_TICK))


def _pt_colors(sub, cnv_sub, chrom):
    """각 데이터 포인트를 해당 CNV 구간 색으로 염색."""
    colors = np.full(len(sub), '#94A3B8', dtype=object)
    for _, seg in cnv_sub.iterrows():
        mask = (sub['pos'] >= seg['start']) & (sub['pos'] <= seg['end'])
        c = seg.get('color')
        if not c or str(c).lower() == 'nan':
            from annotations import infer_event_type
            c = DEFAULT_COLORS[infer_event_type(seg)]
        colors[mask] = c
    return colors.tolist()


def chromosome_figure(signal, chrom, cnv):
    sub     = signal[signal['chrom'] == chrom].copy()
    cnv_sub = cnv[cnv['chrom'] == chrom]
    has_baf = 'baf' in sub.columns and sub['baf'].notna().any()
    nrows   = 2 if has_baf else 1

    titles = [f'Copy Number — chr{chrom}']
    if has_baf:
        titles.append(f'BAF — chr{chrom}')

    fig = make_subplots(
        rows=nrows, cols=1, shared_xaxes=True,
        vertical_spacing=0.10,
        row_heights=[0.55, 0.45] if has_baf else [1.0],
        subplot_titles=titles,
    )

    if sub.empty:
        fig.add_annotation(text=f'chr{chrom}: no signal data',
                           x=0.5, y=0.5, xref='paper', yref='paper', showarrow=False)
    else:
        xmb    = sub['pos'] / 1e6
        span   = xmb.max() - xmb.min()
        msize  = 4.5 if span <= 30 else 3.0
        colors = _pt_colors(sub, cnv_sub, chrom)

        # CN scatter
        fig.add_trace(go.Scattergl(
            x=xmb, y=sub['cn'], mode='markers',
            marker=dict(size=msize, color=colors, opacity=0.82),
            hovertemplate='%{x:.3f} Mb<br>CN: %{y:.3f}<extra></extra>',
            name='CN',
        ), row=1, col=1)

        # BAF scatter
        if has_baf:
            fig.add_trace(go.Scattergl(
                x=xmb, y=sub['baf'], mode='markers',
                marker=dict(size=msize, color=colors, opacity=0.70),
                hovertemplate='%{x:.3f} Mb<br>BAF: %{y:.3f}<extra></extra>',
                name='BAF',
            ), row=2, col=1)

    # CNV 구간 shade
    for _, seg in cnv_sub.iterrows():
        x0, x1 = seg['start'] / 1e6, seg['end'] / 1e6
        c = seg.get('color')
        if not c or str(c).lower() == 'nan':
            from annotations import infer_event_type
            c = DEFAULT_COLORS[infer_event_type(seg)]
        for row in range(1, nrows + 1):
            fig.add_vrect(x0=x0, x1=x1,
                          fillcolor=_hex_to_rgba(c, 0.09), line_color=c,
                          line_width=1.2, row=row, col=1)

    # 기준선
    for v, op in [(1, 0.4), (2, 0.65), (3, 0.4)]:
        fig.add_hline(y=v, line=dict(color='#64748B', width=1, dash='dot'),
                      opacity=op, row=1, col=1)
    if has_baf:
        for v in [0.33, 0.5, 0.67]:
            fig.add_hline(y=v, line=dict(color='#64748B', width=1, dash='dot'),
                          opacity=0.4, row=2, col=1)

    fig.update_yaxes(title_text='CN', range=[0, 5.5], row=1, col=1, **_AXIS)
    if has_baf:
        fig.update_yaxes(
            title_text='BAF', range=[-0.05, 1.05],
            tickvals=[0, 0.25, 0.33, 0.5, 0.67, 0.75, 1],
            ticktext=['0', '.25', '⅓', '.5', '⅔', '.75', '1'],
            row=2, col=1, **_AXIS,
        )
    fig.update_xaxes(title_text='Position (Mb)', row=nrows, col=1, **_AXIS)
    fig.update_annotations(font=dict(size=12, color='#334155'))

    fig.update_layout(
        height=380 if has_baf else 260,
        margin=dict(l=58, r=30, t=48, b=40),
        paper_bgcolor=_PAPER, plot_bgcolor=_BG,
        showlegend=False, hovermode='x unified',
        font=dict(family='Inter, system-ui, sans-serif', size=11),
        hoverlabel=dict(bgcolor='white', bordercolor=_GRID,
                        font=dict(size=11, color='#1E293B')),
    )
    return fig


def build_plot_html(signal, cnv, include_plotlyjs=True):
    order  = {str(i): i for i in range(1, 23)}
    order.update({'X': 23, 'Y': 24})
    chroms = sorted(signal['chrom'].unique(),
                    key=lambda c: order.get(c.upper(), 99))
    out   = {}
    first = True
    for chrom in chroms:
        fig = chromosome_figure(signal, chrom, cnv)
        out[chrom] = to_html(
            fig, full_html=False,
            include_plotlyjs=('cdn' if include_plotlyjs and first else False),
            config={'responsive': True, 'displaylogo': False,
                    'modeBarButtonsToRemove': ['select2d', 'lasso2d']},
        )
        first = False
    return out
