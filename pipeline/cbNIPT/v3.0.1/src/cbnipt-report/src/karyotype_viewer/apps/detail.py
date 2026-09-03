"""
Dash App 2 – Chromosome Detail (report-style layout).

레이아웃: reporter.py의 cbNIPT report 스타일을 Dash로 재현.
  - 상단: 종합 판정 배너 + sample info
  - 좌: Ideogram + brush + syndrome 카드
  - 우: CN/BAF graph + syndrome summary chart
"""

from __future__ import annotations
import json
from typing import Optional

import pandas as pd

from dash import Dash, html, dcc, Input, Output, State, no_update, ctx
import dash_bio as dashbio

from ..core.models import SampleData
from ..core.annot import build_iscn, annotation_options_for_chrom
from ..core.data import demo_dataframe
from ..core.reference import CHROM_SIZES
from ..core.cnv_loader import get_chrom_df
from ..core.nipt_markers import NiptSyndrome
from ..ui.figures import region_fig, syndrome_summary_fig

# ── Report colour tokens (CSS vars via inline style) ──────────────────────
_CSS = """
<style>
:root {
  --bg:#F4F6F9; --surface:#FFFFFF; --surface2:#F8F9FB;
  --border:#E1E5EC; --border-s:#BEC5D1;
  --text:#0E1520; --text-sub:#48536A; --text-muted:#8A94A8;
  --navy:#162040; --navy-l:#EAF0FB; --navy-b:#A8BDD8;
  --teal:#0A6E5C; --teal-l:#E8F7F4; --teal-b:#6EC4B5;
  --red:#9B1B1B;  --red-l:#FDF0F0;  --red-b:#EFA5A5;
  --amber:#7A4800;--amber-l:#FEF8EC;--amber-b:#F0C060;
  --mono:'SF Mono','Fira Code','Consolas',monospace;
  --sans:-apple-system,'Segoe UI','Apple SD Gothic Neo','Malgun Gothic',sans-serif;
  --r:10px;
}
* { box-sizing: border-box; }
body {
  font-family: var(--sans);
  background: var(--bg);
  color: var(--text);
  font-size: 13px;
  line-height: 1.55;
  margin: 0;
}
.rpt-page {
  max-width: 1100px;
  margin: 0 auto;
  padding: 12px 14px 40px;
}

/* Banner */
.rbanner { border-radius: var(--r); padding: .9rem 1.2rem; margin-bottom: 1rem;
           display: flex; align-items: center; gap: 1rem; flex-wrap: wrap; }
.rbanner.lr { background: var(--teal-l); border: 1.5px solid var(--teal-b); }
.rbanner.hr { background: var(--red-l);  border: 1.5px solid var(--red-b); }
.rbanner.mo { background: var(--amber-l);border: 1.5px solid var(--amber-b); }
.ricon { width: 38px; height: 38px; border-radius: 50%;
         display: flex; align-items: center; justify-content: center;
         font-size: 16px; flex-shrink: 0; }
.lr .ricon { background: var(--teal-b); }
.hr .ricon { background: var(--red-b); }
.mo .ricon { background: var(--amber-b); }
.rlbl { font-size: 10px; font-weight: 600; color: var(--text-muted); }
.rres { font-size: 16px; font-weight: 700; }
.lr .rres { color: var(--teal); } .hr .rres { color: var(--red); } .mo .rres { color: var(--amber); }
.rdesc { font-size: 12px; color: var(--text-sub); margin-top: 2px; }

/* Cards */
.rcard { background: var(--surface); border: 1px solid var(--border);
         border-radius: var(--r); margin-bottom: 1rem; overflow: hidden; }
.rch { display: flex; align-items: center; justify-content: space-between;
       padding: 8px 14px; border-bottom: 1px solid var(--border);
       background: var(--surface2); }
.rct { font-size: 10px; font-weight: 700; letter-spacing: .07em;
       text-transform: uppercase; color: var(--text-muted); }
.rcb { padding: 12px 14px; }

/* Info grid */
.igrid { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
.iitem { padding: 8px 12px; border-right: 1px solid var(--border); border-bottom: 1px solid var(--border); }
.iitem .il { font-size: 10px; color: var(--text-muted); margin-bottom: 1px; }
.iitem .iv { font-size: 13px; font-weight: 500; }

/* Syndrome call pills */
.pill { display: inline-block; font-size: 11px; font-weight: 700;
        padding: 3px 9px; border-radius: 99px; border: 1px solid; white-space: nowrap; }
.plr { background: var(--teal-l); color: var(--teal); border-color: var(--teal-b); }
.phr { background: var(--red-l);  color: var(--red);  border-color: var(--red-b); }
.pmo { background: var(--amber-l);color: var(--amber); border-color: var(--amber-b); }
.pnc { background: var(--surface2); color: var(--text-muted); border-color: var(--border); }

/* Syndrome table */
.stbl { width: 100%; border-collapse: collapse; font-size: 12px; }
.stbl thead th { background: var(--surface2); padding: 7px 10px;
  font-size: 10px; font-weight: 700; text-transform: uppercase;
  letter-spacing: .05em; color: var(--text-muted);
  border-bottom: 1px solid var(--border-s); text-align: left; }
.stbl tbody td { padding: 8px 10px; border-bottom: 1px solid var(--border); vertical-align: middle; }
.stbl tbody tr:last-child td { border-bottom: none; }
.catrow td { background: var(--surface2) !important; font-size: 9px; font-weight: 700;
  letter-spacing: .07em; text-transform: uppercase; color: var(--text-muted);
  padding: 5px 10px !important; border-top: 1px solid var(--border) !important; }
.cn-name { font-weight: 600; font-size: 12px; }
.cn-sub  { font-size: 10px; color: var(--text-muted); }

/* Ideogram brush chip */
.brush-chip { font-family: var(--mono); font-size: 11px;
              color: #3182CE; font-weight: 600; }
.brush-hint { font-size: 10px; color: var(--text-muted); }

/* Dropdown override */
.Select-control { font-size: 11px !important; }

/* Two-column layout */
.rpt-cols { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }
@media (max-width: 800px) { .rpt-cols { grid-template-columns: 1fr; } }
</style>
"""

# ---------------------------------------------------------------------------
# postMessage JS
# ---------------------------------------------------------------------------
_DETAIL_MESSAGE_JS = r"""
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
    if ((attempt || 0) < 50)
      setTimeout(function () { setChrom(chrom, (attempt||0)+1); }, 100);
  }
  window.addEventListener('message', function (event) {
    if (event.origin !== window.location.origin) return;
    const d = event.data || {};
    if (d.type !== 'karyotype-chromosome-selected') return;
    setChrom(d.chrom, 0);
  });
})();
</script>
"""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_PILL_MAP = {
    "ABNORMAL":  ("phr", "ABNORMAL"),
    "SUSPICIOUS":("pmo", "SUSPICIOUS"),
    "NORMAL":    ("plr", "NORMAL"),
    "UNKNOWN":   ("pnc", "UNKNOWN"),
}
_GROUP_ORDER = ["Autosome Abnormality", "Sex Chromosome Abnormality", "Micro Deletion"]
_GROUP_COLOR = {
    "Autosome Abnormality":       "#FC8181",
    "Sex Chromosome Abnormality": "#90CDF4",
    "Micro Deletion":             "#F6AD55",
}

def _parse_bp(value, default: int = 0) -> int:
    if value is None or value == "": return default
    if isinstance(value, (int, float)): return int(value)
    try: return int(float(str(value).replace(",","").replace("bp","").strip()))
    except: return default


def _overall_call(syndromes: dict) -> tuple[str, str, str, str]:
    """(css_class, icon, title, desc)"""
    calls = [s.call for s in syndromes.values()]
    abn  = [s.syndrome for s in syndromes.values() if s.call == "ABNORMAL"]
    susp = [s.syndrome for s in syndromes.values() if s.call == "SUSPICIOUS"]
    if "ABNORMAL" in calls:
        desc = f"{', '.join(abn)} — 이상 소견 확인."
        if susp: desc += f" 추가 확인 필요: {', '.join(susp)}."
        return "hr", "⚠", "ABNORMAL — 이상 소견 확인", desc
    if "SUSPICIOUS" in calls:
        return "mo", "〜", "SUSPICIOUS — 추가 확인 필요", f"{', '.join(susp)}"
    return "lr", "✓", "NORMAL — 이상 소견 없음", "검사한 전 항목에서 이상 소견이 관찰되지 않았습니다."


def _syndrome_table(syndromes: dict) -> html.Table:
    rows = []
    for grp in _GROUP_ORDER:
        items = [s for s in syndromes.values() if s.group == grp]
        if not items: continue
        gc = _GROUP_COLOR.get(grp, "#CBD5E0")
        rows.append(html.Tr(className="catrow", children=[
            html.Td(grp, colSpan=4,
                    style={"borderLeft": f"3px solid {gc}"})
        ]))
        for s in sorted(items, key=lambda x: x.syndrome):
            pill_cls, pill_txt = _PILL_MAP.get(s.call, ("pnc", s.call))
            cn_str = f"{s.cn_value:.3f}" if s.cn_value is not None else "—"
            primary = next(
                (f for f in s.features
                 if f.feature_type in ("TargetChromosome","PrimaryTargetRegion")),
                s.features[0] if s.features else None,
            )
            region_str = (
                f"chr{primary.chromosome.replace('chr','')} "
                f"{primary.start/1e6:.1f}–{primary.end/1e6:.1f} Mb"
            ) if primary else "—"
            rows.append(html.Tr([
                html.Td([
                    html.Div(s.syndrome, className="cn-name"),
                    html.Div(s.nipt_id,  className="cn-sub"),
                ]),
                html.Td(html.Span(pill_txt, className=f"pill {pill_cls}")),
                html.Td(cn_str, style={"fontFamily":"monospace","fontSize":"11px"}),
                html.Td(region_str, style={"fontSize":"10px","color":"var(--text-muted)"}),
            ]))
    return html.Table(className="stbl", children=[
        html.Thead(html.Tr([
            html.Th("Syndrome"), html.Th("판정"),
            html.Th("CN"), html.Th("Target Region"),
        ])),
        html.Tbody(rows),
    ])


def _chrom_syndrome_cards(syndromes: dict, chrom: str) -> list:
    chrom = chrom.replace("chr", "")
    items = [s for s in syndromes.values() if s.primary_chrom == chrom]
    if not items:
        return [html.Span("이 염색체에 대한 마커 없음",
                          style={"fontSize":"11px","color":"var(--text-muted)"})]
    cards = []
    for s in items:
        pill_cls, pill_txt = _PILL_MAP.get(s.call, ("pnc", s.call))
        c = s.call_color
        cn_str = f"CN {s.cn_value:.3f}" if s.cn_value is not None else ""
        cards.append(html.Div(style={
            "border": f"1px solid {c}44", "borderLeft": f"4px solid {c}",
            "borderRadius": "6px", "padding": "8px 12px",
            "background": f"{c}0d", "flex": "1", "minWidth": "150px",
        }, children=[
            html.Div([
                html.Span(s.syndrome, style={"fontWeight":"700","fontSize":"12px"}),
                html.Span(f"  {s.nipt_id}", style={"fontSize":"10px","color":"var(--text-muted)","fontFamily":"monospace"}),
            ]),
            html.Div([
                html.Span(pill_txt, className=f"pill {pill_cls}"),
                html.Span(cn_str, style={"fontSize":"11px","color":"var(--text-sub)","marginLeft":"6px","fontFamily":"monospace"}),
            ], style={"marginTop":"5px"}),
        ]))
    return cards


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------
def create_detail_app(
    sample: SampleData,
    initial_chrom: str = "21",
    requests_prefix: str = "/detail/",
    cnv_data: Optional[dict[str, pd.DataFrame]] = None,
    syndromes: Optional[dict[str, NiptSyndrome]] = None,
) -> Dash:
    initial_size = CHROM_SIZES[initial_chrom]
    cnv_data  = cnv_data  or {}
    syndromes = syndromes or {}
    iscn      = build_iscn(sample)
    sex_label = "♀ Female" if sample.sex == "female" else "♂ Male"

    def _get_df(chrom, s, e):
        return get_chrom_df(cnv_data, chrom,
                            fallback_fn=lambda: demo_dataframe(sample, chrom, s, e))

    def _chrom_syn(chrom):
        c = chrom.replace("chr","")
        return [s for s in syndromes.values() if s.primary_chrom == c]

    app = Dash(
        __name__ + "_detail",
        requests_pathname_prefix=requests_prefix,
        routes_pathname_prefix="/",
        suppress_callback_exceptions=True,
        title="Chromosome Detail",
    )

    # overall banner
    bcls, bicon, btitle, bdesc = _overall_call(syndromes)

    # initial data
    _init_df  = _get_df(initial_chrom, 1, initial_size)
    _init_syn = _chrom_syn(initial_chrom)
    _init_fig = region_fig(_init_df, initial_chrom, 1, initial_size, _init_syn)
    _init_fig.update_layout(height=380)
    _sum_fig  = syndrome_summary_fig(syndromes)
    _sum_fig.update_layout(height=max(260, len(syndromes)*26+60))

    app.layout = html.Div([
        dcc.Store(id="detail-selected-chrom", data=initial_chrom),
        dcc.Store(id="brush-region",
                  data={"chrom": initial_chrom, "start": 1, "end": initial_size}),

        html.Div(className="rpt-page", children=[

            # ── Banner ───────────────────────────────────────────────────
            html.Div(className=f"rbanner {bcls}", children=[
                html.Div(bicon, className="ricon"),
                html.Div([
                    html.Div("종합 판정", className="rlbl"),
                    html.Div(btitle, className="rres"),
                    html.Div(bdesc,  className="rdesc"),
                ]),
            ]),

            # ── Sample info ───────────────────────────────────────────────
            html.Div(className="rcard", children=[
                html.Div(className="rch", children=[
                    html.Span("Sample Information", className="rct"),
                ]),
                html.Div(className="igrid", children=[
                    html.Div(className="iitem", children=[
                        html.Div("Sample ID", className="il"),
                        html.Div(sample.id, className="iv"),
                    ]),
                    html.Div(className="iitem", children=[
                        html.Div("Sex", className="il"),
                        html.Div(sex_label, className="iv"),
                    ]),
                    html.Div(className="iitem", children=[
                        html.Div("ISCN", className="il"),
                        html.Div(iscn, className="iv",
                                 style={"fontFamily":"monospace","fontSize":"11px"}),
                    ]),
                    html.Div(className="iitem", children=[
                        html.Div("CNV Events", className="il"),
                        html.Div(str(len(sample.events)), className="iv"),
                    ]),
                ]),
            ]),

            # ── Two-column: left=karyotype / right=CNV+summary ────────────
            html.Div(className="rpt-cols", children=[

                # ── LEFT: Ideogram + brush + syndrome cards ───────────────
                html.Div([

                    # Chromosome header + dropdown
                    html.Div(className="rcard", children=[
                        html.Div(className="rch", children=[
                            html.Div([
                                html.Span("Chromosome",
                                          style={"fontSize":"10px","fontWeight":"700",
                                                 "letterSpacing":".07em","textTransform":"uppercase",
                                                 "color":"var(--text-muted)"}),
                                html.Span(id="detail-chrom-badge",
                                          children=f"chr{initial_chrom}",
                                          style={"fontFamily":"monospace","fontWeight":"700",
                                                 "fontSize":"16px","color":"var(--text)",
                                                 "marginLeft":"8px"}),
                            ], style={"display":"flex","alignItems":"center"}),
                            dcc.Dropdown(
                                id="annotation-select",
                                options=annotation_options_for_chrom(sample, initial_chrom),
                                placeholder="Region 선택 → brush",
                                clearable=True,
                                style={"fontSize":"11px","minWidth":"200px"},
                            ),
                        ]),

                        # Ideogram
                        html.Div(
                            dashbio.Ideogram(
                                id="detail-ideo",
                                organism="human", assembly="GRCh38",
                                orientation="horizontal",
                                chromosomes=[initial_chrom], rows=1,
                                chrHeight=600, chrWidth=22, chrMargin=20,
                                rotatable=False,
                                showBandLabels=True, showChromosomeLabels=True,
                                resolution=850, annotations=None,
                                showAnnotTooltip=False,
                                brush=f"chr{initial_chrom}:1-{initial_size}",
                            ),
                            style={"width":"100%","overflowX":"auto","overflowY":"visible",
                                   "padding":"8px 12px 4px","background":"white"},
                        ),

                        # brush chip
                        html.Div([
                            html.Span(id="brush-chip", className="brush-chip"),
                            html.Span(" ← drag brush or select annotation",
                                      className="brush-hint"),
                        ], style={"padding":"4px 14px 10px","borderTop":"1px solid var(--border)"}),
                    ]),

                    # Syndrome cards for this chrom
                    html.Div(className="rcard", children=[
                        html.Div(className="rch", children=[
                            html.Span("Syndrome · 이 염색체", className="rct"),
                        ]),
                        html.Div(
                            id="syndrome-cards",
                            children=_chrom_syndrome_cards(syndromes, initial_chrom),
                            style={"display":"flex","flexWrap":"wrap","gap":"8px",
                                   "padding":"10px 12px"},
                        ),
                    ]),

                    # Syndrome table (all)
                    html.Div(className="rcard", children=[
                        html.Div(className="rch", children=[
                            html.Span("전체 Syndrome 판정", className="rct"),
                        ]),
                        html.Div(className="rcb", style={"padding":"0"},
                                 children=_syndrome_table(syndromes)),
                    ]),
                ]),

                # ── RIGHT: CN/BAF + summary chart ────────────────────────
                html.Div([

                    # CN/BAF graph
                    html.Div(className="rcard", children=[
                        html.Div(className="rch", children=[
                            html.Span("CN / BAF", className="rct"),
                            html.Span(id="region-title-chip",
                                      children=f"chr{initial_chrom}: full",
                                      style={"fontFamily":"monospace","fontSize":"10px",
                                             "color":"var(--text-muted)"}),
                        ]),
                        dcc.Graph(
                            id="region-graph",
                            figure=_init_fig,
                            config={"displayModeBar":True,"responsive":True,
                                    "modeBarButtonsToRemove":["select2d","lasso2d"],
                                    "toImageButtonOptions":{"format":"png","scale":2}},
                            style={"height":"380px"},
                        ),
                    ]),

                    # Syndrome summary chart
                    html.Div(className="rcard", children=[
                        html.Div(className="rch", children=[
                            html.Span("Syndrome Call Summary", className="rct"),
                        ]),
                        dcc.Graph(
                            id="syndrome-summary-graph",
                            figure=_sum_fig,
                            config={"displayModeBar":False,"responsive":True},
                            style={"height":f"{max(260, len(syndromes)*26+60)}px"},
                        ),
                    ]),
                ]),
            ]),
        ]),
    ])

    # ── inject CSS + postMessage JS ──────────────────────────────────────
    app.index_string = app.index_string.replace(
        "</head>", _CSS + "\n</head>"
    ).replace(
        "</body>", _DETAIL_MESSAGE_JS + "\n</body>"
    )

    # -----------------------------------------------------------------------
    # Callbacks
    # -----------------------------------------------------------------------
    @app.callback(
        Output("detail-ideo", "chromosomes"),
        Output("detail-ideo", "brush"),
        Output("annotation-select", "options"),
        Output("annotation-select", "value"),
        Output("detail-chrom-badge", "children"),
        Output("brush-region", "data"),
        Output("brush-chip", "children"),
        Output("syndrome-cards", "children"),
        Input("detail-selected-chrom", "data"),
    )
    def _change_chrom(chrom):
        chrom = str(chrom or initial_chrom)
        if chrom not in CHROM_SIZES: chrom = initial_chrom
        size = CHROM_SIZES[chrom]
        return (
            [chrom],
            f"chr{chrom}:1-{size}",
            annotation_options_for_chrom(sample, chrom),
            None,
            f"chr{chrom}",
            {"chrom": chrom, "start": 1, "end": size},
            f"chr{chrom}: 0.000 – {size/1e6:.3f} Mb",
            _chrom_syndrome_cards(syndromes, chrom),
        )

    @app.callback(
        Output("detail-ideo", "brush", allow_duplicate=True),
        Output("brush-region", "data", allow_duplicate=True),
        Output("brush-chip", "children", allow_duplicate=True),
        Input("annotation-select", "value"),
        State("detail-selected-chrom", "data"),
        prevent_initial_call=True,
    )
    def _annot_to_brush(value, chrom):
        if not value or not chrom: return no_update, no_update, no_update
        try:
            item = json.loads(value)
            chrom = str(chrom)
            s = max(1, _parse_bp(item.get("start"), 1))
            e = min(CHROM_SIZES[chrom], _parse_bp(item.get("end"), CHROM_SIZES[chrom]))
        except: return no_update, no_update, no_update
        return (f"chr{chrom}:{s}-{e}",
                {"chrom": chrom, "start": s, "end": e},
                f"chr{chrom}: {s/1e6:.3f} – {e/1e6:.3f} Mb")

    @app.callback(
        Output("brush-region", "data", allow_duplicate=True),
        Output("brush-chip", "children", allow_duplicate=True),
        Input("detail-ideo", "brushData"),
        State("detail-selected-chrom", "data"),
        prevent_initial_call=True,
    )
    def _brush_to_region(data, chrom):
        if not data or not chrom: return no_update, no_update
        chrom = str(chrom)
        s = _parse_bp(data.get("start"), 1)
        e = _parse_bp(data.get("end"), CHROM_SIZES[chrom])
        if e < s: s, e = e, s
        s = max(1, s); e = min(CHROM_SIZES[chrom], e)
        if e <= s: return no_update, no_update
        return ({"chrom": chrom, "start": s, "end": e},
                f"chr{chrom}: {s/1e6:.3f} – {e/1e6:.3f} Mb")

    @app.callback(
        Output("region-graph", "figure"),
        Output("region-title-chip", "children"),
        Input("brush-region", "data"),
    )
    def _update_graph(region):
        if not region: return no_update, no_update
        chrom = str(region["chrom"])
        s = _parse_bp(region["start"], 1)
        e = _parse_bp(region["end"], CHROM_SIZES[chrom])
        df  = _get_df(chrom, s, e)
        fig = region_fig(df, chrom, s, e, _chrom_syn(chrom))
        fig.update_layout(height=380)
        return fig, f"chr{chrom}: {s/1e6:.3f}–{e/1e6:.3f} Mb"

    return app
