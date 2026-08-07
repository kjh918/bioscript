"""
Dash App 2 – Single-chromosome Detail view.
CNV 데이터: --cnv 디렉토리에서 사전 로드, 없으면 demo fallback.
Syndrome overlay: NiptSyndrome 리스트를 region_fig에 전달.
"""

from __future__ import annotations

import json
from typing import Optional

import pandas as pd

from dash import Dash, html, dcc, Input, Output, State, no_update, ctx
import dash_bio as dashbio

from ..core.models import SampleData
from ..core.annot import annotation_options_for_chrom
from ..core.data import demo_dataframe
from ..core.reference import CHROM_SIZES
from ..core.cnv_loader import get_chrom_df
from ..core.nipt_markers import NiptSyndrome, CALL_COLORS, load_marker_tsv
from ..ui.styles import PAGE, card, _label
from ..ui.figures import region_fig, syndrome_summary_fig


def _parse_bp(value, default: int = 0) -> int:
    if value is None or value == "":
        return default
    if isinstance(value, (int, float)):
        return int(value)
    try:
        return int(float(str(value).replace(",", "").replace("bp", "").strip()))
    except (TypeError, ValueError):
        return default


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

# call badge styling
_CALL_BADGE = {
    "ABNORMAL":   {"background": "#FED7D7", "color": "#C53030", "border": "1px solid #FC8181"},
    "SUSPICIOUS": {"background": "#FEEBC8", "color": "#C05621", "border": "1px solid #F6AD55"},
    "NORMAL":     {"background": "#C6F6D5", "color": "#276749", "border": "1px solid #68D391"},
    "UNKNOWN":    {"background": "#E2E8F0", "color": "#4A5568", "border": "1px solid #CBD5E0"},
}


def _call_badge(call: str) -> html.Span:
    style = {**_CALL_BADGE.get(call, _CALL_BADGE["UNKNOWN"]),
             "borderRadius": "4px", "padding": "2px 8px",
             "fontSize": "11px", "fontWeight": "700", "fontFamily": "monospace"}
    return html.Span(call, style=style)


def _syndrome_cards(syndromes: dict[str, NiptSyndrome], chrom: str) -> list:
    """해당 염색체와 관련된 syndrome 카드 렌더링."""
    chrom = chrom.replace("chr", "")
    cards = []
    for syn in syndromes.values():
        if syn.primary_chrom != chrom:
            continue
        border_color = syn.call_color
        cards.append(html.Div(style={
            "border":       f"1px solid {border_color}44",
            "borderLeft":   f"4px solid {border_color}",
            "borderRadius": "6px",
            "padding":      "8px 12px",
            "background":   "white",
            "minWidth":     "180px",
            "flex":         "1",
        }, children=[
            html.Div([
                html.Span(syn.syndrome, style={
                    "fontWeight": "700", "fontSize": "12px", "color": "#2D3748",
                }),
                html.Span(f"  {syn.nipt_id}", style={
                    "fontSize": "10px", "color": "#A0AEC0", "fontFamily": "monospace",
                }),
            ]),
            html.Div(style={"display": "flex", "alignItems": "center",
                            "gap": "8px", "marginTop": "5px"}, children=[
                _call_badge(syn.call),
                html.Span(
                    f"CN {syn.cn_value:.2f}" if syn.cn_value is not None else "",
                    style={"fontSize": "11px", "color": "#4A5568", "fontFamily": "monospace"},
                ),
            ]),
            html.Div(style={"marginTop": "5px"}, children=[
                html.Span(syn.group, style={
                    "fontSize": "9px", "color": "#718096",
                    "background": syn.group_color + "33",
                    "border": f"1px solid {syn.group_color}66",
                    "borderRadius": "3px", "padding": "1px 5px",
                }),
            ]),
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

    def _get_df(chrom: str, start_bp: int, end_bp: int) -> pd.DataFrame:
        return get_chrom_df(
            cnv_data, chrom,
            fallback_fn=lambda: demo_dataframe(sample, chrom, start_bp, end_bp),
        )

    def _chrom_syndromes(chrom: str) -> list[NiptSyndrome]:
        chrom = chrom.replace("chr", "")
        return [s for s in syndromes.values() if s.primary_chrom == chrom]

    app = Dash(
        __name__ + "_detail",
        requests_pathname_prefix=requests_prefix,
        routes_pathname_prefix="/",
        suppress_callback_exceptions=True,
        title="Chromosome Detail",
    )

    # initial figures
    _init_df  = _get_df(initial_chrom, 1, initial_size)
    _init_syn = _chrom_syndromes(initial_chrom)

    app.layout = html.Div([
        dcc.Store(id="detail-selected-chrom", data=initial_chrom),
        dcc.Store(id="brush-region",
                  data={"chrom": initial_chrom, "start": 1, "end": initial_size}),

        # ── header ───────────────────────────────────────────────────────
        html.Div([
            html.Div([
                _label("CHROMOSOME DETAIL"),
                html.Span(
                    id="detail-chrom-badge", children=f"chr{initial_chrom}",
                    style={"fontFamily": "monospace", "fontWeight": "700",
                           "fontSize": "17px", "color": "#2D3748", "marginLeft": "8px"},
                ),
            ], style={"display": "flex", "alignItems": "center"}),

            dcc.Dropdown(
                id="annotation-select",
                options=annotation_options_for_chrom(sample, initial_chrom),
                placeholder="CNV / Gene region 선택 → brush 이동",
                clearable=True,
                style={"fontSize": "11px", "minWidth": "280px"},
            ),
        ], style={
            "display": "flex", "justifyContent": "space-between",
            "alignItems": "center", "flexWrap": "wrap", "gap": "10px",
            "padding": "8px 12px 6px",
            "background": "white", "borderBottom": "1px solid #E2E8F0",
        }),

        # ── Ideogram ─────────────────────────────────────────────────────
        html.Div(
            dashbio.Ideogram(
                id="detail-ideo",
                organism="human", assembly="GRCh38",
                orientation="horizontal", chromosomes=[initial_chrom],
                rows=1, chrHeight=600, chrWidth=22, chrMargin=20,
                rotatable=False,
                showBandLabels=True, showChromosomeLabels=True,
                resolution=850, annotations=None, showAnnotTooltip=False,
                brush=f"chr{initial_chrom}:1-{initial_size}",
            ),
            style={
                "width": "100%", "overflowX": "auto", "overflowY": "visible",
                "padding": "8px 12px 4px", "boxSizing": "border-box",
                "background": "white",
            },
        ),

        # brush chip
        html.Div([
            html.Span(id="brush-chip", style={
                "fontFamily": "monospace", "fontSize": "11px",
                "color": "#3182CE", "fontWeight": "600",
            }),
            html.Span("  ← drag brush or select annotation above",
                      style={"fontSize": "10px", "color": "#A0AEC0"}),
        ], style={
            "padding": "4px 12px 8px", "background": "white",
            "marginBottom": "8px", "borderTop": "1px solid #EDF2F7",
        }),

        # ── Syndrome cards for this chrom ─────────────────────────────────
        html.Div(
            id="syndrome-cards",
            children=_syndrome_cards(syndromes, initial_chrom),
            style={
                "display": "flex", "flexWrap": "wrap", "gap": "8px",
                "padding": "8px 12px", "background": "#F7FAFC",
                "borderRadius": "7px", "marginBottom": "8px",
            },
        ),

        # ── CN / BAF graph ────────────────────────────────────────────────
        card(
            "CN / BAF REGION DETAIL",
            dcc.Graph(
                id="region-graph",
                figure=region_fig(_init_df, initial_chrom, 1, initial_size, _init_syn),
                config={
                    "displayModeBar": True, "responsive": True,
                    "modeBarButtonsToRemove": ["select2d", "lasso2d"],
                    "toImageButtonOptions": {"format": "png", "scale": 2},
                },
                style={"height": "400px"},
            ),
            right=html.Span(
                id="region-title-chip", children=f"chr{initial_chrom}: full",
                style={"fontFamily": "monospace", "fontWeight": "400"},
            ),
            body_style={"padding": "6px 8px"},
        ),

        # ── Syndrome summary chart ────────────────────────────────────────
        card(
            "NIPT SYNDROME CALL SUMMARY",
            dcc.Graph(
                id="syndrome-summary-graph",
                figure=syndrome_summary_fig(syndromes),
                config={"displayModeBar": False, "responsive": True},
                style={"height": f"{max(260, len(syndromes)*28+60)}px"},
            ),
            body_style={"padding": "4px 8px"},
        ),

    ], style=PAGE)

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
        if chrom not in CHROM_SIZES:
            chrom = initial_chrom
        size   = CHROM_SIZES[chrom]
        return (
            [chrom],
            f"chr{chrom}:1-{size}",
            annotation_options_for_chrom(sample, chrom),
            None,
            f"chr{chrom}",
            {"chrom": chrom, "start": 1, "end": size},
            f"chr{chrom}: 0.000 – {size/1e6:.3f} Mb",
            _syndrome_cards(syndromes, chrom),
        )

    @app.callback(
        Output("detail-ideo", "brush", allow_duplicate=True),
        Output("brush-region", "data", allow_duplicate=True),
        Output("brush-chip", "children", allow_duplicate=True),
        Input("annotation-select", "value"),
        State("detail-selected-chrom", "data"),
        prevent_initial_call=True,
    )
    def _annotation_to_brush(value, chrom):
        if not value or not chrom:
            return no_update, no_update, no_update
        try:
            item     = json.loads(value)
            chrom    = str(chrom)
            start_bp = max(1, _parse_bp(item.get("start"), 1))
            end_bp   = min(CHROM_SIZES[chrom], _parse_bp(item.get("end"), CHROM_SIZES[chrom]))
        except Exception:
            return no_update, no_update, no_update
        return (
            f"chr{chrom}:{start_bp}-{end_bp}",
            {"chrom": chrom, "start": start_bp, "end": end_bp},
            f"chr{chrom}: {start_bp/1e6:.3f} – {end_bp/1e6:.3f} Mb",
        )

    @app.callback(
        Output("brush-region", "data", allow_duplicate=True),
        Output("brush-chip", "children", allow_duplicate=True),
        Input("detail-ideo", "brushData"),
        State("detail-selected-chrom", "data"),
        prevent_initial_call=True,
    )
    def _brush_to_region(data, chrom):
        if not data or not chrom:
            return no_update, no_update
        chrom    = str(chrom)
        start_bp = _parse_bp(data.get("start"), 1)
        end_bp   = _parse_bp(data.get("end"), CHROM_SIZES[chrom])
        if end_bp < start_bp:
            start_bp, end_bp = end_bp, start_bp
        start_bp = max(1, start_bp)
        end_bp   = min(CHROM_SIZES[chrom], end_bp)
        if end_bp <= start_bp:
            return no_update, no_update
        return (
            {"chrom": chrom, "start": start_bp, "end": end_bp},
            f"chr{chrom}: {start_bp/1e6:.3f} – {end_bp/1e6:.3f} Mb",
        )

    @app.callback(
        Output("region-graph", "figure"),
        Output("region-title-chip", "children"),
        Input("brush-region", "data"),
    )
    def _update_graph(region):
        if not region:
            return no_update, no_update
        chrom    = str(region["chrom"])
        start_bp = _parse_bp(region["start"], 1)
        end_bp   = _parse_bp(region["end"], CHROM_SIZES[chrom])
        df       = _get_df(chrom, start_bp, end_bp)
        syn_list = _chrom_syndromes(chrom)
        title    = f"chr{chrom}: {start_bp/1e6:.3f}–{end_bp/1e6:.3f} Mb"
        return region_fig(df, chrom, start_bp, end_bp, syn_list), title

    app.index_string = app.index_string.replace(
        "</body>", _DETAIL_MESSAGE_JS + "\n</body>"
    )
    return app
