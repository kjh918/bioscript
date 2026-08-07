"""
Dash App 1 – Whole-karyotype Overview.
"""

from __future__ import annotations

from dash import Dash, html, dcc, Input, Output
import dash_bio as dashbio

from ..core.models import SampleData
from ..core.annot import build_iscn, build_raw_annots
from ..ui.styles import PAGE, DASHBAR_STYLE, card, dashbar, _label, badge, _rgba


# ---------------------------------------------------------------------------
# Injected JS: chromosome click → postMessage to host
# ---------------------------------------------------------------------------
_OVERVIEW_CLICK_JS = r"""
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
    const g = target && target.closest ? target.closest('g') : null;
    if (g) {
      for (const t of Array.from(g.querySelectorAll('text'))) {
        const c = normaliseChrom(t.textContent);
        if (c) return c;
      }
    }
    if (target && target.tagName && target.tagName.toLowerCase() === 'text') {
      return normaliseChrom(target.textContent);
    }
    return null;
  }

  function publish(chrom) {
    if (!chrom) return;
    if (window.dash_clientside && window.dash_clientside.set_props) {
      window.dash_clientside.set_props('overview-selected-chrom', {data: chrom});
    }
    window.parent.postMessage(
      {type: 'karyotype-chromosome-selected', chrom: chrom},
      window.location.origin
    );
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

_ANNOTATION_TRACKS = [
    {"id": "cnv",  "displayName": "CNV",  "color": "#E53E3E", "shape": "rectangle"},
    {"id": "gene", "displayName": "Gene", "color": "#6B46C1", "shape": "circle"},
]
_LEGEND = [{"name": "Legend", "rows": [
    {"name": "Trisomy / Gain",  "color": "#FC8181", "shape": "rectangle"},
    {"name": "Monosomy / Loss", "color": "#90CDF4", "shape": "rectangle"},
    {"name": "Partial Gain",    "color": "#FBD38D", "shape": "rectangle"},
    {"name": "Gene",            "color": "#B794F4", "shape": "circle"},
]}]


def _event_cards(sample: SampleData) -> list:
    cards = []
    for ev in sample.events:
        cards.append(html.Div(style={
            "border":       f"1px solid {ev.color}88",
            "borderLeft":   f"4px solid {ev.color}",
            "borderRadius": "5px",
            "padding":      "7px 12px",
            "background":   _rgba(ev.color, 0.06),
            "minWidth":     "140px",
        }, children=[
            html.Div(ev.iscn, style={
                "fontWeight": "700", "fontFamily": "monospace",
                "fontSize": "13px", "color": ev.color,
            }),
            html.Div(
                f"chr{ev.chr}  {ev.type.replace('_', ' ')}",
                style={"fontSize": "11px", "color": "#718096", "marginTop": "2px"},
            ),
            html.Div(
                f"CN = {ev.cn}",
                style={"fontSize": "11px", "color": "#4A5568", "marginTop": "1px"},
            ),
        ]))
    return cards


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------
def create_overview_app(
    sample: SampleData,
    initial_chrom: str = "5",
    requests_prefix: str = "/overview/",
) -> Dash:
    iscn       = build_iscn(sample)
    raw_annots = build_raw_annots(sample)
    sex_label  = "♀ Female" if sample.sex == "female" else "♂ Male"
    sex_color  = "#E53E3E"  if sample.sex == "female" else "#3182CE"

    try:
        total_n = int(iscn.split(",")[0])
    except (ValueError, IndexError):
        total_n = 46
    total_color = "#E53E3E" if total_n != 46 else "#2D3748"

    app = Dash(
        __name__ + "_overview",
        requests_pathname_prefix=requests_prefix,
        routes_pathname_prefix="/",
        suppress_callback_exceptions=True,
        title="Karyotype Overview",
    )

    app.layout = html.Div([
        dcc.Store(id="overview-selected-chrom", data=initial_chrom),

        # ── Sample info (card 사용) ────────────────────────────────────────
        card(
            "SAMPLE INFO",
            html.Div(style={
                "display": "grid",
                "gridTemplateColumns": "repeat(auto-fill, minmax(200px, 1fr))",
                "gap": "16px",
            }, children=[
                html.Div([
                    _label("Sample ID"),
                    html.Div(sample.id, style={
                        "fontWeight": "700", "fontSize": "14px", "marginTop": "2px",
                    }),
                ]),
                html.Div([
                    _label("Sex"),
                    html.Div(sex_label, style={
                        "fontWeight": "700", "fontSize": "14px",
                        "color": sex_color, "marginTop": "2px",
                    }),
                ]),
                html.Div([
                    _label("Karyotype (ISCN)"),
                    html.Div(badge(iscn), style={"marginTop": "4px"}),
                ]),
                html.Div([
                    _label("Total chromosomes"),
                    html.Div(str(total_n), style={
                        "fontWeight": "700", "fontSize": "20px",
                        "color": total_color, "marginTop": "2px",
                    }),
                ]),
                html.Div(style={"gridColumn": "1 / -1"}, children=[
                    _label("Chromosomal Events"),
                    html.Div(style={
                        "display": "flex", "gap": "8px",
                        "flexWrap": "wrap", "marginTop": "6px",
                    }, children=_event_cards(sample)),
                ]),
            ]),
            body_style={"padding": "10px 14px"},
        ),

        # ── Karyotype Ideogram (card 없이 직접 배치) ──────────────────────
        # dashbar만 단독 렌더링 – card로 감싸면 고정 height에 Ideogram이 잘림
        dashbar(
            "KARYOTYPE OVERVIEW",
            right=html.Span(
                "click chromosome → detail",
                style={"fontSize": "10px", "fontWeight": "400",
                       "letterSpacing": "0", "color": "#A0AEC0"},
            ),
        ),

        # Ideogram wrap: height 고정 없음 → Ideogram이 자체 크기 결정
        html.Div(
            dashbio.Ideogram(
                id="overview-ideo",
                organism="human",
                assembly="GRCh38",
                orientation="vertical",
                chromosomes=sample.display_chroms,
                rows=1,
                chrHeight=420,
                chrWidth=15,
                chrMargin=14,
                rotatable=False,
                showBandLabels=True,
                showChromosomeLabels=True,
                resolution=550,
                annotations=raw_annots,
                annotationsLayout="tracks",
                annotationTracks=_ANNOTATION_TRACKS,
                annotationHeight=8,
                barWidth=3,
                histogramScaling="relative",
                showAnnotTooltip=True,
                legend=_LEGEND,
            ),
            id="overview-ideo-wrap",
            style={
                "width": "100%",
                "overflowX": "auto",
                "overflowY": "visible",
                "padding": "8px 12px 14px",
                "boxSizing": "border-box",
                "background": "white",
            },
        ),

        html.Div([
            html.Span("Chromosome을 직접 클릭하세요. ",
                      style={"color": "#718096"}),
            html.Span(id="overview-status", style={
                "fontFamily": "monospace", "color": "#3182CE", "fontWeight": "600",
            }),
        ], style={
            "fontSize": "10px",
            "padding": "4px 12px 8px",
            "background": "white",
            "marginBottom": "8px",
        }),

    ], style=PAGE)

    @app.callback(
        Output("overview-status", "children"),
        Input("overview-selected-chrom", "data"),
    )
    def _status(chrom):
        return f"chr{chrom} selected"

    app.index_string = app.index_string.replace(
        "</body>", _OVERVIEW_CLICK_JS + "\n</body>"
    )
    return app
