"""
Dash App 1 – Whole-karyotype Overview.

Two independent rows=1 Ideogram components (chr1-12 / chr13-22+X/Y).
Layout is fluid: no fixed minWidth, chrWidth scales to container.
"""

from __future__ import annotations

from dash import Dash, html, dcc, Input, Output
import dash_bio as dashbio

from ..core.models import SampleData
from ..core.annot import build_iscn, build_raw_annots
from ..ui.styles import PAGE, card, dashbar, _label, badge, _rgba


# ---------------------------------------------------------------------------
# Click JS – supports both ideogram roots
# ---------------------------------------------------------------------------
_OVERVIEW_CLICK_JS = r"""
<script>
(function () {
  const allowed = new Set([
    '1','2','3','4','5','6','7','8','9','10','11','12',
    '13','14','15','16','17','18','19','20','21','22','X','Y'
  ]);
  const ROOT_IDS = ['overview-ideo-top', 'overview-ideo-bottom'];

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
      if (m) { const c = m[1].toUpperCase(); if (allowed.has(c)) return c; }
    }
    return null;
  }

  function chromFromNode(target, root) {
    let node = target;
    while (node && node !== root && node !== document.body) {
      const cands = [];
      if (node.id) cands.push(node.id);
      if (node.getAttribute) {
        cands.push(node.getAttribute('class'));
        cands.push(node.getAttribute('data-chr'));
        cands.push(node.getAttribute('data-chromosome'));
        cands.push(node.getAttribute('aria-label'));
      }
      for (const v of cands) { const c = normaliseChrom(v); if (c) return c; }
      node = node.parentElement;
    }
    const g = target && target.closest ? target.closest('g') : null;
    if (g) {
      for (const t of Array.from(g.querySelectorAll('text'))) {
        const c = normaliseChrom(t.textContent); if (c) return c;
      }
    }
    if (target && target.tagName && target.tagName.toLowerCase() === 'text')
      return normaliseChrom(target.textContent);
    return null;
  }

  function publish(chrom) {
    if (!chrom) return;
    if (window.dash_clientside && window.dash_clientside.set_props)
      window.dash_clientside.set_props('overview-selected-chrom', {data: chrom});
    window.parent.postMessage(
      {type: 'karyotype-chromosome-selected', chrom: chrom},
      window.location.origin
    );
  }

  function installOne(rootId) {
    const root = document.getElementById(rootId);
    if (!root || root.dataset.chromClickInstalled === '1') return !!root;
    root.dataset.chromClickInstalled = '1';
    root.style.cursor = 'pointer';
    root.addEventListener('click', function (ev) {
      const c = chromFromNode(ev.target, root); if (c) publish(c);
    }, true);
    return true;
  }

  let tries = 0;
  const timer = setInterval(function () {
    tries++;
    if (ROOT_IDS.map(installOne).every(Boolean) || tries > 100) clearInterval(timer);
  }, 100);
})();
</script>
"""

# ---------------------------------------------------------------------------
# Annotation tracks
# ---------------------------------------------------------------------------
_ANNOTATION_TRACKS = [
    {"id": "cnv",  "displayName": "CNV",  "color": "#E53E3E", "shape": "rectangle"},
    {"id": "gene", "displayName": "Gene", "color": "#6B46C1", "shape": "circle"},
]

# ---------------------------------------------------------------------------
# Legend – cytoband staining + CNV/Gene
# ---------------------------------------------------------------------------
_LEGEND = [{
    "name": "Cytoband",
    "rows": [
        {"name": "Negative (Gneg)",      "color": "#f0f0f0", "shape": "rectangle"},
        {"name": "Positive 25% (Gpos25)","color": "#c8c8c8", "shape": "rectangle"},
        {"name": "Positive 50% (Gpos50)","color": "#a0a0a0", "shape": "rectangle"},
        {"name": "Positive 75% (Gpos75)","color": "#686868", "shape": "rectangle"},
        {"name": "Positive 100% (Gpos100)","color": "#1e1e1e","shape": "rectangle"},
        {"name": "Centromere (Acen)",    "color": "#E53E3E", "shape": "triangle"},
        {"name": "Variable (Gvar)",      "color": "#a0a0a0", "shape": "rectangle"},
        {"name": "Stalk (Stalk)",        "color": "#6B46C1", "shape": "rectangle"},
    ],
}, {
    "name": "Annotation",
    "rows": [
        {"name": "Trisomy / Gain",  "color": "#FC8181", "shape": "rectangle"},
        {"name": "Monosomy / Loss", "color": "#90CDF4", "shape": "rectangle"},
        {"name": "Partial Gain",    "color": "#FBD38D", "shape": "rectangle"},
        {"name": "Gene",            "color": "#B794F4", "shape": "circle"},
    ],
}]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
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
            html.Div(f"chr{ev.chr}  {ev.type.replace('_', ' ')}",
                     style={"fontSize": "11px", "color": "#718096", "marginTop": "2px"}),
            html.Div(f"CN = {ev.cn}",
                     style={"fontSize": "11px", "color": "#4A5568", "marginTop": "1px"}),
        ]))
    return cards


def _filter_annots(raw_annots: dict, chroms: list[str]) -> dict:
    wanted = {str(c) for c in chroms}
    return {
        "keys": list(raw_annots.get("keys", [])),
        "annots": [r for r in raw_annots.get("annots", []) if str(r.get("chr")) in wanted],
    }


def _ideo_row(component_id: str, chroms: list[str], raw_annots: dict,
              show_legend: bool = False) -> html.Div:
    """
    One fluid-width ideogram row.
    No minWidth → stretches to container, no horizontal whitespace.
    chrWidth/chrMargin kept small so all chromosomes fit even on narrow screens.
    """
    return html.Div(
        style={
            "width": "100%",
            "overflowX": "auto",   # horizontal scroll only when truly too narrow
            "overflowY": "visible",
            "padding": "6px 8px 10px",
            "boxSizing": "border-box",
            "background": "white",
        },
        children=dashbio.Ideogram(
            id=component_id,
            organism="human",
            assembly="GRCh38",
            orientation="vertical",
            chromosomes=chroms,
            rows=1,
            chrHeight=260,
            chrWidth=13,
            chrMargin=16,
            rotatable=False,
            showBandLabels=True,
            showChromosomeLabels=True,
            resolution=550,
            annotations=_filter_annots(raw_annots, chroms),
            annotationsLayout="tracks",
            annotationTracks=_ANNOTATION_TRACKS,
            annotationHeight=7,
            showAnnotTooltip=True,
            legend=_LEGEND if show_legend else None,
        ),
    )


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

    displayed    = [str(c) for c in sample.display_chroms]
    top_chroms   = [c for c in [str(i) for i in range(1, 13)]       if c in displayed]
    bot_chroms   = [c for c in [str(i) for i in range(13, 23)] + ["X","Y"] if c in displayed]

    app = Dash(
        __name__ + "_overview",
        requests_pathname_prefix=requests_prefix,
        routes_pathname_prefix="/",
        suppress_callback_exceptions=True,
        title="Karyotype Overview",
    )

    app.layout = html.Div([
        dcc.Store(id="overview-selected-chrom", data=initial_chrom),

        # ── Sample info ───────────────────────────────────────────────────
        card("SAMPLE INFO",
            html.Div(style={
                "display": "grid",
                "gridTemplateColumns": "repeat(auto-fill, minmax(190px, 1fr))",
                "gap": "14px",
            }, children=[
                html.Div([_label("Sample ID"),
                          html.Div(sample.id, style={"fontWeight":"700","fontSize":"14px","marginTop":"2px"})]),
                html.Div([_label("Sex"),
                          html.Div(sex_label, style={"fontWeight":"700","fontSize":"14px","color":sex_color,"marginTop":"2px"})]),
                html.Div([_label("Karyotype (ISCN)"),
                          html.Div(badge(iscn), style={"marginTop":"4px"})]),
                html.Div([_label("Total chromosomes"),
                          html.Div(str(total_n), style={"fontWeight":"700","fontSize":"20px","color":total_color,"marginTop":"2px"})]),
                html.Div(style={"gridColumn":"1/-1"}, children=[
                    _label("Chromosomal Events"),
                    html.Div(style={"display":"flex","gap":"8px","flexWrap":"wrap","marginTop":"6px"},
                             children=_event_cards(sample)),
                ]),
            ]),
            body_style={"padding": "10px 14px"},
        ),

        # ── Karyotype header (one dashbar, no per-row header) ─────────────
        dashbar("KARYOTYPE OVERVIEW",
            right=html.Span("click chromosome → detail",
                style={"fontSize":"10px","fontWeight":"400","color":"#A0AEC0"}),
        ),

        # ── Row 1: chr1-12 ────────────────────────────────────────────────
        _ideo_row("overview-ideo-top", top_chroms, raw_annots, show_legend=False),

        html.Div(style={"height":"1px","background":"#EDF2F7","margin":"0 8px"}),

        # ── Row 2: chr13-22 + X/Y  (legend here only) ────────────────────
        _ideo_row("overview-ideo-bottom", bot_chroms, raw_annots, show_legend=True),

        # ── Status line ───────────────────────────────────────────────────
        html.Div([
            html.Span("Chromosome을 직접 클릭하세요. ", style={"color":"#718096"}),
            html.Span(id="overview-status",
                      style={"fontFamily":"monospace","color":"#3182CE","fontWeight":"600"}),
        ], style={"fontSize":"10px","padding":"4px 12px 8px",
                  "background":"white","marginBottom":"8px"}),

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
