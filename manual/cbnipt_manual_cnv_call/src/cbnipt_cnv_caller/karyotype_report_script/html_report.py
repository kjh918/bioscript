from __future__ import annotations
import html, json
from pathlib import Path
from config import ALL_CHROMS, IDEOGRAM_CDN, DEFAULT_COLORS


def _e(v) -> str:
    return html.escape('' if v is None else str(v))


def _cnv_badge(type_: str) -> str:
    color = DEFAULT_COLORS.get(type_, '#718096')
    bg    = {'gain': '#FEF2F2', 'loss': '#EFF6FF'}.get(type_, '#F8FAFC')
    tc    = {'gain': '#DC2626', 'loss': '#2563EB'}.get(type_, '#64748B')
    return (f'<span style="background:{bg};color:{tc};border:1px solid {color}44;'
            f'border-radius:4px;padding:1px 7px;font-size:10px;font-weight:700;'
            f'text-transform:uppercase;letter-spacing:.05em">{_e(type_)}</span>')


def _cnv_table(rows: list[dict]) -> str:
    if not rows:
        return '<div class="empty">CNV 결과 없음</div>'
    trs = []
    for r in rows:
        cn_txt = '' if r['cn'] is None else f'{float(r["cn"]):.2f}'
        t      = r['type']
        bar_w  = min(100, max(0, r['length_mb'] / 250 * 100))
        bar_c  = DEFAULT_COLORS.get(t, '#94A3B8')
        trs.append(
            f'<tr>'
            f'<td><span class="chrom-chip">chr{_e(r["chrom"])}</span></td>'
            f'<td><b>{_e(r["name"])}</b></td>'
            f'<td>{_cnv_badge(t)}</td>'
            f'<td class="mono">{cn_txt}</td>'
            f'<td class="mono">{r["start"]/1e6:.3f}</td>'
            f'<td class="mono">{r["end"]/1e6:.3f}</td>'
            f'<td>'
            f'  <div class="len-bar-wrap">'
            f'    <div class="len-bar" style="width:{bar_w:.1f}%;background:{bar_c}"></div>'
            f'    <span class="mono">{r["length_mb"]:.2f}</span>'
            f'  </div>'
            f'</td>'
            f'</tr>'
        )
    return '''<div class="table-wrap"><table class="cnv-table">
<thead><tr>
  <th>Chr</th><th>Name</th><th>Type</th><th>CN</th>
  <th>Start Mb</th><th>End Mb</th><th>Length Mb</th>
</tr></thead>
<tbody>''' + ''.join(trs) + '</tbody></table></div>'


def render_report(
    output,
    sample_id: str,
    sex: str,
    annotations: list[dict],
    cnv_rows: list[dict],
    plot_html: dict[str, str],
) -> Path:
    output    = Path(output)
    ann_json  = json.dumps(annotations, ensure_ascii=False)
    chroms    = [c for c in ALL_CHROMS if not (sex.lower() == 'female' and c == 'Y')]
    chr_json  = json.dumps(chroms)

    gain_n = sum(1 for r in cnv_rows if r['type'] == 'gain')
    loss_n = sum(1 for r in cnv_rows if r['type'] == 'loss')

    # TOC + plot sections
    toc_items = ''.join(
        f'<a class="toc-chip" href="#chr-{_e(c)}">chr{_e(c)}</a>'
        for c in plot_html
    )
    plot_sections = ''.join(
        f'''<section class="card chr-panel" id="chr-{_e(c)}">
  <div class="card-head">
    <span class="chrom-chip lg">chr{_e(c)}</span>
    <span class="head-title">CN / BAF</span>
    <a class="back-top" href="#top">↑ top</a>
  </div>
  <div class="plot-body">{ph}</div>
</section>'''
        for c, ph in plot_html.items()
    )

    page = f'''<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>{_e(sample_id)} · Karyotype Report</title>
<script src="{IDEOGRAM_CDN}"></script>
<style>
/* ── Reset & Base ─────────────────────────────────────── */
*{{box-sizing:border-box;margin:0;padding:0}}
:root{{
  --bg:#F1F5F9;
  --card:#FFFFFF;
  --border:#E2E8F0;
  --text:#1E293B;
  --muted:#64748B;
  --accent:#2563EB;
  --gain:#DC2626;
  --loss:#2563EB;
  --gene:#7C3AED;
  --radius:10px;
  --shadow:0 1px 4px rgba(0,0,0,.07),0 4px 12px rgba(0,0,0,.04);
  --mono:'JetBrains Mono','Fira Mono','Courier New',monospace;
}}
body{{
  background:var(--bg);color:var(--text);
  font-family:Inter,'Segoe UI',system-ui,sans-serif;
  font-size:13px;line-height:1.6;
}}

/* ── Layout ───────────────────────────────────────────── */
#top{{scroll-margin-top:0}}
.page{{max-width:1400px;margin:0 auto;padding:18px 20px 40px;display:flex;flex-direction:column;gap:14px}}

/* ── Card ─────────────────────────────────────────────── */
.card{{background:var(--card);border:1px solid var(--border);border-radius:var(--radius);box-shadow:var(--shadow);overflow:hidden}}
.card-head{{
  display:flex;align-items:center;gap:10px;
  padding:9px 16px;border-bottom:1px solid var(--border);
  background:#F8FAFC;
  font-size:10px;font-weight:700;color:var(--muted);
  letter-spacing:.1em;text-transform:uppercase;
}}
.head-title{{flex:1}}
.back-top{{color:var(--accent);text-decoration:none;font-size:11px;font-weight:400;text-transform:none;letter-spacing:0}}
.back-top:hover{{text-decoration:underline}}
.card-body{{padding:14px 16px}}

/* ── Header strip ─────────────────────────────────────── */
.site-header{{
  background:var(--text);color:#F8FAFC;
  padding:14px 20px;display:flex;align-items:center;gap:12px;
}}
.site-header h1{{font-size:14px;font-weight:700;letter-spacing:.04em}}
.site-header .sample-badge{{
  font-family:var(--mono);font-size:12px;
  background:rgba(255,255,255,.12);border:1px solid rgba(255,255,255,.2);
  border-radius:5px;padding:2px 10px;color:#e2e8f0;
}}
.sex-badge{{
  font-size:11px;font-weight:600;padding:2px 9px;border-radius:4px;
  background:rgba(255,255,255,.1);color:#CBD5E1;
}}

/* ── Stats grid ───────────────────────────────────────── */
.stats{{display:grid;grid-template-columns:repeat(auto-fill,minmax(130px,1fr));gap:12px}}
.stat{{background:#F8FAFC;border:1px solid var(--border);border-radius:7px;padding:12px 14px}}
.stat-label{{font-size:10px;color:var(--muted);text-transform:uppercase;letter-spacing:.07em;font-weight:700}}
.stat-value{{font-size:20px;font-weight:800;color:var(--text);margin-top:3px;font-variant-numeric:tabular-nums}}
.stat-value.gain{{color:var(--gain)}}
.stat-value.loss{{color:var(--loss)}}
.stat-sub{{font-size:10px;color:var(--muted);margin-top:1px}}

/* ── Event cards ──────────────────────────────────────── */
.events{{display:flex;gap:8px;flex-wrap:wrap;margin-top:4px}}
.event-card{{
  border-radius:7px;padding:9px 13px;min-width:150px;
  border-left-width:4px;border-left-style:solid;
  font-size:12px;
}}
.event-card .ev-name{{font-weight:800;font-family:var(--mono);font-size:13px}}
.event-card .ev-meta{{font-size:10px;color:var(--muted);margin-top:2px}}

/* ── Ideogram ─────────────────────────────────────────── */
.ideo-wrap{{overflow-x:auto;overflow-y:visible;padding:4px 0 12px;min-height:360px}}
.ideo-inner{{min-width:1200px}}
.legend{{display:flex;gap:16px;flex-wrap:wrap;margin-top:6px;font-size:11px;color:var(--muted);align-items:center}}
.sw{{display:inline-block;width:9px;height:9px;margin-right:4px;vertical-align:-1px;border-radius:1px}}
.sw.circle{{border-radius:50%}}

/* ── CNV table ────────────────────────────────────────── */
.table-wrap{{overflow-x:auto}}
.cnv-table{{width:100%;border-collapse:collapse;font-size:12px}}
.cnv-table th{{
  padding:7px 10px;text-align:left;white-space:nowrap;
  font-size:10px;color:var(--muted);font-weight:700;
  text-transform:uppercase;letter-spacing:.07em;
  background:#F8FAFC;border-bottom:2px solid var(--border);
}}
.cnv-table td{{padding:8px 10px;border-bottom:1px solid var(--border);vertical-align:middle}}
.cnv-table tr:last-child td{{border-bottom:none}}
.cnv-table tr:hover td{{background:#F8FAFC}}

.chrom-chip{{
  display:inline-block;
  background:#EFF6FF;color:#1D4ED8;
  border:1px solid #BFDBFE;border-radius:4px;
  padding:1px 7px;font-size:11px;font-weight:700;
  font-family:var(--mono);
}}
.chrom-chip.lg{{font-size:12px;padding:2px 10px}}

.mono{{font-family:var(--mono);font-size:11px}}

.len-bar-wrap{{display:flex;align-items:center;gap:8px;min-width:120px}}
.len-bar{{height:5px;border-radius:2px;min-width:2px;transition:width .3s}}

/* ── TOC ──────────────────────────────────────────────── */
.toc{{display:flex;gap:5px;flex-wrap:wrap}}
.toc-chip{{
  text-decoration:none;font-size:11px;font-weight:600;
  color:var(--accent);background:#EFF6FF;
  border:1px solid #BFDBFE;border-radius:4px;
  padding:3px 8px;font-family:var(--mono);
  transition:background .15s;
}}
.toc-chip:hover{{background:#DBEAFE}}

/* ── Chr panels ───────────────────────────────────────── */
.chr-panel{{scroll-margin-top:16px}}
.plot-body{{padding:6px 10px 10px}}
.empty{{color:var(--muted);padding:20px;text-align:center}}

@media(max-width:700px){{
  .page{{padding:10px}}.stats{{grid-template-columns:1fr 1fr}}.ideo-inner{{min-width:980px}}
}}
</style>
</head>
<body id="top">

<!-- Site Header -->
<header class="site-header">
  <span>🧬</span>
  <h1>Karyotype Report</h1>
  <span class="sample-badge">{_e(sample_id)}</span>
  <span class="sex-badge">{'♀ Female' if sex.lower()=='female' else '♂ Male'}</span>
</header>

<div class="page">

<!-- 1. Sample Info -->
<section class="card">
  <div class="card-head"><span class="head-title">Sample Info</span></div>
  <div class="card-body">
    <div class="stats">
      <div class="stat">
        <div class="stat-label">Sample ID</div>
        <div class="stat-value" style="font-size:14px">{_e(sample_id)}</div>
      </div>
      <div class="stat">
        <div class="stat-label">Sex</div>
        <div class="stat-value" style="font-size:16px">{'♀ Female' if sex.lower()=='female' else '♂ Male'}</div>
      </div>
      <div class="stat">
        <div class="stat-label">CNV Calls</div>
        <div class="stat-value">{len(cnv_rows)}</div>
      </div>
      <div class="stat">
        <div class="stat-label">Gains</div>
        <div class="stat-value gain">{gain_n}</div>
      </div>
      <div class="stat">
        <div class="stat-label">Losses</div>
        <div class="stat-value loss">{loss_n}</div>
      </div>
      <div class="stat">
        <div class="stat-label">Chromosomes</div>
        <div class="stat-value">{len(set(r["chrom"] for r in cnv_rows))}</div>
        <div class="stat-sub">affected</div>
      </div>
    </div>

    <!-- Event cards -->
    <div style="margin-top:14px">
      <div class="stat-label" style="margin-bottom:8px">Detected Events</div>
      <div class="events">
        {_event_cards(cnv_rows)}
      </div>
    </div>
  </div>
</section>

<!-- 2. Karyotype Overview -->
<section class="card">
  <div class="card-head"><span class="head-title">Karyotype Overview</span>
    <span style="font-weight:400;font-size:11px;color:#94A3B8;text-transform:none;letter-spacing:0">ideogram.js · GRCh38</span>
  </div>
  <div class="card-body">
    <div class="ideo-wrap"><div class="ideo-inner" id="overview-ideo"></div></div>
    <div class="legend">
      <span><span class="sw" style="background:var(--gain)"></span>CN Gain</span>
      <span><span class="sw" style="background:var(--loss)"></span>CN Loss</span>
      <span><span class="sw circle" style="background:var(--gene)"></span>Gene</span>
    </div>
  </div>
</section>

<!-- 3. CNV Table -->
<section class="card">
  <div class="card-head"><span class="head-title">CNV Summary</span>
    <span style="font-weight:400;font-size:11px;color:#94A3B8;text-transform:none;letter-spacing:0">{len(cnv_rows)} events</span>
  </div>
  <div class="card-body" style="padding:0">
    {_cnv_table(cnv_rows)}
  </div>
</section>

<!-- 4. Chromosome Plots TOC -->
<section class="card">
  <div class="card-head"><span class="head-title">Chromosome Plots</span></div>
  <div class="card-body"><div class="toc">{toc_items}</div></div>
</section>

<!-- 5. Per-chromosome CN/BAF plots -->
{plot_sections}

</div>

<script>
const annotations  = {ann_json};
const chromosomes  = {chr_json};

window.addEventListener('load', () => {{
  new Ideogram({{
    organism:'human', assembly:'GRCh38',
    container:'#overview-ideo',
    orientation:'vertical',
    chromosomes,
    chrHeight:300, chrWidth:14, chrMargin:20, rows:1,
    showBandLabels:true, showChromosomeLabels:true, resolution:550,
    annotations,
    annotationsLayout:'overlay',
    annotationHeight:6,
    showAnnotTooltip:true,
    onLoad() {{
      document.querySelectorAll('#overview-ideo g.chromosome').forEach(el => {{
        el.style.cursor = 'pointer';
        el.addEventListener('click', () => {{
          const m = (el.id||'').match(/chr([0-9]{{1,2}}|X|Y)/i)
                 || (el.textContent||'').match(/([0-9]{{1,2}}|X|Y)$/);
          if (m) {{
            const t = document.getElementById('chr-' + m[1].toUpperCase());
            if (t) t.scrollIntoView({{behavior:'smooth', block:'start'}});
          }}
        }});
      }});
    }}
  }});
}});
</script>
</body></html>'''

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding='utf-8')
    return output


def _event_cards(cnv_rows):
    if not cnv_rows:
        return '<div class="empty">No events</div>'
    cards = []
    for r in cnv_rows:
        t  = r['type']
        c  = DEFAULT_COLORS.get(t, '#718096')
        bg = {'gain':'#FEF2F2','loss':'#EFF6FF'}.get(t,'#F8FAFC')
        cn_txt = '' if r['cn'] is None else f'  CN={float(r["cn"]):.2f}'
        cards.append(
            f'<div class="event-card" style="background:{bg};border-left-color:{c}">'
            f'  <div class="ev-name" style="color:{c}">{_e(r["name"])}</div>'
            f'  <div class="ev-meta">chr{_e(r["chrom"])} · {t}{cn_txt}</div>'
            f'  <div class="ev-meta">{r["start"]/1e6:.2f}–{r["end"]/1e6:.2f} Mb ({r["length_mb"]:.2f} Mb)</div>'
            f'</div>'
        )
    return ''.join(cards)
