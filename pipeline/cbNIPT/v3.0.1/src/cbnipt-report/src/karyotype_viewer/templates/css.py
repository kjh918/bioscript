"""
CSS 디자인 토큰 + 컴포넌트 스타일.
f-string 없이 순수 문자열 → reporter.py에서 <style> 태그로 삽입.
"""

CSS = """
:root{
  --bg:#F4F6F9; --surface:#fff; --surface2:#F8F9FB;
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
*{box-sizing:border-box;margin:0;padding:0;}
body{font-family:var(--sans);background:var(--bg);color:var(--text);font-size:13px;line-height:1.55;}
@media print{
  body{background:#fff;font-size:11px;}
  .no-print{display:none!important;}
  .page{max-width:100%;padding:0;}
}
.page{max-width:1100px;margin:0 auto;padding:14px 16px 60px;}

/* ── Header ─────────────────────────────────────────────────────── */
.hdr{background:var(--navy);color:#fff;border-radius:12px;padding:1.2rem 1.6rem;
     margin-bottom:1rem;display:flex;justify-content:space-between;
     align-items:flex-start;flex-wrap:wrap;gap:1rem;}
.hdr-brand{font-size:10px;font-weight:700;letter-spacing:.12em;opacity:.5;margin-bottom:4px;}
.hdr h1{font-size:18px;font-weight:700;}
.hdr-sub{font-size:11px;opacity:.45;margin-top:2px;}
.hdr-right{text-align:right;}
.hdr-lbl{font-size:10px;opacity:.4;letter-spacing:.05em;}
.hdr-val{font-size:13px;font-weight:600;font-family:var(--mono);}

/* ── Banner ──────────────────────────────────────────────────────── */
.banner{border-radius:var(--r);padding:.9rem 1.2rem;margin-bottom:1rem;
        display:flex;align-items:center;gap:1rem;flex-wrap:wrap;}
.banner.NORMAL   {background:var(--teal-l);border:1.5px solid var(--teal-b);}
.banner.ABNORMAL {background:var(--red-l); border:1.5px solid var(--red-b);}
.banner.SUSPICIOUS{background:var(--amber-l);border:1.5px solid var(--amber-b);}
.bicon{width:40px;height:40px;border-radius:50%;display:flex;align-items:center;
       justify-content:center;font-size:17px;flex-shrink:0;}
.NORMAL   .bicon{background:var(--teal-b);}
.ABNORMAL .bicon{background:var(--red-b);}
.SUSPICIOUS .bicon{background:var(--amber-b);}
.bres{font-size:16px;font-weight:700;}
.NORMAL .bres{color:var(--teal);}
.ABNORMAL .bres{color:var(--red);}
.SUSPICIOUS .bres{color:var(--amber);}
.bdesc{font-size:12px;color:var(--text-sub);margin-top:2px;}

/* ── Cards ───────────────────────────────────────────────────────── */
.card{background:var(--surface);border:1px solid var(--border);
      border-radius:var(--r);margin-bottom:1rem;overflow:hidden;}
.ch{display:flex;align-items:center;justify-content:space-between;
    padding:8px 14px;border-bottom:1px solid var(--border);background:var(--surface2);}
.ct{font-size:10px;font-weight:700;letter-spacing:.07em;
    text-transform:uppercase;color:var(--text-muted);}
.cb{padding:12px 14px;}

/* ── Info grid ───────────────────────────────────────────────────── */
.igrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));}
.iitem{padding:8px 12px;border-right:1px solid var(--border);border-bottom:1px solid var(--border);}
.iitem .il{font-size:10px;color:var(--text-muted);margin-bottom:1px;}
.iitem .iv{font-size:13px;font-weight:500;}

/* ── Chips ───────────────────────────────────────────────────────── */
.chips{display:flex;gap:6px;flex-wrap:wrap;margin-bottom:1rem;}
.chip{display:inline-flex;align-items:center;gap:5px;padding:5px 11px;
      border-radius:99px;font-size:11px;font-weight:600;border:1px solid;}
.chip-ab{background:var(--red-l);color:var(--red);border-color:var(--red-b);}
.chip-su{background:var(--amber-l);color:var(--amber);border-color:var(--amber-b);}
.chip-ok{background:var(--teal-l);color:var(--teal);border-color:var(--teal-b);}
.chip-ev{background:var(--navy-l);color:var(--navy);border-color:var(--navy-b);}

/* ── Pills ───────────────────────────────────────────────────────── */
.pill{display:inline-block;font-size:10px;font-weight:700;
      padding:2px 8px;border-radius:99px;border:1px solid;white-space:nowrap;}
.plr{background:var(--teal-l);color:var(--teal);border-color:var(--teal-b);}
.phr{background:var(--red-l);color:var(--red);border-color:var(--red-b);}
.pmo{background:var(--amber-l);color:var(--amber);border-color:var(--amber-b);}
.pnc{background:var(--surface2);color:var(--text-muted);border-color:var(--border);}

/* ── Event cards ─────────────────────────────────────────────────── */
.ev-cards{display:flex;gap:8px;flex-wrap:wrap;}
.ev-card{border-radius:6px;padding:8px 12px;min-width:120px;}

/* ── Ideogram ────────────────────────────────────────────────────── */
#ideogram-wrap { cursor:pointer; }
#ideogram-wrap svg { max-width:none; }
.ideo-selected{font-family:var(--mono);color:var(--navy);font-weight:600;}

/* ── Chromosome detail (single) ──────────────────────────────────── */
#chrom-detail-wrap{
  width:100%;overflow-x:auto;padding:6px 4px 12px;
  background:white;min-height:60px;
}
.brush-hint{font-size:10px;color:var(--text-muted);margin-top:4px;padding-left:4px;}

/* ── Syndrome table ──────────────────────────────────────────────── */
.stbl{width:100%;border-collapse:collapse;font-size:12px;}
.stbl thead th{background:var(--surface2);padding:7px 10px;font-size:10px;
               font-weight:700;text-transform:uppercase;letter-spacing:.05em;
               color:var(--text-muted);border-bottom:1px solid var(--border-s);text-align:left;}
.stbl tbody td{padding:8px 10px;border-bottom:1px solid var(--border);vertical-align:middle;}
.stbl tbody tr:last-child td{border-bottom:none;}
.stbl tbody tr:hover{background:var(--surface2);}
.stbl tbody tr{cursor:pointer;}
.catrow td{background:var(--surface2)!important;font-size:9px;font-weight:700;
           letter-spacing:.07em;text-transform:uppercase;color:var(--text-muted);
           padding:5px 10px!important;border-top:1px solid var(--border)!important;
           cursor:default;}
.sname{font-weight:600;font-size:12px;}
.ssub{font-size:10px;color:var(--text-muted);}

/* ── CNV tabs ────────────────────────────────────────────────────── */
.chrom-tabs{display:flex;gap:4px;flex-wrap:wrap;margin-bottom:10px;}
.ctab{padding:4px 10px;border:1px solid var(--border);border-radius:6px;
      font-size:11px;font-weight:600;cursor:pointer;background:var(--surface);
      color:var(--text-muted);transition:all .12s;}
.ctab:hover{border-color:var(--navy-b);color:var(--navy);}
.ctab.active{background:var(--navy);color:#fff;border-color:var(--navy);}
.ctab.abn{border-color:var(--red-b);color:var(--red);}
.ctab.abn.active{background:var(--red);border-color:var(--red);color:#fff;}
.ctab.sus{border-color:var(--amber-b);color:var(--amber);}
.ctab.sus.active{background:var(--amber);border-color:var(--amber);color:#fff;}
.ctab.nml{color:var(--text-muted);}
.cnv-panel{display:none;}
.cnv-panel.active{display:block;}
.syn-chips{display:flex;gap:6px;flex-wrap:wrap;margin-bottom:8px;}
.region-list{font-size:11px;color:var(--text-muted);margin-bottom:6px;
             font-family:var(--mono);}

/* ── Two-col layout ──────────────────────────────────────────────── */
.cols{display:grid;grid-template-columns:380px 1fr;gap:1rem;}
@media(max-width:860px){.cols{grid-template-columns:1fr;}}

/* ── Loading spinner ─────────────────────────────────────────────── */
.loading{display:flex;align-items:center;gap:8px;color:var(--text-muted);
         font-size:12px;padding:20px 0;}
.spinner{width:16px;height:16px;border:2px solid var(--border);
         border-top-color:var(--navy);border-radius:50%;
         animation:spin .6s linear infinite;}
@keyframes spin{to{transform:rotate(360deg);}}

/* ── Footer ──────────────────────────────────────────────────────── */
.footer{margin-top:2rem;padding-top:1rem;border-top:1px solid var(--border);
        font-size:11px;color:var(--text-muted);
        display:flex;justify-content:space-between;flex-wrap:wrap;gap:6px;}
"""
