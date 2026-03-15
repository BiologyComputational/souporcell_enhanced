// ============================================================================
// report/css.rs — All CSS for the white-theme HTML report
// ============================================================================

pub fn report_css() -> &'static str {
    r#"
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
:root{
  --bg:#F8F7F4;--surface:#FFFFFF;--surface2:#F1EFE9;--border:#DDD9D0;--border2:#C8C4BB;
  --teal:#006A6A;--teal-lt:#E0F2F2;--teal-mid:#B2DCDC;
  --amber:#C75B00;--amber-lt:#FEF3E2;
  --red:#B91C1C;--red-lt:#FEE2E2;
  --green:#15803D;--green-lt:#DCFCE7;
  --ink:#0F172A;--ink2:#334155;--muted:#64748B;--muted2:#94A3B8;
  --font-body:'DM Sans',system-ui,sans-serif;
  --font-serif:'Source Serif 4',Georgia,serif;
  --font-mono:'JetBrains Mono',monospace;
  --radius:6px;--radius-lg:10px;
  --shadow:0 1px 4px rgba(0,0,0,.07),0 4px 16px rgba(0,0,0,.05);
  --shadow-lg:0 4px 24px rgba(0,0,0,.10);
}
html{scroll-behavior:smooth}
body{background:var(--bg);color:var(--ink);font-family:var(--font-body);font-size:14px;
     line-height:1.65;-webkit-font-smoothing:antialiased}

/* Layout */
.page{max-width:1380px;margin:0 auto;padding:0 32px 60px}
.two-col{display:grid;grid-template-columns:1fr 1fr;gap:16px}

/* Cover */
.cover{background:var(--teal);color:#fff;padding:48px 40px 40px;margin:0 -32px 40px;
       position:relative;overflow:hidden}
.cover::after{content:'';position:absolute;right:-60px;top:-60px;width:340px;height:340px;
              border-radius:50%;border:60px solid rgba(255,255,255,.06);pointer-events:none}
.cover::before{content:'';position:absolute;right:120px;bottom:-80px;width:220px;height:220px;
               border-radius:50%;border:40px solid rgba(255,255,255,.04);pointer-events:none}
.cover-eyebrow{font-size:11px;font-weight:600;letter-spacing:.12em;text-transform:uppercase;
               opacity:.7;margin-bottom:8px}
.cover h1{font-family:var(--font-serif);font-size:2rem;font-weight:600;
           line-height:1.2;margin-bottom:12px}
.cover h1 em{font-style:italic;opacity:.85}
.cover-meta{display:flex;flex-wrap:wrap;gap:10px;margin-top:20px}
.cover-badge{display:inline-flex;align-items:center;gap:6px;padding:5px 14px;
             border-radius:20px;font-size:11px;font-weight:600;letter-spacing:.04em;
             border:1.5px solid rgba(255,255,255,.3);background:rgba(255,255,255,.12)}
.cover-badge.pass{background:rgba(255,255,255,.2);border-color:rgba(255,255,255,.5)}
.cover-badge.warn{background:rgba(199,91,0,.35);border-color:rgba(199,91,0,.6)}

/* Section headings */
.section-label{font-size:10px;font-weight:700;letter-spacing:.12em;text-transform:uppercase;
               color:var(--teal);margin-bottom:4px}
h2{font-family:var(--font-serif);font-size:1.4rem;font-weight:600;color:var(--ink);
   margin-bottom:6px;line-height:1.3}
h2+.section-desc{font-size:13px;color:var(--muted);margin-bottom:20px;max-width:700px}
h3{font-size:.88rem;font-weight:700;color:var(--ink2);text-transform:uppercase;
   letter-spacing:.06em;margin:22px 0 9px}
.section-divider{border:none;border-top:1px solid var(--border);margin:40px 0}

/* Metric cards */
.metric-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(190px,1fr));
             gap:12px;margin-bottom:32px}
.metric-card{background:var(--surface);border:1px solid var(--border);
             border-radius:var(--radius-lg);padding:18px 20px;
             box-shadow:var(--shadow);position:relative;overflow:hidden}
.metric-card::before{content:'';position:absolute;top:0;left:0;right:0;height:3px;
                     background:var(--teal)}
.metric-card.amber::before{background:var(--amber)}
.metric-card.red::before  {background:var(--red)}
.metric-card.green::before{background:var(--green)}
.mc-label{font-size:10px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;
          color:var(--muted);margin-bottom:6px}
.mc-value{font-size:1.75rem;font-weight:700;color:var(--ink);line-height:1;
          font-family:var(--font-mono)}
.metric-card.amber .mc-value{color:var(--amber)}
.metric-card.red   .mc-value{color:var(--red)}
.metric-card.green .mc-value{color:var(--green)}
.metric-card.teal  .mc-value{color:var(--teal)}
.mc-sub{font-size:11px;color:var(--muted);margin-top:5px}

/* Generic card */
.card{background:var(--surface);border:1px solid var(--border);
      border-radius:var(--radius-lg);padding:22px 24px;
      box-shadow:var(--shadow);margin-bottom:16px}
.card-title{font-size:11px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;
            color:var(--muted);margin-bottom:12px}

/* Tables */
.table-wrap{overflow-x:auto;margin-bottom:16px}
table{width:100%;border-collapse:collapse;font-size:13px}
thead tr{background:var(--surface2);border-bottom:2px solid var(--border2)}
th{text-align:left;padding:9px 14px;font-size:10.5px;font-weight:700;
   text-transform:uppercase;letter-spacing:.07em;color:var(--muted);white-space:nowrap}
td{padding:8px 14px;border-bottom:1px solid var(--border);color:var(--ink2);
   vertical-align:middle}
tbody tr:hover{background:#F9F8F5}
tbody tr:last-child td{border-bottom:none}
tfoot tr{background:var(--surface2)}
.mono{font-family:var(--font-mono);font-size:12px}
.num{text-align:right;font-family:var(--font-mono)}
.td-good{color:var(--green);font-weight:700}
.td-warn{color:var(--amber);font-weight:700}
.td-bad {color:var(--red); font-weight:700}
.td-teal{color:var(--teal);font-weight:600}
.pill{display:inline-block;padding:2px 8px;border-radius:10px;font-size:10px;font-weight:700}
.pill-green{background:var(--green-lt);color:var(--green)}
.pill-amber{background:var(--amber-lt);color:var(--amber)}
.pill-red  {background:var(--red-lt);  color:var(--red)}
.pill-teal {background:var(--teal-lt); color:var(--teal)}

/* Pipeline graph */
#graph-wrap{display:flex;background:var(--surface);border:1px solid var(--border);
            border-radius:var(--radius-lg);overflow:hidden;height:620px;
            margin-bottom:24px;box-shadow:var(--shadow)}
#graph-sidebar{width:230px;flex-shrink:0;border-right:1px solid var(--border);
               display:flex;flex-direction:column;overflow:hidden;background:var(--surface2)}
#graph-stage-list{flex:1;overflow-y:auto;padding:4px 0}
#graph-svg{width:100%;height:100%;display:block}
#detail-panel{display:none;position:absolute;top:12px;right:12px;width:290px;
              background:var(--surface);border:1px solid var(--border2);
              border-radius:var(--radius-lg);padding:16px;box-shadow:var(--shadow-lg);
              z-index:10;max-height:94%;overflow-y:auto}
#step-desc{display:none;position:absolute;bottom:10px;left:10px;right:10px;
           background:rgba(255,255,255,.95);border:1px solid var(--border2);
           border-radius:var(--radius);padding:9px 13px;font-size:12px;
           color:var(--muted);z-index:5;backdrop-filter:blur(4px)}
.sb-btn{flex:1;font-size:11px;font-weight:600;padding:6px;border:1px solid var(--border2);
        border-radius:var(--radius);background:var(--surface);color:var(--ink2);
        cursor:pointer;transition:background .15s}
.sb-btn:hover{background:var(--teal-lt);color:var(--teal)}

/* Interpretation cards */
.interp-grid{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-bottom:24px}
.interp-card{background:var(--surface);border:1px solid var(--border);
             border-radius:var(--radius-lg);padding:18px 20px;box-shadow:var(--shadow)}
.ic-title{font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.09em;
          color:var(--teal);margin-bottom:6px}
.interp-card p{font-size:12.5px;color:var(--muted);line-height:1.65}
.ic-verdict{margin-top:10px;padding:7px 12px;border-radius:var(--radius);
            font-size:11.5px;font-weight:600;border-left:3px solid}
.ic-verdict.good{background:var(--green-lt);color:var(--green);border-color:var(--green)}
.ic-verdict.warn{background:var(--amber-lt);color:var(--amber);border-color:var(--amber)}
.ic-verdict.info{background:var(--teal-lt); color:var(--teal); border-color:var(--teal)}

/* Experiment / hypothesis boxes */
.hypo-box{background:var(--teal-lt);border:1px solid var(--teal-mid);
          border-left:4px solid var(--teal);border-radius:var(--radius);
          padding:14px 18px;margin-bottom:14px}
.hypo-label{font-size:10px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;
            color:var(--teal);margin-bottom:4px}
.hypo-box p{font-size:13px;color:var(--ink2);line-height:1.6;font-family:var(--font-serif)}
.pipeline-steps{display:flex;flex-direction:column;gap:5px;margin:10px 0}
.pipeline-step{display:flex;align-items:flex-start;gap:12px}
.ps-num{width:24px;height:24px;border-radius:50%;background:var(--teal);color:#fff;
        font-size:11px;font-weight:700;display:flex;align-items:center;
        justify-content:center;flex-shrink:0;margin-top:2px}
.ps-name{font-weight:700;font-size:13px;color:var(--ink)}
.ps-desc{font-size:12px;color:var(--muted)}

/* SVG diagnostic plots */
.plots-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(420px,1fr));
            gap:16px;margin-bottom:24px}
.plot-card{background:var(--surface);border:1px solid var(--border);
           border-radius:var(--radius-lg);overflow:hidden;box-shadow:var(--shadow)}
.plot-card .plot-title{padding:10px 16px;font-size:11px;font-weight:700;
                        letter-spacing:.08em;text-transform:uppercase;color:var(--muted);
                        border-bottom:1px solid var(--border);background:var(--surface2)}
.plot-card img{width:100%;display:block}

/* Methods / references */
.ref-list{list-style:none;padding:0}
.ref-list li{padding:6px 0;border-bottom:1px solid var(--border);
             font-size:12px;color:var(--muted);line-height:1.55}
.ref-list li:last-child{border-bottom:none}
.ref-list strong{color:var(--ink2)}

/* Footer */
.footer{border-top:2px solid var(--border);margin-top:40px;padding-top:20px;
        display:flex;justify-content:space-between;align-items:flex-start;
        flex-wrap:wrap;gap:12px}
.footer p{font-size:11px;color:var(--muted);line-height:1.6}

@media(max-width:900px){
  .two-col,.interp-grid{grid-template-columns:1fr}
  .cover{padding:32px 24px 28px}
  .page{padding:0 16px 40px}
}
"#
}
