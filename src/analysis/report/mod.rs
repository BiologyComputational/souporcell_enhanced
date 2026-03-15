// ============================================================================
// analysis/report/mod.rs — HTML report module root
// ============================================================================
//
// Sub-modules (one responsibility each):
//   css             — all CSS for the white-theme report
//   cards           — §2  executive-summary metric card grid
//   experiment      — §3  experiment design, hypotheses, pipeline overview
//   graph           — §4  D3.js graph JSON builder (reads graph_template.json)
//   graph_js        — §4  D3.js <script> block
//   tables          — §5  detailed results tables (5 tables, all live data)
//   svgs            — §6  inline SVG balance + confidence charts
//   plot_embeds     — §7  embedded SVG file cards
//   interpretation  — §8  PhD-level per-metric commentary
//   params          — §9  full run parameters table
//   methods         — §10 algorithmic methods + references
//   utils           — shared colour constants + helper functions
// ============================================================================

pub mod cards;
pub mod css;
pub mod experiment;
pub mod graph;
pub mod graph_js;
pub mod interpretation;
pub mod methods;
pub mod params;
pub mod plot_embeds;
pub mod svgs;
pub mod tables;
pub mod utils;

use std::fs;
use std::time::Duration;

use crate::config::params::Params;
use crate::analysis::plots::PlotData;

// ── Public data carrier ───────────────────────────────────────────────────────

pub struct ReportData {
    pub params:            Params,
    pub total_cells:       usize,
    pub loci_raw:          usize,
    pub loci_used:         usize,
    pub zero_cov_frac:     f32,
    pub load_time:         Duration,
    pub cluster_time:      Duration,
    pub total_time:        Duration,
    pub cells_per_cluster: Vec<usize>,
    pub total_restarts:    usize,
    pub global_best_lp:    f32,
    pub plot_files:        Vec<String>,
    pub plot_dir:          String,
}

// ── Public entry point ────────────────────────────────────────────────────────

pub fn write_report(data: &ReportData, plot_data: &PlotData) -> String {
    let path = format!("{}/souporcell_report.html", data.plot_dir);
    let html = build_html(data, plot_data);
    fs::write(&path, &html).unwrap_or_else(|e| {
        eprintln!("Warning: could not write report: {}", e);
    });
    path
}

// ── HTML orchestrator ─────────────────────────────────────────────────────────

fn build_html(data: &ReportData, plot_data: &PlotData) -> String {
    // ── Pre-compute shared derived quantities ──────────────────────────────
    let k           = data.params.num_clusters;
    let total_cells: usize = data.cells_per_cluster.iter().sum();
    let expected    = if k > 0 { total_cells / k } else { 1 };

    let at_best = plot_data.restart_losses.iter()
        .filter(|&&(_, _, l)| (l - plot_data.global_best).abs() < 1.0)
        .count();
    let conv_pct = if data.total_restarts > 0 {
        at_best * 100 / data.total_restarts
    } else { 0 };
    let conv_ok  = conv_pct >= 70;

    let max_imb = data.cells_per_cluster.iter()
        .map(|&c| ((c as isize - expected as isize).abs() as f32
                   / expected as f32 * 100.0) as usize)
        .max().unwrap_or(0);
    let bal_ok  = max_imb < 50;

    let (median_conf, ambig_pct) =
        utils::compute_conf_stats(plot_data, total_cells);

    let run_id = format!("k{}_r{}_s{}",
                         k, data.total_restarts, data.params.seed);
    let ts     = utils::chrono_now();

    // ── Delegate to sub-modules ────────────────────────────────────────────
    let css_block      = css::report_css();
    let metric_cards   = cards::build_metric_cards(
                             data, plot_data, at_best, conv_pct,
                             max_imb, median_conf, ambig_pct);
    let experiment_sec = experiment::build_experiment_section(
                             data, plot_data, at_best, conv_pct,
                             max_imb, median_conf);
    let graph_json     = graph::build_graph_json(data, plot_data);
    let graph_script   = graph_js::build_graph_script(&graph_json, conv_ok, bal_ok);
    let results_tables = tables::build_results_tables(
                             data, plot_data, at_best, conv_pct,
                             max_imb, median_conf, ambig_pct);
    let balance_svg    = svgs::build_balance_svg(data);
    let confidence_svg = svgs::build_confidence_svg(plot_data);
    let plot_embeds    = plot_embeds::build_plot_embeds(data);
    let interp_guide   = interpretation::build_interpretation_guide(
                             data, plot_data, conv_pct,
                             max_imb, median_conf, ambig_pct);
    let params_table   = params::build_params_table(data);
    let methods_sec    = methods::build_methods_section(data);

    let conv_badge = if conv_ok { " pass" } else { " warn" };
    let conv_flag  = if conv_ok { "" } else { " \u{26a0}" };

    format!(r####"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>souporcell v2.6 &middot; Demultiplexing Report [{run_id}]</title>
<link rel="preconnect" href="https://fonts.googleapis.com"/>
<link href="https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,wght@0,300;0,400;0,600;0,700;1,400&family=JetBrains+Mono:wght@400;600&family=DM+Sans:wght@400;500;600;700&display=swap" rel="stylesheet"/>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
<style>{css_block}</style>
</head>
<body>

<!-- §1 COVER -->
<div class="cover">
  <div class="page" style="padding-top:0;padding-bottom:0">
    <div class="cover-eyebrow">Single-Cell RNA-seq Genotype Demultiplexing &middot; souporcell v2.6</div>
    <h1>Donor Demultiplexing Report<br/><em>Execution &amp; Quality Assessment</em></h1>
    <div style="font-size:13px;opacity:.8;max-width:640px;margin-top:8px">
      Genetic-barcode-free pooled-sample decomposition using sparse allelic read evidence
      at common SNP loci.
      k&nbsp;=&nbsp;{k} donors &middot; {method} algorithm &middot;
      {restarts} restarts &middot; seed {seed}
    </div>
    <div class="cover-meta">
      <span class="cover-badge pass">&#x2713; Run completed</span>
      <span class="cover-badge">Run ID: {run_id}</span>
      <span class="cover-badge">k = {k}</span>
      <span class="cover-badge">{method}</span>
      <span class="cover-badge{conv_badge}">Convergence {conv_pct}%{conv_flag}</span>
      <span class="cover-badge">Generated: {ts}</span>
    </div>
  </div>
</div>

<div class="page">

<!-- §2 METRICS -->
<div class="section-label">Executive Summary</div>
<h2>Run Metrics at a Glance</h2>
<p class="section-desc">
  All values computed dynamically from this run&apos;s output &mdash; no hardcoded numbers.
</p>
<div class="metric-grid">{metric_cards}</div>

<hr class="section-divider"/>

<!-- §3 EXPERIMENT DESIGN -->
{experiment_sec}

<hr class="section-divider"/>

<!-- §4 PIPELINE GRAPH -->
<div class="section-label">Interactive Pipeline Trace</div>
<h2>Pipeline Execution Graph</h2>
<p class="section-desc">
  Click nodes for stage details. Use <strong>Step mode</strong> to walk sequentially.
  <strong>Drag</strong> to rearrange; <strong>scroll</strong> to zoom.
</p>
<div id="graph-wrap">
  <div id="graph-sidebar">
    <div style="padding:12px 14px;border-bottom:1px solid var(--border)">
      <div style="font-size:13px;font-weight:700;color:var(--ink)">Pipeline stages</div>
      <div style="font-size:11px;color:var(--muted);margin-top:2px">
        11 stages &middot; 13 data flows
      </div>
    </div>
    <div style="padding:8px;border-bottom:1px solid var(--border);display:flex;gap:6px">
      <button class="sb-btn" onclick="graphShowAll()">Show all</button>
      <button class="sb-btn" onclick="graphStartStep()">Step mode</button>
      <button class="sb-btn" style="flex:0;padding:6px 10px"
              onclick="graphToggleFS()" title="Fullscreen">&#x26F6;</button>
    </div>
    <div id="graph-step-controls"
         style="display:none;padding:7px 8px;border-bottom:1px solid var(--border);
                align-items:center;gap:6px">
      <button class="sb-btn" style="flex:0;padding:4px 12px"
              onclick="graphPrev()">&#8592;</button>
      <span id="graph-step-label"
            style="flex:1;text-align:center;font-size:11px;color:var(--muted)"></span>
      <button class="sb-btn" style="flex:0;padding:4px 12px"
              onclick="graphNext()">&#8594;</button>
    </div>
    <div id="graph-stage-list"></div>
    <div style="border-top:1px solid var(--border)">
      <div style="padding:9px 14px 3px;font-size:10px;font-weight:700;color:var(--muted);
                  text-transform:uppercase;letter-spacing:.07em">Quality checklist</div>
      <div id="graph-checklist" style="padding:0 14px 10px"></div>
    </div>
  </div>
  <div style="flex:1;position:relative;background:var(--surface)">
    <svg id="graph-svg"></svg>
    <div id="detail-panel">
      <div style="display:flex;justify-content:space-between;align-items:center;
                  margin-bottom:8px">
        <span id="dp-title"
              style="font-weight:700;font-size:13px;color:var(--teal)"></span>
        <button onclick="document.getElementById('detail-panel').style.display='none'"
                style="border:none;background:none;cursor:pointer;
                       color:var(--muted);font-size:16px">&#x2715;</button>
      </div>
      <div id="dp-badge"       style="margin-bottom:8px"></div>
      <div id="dp-desc"        style="color:var(--muted);margin-bottom:10px;
                                      line-height:1.55;font-size:11.5px"></div>
      <div id="dp-rows"        style="border-top:1px solid var(--border);
                                      padding-top:8px"></div>
      <div id="dp-node-checks" style="margin-top:8px;border-top:1px solid var(--border);
                                      padding-top:8px"></div>
    </div>
    <div id="step-desc"></div>
  </div>
</div>

<hr class="section-divider"/>

<!-- §5 RESULTS TABLES -->
{results_tables}

<hr class="section-divider"/>

<!-- §6 QUALITY CHARTS -->
<div class="section-label">Assignment Quality</div>
<h2>Cluster Balance &amp; Posterior Confidence</h2>
<p class="section-desc">
  Left: cell counts per inferred donor cluster vs. ideal equal split.
  Right: per-cell posterior confidence distribution.
  Cells in the red zone (&lt;10 log-units) are assignment-ambiguous.
</p>
<div class="two-col" style="margin-bottom:24px">
  <div class="card" style="padding:0;overflow:hidden">
    <div class="card-title"
         style="padding:10px 16px;border-bottom:1px solid var(--border)">
      Cluster Cell Balance
    </div>
    {balance_svg}
  </div>
  <div class="card" style="padding:0;overflow:hidden">
    <div class="card-title"
         style="padding:10px 16px;border-bottom:1px solid var(--border)">
      Posterior Confidence Distribution
    </div>
    {confidence_svg}
  </div>
</div>

<!-- §7 SVG PLOTS -->
{plot_embeds}

<hr class="section-divider"/>

<!-- §8 INTERPRETATION -->
<div class="section-label">Scientific Interpretation</div>
<h2>PhD-Level Result Assessment</h2>
<p class="section-desc">
  Per-metric commentary grounded in single-cell demultiplexing theory.
  All verdicts derive from this run&apos;s computed values.
</p>
<div class="interp-grid">{interp_guide}</div>

<hr class="section-divider"/>

<!-- §9 PARAMETERS -->
<div class="section-label">Reproducibility</div>
<h2>Full Run Parameters</h2>
<p class="section-desc">
  Complete parameter set sufficient to reproduce this run identically.
</p>
<div class="card">{params_table}</div>

<hr class="section-divider"/>

<!-- §10 METHODS -->
{methods_sec}

<!-- §11 FOOTER -->
<div class="footer">
  <p>
    <strong>souporcell v2.6</strong> &mdash; Enhanced modular rewrite<br/>
    Author: Jahidul Arafat &middot;
    Architecture: domain / config / core / infra / analysis / analysis::report
  </p>
  <p style="text-align:right">
    Original algorithm: Heaton <em>et al.</em>, <em>Nature Methods</em> 2020<br/>
    Report generated: {ts}
  </p>
</div>

</div><!-- /page -->

{graph_script}
</body>
</html>
"####,
        run_id         = run_id,
        k              = k,
        method         = data.params.clustering_method.name().to_uppercase(),
        restarts       = data.total_restarts,
        seed           = data.params.seed,
        ts             = ts,
        conv_pct       = conv_pct,
        conv_badge     = conv_badge,
        conv_flag      = conv_flag,
        css_block      = css_block,
        metric_cards   = metric_cards,
        experiment_sec = experiment_sec,
        results_tables = results_tables,
        balance_svg    = balance_svg,
        confidence_svg = confidence_svg,
        plot_embeds    = plot_embeds,
        interp_guide   = interp_guide,
        params_table   = params_table,
        methods_sec    = methods_sec,
        graph_script   = graph_script,
    )
}
