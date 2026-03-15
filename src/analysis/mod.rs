// ============================================================================
// analysis/mod.rs — Post-clustering analysis layer
// ============================================================================
//
// Everything that runs after the main clustering loop: diagnostic plots
// and the HTML run report.  These modules consume the clustering output
// (log-probabilities, cluster centres) and produce artefacts for human
// inspection — they never feed back into the algorithm.
//
// Dependency rule: analysis may import domain, config, and infra.  It must NOT
// import core (no algorithm logic here).  Modules within this layer are
// independent of each other (plots does not import report, etc.).
//
// Modules:
//   plots          — 6 SVG diagnostic plots accumulated during the run via PlotData
//   report/        — Self-contained HTML run report (multi-file sub-module)
//     mod.rs       — ReportData struct, write_report(), build_html() orchestrator
//     css.rs       — All CSS for the white-theme report
//     cards.rs     — §2 executive metric card grid
//     experiment.rs— §3 experiment design, hypotheses, pipeline overview
//     graph.rs     — §4 D3 graph JSON builder (reads report/graph_template.json)
//     graph_js.rs  — §4 D3 <script> block
//     tables.rs    — §5 detailed results tables (5 tables, all live data)
//     svgs.rs      — §6 inline SVG balance + confidence charts
//     plot_embeds.rs— §7 embedded SVG file cards
//     interpretation.rs — §8 PhD-level per-metric commentary
//     params.rs    — §9 full run parameters table
//     methods.rs   — §10 algorithmic methods + references
//     utils.rs     — shared colour constants + helper functions
// ============================================================================

pub mod plots;
pub mod report;

// Re-export primary public types and entry points.
#[allow(unused_imports)]
pub use plots::{PlotData, PlotFlags, write_plots};
#[allow(unused_imports)]
pub use report::{ReportData, write_report};
