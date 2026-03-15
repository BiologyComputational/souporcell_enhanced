// ============================================================================
// report/graph.rs — D3.js pipeline graph JSON builder (§4)
//
// Loads graph_template.json at compile time via include_str! and substitutes
// {{PLACEHOLDER}} tokens with live run data.
// ============================================================================

use crate::analysis::report::ReportData;
use crate::analysis::plots::PlotData;

pub fn build_graph_json(data: &ReportData, plot_data: &PlotData) -> String {
    let k           = data.params.num_clusters;
    let n_cells     = data.total_cells;
    let n_loci_raw  = data.loci_raw;
    let n_loci_used = data.loci_used;
    let zero_pct    = data.zero_cov_frac * 100.0;
    let restarts    = data.total_restarts;
    let threads     = data.params.threads;
    let best_lp     = data.global_best_lp;
    let method      = data.params.clustering_method.name().to_uppercase();
    let init        = format!("{:?}", data.params.initialization_strategy);

    let at_best = plot_data.restart_losses.iter()
        .filter(|&&(_, _, l)| (l - plot_data.global_best).abs() < 1.0)
        .count();
    let conv_quality = if restarts > 0 {
        (at_best as f32 / restarts as f32 * 100.0) as usize
    } else { 0 };

    let total: usize = data.cells_per_cluster.iter().sum();
    let expected     = if k > 0 { total / k } else { 1 };
    let max_imbalance = data.cells_per_cluster.iter()
        .map(|&c| ((c as isize - expected as isize).abs() as f32
                   / expected as f32 * 100.0) as usize)
        .max().unwrap_or(0);

    let load_ms   = data.load_time.as_millis();
    let clust_ms  = data.cluster_time.as_millis();
    let total_ms  = data.total_time.as_millis();
    let n_ts      = plot_data.temp_step_losses.len();
    let rpt       = (restarts as f32 / threads as f32).ceil() as usize;
    let n_plots   = data.plot_files.len();
    let loci_f    = n_loci_raw.saturating_sub(n_loci_used);
    let flt_pct   = loci_f as f32 / n_loci_raw.max(1) as f32 * 100.0;

    include_str!("graph_template.json")
        .replace("{{K}}",             &k.to_string())
        .replace("{{RESTARTS}}",      &restarts.to_string())
        .replace("{{THREADS}}",       &threads.to_string())
        .replace("{{SEED}}",          &data.params.seed.to_string())
        .replace("{{METHOD}}",        &method)
        .replace("{{INIT}}",          &init)
        .replace("{{LOAD_MS}}",       &load_ms.to_string())
        .replace("{{CLUST_MS}}",      &clust_ms.to_string())
        .replace("{{TOTAL_MS}}",      &total_ms.to_string())
        .replace("{{N_CELLS}}",       &n_cells.to_string())
        .replace("{{N_LOCI_RAW}}",    &n_loci_raw.to_string())
        .replace("{{N_LOCI_USED}}",   &n_loci_used.to_string())
        .replace("{{ZERO_PCT}}",      &format!("{:.1}", zero_pct))
        .replace("{{LOCI_FILTERED}}", &loci_f.to_string())
        .replace("{{FILTER_PCT}}",    &format!("{:.1}", flt_pct))
        .replace("{{MIN_REF}}",       &data.params.min_ref.to_string())
        .replace("{{MIN_ALT}}",       &data.params.min_alt.to_string())
        .replace("{{AT_BEST}}",       &at_best.to_string())
        .replace("{{CONV_QUALITY}}",  &conv_quality.to_string())
        .replace("{{BEST_LP}}",       &format!("{:.0}", best_lp))
        .replace("{{N_TS}}",          &n_ts.to_string())
        .replace("{{RPT}}",           &rpt.to_string())
        .replace("{{MAX_IMBALANCE}}", &max_imbalance.to_string())
        .replace("{{N_PLOTS}}",       &n_plots.to_string())
}
