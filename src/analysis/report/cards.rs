// ============================================================================
// report/cards.rs — Executive-summary metric card grid (§2)
// ============================================================================

use crate::analysis::report::ReportData;
use crate::analysis::plots::PlotData;

pub fn build_metric_cards(
    data:        &ReportData,
    _plot_data:  &PlotData,
    at_best:     usize,
    conv_pct:    usize,
    max_imb:     usize,
    median_conf: f32,
    ambig_pct:   usize,
) -> String {
    let k        = data.params.num_clusters;
    let total: usize = data.cells_per_cluster.iter().sum();
    let loci_pct = data.loci_used as f32 / data.loci_raw.max(1) as f32 * 100.0;

    let total_secs = data.total_time.as_secs_f32();
    let time_str   = if total_secs < 60.0 {
        format!("{:.1}s", total_secs)
    } else {
        format!("{:.0}m {:.0}s", total_secs / 60.0, total_secs % 60.0)
    };

    let conv_cls = if conv_pct >= 70 { "green" } else if conv_pct >= 40 { "amber" } else { "red" };
    let imb_cls  = if max_imb  < 20  { "green" } else if max_imb  < 40 { "teal"  } else { "amber" };
    let conf_cls = if ambig_pct < 5  { "green" } else { "amber" };

    fn card(label: &str, value: &str, sub: &str, cls: &str) -> String {
        format!(
            r#"<div class="metric-card {cls}">
  <div class="mc-label">{label}</div>
  <div class="mc-value">{value}</div>
  <div class="mc-sub">{sub}</div>
</div>"#
        )
    }

    let mut out = String::new();
    out += &card("Cells demultiplexed",
                 &format!("{}", total),
                 "total barcodes assigned", "teal");
    out += &card("SNP loci used",
                 &format!("{}", data.loci_used),
                 &format!("{:.1}% of {} raw loci", loci_pct, data.loci_raw), "");
    out += &card("Donor clusters k",
                 &format!("{}", k),
                 &format!("{} &middot; {} restarts",
                          data.params.clustering_method.name(), data.total_restarts), "");
    out += &card("Global best log-prob",
                 &format!("{:.0}", data.global_best_lp),
                 &format!("{}/{} restarts at optimum", at_best, data.total_restarts), "teal");
    out += &card("Convergence stability",
                 &format!("{}%", conv_pct),
                 &format!("{} of {} restarts agreed", at_best, data.total_restarts), conv_cls);
    out += &card("Max cluster imbalance",
                 &format!("{}%", max_imb),
                 "deviation from equal split", imb_cls);
    out += &card("Median cell confidence",
                 &format!("{:.0}", median_conf),
                 &format!("{}% cells ambiguous (&lt;10 log-units)", ambig_pct), conf_cls);
    out += &card("Wall-clock time",
                 &time_str,
                 &format!("{} threads &middot; seed {}",
                          data.params.threads, data.params.seed), "");
    out
}
