// ============================================================================
// report/tables.rs — Detailed results tables (§5)
//
// Table 1: per-cluster statistics (cells, fraction, imbalance, confidence)
// Table 2: restart landscape summary (global/mean/worst LP, convergence)
// Table 3: annealing trace summary   (temp-steps, LP range)
// Table 4: per-thread efficiency
// Table 5: locus QC accounting
// ============================================================================

use crate::analysis::report::ReportData;
use crate::analysis::plots::PlotData;

pub fn build_results_tables(
    data:        &ReportData,
    plot_data:   &PlotData,
    _at_best:    usize,
    conv_pct:    usize,
    max_imb:     usize,
    median_conf: f32,
    ambig_pct:   usize,
) -> String {
    let k      = data.params.num_clusters;
    let total: usize = data.cells_per_cluster.iter().sum();
    let exp    = if k > 0 { total / k } else { 1 };
    let global_lp = plot_data.global_best;

    // ── Table 1: per-cluster stats ──────────────────────────────────────────
    let cluster_rows: String = data.cells_per_cluster.iter().enumerate()
        .map(|(i, &count)| {
            let frac = count as f32 / total.max(1) as f32 * 100.0;
            let dev  = ((count as isize - exp as isize).abs() as f32
                        / exp as f32 * 100.0) as usize;

            let cl_confs: Vec<f32> = plot_data.cell_assignments.iter()
                .filter(|(cl, lps)| *cl == i && lps.len() >= 2)
                .map(|(_, lps)| {
                    let mut s = lps.clone();
                    s.sort_by(|a, b| b.total_cmp(a));
                    s[0] - s[1]
                })
                .collect();

            let (med_c, min_c, max_c) = if cl_confs.is_empty() {
                (0.0f32, 0.0f32, 0.0f32)
            } else {
                let mut sc = cl_confs.clone();
                sc.sort_by(|a, b| a.total_cmp(b));
                (sc[sc.len() / 2], sc[0], *sc.last().unwrap())
            };

            let cl_ambig  = cl_confs.iter().filter(|&&c| c < 10.0).count();
            let ambig_cl  = if count > 0 { cl_ambig * 100 / count } else { 0 };
            let conf_cls  = if med_c > 50.0 { "td-good" } else if med_c > 15.0 { "td-teal" } else { "td-warn" };
            let dev_cls   = if dev  < 20    { "td-good" } else if dev    < 40   { "td-warn" } else { "td-bad" };

            format!(
                "<tr>\
                 <td class='mono'>Cluster {i}</td>\
                 <td class='num'>{count}</td>\
                 <td class='num'>{frac:.1}%</td>\
                 <td class='num {dev_cls}'>{dev}%</td>\
                 <td class='num {conf_cls}'>{med_c:.0}</td>\
                 <td class='num'>{min_c:.0}</td>\
                 <td class='num'>{max_c:.0}</td>\
                 <td class='num'>{ambig_cl}%</td>\
                 </tr>"
            )
        })
        .collect();

    // ── Table 2: restart landscape ──────────────────────────────────────────
    let n_restarts = data.total_restarts;
    let within_1:  usize = plot_data.restart_losses.iter()
        .filter(|&&(_, _, l)| (l - global_lp).abs() < 1.0).count();
    let within_5:  usize = plot_data.restart_losses.iter()
        .filter(|&&(_, _, l)| (l - global_lp).abs() < 5.0).count();
    let within_10: usize = plot_data.restart_losses.iter()
        .filter(|&&(_, _, l)| (l - global_lp).abs() < 10.0).count();
    let worst_lp  = plot_data.restart_losses.iter()
        .map(|&(_, _, l)| l).fold(f32::INFINITY, f32::min);
    let mean_lp: f32 = if n_restarts > 0 {
        plot_data.restart_losses.iter().map(|&(_, _, l)| l).sum::<f32>()
            / n_restarts as f32
    } else { 0.0 };

    // ── Table 3: annealing trace ────────────────────────────────────────────
    let n_ts    = plot_data.temp_step_losses.len();
    let ts_per_r = if n_restarts > 0 { n_ts / n_restarts } else { 0 };
    let (ts_min, ts_max) = plot_data.temp_step_losses.iter()
        .fold((f32::INFINITY, f32::NEG_INFINITY),
              |(mn, mx), &(_, _, _, _, l)| (mn.min(l), mx.max(l)));
    let ts_min_s  = if ts_min.is_infinite() { 0.0 } else { ts_min };
    let ts_max_s  = if ts_max.is_infinite() { 0.0 } else { ts_max };
    let ts_range  = (ts_max_s - ts_min_s).abs();

    // ── Table 4: per-thread ─────────────────────────────────────────────────
    let thread_rows: String = (0..data.params.threads).map(|tid| {
        let tl: Vec<f32> = plot_data.restart_losses.iter()
            .filter(|&&(t, _, _)| t == tid).map(|&(_, _, l)| l).collect();
        let n_t    = tl.len();
        let best_t = tl.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        let at_g   = tl.iter().filter(|&&l| (l - global_lp).abs() < 1.0).count();
        let eff    = if n_t > 0 { at_g * 100 / n_t } else { 0 };
        let eff_c  = if eff >= 70 { "td-good" } else if eff >= 40 { "td-warn" } else { "td-bad" };
        format!(
            "<tr>\
             <td class='num'>Thread {tid}</td>\
             <td class='num'>{n_t}</td>\
             <td class='num mono'>{best_t:.2}</td>\
             <td class='num {eff_c}'>{eff}%</td>\
             </tr>"
        )
    }).collect();

    let loci_filtered = data.loci_raw.saturating_sub(data.loci_used);
    let loci_pct      = data.loci_used as f32 / data.loci_raw.max(1) as f32 * 100.0;
    let filter_pct    = loci_filtered as f32 / data.loci_raw.max(1) as f32 * 100.0;

    format!(r#"<div class="section-label">Detailed Results</div>
<h2>Comprehensive Output Tables</h2>
<p class="section-desc">
  All values derived dynamically from this run&apos;s output &mdash; no data is hardcoded.
</p>

<h3>Table 1 &middot; Per-Cluster Statistics</h3>
<div class="table-wrap"><table>
  <thead><tr>
    <th>Cluster</th><th class="num">Cells</th><th class="num">Fraction</th>
    <th class="num">Imbalance</th><th class="num">Median conf.</th>
    <th class="num">Min conf.</th><th class="num">Max conf.</th>
    <th class="num">Ambiguous</th>
  </tr></thead>
  <tbody>{cluster_rows}</tbody>
  <tfoot><tr>
    <td>Total / summary</td>
    <td class="num">{total}</td>
    <td class="num">100%</td>
    <td class="num">{max_imb}% max</td>
    <td class="num">{median_conf:.0}</td>
    <td class="num">&mdash;</td>
    <td class="num">&mdash;</td>
    <td class="num">{ambig_pct}%</td>
  </tr></tfoot>
</table></div>

<h3>Table 2 &middot; Restart Landscape Summary</h3>
<div class="table-wrap"><table>
  <thead><tr><th>Metric</th><th class="num">Value</th><th>Interpretation</th></tr></thead>
  <tbody>
    <tr><td>Total restarts executed</td>
        <td class="num">{n_restarts}</td><td>Configured via --restarts</td></tr>
    <tr><td>Global best log-prob</td>
        <td class="num mono">{global_lp:.2}</td>
        <td>Highest observed mixture log-likelihood</td></tr>
    <tr><td>Mean log-prob across restarts</td>
        <td class="num mono">{mean_lp:.2}</td>
        <td>&Delta; = {delta_mean:.2} from global best</td></tr>
    <tr><td>Worst log-prob (deepest local opt.)</td>
        <td class="num mono">{worst_lp:.2}</td>
        <td>Deepest local optimum reached</td></tr>
    <tr><td>Restarts within 1 unit of global</td>
        <td class="num">{within_1} ({pct1}%)</td>
        <td>Convergence stability (target &ge;70%)</td></tr>
    <tr><td>Restarts within 5 units of global</td>
        <td class="num">{within_5} ({pct5}%)</td>
        <td>Near-global solutions</td></tr>
    <tr><td>Restarts within 10 units of global</td>
        <td class="num">{within_10} ({pct10}%)</td>
        <td>Broad convergence basin</td></tr>
    <tr><td>Convergence stability</td>
        <td class="num {conv_cls}">{conv_pct}%</td>
        <td>{conv_verdict}</td></tr>
  </tbody>
</table></div>

<h3>Table 3 &middot; Annealing Trace Summary</h3>
<div class="table-wrap"><table>
  <thead><tr><th>Metric</th><th class="num">Value</th><th>Notes</th></tr></thead>
  <tbody>
    <tr><td>Total temperature steps recorded</td>
        <td class="num">{n_ts}</td><td>Across all restarts</td></tr>
    <tr><td>Temperature steps per restart (avg)</td>
        <td class="num">{ts_per_r}</td><td>Annealing schedule length</td></tr>
    <tr><td>Min log-prob during annealing</td>
        <td class="num mono">{ts_min_s:.2}</td><td>Hottest-temperature state</td></tr>
    <tr><td>Max log-prob during annealing</td>
        <td class="num mono">{ts_max_s:.2}</td><td>Coolest-temperature converged value</td></tr>
    <tr><td>Log-prob improvement (max&minus;min)</td>
        <td class="num">{ts_range:.2}</td><td>Landscape traversal depth</td></tr>
  </tbody>
</table></div>

<h3>Table 4 &middot; Per-Thread Efficiency</h3>
<div class="table-wrap"><table>
  <thead><tr>
    <th>Thread</th><th class="num">Restarts run</th>
    <th class="num">Best log-prob</th><th class="num">At global opt.</th>
  </tr></thead>
  <tbody>{thread_rows}</tbody>
</table></div>

<h3>Table 5 &middot; Locus QC Accounting</h3>
<div class="table-wrap"><table>
  <thead><tr>
    <th>Stage</th><th class="num">Count</th><th class="num">Fraction</th><th>Notes</th>
  </tr></thead>
  <tbody>
    <tr><td>Raw loci loaded</td>
        <td class="num">{loci_raw}</td><td class="num">100%</td>
        <td>All positions in ref + alt matrices</td></tr>
    <tr><td>Loci passing QC filters</td>
        <td class="num">{loci_used}</td><td class="num">{loci_pct:.1}%</td>
        <td>min_ref / min_alt filter applied</td></tr>
    <tr><td>Loci filtered out</td>
        <td class="num td-warn">{loci_filtered}</td>
        <td class="num td-warn">{filter_pct:.1}%</td>
        <td>Too few informative reads at these positions</td></tr>
    <tr><td>Matrix sparsity (post-filter)</td>
        <td class="num">{sparsity:.1}%</td><td class="num">&mdash;</td>
        <td>Zero-coverage locus&times;cell fraction</td></tr>
  </tbody>
</table></div>
"#,
        cluster_rows  = cluster_rows,
        total         = total,
        max_imb       = max_imb,
        median_conf   = median_conf,
        ambig_pct     = ambig_pct,
        n_restarts    = n_restarts,
        global_lp     = global_lp,
        mean_lp       = mean_lp,
        delta_mean    = (global_lp - mean_lp).abs(),
        worst_lp      = worst_lp,
        within_1      = within_1,  pct1  = within_1  * 100 / n_restarts.max(1),
        within_5      = within_5,  pct5  = within_5  * 100 / n_restarts.max(1),
        within_10     = within_10, pct10 = within_10 * 100 / n_restarts.max(1),
        conv_pct      = conv_pct,
        conv_cls      = if conv_pct >= 70 { "td-good" } else { "td-warn" },
        conv_verdict  = if conv_pct >= 70 { "Strong: global optimum reliably found" }
                        else { "Moderate: consider --restarts 200+" },
        n_ts          = n_ts,
        ts_per_r      = ts_per_r,
        ts_min_s      = ts_min_s,
        ts_max_s      = ts_max_s,
        ts_range      = ts_range,
        thread_rows   = thread_rows,
        loci_raw      = data.loci_raw,
        loci_used     = data.loci_used,
        loci_pct      = loci_pct,
        loci_filtered = loci_filtered,
        filter_pct    = filter_pct,
        sparsity      = data.zero_cov_frac * 100.0,
    )
}
