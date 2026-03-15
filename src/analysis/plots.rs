// ============================================================================
// plots.rs — PhD-level convergence and quality diagnostics  (souporcell v2.6)
// ============================================================================
//
// Generates self-contained SVG plots from the data already computed during
// clustering.  No external dependencies — pure Rust string generation.
// All output is written to --plot_dir (default: "souporcell_plots/").
//
// Plots produced (all optional, controlled by CLI):
//
//   1. restart_landscape.svg
//      Final log-loss per restart, sorted ascending.  Horizontal line at
//      global best.  Colour-coded by thread.  Shows whether enough restarts
//      were run to reliably find the global optimum.
//
//   2. convergence_curves.svg
//      Per-restart loss trace across EM/KHM iterations, one curve per restart
//      (sampled from CONV data).  Log-scale Y axis.  Annealing temperature
//      steps marked as vertical bands.
//
//   3. annealing_profile.svg
//      Mean ± SD of loss at each of the 9 temperature steps across all
//      restarts.  Shows how the annealing schedule drives convergence.
//
//   4. cluster_balance.svg
//      Horizontal bar chart of final cell counts per cluster with percentage
//      labels.  Includes an expected-equal-split reference line.
//
//   5. posterior_confidence.svg
//      Histogram of per-cell posterior confidence = best_log_prob −
//      second_best_log_prob.  Cells near zero are ambiguous / potential
//      doublets.  Vertical line at median.
//
//   6. thread_efficiency.svg
//      Per-thread restart count vs best loss found.  Reveals thread imbalance
//      or lucky threads that found the global optimum early.
//
// PlotData struct collects everything during the run; write_plots() is called
// once at the end when params.plot_dir.is_some().
// ============================================================================

use std::fs;

// ── Colour constants (avoids # ambiguity inside format! strings) ─────────────
const C_BLUE:    &str = "#2166ac";
const C_RED:     &str = "#c0392b";
const C_GREEN:   &str = "#1b7837";
const C_GRAY:    &str = "#e0e0e0";
const C_LGRAY:   &str = "#e8e8e8";
const C_MGRAY:   &str = "#999999";
const C_DGRAY:   &str = "#333333";
const C_TEXT:    &str = "#555555";
const C_DIM:     &str = "#666666";
const C_WHITE:   &str = "white";
const C_BG:      &str = "#fafafa";
const C_AMBER:   &str = "#ffcccc";
const C_BGBAND:  &str = "#f0f4f8";
const C_TITLE:   &str = "#1a1a2e";


// ── PlotData — accumulates everything needed for all plots ────────────────────

/// Collected during the clustering run; passed to write_plots() at the end.
pub struct PlotData {
    /// (thread_id, restart_index, final_loss)
    pub restart_losses:   Vec<(usize, usize, f32)>,

    /// CONV data: (thread, epoch/restart, temp_step, iter, loss)
    /// Only populated when params.verbose = true.
    pub conv_trace:       Vec<(usize, usize, usize, usize, f32)>,

    /// EM/KHM per-temp-step summary: (thread, epoch, temp_step, iters, loss)
    /// Always populated (these come from EM/KHM lines, not CONV).
    pub temp_step_losses: Vec<(usize, usize, usize, usize, f32)>,

    /// Final cell cluster assignments with log-probability vectors.
    /// (assigned_cluster, log_probs[k])
    pub cell_assignments: Vec<(usize, Vec<f32>)>,

    /// Number of clusters k.
    pub num_clusters:     usize,

    /// Global best log-loss.
    pub global_best:      f32,
}

impl PlotData {
    pub fn new(num_clusters: usize) -> Self {
        PlotData { num_clusters, global_best: f32::NEG_INFINITY, ..Default::default() }
    }

    pub fn record_restart(&mut self, thread: usize, restart: usize, loss: f32) {
        self.restart_losses.push((thread, restart, loss));
        if loss > self.global_best { self.global_best = loss; }
    }

    #[allow(dead_code)]
    pub fn record_conv(&mut self, thread: usize, epoch: usize, temp_step: usize,
                       iter: usize, loss: f32) {
        self.conv_trace.push((thread, epoch, temp_step, iter, loss));
    }

    pub fn record_temp_step(&mut self, thread: usize, epoch: usize,
                             temp_step: usize, iters: usize, loss: f32) {
        self.temp_step_losses.push((thread, epoch, temp_step, iters, loss));
        self.conv_trace.push((thread, epoch, temp_step, iters, loss));
    }

    pub fn record_cell(&mut self, assigned: usize, log_probs: Vec<f32>) {
        self.cell_assignments.push((assigned, log_probs));
    }

    /// Minimal clone used as fallback when Arc unwrap fails (should never happen).
    pub fn clone_basic(&self) -> PlotData {
        PlotData {
            restart_losses:   self.restart_losses.clone(),
            conv_trace:       self.conv_trace.clone(),
            temp_step_losses: self.temp_step_losses.clone(),
            cell_assignments: Vec::new(),  // large — skip in fallback
            num_clusters:     self.num_clusters,
            global_best:      self.global_best,
        }
    }
}

impl Default for PlotData {
    fn default() -> Self {
        PlotData {
            restart_losses:   Vec::new(),
            conv_trace:       Vec::new(),
            temp_step_losses: Vec::new(),
            cell_assignments: Vec::new(),
            num_clusters:     0,
            global_best:      f32::NEG_INFINITY,  // correct sentinel
        }
    }
}

// ── Main entry point ──────────────────────────────────────────────────────────

/// Write all enabled plots to plot_dir.  Returns the list of files written.
pub fn write_plots(data: &PlotData, plot_dir: &str, flags: &PlotFlags)
    -> Vec<String>
{
    fs::create_dir_all(plot_dir)
        .unwrap_or_else(|e| panic!("Cannot create plot directory '{}': {}", plot_dir, e));

    let mut written = Vec::new();

    if flags.restart_landscape {
        let path = format!("{}/restart_landscape.svg", plot_dir);
        fs::write(&path, plot_restart_landscape(data))
            .unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
        written.push(path);
    }
    if flags.convergence_curves && !data.conv_trace.is_empty() {
        let path = format!("{}/convergence_curves.svg", plot_dir);
        fs::write(&path, plot_convergence_curves(data))
            .unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
        written.push(path);
    }
    if flags.annealing_profile && !data.temp_step_losses.is_empty() {
        let path = format!("{}/annealing_profile.svg", plot_dir);
        fs::write(&path, plot_annealing_profile(data))
            .unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
        written.push(path);
    }
    if flags.cluster_balance && !data.cell_assignments.is_empty() {
        let path = format!("{}/cluster_balance.svg", plot_dir);
        fs::write(&path, plot_cluster_balance(data))
            .unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
        written.push(path);
    }
    if flags.posterior_confidence && !data.cell_assignments.is_empty() {
        let path = format!("{}/posterior_confidence.svg", plot_dir);
        fs::write(&path, plot_posterior_confidence(data))
            .unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
        written.push(path);
    }
    if flags.thread_efficiency && !data.restart_losses.is_empty() {
        let path = format!("{}/thread_efficiency.svg", plot_dir);
        fs::write(&path, plot_thread_efficiency(data))
            .unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
        written.push(path);
    }

    written
}

// ── PlotFlags — which plots to generate ──────────────────────────────────────

pub struct PlotFlags {
    pub restart_landscape:    bool,
    pub convergence_curves:   bool,
    pub annealing_profile:    bool,
    pub cluster_balance:      bool,
    pub posterior_confidence: bool,
    pub thread_efficiency:    bool,
}

impl PlotFlags {
    /// Enable all plots.
    pub fn all() -> Self {
        PlotFlags {
            restart_landscape:    true,
            convergence_curves:   true,
            annealing_profile:    true,
            cluster_balance:      true,
            posterior_confidence: true,
            thread_efficiency:    true,
        }
    }

    /// Parse from a comma-separated CLI string like "restart,balance,confidence"
    /// or the special values "all" / "none".
    pub fn from_str(s: &str) -> Self {
        if s == "all"  { return Self::all(); }
        if s == "none" { return Self::none(); }
        let parts: Vec<&str> = s.split(',').map(str::trim).collect();
        PlotFlags {
            restart_landscape:    parts.iter().any(|p| *p == "restart"    || *p == "restart_landscape"),
            convergence_curves:   parts.iter().any(|p| *p == "convergence"|| *p == "convergence_curves"),
            annealing_profile:    parts.iter().any(|p| *p == "annealing"  || *p == "annealing_profile"),
            cluster_balance:      parts.iter().any(|p| *p == "balance"    || *p == "cluster_balance"),
            posterior_confidence: parts.iter().any(|p| *p == "confidence" || *p == "posterior_confidence"),
            thread_efficiency:    parts.iter().any(|p| *p == "threads"    || *p == "thread_efficiency"),
        }
    }

    pub fn none() -> Self {
        PlotFlags {
            restart_landscape: false, convergence_curves: false,
            annealing_profile: false, cluster_balance: false,
            posterior_confidence: false, thread_efficiency: false,
        }
    }
}

// ── SVG helpers ───────────────────────────────────────────────────────────────

const W: f32 = 900.0;
const H: f32 = 520.0;
const PAD_L: f32 = 80.0;
const PAD_R: f32 = 40.0;
const PAD_T: f32 = 50.0;
const PAD_B: f32 = 70.0;
const PW: f32 = W - PAD_L - PAD_R;   // plot area width
const PH: f32 = H - PAD_T - PAD_B;   // plot area height

// PhD-grade colour palette (colour-blind safe, print-friendly)
const COLOURS: &[&str] = &[
    "#2166ac", "#d6604d", "#1b7837", "#762a83",
    "#b35806", "#4d4d4d", "#878787", "#01665e",
];

fn colour(i: usize) -> &'static str { COLOURS[i % COLOURS.len()] }

fn svg_header(title: &str, subtitle: &str) -> String {
    format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="'Helvetica Neue',Arial,sans-serif">
  <rect width="{W}" height="{H}" fill="{C_BG}" rx="6" />
  <!-- title -->
  <text x="{cx}" y="28" text-anchor="middle" font-size="15" font-weight="bold" fill="{C_TITLE}">{title}</text>
  <text x="{cx}" y="44" text-anchor="middle" font-size="11" fill="{C_TEXT}">{subtitle}</text>
  <!-- plot area border -->
  <rect x="{PAD_L}" y="{PAD_T}" width="{PW}" height="{PH}" fill="{C_WHITE}" stroke="{C_MGRAY}" stroke-width="1" rx="2"/>
"#,
        W = W, H = H, cx = W / 2.0,
        title = title, subtitle = subtitle,
        PAD_L = PAD_L, PAD_T = PAD_T, PW = PW, PH = PH
    )
}

fn svg_footer() -> &'static str { "</svg>\n" }

/// Map a data value to SVG x coordinate within the plot area.
fn px(val: f32, min: f32, max: f32) -> f32 {
    PAD_L + (val - min) / (max - min).max(1e-9) * PW
}

/// Map a data value to SVG y coordinate within the plot area (Y inverted).
fn py(val: f32, min: f32, max: f32) -> f32 {
    PAD_T + PH - (val - min) / (max - min).max(1e-9) * PH
}

fn x_axis_label(label: &str) -> String {
    format!(
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" fill="{C_DGRAY}">{}</text>"#,
        PAD_L + PW / 2.0,
        H - 10.0,
        label
    )
}

fn y_axis_label(label: &str) -> String {
    format!(
        r#"<text transform="rotate(-90)" x="{}" y="18" text-anchor="middle" font-size="12" fill="{C_DGRAY}">{}</text>"#,
        -(PAD_T + PH / 2.0),
        label
    )
}

fn tick_y(y_svg: f32, label: &str) -> String {
    let xt  = PAD_L - 4.0;
    let yt  = y_svg + 3.5;
    format!(
        r#"<line x1="{PAD_L}" y1="{y_svg}" x2="{xt}" y2="{y_svg}" stroke="{C_GRAY}" stroke-width="0.8"/>
           <text x="{xt}" y="{yt}" text-anchor="end" font-size="9" fill="{C_DIM}">{label}</text>"#,
        y_svg = y_svg, xt = xt, yt = yt, label = label,
    )
}

fn tick_x(x_svg: f32, label: &str) -> String {
    format!(
        r#"<line x1="{x}" y1="{}" x2="{x}" y2="{}" stroke="{C_GRAY}" stroke-width="0.8"/>
           <text x="{x}" y="{}" text-anchor="middle" font-size="9" fill="{C_DIM}">{}</text>"#,
        PAD_T, PAD_T + PH + 4.0, PAD_T + PH + 14.0, label,
        x = x_svg
    )
}

fn nice_ticks(min: f32, max: f32, count: usize) -> Vec<f32> {
    let range = max - min;
    if range <= 0.0 { return vec![min]; }
    let raw_step = range / count as f32;
    let mag = raw_step.log10().floor();
    let step = (raw_step / 10f32.powf(mag)).ceil() * 10f32.powf(mag);
    let start = (min / step).ceil() * step;
    let mut ticks = Vec::new();
    let mut t = start;
    while t <= max + step * 0.01 {
        ticks.push(t);
        t += step;
    }
    ticks
}

fn fmt_loss(v: f32) -> String {
    if v.abs() >= 1_000_000.0 { format!("{:.1}M", v / 1_000_000.0) }
    else if v.abs() >= 1_000.0 { format!("{:.0}k", v / 1_000.0) }
    else { format!("{:.1}", v) }
}

// ── Plot 1: Restart landscape ─────────────────────────────────────────────────

fn plot_restart_landscape(data: &PlotData) -> String {
    let mut losses_sorted: Vec<(usize, f32)> = data.restart_losses
        .iter()
        .enumerate()
        .map(|(_i, &(thread, _, loss))| (thread, loss))
        .collect();
    // Sort by loss value ascending (worst → best)
    losses_sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let min_loss = losses_sorted.first().map(|x| x.1).unwrap_or(-1.0);
    let max_loss = losses_sorted.last().map(|x| x.1).unwrap_or(0.0);
    let n = losses_sorted.len();

    let mut svg = svg_header(
        "Restart Landscape",
        &format!("{} restarts — sorted by final log-likelihood. Horizontal line = global optimum.", n),
    );

    // Y grid ticks
    let yticks = nice_ticks(min_loss, max_loss, 6);
    for &t in &yticks {
        svg += &tick_y(py(t, min_loss, max_loss), &fmt_loss(t));
    }
    // X grid ticks (every 10 restarts)
    for i in (0..=n).step_by(10) {
        let x = px(i as f32, 0.0, n as f32);
        svg += &tick_x(x, &i.to_string());
    }

    // Global best line
    let yb = py(data.global_best, min_loss, max_loss);
    svg += &format!(
        r#"<line x1="{}" y1="{yb}" x2="{}" y2="{yb}" stroke="{C_RED}" stroke-width="1.5" stroke-dasharray="6,3"/>
           <text x="{}" y="{}" font-size="9" fill="{C_RED}">global best {}</text>"#,
        PAD_L, PAD_L + PW,
        PAD_L + PW - 2.0, yb - 4.0,
        fmt_loss(data.global_best)
    );

    // Bars (coloured by thread)
    let bar_w = (PW / n as f32 * 0.7).max(1.5);
    for (i, &(thread, loss)) in losses_sorted.iter().enumerate() {
        let x  = px(i as f32 + 0.5, 0.0, n as f32) - bar_w / 2.0;
        let y1 = py(loss, min_loss, max_loss);
        let y0 = py(min_loss, min_loss, max_loss);
        let h  = (y0 - y1).abs().max(1.0);
        let col = colour(thread);
        // highlight global-best restarts
        let stroke = if (loss - data.global_best).abs() < 0.1 { "{C_RED}" } else { col };
        svg += &format!(
            r#"<rect x="{x:.1}" y="{y1:.1}" width="{bar_w:.1}" height="{h:.1}" fill="{col}" fill-opacity="0.75" stroke="{stroke}" stroke-width="0.4"/>"#,
        );
    }

    // Legend: one swatch per thread
    let threads: Vec<usize> = {
        let mut v: Vec<usize> = data.restart_losses.iter().map(|&(t,_,_)| t).collect();
        v.sort_unstable(); v.dedup(); v
    };
    let lx = PAD_L + PW - (threads.len() as f32) * 65.0;
    for (j, &t) in threads.iter().enumerate() {
        let x = lx + j as f32 * 65.0;
        svg += &format!(
            r#"<rect x="{x}" y="{}" width="10" height="10" fill="{}"/>
               <text x="{}" y="{}" font-size="9" fill="{C_DGRAY}">thread {t}</text>"#,
            PAD_T + 8.0, colour(t),
            x + 13.0, PAD_T + 17.0
        );
    }

    svg += &x_axis_label("Restart index (sorted by loss)");
    svg += &y_axis_label("Final log-likelihood");
    svg += svg_footer();
    svg
}

// ── Plot 2: Convergence curves ────────────────────────────────────────────────

fn plot_convergence_curves(data: &PlotData) -> String {
    // Group CONV trace by (thread, epoch) → list of (iter, loss)
    use std::collections::BTreeMap;
    let mut by_restart: BTreeMap<(usize,usize), Vec<(usize,f32)>> = BTreeMap::new();
    for &(thread, epoch, _ts, iter, loss) in &data.conv_trace {
        if loss.is_infinite() || loss.is_nan() { continue; }
        by_restart.entry((thread, epoch)).or_default().push((iter, loss));
    }
    // Each restart: sort by iter, compute cumulative iter offset across temp steps
    let mut restart_traces: Vec<Vec<(f32,f32)>> = Vec::new();
    let mut thread_for_restart: Vec<usize> = Vec::new();
    for (&(thread, _epoch), pts) in &by_restart {
        let mut sorted = pts.clone();
        sorted.sort_by_key(|&(i,_)| i);
        restart_traces.push(sorted.iter().enumerate().map(|(i,&(_,l))| (i as f32, l)).collect());
        thread_for_restart.push(thread);
    }

    if restart_traces.is_empty() {
        return String::from("<!-- no CONV data: run with --verbose to enable convergence curves -->");
    }

    let max_iter = restart_traces.iter().map(|t| t.len()).max().unwrap_or(1) as f32;
    let all_losses: Vec<f32> = restart_traces.iter().flatten().map(|&(_,l)| l).collect();
    // exclude extreme outliers (first iteration always has huge loss)
    let mut sorted_l = all_losses.clone(); sorted_l.sort_by(|a,b| a.partial_cmp(b).unwrap());
    let p5  = sorted_l[(sorted_l.len() as f32 * 0.02) as usize];
    let p95 = sorted_l[((sorted_l.len() as f32 * 0.98) as usize).min(sorted_l.len()-1)];

    let mut svg = svg_header(
        "EM Convergence Curves",
        &format!("{} restarts × all iterations. Y-axis clipped to 2nd–98th percentile of loss.", restart_traces.len()),
    );

    let yticks = nice_ticks(p5, p95, 6);
    for &t in &yticks {
        svg += &tick_y(py(t, p5, p95), &fmt_loss(t));
    }
    let xticks = nice_ticks(0.0, max_iter, 8);
    for &t in &xticks {
        svg += &tick_x(px(t, 0.0, max_iter), &format!("{:.0}", t));
    }

    // Draw curves (thin, semi-transparent — many will overlap)
    for (trace, &thread) in restart_traces.iter().zip(thread_for_restart.iter()) {
        if trace.len() < 2 { continue; }
        let pts: String = trace.iter()
            .map(|&(i,l)| {
                let x = px(i, 0.0, max_iter);
                let y = py(l.clamp(p5, p95), p5, p95);
                format!("{:.1},{:.1}", x, y)
            })
            .collect::<Vec<_>>()
            .join(" ");
        svg += &format!(
            r#"<polyline points="{pts}" fill="none" stroke="{}" stroke-width="0.7" stroke-opacity="0.35"/>"#,
            colour(thread)
        );
    }

    // Overlay the best-restart curve thicker
    if let Some((_best_idx, _)) = data.restart_losses.iter()
        .enumerate()
        .max_by(|a, b| a.1.2.partial_cmp(&b.1.2).unwrap())
    {
        // find the trace for the best restart — identify by global_best loss value
        for (trace, &_thread) in restart_traces.iter().zip(thread_for_restart.iter()) {
            if trace.last().map(|&(_,l)| (l - data.global_best).abs() < 1.0).unwrap_or(false) {
                let pts: String = trace.iter()
                    .map(|&(i,l)| {
                        let x = px(i, 0.0, max_iter);
                        let y = py(l.clamp(p5, p95), p5, p95);
                        format!("{:.1},{:.1}", x, y)
                    })
                    .collect::<Vec<_>>()
                    .join(" ");
                svg += &format!(
                    r#"<polyline points="{pts}" fill="none" stroke="{C_RED}" stroke-width="2.0"/>"#
                );
                break;
            }
        }
    }

    // Annotation: best run label
    svg += &format!(
        r#"<text x="{}" y="{}" font-size="10" fill="{C_RED}">— best restart</text>"#,
        PAD_L + PW - 100.0, PAD_T + 16.0
    );

    svg += &x_axis_label("Iteration (cumulative across temperature steps)");
    svg += &y_axis_label("Log-likelihood");
    svg += svg_footer();
    svg
}

// ── Plot 3: Annealing profile ─────────────────────────────────────────────────

fn plot_annealing_profile(data: &PlotData) -> String {
    // Aggregate: for each temp_step (0..8), collect all final losses
    let mut by_step: Vec<Vec<f32>> = vec![Vec::new(); 9];
    for &(_thread, _epoch, step, _iters, loss) in &data.temp_step_losses {
        if step < 9 && !loss.is_nan() && !loss.is_infinite() {
            by_step[step].push(loss);
        }
    }

    // Compute mean ± sd per step
    let stats: Vec<(f32, f32, f32, f32)> = by_step.iter().map(|vals| {
        if vals.is_empty() { return (0.0, 0.0, 0.0, 0.0); }
        let n = vals.len() as f32;
        let mean = vals.iter().sum::<f32>() / n;
        let sd = (vals.iter().map(|v| (v - mean).powi(2)).sum::<f32>() / n).sqrt();
        let min = vals.iter().cloned().fold(f32::INFINITY, f32::min);
        let max = vals.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        (mean, sd, min, max)
    }).collect();

    let all_vals: Vec<f32> = data.temp_step_losses.iter()
        .filter(|&&(_,_,_,_,l)| l.is_finite())
        .map(|&(_,_,_,_,l)| l).collect();
    let min_v = all_vals.iter().cloned().fold(f32::INFINITY, f32::min);
    let max_v = all_vals.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

    let mut svg = svg_header(
        "Annealing Profile — Loss by Temperature Step",
        "Mean ± 1 SD of log-likelihood at each of the 9 annealing steps. Lower is less converged.",
    );

    let yticks = nice_ticks(min_v, max_v, 6);
    for &t in &yticks {
        svg += &tick_y(py(t, min_v, max_v), &fmt_loss(t));
    }
    for step in 0..9usize {
        let x = px(step as f32, 0.0, 8.0);
        svg += &tick_x(x, &format!("T{}", step));
    }

    // Temperature step background bands (alternating)
    for step in 0..9usize {
        if step % 2 == 0 {
            let x1 = px(step as f32 - 0.5, 0.0, 8.0).max(PAD_L);
            let x2 = px(step as f32 + 0.5, 0.0, 8.0).min(PAD_L + PW);
            svg += &format!(
                r#"<rect x="{x1:.1}" y="{PAD_T}" width="{:.1}" height="{PH}" fill="{C_BGBAND}" opacity="0.5"/>"#,
                x2 - x1
            );
        }
    }

    // Min-max range band
    let band_pts: String = (0..9usize).map(|s| {
        let x = px(s as f32, 0.0, 8.0);
        format!("{:.1},{:.1}", x, py(stats[s].3, min_v, max_v))
    }).chain((0..9usize).rev().map(|s| {
        let x = px(s as f32, 0.0, 8.0);
        format!("{:.1},{:.1}", x, py(stats[s].2, min_v, max_v))
    })).collect::<Vec<_>>().join(" ");
    svg += &format!(r#"<polygon points="{band_pts}" fill="{C_BLUE}" fill-opacity="0.08"/>"#);

    // SD band
    let sd_pts: String = (0..9usize).map(|s| {
        let x = px(s as f32, 0.0, 8.0);
        let y = py((stats[s].0 + stats[s].1).min(max_v), min_v, max_v);
        format!("{:.1},{:.1}", x, y)
    }).chain((0..9usize).rev().map(|s| {
        let x = px(s as f32, 0.0, 8.0);
        let y = py((stats[s].0 - stats[s].1).max(min_v), min_v, max_v);
        format!("{:.1},{:.1}", x, y)
    })).collect::<Vec<_>>().join(" ");
    svg += &format!(r#"<polygon points="{sd_pts}" fill="{C_BLUE}" fill-opacity="0.22"/>"#);

    // Mean line
    let mean_pts: String = (0..9usize).map(|s| {
        let x = px(s as f32, 0.0, 8.0);
        let y = py(stats[s].0, min_v, max_v);
        format!("{:.1},{:.1}", x, y)
    }).collect::<Vec<_>>().join(" ");
    svg += &format!(r#"<polyline points="{mean_pts}" fill="none" stroke="{C_BLUE}" stroke-width="2.5"/>"#);

    // Mean dots + annotations
    for step in 0..9usize {
        let x = px(step as f32, 0.0, 8.0);
        let y = py(stats[step].0, min_v, max_v);
        svg += &format!(
            r#"<circle cx="{x:.1}" cy="{y:.1}" r="4" fill="{C_BLUE}"/>
               <text x="{x:.1}" y="{:.1}" text-anchor="middle" font-size="8" fill="{C_BLUE}">{}</text>"#,
            y - 7.0, fmt_loss(stats[step].0)
        );
    }

    // Legend
    svg += &format!(
        r#"<rect x="{}" y="{}" width="10" height="10" fill="{C_BLUE}" fill-opacity="0.22"/>
           <text x="{}" y="{}" font-size="9" fill="{C_DGRAY}">± 1 SD</text>
           <rect x="{}" y="{}" width="10" height="10" fill="{C_BLUE}" fill-opacity="0.08"/>
           <text x="{}" y="{}" font-size="9" fill="{C_DGRAY}">min–max range</text>"#,
        PAD_L + 8.0, PAD_T + 8.0,
        PAD_L + 22.0, PAD_T + 17.0,
        PAD_L + 90.0, PAD_T + 8.0,
        PAD_L + 104.0, PAD_T + 17.0
    );

    svg += &x_axis_label("Temperature step (T0 = highest temp, T8 = final)");
    svg += &y_axis_label("Log-likelihood");
    svg += svg_footer();
    svg
}

// ── Plot 4: Cluster balance ───────────────────────────────────────────────────

fn plot_cluster_balance(data: &PlotData) -> String {
    let mut counts = vec![0usize; data.num_clusters];
    for &(cluster, _) in &data.cell_assignments {
        if cluster < data.num_clusters { counts[cluster] += 1; }
    }
    let total: usize = counts.iter().sum();
    let expected = total as f32 / data.num_clusters as f32;
    let max_count = *counts.iter().max().unwrap_or(&1) as f32;

    let mut svg = svg_header(
        "Cluster Balance",
        &format!("{} cells across {} clusters. Dashed line = expected if equal split.", total, data.num_clusters),
    );

    let bar_h = PH / data.num_clusters as f32 * 0.65;
    let gap   = PH / data.num_clusters as f32;

    // Expected line
    let xe = px(expected, 0.0, max_count * 1.05);
    svg += &format!(
        r#"<line x1="{xe:.1}" y1="{PAD_T}" x2="{xe:.1}" y2="{}" stroke="{C_MGRAY}" stroke-width="1.2" stroke-dasharray="5,3"/>
           <text x="{xe:.1}" y="{}" text-anchor="middle" font-size="9" fill="{C_MGRAY}">expected ({:.0})</text>"#,
        PAD_T + PH,
        PAD_T - 4.0,
        expected
    );

    // X ticks
    let xticks = nice_ticks(0.0, max_count * 1.05, 6);
    for &t in &xticks {
        svg += &tick_x(px(t, 0.0, max_count * 1.05), &format!("{:.0}", t));
    }

    for (i, &count) in counts.iter().enumerate() {
        let y    = PAD_T + i as f32 * gap + (gap - bar_h) / 2.0;
        let w    = px(count as f32, 0.0, max_count * 1.05) - PAD_L;
        let pct  = 100.0 * count as f32 / total as f32;
        let col  = colour(i);

        // Bar
        svg += &format!(
            r#"<rect x="{PAD_L}" y="{y:.1}" width="{w:.1}" height="{bar_h:.1}" fill="{col}" rx="2"/>"#,
        );

        // Label inside or outside
        let label_x = if w > 80.0 { PAD_L + w - 5.0 } else { PAD_L + w + 5.0 };
        let anchor   = if w > 80.0 { "end" } else { "start" };
        svg += &format!(
            r#"<text x="{label_x:.1}" y="{:.1}" text-anchor="{anchor}" font-size="11" font-weight="bold" fill="{}">
                   {count} ({pct:.1}%)
               </text>"#,
            y + bar_h * 0.65,
            if w > 80.0 { "white" } else { col }
        );

        // Cluster label on left
        svg += &format!(
            r#"<text x="{}" y="{:.1}" text-anchor="end" font-size="11" fill="{C_DGRAY}">Cluster {i}</text>"#,
            PAD_L - 6.0, y + bar_h * 0.65
        );
    }

    svg += &x_axis_label("Cell count");
    svg += svg_footer();
    svg
}

// ── Plot 5: Posterior confidence histogram ────────────────────────────────────

fn plot_posterior_confidence(data: &PlotData) -> String {
    // Confidence = best_log_prob − second_best_log_prob (always ≥ 0)
    let confidences: Vec<f32> = data.cell_assignments.iter().map(|(_, lps)| {
        if lps.len() < 2 { return 0.0; }
        let mut sorted = lps.clone();
        sorted.sort_by(|a, b| b.partial_cmp(a).unwrap()); // descending
        sorted[0] - sorted[1]                              // best - second
    }).collect();

    let _max_conf = confidences.iter().cloned().fold(0.0f32, f32::max);
    let n_cells  = confidences.len();

    // Histogram with 60 bins, cap display at 95th percentile for readability
    let n_bins = 60usize;
    let mut sorted_c = confidences.clone();
    sorted_c.sort_by(|a,b| a.partial_cmp(b).unwrap());
    let p95 = sorted_c[((sorted_c.len() as f32 * 0.95) as usize).min(sorted_c.len()-1)];
    let p50 = sorted_c[sorted_c.len() / 2];
    let bin_w = p95 / n_bins as f32;
    let mut bins = vec![0usize; n_bins + 1];
    for &c in &confidences {
        let bi = ((c / bin_w) as usize).min(n_bins);
        bins[bi] += 1;
    }
    let max_bin = *bins.iter().max().unwrap_or(&1) as f32;

    // Ambiguous cell threshold: confidence < 10 log units
    let ambig_threshold = 10.0f32;
    let n_ambig = confidences.iter().filter(|&&c| c < ambig_threshold).count();
    let pct_ambig = 100.0 * n_ambig as f32 / n_cells as f32;

    let mut svg = svg_header(
        "Posterior Confidence Distribution",
        &format!(
            "{n_cells} cells. Confidence = best − 2nd-best log-prob. Low confidence → potential doublets.",
        ),
    );

    let yticks = nice_ticks(0.0, max_bin, 5);
    for &t in &yticks {
        svg += &tick_y(py(t, 0.0, max_bin), &format!("{:.0}", t));
    }
    let xticks = nice_ticks(0.0, p95, 7);
    for &t in &xticks {
        svg += &tick_x(px(t, 0.0, p95), &format!("{:.0}", t));
    }

    // Ambiguous zone shading
    let x_amb = px(ambig_threshold, 0.0, p95);
    svg += &format!(
        r#"<rect x="{PAD_L}" y="{PAD_T}" width="{:.1}" height="{PH}" fill="{C_AMBER}" opacity="0.35"/>
           <text x="{}" y="{}" font-size="9" fill="{C_RED}">ambiguous &lt;{:.0} ({:.1}% of cells)</text>"#,
        (x_amb - PAD_L).max(0.0),
        PAD_L + 3.0, PAD_T + 13.0,
        ambig_threshold, pct_ambig
    );

    // Histogram bars
    for (bi, &count) in bins.iter().enumerate() {
        if count == 0 { continue; }
        let x  = px(bi as f32 * bin_w, 0.0, p95);
        let x2 = px((bi + 1) as f32 * bin_w, 0.0, p95);
        let y  = py(count as f32, 0.0, max_bin);
        let y0 = PAD_T + PH;
        let col = if (bi as f32 * bin_w) < ambig_threshold { "{C_RED}" } else { "{C_BLUE}" };
        svg += &format!(
            r#"<rect x="{x:.1}" y="{y:.1}" width="{:.1}" height="{:.1}" fill="{col}" fill-opacity="0.7"/>"#,
            (x2 - x).max(1.0), y0 - y
        );
    }

    // Median line
    let xm = px(p50, 0.0, p95);
    svg += &format!(
        r#"<line x1="{xm:.1}" y1="{PAD_T}" x2="{xm:.1}" y2="{}" stroke="{C_GREEN}" stroke-width="1.5" stroke-dasharray="4,3"/>
           <text x="{:.1}" y="{}" text-anchor="middle" font-size="9" fill="{C_GREEN}">median {:.0}</text>"#,
        PAD_T + PH,
        xm, PAD_T - 4.0, p50
    );

    // Stats box
    svg += &format!(
        r#"<rect x="{}" y="{}" width="190" height="50" fill="{C_WHITE}" stroke="{C_MGRAY}" rx="3"/>
           <text x="{}" y="{}" font-size="9" fill="{C_DGRAY}">n cells = {n_cells}</text>
           <text x="{}" y="{}" font-size="9" fill="{C_DGRAY}">median confidence = {p50:.1} log-units</text>
           <text x="{}" y="{}" font-size="9" fill="{C_RED}">ambiguous (&lt;{:.0}) = {n_ambig} ({pct_ambig:.1}%)</text>"#,
        PAD_L + PW - 200.0, PAD_T + 10.0,
        PAD_L + PW - 196.0, PAD_T + 24.0,
        PAD_L + PW - 196.0, PAD_T + 36.0,
        PAD_L + PW - 196.0, PAD_T + 48.0,
        ambig_threshold
    );

    svg += &x_axis_label("Posterior confidence (log-probability units, clipped at 95th percentile)");
    svg += &y_axis_label("Cell count");
    svg += svg_footer();
    svg
}

// ── Plot 6: Thread efficiency ─────────────────────────────────────────────────

fn plot_thread_efficiency(data: &PlotData) -> String {
    // Per thread: restart count and best loss
    let mut thread_restarts: std::collections::BTreeMap<usize, (usize, f32)> =
        std::collections::BTreeMap::new();
    for &(thread, _, loss) in &data.restart_losses {
        let e = thread_restarts.entry(thread).or_insert((0, f32::NEG_INFINITY));
        e.0 += 1;
        if loss > e.1 { e.1 = loss; }
    }

    let threads: Vec<usize> = thread_restarts.keys().cloned().collect();
    let n_threads = threads.len();

    let losses: Vec<f32> = threads.iter().map(|t| thread_restarts[t].1).collect();
    let counts: Vec<usize> = threads.iter().map(|t| thread_restarts[t].0).collect();
    let max_count = *counts.iter().max().unwrap_or(&1) as f32;
    let min_loss = losses.iter().cloned().fold(f32::INFINITY, f32::min);
    let max_loss = losses.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

    // Dual-axis plot: bar = restarts per thread, line = best loss per thread
    let mut svg = svg_header(
        "Thread Efficiency",
        "Bars = restarts per thread. Line = best log-likelihood found by each thread.",
    );

    let bar_w  = PW / n_threads as f32 * 0.6;
    let x_step = PW / n_threads as f32;

    // X ticks (thread labels)
    for (j, &t) in threads.iter().enumerate() {
        let x = PAD_L + (j as f32 + 0.5) * x_step;
        svg += &format!(
            r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="10" fill="{C_DGRAY}">T{t}</text>"#,
            PAD_T + PH + 14.0
        );
    }

    // Y left (restart count)
    let yticks_l = nice_ticks(0.0, max_count, 5);
    for &t in &yticks_l {
        let y = py(t, 0.0, max_count);
        svg += &format!(
            r#"<line x1="{PAD_L}" y1="{y:.1}" x2="{}" y2="{y:.1}" stroke="{C_LGRAY}" stroke-width="0.8"/>
               <text x="{}" y="{:.1}" text-anchor="end" font-size="9" fill="{C_BLUE}">{:.0}</text>"#,
            PAD_L + PW, PAD_L - 4.0, y + 3.0, t
        );
    }

    // Bars (restart count, blue)
    for (j, (&t, &count)) in threads.iter().zip(counts.iter()).enumerate() {
        let cx = PAD_L + (j as f32 + 0.5) * x_step;
        let x  = cx - bar_w / 2.0;
        let y  = py(count as f32, 0.0, max_count);
        let h  = PAD_T + PH - y;
        svg += &format!(
            r#"<rect x="{x:.1}" y="{y:.1}" width="{bar_w:.1}" height="{h:.1}" fill="{}" fill-opacity="0.55" rx="2"/>
               <text x="{cx:.1}" y="{:.1}" text-anchor="middle" font-size="10" fill="{C_BLUE}" font-weight="bold">{count}</text>"#,
            colour(t), y - 4.0
        );
    }

    // Y right axis (loss) — map to same pixel space but different scale
    let yr_loss_range = (max_loss - min_loss).max(1.0);
    let loss_to_y = |l: f32| -> f32 {
        PAD_T + PH - (l - min_loss) / yr_loss_range * PH
    };

    // Right axis ticks
    let yticks_r = nice_ticks(min_loss, max_loss, 5);
    let rx = PAD_L + PW;
    for &t in &yticks_r {
        let y = loss_to_y(t);
        svg += &format!(
            r#"<text x="{}" y="{:.1}" text-anchor="start" font-size="9" fill="{C_RED}">{}</text>"#,
            rx + 3.0, y + 3.0, fmt_loss(t)
        );
    }

    // Loss line (red, right axis)
    let loss_pts: String = threads.iter().zip(losses.iter()).enumerate().map(|(j, (_, &l))| {
        let x = PAD_L + (j as f32 + 0.5) * x_step;
        let y = loss_to_y(l);
        format!("{:.1},{:.1}", x, y)
    }).collect::<Vec<_>>().join(" ");
    svg += &format!(
        r#"<polyline points="{loss_pts}" fill="none" stroke="{C_RED}" stroke-width="2"/>
           <text x="{}" y="{}" font-size="9" fill="{C_RED}">→ best loss (right axis)</text>"#,
        PAD_L + PW - 170.0, PAD_T + 16.0
    );
    // Dots on loss line
    for (j, &l) in losses.iter().enumerate() {
        let x = PAD_L + (j as f32 + 0.5) * x_step;
        let y = loss_to_y(l);
        svg += &format!(r#"<circle cx="{x:.1}" cy="{y:.1}" r="4" fill="{C_RED}"/>"#);
    }

    // Global best marker
    let yg = loss_to_y(data.global_best);
    svg += &format!(
        r#"<line x1="{PAD_L}" y1="{yg:.1}" x2="{rx:.1}" y2="{yg:.1}" stroke="{C_RED}" stroke-width="1" stroke-dasharray="4,3"/>
           <text x="{}" y="{:.1}" font-size="9" fill="{C_RED}">global best</text>"#,
        rx + 3.0, yg - 2.0
    );

    // Legend
    svg += &format!(
        r#"<rect x="{}" y="{}" width="10" height="10" fill="{C_BLUE}" fill-opacity="0.55"/>
           <text x="{}" y="{}" font-size="9" fill="{C_DGRAY}">restarts / thread</text>"#,
        PAD_L + 8.0, PAD_T + 8.0,
        PAD_L + 22.0, PAD_T + 17.0
    );

    svg += &x_axis_label("Thread");
    svg += &y_axis_label("Restarts per thread");
    svg += svg_footer();
    svg
}
