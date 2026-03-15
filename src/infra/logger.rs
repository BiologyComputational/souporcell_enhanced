// ============================================================================
// logger.rs — Structured stderr logging, timing, progress, convergence summary
// ============================================================================
// All terminal output from souporcell goes through this module so formatting
// is consistent, filterable, and machine-parseable if piped to a log file.
//
// Design principles:
//   - Every timed section uses the Logger struct to track wall-clock elapsed
//   - Progress lines are prefixed with a fixed tag so `grep` can extract them
//   - Convergence lines carry structured fields (thread, epoch, iter, loss)
//   - Summary at end is human-readable but also tab-delimited for scripting
// ============================================================================

use std::time::{Duration, Instant};

// ── ANSI colour helpers (only if stderr is a tty) ────────────────────────────

fn is_tty() -> bool {
    // simple heuristic: check if TERM is set
    std::env::var("TERM").map(|t| t != "dumb").unwrap_or(false)
}

fn cyan()    -> &'static str { if is_tty() { "\x1b[36m"    } else { "" } }
fn yellow()  -> &'static str { if is_tty() { "\x1b[33m"    } else { "" } }
fn green()   -> &'static str { if is_tty() { "\x1b[32m"    } else { "" } }
fn red()     -> &'static str { if is_tty() { "\x1b[31m"    } else { "" } }
fn bold()    -> &'static str { if is_tty() { "\x1b[1m"     } else { "" } }
fn dim()     -> &'static str { if is_tty() { "\x1b[2m"     } else { "" } }
fn reset()   -> &'static str { if is_tty() { "\x1b[0m"     } else { "" } }

// ── Top-level pipeline logger ────────────────────────────────────────────────

pub struct Logger {
    start: Instant,
}

impl Logger {
    pub fn new() -> Self {
        Logger { start: Instant::now() }
    }

    /// Print a section banner  ─────────────────────────────────────────────
    pub fn section(&self, title: &str) {
        let elapsed = fmt_duration(self.start.elapsed());
        eprintln!(
            "\n{}{}══ {} {}[+{}]{}",
            bold(), cyan(), title, dim(), elapsed, reset()
        );
        eprintln!("{}{}{}",
            cyan(),
            "─".repeat(78),
            reset()
        );
    }

    /// Print a key-value info line  ─────────────────────────────────────────
    pub fn info(&self, key: &str, val: &str) {
        eprintln!("{}  │{} {:.<40} {}{}", dim(), reset(), key, yellow(), val);
        // the reset at end is intentional — no trailing escape on plain pipes
        eprint!("{}", reset());
    }

    /// Print a success / completion line with elapsed time ──────────────────
    pub fn ok(&self, msg: &str, since: Option<Instant>) {
        let timing = match since {
            Some(t) => format!("  [+{}]", fmt_duration(t.elapsed())),
            None    => String::new(),
        };
        eprintln!("{}  ✔  {}{}{}{}",
            green(), msg, dim(), timing, reset());
    }

    /// Print a warning  ──────────────────────────────────────────────────────
    pub fn warn(&self, msg: &str) {
        eprintln!("{}  ⚠  WARNING: {}{}", yellow(), msg, reset());
    }

    /// Print a fatal error and exit  ─────────────────────────────────────────
    #[allow(dead_code)]
    pub fn fatal(&self, msg: &str) -> ! {
        eprintln!("{}{}  ✖  FATAL: {}{}", bold(), red(), msg, reset());
        std::process::exit(1);
    }

    /// Total elapsed since Logger::new()  ────────────────────────────────────
    pub fn total_elapsed(&self) -> String {
        fmt_duration(self.start.elapsed())
    }
}

// ── EM / KHM convergence line ─────────────────────────────────────────────────
//
// Format (tab-separated, machine-parseable):
//   CONV  thread  epoch  iter  temp_step  log_loss  delta  clusters_populated
//
// The CONV prefix allows:  grep '^CONV' clustering.err | cut -f5  → loss series
pub fn log_convergence(
    thread_num:  usize,
    epoch:       usize,
    iteration:   usize,
    temp_step:   usize,
    log_loss:    f32,
    delta:       f32,
    populated:   usize,
    method:      &str,
) {
    eprintln!(
        "CONV\t{method}\t{thread}\t{epoch}\t{iter}\t{ts}\t{loss:.4}\t{delta:.4}\t{pop}",
        method = method,
        thread = thread_num,
        epoch  = epoch,
        iter   = iteration,
        ts     = temp_step,
        loss   = log_loss,
        delta  = delta,
        pop    = populated,
    );
}

// ── Thread restart progress ───────────────────────────────────────────────────
//
// Format: RESTART  thread  iteration  current_loss  best_so_far
pub fn log_restart(thread_num: usize, iteration: usize, loss: f32, best: f32) {
    eprintln!(
        "RESTART\t{thread}\t{iter}\t{loss:.4}\t{best:.4}",
        thread = thread_num,
        iter   = iteration,
        loss   = loss,
        best   = best,
    );
}

// ── Data loading progress ─────────────────────────────────────────────────────
pub fn log_loading_stats(
    total_loci:    usize,
    used_loci:     usize,
    total_cells:   usize,
    zero_coverage: usize,
) {
    let zero_pct = 100.0 * zero_coverage as f64 / (used_loci * total_cells).max(1) as f64;
    eprintln!(
        "{}  │  Loci in matrix:    {}{}",
        dim(), total_loci, reset()
    );
    eprintln!(
        "{}  │  Loci passing QC:   {}{}",
        dim(), used_loci, reset()
    );
    eprintln!(
        "{}  │  Cells loaded:      {}{}",
        dim(), total_cells, reset()
    );
    eprintln!(
        "{}  │  Zero-cov entries:  {} ({:.1}% of locus×cell pairs — expected >95% for sparse scRNA){}",
        dim(), zero_coverage, zero_pct, reset()
    );
    // Warn if unexpectedly low sparsity (suggests data issue)
    if zero_pct < 80.0 {
        eprintln!(
            "{}  ⚠  Sparsity {:.1}% is lower than expected for scRNA-seq — check input matrices{}",
            yellow(), zero_pct, reset()
        );
    }
}

// ── Bad cluster detection report ─────────────────────────────────────────────
pub fn log_cluster_analysis(
    run:              usize,
    cluster_idx:      usize,
    cell_count:       usize,
    is_outlier:       bool,
) {
    let tag = if is_outlier {
        format!("{}OUTLIER — will be reinitialized{}", red(), reset())
    } else {
        format!("{}ok{}", green(), reset())
    };
    eprintln!(
        "  CLUSTER_ANALYSIS\trun={}\tcluster={}\tcells={}\tstatus={}",
        run, cluster_idx, cell_count, tag
    );
}

// ── Final run summary ─────────────────────────────────────────────────────────
pub fn log_run_summary(
    best_log_probability: f32,
    num_clusters:         usize,
    cells_per_cluster:    &[usize],
    total_cells:          usize,
    restarts_run:         usize,
    total_elapsed:        &str,
) {
    eprintln!("\n{}{}══ CLUSTERING COMPLETE {}══{}",
        bold(), green(), "─".repeat(40), reset());
    eprintln!("{}  best total log probability : {:.4}{}",
        yellow(), best_log_probability, reset());
    eprintln!("{}  restarts completed         : {}{}",
        dim(), restarts_run, reset());
    eprintln!("{}  total runtime              : {}{}",
        dim(), total_elapsed, reset());
    eprintln!("{}  cluster breakdown:{}",
        dim(), reset());
    for (i, &count) in cells_per_cluster.iter().enumerate() {
        let pct = 100.0 * count as f64 / total_cells.max(1) as f64;
        let bar = "#".repeat((pct / 2.0) as usize);  // 50-char max bar
        eprintln!(
            "    cluster {:>2}  {:>5} cells  ({:>5.1}%)  {}{}{}",
            i, count, pct,
            cyan(), bar, reset()
        );
    }

    // Warn if any cluster is severely underpopulated (< 1% of cells)
    for (i, &count) in cells_per_cluster.iter().enumerate() {
        let pct = 100.0 * count as f64 / total_cells.max(1) as f64;
        if pct < 1.0 {
            eprintln!(
                "{}  ⚠  Cluster {} has only {} cells ({:.1}%) — consider reducing k{}",
                yellow(), i, count, pct, reset()
            );
        }
    }

    // Warn if num_clusters > 4 (doublet search cap in troublet)
    if num_clusters > 4 {
        eprintln!(
            "{}  ⚠  k={} > 4: troublet doublet search is capped at top-4 clusters by posterior.\n\
             {}     Some inter-donor doublets involving lower-ranked clusters may be missed.{}",
            yellow(), num_clusters,
            yellow(), reset()
        );
    }
}

// ── Utilities ─────────────────────────────────────────────────────────────────

/// Format a Duration as "Xh Ym Zs" or "Ym Zs" or "Zs" depending on magnitude
pub fn fmt_duration(d: Duration) -> String {
    let secs  = d.as_secs();
    let hours = secs / 3600;
    let mins  = (secs % 3600) / 60;
    let s     = secs % 60;
    if hours > 0 {
        format!("{}h {:02}m {:02}s", hours, mins, s)
    } else if mins > 0 {
        format!("{}m {:02}s", mins, s)
    } else {
        format!("{}.{:03}s", s, d.subsec_millis())
    }
}
