// ============================================================================
// main.rs — souporcell v2 entry point and pipeline orchestration
// ============================================================================
//
// Module map:
//   params.rs       — CLI parsing, Params struct, algorithm enums
//   types.rs        — CellData, ThreadData, ConvergenceStats
//   io.rs           — matrix loading, barcode loading, result writing
//   math.rs         — log_sum_exp, normalize, binomial_loss
//   em.rs           — EM algorithm + shared centre-update helpers
//   khm.rs          — KHM algorithm
//   cluster_init.rs — all cluster centre initialisation strategies
//   logger.rs       — structured stderr logging, timing, summaries
//
// main() responsibilities:
//   1. Parse params (params.rs)
//   2. Load data and log QC statistics (io.rs, logger.rs)
//   3. Run parallel restarts with EM or KHM (rayon)
//   4. Optionally run souporcell3 bad-cluster detection loop
//   5. Write results to stdout (io.rs)
//   6. Print final summary to stderr (logger.rs)
// ============================================================================

#[macro_use]
extern crate clap;

// ── Module tree ──────────────────────────────────────────────────────────────
//
//   domain/     Pure types and math (no I/O, no CLI)
//   config/     CLI parsing and Params struct
//   core/       Clustering algorithms (EM, KHM, init)
//   infra/      File I/O and terminal logging
//   analysis/   Post-clustering: diagnostic plots and HTML report
//
mod domain;
mod config;
mod core;
mod infra;
mod analysis;

use std::time::Instant;

use hashbrown::HashMap;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use std::sync::{Arc, Mutex};
use crate::core::cluster_init::init_cluster_centers;
use crate::analysis::plots::{PlotData, PlotFlags, write_plots};
use crate::analysis::report::{ReportData, write_report};
use crate::core::em::em;
use crate::infra::io::{load_barcodes, load_cell_data, write_cluster_assignments};
use crate::infra::preflight::{show_preflight, ApprovalResult};
use crate::infra::gui_server;
use crate::core::khm::khm;
use crate::infra::logger::{log_cluster_analysis, log_loading_stats, log_run_summary,
             log_restart, Logger};
use crate::config::params::{load_params, ClusterMethod, Params};
use crate::domain::types::{CellData, ThreadData};

// ── main ─────────────────────────────────────────────────────────────────────

fn main() {
    // ── GUI mode — check raw args before clap ─────────────────────────────────
    // Intercept --gui before load_params() so it works without matrix paths.
    let raw_args: Vec<String> = std::env::args().collect();
    if raw_args.iter().any(|a| a == "--gui") {
        let port: u16 = raw_args.iter()
            .position(|a| a == "--gui_port")
            .and_then(|i| raw_args.get(i + 1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(7979);
        gui_server::launch(port);
        return;
    }

    let log    = Logger::new();
    let params = load_params();
    let t_start = Instant::now();

    // ── Pre-flight approval gate (--dry_run) ──────────────────────────────────
    if params.dry_run {
        match show_preflight(&params) {
            ApprovalResult::Abort => std::process::exit(0),
            ApprovalResult::Proceed => {}
        }
    }

    // ── Data loading ──────────────────────────────────────────────────────────
    log.section("DATA LOADING");
    let t_load = Instant::now();
    log.info("ref matrix", &params.ref_mtx);
    log.info("alt matrix", &params.alt_mtx);
    log.info("barcodes",   &params.barcodes);
    log.info("k (clusters)", &params.num_clusters.to_string());
    log.info("restarts",   &params.restarts.to_string());
    log.info("threads",    &params.threads.to_string());
    log.info("seed",       &params.seed.to_string());
    log.info("method",     params.clustering_method.name());

    let barcodes = load_barcodes(&params);
    log.info("barcodes loaded", &barcodes.len().to_string());

    let loaded = load_cell_data(&params);
    log_loading_stats(
        loaded.total_loci_raw,
        loaded.loci_used,
        loaded.total_cells,
        loaded.zero_cov_entries,
    );

    // Warn if k > 4 (doublet search limitation in troublet)
    if params.num_clusters > 4 {
        log.warn(&format!(
            "k={} > 4: troublet doublet search is capped at the top-4 clusters \
             by posterior — some inter-donor doublets may be missed for larger k.",
            params.num_clusters
        ));
    }

    let load_elapsed = t_load.elapsed();
    log.ok("Data loading complete", Some(t_load));

    // ── Clustering ────────────────────────────────────────────────────────────
    log.section("CLUSTERING");
    let t_clust = Instant::now();

    let (best_log_probs, plot_data) = run_clustering(
        &log, &params,
        loaded.loci_used,
        &loaded.cell_data,
        &loaded.locus_to_index,
    );
    let clust_elapsed = t_clust.elapsed();

    // ── Compute cluster sizes ─────────────────────────────────────────────────
    let cells_per_cluster = compute_cells_per_cluster(&best_log_probs, params.num_clusters);

    // ── Output ────────────────────────────────────────────────────────────────
    log.section("OUTPUT");
    write_cluster_assignments(&barcodes, &best_log_probs);
    log.ok("Cluster assignments written to stdout", None);


    // ── Diagnostic plots (optional) ────────────────────────────────────────
    let mut plot_files: Vec<String> = Vec::new();
    if let Some(ref plot_dir) = params.plot_dir {
        log.section("DIAGNOSTIC PLOTS");
        let flags = PlotFlags::from_str(&params.plots);
        plot_files = write_plots(&plot_data, plot_dir, &flags);
        for f in &plot_files {
            log.ok(&format!("Wrote {}", f), None);
        }
        if plot_files.is_empty() {
            log.warn("No plots written — check --plots flag or run with --verbose for convergence_curves");
        }
    }

    // ── HTML report ───────────────────────────────────────────────────────────
    if let Some(ref plot_dir) = params.plot_dir {
        log.section("HTML REPORT");
        let zero_frac = loaded.zero_cov_entries as f32
            / (loaded.total_cells.max(1) * loaded.loci_used.max(1)) as f32;
        let report_data = ReportData {
            params:            params.clone(),
            total_cells:       loaded.total_cells,
            loci_raw:          loaded.total_loci_raw,
            loci_used:         loaded.loci_used,
            zero_cov_frac:     zero_frac,
            load_time:         load_elapsed,
            cluster_time:      clust_elapsed,
            total_time:        t_start.elapsed(),
            cells_per_cluster: cells_per_cluster.clone(),
            total_restarts:    plot_data.restart_losses.len(),
            global_best_lp:    plot_data.global_best,
            plot_files:        plot_files.clone(),
            plot_dir:          plot_dir.clone(),
        };
        let report_path = write_report(&report_data, &plot_data);
        log.ok(&format!("Report → {}", report_path), None);
    }

    // ── Cleanup (optional) ────────────────────────────────────────────────────
    if params.cleanup {
        log.section("CLEANUP");
        // In the Rust binary context, the main intermediate files are the
        // BAM/pileup temp files created by the Python pipeline wrapper.
        // The binary itself doesn't create temp files beyond stdout/stderr,
        // so we report the plot_dir disk usage as a courtesy and note that
        // full cleanup of upstream intermediates must be done by the pipeline.
        if let Some(ref plot_dir) = params.plot_dir {
            match std::fs::read_dir(plot_dir) {
                Ok(entries) => {
                    let mut total_bytes = 0u64;
                    let mut n_files     = 0usize;
                    for entry in entries.flatten() {
                        if let Ok(meta) = entry.metadata() {
                            if meta.is_file() {
                                total_bytes += meta.len();
                                n_files     += 1;
                            }
                        }
                    }
                    log.info("plot_dir",    plot_dir);
                    log.info("files",       &n_files.to_string());
                    log.info("total_size",  &format!("{:.1} KB", total_bytes as f64 / 1024.0));
                    log.ok("Cleanup complete — intermediate BAM/pileup files must be removed by the pipeline wrapper", None);
                }
                Err(e) => log.warn(&format!("Could not read plot_dir for cleanup report: {}", e)),
            }
        } else {
            log.info("cleanup", "no plot_dir set — nothing to report");
        }
    }

    eprintln!("\n{}  Total runtime: {}{}", "\x1b[2m", log.total_elapsed(), "\x1b[0m");
}

// ── run_clustering ────────────────────────────────────────────────────────────
//
// Manages the outer restart loop (parallel over threads), the optional
// souporcell3 multi-run logic, and the global-best bookkeeping.
//
// Returns the best per-cell log-probability matrix found across all restarts.
fn run_clustering(
    log:           &Logger,
    params:        &Params,
    loci_used:     usize,
    cell_data:     &[CellData],
    locus_to_index: &HashMap<usize, usize>,
) -> (Vec<Vec<f32>>, PlotData) {
    let seed_array: [u8; 32] = {
        let mut a = [0u8; 32];
        // Spread the u64 seed across all 32 bytes for StdRng
        let bytes = params.seed.to_le_bytes();
        for i in 0..32 { a[i] = bytes[i % 8]; }
        a
    };
    let mut global_rng: StdRng = SeedableRng::from_seed(seed_array);

    let solves_per_thread =
        ((params.restarts as f32) / (params.threads as f32)).ceil() as usize;

    let mut global_best_lp        = f32::NEG_INFINITY;
    let mut global_best_log_probs: Vec<Vec<f32>> = Vec::new();
    let mut global_best_centers:   Vec<Vec<f32>> = Vec::new();
    let mut total_restarts_run     = 0usize;

    // Plot data collector — shared across threads via Arc<Mutex>
    let shared_plot = Arc::new(Mutex::new(PlotData::new(params.num_clusters)));

    // souporcell3 runs up to 3 passes; default runs exactly once (run 0 only)
    let max_runs = if params.souporcell3 { 3 } else { 1 };

    for run in 0..max_runs {
        let t_run = Instant::now();

        // Bad-cluster detection for runs > 0
        let bad_clusters = bad_cluster_detection(run, params.num_clusters, &global_best_log_probs);
        if run > 0 && bad_clusters.is_empty() {
            log.ok(
                &format!("No outlier clusters detected on run {} — converged early", run + 1),
                None
            );
            break;
        }

        if run > 0 {
            log.section(&format!("SOUPORCELL3 — RUN {}/3  (reinitialising {} clusters)",
                run + 1, bad_clusters.len()));
        }

        // Build per-thread seed array before spawning
        let thread_seeds: Vec<[u8; 32]> = (0..params.threads)
            .map(|_| new_seed(&mut global_rng))
            .collect();

        let mut threads: Vec<ThreadData> = thread_seeds.into_iter().enumerate()
            .map(|(i, seed)| ThreadData::from_seed(seed, solves_per_thread, i))
            .collect();

        // Parallel restart loop — shared_plot_clone is captured by each thread
        let shared_plot_clone = Arc::clone(&shared_plot);
        threads.par_iter_mut().for_each(|td| {
            let plot_ref = Arc::clone(&shared_plot_clone);
            for iteration in 0..td.solves_per_thread {
                let (cluster_centers, locked) = prepare_centers(
                    run, loci_used, cell_data, params, &mut td.rng,
                    locus_to_index, &bad_clusters, &global_best_centers,
                );
                let mut cc = cluster_centers;  // need mut for em/khm

                let (loss, log_probs, _stats, ts_losses) = match params.clustering_method {
                    ClusterMethod::EM  => em(loci_used, &mut cc, cell_data, params,
                                             iteration, td.thread_num, locked),
                    ClusterMethod::KHM => khm(loci_used, &mut cc, cell_data, params,
                                              iteration, td.thread_num, locked),
                };
                // record temp-step data for annealing profile plot
                if let Ok(mut pd) = plot_ref.lock() {
                    for rec in ts_losses { pd.record_temp_step(rec.0, rec.1, rec.2, rec.3, rec.4); }
                }

                if loss > td.best_total_log_probability {
                    td.best_total_log_probability = loss;
                    td.best_log_probabilities     = log_probs;
                    td.cluster_centers            = cc;
                }

                log_restart(td.thread_num, iteration, loss, td.best_total_log_probability);
                td.restarts_completed += 1;
                // record for plot 1 (restart landscape) and plot 6 (thread efficiency)
                if let Ok(mut pd) = plot_ref.lock() {
                    pd.record_restart(td.thread_num, iteration, loss);
                }
            }
        });

        // Merge thread results into global best
        for td in &threads {
            total_restarts_run += td.restarts_completed;
            if td.best_total_log_probability > global_best_lp {
                global_best_lp        = td.best_total_log_probability;
                global_best_log_probs = td.best_log_probabilities.clone();
                global_best_centers   = td.cluster_centers.clone();
            }
        }

        log.ok(
            &format!("Run {} complete — best log-prob: {:.4}", run + 1, global_best_lp),
            Some(t_run)
        );
    }

    // Populate cell assignment data for plots 4 & 5
    {
        let mut pd = shared_plot.lock().unwrap();
        pd.global_best = global_best_lp;
        for lps in &global_best_log_probs {
            if lps.is_empty() { continue; }
            let best = lps.iter().enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(i, _)| i).unwrap_or(0);
            pd.record_cell(best, lps.clone());
        }
    }

    // Unwrap Arc → Mutex → PlotData (must happen before log_run_summary)
    let plot_data = match Arc::try_unwrap(shared_plot) {
        Ok(mutex)  => mutex.into_inner().unwrap_or_else(|_| PlotData::new(params.num_clusters)),
        Err(arc)   => {
            eprintln!("Warning: plot data Arc has multiple owners — cloning");
            arc.lock().unwrap().clone_basic()
        }
    };

    // ── Final summary ──────────────────────────────────────────────────────────
    let cells_per_cluster_pre = compute_cells_per_cluster(
        &global_best_log_probs, params.num_clusters
    );
    log_run_summary(
        plot_data.global_best,
        params.num_clusters,
        &cells_per_cluster_pre,
        cell_data.len(),
        total_restarts_run,
        &log.total_elapsed(),
    );

    (global_best_log_probs, plot_data)
}

// ── prepare_centers ───────────────────────────────────────────────────────────
//
// For run 0: initialise all cluster centres fresh.
// For run > 0 (souporcell3): keep good centres locked, reinitialise bad ones.
fn prepare_centers(
    run:              usize,
    loci_used:        usize,
    cell_data:        &[CellData],
    params:           &Params,
    rng:              &mut StdRng,
    locus_to_index:   &HashMap<usize, usize>,
    bad_clusters:     &[usize],
    best_centers:     &[Vec<f32>],
) -> (Vec<Vec<f32>>, Vec<usize>) {
    if run == 0 {
        let centers = init_cluster_centers(
            loci_used, cell_data, params, rng, locus_to_index, params.num_clusters
        );
        (centers, vec![])
    } else {
        let reinit = init_cluster_centers(
            loci_used, cell_data, params, rng, locus_to_index, bad_clusters.len()
        );
        let mut centers = best_centers.to_vec();
        for (idx, &cluster_idx) in bad_clusters.iter().enumerate() {
            centers[cluster_idx] = reinit[idx].clone();
        }
        let locked: Vec<usize> = (0..params.num_clusters)
            .filter(|c| !bad_clusters.contains(c))
            .collect();
        (centers, locked)
    }
}

// ── bad_cluster_detection ─────────────────────────────────────────────────────
//
// Identifies under- or over-populated clusters for souporcell3 reinitialization.
// Returns empty vec on run 0 or when num_clusters < 16 (no detection needed).
//
// Enhancement: now logs each cluster's status via log_cluster_analysis()
// instead of bare eprint! calls.
fn bad_cluster_detection(
    run:                  usize,
    num_clusters:         usize,
    best_log_probabilities: &[Vec<f32>],
) -> Vec<usize> {
    if num_clusters < 16 || run == 0 {
        return vec![];
    }

    // Count cells assigned to each cluster
    let mut assigned: Vec<(usize, usize)> = (0..num_clusters).map(|i| (i, 0)).collect();
    for lps in best_log_probabilities {
        let best = lps.iter().enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(i, _)| i)
            .unwrap_or(0);
        assigned[best].1 += 1;
    }
    assigned.sort_by_key(|(_, count)| *count);

    // Reinitialisation window shrinks with each run
    let run_value      = 3usize.saturating_sub(run);
    let lo_cut         = run_value * (assigned.len() / 16);
    let hi_cut         = (16 - run_value) * (assigned.len() / 16);

    let mut replace = Vec::new();
    eprintln!("CLUSTER_ANALYSIS\trun={}", run);
    for (rank, &(cluster, count)) in assigned.iter().enumerate() {
        let is_outlier = rank < lo_cut || rank > hi_cut;
        log_cluster_analysis(run, cluster, count, is_outlier);
        if is_outlier && !replace.contains(&cluster) {
            replace.push(cluster);
        }
    }
    replace
}

// ── Utilities ─────────────────────────────────────────────────────────────────

fn compute_cells_per_cluster(
    best_log_probs: &[Vec<f32>],
    num_clusters:   usize,
) -> Vec<usize> {
    let mut counts = vec![0usize; num_clusters];
    for lps in best_log_probs {
        if lps.is_empty() { continue; }
        let best = lps.iter().enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(i, _)| i)
            .unwrap_or(0);
        counts[best] += 1;
    }
    counts
}

fn new_seed(rng: &mut StdRng) -> [u8; 32] {
    let mut seed = [0u8; 32];
    for b in &mut seed { *b = rng.gen::<u8>(); }
    seed
}

