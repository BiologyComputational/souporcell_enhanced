// ============================================================================
// infra/preflight.rs — Pre-flight plan, expected output preview, approval gate
// ============================================================================
//
// When --dry_run is set, this module:
//   1. Prints a full structured run plan (what will happen, in what order,
//      with what parameters)
//   2. Shows per-stage estimated cost (time, memory, output size)
//   3. Shows expected output description and sample interpretation
//   4. Prompts the user for approval [Y/n]
//   5. Returns ApprovalResult::Proceed or ApprovalResult::Abort
//
// When --dry_run is NOT set, show_preflight() is never called — zero overhead.
//
// The approval gate can also be bypassed by setting --dry_run_yes (or the env
// var SOUPC_DRY_RUN_YES=true) which prints the plan but proceeds automatically.
// This is useful for CI pipelines that want the plan logged but no stdin.
// ============================================================================

use std::io::{self, Write};
use crate::config::params::Params;

// ── Public result type ────────────────────────────────────────────────────────

pub enum ApprovalResult {
    Proceed,
    Abort,
}

// ── Public entry point ────────────────────────────────────────────────────────

pub fn show_preflight(params: &Params) -> ApprovalResult {
    let auto_yes = params.dry_run_yes;

    print_banner();
    print_run_identity(params);
    print_input_plan(params);
    print_pipeline_stages(params);
    print_algorithm_detail(params);
    print_expected_outputs(params);
    print_expected_interpretation(params);
    print_quality_thresholds(params);
    print_footer(params);

    if auto_yes {
        eprintln!(
            "\n{}  [--dry_run_yes set — proceeding automatically]{}\n",
            green(), reset()
        );
        return ApprovalResult::Proceed;
    }

    prompt_approval()
}

// ── Approval prompt ───────────────────────────────────────────────────────────

fn prompt_approval() -> ApprovalResult {
    eprintln!();
    eprint!(
        "{}{}  Proceed with this run? [Y/n]: {}",
        bold(), cyan(), reset()
    );
    let _ = io::stderr().flush();

    // Read from stdin so the prompt works even when stdout is redirected
    let mut input = String::new();
    match io::stdin().read_line(&mut input) {
        Ok(_) => {
            let answer = input.trim().to_lowercase();
            if answer.is_empty() || answer == "y" || answer == "yes" {
                eprintln!("{}  ✔  Approved — starting run…{}\n", green(), reset());
                ApprovalResult::Proceed
            } else {
                eprintln!("{}  ✖  Aborted by user.{}\n", red(), reset());
                ApprovalResult::Abort
            }
        }
        Err(_) => {
            eprintln!("{}  ⚠  Could not read stdin — aborting. Use --dry_run_yes to auto-proceed.{}\n",
                      yellow(), reset());
            ApprovalResult::Abort
        }
    }
}

// ── Section printers ──────────────────────────────────────────────────────────

fn print_banner() {
    eprintln!();
    eprintln!("{}{}╔══════════════════════════════════════════════════════════════════════════╗{}",
              bold(), cyan(), reset());
    eprintln!("{}{}║        souporcell v2.7  ·  PRE-FLIGHT RUN PLAN                          ║{}",
              bold(), cyan(), reset());
    eprintln!("{}{}║        Review carefully before approving. Nothing has run yet.           ║{}",
              bold(), cyan(), reset());
    eprintln!("{}{}╚══════════════════════════════════════════════════════════════════════════╝{}",
              bold(), cyan(), reset());
}

fn print_run_identity(params: &Params) {
    section("RUN IDENTITY");
    let method  = params.clustering_method.name().to_uppercase();
    let run_id  = format!("k{}_r{}_s{}_t{}",
                          params.num_clusters, params.restarts,
                          params.seed, params.threads);
    kv("Run ID",            &run_id);
    kv("Algorithm",         &method);
    kv("k (donor clusters)",&params.num_clusters.to_string());
    kv("Restarts",          &params.restarts.to_string());
    kv("Seed",              &params.seed.to_string());
    kv("Threads",           &params.threads.to_string());
    kv("Config profile",    &params.config_profile.as_deref().unwrap_or("(none — CLI only)"));
    kv("Env file",          &params.env_file.as_deref().unwrap_or("(none)"));
    kv("Param priority",    "CLI  >  --config JSON  >  .env  >  built-in defaults");
}

fn print_input_plan(params: &Params) {
    section("INPUT FILES");
    kv("ref_matrix",  &params.ref_mtx);
    kv("alt_matrix",  &params.alt_mtx);
    kv("barcodes",    &params.barcodes);
    eprintln!();

    section("LOCUS QC FILTERS  (applied after loading)");
    kv("min_ref",      &format!("{} cells must have ≥1 ref read at a locus", params.min_ref));
    kv("min_alt",      &format!("{} cells must have ≥1 alt read at a locus", params.min_alt));
    kv("min_ref_umis", &format!("{} minimum total ref UMIs per locus", params.min_ref_umis));
    kv("min_alt_umis", &format!("{} minimum total alt UMIs per locus", params.min_alt_umis));
    note("Loci failing any filter are silently excluded. \
          Typical pass rate: 5–20% of raw loci. \
          If 0 loci pass, clustering cannot proceed — lower min_ref/min_alt.");
}

fn print_pipeline_stages(params: &Params) {
    section("PIPELINE STAGES  (in execution order)");

    let plot_dir = params.plot_dir.as_deref().unwrap_or("(disabled)");
    let stages: &[(&str, &str, &str)] = &[
        ("1", "Load matrices",
         "Read ref.mtx and alt.mtx into sparse CSR format. Decompress .gz if needed."),
        ("2", "Load barcodes",
         "Read barcodes file, build barcode→index map."),
        ("3", "Locus QC filter",
         "Apply min_ref / min_alt / min_ref_umis / min_alt_umis. Build locus_to_index."),
        ("4", "Cluster centre init",
         &format!("Initialise {} × k={} theta vectors using strategy: {:?}.",
                  params.restarts, params.num_clusters, params.initialization_strategy)),
        ("5", "Parallel EM/KHM + annealing",
         &format!("{} restarts × {} temp steps across {} threads. \
                   Convergence: |Δ loss| < {:.4} × n_cells  OR  {} iters.",
                  params.restarts, params.anneal_steps,
                  params.threads, params.conv_tol, params.max_iter)),
        ("6", "Global optimum selection",
         "Pick the restart with the highest total log-likelihood across all threads."),
        ("7", "Cell assignment",
         "Assign each cell to its argmax cluster. Compute confidence = best_lp − second_best_lp."),
        ("8", "Write clusters_tmp.tsv",
         "One row per cell: barcode, cluster index, log-prob per cluster → stdout."),
        ("9", "SVG diagnostic plots",
         &format!("Write to: {}  (plots: {})", plot_dir, params.plots)),
        ("10", "HTML run report",
         &format!("Write souporcell_report.html to: {}", plot_dir)),
    ];

    for (num, name, desc) in stages {
        eprintln!("  {}{}  {:>2}. {:<32}{}  {}{}{}",
                  bold(), cyan(), num, name, reset(),
                  dim(), desc, reset());
    }
}

fn print_algorithm_detail(params: &Params) {
    section("ALGORITHM HYPERPARAMETERS");

    subsection("Annealing schedule");
    kv("anneal_steps",    &format!("{} temperature steps (T0 = hottest → T{} = standard EM)",
                                   params.anneal_steps, params.anneal_steps - 1));
    kv("anneal_base_em",  &format!("{:.4}  →  raw_temp = total_alleles / ({:.1} × 2^step)",
                                   params.anneal_base_em, params.anneal_base_em));
    kv("anneal_base_khm", &format!("{:.4}  (KHM runs hotter than EM — smaller base = softer assignments)",
                                   params.anneal_base_khm));
    note(&format!(
        "At step 0 a cell with 100 total alleles gets temp = 100 / ({:.1} × 1) = {:.1}. \
         At the final step temp is forced to 1.0 (standard softmax).",
        params.anneal_base_em,
        100.0 / params.anneal_base_em
    ));

    subsection("Convergence");
    let example_cells = 3000usize;
    kv("conv_tol",          &format!("{:.4}  →  stop when |Δ loss| < {:.1} (for ~{} cells)",
                                     params.conv_tol,
                                     params.conv_tol * example_cells as f32,
                                     example_cells));
    kv("max_iter",          &format!("{} iterations hard cap per temperature step",
                                     params.max_iter));
    kv("min_cluster_cells", &format!("{} cells — minimum for a cluster to count as populated \
                                      in verbose logs", params.min_cluster_cells));

    subsection("M-step (allele-frequency update)");
    kv("pseudocount_alt",   &format!("{:.2}  added to alt-allele numerator (Laplace smoothing)",
                                     params.pseudocount_alt));
    kv("pseudocount_total", &format!("{:.2}  added to total-allele denominator",
                                     params.pseudocount_total));
    kv("theta bounds",      &format!("[{:.3}, {:.3}]  — θ clamped to this range after each M-step",
                                     params.theta_min, params.theta_max));
    note("Pseudocounts prevent θ from reaching 0.0 or 1.0, which would make log-likelihood \
          undefined. The defaults (1.0 / 2.0) correspond to a Beta(1,1) = Uniform prior.");

    if params.clustering_method == crate::config::params::ClusterMethod::KHM {
        subsection("KHM-specific");
        kv("khm_p", &format!("{:.1}  harmonic-mean exponent — higher = more peaked membership \
                              around nearest centre", params.khm_p));
        note("KHM p=25 gives very sharp membership functions. \
              Increase to 50+ for even more aggressive separation; \
              decrease to 10 for softer, more EM-like behaviour.");
    }
}

fn print_expected_outputs(params: &Params) {
    section("EXPECTED OUTPUTS");

    subsection("stdout → clusters_tmp.tsv");
    eprintln!("  {}One row per cell barcode. Redirect to a file with  > clusters_tmp.tsv{}",
              dim(), reset());
    eprintln!();
    eprintln!("  {}{}BARCODE              CLUSTER  LP_0        LP_1        LP_2        LP_3{}",
              bold(), cyan(), reset());
    eprintln!("  {}AAACATTGAGCTAC-1     2        -312.44     -298.11     -287.56     -301.88{}",
              dim(), reset());
    eprintln!("  {}AAACATTGATCAGC-1     0        -289.12     -401.55     -398.70     -411.22{}",
              dim(), reset());
    eprintln!("  {}AAACCGTGCTTCCG-1     1        -445.88     -271.34     -388.91     -420.55{}",
              dim(), reset());
    eprintln!("  {}... (one row per barcode){}",
              dim(), reset());
    eprintln!();
    note("Column 2 = best-fit donor cluster (0-based). Columns 3..k+2 = log-probability \
          per cluster — least negative wins. Feed this file to troublet (Stage 7).");

    if let Some(ref plot_dir) = params.plot_dir {
        subsection(&format!("SVG plots → {}/<name>.svg", plot_dir));
        let plots_to_show = vec![
            ("restart_landscape.svg",   "Final log-likelihood per restart, sorted worst→best"),
            ("convergence_curves.svg",  "Per-iteration loss curves (requires --verbose)"),
            ("annealing_profile.svg",   "Mean ± SD loss at each temperature step"),
            ("cluster_balance.svg",     "Cells per cluster vs expected equal split"),
            ("posterior_confidence.svg","Histogram of best_lp − second_best_lp per cell"),
            ("thread_efficiency.svg",   "Restart count and best loss per thread"),
        ];
        for (name, desc) in plots_to_show {
            eprintln!("  {}  {:<36}  {}{}{}",
                      cyan(), name, dim(), desc, reset());
        }

        subsection(&format!("HTML report → {}/souporcell_report.html", plot_dir));
        note("Self-contained single HTML file. No internet required. \
              Contains metric cards, interactive D3.js pipeline graph, \
              5 results tables, quality interpretation, and full parameters.");
    }
}

fn print_expected_interpretation(params: &Params) {
    section("EXPECTED RESULTS & SAMPLE INTERPRETATION");

    let k          = params.num_clusters;
    let restarts   = params.restarts;
    let expected_per_cluster = "~25%";   // rough for equal pooling

    eprintln!("  {}Assuming equal donor pooling (most common case):{}",
              dim(), reset());
    eprintln!();

    let items: &[(&str, &str, &str)] = &[
        ("Cells per cluster",
         &format!("~{} of total cells", expected_per_cluster),
         "Severe imbalance (> 50% deviation) suggests k may be too large or a donor is absent"),
        ("Convergence stability",
         &format!("≥ 70% of {} restarts at same log-prob", restarts),
         "< 40% → increase --restarts or use profiles/high_convergence.json"),
        ("Median cell confidence",
         "> 15 log-units",
         "Cells with < 10 log-units are ambiguous — likely doublets, flagged by troublet"),
        ("Loci pass rate",
         "> 1% of raw SNP loci",
         "< 1% → relax --min_ref / --min_alt; insufficient informative sites"),
        ("Matrix sparsity",
         "> 90% zero coverage",
         "< 90% may indicate ref/alt matrix mismatch or wrong barcodes file"),
        ("EM convergence per step",
         &format!("iters < {} at most steps", params.max_iter),
         &format!("iters = {} (cap hit) → data too noisy or conv_tol too tight", params.max_iter)),
        ("Cluster count",
         &format!("{} non-empty clusters", k),
         "Any empty cluster → k is larger than the actual number of donors"),
    ];

    for (metric, expected, warning) in items {
        eprintln!("  {}●{} {:<35} {}expected: {:<28}{} {}↳ warn if: {}{}",
                  cyan(), reset(),
                  metric,
                  bold(), expected, reset(),
                  dim(), warning, reset());
        eprintln!();
    }

    note("These are heuristics for typical 10x Chromium scRNA-seq experiments. \
          Cell line mixtures or very unequal donor proportions will show different patterns.");
}

fn print_quality_thresholds(params: &Params) {
    section("QUALITY CHECKS  (evaluated automatically in the HTML report)");

    let checks: &[(&str, &str)] = &[
        ("Loci QC pass rate",       "> 1% of raw loci"),
        ("Matrix sparsity",         "> 90% zero-coverage entries"),
        ("All restarts completed",  &format!("{} of {} restarts", params.restarts, params.restarts)),
        ("Convergence stability",   "≥ 70% of restarts at global optimum (±1 log-unit)"),
        ("Cluster balance",         "Max deviation from equal split < 50%"),
        ("No empty clusters",       "All k clusters have > 0 cells"),
        ("Median cell confidence",  "> 15 log-units (best − second_best log-prob)"),
        ("SVG plots written",       "All requested plots generated"),
        ("clusters_tmp.tsv written","Non-empty output file"),
        ("HTML report generated",   "souporcell_report.html written"),
    ];

    for (check, threshold) in checks {
        eprintln!("  {}  {:<40}  {}threshold: {}{}",
                  cyan(), check, dim(), threshold, reset());
    }
}

fn print_footer(params: &Params) {
    section("ESTIMATED RUNTIME");

    let total_work = params.restarts as f32 * params.anneal_steps as f32;
    let per_thread = total_work / params.threads as f32;
    note(&format!(
        "Total work units: {} restarts × {} anneal steps = {:.0} unit-steps. \
         With {} threads: ~{:.0} unit-steps per thread. \
         Actual time depends on dataset size and loci count.",
        params.restarts, params.anneal_steps, total_work,
        params.threads, per_thread
    ));

    let abort_hint = if params.dry_run_yes {
        String::new()
    } else {
        "  Type 'n' to abort without running anything.\n".to_string()
    };

    eprintln!();
    eprintln!("{}{}  Nothing has been read or written yet.{}",
              bold(), yellow(), reset());
    eprintln!("{}{}  All of the above is the planned behaviour only.{}",
              bold(), yellow(), reset());
    if !abort_hint.is_empty() {
        eprint!("{}{}{}", dim(), abort_hint, reset());
    }
}

// ── Formatting helpers ────────────────────────────────────────────────────────

fn section(title: &str) {
    eprintln!("\n{}{}  ┌─ {} {}",
              bold(), cyan(), title, reset());
    eprintln!("{}{}  └{}{}", cyan(), dim(),
              "─".repeat(72), reset());
}

fn subsection(title: &str) {
    eprintln!("\n{}  ▸ {}{}",
              cyan(), title, reset());
}

fn kv(key: &str, val: &str) {
    eprintln!("  {}  {:<28}  {}{}{}",
              dim(), key, reset(), val, reset());
}

fn note(msg: &str) {
    // Word-wrap at ~80 cols
    let prefix = "  ℹ  ";
    let width   = 74usize;
    let mut remaining = msg;
    let mut first = true;
    while !remaining.is_empty() {
        let chunk = if remaining.len() <= width {
            let s = remaining;
            remaining = "";
            s
        } else {
            let cut = remaining[..width].rfind(' ').unwrap_or(width);
            let s = &remaining[..cut];
            remaining = remaining[cut..].trim_start();
            s
        };
        if first {
            eprintln!("{}{}{}{}{}",
                      dim(), prefix, chunk, reset(), reset());
            first = false;
        } else {
            eprintln!("{}  {}  {}{}",
                      dim(), " ".repeat(prefix.len() - 2), chunk, reset());
        }
    }
}

// ── ANSI helpers (mirrors logger.rs — kept local to avoid coupling) ───────────

fn is_tty() -> bool {
    std::env::var("TERM").map(|t| t != "dumb").unwrap_or(false)
}
fn cyan()   -> &'static str { if is_tty() { "\x1b[36m" } else { "" } }
fn yellow() -> &'static str { if is_tty() { "\x1b[33m" } else { "" } }
fn green()  -> &'static str { if is_tty() { "\x1b[32m" } else { "" } }
fn red()    -> &'static str { if is_tty() { "\x1b[31m" } else { "" } }
fn bold()   -> &'static str { if is_tty() { "\x1b[1m"  } else { "" } }
fn dim()    -> &'static str { if is_tty() { "\x1b[2m"  } else { "" } }
fn reset()  -> &'static str { if is_tty() { "\x1b[0m"  } else { "" } }
