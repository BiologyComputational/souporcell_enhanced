// ============================================================================
// config/params.rs — CLI parsing, config profile merging, Params struct
// ============================================================================
//
// Parameter resolution order (highest priority first):
//   1. CLI flags          --restarts 500
//   2. --config JSON      --config profiles/k4.json
//   3. .souporcell.env    auto-discovered in CWD or ~/.souporcell.env
//   4. Built-in defaults  (same values as previous hardcoded constants)
//
// Nothing in core/, domain/, or analysis/ contains a hardcoded algorithmic
// constant — every tunable value flows through Params.
// ============================================================================

use clap::App;
use crate::config::env_loader::{load_env, EnvMap};
use crate::config::config_loader::{load_config, ConfigProfile};

// ── Algorithm strategy enums ─────────────────────────────────────────────────

#[derive(Clone, Debug, PartialEq)]
pub enum ClusterInit {
    KmeansPP,
    RandomUniform,
    RandomAssignment,
    MiddleVariance,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ClusterMethod {
    EM,
    KHM,
}

impl ClusterMethod {
    pub fn name(&self) -> &'static str {
        match self {
            ClusterMethod::EM  => "em",
            ClusterMethod::KHM => "khm",
        }
    }
}

// ── Params ────────────────────────────────────────────────────────────────────

#[derive(Clone)]
pub struct Params {
    // ── I/O paths ────────────────────────────────────────────────────────────
    pub ref_mtx:                      String,
    pub alt_mtx:                      String,
    pub barcodes:                     String,

    // ── Clustering ───────────────────────────────────────────────────────────
    pub num_clusters:                 usize,
    pub restarts:                     u32,
    pub seed:                         u64,
    pub clustering_method:            ClusterMethod,
    pub initialization_strategy:      ClusterInit,
    pub souporcell3:                  bool,

    // ── Locus quality filters ────────────────────────────────────────────────
    pub min_alt:                      u32,
    pub min_ref:                      u32,
    pub min_alt_umis:                 u32,
    pub min_ref_umis:                 u32,

    // ── Parallelism ──────────────────────────────────────────────────────────
    pub threads:                      usize,

    // ── Optional priors / genotype inputs ────────────────────────────────────
    pub known_cell_assignments:       Option<String>,
    pub known_genotypes:              Option<String>,
    pub known_genotypes_sample_names: Vec<String>,

    // ── Logging / verbosity ──────────────────────────────────────────────────
    pub verbose:                      bool,
    #[allow(dead_code)]
    pub progress_interval:            usize,

    // ── Diagnostics / plots ──────────────────────────────────────────────────
    pub plot_dir:                     Option<String>,
    pub plots:                        String,
    pub cleanup:                      bool,

    // ── Annealing schedule ────────────────────────────────────────────────────
    /// Number of temperature steps in the annealing schedule (default: 9)
    pub anneal_steps:                 usize,
    /// EM temperature divisor base: raw_temp = total_alleles / (base × 2^step) (default: 20.0)
    pub anneal_base_em:               f32,
    /// KHM temperature divisor base (default: 0.5)
    pub anneal_base_khm:              f32,

    // ── EM / KHM convergence ──────────────────────────────────────────────────
    /// Convergence tolerance: stop when |Δ loss| < tol × n_cells (default: 0.01)
    pub conv_tol:                     f32,
    /// Maximum EM/KHM iterations per temperature step (default: 1000)
    pub max_iter:                     usize,
    /// Minimum cells assigned to a cluster for it to count as "populated" (default: 200)
    pub min_cluster_cells:            usize,

    // ── M-step pseudocounts ───────────────────────────────────────────────────
    /// Pseudocount added to alt-allele sum in M-step (default: 1.0)
    pub pseudocount_alt:              f32,
    /// Pseudocount added to total-allele denominator in M-step (default: 2.0)
    pub pseudocount_total:            f32,

    // ── Allele-frequency bounds ───────────────────────────────────────────────
    /// Lower bound for θ after M-step clamp (default: 0.01)
    pub theta_min:                    f32,
    /// Upper bound for θ after M-step clamp (default: 0.99)
    pub theta_max:                    f32,

    // ── KHM-specific ──────────────────────────────────────────────────────────
    /// KHM harmonic-mean exponent p (default: 25.0)
    pub khm_p:                        f32,

    // ── Config provenance (informational) ────────────────────────────────────
    pub config_profile:               Option<String>,
    pub env_file:                     Option<String>,

    // ── Pre-flight / dry-run ──────────────────────────────────────────────────
    /// Print the full run plan and wait for [Y/n] approval before proceeding.
    pub dry_run:                      bool,
    /// Requires dry_run=true. Print the plan and proceed without prompting.
    /// For CI pipelines. Also settable via SOUPC_DRY_RUN_YES=true.
    pub dry_run_yes:                  bool,
}

// ── load_params ───────────────────────────────────────────────────────────────

pub fn load_params() -> Params {
    let yaml    = load_yaml!("params.yml");
    let matches = App::from_yaml(yaml).get_matches();

    // ── Load env + config layers ──────────────────────────────────────────────
    let env_path    = matches.value_of("env");
    let config_path = matches.value_of("config");

    let env_raw  = load_env(env_path);
    let env_file = env_path.map(str::to_string).or_else(find_env_path);
    let env      = EnvMap(env_raw);

    let cfg: Option<ConfigProfile> = config_path.map(load_config);

    // ── Resolve macros (CLI > config > env > default) ────────────────────────

    macro_rules! resolve_str {
        ($cli:expr, $ck:expr, $ek:expr, $def:expr) => {{
            if let Some(v) = matches.value_of($cli) {
                v.to_string()
            } else if let Some(ref c) = cfg {
                c.str_val($ck).map(str::to_string)
                    .or_else(|| env.str($ek).map(str::to_string))
                    .unwrap_or_else(|| $def.to_string())
            } else {
                env.str($ek).map(str::to_string)
                    .unwrap_or_else(|| $def.to_string())
            }
        }};
    }

    macro_rules! resolve_parse {
        ($cli:expr, $ck:expr, $ek:expr, $T:ty, $def:expr) => {{
            if let Some(raw) = matches.value_of($cli) {
                raw.parse::<$T>()
                   .unwrap_or_else(|_| panic!("--{} must be a valid {}", $cli, stringify!($T)))
            } else if let Some(ref c) = cfg {
                c.parse::<$T>($ck)
                    .or_else(|| env.parse::<$T>($ek))
                    .unwrap_or($def)
            } else {
                env.parse::<$T>($ek).unwrap_or($def)
            }
        }};
    }

    macro_rules! resolve_bool_flag {
        ($cli:expr, $ck:expr, $ek:expr, $def:expr) => {{
            if matches.is_present($cli) { true }
            else if let Some(ref c) = cfg {
                c.bool_val($ck).or_else(|| env.bool($ek)).unwrap_or($def)
            } else {
                env.bool($ek).unwrap_or($def)
            }
        }};
    }

    macro_rules! resolve_bool_val {
        ($cli:expr, $ck:expr, $ek:expr, $def:expr) => {{
            if let Some(raw) = matches.value_of($cli) {
                raw.parse::<bool>()
                   .unwrap_or_else(|_| panic!("--{} must be true or false", $cli))
            } else if let Some(ref c) = cfg {
                c.bool_val($ck).or_else(|| env.bool($ek)).unwrap_or($def)
            } else {
                env.bool($ek).unwrap_or($def)
            }
        }};
    }

    macro_rules! resolve_opt_str {
        ($cli:expr, $ck:expr, $ek:expr) => {{
            if let Some(v) = matches.value_of($cli) {
                Some(v.to_string())
            } else if let Some(ref c) = cfg {
                c.str_val($ck).map(str::to_string)
                    .or_else(|| env.str($ek).map(str::to_string))
            } else {
                env.str($ek).map(str::to_string)
            }
        }};
    }

    // ── I/O paths ─────────────────────────────────────────────────────────────
    let ref_mtx  = expand_tilde(&resolve_str!("ref_matrix", "ref_matrix", "SOUPC_REF_MATRIX", ""));
    let alt_mtx  = expand_tilde(&resolve_str!("alt_matrix", "alt_matrix", "SOUPC_ALT_MATRIX", ""));
    let barcodes = expand_tilde(&resolve_str!("barcodes",   "barcodes",   "SOUPC_BARCODES",   ""));

    if ref_mtx.is_empty()  { eprintln!("Error: ref_matrix required (-r / SOUPC_REF_MATRIX)");  std::process::exit(1); }
    if alt_mtx.is_empty()  { eprintln!("Error: alt_matrix required (-a / SOUPC_ALT_MATRIX)");  std::process::exit(1); }
    if barcodes.is_empty() { eprintln!("Error: barcodes required   (-b / SOUPC_BARCODES)");     std::process::exit(1); }

    // ── Clustering ────────────────────────────────────────────────────────────
    let num_clusters = resolve_parse!("num_clusters", "num_clusters", "SOUPC_NUM_CLUSTERS", usize, 0);
    if num_clusters < 2 {
        eprintln!("Error: num_clusters must be >= 2 (got {}). Set via -k / SOUPC_NUM_CLUSTERS.", num_clusters);
        std::process::exit(1);
    }

    let restarts = resolve_parse!("restarts", "restarts", "SOUPC_RESTARTS", u32,   100);
    let seed     = resolve_parse!("seed",     "seed",     "SOUPC_SEED",     u64,   4);
    let threads  = resolve_parse!("threads",  "threads",  "SOUPC_THREADS",  usize, 1);
    assert!(threads >= 1, "threads must be >= 1");

    let method_str = resolve_str!("clustering_method", "clustering_method",
                                  "SOUPC_CLUSTERING_METHOD", "em");
    let clustering_method = match method_str.as_str() {
        "em"  => ClusterMethod::EM,
        "khm" => ClusterMethod::KHM,
        other => panic!("Unknown clustering_method '{}'. Use 'em' or 'khm'.", other),
    };

    let init_str = resolve_str!("initialization_strategy", "initialization_strategy",
                                "SOUPC_INIT_STRATEGY", "random_uniform");
    let initialization_strategy = match init_str.as_str() {
        "kmeans_pp" | "kmeans++"  => ClusterInit::KmeansPP,
        "random_uniform"          => ClusterInit::RandomUniform,
        "random_cell_assignment"  => ClusterInit::RandomAssignment,
        "middle_variance"         => ClusterInit::MiddleVariance,
        other => panic!("Unknown initialization_strategy '{}'.", other),
    };

    let souporcell3 = resolve_bool_val!("souporcell3", "souporcell3", "SOUPC_SOUPORCELL3", false);

    // ── Locus QC ──────────────────────────────────────────────────────────────
    let min_alt      = resolve_parse!("min_alt",      "min_alt",      "SOUPC_MIN_ALT",      u32, 4);
    let min_ref      = resolve_parse!("min_ref",      "min_ref",      "SOUPC_MIN_REF",      u32, 4);
    let min_alt_umis = resolve_parse!("min_alt_umis", "min_alt_umis", "SOUPC_MIN_ALT_UMIS", u32, 0);
    let min_ref_umis = resolve_parse!("min_ref_umis", "min_ref_umis", "SOUPC_MIN_REF_UMIS", u32, 0);

    // ── Annealing schedule ────────────────────────────────────────────────────
    let anneal_steps   = resolve_parse!("anneal_steps",   "anneal_steps",   "SOUPC_ANNEAL_STEPS",   usize, 9);
    let anneal_base_em = resolve_parse!("anneal_base_em", "anneal_base_em", "SOUPC_ANNEAL_BASE_EM", f32,   20.0);
    let anneal_base_khm= resolve_parse!("anneal_base_khm","anneal_base_khm","SOUPC_ANNEAL_BASE_KHM",f32,   0.5);

    // ── EM / KHM convergence ──────────────────────────────────────────────────
    let conv_tol         = resolve_parse!("conv_tol",         "conv_tol",         "SOUPC_CONV_TOL",         f32,   0.01);
    let max_iter         = resolve_parse!("max_iter",         "max_iter",         "SOUPC_MAX_ITER",         usize, 1000);
    let min_cluster_cells= resolve_parse!("min_cluster_cells","min_cluster_cells","SOUPC_MIN_CLUSTER_CELLS",usize, 200);

    // ── M-step pseudocounts ───────────────────────────────────────────────────
    let pseudocount_alt  = resolve_parse!("pseudocount_alt",  "pseudocount_alt",  "SOUPC_PSEUDOCOUNT_ALT",  f32,   1.0);
    let pseudocount_total= resolve_parse!("pseudocount_total","pseudocount_total","SOUPC_PSEUDOCOUNT_TOTAL", f32,   2.0);

    // ── Allele-frequency bounds ───────────────────────────────────────────────
    let theta_min = resolve_parse!("theta_min", "theta_min", "SOUPC_THETA_MIN", f32, 0.01);
    let theta_max = resolve_parse!("theta_max", "theta_max", "SOUPC_THETA_MAX", f32, 0.99);
    assert!(theta_min < theta_max, "theta_min ({}) must be < theta_max ({})", theta_min, theta_max);

    // ── KHM-specific ──────────────────────────────────────────────────────────
    let khm_p = resolve_parse!("khm_p", "khm_p", "SOUPC_KHM_P", f32, 25.0);

    // ── Logging ───────────────────────────────────────────────────────────────
    let verbose           = resolve_bool_flag!("verbose",           "verbose",           "SOUPC_VERBOSE",           false);
    let progress_interval = resolve_parse!("progress_interval", "progress_interval", "SOUPC_PROGRESS_INTERVAL", usize, 500);

    // ── Optional priors ───────────────────────────────────────────────────────
    let known_cell_assignments = resolve_opt_str!("known_cell_assignments",
                                                  "known_cell_assignments",
                                                  "SOUPC_KNOWN_CELL_ASSIGNMENTS");
    let known_genotypes        = resolve_opt_str!("known_genotypes",
                                                  "known_genotypes",
                                                  "SOUPC_KNOWN_GENOTYPES");
    let known_genotypes_sample_names: Vec<String> = {
        let cli_vals: Vec<String> = matches
            .values_of("known_genotypes_sample_names")
            .map(|vs| vs.map(str::to_string).collect())
            .unwrap_or_default();
        if !cli_vals.is_empty() {
            cli_vals
        } else if let Some(ref c) = cfg {
            let from_cfg = c.str_vec("known_genotypes_sample_names");
            if !from_cfg.is_empty() { from_cfg }
            else {
                env.str("SOUPC_KNOWN_GENOTYPES_SAMPLE_NAMES")
                    .map(|s| s.split(',').map(|x| x.trim().to_string()).collect())
                    .unwrap_or_default()
            }
        } else {
            env.str("SOUPC_KNOWN_GENOTYPES_SAMPLE_NAMES")
                .map(|s| s.split(',').map(|x| x.trim().to_string()).collect())
                .unwrap_or_default()
        }
    };

    // ── Diagnostics ───────────────────────────────────────────────────────────
    let plot_dir = resolve_opt_str!("plot_dir", "plot_dir", "SOUPC_PLOT_DIR")
                       .map(|p| expand_tilde(&p));
    let plots          = resolve_str!("plots", "plots", "SOUPC_PLOTS", "all");
    let cleanup        = resolve_bool_flag!("cleanup", "cleanup", "SOUPC_CLEANUP", false);
    let dry_run        = resolve_bool_flag!("dry_run",     "dry_run",     "SOUPC_DRY_RUN",     false);
    let dry_run_yes    = resolve_bool_flag!("dry_run_yes", "dry_run_yes", "SOUPC_DRY_RUN_YES", false);
    let config_profile = config_path.map(str::to_string);

    log_param_sources(config_path, env_file.as_deref(), &ref_mtx, &alt_mtx,
                      &barcodes, num_clusters, restarts, seed, threads,
                      anneal_steps, anneal_base_em, anneal_base_khm,
                      conv_tol, max_iter, pseudocount_alt, pseudocount_total,
                      theta_min, theta_max, khm_p, min_cluster_cells);

    Params {
        ref_mtx, alt_mtx, barcodes,
        num_clusters, restarts, seed,
        clustering_method, initialization_strategy, souporcell3,
        min_alt, min_ref, min_alt_umis, min_ref_umis,
        threads,
        known_cell_assignments, known_genotypes, known_genotypes_sample_names,
        verbose, progress_interval,
        plot_dir, plots, cleanup,
        anneal_steps, anneal_base_em, anneal_base_khm,
        conv_tol, max_iter, min_cluster_cells,
        pseudocount_alt, pseudocount_total,
        theta_min, theta_max,
        khm_p,
        config_profile, env_file,
        dry_run, dry_run_yes,
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn find_env_path() -> Option<String> {
    let cwd = std::env::current_dir().unwrap_or_else(|_| std::path::PathBuf::from("."));
    let p = cwd.join(".souporcell.env");
    if p.exists() { return Some(p.to_string_lossy().to_string()); }
    if let Ok(home) = std::env::var("HOME") {
        let hp = std::path::PathBuf::from(home).join(".souporcell.env");
        if hp.exists() { return Some(hp.to_string_lossy().to_string()); }
    }
    None
}

/// Expand a leading ~ to $HOME so paths work from .env, config JSON, and CLI.
fn expand_tilde(path: &str) -> String {
    if path == "~" {
        return std::env::var("HOME").unwrap_or_else(|_| "~".to_string());
    }
    if path.starts_with("~/") || path.starts_with("~\\") {
        if let Ok(home) = std::env::var("HOME") {
            return format!("{}/{}", home.trim_end_matches('/'), &path[2..]);
        }
    }
    path.to_string()
}

#[allow(clippy::too_many_arguments)]
fn log_param_sources(
    config_path:       Option<&str>,
    env_file:          Option<&str>,
    ref_mtx:           &str,
    alt_mtx:           &str,
    barcodes:          &str,
    k:                 usize,
    restarts:          u32,
    seed:              u64,
    threads:           usize,
    anneal_steps:      usize,
    anneal_base_em:    f32,
    anneal_base_khm:   f32,
    conv_tol:          f32,
    max_iter:          usize,
    pseudocount_alt:   f32,
    pseudocount_total: f32,
    theta_min:         f32,
    theta_max:         f32,
    khm_p:             f32,
    min_cluster_cells: usize,
) {
    eprintln!("[PARAMS] ──────────────────────────────────────────────────────────");
    match (config_path, env_file) {
        (Some(c), Some(e)) => eprintln!("[PARAMS] Sources         : CLI > config({}) > env({})", c, e),
        (Some(c), None)    => eprintln!("[PARAMS] Sources         : CLI > config({})", c),
        (None,    Some(e)) => eprintln!("[PARAMS] Sources         : CLI > env({})", e),
        (None,    None)    => eprintln!("[PARAMS] Sources         : CLI > built-in defaults"),
    }
    eprintln!("[PARAMS] ref_matrix      : {}", ref_mtx);
    eprintln!("[PARAMS] alt_matrix      : {}", alt_mtx);
    eprintln!("[PARAMS] barcodes        : {}", barcodes);
    eprintln!("[PARAMS] k={} restarts={} seed={} threads={}", k, restarts, seed, threads);
    eprintln!("[PARAMS] ── Annealing ──────────────────────────────────────────────");
    eprintln!("[PARAMS] anneal_steps    : {}   base_em: {}   base_khm: {}",
              anneal_steps, anneal_base_em, anneal_base_khm);
    eprintln!("[PARAMS] ── Convergence ─────────────────────────────────────────────");
    eprintln!("[PARAMS] conv_tol        : {}   max_iter: {}   min_cluster_cells: {}",
              conv_tol, max_iter, min_cluster_cells);
    eprintln!("[PARAMS] ── M-step / Theta ─────────────────────────────────────────");
    eprintln!("[PARAMS] pseudocount     : alt={} total={}   theta: [{}, {}]",
              pseudocount_alt, pseudocount_total, theta_min, theta_max);
    eprintln!("[PARAMS] khm_p           : {}", khm_p);
    eprintln!("[PARAMS] ──────────────────────────────────────────────────────────");
}
