// ============================================================================
// types.rs — Core data structures shared across modules
// ============================================================================

use rand::rngs::StdRng;
use rand::SeedableRng;

// ── CellData ──────────────────────────────────────────────────────────────────
//
// Holds all allele-level information for a single cell barcode.
// Populated by io::load_cell_data() and consumed by em / khm / math modules.
//
// Memory layout note: all Vec fields are co-indexed — element i of
// `loci`, `alt_counts`, `ref_counts`, `allele_fractions`, and
// `log_binomial_coefficient` all refer to the same locus observation.
pub struct CellData {
    /// Index into the compressed locus table (0..loci_used)
    pub loci: Vec<usize>,
    /// Alt allele UMI count at each locus
    pub alt_counts: Vec<u32>,
    /// Ref allele UMI count at each locus
    pub ref_counts: Vec<u32>,
    /// alt / (alt + ref) at each locus — precomputed for EM speed
    pub allele_fractions: Vec<f32>,
    /// ln C(n, k) precomputed once — used in binomial log-likelihood
    pub log_binomial_coefficient: Vec<f32>,
    /// Sum of (alt + ref) across all loci — used for temperature annealing
    pub total_alleles: f32,
}

impl CellData {
    pub fn new() -> CellData {
        CellData {
            loci:                    Vec::new(),
            alt_counts:              Vec::new(),
            ref_counts:              Vec::new(),
            allele_fractions:        Vec::new(),
            log_binomial_coefficient: Vec::new(),
            total_alleles:           0.0,
        }
    }

    /// Number of covered loci for this cell
    #[allow(dead_code)]
    pub fn num_loci(&self) -> usize {
        self.loci.len()
    }

    /// Total read depth across all loci (same as total_alleles but usize)
    #[allow(dead_code)]
    pub fn total_depth(&self) -> u32 {
        self.alt_counts.iter().zip(self.ref_counts.iter())
            .map(|(a, r)| a + r)
            .sum()
    }
}

// ── ThreadData ────────────────────────────────────────────────────────────────
//
// Per-thread state for the parallel restart loop in souporcell_main().
// Each thread owns its RNG independently so results are deterministic
// given a fixed global seed regardless of thread scheduling.
pub struct ThreadData {
    /// Best per-cell log-probabilities (shape: cells × clusters) from all
    /// restarts this thread has completed so far
    pub best_log_probabilities:    Vec<Vec<f32>>,
    /// Total log-loss of the best solution found so far
    pub best_total_log_probability: f32,
    /// Independent RNG seeded from the global seed + thread index
    pub rng:                       StdRng,
    /// How many EM/KHM restarts this thread is responsible for
    pub solves_per_thread:         usize,
    /// Zero-based thread index — used in CONV log lines
    pub thread_num:                usize,
    /// Cluster centre vectors corresponding to best_log_probabilities
    pub cluster_centers:           Vec<Vec<f32>>,
    /// Count of restarts actually completed (may be less than solves_per_thread
    /// if early termination is added in future)
    pub restarts_completed:        usize,
}

impl ThreadData {
    pub fn from_seed(seed: [u8; 32], solves_per_thread: usize, thread_num: usize) -> ThreadData {
        ThreadData {
            best_log_probabilities:     Vec::new(),
            best_total_log_probability: f32::NEG_INFINITY,
            rng:                        SeedableRng::from_seed(seed),
            solves_per_thread,
            thread_num,
            cluster_centers:            Vec::new(),
            restarts_completed:         0,
        }
    }
}

// ── ConvergenceStats ─────────────────────────────────────────────────────────
//
// Returned by em() and khm() so callers (main / tests) can inspect
// convergence behaviour without parsing stderr.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConvergenceStats {
    pub total_iterations: usize,
    pub final_log_loss:   f32,
    pub temp_steps_used:  usize,
}
