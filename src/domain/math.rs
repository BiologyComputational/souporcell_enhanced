// ============================================================================
// math.rs — Numerically-stable math primitives used by EM and KHM
// ============================================================================
//
// All functions here are pure (no I/O, no side effects) and operate on
// f32 to match the original codebase's memory / speed trade-off.
//
// Functions:
//   log_sum_exp              — numerically stable log(Σ exp(xᵢ))
//   normalize_in_log_with_temp — softmax with temperature
//   binomial_loss_with_min_index — per-cell log-likelihood across all clusters
//
// The `_` prefixed originals (_sum_of_squares_loss, _normalize_in_log,
// _update_centers_flat) are retained as dead-code stubs so the git diff
// remains clean. They are not pub.
// ============================================================================

use crate::domain::types::CellData;

// ── log_sum_exp ───────────────────────────────────────────────────────────────
//
// Computes log(Σᵢ exp(pᵢ)) in a numerically stable way by factoring out
// the maximum value:  log(exp(max) · Σᵢ exp(pᵢ - max)) = max + log(Σᵢ exp(pᵢ - max))
//
// Returns NEG_INFINITY for an empty slice.
#[inline]
pub fn log_sum_exp(p: &[f32]) -> f32 {
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    if max_p == f32::NEG_INFINITY {
        return f32::NEG_INFINITY;
    }
    let sum: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum.ln()
}

// ── normalize_in_log_with_temp ────────────────────────────────────────────────
//
// Softmax with temperature:  p̃ᵢ = exp((log_pᵢ / temp) - log_Z)
//
// temp = 1.0   → standard softmax (used in final annealing step)
// temp > 1.0   → softer / more uniform (early annealing steps)
// temp → 0     → hard argmax
//
// The original used cell.total_alleles as the temperature divisor, which
// means cells with more reads experience less smoothing (correct behaviour).
pub fn normalize_in_log_with_temp(log_probs: &[f32], temp: f32) -> Vec<f32> {
    let tempered: Vec<f32> = log_probs.iter().map(|lp| lp / temp).collect();
    let z = log_sum_exp(&tempered);
    tempered.iter().map(|lp| (lp - z).exp()).collect()
}

// ── binomial_loss_with_min_index ─────────────────────────────────────────────
//
// For a single cell, computes the log-likelihood under each cluster's
// current allele-fraction estimates using a binomial model:
//
//   log P(cell | cluster c) = log_prior
//       + Σ_locus [ ln C(n,k) + k·ln(θ_c) + (n-k)·ln(1-θ_c) ]
//
// where θ_c = cluster_centers[c][locus], k = alt count, n = total count.
//
// Returns (log_probabilities, index_of_best_cluster).
pub fn binomial_loss_with_min_index(
    cell_data:       &CellData,
    cluster_centers: &[Vec<f32>],
    log_prior:       f32,
) -> (Vec<f32>, usize) {
    let mut log_probabilities = Vec::with_capacity(cluster_centers.len());
    let mut best_log = f32::MIN;
    let mut best_idx = 0;

    for (cluster, center) in cluster_centers.iter().enumerate() {
        let mut lp = log_prior;
        for (locus_idx, &locus) in cell_data.loci.iter().enumerate() {
            let theta     = center[locus];
            let alt_c     = cell_data.alt_counts[locus_idx] as f32;
            let ref_c     = cell_data.ref_counts[locus_idx] as f32;
            let log_binom = cell_data.log_binomial_coefficient[locus_idx];
            lp += log_binom
                + alt_c * theta.ln()
                + ref_c * (1.0 - theta).ln();
        }
        if lp > best_log {
            best_log = lp;
            best_idx = cluster;
        }
        log_probabilities.push(lp);
    }
    (log_probabilities, best_idx)
}

// ── Dead-code stubs (kept for reference / future use) ────────────────────────

#[allow(dead_code)]
fn sum_of_squares_loss(
    cell_data:       &CellData,
    cluster_centers: &[Vec<f32>],
    log_prior:       f32,
) -> Vec<f32> {
    let mut lps = Vec::with_capacity(cluster_centers.len());
    for (cluster, center) in cluster_centers.iter().enumerate() {
        lps.push(log_prior);
        for (locus_idx, &locus) in cell_data.loci.iter().enumerate() {
            lps[cluster] -= (cell_data.allele_fractions[locus_idx] - center[locus]).powf(2.0);
        }
    }
    lps
}

#[allow(dead_code)]
fn normalize_in_log(log_probs: &[f32]) -> Vec<f32> {
    let z = log_sum_exp(log_probs);
    log_probs.iter().map(|lp| (lp - z).exp()).collect()
}
