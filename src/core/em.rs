// ============================================================================
// em.rs — Expectation-Maximisation clustering with deterministic annealing
// ============================================================================
//
// All formerly hardcoded constants are now read from Params:
//   params.anneal_steps      (was: 9)
//   params.anneal_base_em    (was: 20.0)
//   params.conv_tol          (was: 0.01)
//   params.max_iter          (was: 1000)
//   params.min_cluster_cells (was: 200)
//   params.pseudocount_alt   (was: 1.0)
//   params.pseudocount_total (was: 2.0)
//   params.theta_min         (was: 0.01)
//   params.theta_max         (was: 0.99)
//
// Algorithm:
//   E-step: for each cell, compute posterior P(cluster | allele counts, θ)
//           with temperature-scaled softmax to prevent early convergence
//   M-step: update θ_c = (Σ prob_i × alt_i + pseudocount_alt)
//                       / (Σ prob_i × total_i + pseudocount_total)
//           clamped to [theta_min, theta_max]
//
// Convergence: |Δ total log-loss| < conv_tol × n_cells  OR  max_iter iters
// ============================================================================

use crate::infra::logger;
use crate::domain::math::{binomial_loss_with_min_index, log_sum_exp, normalize_in_log_with_temp};
use crate::config::params::Params;
use crate::domain::types::{CellData, ConvergenceStats};

// ── em ────────────────────────────────────────────────────────────────────────

pub fn em(
    loci:                   usize,
    cluster_centers:        &mut Vec<Vec<f32>>,
    cell_data:              &[CellData],
    params:                 &Params,
    epoch:                  usize,
    thread_num:             usize,
    locked_cluster_centers: Vec<usize>,
) -> (f32, Vec<Vec<f32>>, ConvergenceStats, Vec<(usize,usize,usize,usize,f32)>) {

    let num_clusters  = params.num_clusters;
    let n_cells       = cell_data.len();
    let log_prior:f32 = (1.0 / num_clusters as f32).ln();

    // ── All tunable constants come from Params ───────────────────────────────
    let temp_steps         = params.anneal_steps;
    let anneal_base        = params.anneal_base_em;
    let log_loss_chg_limit = params.conv_tol * n_cells as f32;
    let max_iterations     = params.max_iter;
    let min_cells_pop      = params.min_cluster_cells;

    let mut ts_losses: Vec<(usize,usize,usize,usize,f32)> = Vec::new();
    let mut sums   = alloc_sums(loci, num_clusters, params.pseudocount_alt);
    let mut denoms = alloc_denoms(loci, num_clusters, params.pseudocount_total);

    let mut total_log_loss = f32::NEG_INFINITY;
    let mut final_log_probabilities: Vec<Vec<f32>> = vec![Vec::new(); n_cells];

    let mut last_log_loss    = f32::NEG_INFINITY;
    let mut total_iterations = 0usize;

    for temp_step in 0..temp_steps {
        let mut log_loss_change = 10_000.0f32;
        let mut iter_in_step    = 0usize;

        while log_loss_change > log_loss_chg_limit && total_iterations < max_iterations {
            let mut log_binom_loss = 0.0f32;
            reset_sums_denoms(loci, &mut sums, &mut denoms, num_clusters,
                              params.pseudocount_alt, params.pseudocount_total);

            for (celldex, cell) in cell_data.iter().enumerate() {
                let (log_binoms, _) =
                    binomial_loss_with_min_index(cell, cluster_centers, log_prior);
                log_binom_loss += log_sum_exp(&log_binoms);

                // Temperature: high at step 0 (soft) → 1.0 at last step (hard)
                let raw_temp = cell.total_alleles / (anneal_base * 2.0f32.powf(temp_step as f32));
                let temp     = if temp_step == temp_steps - 1 { 1.0 } else { raw_temp.max(1.0) };

                let probs = normalize_in_log_with_temp(&log_binoms, temp);
                update_centers_average(&mut sums, &mut denoms, cell, &probs);
                final_log_probabilities[celldex] = log_binoms;
            }

            total_log_loss  = log_binom_loss;
            log_loss_change = log_binom_loss - last_log_loss;
            last_log_loss   = log_binom_loss;
            update_final_with_lock(loci, &sums, &denoms, cluster_centers,
                                   &locked_cluster_centers,
                                   params.theta_min, params.theta_max);

            total_iterations += 1;
            iter_in_step     += 1;

            let populated = count_populated_clusters(
                &final_log_probabilities, num_clusters, min_cells_pop);

            if params.verbose {
                logger::log_convergence(
                    thread_num, epoch, total_iterations, temp_step,
                    log_binom_loss, log_loss_change, populated, "em",
                );
            }
        }

        eprintln!(
            "EM\tthread={}\tepoch={}\ttemp_step={}/{}\titers={}\tloss={:.4}\tdelta={:.4}",
            thread_num, epoch, temp_step, temp_steps - 1,
            iter_in_step, last_log_loss, log_loss_change
        );
        ts_losses.push((thread_num, epoch, temp_step, iter_in_step, last_log_loss));
    }

    let stats = ConvergenceStats {
        total_iterations,
        final_log_loss: total_log_loss,
        temp_steps_used: temp_steps,
    };
    (total_log_loss, final_log_probabilities, stats, ts_losses)
}

// ── Shared centre-update helpers (pub(crate) — also used by khm.rs) ──────────

/// Allocate sums matrix initialised to pseudocount_alt
pub(crate) fn alloc_sums(loci: usize, num_clusters: usize, pseudocount: f32) -> Vec<Vec<f32>> {
    vec![vec![pseudocount; loci]; num_clusters]
}

/// Allocate denoms matrix initialised to pseudocount_total
pub(crate) fn alloc_denoms(loci: usize, num_clusters: usize, pseudocount: f32) -> Vec<Vec<f32>> {
    vec![vec![pseudocount; loci]; num_clusters]
}

/// Reset sums/denoms to pseudocounts before each E-step
pub(crate) fn reset_sums_denoms(
    loci:              usize,
    sums:              &mut Vec<Vec<f32>>,
    denoms:            &mut Vec<Vec<f32>>,
    num_clusters:      usize,
    pseudocount_alt:   f32,
    pseudocount_total: f32,
) {
    for cluster in 0..num_clusters {
        for index in 0..loci {
            sums[cluster][index]   = pseudocount_alt;
            denoms[cluster][index] = pseudocount_total;
        }
    }
}

/// M-step: accumulate weighted alt counts (soft assignment)
pub(crate) fn update_centers_average(
    sums:          &mut Vec<Vec<f32>>,
    denoms:        &mut Vec<Vec<f32>>,
    cell:          &CellData,
    probabilities: &[f32],
) {
    for locus_idx in 0..cell.loci.len() {
        let locus = cell.loci[locus_idx];
        let alt   = cell.alt_counts[locus_idx] as f32;
        let total = alt + cell.ref_counts[locus_idx] as f32;
        for (cluster, &prob) in probabilities.iter().enumerate() {
            sums[cluster][locus]   += prob * alt;
            denoms[cluster][locus] += prob * total;
        }
    }
}

/// Apply updated sums/denoms to cluster_centers, skipping locked clusters.
/// theta clamped to [theta_min, theta_max].
pub(crate) fn update_final_with_lock(
    loci:                   usize,
    sums:                   &Vec<Vec<f32>>,
    denoms:                 &Vec<Vec<f32>>,
    cluster_centers:        &mut Vec<Vec<f32>>,
    locked_cluster_centers: &[usize],
    theta_min:              f32,
    theta_max:              f32,
) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            if locked_cluster_centers.contains(&cluster) { continue; }
            cluster_centers[cluster][locus] =
                (sums[cluster][locus] / denoms[cluster][locus]).clamp(theta_min, theta_max);
        }
    }
}

/// Count clusters with > threshold cells assigned
pub(crate) fn count_populated_clusters(
    final_log_probs: &[Vec<f32>],
    num_clusters:    usize,
    threshold:       usize,
) -> usize {
    let mut cells_per = vec![0usize; num_clusters];
    for lps in final_log_probs {
        if lps.is_empty() { continue; }
        let best = lps.iter().enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(i, _)| i).unwrap_or(0);
        cells_per[best] += 1;
    }
    cells_per.iter().filter(|&&c| c > threshold).count()
}
