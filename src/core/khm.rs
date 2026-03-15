// ============================================================================
// khm.rs — K-Harmonic Means clustering with deterministic annealing
// ============================================================================
//
// All formerly hardcoded constants are now read from Params:
//   params.anneal_steps      (was: 9)
//   params.anneal_base_khm   (was: 0.5)
//   params.conv_tol          (was: 0.01)
//   params.max_iter          (was: 1000)
//   params.min_cluster_cells (was: 200)
//   params.pseudocount_alt   (was: 1.0)
//   params.pseudocount_total (was: 2.0)
//   params.theta_min         (was: 0.01)
//   params.theta_max         (was: 0.99)
//   params.khm_p             (was: 25.0 const KHM_P)
//
// KHM membership function uses a harmonic-mean-based weighting q(cell, cluster)
// that is more robust to outlier loci than EM's posterior.
// ============================================================================

use crate::core::em::{
    alloc_denoms, alloc_sums, count_populated_clusters, reset_sums_denoms,
    update_centers_average, update_final_with_lock,
};
use crate::infra::logger;
use crate::domain::math::{binomial_loss_with_min_index, log_sum_exp, normalize_in_log_with_temp};
use crate::config::params::Params;
use crate::domain::types::{CellData, ConvergenceStats};

// ── khm ───────────────────────────────────────────────────────────────────────

pub fn khm(
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
    let anneal_base        = params.anneal_base_khm;
    let log_loss_chg_limit = params.conv_tol * n_cells as f32;
    let max_iterations     = params.max_iter;
    let min_cells_pop      = params.min_cluster_cells;
    let khm_p              = params.khm_p;

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
                let (log_binoms, min_clus) =
                    binomial_loss_with_min_index(cell, cluster_centers, log_prior);
                log_binom_loss += log_sum_exp(&log_binoms);

                let (q_vec, q_sum) = calculate_q(&log_binoms, min_clus, khm_p);
                let khm_log_probs: Vec<f32> = q_vec.iter().map(|&q| q - q_sum).collect();

                // Temperature: high at step 0 (soft) → 1.0 at last step (hard)
                let raw_temp = cell.total_alleles / (anneal_base * 2.0f32.powf(temp_step as f32));
                let temp     = if temp_step == temp_steps - 1 { 1.0 } else { raw_temp.max(1.0) };
                let probs    = normalize_in_log_with_temp(&khm_log_probs, temp);

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
                    log_binom_loss, log_loss_change, populated, "khm",
                );
            }
        }

        eprintln!(
            "KHM\tthread={}\tepoch={}\ttemp_step={}/{}\titers={}\tloss={:.4}\tdelta={:.4}",
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

// ── calculate_q ───────────────────────────────────────────────────────────────
//
// KHM membership function. khm_p now comes from params instead of const KHM_P.
//
//   log_q_denom = log_sum_exp([ p*(log_winner + log_loss_i)  for i≠winner,
//                                0.0  for winner ])
//   log_q_i     = (2*p*log_winner) - ((p+2)*(-log_loss_i)) - 2*log_q_denom
//   q_sum       = log_sum_exp(log_q_vec)
fn calculate_q(log_loss_vec: &[f32], min_clus: usize, khm_p: f32) -> (Vec<f32>, f32) {
    let log_winner = -log_loss_vec[min_clus];

    let log_denom_terms: Vec<f32> = log_loss_vec.iter().enumerate().map(|(i, &lp)| {
        if i == min_clus { 0.0 }
        else             { khm_p * (log_winner + lp) }
    }).collect();
    let q_denom = log_sum_exp(&log_denom_terms);

    let q_vec: Vec<f32> = log_loss_vec.iter().map(|&lp| {
        (2.0 * khm_p * log_winner) - ((khm_p + 2.0) * (-lp)) - (2.0 * q_denom)
    }).collect();

    let q_sum = log_sum_exp(&q_vec);
    (q_vec, q_sum)
}
