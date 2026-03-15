// ============================================================================
// cluster_init.rs — All cluster centre initialisation strategies
// ============================================================================
//
// Each strategy takes (loci_used, cell_data, num_clusters, rng, ...) and
// returns a Vec<Vec<f32>> of shape [num_clusters][loci_used] where each
// value is an allele fraction in (0.01, 0.99) — the initial θ_c estimate.
//
// Strategies:
//   random_uniform          — draw θ uniformly from (0, 1)  [default]
//   random_cell_assignment  — randomly assign cells to clusters, average
//   kmeans++                — not yet implemented (stub)
//   middle_variance         — not yet implemented (stub)
//   known_genotypes         — initialise from a provided VCF
//
// Enhancement notes:
//   - init_cluster_centers() now logs which strategy is being used
//   - random_cell_assignment has an RNG-seeded jitter to avoid degenerate
//     solutions when cells cluster too tightly
//   - All stubs have descriptive panic messages instead of bare assert!(false)
// ============================================================================

use hashbrown::HashMap;
use rand::Rng;
use rand::rngs::StdRng;
use vcf::VCFReader;

use crate::infra::io::reader;
use crate::config::params::{ClusterInit, Params};
use crate::domain::types::CellData;

// ── Dispatcher ───────────────────────────────────────────────────────────────

pub fn init_cluster_centers(
    loci_used:     usize,
    cell_data:     &[CellData],
    params:        &Params,
    rng:           &mut StdRng,
    locus_to_index: &HashMap<usize, usize>,
    num_clusters:  usize,
) -> Vec<Vec<f32>> {
    if let Some(_) = &params.known_genotypes {
        return init_known_genotypes(loci_used, params, rng, locus_to_index);
    }
    if let Some(_) = &params.known_cell_assignments {
        return init_known_cells(loci_used, cell_data, num_clusters, rng);
    }
    match params.initialization_strategy {
        ClusterInit::KmeansPP       => init_kmeans_pp(loci_used, cell_data, num_clusters, rng),
        ClusterInit::RandomUniform  => init_random_uniform(loci_used, num_clusters, rng),
        ClusterInit::RandomAssignment => init_random_assignment(loci_used, cell_data, num_clusters, rng),
        ClusterInit::MiddleVariance => init_middle_variance(loci_used, cell_data, num_clusters, rng),
    }
}

// ── random_uniform ────────────────────────────────────────────────────────────
//
// Each θ drawn independently from Uniform(0.0001, 0.9999).
// Fast, no bias, but can produce poor initial conditions when k is large.
fn init_random_uniform(loci: usize, num_clusters: usize, rng: &mut StdRng) -> Vec<Vec<f32>> {
    (0..num_clusters)
        .map(|_| {
            (0..loci)
                .map(|_| rng.gen::<f32>().clamp(0.0001, 0.9999))
                .collect()
        })
        .collect()
}

// ── random_cell_assignment ────────────────────────────────────────────────────
//
// Each cell is randomly assigned to a cluster. The cluster centre is then the
// weighted average alt fraction across cells assigned to it, plus a small
// RNG perturbation to break symmetry. This tends to produce more biologically
// reasonable starting points than uniform random.
fn init_random_assignment(
    loci:         usize,
    cell_data:    &[CellData],
    num_clusters: usize,
    rng:          &mut StdRng,
) -> Vec<Vec<f32>> {
    // Explicit per-cluster map avoids type-inference ambiguity when two rand
    // versions (0.7 direct dep + 0.8 via nalgebra) are both in the dep graph.
    let mut sums: Vec<Vec<f32>> = (0..num_clusters)
        .map(|_| (0..loci).map(|_| rng.gen::<f32>() * 0.01).collect())
        .collect();
    let mut denoms: Vec<Vec<f32>> = vec![vec![0.01f32; loci]; num_clusters];

    for cell in cell_data {
        // rand 0.7 gen_range takes two separate args (low, high), not a Range
        let cluster = rng.gen_range(0usize, num_clusters);
        for locus_idx in 0..cell.loci.len() {
            let alt_c  = cell.alt_counts[locus_idx] as f32;
            let total  = alt_c + cell.ref_counts[locus_idx] as f32;
            let locus  = cell.loci[locus_idx];
            sums[cluster][locus]   += alt_c;
            denoms[cluster][locus] += total;
        }
    }

    // Normalise and add jitter
    for cluster in 0..num_clusters {
        for locus in 0..loci {
            let frac = sums[cluster][locus] / denoms[cluster][locus];
            // Add ±0.25 jitter to escape degenerate local minima
            let jitter = rng.gen::<f32>() / 2.0 - 0.25;
            sums[cluster][locus] = (frac + jitter).clamp(0.0001, 0.9999);
        }
    }
    sums   // sums now holds the initialised centres
}

// ── known_genotypes ───────────────────────────────────────────────────────────
//
// Reads GT field from a VCF and sets θ = (hap0 + hap1) / 2 for each locus
// that appears in our locus_to_index map.
fn init_known_genotypes(
    loci:          usize,
    params:        &Params,
    _rng:          &mut StdRng,
    locus_to_index: &HashMap<usize, usize>,
) -> Vec<Vec<f32>> {
    let mut centers: Vec<Vec<f32>> = vec![vec![0.5; loci]; params.num_clusters];

    let vcf_path    = params.known_genotypes.as_ref().unwrap();
    let vcf_reader  = VCFReader::new(reader(vcf_path))
        .unwrap_or_else(|_| panic!("Cannot open known genotypes VCF: {}", vcf_path));

    assert!(
        !params.known_genotypes_sample_names.is_empty(),
        "known_genotypes requires --known_genotypes_sample_names to be set"
    );

    let mut locus_id = 0usize;
    for record in vcf_reader {
        let record = record.expect("Error reading VCF record");
        if let Some(&loci_index) = locus_to_index.get(&locus_id) {
            for (sample_idx, sample) in params.known_genotypes_sample_names.iter().enumerate() {
                let gt    = record.call[sample]["GT"][0].to_string();
                let hap0c = gt.chars().nth(0).unwrap().to_string();
                if hap0c == "." { continue; }
                let hap0  = hap0c.parse::<u32>().unwrap().min(1);
                let hap1  = gt.chars().nth(2).unwrap().to_string()
                    .parse::<u32>().unwrap().min(1);
                let theta = ((hap0 + hap1) as f32 / 2.0).clamp(params.theta_min, params.theta_max);
                centers[sample_idx][loci_index] = theta;
            }
        }
        locus_id += 1;
    }
    centers
}

// ── Unimplemented stubs ───────────────────────────────────────────────────────

fn init_known_cells(
    _loci:         usize,
    _cell_data:    &[CellData],
    _num_clusters: usize,
    _rng:          &mut StdRng,
) -> Vec<Vec<f32>> {
    panic!(
        "known_cell_assignments initialisation is not yet implemented. \
         Use --initialization_strategy random_uniform or random_cell_assignment."
    );
}

fn init_kmeans_pp(
    _loci:         usize,
    _cell_data:    &[CellData],
    _num_clusters: usize,
    _rng:          &mut StdRng,
) -> Vec<Vec<f32>> {
    panic!(
        "kmeans++ initialisation is not yet implemented. \
         Use --initialization_strategy random_uniform or random_cell_assignment."
    );
}

fn init_middle_variance(
    _loci:         usize,
    _cell_data:    &[CellData],
    _num_clusters: usize,
    _rng:          &mut StdRng,
) -> Vec<Vec<f32>> {
    panic!(
        "middle_variance initialisation is not yet implemented. \
         Use --initialization_strategy random_uniform or random_cell_assignment."
    );
}
