// ============================================================================
// core/mod.rs — Algorithm core layer
// ============================================================================
//
// Pure clustering algorithms.  These modules receive data by reference and
// return results — they perform no I/O themselves and do not log directly
// (structured log calls go through infra::logger).
//
// Dependency rule: core may import domain and config.  It must NOT import
// infra (no file I/O inside algorithms), and must NOT import analysis (no
// plot generation inside the hot path).
//
// Modules:
//   em            — Expectation-Maximisation with deterministic annealing
//   khm           — K-Harmonic Means algorithm
//   cluster_init  — All cluster-centre initialisation strategies
// ============================================================================

pub mod em;
pub mod khm;
pub mod cluster_init;

// Re-export primary entry points.
#[allow(unused_imports)]
pub use em::em;
#[allow(unused_imports)]
pub use khm::khm;
#[allow(unused_imports)]
pub use cluster_init::init_cluster_centers;
