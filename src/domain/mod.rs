// ============================================================================
// domain/mod.rs — Pure domain layer
// ============================================================================
//
// This layer contains zero I/O, zero CLI dependencies, and zero knowledge of
// the outside world.  It defines the core data structures and pure mathematical
// operations that every other layer builds on.
//
// Dependency rule: domain modules may only depend on each other and the Rust
// standard library.  They must NOT import from config, core, infra, or analysis.
//
// Modules:
//   types  — CellData, ThreadData, ConvergenceStats (shared data structures)
//   math   — log_sum_exp, normalize, binomial_loss (pure f32 functions)
// ============================================================================

pub mod types;
pub mod math;

// Re-export the most commonly used items so callers can write
// `use crate::domain::{CellData, log_sum_exp}` instead of the full path.
#[allow(unused_imports)]
pub use types::{CellData, ThreadData, ConvergenceStats};
#[allow(unused_imports)]
pub use math::{log_sum_exp, normalize_in_log_with_temp, binomial_loss_with_min_index};
