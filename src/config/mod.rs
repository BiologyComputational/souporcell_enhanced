// ============================================================================
// config/mod.rs — Configuration layer
// ============================================================================
//
// Modules:
//   params         — Params struct, CLI parsing, 3-layer merge
//   env_loader     — .souporcell.env file reader
//   config_loader  — JSON config profile reader
// ============================================================================

pub mod env_loader;
pub mod config_loader;
pub mod params;

#[allow(unused_imports)]
pub use params::{load_params, Params, ClusterMethod, ClusterInit};
