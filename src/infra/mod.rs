// ============================================================================
// infra/mod.rs — Infrastructure layer
// ============================================================================
//
// Everything that touches the outside world: filesystem and terminal.
// Modules here are the only ones allowed to perform I/O.
//
// Dependency rule: infra may import domain and config.  It must NOT import
// core or analysis — infrastructure serves all layers but depends on none
// of them except domain types and config.
//
// Modules:
//   io          — Matrix loading (.mtx), barcode loading (.tsv/.gz), result writing
//   logger      — Structured stderr logging, ANSI colour, section timing
//   preflight   — Pre-flight run plan printer and user approval gate (--dry_run)
//   gui_server  — Zero-dependency HTTP server for --gui mode
// ============================================================================

pub mod io;
pub mod logger;
pub mod preflight;
pub mod gui_server;

// Re-export the most-used items.
#[allow(unused_imports)]
pub use io::{LoadResult, load_barcodes, load_cell_data, write_cluster_assignments, reader};
#[allow(unused_imports)]
pub use logger::{Logger, log_convergence, log_restart, log_loading_stats,
                 log_cluster_analysis, log_run_summary};
#[allow(unused_imports)]
pub use preflight::{show_preflight, ApprovalResult};
#[allow(unused_imports)]
pub use gui_server::launch as launch_gui;
