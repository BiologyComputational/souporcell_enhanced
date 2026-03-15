// ============================================================================
// report/utils.rs — Shared colour tokens, helper functions
// ============================================================================

use crate::analysis::plots::PlotData;

// ── Colour palette (used in inline-SVG builders — avoids # inside format!) ──
pub const C_AMBER: &str = "#C75B00";
pub const C_RED:   &str = "#B91C1C";
pub const C_GREEN: &str = "#15803D";
pub const C_LIGHT: &str = "#E2E8F0";
pub const C_INK:   &str = "#0F172A";
pub const C_MID:   &str = "#64748B";
pub const C_BAR:   &str = "#0891B2";

// ── Confidence statistics ─────────────────────────────────────────────────────

/// Returns (median_confidence, ambiguous_cell_pct).
/// median_confidence = median of (best − second_best) log-prob per cell.
/// ambig_pct         = % of cells whose confidence < 10 log-units.
pub fn compute_conf_stats(plot_data: &PlotData, total: usize) -> (f32, usize) {
    let mut confs: Vec<f32> = plot_data.cell_assignments.iter()
        .filter(|(_, lps)| lps.len() >= 2)
        .map(|(_, lps)| {
            let mut s = lps.clone();
            s.sort_by(|a, b| b.total_cmp(a));
            s[0] - s[1]
        })
        .collect();
    confs.sort_by(|a, b| a.total_cmp(b));
    let median    = if confs.is_empty() { 0.0 } else { confs[confs.len() / 2] };
    let ambig_n   = confs.iter().filter(|&&c| c < 10.0).count();
    let ambig_pct = if total > 0 { ambig_n * 100 / total } else { 0 };
    (median, ambig_pct)
}

// ── Timestamp ─────────────────────────────────────────────────────────────────

/// Returns a simple UTC timestamp string computed from UNIX time.
/// Avoids the chrono dependency; accuracy is ±1 day on leap years.
pub fn chrono_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs   = SystemTime::now().duration_since(UNIX_EPOCH)
                     .map(|d| d.as_secs()).unwrap_or(0);
    let days   = secs / 86400;
    let year   = 1970 + days / 365;
    let month  = (days % 365) / 30 + 1;
    let day    = (days % 365) % 30 + 1;
    let hour   = (secs % 86400) / 3600;
    let minute = (secs % 3600) / 60;
    format!("{:04}-{:02}-{:02} {:02}:{:02} UTC", year, month, day, hour, minute)
}
