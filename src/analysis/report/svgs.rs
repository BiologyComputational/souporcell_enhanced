// ============================================================================
// report/svgs.rs — Inline SVG chart builders (§6)
//
// build_balance_svg()    — horizontal bar chart of cells per cluster
// build_confidence_svg() — histogram of per-cell posterior confidence
// ============================================================================

use crate::analysis::report::ReportData;
use crate::analysis::plots::PlotData;
use crate::analysis::report::utils::{C_AMBER, C_RED, C_GREEN, C_LIGHT, C_INK, C_MID, C_BAR};

// ── Cluster balance chart ─────────────────────────────────────────────────────

pub fn build_balance_svg(data: &ReportData) -> String {
    let k     = data.params.num_clusters;
    let cells = &data.cells_per_cluster;
    let total: usize = cells.iter().sum();
    if total == 0 {
        return "<svg width='420' height='220'></svg>".to_string();
    }
    let expected  = total / k.max(1);
    let max_count = cells.iter().copied().max().unwrap_or(1).max(expected);

    let (w, h, pl, pr, pt, pb) = (420usize, 220usize, 70usize, 24usize, 20usize, 36usize);
    let bar_h = (h - pt - pb) / k.max(1);
    let pw    = w - pl - pr;

    let mut svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" style="background:#fff">"##
    );

    // Percentage grid lines
    for pct in [25, 50, 75, 100usize] {
        let x = pl + (pct as f32 / 100.0 * pw as f32) as usize;
        svg += &format!(
            r##"<line x1="{x}" y1="{pt}" x2="{x}" y2="{}" stroke="{C_LIGHT}" stroke-width="1"/>
<text x="{x}" y="{}" font-size="9" fill="{C_MID}" text-anchor="middle">{pct}%</text>"##,
            h - pb, h - pb + 12
        );
    }

    // Ideal equal-split marker
    let xe = pl + (expected as f32 / max_count as f32 * pw as f32) as usize;
    svg += &format!(
        r##"<line x1="{xe}" y1="{pt}" x2="{xe}" y2="{}" stroke="{C_AMBER}" stroke-width="1.5" stroke-dasharray="5,3"/>"##,
        h - pb
    );

    // Bars
    for (i, &count) in cells.iter().enumerate() {
        let bw  = (count as f32 / max_count as f32 * pw as f32) as usize;
        let y   = pt + i * bar_h;
        let pct = count * 100 / total.max(1);
        let bh  = bar_h.saturating_sub(4);
        svg += &format!(
            r##"<rect x="{pl}" y="{y}" width="{bw}" height="{bh}" fill="{C_BAR}" opacity=".8" rx="3"/>
<text x="{}" y="{}" font-size="11" fill="{C_INK}" dominant-baseline="middle">{count} ({pct}%)</text>
<text x="{}" y="{}" font-size="10" fill="{C_MID}" text-anchor="end" dominant-baseline="middle">Cl.{i}</text>"##,
            pl + bw + 7, y + bh / 2,
            pl - 5,      y + bh / 2,
        );
    }

    svg += &format!(
        r##"<text x="{xe}" y="{}" font-size="9" fill="{C_AMBER}" text-anchor="middle">ideal ({expected})</text>"##,
        h - pb + 24
    );
    svg += "</svg>";
    svg
}

// ── Posterior confidence histogram ────────────────────────────────────────────

pub fn build_confidence_svg(plot_data: &PlotData) -> String {
    let mut confs: Vec<f32> = plot_data.cell_assignments.iter()
        .filter(|(_, lps)| lps.len() >= 2)
        .map(|(_, lps)| {
            let mut s = lps.clone();
            s.sort_by(|a, b| b.total_cmp(a));
            (s[0] - s[1]).min(200.0)
        })
        .collect();

    if confs.is_empty() {
        return "<svg width='420' height='220'>\
                <text x='210' y='110' fill='#94A3B8' text-anchor='middle' font-size='13'>\
                No confidence data</text></svg>".to_string();
    }

    confs.sort_by(|a, b| a.total_cmp(b));
    let max_c    = confs.last().copied().unwrap_or(1.0).max(1.0);
    let n_bins   = 40usize;
    let mut bins = vec![0usize; n_bins];
    for &c in &confs {
        let b = ((c / max_c * (n_bins - 1) as f32) as usize).min(n_bins - 1);
        bins[b] += 1;
    }
    let max_bin = *bins.iter().max().unwrap_or(&1).max(&1);

    let (w, h, pl, pr, pt, pb) = (420usize, 220usize, 12usize, 12usize, 14usize, 30usize);
    let pw = w - pl - pr;
    let ph = h - pt - pb;
    let bw = pw / n_bins;

    let ambig_x = pl + (10.0 / max_c * pw as f32) as usize;
    let median  = confs[confs.len() / 2];
    let med_x   = pl + (median / max_c * pw as f32) as usize;

    let mut svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" style="background:#fff">"##
    );

    // Ambiguous zone shading
    svg += &format!(
        r##"<rect x="{pl}" y="{pt}" width="{ambig_x}" height="{ph}" fill="#FEE2E2" rx="2"/>"##
    );

    // Histogram bars
    for (i, &cnt) in bins.iter().enumerate() {
        let bh  = (cnt as f32 / max_bin as f32 * ph as f32) as usize;
        let x   = pl + i * bw;
        let y   = pt + ph - bh;
        let col = if x < ambig_x { C_RED } else { C_BAR };
        svg += &format!(
            r##"<rect x="{x}" y="{y}" width="{}" height="{bh}" fill="{col}" opacity=".75" rx="2"/>"##,
            bw.saturating_sub(1)
        );
    }

    // Median line
    svg += &format!(
        r##"<line x1="{med_x}" y1="{pt}" x2="{med_x}" y2="{}" stroke="{C_GREEN}" stroke-width="2" stroke-dasharray="5,3"/>"##,
        h - pb
    );

    // Labels
    svg += &format!(r##"<text x="{pl}" y="{}" font-size="9" fill="{C_MID}">0</text>"##, h - pb + 12);
    svg += &format!(
        r##"<text x="{}" y="{}" font-size="9" fill="{C_MID}" text-anchor="end">{:.0}</text>"##,
        w - pr, h - pb + 12, max_c
    );
    svg += &format!(
        r##"<text x="{med_x}" y="{}" font-size="9" fill="{C_GREEN}" text-anchor="middle">median {:.0}</text>"##,
        h - pb + 12, median
    );
    svg += &format!(
        r##"<text x="{}" y="{}" font-size="9" fill="{C_RED}">ambiguous</text>"##,
        pl + 3, pt + 10
    );
    svg += "</svg>";
    svg
}
