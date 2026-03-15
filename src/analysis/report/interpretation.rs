// ============================================================================
// report/interpretation.rs — PhD-level per-metric commentary (§8)
// ============================================================================

use crate::analysis::report::ReportData;
use crate::analysis::plots::PlotData;

pub fn build_interpretation_guide(
    data:        &ReportData,
    _plot_data:  &PlotData,
    conv_pct:    usize,
    max_imb:     usize,
    median_conf: f32,
    ambig_pct:   usize,
) -> String {
    let loci_pct = data.loci_used as f32 / data.loci_raw.max(1) as f32 * 100.0;
    let k        = data.params.num_clusters;

    struct Card {
        title:   &'static str,
        cls:     &'static str,
        verdict: String,
        desc:    &'static str,
    }

    let cards = vec![
        Card {
            title: "Convergence stability",
            cls:   if conv_pct >= 70 { "good" } else { "warn" },
            verdict: if conv_pct >= 70 {
                format!(
                    "Strong: {}% of restarts converged to the global optimum. \
                     Posterior landscape is well-defined for k={}.",
                    conv_pct, k
                )
            } else {
                format!(
                    "Moderate at {}%. Sparse SNP signal ({:.1}% of loci passed QC) \
                     creates shallow local optima. Recommend --restarts 200\u{2013}500.",
                    conv_pct, loci_pct
                )
            },
            desc: "Fraction of independent restarts converging to the same \
                   global-optimum log-probability (\u{00b1}1 log-unit). \
                   \u{2265}70% confirms a well-identifiable solution. \
                   <40% suggests k may be too large for the available SNP evidence.",
        },
        Card {
            title: "Cluster balance",
            cls:   if max_imb < 20 { "good" } else { "warn" },
            verdict: if max_imb < 20 {
                format!(
                    "Excellent: max deviation {}% \u{2014} clusters near equal size, \
                     consistent with uniform donor pooling.",
                    max_imb
                )
            } else {
                format!(
                    "Imbalance {}%. Possible causes: unequal donor contributions, \
                     one donor with lower RNA quality, or k > actual donors.",
                    max_imb
                )
            },
            desc: "Maximum fractional deviation of any cluster cell count from the \
                   expected equal-split. Severe imbalance (>50%) can indicate k is \
                   larger than the true number of donors or that a donor\u{2019}s cells \
                   were lost in library preparation.",
        },
        Card {
            title: "Posterior confidence",
            cls:   if ambig_pct < 5 { "good" } else { "warn" },
            verdict: format!(
                "{}% ambiguous cells (gap < 10 log-units). Median = {:.0} log-units. {}",
                ambig_pct, median_conf,
                if ambig_pct < 5 {
                    "Excellent assignment certainty."
                } else {
                    "High ambiguity may indicate doublets or low-coverage cells."
                }
            ),
            desc: "Per-cell confidence = best_log_prob \u{2212} second_best_log_prob. \
                   Values near zero indicate allele patterns consistent with multiple \
                   donors (likely doublets). The median is the single most informative \
                   quality indicator.",
        },
        Card {
            title: "Matrix sparsity",
            cls:   if data.zero_cov_frac > 0.90 { "good" } else { "info" },
            verdict: format!(
                "{:.1}% sparsity \u{2014} {} for scRNA-seq allele data.",
                data.zero_cov_frac * 100.0,
                if data.zero_cov_frac > 0.90 { "normal" } else { "unusually dense" }
            ),
            desc: "Fraction of locus\u{00d7}cell pairs with zero read coverage. \
                   scRNA-seq allele data is inherently sparse (typically 92\u{2013}98%). \
                   Unusually low sparsity may indicate a genome reference mismatch or \
                   overly permissive locus QC filters.",
        },
        Card {
            title: "Locus QC pass rate",
            cls:   if loci_pct > 5.0 { "good" } else { "warn" },
            verdict: format!(
                "{:.1}% of loci pass QC ({} of {} raw). {}",
                loci_pct, data.loci_used, data.loci_raw,
                if loci_pct > 5.0 {
                    "Sufficient informative sites."
                } else {
                    "Very aggressive filtering \u{2014} consider relaxing min_ref/min_alt."
                }
            ),
            desc: "Fraction of raw SNP loci retaining \u{2265} min_ref ref reads and \
                   \u{2265} min_alt alt reads. Too few loci reduce statistical power; \
                   too many add noise. 5\u{2013}15% is typical for stringent thresholds.",
        },
        Card {
            title: "Global log-prob",
            cls:   "info",
            verdict: format!(
                "Global best = {:.0}. Compare runs with identical k and seed to assess \
                 model fit. Higher (less negative) = better fit.",
                data.global_best_lp
            ),
            desc: "The global best log-probability is the mixture model objective value. \
                   Its absolute value is not interpretable in isolation \u{2014} only \
                   relative comparisons (across k values with same data) are meaningful.",
        },
    ];

    cards.iter().map(|c| {
        format!(
            r#"<div class="interp-card">
  <div class="ic-title">{}</div>
  <p>{}</p>
  <div class="ic-verdict {}">{}</div>
</div>
"#,
            c.title, c.desc, c.cls, c.verdict
        )
    }).collect()
}
