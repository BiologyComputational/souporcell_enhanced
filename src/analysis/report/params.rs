// ============================================================================
// report/params.rs — Full run parameters table (§9)
// ============================================================================

use crate::analysis::report::ReportData;

pub fn build_params_table(data: &ReportData) -> String {
    let p = &data.params;

    let rows: &[(&str, String, &str)] = &[
        ("ref_matrix",       p.ref_mtx.clone(),
         "Sparse ref-allele UMI count matrix (.mtx)"),
        ("alt_matrix",       p.alt_mtx.clone(),
         "Sparse alt-allele UMI count matrix (.mtx)"),
        ("barcodes",         p.barcodes.clone(),
         "Cell barcode whitelist (one per line)"),
        ("k (num_clusters)", p.num_clusters.to_string(),
         "Number of donors to demultiplex"),
        ("restarts",         p.restarts.to_string(),
         "EM random restarts for global optimum search"),
        ("seed",             p.seed.to_string(),
         "RNG seed for reproducibility"),
        ("clustering_method",p.clustering_method.name().to_string(),
         "EM variant (em / khm)"),
        ("init_strategy",    format!("{:?}", p.initialization_strategy),
         "Cluster centre initialisation method"),
        ("threads",          p.threads.to_string(),
         "Parallel threads for restarts"),
        ("min_ref",          p.min_ref.to_string(),
         "Min ref-allele reads per locus to include"),
        ("min_alt",          p.min_alt.to_string(),
         "Min alt-allele reads per locus to include"),
        ("min_ref_umis",     p.min_ref_umis.to_string(),
         "Min ref UMIs per cell"),
        ("min_alt_umis",     p.min_alt_umis.to_string(),
         "Min alt UMIs per cell"),
        ("verbose",          p.verbose.to_string(),
         "Emit per-iteration convergence logs to stderr"),
        ("souporcell3",      p.souporcell3.to_string(),
         "Use v3 compatibility mode"),
        ("plot_dir",         p.plot_dir.clone()
                              .unwrap_or_else(|| "(disabled)".to_string()),
         "Output directory for SVGs and HTML report"),
        ("plots",            p.plots.clone(),
         "Comma-separated plot list, or 'all' / 'none'"),
        ("cleanup",          p.cleanup.to_string(),
         "Remove intermediate files after run"),
        ("known_genotypes",  p.known_genotypes.clone()
                              .unwrap_or_else(|| "(none)".to_string()),
         "Optional known-genotype VCF for supervised mode"),
    ];

    let mut t = "\
        <div class='table-wrap'><table>\
        <thead><tr>\
          <th>Parameter</th><th>Value</th><th>Purpose</th>\
        </tr></thead>\
        <tbody>".to_string();

    for (key, val, desc) in rows {
        t += &format!(
            "<tr>\
             <td class='mono'>{key}</td>\
             <td class='mono'>{val}</td>\
             <td style='color:var(--muted);font-size:12px'>{desc}</td>\
             </tr>"
        );
    }
    t += "</tbody></table></div>";
    t
}
