// ============================================================================
// report/methods.rs — Algorithmic methods and references section (§10)
// ============================================================================

use crate::analysis::report::ReportData;

pub fn build_methods_section(data: &ReportData) -> String {
    let method   = data.params.clustering_method.name().to_uppercase();
    let restarts = data.total_restarts;

    format!(r#"<div class="section-label">Methods &amp; References</div>
<h2>Algorithmic Methods &amp; Citations</h2>
<div class="two-col">
<div>
  <h3>Algorithm Summary</h3>
  <div class="card">
    <p style="font-size:13px;color:var(--ink2);line-height:1.8">
      souporcell v2.6 performs <em>unsupervised genotype demultiplexing</em> of pooled
      scRNA-seq libraries. Each cell barcode is modelled as drawn from one of
      <strong>k</strong> donor genotype clusters. The allele observation model is
      binomial: given donor allele frequency &theta;, the probability of observing
      <em>r</em> ref reads out of <em>r+a</em> total is Binomial(<em>r+a</em>,&nbsp;&theta;).
      The mixture log-likelihood is maximised via <strong>{method}</strong> with
      simulated annealing and <strong>{restarts}</strong> random restarts. The global
      optimum is selected by highest log-likelihood; cells are assigned by argmax over
      posterior donor probabilities.
    </p>
  </div>

  <h3>Quality Thresholds Used</h3>
  <div class="table-wrap"><table>
    <thead><tr><th>Metric</th><th>Threshold</th><th>Rationale</th></tr></thead>
    <tbody>
      <tr><td>Loci QC pass rate</td><td>&gt;1%</td>
          <td>Minimum informative sites for clustering</td></tr>
      <tr><td>Matrix sparsity</td><td>&gt;90%</td>
          <td>Expected for scRNA-seq allele data</td></tr>
      <tr><td>Convergence stability</td><td>&ge;70%</td>
          <td>Reliable global optimum identification</td></tr>
      <tr><td>Cluster imbalance</td><td>&lt;50%</td>
          <td>No cluster should dominate 3&times; over others</td></tr>
      <tr><td>Median cell confidence</td><td>&gt;15 log-units</td>
          <td>Assignment certainty above noise floor</td></tr>
      <tr><td>Ambiguous cells</td><td>&lt;5%</td>
          <td>Fraction with confidence &lt;10 log-units</td></tr>
    </tbody>
  </table></div>
</div>

<div>
  <h3>Primary References</h3>
  <ul class="ref-list">
    <li><strong>Heaton H, et al.</strong> (2020). Souporcell: robust clustering of
        single-cell RNA-seq data by genotype without reference genotypes.
        <em>Nature Methods</em> 17, 615&ndash;620.
        doi:10.1038/s41592-020-0820-1</li>
    <li><strong>Dempster AP, Laird NM, Rubin DB.</strong> (1977). Maximum likelihood
        from incomplete data via the EM algorithm.
        <em>J R Stat Soc B</em> 39(1):1&ndash;38.</li>
    <li><strong>Kang HM, et al.</strong> (2018). Multiplexed droplet single-cell
        RNA-sequencing using natural genetic variation.
        <em>Nature Biotechnology</em> 36, 89&ndash;94.</li>
    <li><strong>McCarthy DJ, et al.</strong> (2020). Cardelino: computational
        integration of somatic clonal inference for single cells.
        <em>Nature Methods</em>.</li>
    <li><strong>Stoeckius M, et al.</strong> (2018). Cell Hashing with barcoded
        antibodies enables multiplexing and doublet detection for single cell
        genomics. <em>Genome Biology</em> 19:224.</li>
  </ul>

  <h3>Software Dependencies</h3>
  <div class="table-wrap"><table>
    <thead><tr><th>Tool</th><th>Version</th><th>Role</th></tr></thead>
    <tbody>
      <tr><td>Rust</td><td class="mono">&ge;1.75</td>
          <td>Core algorithm runtime</td></tr>
      <tr><td>D3.js</td><td class="mono">7.8.5</td>
          <td>Interactive pipeline graph</td></tr>
      <tr><td>rayon</td><td class="mono">1.x</td>
          <td>Parallel restart execution</td></tr>
      <tr><td>clap</td><td class="mono">3.x</td>
          <td>CLI argument parsing</td></tr>
      <tr><td>serde_json</td><td class="mono">1.x</td>
          <td>Graph template deserialisation</td></tr>
    </tbody>
  </table></div>
</div>
</div>
"#,
        method   = method,
        restarts = restarts,
    )
}
