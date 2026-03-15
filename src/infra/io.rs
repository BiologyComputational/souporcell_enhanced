// ============================================================================
// io.rs — All file I/O: matrix loading, barcode loading, result writing
// ============================================================================
//
// Enhancements vs original:
//   - load_cell_data returns a LoadResult struct instead of a 5-tuple
//   - Zero-coverage count is returned to logger (was an ambiguous eprintln!)
//   - Both .gz and plain files handled via the reader() helper
//   - write_cluster_assignments() extracted from main so output logic is
//     testable independently of the full pipeline
//   - Per-locus coverage statistics computed and returned for the summary
// ============================================================================

use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use hashbrown::{HashMap, HashSet};
use statrs::function::factorial;

use crate::config::params::Params;
use crate::domain::types::CellData;

// ── Returned by load_cell_data ────────────────────────────────────────────────
pub struct LoadResult {
    pub loci_used:        usize,
    pub total_cells:      usize,
    pub cell_data:        Vec<CellData>,
    /// Maps compressed index → original locus id (for VCF reconstruction)
    #[allow(dead_code)]
    pub index_to_locus:   Vec<usize>,
    /// Maps original locus id → compressed index (for known-genotype init)
    pub locus_to_index:   HashMap<usize, usize>,
    /// Number of (locus, cell) pairs with zero coverage — passed to logger
    pub zero_cov_entries: usize,
    /// Total loci in the raw matrix before QC filters
    pub total_loci_raw:   usize,
}

// ── reader ────────────────────────────────────────────────────────────────────
//
// Returns a boxed BufRead that transparently handles .gz files.
// Also expands a leading ~ to $HOME so paths from .env / config profiles work.
pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let expanded = expand_tilde(filename);
    let path = Path::new(&expanded);
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("Could not open file '{}': {}", expanded, e));
    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(128 * 1024, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

// ── expand_tilde ──────────────────────────────────────────────────────────────
//
// Replaces a leading ~ with the value of $HOME (or $USERPROFILE on Windows).
// Handles ~/path and ~ alone. Does NOT expand ~username.
fn expand_tilde(path: &str) -> String {
    if path == "~" {
        return home_dir().unwrap_or_else(|| "~".to_string());
    }
    if path.starts_with("~/") || path.starts_with("~\\") {
        if let Some(home) = home_dir() {
            return format!("{}/{}", home.trim_end_matches('/'), &path[2..]);
        }
    }
    path.to_string()
}

fn home_dir() -> Option<String> {
    std::env::var("HOME").ok()
        .or_else(|| std::env::var("USERPROFILE").ok())
}

// ── load_barcodes ─────────────────────────────────────────────────────────────
//
// Reads one barcode per line (supports .gz).
// The order here defines the cell column order — must match the mtx files.
pub fn load_barcodes(params: &Params) -> Vec<String> {
    reader(&params.barcodes)
        .lines()
        .map(|l| l.expect("Could not read barcode line"))
        .filter(|l| !l.is_empty())
        .collect()
}

// ── load_cell_data ────────────────────────────────────────────────────────────
//
// Parses ref.mtx and alt.mtx in lockstep (they must be identical in
// row/column ordering, which vartrix guarantees).
//
// Two passes:
//   Pass 1 — aggregate per-locus cell counts and UMI counts to apply
//             min_ref / min_alt / min_ref_umis / min_alt_umis filters
//   Pass 2 — build CellData entries only for loci that pass QC
pub fn load_cell_data(params: &Params) -> LoadResult {
    // ── Pass 1: collect per-locus statistics ─────────────────────────────────
    let alt_reader = File::open(expand_tilde(&params.alt_mtx))
        .unwrap_or_else(|e| panic!("Cannot open alt matrix '{}': {}", params.alt_mtx, e));
    let ref_reader = File::open(expand_tilde(&params.ref_mtx))
        .unwrap_or_else(|e| panic!("Cannot open ref matrix '{}': {}", params.ref_mtx, e));
    let mut alt_lines = BufReader::new(alt_reader).lines();
    let mut ref_lines = BufReader::new(ref_reader).lines();

    let mut total_loci  = 0usize;
    let mut total_cells = 0usize;

    // per-locus: [cells_with_ref, cells_with_alt]
    let mut locus_cell_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    // per-locus: [total_ref_umis, total_alt_umis]
    let mut locus_umi_counts:  HashMap<usize, [u32; 2]> = HashMap::new();
    // per-locus per-cell: [ref_count, alt_count]
    let mut locus_counts: HashMap<usize, HashMap<usize, [u32; 2]>> = HashMap::new();
    let mut all_loci: HashSet<usize> = HashSet::new();

    // Skip the two header lines (% comment and blank), then the dimension line
    for line_number in 0.. {
        let alt_line = match alt_lines.next() {
            Some(l) => l.expect("alt mtx read error"),
            None    => break,
        };
        let ref_line = ref_lines.next()
            .expect("ref mtx ended before alt mtx")
            .expect("ref mtx read error");

        if line_number < 2 {
            // MatrixMarket comment / format lines
            continue;
        }
        if line_number == 2 {
            // Dimension line: rows cols nnz
            let tokens: Vec<&str> = alt_line.split_whitespace().collect();
            total_loci  = tokens[0].parse::<usize>().unwrap();
            total_cells = tokens[1].parse::<usize>().unwrap();
            continue;
        }
        // Data lines: row col value  (1-indexed)
        let at: Vec<&str> = alt_line.split_whitespace().collect();
        let rt: Vec<&str> = ref_line.split_whitespace().collect();
        let locus = at[0].parse::<usize>().unwrap() - 1;
        let cell  = at[1].parse::<usize>().unwrap() - 1;
        let ref_count = rt[2].parse::<u32>().unwrap();
        let alt_count = at[2].parse::<u32>().unwrap();

        all_loci.insert(locus);

        let cc = locus_cell_counts.entry(locus).or_insert([0u32; 2]);
        let uc = locus_umi_counts.entry(locus).or_insert([0u32; 2]);
        if ref_count > 0 { cc[0] += 1; uc[0] += ref_count; }
        if alt_count > 0 { cc[1] += 1; uc[1] += alt_count; }

        locus_counts
            .entry(locus)
            .or_insert_with(HashMap::new)
            .insert(cell, [ref_count, alt_count]);
    }

    // ── Sort loci for deterministic index assignment ──────────────────────────
    let mut sorted_loci: Vec<usize> = all_loci.into_iter().collect();
    sorted_loci.sort_unstable();

    // ── Pass 2: build CellData for QC-passing loci ───────────────────────────
    let mut used_loci:     HashSet<usize>           = HashSet::new();
    let mut index_to_locus: Vec<usize>              = Vec::new();
    let mut locus_to_index: HashMap<usize, usize>   = HashMap::new();
    let mut cell_data:      Vec<CellData>           = (0..total_cells).map(|_| CellData::new()).collect();
    let mut zero_cov_entries = 0usize;

    let mut locus_index = 0usize;
    for locus in sorted_loci {
        let cc = locus_cell_counts.get(&locus).unwrap();
        let uc = locus_umi_counts.get(&locus).unwrap();

        if cc[0] >= params.min_ref
            && cc[1] >= params.min_alt
            && uc[0] >= params.min_ref_umis
            && uc[1] >= params.min_alt_umis
        {
            used_loci.insert(locus);
            index_to_locus.push(locus);
            locus_to_index.insert(locus, locus_index);

            let cell_map = locus_counts.get(&locus).unwrap();
            for cell in 0..total_cells {
                match cell_map.get(&cell) {
                    None => {
                        // Cell has no reads at this locus — contributes to
                        // the zero-coverage count but does not add data
                        zero_cov_entries += 1;
                    }
                    Some(&[ref_c, alt_c]) => {
                        if ref_c + alt_c == 0 {
                            zero_cov_entries += 1;
                            continue;
                        }
                        let total = (ref_c + alt_c) as u64;
                        let cd = &mut cell_data[cell];
                        cd.alt_counts.push(alt_c);
                        cd.ref_counts.push(ref_c);
                        cd.loci.push(locus_index);
                        cd.allele_fractions.push(alt_c as f32 / (ref_c + alt_c) as f32);
                        cd.log_binomial_coefficient.push(
                            factorial::ln_binomial(total, alt_c as u64) as f32
                        );
                        cd.total_alleles += (ref_c + alt_c) as f32;
                    }
                }
            }
            locus_index += 1;
        }
    }

    LoadResult {
        loci_used:        used_loci.len(),
        total_cells,
        cell_data,
        index_to_locus,
        locus_to_index,
        zero_cov_entries,
        total_loci_raw:   total_loci,
    }
}

// ── write_cluster_assignments ────────────────────────────────────────────────
//
// Writes the clusters_tmp.tsv to stdout (captured by the pipeline).
// Format (one row per barcode, tab-separated):
//   barcode   best_cluster   lp_cluster0   lp_cluster1 ...
//
// This is the format that troublet and consensus.py expect as input.
pub fn write_cluster_assignments(
    barcodes:            &[String],
    best_log_probabilities: &[Vec<f32>],
) {
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    for (bc, log_probs) in barcodes.iter().zip(best_log_probabilities.iter()) {
        let best = log_probs.iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(i, _)| i)
            .unwrap_or(0);
        write!(out, "{}\t{}", bc, best).unwrap();
        for lp in log_probs {
            write!(out, "\t{}", lp).unwrap();
        }
        writeln!(out).unwrap();
    }
}
