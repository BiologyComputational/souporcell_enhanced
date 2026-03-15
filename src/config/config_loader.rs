// ============================================================================
// config/config_loader.rs — JSON config profile loader
// ============================================================================
//
// Reads a named JSON profile passed via --config path/to/profile.json.
// The profile sets any subset of params — fields not present fall through to
// .env values or built-in defaults.
//
// JSON schema (all fields optional except for documentation):
// {
//   "profile":     "k4_donors",             // human label, informational only
//   "description": "4-donor GRCh38 run",    // informational only
//   "params": {
//     "ref_matrix":              "/path/to/ref.mtx",
//     "alt_matrix":              "/path/to/alt.mtx",
//     "barcodes":                "/path/to/barcodes.tsv",
//     "num_clusters":            4,
//     "restarts":                200,
//     "seed":                    42,
//     "threads":                 8,
//     "clustering_method":       "em",
//     "initialization_strategy": "random_uniform",
//     "souporcell3":             false,
//     "min_alt":                 4,
//     "min_ref":                 4,
//     "min_alt_umis":            0,
//     "min_ref_umis":            0,
//     "verbose":                 false,
//     "progress_interval":       500,
//     "plot_dir":                "./output",
//     "plots":                   "all",
//     "cleanup":                 false,
//     "known_genotypes":         null,
//     "known_genotypes_sample_names": []
//   }
// }
//
// Parsing is done manually (no serde dependency) using a minimal JSON
// value extractor so we don't need to add serde_json to Cargo.toml.
// ============================================================================

use std::collections::HashMap;
use std::fs;
use std::path::Path;

// ── Public profile struct ─────────────────────────────────────────────────────

pub struct ConfigProfile {
    pub profile_name: String,
    pub description:  String,
    /// Flat string map of param_name → raw_string_value.
    /// Consumers cast via typed accessors below.
    pub values:       HashMap<String, String>,
}

impl ConfigProfile {
    pub fn str_val(&self, key: &str) -> Option<&str> {
        self.values.get(key).map(String::as_str)
    }

    pub fn parse<T>(&self, key: &str) -> Option<T>
    where
        T: std::str::FromStr,
        T::Err: std::fmt::Debug,
    {
        self.values.get(key).and_then(|v| {
            if v == "null" { return None; }
            v.parse::<T>().map_err(|e| {
                eprintln!("[CONFIG] Warning: could not parse {}={:?}: {:?}", key, v, e);
            }).ok()
        })
    }

    pub fn bool_val(&self, key: &str) -> Option<bool> {
        self.values.get(key).and_then(|v| match v.as_str() {
            "true"  => Some(true),
            "false" => Some(false),
            "null"  => None,
            other   => {
                eprintln!("[CONFIG] Warning: {}={:?} is not a valid bool.", key, other);
                None
            }
        })
    }

    pub fn str_vec(&self, key: &str) -> Vec<String> {
        // Stored as comma-separated in the flat map after array extraction
        self.values.get(key)
            .map(|v| {
                if v.is_empty() || v == "null" { return vec![]; }
                v.split(',').map(|s| s.trim().to_string()).collect()
            })
            .unwrap_or_default()
    }
}

// ── Public entry point ────────────────────────────────────────────────────────

pub fn load_config(path: &str) -> ConfigProfile {
    let p = Path::new(path);
    if !p.exists() {
        eprintln!("[CONFIG] Error: --config file '{}' not found.", path);
        std::process::exit(1);
    }

    let content = fs::read_to_string(p).unwrap_or_else(|e| {
        eprintln!("[CONFIG] Error: cannot read '{}': {}", path, e);
        std::process::exit(1);
    });

    eprintln!("[CONFIG] Loading profile from: {}", path);
    let profile = parse_profile(&content);
    eprintln!(
        "[CONFIG] Profile '{}' loaded — {} param(s) set.",
        profile.profile_name,
        profile.values.len()
    );
    profile
}

// ── Minimal JSON parser (no serde dependency) ─────────────────────────────────
//
// Supports only the flat subset we need:
//   - top-level "profile" and "description" strings
//   - "params" object containing string / number / bool / null / string-array values
//
// This is intentionally NOT a general-purpose JSON parser.

fn parse_profile(json: &str) -> ConfigProfile {
    let mut profile_name = String::from("(unnamed)");
    let mut description  = String::new();
    let mut values       = HashMap::new();

    // Extract top-level "profile" string
    if let Some(v) = extract_string(json, "profile") {
        profile_name = v;
    }
    if let Some(v) = extract_string(json, "description") {
        description = v;
    }

    // Find the "params": { ... } block
    if let Some(params_block) = extract_object_block(json, "params") {
        // String fields
        for key in &[
            "ref_matrix", "alt_matrix", "barcodes",
            "clustering_method", "initialization_strategy",
            "plot_dir", "plots", "known_genotypes",
        ] {
            if let Some(v) = extract_string(&params_block, key) {
                values.insert(key.to_string(), v);
            }
        }

        // Numeric fields (stored as raw string, cast later)
        for key in &[
            // Core clustering
            "num_clusters", "restarts", "seed", "threads",
            // Locus QC
            "min_alt", "min_ref", "min_alt_umis", "min_ref_umis",
            // Annealing schedule
            "anneal_steps", "anneal_base_em", "anneal_base_khm",
            // EM/KHM convergence
            "conv_tol", "max_iter", "min_cluster_cells",
            // M-step pseudocounts
            "pseudocount_alt", "pseudocount_total",
            // Allele-frequency bounds
            "theta_min", "theta_max",
            // KHM
            "khm_p",
            // Logging
            "progress_interval",
        ] {
            if let Some(v) = extract_number(&params_block, key) {
                values.insert(key.to_string(), v);
            }
        }

        // Bool fields
        for key in &["souporcell3", "verbose", "cleanup"] {
            if let Some(v) = extract_bool(&params_block, key) {
                values.insert(key.to_string(), v.to_string());
            }
        }

        // String-array field
        if let Some(arr) = extract_string_array(&params_block, "known_genotypes_sample_names") {
            values.insert(
                "known_genotypes_sample_names".to_string(),
                arr.join(","),
            );
        }
    }

    ConfigProfile { profile_name, description, values }
}

// ── JSON value extractors ─────────────────────────────────────────────────────

/// Extract `"key": "value"` → value string (quotes stripped).
fn extract_string(json: &str, key: &str) -> Option<String> {
    let needle = format!("\"{}\"", key);
    let start  = json.find(&needle)?;
    let after_key = &json[start + needle.len()..];
    let colon_pos = after_key.find(':')? + 1;
    let after_colon = after_key[colon_pos..].trim_start();

    if after_colon.starts_with("null") {
        return None;
    }
    if !after_colon.starts_with('"') {
        return None; // it's a number or bool, not a string
    }
    // Find closing quote (simple — doesn't handle escaped quotes in values)
    let inner = &after_colon[1..];
    let end = inner.find('"')?;
    Some(inner[..end].to_string())
}

/// Extract `"key": 123` → "123" as String.
fn extract_number(json: &str, key: &str) -> Option<String> {
    let needle = format!("\"{}\"", key);
    let start  = json.find(&needle)?;
    let after_key = &json[start + needle.len()..];
    let colon_pos = after_key.find(':')? + 1;
    let after_colon = after_key[colon_pos..].trim_start();

    if after_colon.starts_with("null") || after_colon.starts_with('"') {
        return None;
    }
    // Read digits (including - and .)
    let num: String = after_colon.chars()
        .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '.')
        .collect();
    if num.is_empty() { None } else { Some(num) }
}

/// Extract `"key": true/false` → bool.
fn extract_bool(json: &str, key: &str) -> Option<bool> {
    let needle = format!("\"{}\"", key);
    let start  = json.find(&needle)?;
    let after_key = &json[start + needle.len()..];
    let colon_pos = after_key.find(':')? + 1;
    let after_colon = after_key[colon_pos..].trim_start();

    if after_colon.starts_with("true")  { return Some(true);  }
    if after_colon.starts_with("false") { return Some(false); }
    None
}

/// Extract `"key": ["a","b","c"]` → Vec<String>.
fn extract_string_array(json: &str, key: &str) -> Option<Vec<String>> {
    let needle = format!("\"{}\"", key);
    let start  = json.find(&needle)?;
    let after_key = &json[start + needle.len()..];
    let colon_pos = after_key.find(':')? + 1;
    let after_colon = after_key[colon_pos..].trim_start();

    if !after_colon.starts_with('[') { return None; }
    let end = after_colon.find(']')?;
    let inner = &after_colon[1..end];

    let items: Vec<String> = inner.split(',')
        .filter_map(|s| {
            let s = s.trim();
            if s.starts_with('"') && s.ends_with('"') {
                Some(s[1..s.len()-1].to_string())
            } else if s.is_empty() || s == "null" {
                None
            } else {
                Some(s.to_string())
            }
        })
        .collect();

    Some(items)
}

/// Extract the body of `"key": { ... }` as a string slice.
fn extract_object_block(json: &str, key: &str) -> Option<String> {
    let needle = format!("\"{}\"", key);
    let start  = json.find(&needle)?;
    let after_key = &json[start + needle.len()..];
    let colon_pos = after_key.find(':')? + 1;
    let after_colon = after_key[colon_pos..].trim_start();

    if !after_colon.starts_with('{') { return None; }

    // Walk forward counting braces to find the matching }
    let mut depth = 0usize;
    let mut end   = 0usize;
    for (i, ch) in after_colon.char_indices() {
        match ch {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 { end = i; break; }
            }
            _ => {}
        }
    }
    Some(after_colon[1..end].to_string())
}
