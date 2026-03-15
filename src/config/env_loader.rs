// ============================================================================
// config/env_loader.rs — .souporcell.env file loader
// ============================================================================
//
// Reads a simple KEY=VALUE env file (shell-style, no export, no $expansion).
// Returns a HashMap<String,String> of every key found.
//
// Auto-discovery order (first found wins):
//   1. Path passed via --env flag (explicit)
//   2. .souporcell.env in the current working directory
//   3. ~/.souporcell.env in the user home directory
//   4. Nothing — silently returns empty map, all defaults apply
//
// Syntax rules (intentionally minimal):
//   - Lines starting with # are comments (whole line)
//   - Blank lines are ignored
//   - KEY=VALUE  (no spaces around =, value can contain spaces)
//   - Values may be optionally quoted with " or ' (quotes stripped)
//   - Inline comments (trailing # ...) are NOT supported — keep it simple
//
// All recognised keys are prefixed SOUPC_ to avoid collision with system env.
// ============================================================================

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

// ── Public entry point ────────────────────────────────────────────────────────

/// Load env vars from `explicit_path` (if given) or via auto-discovery.
/// Returns an empty map — never panics — if no file is found or readable.
pub fn load_env(explicit_path: Option<&str>) -> HashMap<String, String> {
    let path = match explicit_path {
        Some(p) => {
            let pb = PathBuf::from(p);
            if pb.exists() {
                Some(pb)
            } else {
                eprintln!("[ENV] Warning: --env path '{}' not found, skipping.", p);
                None
            }
        }
        None => discover_env_file(),
    };

    match path {
        Some(p) => {
            eprintln!("[ENV] Loading environment from: {}", p.display());
            parse_env_file(&p)
        }
        None => {
            eprintln!("[ENV] No .souporcell.env file found — using CLI/config/defaults.");
            HashMap::new()
        }
    }
}

// ── File discovery ────────────────────────────────────────────────────────────

fn discover_env_file() -> Option<PathBuf> {
    // 1. CWD
    let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    let cwd_env = cwd.join(".souporcell.env");
    if cwd_env.exists() {
        return Some(cwd_env);
    }
    // 2. Home directory
    if let Some(home) = home_dir() {
        let home_env = home.join(".souporcell.env");
        if home_env.exists() {
            return Some(home_env);
        }
    }
    None
}

fn home_dir() -> Option<PathBuf> {
    std::env::var("HOME").ok().map(PathBuf::from)
        .or_else(|| std::env::var("USERPROFILE").ok().map(PathBuf::from))
}

// ── File parser ───────────────────────────────────────────────────────────────

fn parse_env_file(path: &Path) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let content = match fs::read_to_string(path) {
        Ok(c)  => c,
        Err(e) => {
            eprintln!("[ENV] Warning: could not read '{}': {}", path.display(), e);
            return map;
        }
    };

    for (lineno, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        match line.split_once('=') {
            Some((key, val)) => {
                let k = key.trim().to_uppercase();
                let v = expand_tilde_str(strip_quotes(val.trim()));
                if k.starts_with("SOUPC_") {
                    map.insert(k, v);
                } else {
                    eprintln!(
                        "[ENV] Line {}: key '{}' ignored (must start with SOUPC_)",
                        lineno + 1, k
                    );
                }
            }
            None => {
                eprintln!("[ENV] Line {}: malformed (no '='), skipped: {}", lineno + 1, line);
            }
        }
    }

    eprintln!("[ENV] Loaded {} key(s) from env file.", map.len());
    map
}

fn strip_quotes(s: &str) -> &str {
    if (s.starts_with('"') && s.ends_with('"'))
        || (s.starts_with('\'') && s.ends_with('\''))
    {
        &s[1..s.len() - 1]
    } else {
        s
    }
}

/// Expand a leading ~ to $HOME in a string value from the env file.
fn expand_tilde_str(s: &str) -> String {
    if s == "~" {
        return home_dir().map(|p| p.to_string_lossy().to_string())
                         .unwrap_or_else(|| "~".to_string());
    }
    if s.starts_with("~/") || s.starts_with("~\\") {
        if let Some(home) = home_dir() {
            let home_str = home.to_string_lossy();
            return format!("{}/{}", home_str.trim_end_matches('/'), &s[2..]);
        }
    }
    s.to_string()
}

// ── Typed accessors (used by params.rs merge logic) ──────────────────────────

pub struct EnvMap(pub HashMap<String, String>);

impl EnvMap {
    pub fn str(&self, key: &str) -> Option<&str> {
        self.0.get(key).map(String::as_str)
    }

    pub fn parse<T>(&self, key: &str) -> Option<T>
    where
        T: std::str::FromStr,
        T::Err: std::fmt::Debug,
    {
        self.0.get(key).and_then(|v| {
            v.parse::<T>().map_err(|e| {
                eprintln!("[ENV] Warning: could not parse {}={:?}: {:?}", key, v, e);
            }).ok()
        })
    }

    pub fn bool(&self, key: &str) -> Option<bool> {
        self.0.get(key).and_then(|v| match v.to_lowercase().as_str() {
            "true" | "1" | "yes" => Some(true),
            "false" | "0" | "no"  => Some(false),
            other => {
                eprintln!("[ENV] Warning: {}={:?} is not a valid bool, ignoring.", key, other);
                None
            }
        })
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}
