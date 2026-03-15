// ============================================================================
// infra/gui_server.rs — Web GUI server for --gui mode
// ============================================================================
//
// Launched when the user runs:  souporcell --gui  [--gui_port 7979]
//
// Uses ONLY std::net — zero new Cargo.toml dependencies.
//
// Endpoints:
//   GET  /                   → serves the embedded souporcell_gui.html
//   GET  /api/version        → {"version":"2.7.0","binary":"..."}
//   GET  /api/status         → {"status":"ready"}
//   POST /api/run            → spawns souporcell with POSTed JSON config,
//                              streams stdout+stderr back as Server-Sent Events
//   GET  /api/files?dir=...  → lists files in the given directory
//   GET  /files/...          → serves a file from the given absolute path
//
// The GUI HTML is embedded at compile time so the binary is fully self-contained.
// ============================================================================

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;

// ── Embedded GUI ──────────────────────────────────────────────────────────────
// The GUI HTML is embedded at compile time from the file next to this source.
// The binary path is resolved at runtime.
static GUI_HTML: &str = include_str!("souporcell_gui.html");

// ── Embedded analysis page ────────────────────────────────────────────────────
static ALGO_HTML: &str = include_str!("souporcell_algo_analysis.html");


pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");

// ── Active run state ──────────────────────────────────────────────────────────
struct RunState {
    running: bool,
    log:     Vec<String>,
}

impl RunState {
    fn new() -> Self { RunState { running: false, log: Vec::new() } }
}

// ── Public entry point ────────────────────────────────────────────────────────

pub fn launch(port: u16) {
    let addr = format!("127.0.0.1:{}", port);
    let listener = TcpListener::bind(&addr).unwrap_or_else(|e| {
        eprintln!("GUI: cannot bind to {} — {}", addr, e);
        std::process::exit(1);
    });

    // Resolve the path to our own binary
    let binary_path = std::env::current_exe()
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_else(|_| "souporcell".to_string());

    eprintln!("\n╔══════════════════════════════════════════════════════════╗");
    eprintln!("║   souporcell v{} — GUI Mode                          ║", VERSION);
    eprintln!("╠══════════════════════════════════════════════════════════╣");
    eprintln!("║                                                          ║");
    eprintln!("║   Open in your browser:                                  ║");
    eprintln!("║   \x1b[36mhttp://localhost:{}\x1b[0m                              ║", port);
    eprintln!("║                                                          ║");
    eprintln!("║   Press  Ctrl+C  to stop the server.                     ║");
    eprintln!("╚══════════════════════════════════════════════════════════╝\n");

    // Try to auto-open the browser
    let url = format!("http://localhost:{}", port);
    let _ = open_browser(&url);

    let run_state = Arc::new(Mutex::new(RunState::new()));
    let binary    = Arc::new(binary_path);

    for stream in listener.incoming() {
        match stream {
            Ok(s) => {
                let rs  = Arc::clone(&run_state);
                let bin = Arc::clone(&binary);
                thread::spawn(move || handle_connection(s, rs, bin));
            }
            Err(_) => {}
        }
    }
}

// ── Connection handler ────────────────────────────────────────────────────────

fn handle_connection(
    mut stream: TcpStream,
    run_state:  Arc<Mutex<RunState>>,
    binary:     Arc<String>,
) {
    let mut buf = [0u8; 16384];
    let n = match stream.read(&mut buf) {
        Ok(n) if n > 0 => n,
        _ => return,
    };
    let raw = String::from_utf8_lossy(&buf[..n]);
    let (method, path, body) = parse_request(&raw);

    match (method.as_str(), path.as_str()) {
        // ── Serve GUI HTML ────────────────────────────────────────────────────
        ("GET", "/") | ("GET", "/index.html") => {
            let html = inject_version(GUI_HTML);
            respond_html(&mut stream, 200, &html);
        }

        // ── Algorithm analysis page ───────────────────────────────────────────
        ("GET", "/algo-analysis") | ("GET", "/algo-analysis.html") => {
            respond_html(&mut stream, 200, ALGO_HTML);
        }

        // ── Version JSON ──────────────────────────────────────────────────────
        ("GET", "/api/version") => {
            let json = format!(
                r#"{{"version":"{}","binary":"{}","pkg":"{}"}}"#,
                VERSION, &*binary, PKG_NAME
            );
            respond_json(&mut stream, 200, &json);
        }

        // ── Status ───────────────────────────────────────────────────────────
        ("GET", "/api/status") => {
            let running = run_state.lock().map(|s| s.running).unwrap_or(false);
            respond_json(&mut stream,200,
                &format!(r#"{{"status":"{}"}}"#, if running {"running"} else {"ready"}));
        }

        // ── Run souporcell ───────────────────────────────────────────────────
        ("POST", "/api/run") => {
            let args = parse_config_to_args(&body);
            {
                let mut state = run_state.lock().unwrap();
                if state.running {
                    respond_json(&mut stream, 409,
                        r#"{"error":"A run is already in progress"}"#);
                    return;
                }
                state.running = true;
                state.log.clear();
            }

            // Send SSE headers immediately
            let sse_header = "HTTP/1.1 200 OK\r\n\
                Content-Type: text/event-stream\r\n\
                Cache-Control: no-cache\r\n\
                Access-Control-Allow-Origin: *\r\n\
                Connection: keep-alive\r\n\r\n";
            if stream.write_all(sse_header.as_bytes()).is_err() {
                return;
            }

            let bin_path = (*binary).clone();
            let rs = Arc::clone(&run_state);

            // ── Determine output TSV path ─────────────────────────────────────
            // stdout of the souporcell binary is the clusters TSV.
            // We must NOT use Stdio::piped() for stdout — if the pipe buffer
            // fills up (large datasets) the child blocks and never exits,
            // causing child.wait() to hang forever.
            // Instead we redirect stdout directly to a file on disk.
            fn str_val_local(body: &str, key: &str) -> Option<String> {
                let needle = format!("\"{}\":", key);
                let start = body.find(&needle)?;
                let after = body[start + needle.len()..].trim_start();
                if after.starts_with('"') {
                    let inner = &after[1..];
                    let end = inner.find('"')?;
                    let val = inner[..end].trim().to_string();
                    if val.is_empty() { None } else { Some(val) }
                } else { None }
            }
            let plot_dir = str_val_local(&body, "plot_dir")
                .unwrap_or_else(|| ".".to_string());
            // Create the output directory if needed
            let _ = std::fs::create_dir_all(&plot_dir);
            let tsv_path = format!("{}/clusters_tmp.tsv", plot_dir.trim_end_matches('/'));
            let err_path = format!("{}/clusters.err",     plot_dir.trim_end_matches('/'));

            // Open stdout and stderr redirect files
            let tsv_file = std::fs::File::create(&tsv_path);
            let err_file = std::fs::File::create(&err_path);

            let stdout_redirect = match tsv_file {
                Ok(f)  => Stdio::from(f),
                Err(_) => Stdio::null(),
            };
            // We still need stderr as a pipe for SSE streaming
            // (err_file is also written via tee-style below)

            // Spawn the real souporcell binary
            // stdout → file  (never blocks the child regardless of output size)
            // stderr → pipe  (we stream line-by-line to the browser)
            let child = Command::new(&bin_path)
                .args(&args)
                .stdout(stdout_redirect)
                .stderr(Stdio::piped())
                .spawn();

            match child {
                Err(e) => {
                    let msg = format!("data: {{\"type\":\"error\",\"msg\":\"Failed to spawn: {}\"}}\n\n", e);
                    let _ = stream.write_all(msg.as_bytes());
                    let mut state = rs.lock().unwrap();
                    state.running = false;
                }
                Ok(mut child) => {
                    // Send start event with TSV destination
                    let start = format!(
                        "data: {{\"type\":\"start\",\"cmd\":\"{}\",\"tsv_path\":\"{}\"}}\n\n",
                        json_escape(&(shell_escape(&bin_path) + " " + &args.join(" "))),
                        json_escape(&tsv_path)
                    );
                    let _ = stream.write_all(start.as_bytes());

                    // Open err file for mirroring stderr to disk as well
                    let mut err_writer = err_file.ok().map(std::io::BufWriter::new);

                    // Stream stderr line by line (all souporcell progress goes to stderr)
                    if let Some(stderr) = child.stderr.take() {
                        let reader = BufReader::new(stderr);
                        for line in reader.lines().flatten() {
                            let escaped = json_escape(&line);
                            // Detect stage transitions
                            let ev_type = if line.starts_with("[PARAMS]") {
                                "params"
                            } else if line.starts_with("══") || line.starts_with("──") {
                                "section"
                            } else if line.starts_with("  ✔") || line.starts_with("  ✓") {
                                "ok"
                            } else if line.starts_with("RESTART") {
                                "restart"
                            } else if line.starts_with("EM\t") || line.starts_with("KHM\t") {
                                "tempstep"
                            } else if line.contains("⚠") || line.starts_with("  ⚠") {
                                "warn"
                            } else {
                                "log"
                            };
                            let event = format!(
                                "data: {{\"type\":\"{}\",\"msg\":\"{}\"}}\n\n",
                                ev_type, escaped
                            );
                            // Mirror every line to the server terminal and err file
                            println!("[GUI-RUN] {}", &line);
                            if let Some(ref mut w) = err_writer {
                                use std::io::Write as _;
                                let _ = w.write_all(line.as_bytes());
                                let _ = w.write_all(b"\n");
                            }
                            if stream.write_all(event.as_bytes()).is_err() { break; }
                            rs.lock().map(|mut s| s.log.push(line)).ok();
                        }
                    }

                    // Wait for exit
                    let exit_code = child.wait()
                        .map(|s| s.code().unwrap_or(-1))
                        .unwrap_or(-1);

                    // Flush err writer
                    drop(err_writer);

                    let tsv_size = std::fs::metadata(&tsv_path)
                        .map(|m| m.len()).unwrap_or(0);
                    let done = format!(
                        "data: {{\"type\":\"done\",\"exit_code\":{},\"tsv_path\":\"{}\",\"tsv_bytes\":{}}}\n\n",
                        exit_code,
                        json_escape(&tsv_path),
                        tsv_size
                    );
                    println!("[GUI-RUN] [DONE] exit={} tsv={} ({} bytes)",
                             exit_code, &tsv_path, tsv_size);
                    let _ = stream.write_all(done.as_bytes());
                    let mut state = rs.lock().unwrap();
                    state.running = false;
                }
            }
        }

        // ── List output files ─────────────────────────────────────────────────
        ("GET", path) if path.starts_with("/api/files") => {
            let dir = extract_query_param(path, "dir")
                .unwrap_or_else(|| ".".to_string());
            let files = list_files(&dir);
            respond_json(&mut stream, 200, &files);
        }

        // ── Serve an output file ──────────────────────────────────────────────
        ("GET", path) if path.starts_with("/files/") => {
            let file_path = percent_decode(&path[7..]);
            serve_file(&mut stream, &file_path);
        }

        // ── CORS preflight ────────────────────────────────────────────────────
        ("OPTIONS", _) => {
            let resp = "HTTP/1.1 204 No Content\r\n\
                Access-Control-Allow-Origin: *\r\n\
                Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n\
                Access-Control-Allow-Headers: Content-Type\r\n\r\n";
            let _ = stream.write_all(resp.as_bytes());
        }

        _ => {
            respond_json(&mut stream, 404, r#"{"error":"not found"}"#);
        }
    }
}

// ── Config JSON → CLI args ────────────────────────────────────────────────────

fn parse_config_to_args(body: &str) -> Vec<String> {
    let mut args: Vec<String> = Vec::new();

    fn str_val(body: &str, key: &str) -> Option<String> {
        let needle = format!("\"{}\":", key);
        let start  = body.find(&needle)?;
        let after  = body[start + needle.len()..].trim_start();
        if after.starts_with('"') {
            let inner = &after[1..];
            let end   = inner.find('"')?;
            let val   = inner[..end].trim().to_string();
            if val.is_empty() { None } else { Some(val) }
        } else {
            // number or bool
            let end = after.find(|c: char| c == ',' || c == '}' || c == '\n')
                           .unwrap_or(after.len());
            let val = after[..end].trim().to_string();
            if val == "null" || val.is_empty() { None } else { Some(val) }
        }
    }

    fn bool_val(body: &str, key: &str) -> bool {
        str_val(body, key).map(|v| v == "true").unwrap_or(false)
    }

    macro_rules! flag {
        ($flag:expr, $key:expr) => {
            if let Some(v) = str_val(body, $key) {
                args.push(format!("--{}", $flag));
                args.push(v);
            }
        };
    }

    flag!("ref_matrix",            "ref_matrix");
    flag!("alt_matrix",            "alt_matrix");
    flag!("barcodes",              "barcodes");
    flag!("num_clusters",          "num_clusters");
    flag!("restarts",              "restarts");
    flag!("seed",                  "seed");
    flag!("threads",               "threads");
    flag!("clustering_method",     "clustering_method");
    flag!("initialization_strategy","initialization_strategy");
    flag!("min_ref",               "min_ref");
    flag!("min_alt",               "min_alt");
    flag!("min_ref_umis",          "min_ref_umis");
    flag!("min_alt_umis",          "min_alt_umis");
    flag!("anneal_steps",          "anneal_steps");
    flag!("anneal_base_em",        "anneal_base_em");
    flag!("anneal_base_khm",       "anneal_base_khm");
    flag!("conv_tol",              "conv_tol");
    flag!("max_iter",              "max_iter");
    flag!("min_cluster_cells",     "min_cluster_cells");
    flag!("pseudocount_alt",       "pseudocount_alt");
    flag!("pseudocount_total",     "pseudocount_total");
    flag!("theta_min",             "theta_min");
    flag!("theta_max",             "theta_max");
    flag!("khm_p",                 "khm_p");
    flag!("plot_dir",              "plot_dir");
    flag!("plots",                 "plots");

    if bool_val(body, "verbose")      { args.push("--verbose".to_string()); }
    if bool_val(body, "souporcell3")  { args.push("--souporcell3".to_string()); args.push("true".to_string()); }

    args
}

// ── Version injection ─────────────────────────────────────────────────────────

fn inject_version(html: &str) -> String {
    // Replace the hardcoded version string in the GUI HTML with the real one
    html.replace("v2.7.0", &format!("v{}", VERSION))
        .replace("\"2.7.0\"", &format!("\"{}\"", VERSION))
        .replace("v2.7", &format!("v{}", VERSION))
        // Inject a JS variable the GUI can read
        .replace("</head>",
            &format!("<script>window.SOUPC_VERSION=\"{}\";window.SOUPC_SERVER=true;</script></head>",
                     VERSION))
}

// ── File listing ──────────────────────────────────────────────────────────────

fn list_files(dir: &str) -> String {
    let path = std::path::Path::new(dir);
    if !path.exists() || !path.is_dir() {
        return format!(r#"{{"error":"directory not found: {}","files":[]}}"#, json_escape(dir));
    }
    let mut entries: Vec<String> = Vec::new();
    if let Ok(read_dir) = std::fs::read_dir(path) {
        for entry in read_dir.flatten() {
            if let Ok(meta) = entry.metadata() {
                if meta.is_file() {
                    let name  = entry.file_name().to_string_lossy().to_string();
                    let size  = meta.len();
                    let fpath = entry.path().to_string_lossy().to_string();
                    entries.push(format!(
                        r#"{{"name":"{}","size":{},"path":"{}"}}"#,
                        json_escape(&name), size, json_escape(&fpath)
                    ));
                }
            }
        }
    }
    entries.sort();
    format!(r#"{{"dir":"{}","files":[{}]}}"#, json_escape(dir), entries.join(","))
}

// ── File serving ──────────────────────────────────────────────────────────────

fn serve_file(stream: &mut TcpStream, path: &str) {
    match std::fs::read(path) {
        Ok(bytes) => {
            let mime = if path.ends_with(".svg")  { "image/svg+xml" }
                  else if path.ends_with(".html") { "text/html" }
                  else if path.ends_with(".tsv") || path.ends_with(".err") { "text/plain" }
                  else { "application/octet-stream" };
            let fname = std::path::Path::new(path)
                .file_name().map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| "file".to_string());
            let header = format!(
                "HTTP/1.1 200 OK\r\n\
                 Content-Type: {}\r\n\
                 Content-Length: {}\r\n\
                 Content-Disposition: attachment; filename=\"{}\"\r\n\
                 Access-Control-Allow-Origin: *\r\n\r\n",
                mime, bytes.len(), fname
            );
            let _ = stream.write_all(header.as_bytes());
            let _ = stream.write_all(&bytes);
        }
        Err(_) => {
            respond_json(stream, 404, r#"{"error":"file not found"}"#);
        }
    }
}

// ── HTTP helpers ──────────────────────────────────────────────────────────────

fn respond_html(stream: &mut TcpStream, code: u16, body: &str) {
    let resp = format!(
        "HTTP/1.1 {} {}\r\n\
         Content-Type: text/html; charset=utf-8\r\n\
         Content-Length: {}\r\n\
         Access-Control-Allow-Origin: *\r\n\r\n{}",
        code, reason(code), body.len(), body
    );
    let _ = stream.write_all(resp.as_bytes());
}

fn respond_json(stream: &mut TcpStream, code: u16, body: &str) {
    let resp = format!(
        "HTTP/1.1 {} {}\r\n\
         Content-Type: application/json\r\n\
         Content-Length: {}\r\n\
         Access-Control-Allow-Origin: *\r\n\r\n{}",
        code, reason(code), body.len(), body
    );
    let _ = stream.write_all(resp.as_bytes());
}

fn reason(code: u16) -> &'static str {
    match code {
        200 => "OK", 204 => "No Content", 404 => "Not Found",
        409 => "Conflict", _ => "Unknown"
    }
}

fn parse_request(raw: &str) -> (String, String, String) {
    let mut lines = raw.lines();
    let first = lines.next().unwrap_or("");
    let parts: Vec<&str> = first.splitn(3, ' ').collect();
    let method   = parts.get(0).unwrap_or(&"GET").to_string();
    // Strip query string (?...) so routes match regardless of browser cache-busting params
    let raw_path = parts.get(1).unwrap_or(&"/");
    let path     = raw_path.splitn(2, '?').next().unwrap_or("/").to_string();
    // Body is after the blank line
    let body = if let Some(idx) = raw.find("\r\n\r\n") {
        raw[idx+4..].to_string()
    } else if let Some(idx) = raw.find("\n\n") {
        raw[idx+2..].to_string()
    } else {
        String::new()
    };
    (method, path, body)
}

fn extract_query_param(path: &str, key: &str) -> Option<String> {
    let q = path.find('?')
        .map(|i| &path[i+1..])
        .unwrap_or("");
    for pair in q.split('&') {
        if let Some((k, v)) = pair.split_once('=') {
            if k == key {
                return Some(percent_decode(v));
            }
        }
    }
    None
}

fn percent_decode(s: &str) -> String {
    let mut out = String::new();
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '%' {
            let h1 = chars.next().unwrap_or('0');
            let h2 = chars.next().unwrap_or('0');
            let hex = format!("{}{}", h1, h2);
            if let Ok(b) = u8::from_str_radix(&hex, 16) {
                out.push(b as char);
            }
        } else if c == '+' {
            out.push(' ');
        } else {
            out.push(c);
        }
    }
    out
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\")
     .replace('"',  "\\\"")
     .replace('\n', "\\n")
     .replace('\r', "\\r")
     .replace('\t', "\\t")
}

fn shell_escape(s: &str) -> String {
    if s.contains(' ') { format!("\"{}\"", s) } else { s.to_string() }
}

// ── Browser launcher ──────────────────────────────────────────────────────────

fn open_browser(url: &str) -> Result<(), ()> {
    #[cfg(target_os = "macos")]
    { std::process::Command::new("open").arg(url).spawn().map(|_| ()).map_err(|_| ()) }
    #[cfg(target_os = "linux")]
    { std::process::Command::new("xdg-open").arg(url).spawn().map(|_| ()).map_err(|_| ()) }
    #[cfg(target_os = "windows")]
    { std::process::Command::new("cmd").args(["/C","start",url]).spawn().map(|_| ()).map_err(|_| ()) }
    #[cfg(not(any(target_os="macos",target_os="linux",target_os="windows")))]
    { Err(()) }
}
