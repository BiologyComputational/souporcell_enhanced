# souporcell v3.0 — Enhanced Modular Rewrite

> **Genotype-based demultiplexing of pooled scRNA-seq experiments**
>
> Original algorithm: Heaton et al., *Nature Methods* 2020 · https://doi.org/10.1038/s41592-020-0820-1
> Original code: [wheaton5/souporcell](https://github.com/wheaton5/souporcell)
> Enhanced rewrite: Jahidul Arafat — clean architecture, structured logging, reproducible seeding,
> diagnostic SVG plots, HTML run report, config profiles, .env file support, fully configurable hyperparameters

---

## Table of Contents

1. [What this binary does](#1-what-this-binary-does)
2. [Repository layout](#2-repository-layout)
3. [Prerequisites and building](#3-prerequisites-and-building)
4. [Input files](#4-input-files)
5. [Running — all four modes](#5-running--all-four-modes)
6. [Complete parameter reference](#6-complete-parameter-reference)
7. [Config profiles (.json)](#7-config-profiles-json)
8. [Environment file (.souporcell.env)](#8-environment-file-souporcellenv)
9. [Usage examples](#9-usage-examples)
10. [Expected outputs](#10-expected-outputs)
11. [Reading the logs](#11-reading-the-logs)
12. [Diagnostic plots](#12-diagnostic-plots)
13. [HTML run report](#13-html-run-report)
14. [What changed from v2.4 (original)](#14-what-changed-from-v24-original)
15. [Troubleshooting](#15-troubleshooting)

---

## 1. What this binary does

### The biological problem

Modern scRNA-seq experiments pool cells from **multiple donors** into a single sequencing run to reduce cost. After sequencing, all cells arrive unlabelled. **Demultiplexing** recovers which cell came from which donor.

`souporcell` solves this **without requiring prior knowledge of the donors' genotypes**. It discovers donor identity purely from the SNP (single nucleotide polymorphism) signals buried in each cell's RNA reads — different donors carry different DNA variants, and even a small number of informative SNP positions is enough to cluster cells into donor groups with high accuracy.

### Where this binary fits

The full souporcell pipeline has 8 stages. This binary is **Stage 6 — the clustering engine**:

```
Stage 1   samtools merge        Merge per-sample BAM files into one pooled BAM
Stage 2   freebayes             Call candidate SNP variants
Stage 3   vcftools              Filter to high-quality biallelic SNPs
Stage 4   vartrix               Count ref/alt allele UMIs per cell per SNP locus
Stage 5   (matrix prep)         Produce ref.mtx and alt.mtx sparse count matrices
──────────────────────────────────────────────────────────────────────────────
Stage 6 ► souporcell (this)     Cluster cells by genotype → clusters_tmp.tsv     ◄
──────────────────────────────────────────────────────────────────────────────
Stage 7   troublet              Detect doublets
Stage 8   consensus.py          Assign donor labels, estimate ambient RNA %
```

### The algorithm in brief

1. **Load** `ref.mtx` and `alt.mtx`, apply locus QC filters (`--min_ref`, `--min_alt`)
2. **Initialise** k cluster centre vectors (random allele fractions θ per locus)
3. **EM/KHM + deterministic annealing** — `anneal_steps` temperature steps (default 9), gradually hardening soft assignments into hard ones
4. **Repeat** for `--restarts` random starts across `--threads` parallel threads; keep the solution with the highest total log-likelihood
5. **Write** one TSV row per cell barcode to stdout with cluster assignment and per-cluster log-probabilities

---

## 2. Repository layout

```
souporcell/
├── Cargo.toml
├── .souporcell.env              ← example project-level env file
├── profiles/                    ← ready-to-use JSON config profiles
│   ├── default.json
│   ├── fast_test.json
│   ├── high_convergence.json    ← fixes low convergence stability
│   └── k4_grch38.json           ← pre-filled for the GSM2560245 test dataset
│
└── src/
    ├── main.rs                  ← composition root — wires all layers together
    │
    ├── domain/                  ① Pure domain (no I/O, no CLI)
    │   ├── types.rs             CellData, ThreadData, ConvergenceStats
    │   └── math.rs              log_sum_exp, normalize_in_log_with_temp, binomial_loss
    │
    ├── config/                  ② Configuration layer
    │   ├── params.rs            Params struct, load_params(), 3-layer merge
    │   ├── params.yml           clap YAML argument definitions
    │   ├── env_loader.rs        .souporcell.env reader
    │   └── config_loader.rs     JSON profile reader
    │
    ├── core/                    ③ Algorithm core (pure computation, no I/O)
    │   ├── em.rs                EM with deterministic annealing — all constants from Params
    │   ├── khm.rs               KHM algorithm — all constants from Params
    │   └── cluster_init.rs      All initialisation strategies
    │
    ├── infra/                   ④ Infrastructure (filesystem + terminal only)
    │   ├── io.rs                Matrix/barcode loading, result writing
    │   └── logger.rs            Structured stderr logging, ANSI colour, timing
    │
    └── analysis/                ⑤ Post-clustering analysis
        ├── plots.rs             6 SVG diagnostic plots
        └── report/              HTML run report (multi-file sub-module)
            ├── mod.rs           ReportData, write_report(), build_html() orchestrator
            ├── css.rs           White-theme CSS
            ├── cards.rs         §2 executive metric cards
            ├── experiment.rs    §3 hypotheses + pipeline walkthrough
            ├── graph.rs         §4 D3 graph JSON builder
            ├── graph_js.rs      §4 D3 <script> block
            ├── graph_template.json
            ├── tables.rs        §5 five results tables
            ├── svgs.rs          §6 balance + confidence SVGs
            ├── plot_embeds.rs   §7 SVG file cards
            ├── interpretation.rs §8 quality commentary cards
            ├── params.rs        §9 parameters table
            ├── methods.rs       §10 methods + references
            └── utils.rs         shared colour constants + chrono_now()
```

### Architecture — dependency rule

```
domain  ←  config  ←  core  ←  infra  ←  analysis  ←  main.rs
  ①          ②         ③         ④           ⑤
```

No inner layer may import from an outer layer. Algorithms never write to disk. I/O never knows about plot generation.

---

## 3. Prerequisites and building

| Tool | Version | Purpose |
|---|---|---|
| Rust + Cargo | ≥ 1.75 (edition 2021) | Compile the binary |
| C linker | system default | Linked automatically |

```bash
# Build optimised binary
cd souporcell_v2.7
cargo build --release

# Verify
./target/release/souporcell --help
```

Build time: 2–4 minutes on first build (downloads all dependencies). Subsequent builds are incremental.

**Apple Silicon (M1/M2/M3):** no special flags needed. Rust targets `aarch64-apple-darwin` natively.

---

## 4. Input files

All three are **required** (or must be set via config/env):

| File | Flag | Description |
|---|---|---|
| `ref.mtx` | `-r` | MatrixMarket sparse matrix — ref allele UMI counts, from `vartrix --ref-matrix` |
| `alt.mtx` | `-a` | MatrixMarket sparse matrix — alt allele UMI counts, from `vartrix --out-matrix` |
| `barcodes.tsv` | `-b` | One barcode per line, plain text or `.gz`. Must match the column order in the matrices. |

---

## 5. Running — all four modes

v2.7 supports four ways to supply parameters. **All four can be combined** — CLI always wins.

### Mode 1 — Pure CLI (original behaviour, unchanged)

```bash
./target/release/souporcell \
  -r /path/to/ref.mtx \
  -a /path/to/alt.mtx \
  -b /path/to/barcodes.tsv \
  -k 4 --threads 8 --restarts 100 --seed 42 \
  --plot_dir ./output \
  > output/clusters_tmp.tsv \
  2> >(tee output/clusters.err)
```

### Mode 2 — JSON config profile (all params in one file)

```bash
# All params — including paths — come from the profile
./target/release/souporcell \
  --config profiles/k4_grch38.json \
  > output/clusters_tmp.tsv \
  2> >(tee output/clusters.err)

# Profile sets defaults; CLI overrides specific values
./target/release/souporcell \
  --config profiles/k4_grch38.json \
  --restarts 500 --seed 99 \
  > output/clusters_tmp.tsv 2>&1
```

### Mode 3 — .env file (auto-discovered)

Place `.souporcell.env` in your working directory (or `~/.souporcell.env`). It is loaded automatically:

```bash
cd /my/project
cp /path/to/souporcell/.souporcell.env .
# Edit SOUPC_* keys to your paths
./target/release/souporcell -k 4   # paths and most params from .env
```

Or specify an explicit env file:

```bash
./target/release/souporcell --env /data/experiment.env -k 4
```

### Mode 4 — Profile + .env + CLI overrides (full stack)

```bash
./target/release/souporcell \
  --env /data/project.env \
  --config profiles/high_convergence.json \
  --restarts 500 \       # CLI overrides the profile's restarts
  > output/clusters_tmp.tsv
```

**Priority chain (highest → lowest):**

```
CLI flags  >  --config JSON  >  .souporcell.env  >  built-in defaults
```

---

## 6. Complete parameter reference

### 6.1 Required (or set via config/env)

| Flag | Short | Env key | Description |
|---|---|---|---|
| `--ref_matrix` | `-r` | `SOUPC_REF_MATRIX` | Ref allele count matrix (`ref.mtx`) |
| `--alt_matrix` | `-a` | `SOUPC_ALT_MATRIX` | Alt allele count matrix (`alt.mtx`) |
| `--barcodes` | `-b` | `SOUPC_BARCODES` | Cell barcodes file (one per line, `.gz` supported) |
| `--num_clusters` | `-k` | `SOUPC_NUM_CLUSTERS` | Number of donor clusters (≥ 2) |

### 6.2 Config / env file selection (v2.7 new)

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--config <path>` | — | — | Load a JSON parameter profile. Any field not present falls through to .env or defaults. |
| `--env <path>` | — | auto | Path to `.souporcell.env`. Auto-discovered in CWD then `~/.souporcell.env` if omitted. |

### 6.3 Core clustering

| Flag | Short | Env key | Default | Description |
|---|---|---|---|---|
| `--restarts` | — | `SOUPC_RESTARTS` | `100` | Random EM/KHM restarts. More = more likely to find the global optimum. |
| `--seed` | — | `SOUPC_SEED` | `4` | RNG seed (u64). Any value 0–2⁶⁴−1. Reproducibility requires the same seed **and** thread count. |
| `--clustering_method` | `-m` | `SOUPC_CLUSTERING_METHOD` | `em` | `em` or `khm`. EM is recommended; KHM is more robust to poor initialisations for large k. |
| `--initialization_strategy` | — | `SOUPC_INIT_STRATEGY` | `random_uniform` | `random_uniform` · `random_cell_assignment` · `kmeans_pp` · `middle_variance` |
| `--souporcell3` | `-s` | `SOUPC_SOUPORCELL3` | `false` | Multi-pass bad-cluster reinitialisation. Only useful for k ≥ 16. |
| `--threads` | `-t` | `SOUPC_THREADS` | `1` | Parallel threads. Set to your physical CPU core count. |

### 6.4 Locus quality filters

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--min_ref` | `SOUPC_MIN_REF` | `4` | Min cells with ≥1 ref read at a locus. |
| `--min_alt` | `SOUPC_MIN_ALT` | `4` | Min cells with ≥1 alt read at a locus. |
| `--min_ref_umis` | `SOUPC_MIN_REF_UMIS` | `0` | Min total ref UMIs across all cells at a locus. |
| `--min_alt_umis` | `SOUPC_MIN_ALT_UMIS` | `0` | Min total alt UMIs across all cells at a locus. |

### 6.5 Annealing schedule (v2.7 new)

All were previously hardcoded constants. Now fully configurable.

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--anneal_steps` | `SOUPC_ANNEAL_STEPS` | `9` | Number of temperature steps in the annealing schedule. T0 = hottest (soft assignments); T(N-1) = 1.0 (standard EM). |
| `--anneal_base_em` | `SOUPC_ANNEAL_BASE_EM` | `20.0` | EM temperature divisor base. `raw_temp = total_alleles / (base × 2^step)`. Higher = softer assignments early on. |
| `--anneal_base_khm` | `SOUPC_ANNEAL_BASE_KHM` | `0.5` | KHM temperature divisor base. KHM runs hotter than EM by default (0.5 vs 20.0). |

**How temperature works:** At each step, for each cell:
```
raw_temp = cell.total_alleles / (anneal_base × 2^temp_step)
temp     = if last_step then 1.0 else max(raw_temp, 1.0)

# Applied as softmax: p̃ᵢ = exp(log_pᵢ / temp) / Z
# temp >> 1 → uniform over clusters (exploratory)
# temp = 1  → standard EM softmax (exploitative)
```

### 6.6 EM/KHM convergence (v2.7 new)

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--conv_tol` | `SOUPC_CONV_TOL` | `0.01` | Stop iterating when `|Δ loss| < conv_tol × n_cells`. Decrease for tighter convergence; increase for speed. |
| `--max_iter` | `SOUPC_MAX_ITER` | `1000` | Hard iteration cap per temperature step. Guards against non-convergent runs. |
| `--min_cluster_cells` | `SOUPC_MIN_CLUSTER_CELLS` | `200` | Min cells for a cluster to be counted as "populated" in verbose convergence logs. |

### 6.7 M-step pseudocounts (v2.7 new)

Laplace smoothing applied in the M-step to keep allele fraction estimates away from 0 and 1:

```
θ_c[locus] = (Σ prob_i × alt_i + pseudocount_alt)
           / (Σ prob_i × total_i + pseudocount_total)
```

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--pseudocount_alt` | `SOUPC_PSEUDOCOUNT_ALT` | `1.0` | Added to the alt-allele numerator. |
| `--pseudocount_total` | `SOUPC_PSEUDOCOUNT_TOTAL` | `2.0` | Added to the total-allele denominator. Must be ≥ pseudocount_alt. |

### 6.8 Allele-frequency bounds (v2.7 new)

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--theta_min` | `SOUPC_THETA_MIN` | `0.01` | Lower bound for θ after M-step. Keeps log-likelihood finite. |
| `--theta_max` | `SOUPC_THETA_MAX` | `0.99` | Upper bound for θ after M-step. Must be > theta_min. |

### 6.9 KHM-specific (v2.7 new)

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--khm_p` | `SOUPC_KHM_P` | `25.0` | KHM harmonic-mean exponent. Higher = more sharply peaked membership around the nearest cluster centre. |

### 6.10 Optional priors

| Flag | Short | Env key | Description |
|---|---|---|---|
| `--known_genotypes` | `-g` | `SOUPC_KNOWN_GENOTYPES` | Population VCF/BCF with GT fields for each donor. Skips random init if provided. |
| `--known_genotypes_sample_names` | — | `SOUPC_KNOWN_GENOTYPES_SAMPLE_NAMES` | Space-separated (CLI) or comma-separated (env) sample names from the VCF. |

### 6.11 Logging

| Flag | Short | Env key | Default | Description |
|---|---|---|---|---|
| `--verbose` | `-v` | `SOUPC_VERBOSE` | off | Emit a `CONV` line for every single EM/KHM iteration. Warning: with restarts=100 this produces tens of thousands of lines. |
| `--progress_interval` | — | `SOUPC_PROGRESS_INTERVAL` | `500` | Print a progress line every N cells. Set 0 to suppress. |

### 6.12 Diagnostics and output

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--plot_dir` | `SOUPC_PLOT_DIR` | disabled | Directory for SVG plots and HTML report. Zero overhead when omitted. |
| `--plots` | `SOUPC_PLOTS` | `all` | Which plots to generate: `all`, `none`, or comma-separated: `restart,convergence,annealing,balance,confidence,threads` |
| `--cleanup` | `SOUPC_CLEANUP` | off | Report file count and disk usage of plot_dir after a run. |

### 6.13 Pre-flight / dry-run (v2.7 new)

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--dry_run` | `SOUPC_DRY_RUN` | off | Print the full run plan and wait for `[Y/n]` approval **before reading or writing anything**. Shows inputs, pipeline stages, all resolved hyperparameters, expected outputs, and a sample interpretation guide. Type `n` to abort cleanly. |
| `--dry_run_yes` | `SOUPC_DRY_RUN_YES` | off | Requires `--dry_run`. Print the plan and proceed automatically without prompting. For CI pipelines. |

**What the plan shows:**

- **Run identity** — resolved run ID, algorithm, k, restarts, seed, threads, config/env sources
- **Input files** — all three paths, QC filters with plain-English explanation
- **Pipeline stages** — all 10 stages in order with per-stage description populated from your actual params
- **Algorithm hyperparameters** — full annealing schedule (with a worked temperature example), convergence, pseudocounts, theta bounds, KHM exponent
- **Expected outputs** — sample `clusters_tmp.tsv` rows with column annotations, SVG plot list, HTML report
- **Expected interpretation** — 7 quality metrics with expected values and warning thresholds
- **Quality checklist** — the 10 automated checks the HTML report will evaluate
- **Estimated runtime** — work units breakdown across threads

---

## 7. Config profiles (.json)

A profile is a JSON file that sets any subset of parameters. Fields not present fall through to env or defaults. The `--config` flag points to the file.

### Schema

```json
{
  "profile":     "my_experiment",
  "description": "Human-readable description — informational only",
  "params": {
    "ref_matrix":              "/path/to/ref.mtx",
    "alt_matrix":              "/path/to/alt.mtx",
    "barcodes":                "/path/to/barcodes.tsv",
    "num_clusters":            4,
    "restarts":                200,
    "seed":                    42,
    "threads":                 8,
    "clustering_method":       "em",
    "initialization_strategy": "random_uniform",
    "souporcell3":             false,
    "min_ref":                 4,
    "min_alt":                 4,
    "min_ref_umis":            0,
    "min_alt_umis":            0,
    "anneal_steps":            9,
    "anneal_base_em":          20.0,
    "anneal_base_khm":         0.5,
    "conv_tol":                0.01,
    "max_iter":                1000,
    "min_cluster_cells":       200,
    "pseudocount_alt":         1.0,
    "pseudocount_total":       2.0,
    "theta_min":               0.01,
    "theta_max":               0.99,
    "khm_p":                   25.0,
    "verbose":                 false,
    "progress_interval":       500,
    "plot_dir":                "./output",
    "plots":                   "all",
    "cleanup":                 false
  }
}
```

All fields are optional. Omit any you don't want to fix.

### Bundled profiles

| Profile | When to use | Key differences from defaults |
|---|---|---|
| `profiles/default.json` | General-purpose starting point | All defaults |
| `profiles/fast_test.json` | Smoke test — not for real analysis | `anneal_steps=3`, `restarts=10`, `conv_tol=0.05`, `max_iter=100`, `plots=none` |
| `profiles/high_convergence.json` | Sparse SNP coverage, convergence < 70% | `restarts=500`, `min_ref/alt=2`, `conv_tol=0.005`, `max_iter=2000`, `min_cluster_cells=100` |
| `profiles/k4_grch38.json` | Pre-filled for the GSM2560245 test dataset | All paths + `k=4`, `restarts=200`, `threads=8` |

---

## 8. Environment file (.souporcell.env)

A shell-style `KEY=VALUE` file loaded automatically from the current directory or `~/.souporcell.env`.

### Rules

- All keys must start with `SOUPC_`
- Lines starting with `#` are comments
- Values may be optionally quoted with `"` or `'` (quotes are stripped)
- No `export`, no `$VARIABLE` expansion

### Full example

```ini
# .souporcell.env — project-level defaults

# ── I/O paths ──────────────────────────────────────────────────────────
SOUPC_REF_MATRIX=/Volumes/extra_space/demux_test/ref.mtx
SOUPC_ALT_MATRIX=/Volumes/extra_space/demux_test/alt.mtx
SOUPC_BARCODES=~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv

# ── Core clustering ────────────────────────────────────────────────────
SOUPC_NUM_CLUSTERS=4
SOUPC_RESTARTS=100
SOUPC_SEED=42
SOUPC_THREADS=8
SOUPC_CLUSTERING_METHOD=em
SOUPC_INIT_STRATEGY=random_uniform

# ── Locus QC ───────────────────────────────────────────────────────────
SOUPC_MIN_ALT=4
SOUPC_MIN_REF=4
SOUPC_MIN_ALT_UMIS=0
SOUPC_MIN_REF_UMIS=0

# ── Annealing schedule ─────────────────────────────────────────────────
SOUPC_ANNEAL_STEPS=9
SOUPC_ANNEAL_BASE_EM=20.0
SOUPC_ANNEAL_BASE_KHM=0.5

# ── Convergence ────────────────────────────────────────────────────────
SOUPC_CONV_TOL=0.01
SOUPC_MAX_ITER=1000
SOUPC_MIN_CLUSTER_CELLS=200

# ── M-step pseudocounts ────────────────────────────────────────────────
SOUPC_PSEUDOCOUNT_ALT=1.0
SOUPC_PSEUDOCOUNT_TOTAL=2.0

# ── Allele-frequency bounds ────────────────────────────────────────────
SOUPC_THETA_MIN=0.01
SOUPC_THETA_MAX=0.99

# ── KHM ────────────────────────────────────────────────────────────────
SOUPC_KHM_P=25.0

# ── Output ─────────────────────────────────────────────────────────────
SOUPC_PLOT_DIR=~/Setups/ambientRNA_research/souporcell/playground/output
SOUPC_PLOTS=all
SOUPC_CLEANUP=false

# ── Logging ────────────────────────────────────────────────────────────
SOUPC_VERBOSE=false
SOUPC_PROGRESS_INTERVAL=500
```

---

## 9. Usage examples

### Standard run — 4 donors, reproducible

```bash
BIN=./target/release/souporcell
OUT=~/souporcell/playground/output
mkdir -p $OUT

$BIN \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --restarts 100 --seed 42 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

### Using a config profile (all params from file)

```bash
$BIN --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

### Profile + CLI override (restarts and seed overridden on CLI)

```bash
$BIN --config profiles/k4_grch38.json \
  --restarts 500 --seed 99 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Fix low convergence (< 70%) — high_convergence profile

```bash
$BIN --config profiles/high_convergence.json \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Tune the annealing schedule directly

```bash
# More temperature steps, tighter convergence, higher iteration cap
$BIN -r ref.mtx -a alt.mtx -b barcodes.tsv -k 4 \
  --anneal_steps 12 \
  --anneal_base_em 10.0 \
  --conv_tol 0.005 \
  --max_iter 2000 \
  --threads 8 --restarts 200 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Adjust M-step pseudocounts (more aggressive smoothing)

```bash
$BIN -r ref.mtx -a alt.mtx -b barcodes.tsv -k 4 \
  --pseudocount_alt 2.0 \
  --pseudocount_total 4.0 \
  --theta_min 0.005 \
  --theta_max 0.995 \
  --threads 8 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### KHM with custom harmonic exponent

```bash
$BIN -r ref.mtx -a alt.mtx -b barcodes.tsv -k 4 \
  --clustering_method khm \
  --khm_p 50.0 \
  --anneal_base_khm 1.0 \
  --threads 8 --restarts 200 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### High k — 8 donors with stricter filters

```bash
$BIN -r ref.mtx -a alt.mtx -b barcodes.tsv -k 8 \
  --threads 8 --restarts 200 --seed 42 \
  --min_ref 10 --min_alt 10 \
  --min_ref_umis 5 --min_alt_umis 5 \
  --anneal_steps 12 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Preview run plan before executing (--dry_run)

```bash
# Show full plan, wait for Y/n approval, then run
$BIN --dry_run \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --seed 42 --restarts 100 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

```bash
# Show plan with profile, then prompt
$BIN --dry_run --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv 2>&1
```

```bash
# CI mode — print plan, proceed automatically (no stdin needed)
$BIN --dry_run --dry_run_yes --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Quick smoke test

```bash
$BIN --config profiles/fast_test.json \
  -r ref.mtx -a alt.mtx -b barcodes.tsv -k 4 \
  > /dev/null
```

### Verbose convergence debugging

```bash
$BIN -r ref.mtx -a alt.mtx -b barcodes.tsv -k 4 \
  --threads 1 --restarts 5 --verbose \
  > /dev/null \
  2> detailed.err

# Extract loss series per iteration:
grep '^CONV' detailed.err | cut -f7

# Watch convergence per temperature step:
grep '^EM' detailed.err
```

### All 6 SVG diagnostic plots

```bash
$BIN -r ref.mtx -a alt.mtx -b barcodes.tsv -k 4 \
  --threads 8 --restarts 100 --seed 42 --verbose \
  --plot_dir $OUT --plots all \
  > $OUT/clusters_tmp.tsv 2>&1

ls $OUT/
# restart_landscape.svg  convergence_curves.svg  annealing_profile.svg
# cluster_balance.svg    posterior_confidence.svg  thread_efficiency.svg
# souporcell_report.html
```

### Known donor genotypes

```bash
$BIN -r ref.mtx -a alt.mtx -b barcodes.tsv -k 3 \
  --known_genotypes donors.vcf.gz \
  --known_genotypes_sample_names DONOR_A DONOR_B DONOR_C \
  --threads 4 \
  > $OUT/clusters_tmp.tsv 2>&1
```

---

## 10. Expected outputs

### stdout → `clusters_tmp.tsv`

One row per cell barcode:

```
AAACATTGAGCTAC-1    2    -312.44    -298.11    -287.56    -301.88
AAACATTGATCAGC-1    0    -289.12    -401.55    -398.70    -411.22
```

| Column | Content |
|---|---|
| 1 | Cell barcode |
| 2 | Best-fit cluster (0-based, 0 to k-1) |
| 3..k+2 | Log-probability per cluster — least negative = best match |

### SVG plots → `<plot_dir>/` (when `--plot_dir` set)

| File | What it shows |
|---|---|
| `restart_landscape.svg` | Final log-likelihood per restart, sorted worst→best. Healthy: ≥70% of bars at global-best level. |
| `convergence_curves.svg` | Per-iteration loss curves. Requires `--verbose`. |
| `annealing_profile.svg` | Mean ± SD loss at each temperature step T0–T(N-1). |
| `cluster_balance.svg` | Cells per cluster vs expected equal split. |
| `posterior_confidence.svg` | Histogram of `best_lp − second_best_lp`. Red zone = ambiguous (< 10 log-units). |
| `thread_efficiency.svg` | Restarts and best loss per thread. |

### HTML report → `<plot_dir>/souporcell_report.html` (when `--plot_dir` set)

Self-contained, no internet required. Contains 11 sections including metric cards, interactive D3.js pipeline graph, 5 results tables, quality interpretation cards, full parameters table, and algorithmic methods with citations.

---

## 11. Reading the logs

### stderr structure

```
[PARAMS] ──────────────────────────────────────────────────────────
[PARAMS] Sources    : CLI > config(profiles/k4_grch38.json) > env(.souporcell.env)
[PARAMS] ref_matrix : /Volumes/extra_space/demux_test/ref.mtx
[PARAMS] k=4  restarts=200  seed=42  threads=8
[PARAMS] Annealing  : steps=9  base_em=20.0  base_khm=0.5
[PARAMS] Convergence: tol=0.01  max_iter=1000  min_cluster_cells=200
[PARAMS] M-step     : pseudocount alt=1.0 total=2.0  theta=[0.01, 0.99]
[PARAMS] ──────────────────────────────────────────────────────────

EM  thread=0  epoch=0  temp_step=0/8  iters=8   loss=-45221.3  delta=0.0042
EM  thread=0  epoch=0  temp_step=1/8  iters=5   loss=-32118.7  delta=0.0088
...
RESTART  0  0  -13402.12  -13402.12
RESTART  0  1  -13388.44  -13388.44
...

cluster  0    912 cells ( 25.1%)  ################
cluster  1    891 cells ( 24.5%)  ###############
cluster  2    876 cells ( 24.1%)  ###############
cluster  3    960 cells ( 26.4%)  #################
```

### Greppable line formats

| Prefix | Always/Verbose | Fields |
|---|---|---|
| `[PARAMS]` | Always | Parameter source and resolved values |
| `EM` / `KHM` | Always | `thread epoch temp_step iters loss delta` — one per annealing step |
| `RESTART` | Always | `thread iteration loss best_so_far` |
| `CONV` | `--verbose` only | `method thread epoch iter temp_step loss delta clusters_populated` |

```bash
# Track global-best improvement across restarts:
grep '^RESTART' clusters.err | awk '{print $5}' | uniq

# Extract annealing profile for thread 0, restart 0:
grep '^EM' clusters.err | grep 'thread=0' | grep 'epoch=0'

# Check convergence was reached (iters < max_iter):
grep '^EM' clusters.err | awk -F'iters=' '{print $2}' | cut -d$'\t' -f1 | sort -n | tail -5
```

---

## 12. Diagnostic plots

Six SVG plots are generated from data accumulated during the run — no extra passes through the matrices.

| Plot | Key indicator | Healthy | Warning |
|---|---|---|---|
| `restart_landscape.svg` | % restarts at global optimum | ≥ 70% bars at dashed line | < 40%: increase `--restarts` or try `high_convergence` profile |
| `convergence_curves.svg` | Curve smoothness, spread at final iter | Tight band, rapid flattening | Widely spread final values: increase `--restarts` |
| `annealing_profile.svg` | Loss improvement T0→T(N-1) | Steady increase, narrowing SD | Plateau at T4–T5: try increasing `--anneal_steps` |
| `cluster_balance.svg` | Bar heights vs dashed ideal | All bars within ±30% | One bar < 5%: k may be too large |
| `posterior_confidence.svg` | Bimodal distribution | Clear high-confidence peak | Unimodal flat: insufficient informative loci |
| `thread_efficiency.svg` | Red line flatness | All threads reach same best-loss | One thread much worse: per-thread restart budget too low |

---

## 13. HTML run report

Generated at `<plot_dir>/souporcell_report.html`. Single self-contained file — no internet required.

| Section | Content |
|---|---|
| §1 Cover | Run ID, method, convergence badge, timestamp |
| §2 Metric cards | 8 live statistics: cells, loci, k, best log-prob, convergence %, imbalance, confidence, runtime |
| §3 Experiment design | Hypotheses, input data table, 11-step pipeline walkthrough, expected vs observed outcomes table |
| §4 Pipeline graph | Interactive D3.js directed graph — click nodes for detail panels, step-by-step mode, fullscreen |
| §5 Results tables | 5 tables: per-cluster stats, restart landscape, annealing trace, thread efficiency, locus QC |
| §6 Quality charts | Inline cluster balance SVG + posterior confidence histogram |
| §7 SVG plots | Linked diagnostic plots |
| §8 Interpretation | 6 auto-evaluated quality cards (convergence, balance, confidence, sparsity, locus QC, log-prob) |
| §9 Parameters | Complete record of every resolved parameter with source provenance |
| §10 Methods | Algorithm prose, quality threshold table, 5 primary citations, software dependencies |

---

## 14. What changed from v2.4 (original)

The original binary shipped as a single `main.rs` of 856 lines — 27 functions with no module boundaries, no timing, bare unstructured `eprintln!` calls, a developer uncertainty comment printed verbatim in production (`"is this a problem?"`), 6 `assert!(false, ...)` stubs, a `seed` typed as `u8` limited to 256 values, and every algorithmic constant hardcoded with no way to tune them without recompiling.

v2.7 is a complete rewrite across three stages (v2.5 → v2.6 → v2.7).

---

### Architecture (v2.6)

| # | What changed | v2.4 | v2.7 |
|---|---|---|---|
| 1 | File structure | 1 file, 856 lines | 5 dependency-ordered layers, 20+ source files |
| 2 | Module boundaries | None — algorithm, I/O, logging all mixed | `domain / config / core / infra / analysis` with enforced import rules |
| 3 | Algorithm isolation | `em()` and `khm()` called `eprintln!` directly | Core layer has no I/O — all output via `infra/logger.rs` |
| 4 | M-step helpers | Duplicated between EM and KHM | `pub(crate)` in `em.rs`, reused by `khm.rs` |
| 5 | Report module | Single `report.rs` (1290 lines) | `analysis/report/` — 13 files, one responsibility each |

### Reproducibility and correctness (v2.5)

| # | What changed | v2.4 | v2.7 |
|---|---|---|---|
| 6 | Seed type | `u8` — only 256 possible seeds; values > 255 silently overflow to 0 | `u64` — full 64-bit RNG space |
| 7 | Seed spreading | `[params.seed; 32]` — repeats one byte 32 times | `seed.to_le_bytes()` distributed across the 32-byte array |
| 8 | Unimplemented stubs | 6× `assert!(false, "X not yet implemented")` | `panic!("... Use --initialization_strategy random_uniform")` |
| 9 | Known genotypes assertion | `assert!(false, "currently requiring...")` at wrong call site | Positive assertion at correct call site with actionable message |

### Logging and observability (v2.5)

| # | What changed | v2.4 | v2.7 |
|---|---|---|---|
| 10 | Zero-count message | `"N loaded 0 counts, is this a problem?"` printed verbatim | Computed as percentage; warns only if below expected 80% sparsity |
| 11 | CONV lines | 7-column positional, inconsistent names | `CONV method thread epoch iter temp_step loss delta clusters_populated` — 9 named tab-separated fields |
| 12 | RESTART lines | `"thread N iteration N done with X, best so far Y"` | `RESTART thread iter loss best_so_far` — parseable |
| 13 | Per-temp-step summary | Buried in convergence stream | `EM thread=N epoch=N temp_step=N/8 iters=N loss=N delta=N` — always emitted |
| 14 | Final summary | `"best total log probability = X"` only | Per-cluster cell count, percentage, ASCII bar chart, underpopulation warning |
| 15 | Timing | None | `Logger { start: Instant }` — every section prints `[+Xm Ys]` elapsed |
| 16 | ANSI colour | None | `is_tty()` gated — colours in terminal, plain text in pipes/files |
| 17 | Debug eprintlns | 9 commented-out lines left in production code | Removed; replaced by structured always-on or `--verbose`-gated equivalents |

### New CLI flags (v2.5–v2.7)

| # | Flag | Version | Purpose |
|---|---|---|---|
| 18 | `--verbose` / `-v` | v2.5 | CONV line per iteration. Off by default. |
| 19 | `--progress_interval` | v2.5 | Progress every N cells. Default 500; 0 = suppress. |
| 20 | `--plot_dir` / `--plots` | v2.6 | 6 SVG diagnostic plots. Zero overhead when omitted. |
| 21 | `--cleanup` | v2.6 | Report disk usage of plot_dir after run. |
| 22 | `--config` | v2.7 | Load a JSON parameter profile. |
| 23 | `--env` | v2.7 | Explicit path to a `.souporcell.env` file. |
| 24 | `--dry_run` | v2.7 | Print full run plan and prompt `[Y/n]` before executing. Nothing is read or written until approved. |
| 25 | `--dry_run_yes` | v2.7 | Requires `--dry_run`. Auto-approve — for CI pipelines. Also set via `SOUPC_DRY_RUN_YES=true`. |
| 24 | `--anneal_steps` | v2.7 | Number of annealing temperature steps (was hardcoded 9). |
| 25 | `--anneal_base_em` | v2.7 | EM temperature divisor base (was hardcoded 20.0). |
| 26 | `--anneal_base_khm` | v2.7 | KHM temperature divisor base (was hardcoded 0.5). |
| 27 | `--conv_tol` | v2.7 | Convergence tolerance (was hardcoded 0.01). |
| 28 | `--max_iter` | v2.7 | Max iterations per temp step (was hardcoded 1000). |
| 29 | `--min_cluster_cells` | v2.7 | Populated-cluster threshold (was hardcoded 200). |
| 30 | `--pseudocount_alt` | v2.7 | M-step alt pseudocount (was hardcoded 1.0). |
| 31 | `--pseudocount_total` | v2.7 | M-step total pseudocount (was hardcoded 2.0). |
| 32 | `--theta_min` | v2.7 | θ lower bound (was hardcoded 0.01). |
| 33 | `--theta_max` | v2.7 | θ upper bound (was hardcoded 0.99). |
| 34 | `--khm_p` | v2.7 | KHM harmonic exponent (was `const KHM_P: f32 = 25.0`). |

### Previously hardcoded constants — now all configurable (v2.7)

Every algorithmic constant that was hardcoded in the original is now a first-class parameter flowing through the full priority chain (CLI > config > env > default):

| Constant | Was | Flag | Default | Where used |
|---|---|---|---|---|
| Annealing steps | `9` (literal) | `--anneal_steps` | `9` | `em.rs`, `khm.rs` |
| EM temperature base | `20.0` (literal) | `--anneal_base_em` | `20.0` | `em.rs` |
| KHM temperature base | `0.5` (literal) | `--anneal_base_khm` | `0.5` | `khm.rs` |
| Convergence tolerance | `0.01` (literal) | `--conv_tol` | `0.01` | `em.rs`, `khm.rs` |
| Max iterations | `1000` (literal) | `--max_iter` | `1000` | `em.rs`, `khm.rs` |
| Min cluster cells | `200` (literal) | `--min_cluster_cells` | `200` | `em.rs`, `khm.rs` |
| Alt pseudocount | `1.0` (literal) | `--pseudocount_alt` | `1.0` | `em.rs`, `khm.rs` |
| Total pseudocount | `2.0` (literal) | `--pseudocount_total` | `2.0` | `em.rs`, `khm.rs` |
| θ lower bound | `0.01` (literal) | `--theta_min` | `0.01` | `em.rs`, `khm.rs`, `cluster_init.rs` |
| θ upper bound | `0.99` (literal) | `--theta_max` | `0.99` | `em.rs`, `khm.rs`, `cluster_init.rs` |
| KHM exponent p | `const KHM_P = 25.0` | `--khm_p` | `25.0` | `khm.rs` |

### New outputs (v2.6–v2.7)

| # | Output | Description |
|---|---|---|
| 35 | `restart_landscape.svg` | Final log-likelihood per restart, sorted worst→best, colour-coded by thread |
| 36 | `convergence_curves.svg` | Per-iteration loss curve for every restart. Best restart highlighted. Requires `--verbose`. |
| 37 | `annealing_profile.svg` | Mean ± SD of final loss at each of the N temperature steps |
| 38 | `cluster_balance.svg` | Cells per cluster vs expected equal split |
| 39 | `posterior_confidence.svg` | Histogram of `best_lp − second_best_lp` per cell, red-shaded ambiguous zone |
| 40 | `thread_efficiency.svg` | Restart count per thread (bars) and best log-likelihood per thread (line) |
| 41 | `souporcell_report.html` | Self-contained HTML report with D3.js pipeline graph, 5 results tables, quality cards, parameters |

### New data structures (v2.5–v2.7)

| # | What | v2.4 | v2.7 |
|---|---|---|---|
| 42 | `load_cell_data` return | 5-tuple — positional, no field names | `LoadResult` struct with named fields |
| 43 | `em()` / `khm()` return | `(f32, Vec<Vec<f32>>)` | `(f32, Vec<Vec<f32>>, ConvergenceStats, Vec<ts_losses>)` |
| 44 | `Logger` | None | `Logger { start: Instant }` — timed section headers |
| 45 | `PlotData` | None | Accumulates diagnostics during run via `Arc<Mutex<PlotData>>` |
| 46 | `ReportData` | None | Bundles run statistics for HTML report |
| 47 | `ConfigProfile` | None | Holds parsed JSON profile fields with typed accessors |
| 48 | `EnvMap` | None | Wraps loaded env map with typed `parse<T>()`, `bool()`, `str()` accessors |

---

## 15. Troubleshooting

**`cargo build` fails with linker error on macOS**
```bash
xcode-select --install
```

**All cells assigned to one cluster**

Check the log for `Loci passing QC: 0`. If zero, the locus filters are too strict. Try `--min_ref 2 --min_alt 2`.

**Convergence stability < 70%**

This is common when SNP coverage is sparse (< 5% of loci passing QC). Use the bundled profile:
```bash
./target/release/souporcell --config profiles/high_convergence.json \
  -r ref.mtx -a alt.mtx -b barcodes.tsv -k 4
```
Or tune directly: `--restarts 500 --min_ref 2 --min_alt 2 --conv_tol 0.005 --max_iter 2000`

**`iters=1000` appearing in EM lines**

The hard iteration cap was hit. Either the data has noisy loci (try stricter `--min_ref`/`--min_alt`) or the convergence tolerance is too tight (try `--conv_tol 0.02`).

**Cluster sizes very unequal**

Try more restarts: `--restarts 200`. Also try `--initialization_strategy random_cell_assignment`. If k > 4, consider whether k is larger than the true number of donors.

**`--seed` does not give identical results across runs**

Results are reproducible given the same seed **and** thread count. Changing `--threads` distributes restarts differently and will produce a different (equally valid) best solution.

**Profile not found**

```bash
# Check the path is relative to where you run the binary from:
./target/release/souporcell --config ./profiles/k4_grch38.json ...
```

**`.env` file not loaded**

The env file must be named `.souporcell.env` (leading dot) and placed in the current working directory, or in `~/.souporcell.env`. Use `--env /explicit/path.env` to override discovery.

**HTML report is blank or missing plots**

The report is written to `<plot_dir>/souporcell_report.html`. Make sure `--plot_dir` is set and the directory is writable. The SVG plot links are relative paths — keep the SVGs and the HTML in the same directory.

---

*souporcell v2.7 · Original algorithm: Heaton et al., Nature Methods 2020 · Enhanced rewrite: Jahidul Arafat*
