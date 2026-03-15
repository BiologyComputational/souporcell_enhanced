# souporcell v3.0.0 — Enhanced Modular Rewrite

<p align="center">
  <strong>Genotype-based demultiplexing of pooled single-cell RNA-seq experiments</strong><br/>
  Zero prior knowledge of donor genotypes required
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-3.0.0-8B5CF6" alt="version"/>
  <img src="https://img.shields.io/badge/rust-≥1.75-orange" alt="rust"/>
  <img src="https://img.shields.io/badge/edition-2021-blue" alt="edition"/>
  <img src="https://img.shields.io/badge/algorithm-EM%20%2B%20KHM-teal" alt="algorithm"/>
</p>

---

> **Original algorithm:** Heaton et al., *Nature Methods* 2020 · [doi:10.1038/s41592-020-0820-1](https://doi.org/10.1038/s41592-020-0820-1)
>
> **Original code:** [wheaton5/souporcell](https://github.com/wheaton5/souporcell)
>
> **Enhanced rewrite (v3.0.0):** Jahidul Arafat — clean 5-layer architecture · structured logging · reproducible seeding · 6 diagnostic SVG plots · interactive HTML run report · config profiles · `.env` file support · fully configurable hyperparameters · web GUI · pre-flight dry-run · algorithm explorer

---

## Table of Contents

1. [What souporcell does](#1-what-souporcell-does)
2. [Where this binary fits in the pipeline](#2-where-this-binary-fits-in-the-pipeline)
3. [Repository layout](#3-repository-layout)
4. [Prerequisites and building](#4-prerequisites-and-building)
5. [Quick start](#5-quick-start)
6. [All run modes](#6-all-run-modes)
7. [Complete parameter reference](#7-complete-parameter-reference)
8. [Config profiles (.json)](#8-config-profiles-json)
9. [Environment file (.souporcell.env)](#9-environment-file-souporcellenv)
10. [All usage examples](#10-all-usage-examples)
11. [Expected outputs](#11-expected-outputs)
12. [Reading the logs](#12-reading-the-logs)
13. [Diagnostic plots](#13-diagnostic-plots)
14. [HTML run report](#14-html-run-report)
15. [Web GUI](#15-web-gui)
16. [Algorithm explorer](#16-algorithm-explorer)
17. [What changed from v2.4 (original)](#17-what-changed-from-v24-original)
18. [Troubleshooting](#18-troubleshooting)
19. [Citation](#19-citation)

---

## 1. What souporcell does

### The biological problem

Modern scRNA-seq experiments pool cells from **multiple donors** into a single sequencing run to reduce cost and eliminate batch effects. After sequencing, all cells arrive unlabelled — there is no record of which cell came from which person. **Demultiplexing** is the process of recovering that donor identity.

`souporcell` solves this **without requiring prior knowledge of the donors' genotypes**. It discovers donor identity purely from the single-nucleotide polymorphism (SNP) signals embedded in each cell's RNA reads. Different donors carry different DNA variants, and even a small fraction of covered SNP positions is enough to cluster cells into donor groups with high accuracy.

### The algorithm in brief

```
1.  LOAD         ref.mtx + alt.mtx → sparse CellData structures
                 Apply locus QC filters (min_ref, min_alt, ...)
                 Pre-compute ln C(total, alt) for each cell-locus pair

2.  INITIALISE   K cluster centre vectors θ ∈ ℝᴸ (one allele frequency per locus)
                 Strategy: random_uniform / random_cell_assignment / known_genotypes

3.  CLUSTER      For each of R restarts (parallelised over T threads):
                   For each of S annealing temperature steps:
                     EM or KHM iteration until convergence
                     → soft cell-to-cluster assignments p̃_nk
                     → updated θ via weighted M-step with pseudocounts
                 Select restart with highest total log-likelihood

4.  OUTPUT       clusters_tmp.tsv → stdout
                 6 SVG diagnostic plots + HTML report → plot_dir
```

---

## 2. Where this binary fits in the pipeline

souporcell is **Stage 6** of the full demultiplexing pipeline.

```
Stage 1   samtools merge       Merge per-sample BAM files into one pooled BAM
Stage 2   freebayes            Call candidate SNP variants from the merged BAM
Stage 3   vcftools             Filter to high-quality biallelic SNPs only
Stage 4   vartrix              Count ref/alt allele UMIs per (cell × SNP locus)
Stage 5   (matrix prep)        Produce ref.mtx and alt.mtx sparse count matrices
──────────────────────────────────────────────────────────────────────────────────
Stage 6 ► souporcell (this)    Cluster cells by genotype → clusters_tmp.tsv    ◄
──────────────────────────────────────────────────────────────────────────────────
Stage 7   troublet             Detect doublets using posterior log-probabilities
Stage 8   consensus.py         Assign donor labels, estimate ambient RNA fraction
```

**Inputs:** `ref.mtx`, `alt.mtx`, `barcodes.tsv`
**Outputs:** `clusters_tmp.tsv` (stdout), SVG plots, HTML report (when `--plot_dir` is set)

---

## 3. Repository layout

```
souporcell_v3.0.0/
├── Cargo.toml
├── .souporcell.env                  ← project-level env file (edit your paths here)
├── profiles/
│   ├── default.json
│   ├── fast_test.json               ← smoke test only
│   ├── high_convergence.json        ← fixes convergence stability < 70%
│   └── k4_grch38.json               ← pre-filled for the GSM2560245 test dataset
└── src/
    ├── main.rs                      ← composition root
    ├── domain/                      ① Pure domain — no I/O, no CLI
    │   ├── types.rs                   CellData, ThreadData, ConvergenceStats
    │   └── math.rs                    log_sum_exp, normalize_in_log_with_temp,
    │                                  binomial_loss_with_min_index
    ├── config/                      ② Configuration — 3-layer param resolution
    │   ├── params.rs                  Params struct, load_params(), resolve macros
    │   ├── params.yml                 clap YAML argument definitions (34+ flags)
    │   ├── env_loader.rs              .souporcell.env reader with tilde expansion
    │   └── config_loader.rs           JSON profile reader
    ├── core/                        ③ Algorithm core — pure computation, no I/O
    │   ├── em.rs                      EM with deterministic annealing
    │   ├── khm.rs                     K-Harmonic Means algorithm
    │   └── cluster_init.rs            All initialisation strategies
    ├── infra/                       ④ Infrastructure — filesystem + terminal
    │   ├── io.rs                      Matrix/barcode loading, result writing
    │   ├── logger.rs                  Structured stderr logging, ANSI colour, timing
    │   ├── preflight.rs               --dry_run plan printer + [Y/n] gate
    │   ├── gui_server.rs              --gui HTTP server (pure std::net, zero deps)
    │   ├── souporcell_gui.html        Web GUI (embedded at compile time)
    │   └── souporcell_algo_analysis.html  Algorithm Explorer (embedded at compile time)
    └── analysis/                    ⑤ Post-clustering analysis
        ├── plots.rs                   6 SVG diagnostic plots
        └── report/                    HTML run report (13-file sub-module)
            ├── mod.rs  css.rs  cards.rs  experiment.rs
            ├── graph.rs  graph_js.rs  graph_template.json
            ├── tables.rs  svgs.rs  plot_embeds.rs
            └── interpretation.rs  params.rs  methods.rs  utils.rs
```

### Dependency rule

```
domain  ←  config  ←  core  ←  infra  ←  analysis  ←  main.rs
  ①          ②         ③         ④           ⑤
```

No inner layer may import from an outer layer. Algorithms never write to disk. I/O never knows about plot generation.

---

## 4. Prerequisites and building

| Tool | Version | Notes |
|---|---|---|
| Rust + Cargo | ≥ 1.75 (edition 2021) | Install via [rustup.rs](https://rustup.rs) |
| C linker | system default | Linked automatically |
| macOS Xcode tools | any | Only needed on macOS if linker is missing |

```bash
# Unzip and build
cd ~/Setups/ambientRNA_research/souporcell
unzip ~/Downloads/souporcell_v3.0.0.zip
mv souporcell_v2.6_new souporcell_v3.0.0
cd souporcell_v3.0.0

# Build optimised release binary (2–4 min first time)
cargo build --release 2>&1 | tee build.log

# Verify
./target/release/souporcell --version
./target/release/souporcell --help
```

```bash
# Convenience alias — add to ~/.zshrc or ~/.bashrc
export SOUPC_BIN=~/Setups/ambientRNA_research/souporcell/souporcell_v3.0.0/target/release/souporcell
alias souporcell=$SOUPC_BIN
```

> **Apple Silicon (M1/M2/M3/M4):** no special flags needed. Rust targets `aarch64-apple-darwin` natively.
> **macOS linker error:** run `xcode-select --install` then rebuild.

---

## 5. Quick start

```bash
BIN=./target/release/souporcell
OUT=~/souporcell_output
mkdir -p $OUT

# Run
$BIN \
  -r /path/to/ref.mtx \
  -a /path/to/alt.mtx \
  -b /path/to/barcodes.tsv \
  -k 4 --threads 8 --seed 42 --restarts 100 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)

# Check results
tail -5 $OUT/clusters_tmp.tsv
grep 'cluster' $OUT/clusters.err | tail -10
open $OUT/souporcell_report.html
```

Or launch the **web GUI** (no command knowledge required):

```bash
$BIN --gui
# Opens http://localhost:7979 automatically
```

---

## 6. All run modes

### Priority chain (highest → lowest)

```
CLI flags  >  --config JSON profile  >  .souporcell.env  >  built-in defaults
```

### Mode 1 — Pure CLI

```bash
./target/release/souporcell \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --seed 42 --restarts 100 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

### Mode 2 — JSON config profile

```bash
# All params from profile (paths + algorithm settings)
./target/release/souporcell \
  --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)

# Profile sets base; CLI overrides specific values
./target/release/souporcell \
  --config profiles/k4_grch38.json \
  --restarts 500 --seed 99 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Mode 3 — .env file

```bash
# Auto-discovery: CWD then ~/.souporcell.env
cd /my/experiment
cp /path/to/souporcell_v3.0.0/.souporcell.env .
# Edit SOUPC_* keys to your paths
./target/release/souporcell -k 4

# Explicit path
./target/release/souporcell --env /data/project.env -k 4
```

### Mode 4 — Full combined stack

```bash
./target/release/souporcell \
  --env /data/project.env \
  --config profiles/high_convergence.json \
  --restarts 500 --seed 123 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Mode 5 — Web GUI

```bash
./target/release/souporcell --gui              # port 7979
./target/release/souporcell --gui --gui_port 8080
```

Navigate to [http://localhost:7979](http://localhost:7979). Four phases: Configure → Pre-flight → Execute → Outputs. Press `Ctrl+C` to stop the server.

### Mode 6 — Dry run pre-flight

```bash
# Show plan → prompt [Y/n] → run if Y
./target/release/souporcell --dry_run \
  -r $REF -a $ALT -b $BAR \
  -k 4 --threads 8 --seed 42 --restarts 100 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)

# CI mode — auto-approve
./target/release/souporcell --dry_run --dry_run_yes \
  --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv 2>&1

# View plan only (pipe to less, press q, type n at prompt)
./target/release/souporcell --dry_run \
  -r $REF -a $ALT -b $BAR -k 4 \
  2>&1 | less
```

---

## 7. Complete parameter reference

### 7.1 Required inputs (or set via config / env)

| Flag | Short | Env key | Description |
|---|---|---|---|
| `--ref_matrix <path>` | `-r` | `SOUPC_REF_MATRIX` | Ref allele UMI count matrix (`ref.mtx`) from `vartrix --ref-matrix` |
| `--alt_matrix <path>` | `-a` | `SOUPC_ALT_MATRIX` | Alt allele UMI count matrix (`alt.mtx`) from `vartrix --out-matrix` |
| `--barcodes <path>` | `-b` | `SOUPC_BARCODES` | Cell barcodes — one per line, plain text or `.gz`. Column order must match matrices. |
| `--num_clusters <int>` | `-k` | `SOUPC_NUM_CLUSTERS` | Number of donor clusters (≥ 2) |

> Tilde (`~`) is expanded in all path fields regardless of source (CLI / JSON / .env).

### 7.2 Configuration and environment file selection

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--config <path>` | — | — | Load a JSON parameter profile |
| `--env <path>` | — | auto | Explicit path to `.souporcell.env`. Auto-discovered in CWD then `~/.souporcell.env` if omitted. |

### 7.3 Core clustering

| Flag | Short | Env key | Default | Description |
|---|---|---|---|---|
| `--restarts <int>` | — | `SOUPC_RESTARTS` | `100` | Independent EM/KHM restarts. More = better chance of finding the global optimum. |
| `--seed <u64>` | — | `SOUPC_SEED` | `4` | Master RNG seed. Same seed + same threads = identical results. |
| `--threads <int>` | `-t` | `SOUPC_THREADS` | `1` | Parallel threads. Set to physical core count. |
| `--clustering_method <str>` | `-m` | `SOUPC_CLUSTERING_METHOD` | `em` | `em` (recommended) or `khm` (more robust for K ≥ 8) |
| `--initialization_strategy <str>` | — | `SOUPC_INIT_STRATEGY` | `random_uniform` | `random_uniform` · `random_cell_assignment` · `kmeans_pp` · `middle_variance` |
| `--souporcell3` | `-s` | `SOUPC_SOUPORCELL3` | off | Multi-pass bad-cluster reinitialisation. Only useful for K ≥ 16. |

### 7.4 Locus quality filters

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--min_ref <int>` | `SOUPC_MIN_REF` | `4` | Min cells with ≥ 1 ref read at a locus |
| `--min_alt <int>` | `SOUPC_MIN_ALT` | `4` | Min cells with ≥ 1 alt read at a locus |
| `--min_ref_umis <int>` | `SOUPC_MIN_REF_UMIS` | `0` | Min total ref UMIs across all cells at a locus |
| `--min_alt_umis <int>` | `SOUPC_MIN_ALT_UMIS` | `0` | Min total alt UMIs across all cells at a locus |

### 7.5 Annealing schedule *(all were hardcoded in v2.4)*

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--anneal_steps <int>` | `SOUPC_ANNEAL_STEPS` | `9` | Temperature steps. T0 = hottest; T(S-1) = 1.0 (standard EM). |
| `--anneal_base_em <float>` | `SOUPC_ANNEAL_BASE_EM` | `20.0` | EM temperature divisor base. Higher = more exploratory early assignments. |
| `--anneal_base_khm <float>` | `SOUPC_ANNEAL_BASE_KHM` | `0.5` | KHM temperature divisor base. KHM runs much hotter than EM by default. |

**Temperature formula:**
```
raw_temp = cell.total_alleles / (anneal_base × 2^step)
temp     = if last_step then 1.0 else max(raw_temp, 1.0)

# Example: anneal_base_em=20, cell with 100 total alleles
step 0: temp = 100/(20×1)   = 5.0   ← hot, exploratory
step 4: temp = 100/(20×16)  = 0.31  → clamped to 1.0
step 8: temp = 1.0 (forced) ← cold, standard EM
```

### 7.6 Convergence *(all were hardcoded in v2.4)*

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--conv_tol <float>` | `SOUPC_CONV_TOL` | `0.01` | Stop when `|Δloss| < conv_tol × N`. Smaller = tighter but slower. |
| `--max_iter <int>` | `SOUPC_MAX_ITER` | `1000` | Hard iteration cap per temperature step |
| `--min_cluster_cells <int>` | `SOUPC_MIN_CLUSTER_CELLS` | `200` | Min cells for a cluster to be counted as "populated" in logs |

### 7.7 M-step pseudocounts *(were hardcoded in v2.4)*

```
θ_kl = (Σ_n p̃_nk · alt_nl  + pseudocount_alt)
     / (Σ_n p̃_nk · total_nl + pseudocount_total)
```

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--pseudocount_alt <float>` | `SOUPC_PSEUDOCOUNT_ALT` | `1.0` | Alt-allele numerator addend |
| `--pseudocount_total <float>` | `SOUPC_PSEUDOCOUNT_TOTAL` | `2.0` | Total-allele denominator addend (≥ pseudocount_alt) |

### 7.8 Allele-frequency bounds *(were hardcoded in v2.4)*

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--theta_min <float>` | `SOUPC_THETA_MIN` | `0.01` | Lower clamp on θ after M-step. Keeps ln(θ) finite. |
| `--theta_max <float>` | `SOUPC_THETA_MAX` | `0.99` | Upper clamp on θ after M-step |

### 7.9 KHM-specific *(was hardcoded in v2.4)*

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--khm_p <float>` | `SOUPC_KHM_P` | `25.0` | Harmonic-mean exponent. Higher = sharper membership around nearest cluster. |

```
w_nk = (1/d_nk)^p  /  Σ_j(1/d_nj)^p    where d_nk = −log P(cell_n | cluster_k)
```

### 7.10 Optional donor genotype priors

| Flag | Short | Env key | Description |
|---|---|---|---|
| `--known_genotypes <path>` | `-g` | `SOUPC_KNOWN_GENOTYPES` | VCF/BCF with GT fields. Initialises θ from genotype data — reduces restarts needed dramatically. |
| `--known_genotypes_sample_names <names>` | — | `SOUPC_KNOWN_GENOTYPES_SAMPLE_NAMES` | Space-separated sample names matching VCF header |

### 7.11 Logging

| Flag | Short | Env key | Default | Description |
|---|---|---|---|---|
| `--verbose` | `-v` | `SOUPC_VERBOSE` | off | Emit a `CONV` line for every EM/KHM iteration |
| `--progress_interval <int>` | — | `SOUPC_PROGRESS_INTERVAL` | `500` | Progress line every N cells. Set to `0` to suppress. |

### 7.12 Output and diagnostics

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--plot_dir <path>` | `SOUPC_PLOT_DIR` | disabled | Directory for SVG plots and HTML report. Zero overhead when omitted. |
| `--plots <str>` | `SOUPC_PLOTS` | `all` | `all`, `none`, or comma-separated: `restart,convergence,annealing,balance,confidence,threads` |
| `--cleanup` | `SOUPC_CLEANUP` | off | Report file count and disk usage of plot_dir after the run |

### 7.13 Pre-flight / dry-run

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--dry_run` | `SOUPC_DRY_RUN` | off | Print full run plan to stderr, wait for `[Y/n]` before reading or writing anything |
| `--dry_run_yes` | `SOUPC_DRY_RUN_YES` | off | Requires `--dry_run`. Auto-approve without prompting. For CI. |

### 7.14 Web GUI

| Flag | Env key | Default | Description |
|---|---|---|---|
| `--gui` | — | — | Launch the built-in HTTP server and open the web GUI |
| `--gui_port <int>` | — | `7979` | Port for the GUI server |

---

## 8. Config profiles (.json)

### Schema (all fields optional)

```json
{
  "profile":     "my_experiment",
  "description": "Human-readable description",
  "params": {
    "ref_matrix":              "/path/to/ref.mtx",
    "alt_matrix":              "/path/to/alt.mtx",
    "barcodes":                "/path/to/barcodes.tsv",
    "num_clusters":            4,
    "restarts":                100,
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

### Bundled profiles

| Profile | When to use | Key differences |
|---|---|---|
| `profiles/default.json` | General-purpose starting point | All built-in defaults |
| `profiles/fast_test.json` | Smoke test — verify build works | `anneal_steps=3`, `restarts=10`, `conv_tol=0.05`, `max_iter=100`, `plots=none` |
| `profiles/high_convergence.json` | Sparse coverage, stability < 70% | `restarts=500`, `min_ref=2`, `min_alt=2`, `conv_tol=0.005`, `max_iter=2000`, `min_cluster_cells=100` |
| `profiles/k4_grch38.json` | GSM2560245 test dataset | All paths pre-filled, `k=4`, `restarts=200`, `threads=8` |

---

## 9. Environment file (.souporcell.env)

Auto-discovered in CWD then `~/.souporcell.env`. All keys must start with `SOUPC_`.

```ini
# .souporcell.env — copy into your working directory and edit paths

# ── Required inputs ────────────────────────────────────────────────────
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

# ── Locus QC filters ───────────────────────────────────────────────────
SOUPC_MIN_REF=4
SOUPC_MIN_ALT=4
SOUPC_MIN_REF_UMIS=0
SOUPC_MIN_ALT_UMIS=0

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
SOUPC_PLOT_DIR=~/souporcell_output
SOUPC_PLOTS=all
SOUPC_CLEANUP=false

# ── Logging ────────────────────────────────────────────────────────────
SOUPC_VERBOSE=false
SOUPC_PROGRESS_INTERVAL=500
```

---

## 10. All usage examples

```bash
# Set up once
BIN=~/Setups/ambientRNA_research/souporcell/souporcell_v3.0.0/target/release/souporcell
REF=/Volumes/extra_space/demux_test/ref.mtx
ALT=/Volumes/extra_space/demux_test/alt.mtx
BAR=~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv
OUT=~/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output
mkdir -p $OUT
```

### Standard run — 4 donors, reproducible

```bash
$BIN -r $REF -a $ALT -b $BAR \
  -k 4 --threads 8 --seed 42 --restarts 100 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

### Profile only (all params from file, including paths)

```bash
$BIN --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

### Profile + CLI override

```bash
$BIN --config profiles/k4_grch38.json \
  --restarts 500 --seed 99 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Fix low convergence stability (< 70%)

```bash
# Quick: use bundled profile
$BIN --config profiles/high_convergence.json \
  -r $REF -a $ALT -b $BAR -k 4 --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv 2>&1

# Manual tuning
$BIN -r $REF -a $ALT -b $BAR -k 4 \
  --restarts 500 --min_ref 2 --min_alt 2 \
  --conv_tol 0.005 --max_iter 2000 --threads 8 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Tune the annealing schedule

```bash
$BIN -r $REF -a $ALT -b $BAR -k 4 \
  --anneal_steps 12 \
  --anneal_base_em 10.0 \
  --conv_tol 0.005 --max_iter 2000 \
  --restarts 200 --threads 8 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### KHM with custom harmonic exponent

```bash
$BIN -r $REF -a $ALT -b $BAR -k 4 \
  --clustering_method khm \
  --khm_p 50.0 --anneal_base_khm 1.0 \
  --threads 8 --restarts 200 \
  > $OUT/clusters_tmp.tsv 2>&1
```

### High K — 8 donors, stricter QC

```bash
$BIN -r $REF -a $ALT -b $BAR -k 8 \
  --threads 8 --restarts 200 --seed 42 \
  --min_ref 10 --min_alt 10 \
  --min_ref_umis 5 --min_alt_umis 5 \
  --anneal_steps 12 --clustering_method khm \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Adjust pseudocounts (stronger smoothing)

```bash
$BIN -r $REF -a $ALT -b $BAR -k 4 \
  --pseudocount_alt 2.0 --pseudocount_total 4.0 \
  --theta_min 0.005 --theta_max 0.995 \
  --threads 8 > $OUT/clusters_tmp.tsv 2>&1
```

### Known donor genotypes (VCF)

```bash
$BIN -r $REF -a $ALT -b $BAR -k 3 \
  --known_genotypes donors.vcf.gz \
  --known_genotypes_sample_names DONOR_A DONOR_B DONOR_C \
  --threads 4 > $OUT/clusters_tmp.tsv 2>&1
```

### Dry run — preview plan, then prompt

```bash
$BIN --dry_run \
  -r $REF -a $ALT -b $BAR \
  -k 4 --threads 8 --seed 42 --restarts 100 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

### Dry run with profile

```bash
$BIN --dry_run --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv 2>&1
```

### CI mode (auto-approve)

```bash
$BIN --dry_run --dry_run_yes --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv 2>&1
```

### Just view the plan (always abort)

```bash
$BIN --dry_run -r $REF -a $ALT -b $BAR -k 4 \
  2>&1 | less
```

### Launch the web GUI

```bash
$BIN --gui
$BIN --gui --gui_port 8080
```

### Quick smoke test

```bash
$BIN --config profiles/fast_test.json \
  -r $REF -a $ALT -b $BAR -k 4 > /dev/null
```

### Verbose convergence debugging

```bash
$BIN -r $REF -a $ALT -b $BAR -k 4 \
  --threads 1 --restarts 5 --verbose \
  > /dev/null 2> detailed.err

grep '^CONV' detailed.err | cut -f7        # loss series per iteration
grep '^EM' detailed.err                    # per-temperature-step summary
grep '^RESTART' detailed.err | awk '{print $4, $5}'  # global-best progression
```

### Check results after any run

```bash
head -5 $OUT/clusters_tmp.tsv
grep 'cluster' $OUT/clusters.err | tail -10
grep '^RESTART' $OUT/clusters.err | awk '{print $5}' | sort -u   # convergence stability
ls -lh $OUT/
open $OUT/souporcell_report.html
```

---

## 11. Expected outputs

### stdout → `clusters_tmp.tsv`

```
AAACATTGAGCTAC-1    2    -312.44    -298.11    -287.56    -301.88
AAACATTGATCAGC-1    0    -289.12    -401.55    -398.70    -411.22
```

| Column | Content |
|---|---|
| 1 | Cell barcode |
| 2 | Best-fit cluster (0-based, 0 to K-1) |
| 3 … K+2 | Log-probability per cluster — least negative = best match |

The difference between the best and second-best log-probability is the **confidence score** used by doublet detectors.

### SVG diagnostic plots → `<plot_dir>/`

| File | What it shows |
|---|---|
| `restart_landscape.svg` | Final log-likelihood per restart sorted worst→best, coloured by thread |
| `convergence_curves.svg` | Per-iteration loss curve per restart. **Requires `--verbose`.** |
| `annealing_profile.svg` | Mean ± SD of final loss at each temperature step T0–T(S-1) |
| `cluster_balance.svg` | Cells per cluster vs expected equal split |
| `posterior_confidence.svg` | Histogram of best_lp − second_best_lp per cell |
| `thread_efficiency.svg` | Restart count and best log-likelihood per thread |

### HTML report → `<plot_dir>/souporcell_report.html`

Self-contained, no internet required. See [§14](#14-html-run-report).

---

## 12. Reading the logs

### stderr structure

```
[PARAMS] ──────────────────────────────────────────────────────────────────────
[PARAMS] Sources        : CLI > config(profiles/k4_grch38.json) > env(.souporcell.env)
[PARAMS] ref_matrix     : /Volumes/extra_space/demux_test/ref.mtx
[PARAMS] k=4  restarts=200  seed=42  threads=8
[PARAMS] ── Annealing ──────────────────
[PARAMS] anneal_steps=9  base_em=20  base_khm=0.5
[PARAMS] ── Convergence ─────────────────
[PARAMS] conv_tol=0.01  max_iter=1000  min_cluster_cells=200
[PARAMS] ──────────────────────────────────────────────────────────────────────

══ DATA LOADING [+0.001s]
  │ barcodes loaded ........... 3639
  │ Loci in matrix:  80712
  │ Loci passing QC: 8219
  │ Cells loaded:    3639
  ✔  Data loading complete  [+0.390s]

══ CLUSTERING [+0.391s]
EM    thread=0  epoch=0  temp_step=0/8  iters=9   loss=-255987.3  delta=17.7
EM    thread=0  epoch=0  temp_step=1/8  iters=1   loss=-255981.9  delta=5.3
RESTART    0    0    -255661.98    -255661.98
  ✔  Run 1 complete — best log-prob: -255653.19  [+5.983s]

  cluster  0    827 cells  ( 22.7%)  ###########
  cluster  1    843 cells  ( 23.2%)  ###########
  cluster  2    927 cells  ( 25.5%)  ############
  cluster  3   1042 cells  ( 28.6%)  ##############
```

### Greppable line formats

| Prefix | When | Tab-separated fields |
|---|---|---|
| `[PARAMS]` | Always | Parameter source and resolved values |
| `EM` / `KHM` | Always (one per temp step) | `thread epoch temp_step iters loss delta` |
| `RESTART` | Always | `thread iteration loss best_so_far` |
| `CONV` | `--verbose` only | `method thread epoch iter temp_step loss delta n_populated` |

### Useful grep recipes

```bash
# Convergence stability (unique best_so_far values — ideally 1–2)
grep '^RESTART' $OUT/clusters.err | awk '{print $5}' | sort -u

# Global-best improvement trace across all restarts
grep '^RESTART' $OUT/clusters.err | awk '{print $4, $5}'

# Check for max_iter hits (1000 = cap was reached = problem)
grep '^EM' $OUT/clusters.err | awk -F'iters=' '{print $2}' | cut -f1 | sort -n | tail -5

# Per-temp-step profile for thread 0, restart 0
grep '^EM' $OUT/clusters.err | awk '$2==0 && $3==0'

# Extract CONV loss series for plotting
grep '^CONV' $OUT/clusters.err | awk '{print $6}' > loss_series.txt
```

---

## 13. Diagnostic plots

| Plot | Key indicator | Healthy | Warning |
|---|---|---|---|
| `restart_landscape.svg` | % restarts at global optimum | ≥ 70% at dashed line | < 40%: increase `--restarts` or use `high_convergence` profile |
| `convergence_curves.svg` | Spread at final iteration | Tight cluster, rapid flattening | Wide spread: increase `--restarts` or `--max_iter` |
| `annealing_profile.svg` | Loss improvement T0→T(S-1) | Steady increase, narrowing SD | Plateau at T4–T5: try `--anneal_steps 12` |
| `cluster_balance.svg` | Bar heights vs ideal | All within ±30% | One bar < 5%: K may exceed true donor count |
| `posterior_confidence.svg` | Distribution shape | Clear high-confidence peak | Flat unimodal: insufficient informative loci |
| `thread_efficiency.svg` | Red line flatness | All threads reach same best | One thread much worse: per-thread restart budget too low |

---

## 14. HTML run report

Generated at `<plot_dir>/souporcell_report.html`. Single self-contained file, no internet required.

| Section | Content |
|---|---|
| §1 Cover | Run ID, algorithm, convergence badge, timestamp |
| §2 Metric cards | 8 statistics: cells, loci, K, best log-prob, convergence %, imbalance, confidence, runtime |
| §3 Experiment design | Hypotheses, input data table, 11-step pipeline walkthrough |
| §4 Pipeline graph | Interactive D3.js directed graph — click nodes for detail panels |
| §5 Results tables | Per-cluster stats, restart landscape, annealing trace, thread efficiency, locus QC |
| §6 Quality charts | Cluster balance SVG + posterior confidence histogram |
| §7 SVG plots | Linked diagnostic plots with interpretation notes |
| §8 Interpretation | 6 auto-evaluated quality cards |
| §9 Parameters | Complete record of every resolved parameter with provenance |
| §10 Methods | Algorithm prose, quality thresholds, 5 citations, software dependencies |

---

## 15. Web GUI

```bash
./target/release/souporcell --gui              # http://localhost:7979
./target/release/souporcell --gui --gui_port 8080
```

| Phase | What you do |
|---|---|
| 1 — Configure | Select a profile, edit all 34 parameters, review the live CLI preview |
| 2 — Pre-flight | Review the structured dry-run plan |
| 3 — Execute | Click Run. Watch live streaming log with stage checkmarks and restart counter |
| 4 — Outputs | Browse and download all output files |

**Server endpoints:**

| Endpoint | Description |
|---|---|
| `GET /` | Main GUI |
| `GET /algo-analysis` | Algorithm Explorer |
| `GET /api/version` | `{"version":"3.0.0","binary":"..."}` |
| `GET /api/status` | `{"status":"ready"\|"running"}` |
| `POST /api/run` | Start a run, streams SSE events |
| `GET /api/files?dir=...` | List output files |
| `GET /files/<path>` | Download a file |

---

## 16. Algorithm explorer

Available at `http://localhost:7979/algo-analysis` (GUI mode) or standalone at `src/infra/souporcell_algo_analysis.html`.

Designed for readers with zero biology background. Covers every algorithmic stage interactively.

| Feature | Content |
|---|---|
| 9 pipeline stages | Each with WHY / WHAT / HOW / WHEN steps, navigated with Next/Prev or arrow keys |
| Live visualisations | EM convergence simulation · sparse matrix density · annealing temperature curve · Amdahl speedup chart · KHM membership sharpness |
| Complexity table | 42 operations with time/space complexity, calls per run, runtime dominance %. Filterable + searchable. |
| Glossary | 31 plain-English definitions — Biology · Genomics · Statistics · Computer Science · Mathematics — each with an everyday analogy |

---

## 17. What changed from v2.4 (original)

### Architecture

| # | What | v2.4 | v3.0.0 |
|---|---|---|---|
| 1 | Source files | 1 file, 856 lines | 5 dependency-ordered layers, 25+ source files |
| 2 | Module boundaries | None | `domain / config / core / infra / analysis` with enforced import rules |
| 3 | Algorithm isolation | `em()` called `eprintln!` directly | Core layer has zero I/O |
| 4 | Report module | Single 1290-line `report.rs` | `analysis/report/` — 13 focused files |

### Reproducibility

| # | What | v2.4 | v3.0.0 |
|---|---|---|---|
| 5 | Seed type | `u8` — 256 possible values, silently wraps | `u64` — full 64-bit RNG space |
| 6 | Seed spreading | `[params.seed; 32]` — repeats one byte | `seed.to_le_bytes()` distributed across 32-byte array |
| 7 | Unimplemented stubs | 6× `assert!(false, "X not yet implemented")` | Proper `panic!` with actionable messages |

### Logging and observability

| # | What | v2.4 | v3.0.0 |
|---|---|---|---|
| 8  | Zero-count message | `"N loaded 0 counts, is this a problem?"` | Computed as sparsity %; warns only if below 80% |
| 9  | CONV lines | 7-column positional | 9 named tab-separated fields |
| 10 | RESTART lines | Unstructured sentence | `RESTART thread iter loss best_so_far` — parseable |
| 11 | Timing | None | Every section prints `[+Xm Ys]` elapsed |
| 12 | ANSI colour | None | `is_tty()` gated — colour in terminal, plain in pipes |
| 13 | Debug eprintlns | 9 left in production code | Removed; replaced by structured equivalents |

### New CLI flags (v2.5 → v3.0.0)

| Flag | Added | Previously |
|---|---|---|
| `--verbose` / `-v` | v2.5 | No per-iteration logging |
| `--progress_interval` | v2.5 | No progress reporting |
| `--plot_dir` / `--plots` | v2.6 | No diagnostic plots |
| `--config` | v2.7 | No profile support |
| `--env` | v2.7 | No .env support |
| `--dry_run` / `--dry_run_yes` | v2.7 | No pre-flight check |
| `--anneal_steps` | v2.7 | Hardcoded `9` |
| `--anneal_base_em` | v2.7 | Hardcoded `20.0` |
| `--anneal_base_khm` | v2.7 | Hardcoded `0.5` |
| `--conv_tol` | v2.7 | Hardcoded `0.01` |
| `--max_iter` | v2.7 | Hardcoded `1000` |
| `--min_cluster_cells` | v2.7 | Hardcoded `200` |
| `--pseudocount_alt` | v2.7 | Hardcoded `1.0` |
| `--pseudocount_total` | v2.7 | Hardcoded `2.0` |
| `--theta_min` | v2.7 | Hardcoded `0.01` |
| `--theta_max` | v2.7 | Hardcoded `0.99` |
| `--khm_p` | v2.7 | `const KHM_P: f32 = 25.0` |
| `--gui` / `--gui_port` | v3.0 | No GUI |

### New outputs (v2.6 → v3.0.0)

| Output | Added | Description |
|---|---|---|
| `restart_landscape.svg` | v2.6 | Final log-likelihood per restart, sorted, coloured by thread |
| `convergence_curves.svg` | v2.6 | Per-iteration loss curves. Requires `--verbose`. |
| `annealing_profile.svg` | v2.6 | Mean ± SD loss per temperature step |
| `cluster_balance.svg` | v2.6 | Cells per cluster vs equal-split ideal |
| `posterior_confidence.svg` | v2.6 | Histogram of confidence scores |
| `thread_efficiency.svg` | v2.6 | Restart count and best loss per thread |
| `souporcell_report.html` | v2.6 | Self-contained interactive HTML report |
| Web GUI (`--gui`) | v3.0 | Browser-based configuration, execution, and file download |
| Algorithm Explorer | v3.0 | Interactive educational deep-dive + complexity table + glossary |

---

## 18. Troubleshooting

### `cargo build` fails with linker error on macOS
```bash
xcode-select --install
cargo build --release
```

### All cells assigned to one cluster
Check for `Loci passing QC: 0` — QC filters are too strict:
```bash
$BIN -r $REF -a $ALT -b $BAR -k 4 --min_ref 2 --min_alt 2
```

### Convergence stability < 70%
```bash
$BIN --config profiles/high_convergence.json -r $REF -a $ALT -b $BAR -k 4
# or manually:
$BIN -r $REF -a $ALT -b $BAR -k 4 \
  --restarts 500 --min_ref 2 --min_alt 2 --conv_tol 0.005 --max_iter 2000
```

### `iters=1000` in EM lines (max-iter cap hit)
```bash
--conv_tol 0.02        # looser convergence
--max_iter 3000        # higher cap
--min_ref 8 --min_alt 8  # stricter locus filters to reduce noise
```

### Cluster sizes very unequal
```bash
--restarts 200
--initialization_strategy random_cell_assignment
# Also verify K matches the actual number of donors
```

### `--seed` does not give identical results
Same seed + same `--threads` = identical results. Changing `--threads` changes per-thread sub-seeds.

### Profile not found
```bash
# Paths are relative to CWD
cd ~/souporcell_v3.0.0
./target/release/souporcell --config ./profiles/k4_grch38.json ...
# Or use absolute path
$BIN --config ~/souporcell_v3.0.0/profiles/k4_grch38.json ...
```

### `.env` file not loaded
Must be named `.souporcell.env` (leading dot) in CWD or `~/.souporcell.env`:
```bash
$BIN --env /explicit/path/.souporcell.env -k 4
```

### HTML report blank or missing plots
- Ensure `--plot_dir` is set and writable
- Keep SVG files and HTML in the same directory (SVG links are relative)
- Open with a modern browser — not from inside a zip file

### GUI shows `{"error":"not found"}` for `/algo-analysis`
Navigate directly to `http://localhost:7979/algo-analysis` — the route is registered and the query-string stripping bug is fixed in v3.0.0.

---

## 19. Citation

If you use souporcell in published research, please cite:

> Heaton, H., Talman, A.M., Smith, G. *et al.* Souporcell: robust clustering of single-cell RNA-seq data by genotype without reference genotypes. *Nat Methods* **17**, 615–620 (2020). https://doi.org/10.1038/s41592-020-0820-1

---

*souporcell v3.0.0 · Original algorithm: Heaton et al., Nature Methods 2020 · Enhanced rewrite: Jahidul Arafat*
