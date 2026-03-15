Based on your paths from previous sessions:

```bash

# 1. Unzip and build
cd ~/Setups/ambientRNA_research/souporcell
unzip ~/Downloads/souporcell_v2.7.zip
mv souporcell_v2.6_new souporcell_v2.7
cd souporcell_v2.7
cargo build --release 2>&1 | tee build.log
```

```bash
# gui mode
# Launch GUI — opens browser automatically at http://localhost:7979
./target/release/souporcell --gui

# Custom port
./target/release/souporcell --gui --gui_port 8080
```

```bash
# 2. Set your output dir
BIN=souporcell
OUT=/Users/jarotball/Setups/ambientRNA_research/souporcell_pipeline_components/playground/output
mkdir -p $OUT
```

Then pick **one** of these run modes:

```bash
# Mode A — pure CLI (same as before)
$BIN \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b /Users/jarotball/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --seed 42 --restarts 100 --verbose \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

```bash
# Mode B — profile only (all params from file, including paths)
$BIN --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

```bash
# Mode C — fix the convergence failure (was 21%, needs > 70%)
$BIN --config profiles/high_convergence.json \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

```bash
# Mode D — tune annealing directly on CLI
$BIN \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --seed 42 \
  --restarts 500 \
  --anneal_steps 12 \
  --anneal_base_em 10.0 \
  --conv_tol 0.005 \
  --max_iter 2000 \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

```bash
# Mode E — .env file (auto-discovered, no paths on CLI)
cp profiles/.souporcell.env $OUT/.souporcell.env
# edit SOUPC_* paths inside it if needed, then:
cd $OUT
$BIN > clusters_tmp.tsv 2> >(tee clusters.err)
```

**To check results after any run:**

```bash
# Convergence stability — how many restarts found the global optimum
grep '^RESTART' $OUT/clusters.err | awk '{print $5}' | sort -u

# Final cluster breakdown
tail -20 $OUT/clusters.err | grep 'cluster'

# Open the HTML report
open $OUT/souporcell_report.html
```


### User consent edition
**108 KB · 49 files** (one new: `src/infra/preflight.rs`).

---

**Build and run commands:**

```bash
# Unzip and build
cd ~/Setups/ambientRNA_research/souporcell
unzip ~/Downloads/souporcell_v2.7.zip
mv souporcell_v2.6_new souporcell_v2.7
cd souporcell_v2.7
cargo build --release 2>&1 | tee build.log
```

```bash
BIN=~/Setups/ambientRNA_research/souporcell/souporcell_v2.7/target/release/souporcell
OUT=~/Setups/ambientRNA_research/souporcell/playground/output
mkdir -p $OUT
```

**The new dry_run flow — shows the plan first, asks Y/n, then runs:**

```bash
$BIN --dry_run \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --seed 42 --restarts 100 --verbose \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

You'll see the full plan printed to stderr — inputs, all 10 pipeline stages, every resolved hyperparameter with an explanation, a worked temperature example, expected output format with sample rows, 7 interpretation benchmarks, and the 10 quality checklist items. Then:

```
  Proceed with this run? [Y/n]:
```

Type `Y` (or just Enter) to run. Type `n` to abort cleanly — nothing was read or written.

**Using a profile with dry_run:**

```bash
$BIN --dry_run --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

**CI mode — plan logged, no prompt:**

```bash
$BIN --dry_run --dry_run_yes --config profiles/k4_grch38.json \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```

**Dry run without actually executing (just see the plan, always abort):**

```bash
$BIN --dry_run \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --restarts 100 \
  2>&1 | less    # pipe to less to read the full plan comfortably, then type n
```

**Turn off dry_run (original behaviour, no change):**

```bash
# Just omit --dry_run entirely — runs immediately as before
$BIN \
  -r /Volumes/extra_space/demux_test/ref.mtx \
  -a /Volumes/extra_space/demux_test/alt.mtx \
  -b ~/Setups/ambientRNA_research/souporcell_test/GSM2560245_barcodes.tsv \
  -k 4 --threads 8 --seed 42 --restarts 100 --verbose \
  --plot_dir $OUT \
  > $OUT/clusters_tmp.tsv \
  2> >(tee $OUT/clusters.err)
```