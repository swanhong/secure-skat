# Secure RVAS

Run privacy-preserving Burden/SKAT with MVP (Party 1), AoU (Party 2), and a helper without genomic inputs (Party 0).
For a public-data local test, see the [1000 Genomes workflow](python/datasets/onekg/README.md).

Run all commands from the repository root using the same Python environment.

0. Install dependencies once (Linux x86-64; requires Go, Python 3.11+ with pip in a writable conda/venv environment, and R).

```bash
bash scripts/setup/install.sh
source scripts/setup/env.sh
go mod vendor
```

On macOS, install PLINK 2 and R::SKAT separately, then run `python3 -m pip install -r requirements.txt` in a conda/venv environment.
Load `scripts/setup/env.sh` in each new terminal.

1. Configure each site's inputs and output directory.

Edit `config/mvp/configPrepare.toml` and `config/aou/configPrepare.toml`: set `run_dir`, genotype/PGEN prefix, gene panel, annotation, phenotype, covariate, and ancestry paths. Match `masks` to the annotation columns/values; `samples_per_cohort = 0` uses all eligible samples.
Set the helper's `run_dir` in `config/cp0/configPrepare.toml`. Normally, leave `configLocal.toml` unchanged.
Keep `configGlobal.toml` protocol/network settings identical across sites; use matching genome builds and gene IDs.

2. Prepare MVP inputs and create the combined public variant list.

```bash
go run -mod=vendor secure-rvas.go prepare --config config/mvp --party 1
```

Output: `<MVP run_dir>/prepared/` and `<MVP run_dir>/public_variants.tsv`.

3. Copy that one TSV to AoU, set its local path in `config/aou/configPrepare.toml`, then prepare AoU.

```toml
public_var_list = "/local/path/public_variants.tsv"
```

```bash
go run -mod=vendor secure-rvas.go prepare --config config/aou --party 2
```

AoU's `chromosomes` must match the file. Output: `<AoU run_dir>/prepared/`.

4. Generate matching keys for all three parties together.

```bash
go run -mod=vendor secure-rvas.go keygen \
  --config-party0 config/cp0 --config-party1 config/mvp --config-party2 config/aou
```

For separate machines, place each party's own key subset in its `<run_dir>/secure/keys/<ancestry>/` directory.

5. Start each party in a separate terminal or on its corresponding machine.

```bash
# Helper
source scripts/setup/env.sh
go run -mod=vendor secure-rvas.go party --config config/cp0 --party 0
```

```bash
# MVP
source scripts/setup/env.sh
go run -mod=vendor secure-rvas.go party --config config/mvp --party 1
```

```bash
# AoU
source scripts/setup/env.sh
go run -mod=vendor secure-rvas.go party --config config/aou --party 2
```

Results: each data party writes `<run_dir>/secure/<ancestry>/all_secure_results.tsv`.

Set `print_logs = true` in each site's `configLocal.toml` to show its stage timings
and timing tree. Set it to `false` for quiet metrics output; result and metrics
files are still written. If omitted, only party 1 prints. The 1KG configs enable
it only for A; the MVP and AoU configs both enable it.
After each ancestry finishes, enabled data parties print a secure-result check:
Burden/SKAT p-value ranges, Q1/median/Q3, non-finite/out-of-range counts, and exact
zero/one counts across all phenotypes/chromosomes. Quantiles use linear
interpolation over finite values, including out-of-range values without clipping.
This needs no plaintext reference.

6. After all parties finish, optionally run the reference, comparison, plots, and metrics in an authorized environment with both sites' prepared inputs available.

```bash
python3 python/analysis/run_reference.py --config config/mvp --config-b config/aou --engine python
python3 python/analysis/compare_secure_to_reference.py --config config/mvp --config-b config/aou
python3 python/analysis/plot_secure_vs_reference.py --config config/mvp --config-b config/aou
python3 python/analysis/summarize_metrics.py --config config/mvp --config-b config/aou
```

Results: both sites receive the joint plain results in `<run_dir>/reference/`
and their secure/plain comparisons and plots in `<run_dir>/comparison/`.
The joint plain reference is computed once. Each party writes its own `metrics/`.
Each completed party run saves `metrics/process_partyN.csv` with its peak RSS in
bytes. The metrics summary displays this in GiB for the entire party process
across all ancestries and chromosomes; older runs without this file show `NA`.
The comparison command honors each site's `print_logs`; the metrics summary uses
the primary config (`--config-b` adds B to that one combined summary). For an
independent party 2 summary, use:

```bash
python3 python/analysis/compare_secure_to_reference.py --config config/aou --party-id 2
python3 python/analysis/summarize_metrics.py --config config/aou --party-id 2
```

R² requires the matching joint reference under that site's `run_dir/reference/`
(or an existing comparison under `run_dir/comparison/` for the summary alone).
The `party` command does not compute R², and a reference using only one site's
samples is not a reference for the joint secure result.

Caches use `<run_dir>/cache/prepared/<ancestry>/chrN/`. Matching `config.json` reuses the cache; changed preparation settings or old hash-folder caches require a new `run_dir`.

For a local test using AoU-derived inputs for both sites, first provide local VAT annotations at `vat_annotation` in `config/aou/configPrepare.toml`, then generate the split in AoU Workbench:

```bash
python3 -m python.datasets.aou.prepare_aou --config config/aou --config-mvp config/mvp
```

Alternatively, after configuring both sites' inputs on one machine, run preparation through comparison automatically (logs: `output/aou-logs/`):

```bash
INPUT_MODE=sites ./scripts/run_aou_workflow.sh
```
