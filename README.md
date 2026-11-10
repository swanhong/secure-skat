# Secure RVAS

Joint Burden/SKAT with MVP (Party 1), AoU (Party 2), and a helper (Party 0).
Run commands from the repository root. For a public-data test, see the [1000 Genomes workflow](python/datasets/onekg/README.md).

## 1. Install

Linux x86-64: requires Go, Python 3.11+ in a conda/venv environment, and R.

```bash
bash scripts/setup/install.sh
source scripts/setup/env.sh
go mod vendor
```

On macOS, install PLINK 2 and R::SKAT separately, then run `python3 -m pip install -r requirements.txt`.
Load `scripts/setup/env.sh` in each new terminal.

## 2. Configure

- In `config/{mvp,aou}/configPrepare.toml`, set `run_dir` and genotype, gene panel, annotation, phenotype, covariate, and ancestry paths. Use `samples_per_cohort = 0` for all eligible samples.
- Supply annotation TSVs with exactly `variant_key`, `gene_id`, and `gene_symbol`, filtered to **MANE Select + HC LoF**. Exclude variants without MANE Select. For AoU, use `python/datasets/aou/vat_simplify.py`.
- Set `local_max_maf` per site (e.g. MVP `0.05`, AoU `0.02`). Supply ALL covariates even when running only ancestry-specific analyses.
- Set the helper's `run_dir` in `config/cp0/configPrepare.toml`.
- Match `configGlobal.toml` across all three parties: chromosomes, phenotypes, ancestries, protocol parameters, and network addresses. Use matching genome builds and gene IDs.
- Set `global_max_maf = [0.01, 0.001]` for both thresholds, computed over the joint ALL cohort for every ancestry analysis.

Use a new `run_dir` when changing preparation settings.

## 3. Prepare MVP, then AoU

```bash
go run -mod=vendor secure-rvas.go prepare --config config/mvp --party 1
```

Copy `<MVP run_dir>/public_variants.tsv` and `local_cut_variants.tsv` into the same directory on AoU. Set the following in `config/aou/configPrepare.toml`:

```toml
public_var_list = "/local/path/public_variants.tsv"
```

```bash
go run -mod=vendor secure-rvas.go prepare --config config/aou --party 2
```

Prepared inputs: `<run_dir>/prepared/`.

## 4. Generate keys

```bash
go run -mod=vendor secure-rvas.go keygen \
  --config-party0 config/cp0 --config-party1 config/mvp --config-party2 config/aou
```

For separate machines, place each party's own key subset in its `<run_dir>/secure/keys/<ancestry>/` directory.

## 5. Run all three parties

Run each command in a separate terminal or on its corresponding machine:

```bash
# Helper
go run -mod=vendor secure-rvas.go party --config config/cp0 --party 0

# MVP
go run -mod=vendor secure-rvas.go party --config config/mvp --party 1

# AoU
go run -mod=vendor secure-rvas.go party --config config/aou --party 2
```

Results: `<run_dir>/secure/<ancestry>/maf<threshold>/all_secure_results.tsv`.
Without `global_max_maf`, the `maf<threshold>/` directory is omitted.

## 6. Compare with plaintext (optional)

After the run finishes, execute where both sites' prepared inputs are available:

```bash
python3 python/analysis/run_reference.py --config config/mvp --config-b config/aou --engine python
python3 python/analysis/compare_secure_to_reference.py --config config/mvp --config-b config/aou
python3 python/analysis/plot_secure_vs_reference.py --config config/mvp --config-b config/aou
python3 python/analysis/summarize_metrics.py --config config/mvp --config-b config/aou
```

Outputs: `<run_dir>/reference/`, `<run_dir>/comparison/`, and `<run_dir>/metrics/`.

For both sites on one machine, run steps 3–6 together:

```bash
INPUT_MODE=sites ./scripts/run_aou_workflow.sh
```
