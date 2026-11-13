# Secure RVAS

Joint Burden/SKAT with MVP (Party 1), AoU (Party 2), and a helper (Party 0).
Run commands from the repository root. For a public-data test, see the [1000 Genomes workflow](python/datasets/onekg/README.md).

## 1. Install

Linux x86-64: requires Go, Python 3.11+ in a conda/venv environment, and R.

```bash
bash scripts/setup/install.sh
source scripts/setup/env.sh
go mod download all && go mod vendor
```

On macOS, install PLINK 2 and R::SKAT separately, then run `python3 -m pip install -r requirements.txt`.
Load `scripts/setup/env.sh` in each new terminal.

## 2. Configure

- In `config/{mvp,aou}/configPrepare.toml`, set `run_dir` and genotype, gene panel, annotation, phenotype, covariate, and ancestry paths. Use `samples_per_cohort = 0` for all eligible samples.
- Supply annotation TSVs containing `variant_key`, `gene_id`, and `gene_symbol`, filtered to **MANE Select + HC LoF**. Column order is flexible; extra columns are ignored and do not filter rows. Exclude variants without MANE Select. For AoU, use `python/datasets/aou/vat_simplify.py`.
- Set `local_max_maf` per site (e.g. MVP `0.05`, AoU `0.02`). Supply ALL covariates even when running only ancestry-specific analyses.
- Set the helper's `run_dir` in `config/cp0/configPrepare.toml`; it does not need `covariate_columns`. At startup, Party 0 receives the covariate counts from Parties 1 and 2, and all parties stop if they differ. Keep the selected covariates in matching semantic order on the two data holders; this exchange checks counts, not column meanings. Run the updated executable on all three parties.
- Match `configGlobal.toml` across all three parties: chromosomes, phenotypes, ancestries, protocol parameters, and network addresses. Use matching genome builds and gene IDs.
- Set `global_max_maf = [0.01, 0.001]` for both thresholds. Joint analysis uses the combined ALL cohort; independent local analysis uses the site's own ALL cohort. Both apply the same threshold across ancestries.

Use a new `run_dir` when changing preparation settings.

Full analysis uses `config/{mvp,aou,cp0}`: chr1–22, five phenotypes, AFR/EUR/AMR/ALL, and all eligible samples and genes.
Local tests use `config/{mvp-test,aou-test,cp0-test}`: chr21–22, two phenotypes, EUR/AFR, up to 100 samples per ancestry per site, and 50 genes per chromosome. Both split sites use the same 15 AoU covariates, including the original BIOBANK values, from the workspace joint files. Regenerate older 10-PC split inputs and prepared outputs in fresh directories before running this configuration.

`shared_rate` defaults to `1.0` and controls only synthetic variant splitting. The MVP/AoU `-test` configs explicitly use `0.6` (approximately 60% shared, 20% A-only, 20% B-only). Real site preparation does not use this rate: it retains the site's eligible inputs subject to QC, MAF, and gene selection, and determines overlap from the actual variant lists. Use `INPUT_MODE=sites` for real MVP/AoU runs; the workflow script permits split mode only with `-test` config paths.

## 3. Prepare MVP, then AoU

```bash
go run -mod=vendor secure-rvas.go prepare --config config/mvp --party 1
```

Copy `<MVP run_dir>/public_variants.tsv` and `local_cut_variants.tsv` into the same directory on AoU. Set the following in `config/aou/configPrepare.toml`:

```toml
public_var_list = "/local/path/public_variants.tsv"
```

Prepare AoU source inputs before running the site preparation:

```bash
python3 -m python.datasets.aou.prepare_aou --config config/aou
```

AoU preparation reads covariates only from `~/workspace/gwas-data-wgs/results/joint/{anc}/covariates.csv` for `afr`, `amr`, `eur`, and `all`, including ALL when only ancestry-specific analyses are selected. Missing local files stop preparation; there is no covariate download or fallback. The importer accepts comma/tab-separated files and joint exports with a space-separated tail, and writes the TSVs configured by `covariate`. It maps `#IID` to `IID`, `AGE` to `age_at_profile`, `GENDER` to `GENDER_REGULAR`, `GENDER_OTHER` to `OTHER_GENDER`, `AGESQ` to `age_sq`, and `BIOBANK` to `affiliation`, preserving the source values. Samples are subsequently joined by IID and emitted in genotype sample order. No intercept is stored in these 15 covariates; the analysis adds it.

```bash
go run -mod=vendor secure-rvas.go prepare --config config/aou --party 2
```

By default, prepare writes both joint inputs to `<run_dir>/prepared/<ANC>/` and independent local inputs to `<run_dir>/prepared/local/<ANC>/`. Use `--mode joint` or `--mode local` to prepare only one. Local-only preparation does not need the other site's public or local-cut lists. `--clear` removes outputs and prepared caches for the selected mode.

Local preparation retains variants up to `max(global_max_maf)` using the site's ALL cohort. `local_max_maf` applies only to joint preparation. Use a new run directory or regenerate local inputs before increasing the maximum threshold.

Run independent local Burden/SKAT outside Docker/sfkit, using the same command on each site with its own config:

```bash
python3 python/analysis/run_reference.py --config config/mvp --mode local
# On AoU:
python3 python/analysis/run_reference.py --config config/aou --mode local
```

This uses the existing Python plain analysis, including R::SKAT's Davies library. Results: `<run_dir>/local/<ANC>/maf<threshold>/all_local_results.csv`, with the same columns as the joint reference. No other site's inputs are required.

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
`global_max_maf` is required; use `[0.5]` to retain all prepared variants.

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

For a small test using AoU data split into two sites:

```bash
./scripts/run_aou_workflow.sh
```
