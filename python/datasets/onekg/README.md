# Local 1000 Genomes test

Test secure Burden/SKAT against a plaintext reference using public 1000 Genomes data and synthetic annotations/phenotypes.
Synthetic HC selection uses seed 42: one row is selected from each complete block of 100 variant–gene annotations; the final incomplete block is discarded. These are random test labels, not LOFTEE/MANE annotations. The generator preserves the legacy random draws and writes only the previously HC-labeled rows as `variant_key`, `gene_id`, `gene_symbol`. Preparation no longer applies an annotation mask or reads annotation frequencies.

Old six-column annotations are rejected. Preserve old outputs separately and regenerate the gene panel/annotation files and prepared caches before using this version. Genotype downloads do not need to be repeated.

The fixture uses the chromosomes and analyses in `configGlobal.toml` and keeps
A/B samples disjoint with `debug_split_num_parties = 2` in both preparation configs.

`print_logs = true` in A's local config enables stage timings and the timing tree;
B and the helper use `false`. `run_1kg_workflow.sh` displays A's output live while
also saving it to `output/1kg-logs/party1.log` (or `$LOG_DIR/party1.log`).

Run all commands from the repository root using the same Python environment.

0. Install dependencies once (Linux x86-64; requires Go, Python 3.11+ with pip, and R).

```bash
bash scripts/setup/install.sh
```

On macOS, install PLINK 2 and R::SKAT separately, then run:

```bash
python3 -m pip install -r requirements.txt
```

1. Download inputs and generate test data (`data/1kg/generated/`).

```bash
source scripts/setup/env.sh
python3 python/datasets/onekg/prepare_1kgenome.py --config config/1kg-A --num-pheno 2
```

2. Prepare A, then B using A's combined `<A run_dir>/public_variants.tsv`.

`debug_split_num_parties` is prepare-only: `0` disables it (the default), and `2`
partitions the full source PSAM IDs by original ancestry. IDs are sorted and
assigned alternately to A/B, with odd leftovers balanced across ancestries.
`--party 1` and `--party 2` select their respective partitions, consistently
across EUR/ALL analyses and chromosomes. `samples_per_cohort` caps each site's
eligible samples after partitioning. Both configs must use the same full source
PSAM ID set and ancestry mapping. Do not enable this for already disjoint real
sites or AoU artificial split inputs.

If an existing preparation used unsplit samples, choose new `run_dir` values
before enabling the option, and update B's `public_var_list` to A's new path.

```bash
go run -mod=vendor secure-rvas.go prepare --config config/1kg-A --party 1
go run -mod=vendor secure-rvas.go prepare --config config/1kg-B --party 2
```

3. Generate matching keys for all three parties together.

```bash
go run -mod=vendor secure-rvas.go keygen \
  --config-party0 config/1kg-cp0 --config-party1 config/1kg-A --config-party2 config/1kg-B
```

4. Start each party in a separate terminal, from the repository root.

```bash
# Terminal 1
source scripts/setup/env.sh
go run -mod=vendor secure-rvas.go party --config config/1kg-cp0 --party 0
```

```bash
# Terminal 2
source scripts/setup/env.sh
go run -mod=vendor secure-rvas.go party --config config/1kg-A --party 1
```

```bash
# Terminal 3
source scripts/setup/env.sh
go run -mod=vendor secure-rvas.go party --config config/1kg-B --party 2
```

5. After all parties finish, run the reference, compare results, and create plots/metrics.

```bash
python3 python/analysis/run_reference.py --config config/1kg-A --config-b config/1kg-B --engine python
python3 python/analysis/compare_secure_to_reference.py --config config/1kg-A
python3 python/analysis/plot_secure_vs_reference.py --config config/1kg-A
python3 python/analysis/summarize_metrics.py --config config/1kg-A --config-b config/1kg-B
```

Results: `<A run_dir>/{secure,reference,comparison,metrics}/`.
Set each party's `run_dir` in its `configPrepare.toml`; caches use `<run_dir>/cache/prepared/<ancestry>/chrN/`.
Matching `config.json` reuses the cache; changed settings or old hash-folder caches require a new `run_dir`.
If A's `run_dir` changes, update B's `public_var_list` to match.

Alternatively, run every step automatically (logs: `output/1kg-logs/`):

```bash
./scripts/run_1kg_workflow.sh
```
