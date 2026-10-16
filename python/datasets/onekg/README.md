# Local 1000 Genomes test

Test secure Burden/SKAT against a plaintext reference using public 1000 Genomes data and synthetic annotations/phenotypes.
The default fixture uses chr22/EUR with overlapping A/B samples; it is for correctness testing.

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

Results: `output/1kg-A/{secure,reference,comparison,metrics}/`.
Set each party's `run_dir` in its `configPrepare.toml`; caches use `<run_dir>/cache/prepared/<ancestry>/chrN/`.
Matching `config.json` reuses the cache; changed settings or old hash-folder caches require a new `run_dir`.
If A's `run_dir` changes, update B's `public_var_list` to match.

Alternatively, run every step automatically (logs: `output/1kg-logs/`):

```bash
./scripts/run_1kg_workflow.sh
```
