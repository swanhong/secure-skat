# Local 1000 Genomes workflow

This workflow runs the secure Burden/SKAT pipeline locally with public 1000 Genomes data and synthetic annotations / phenotypes. The generated fixture is for correctness testing, not biological interpretation.

Run commands from the repository root.

## Download and generate the test data

```bash
source scripts/setup/env.sh
python3 python/datasets/onekg/prepare_1kgenome.py \
  --config config/1kg-A \
  --num-pheno 2
```

The generator reads the chromosomes from `config/1kg-A/configGlobal.toml`. It downloads missing
1000 Genomes VCFs and sample metadata plus GENCODE v50, then writes the local
fixture under:

```text
data/1kg/generated/
```

Downloaded files are reused on later runs.

## Test separate A and B preparation

`config/1kg-A`, `config/1kg-B`, and `config/1kg-cp0` contain only their own local
party configuration. This smoke test uses chr22, EUR, 100 samples, and three
A-selected genes. A uses MAF <= 0.001; B uses MAF <= 0.01 to allow private
variants. Both sites use the same input and sample seed, so samples overlap.
This fixture tests preparation and connectivity, not independent cohorts.

Run these steps from the repository root:

```bash
python3 python/datasets/onekg/prepare_1kgenome.py --config config/1kg-A --num-pheno 2
go run -mod=vendor secure-rvas.go prepare --config config/1kg-A --party 1
go run -mod=vendor secure-rvas.go prepare --config config/1kg-B --party 2
go run -mod=vendor secure-rvas.go keygen --config-party0 config/1kg-cp0 --config-party1 config/1kg-A --config-party2 config/1kg-B
```

B reads the public lists from `public_var_list` in its prepare configuration.
Keys are written to each party's `<run_dir>/secure/keys/EUR/` directory.
Generate all keys in one invocation; separate invocations create different keys.

Start the following in three terminals at the same time:

```bash
# Terminal 1
go run -mod=vendor secure-rvas.go party --config config/1kg-cp0 --party 0
# Terminal 2
go run -mod=vendor secure-rvas.go party --config config/1kg-A --party 1
# Terminal 3
go run -mod=vendor secure-rvas.go party --config config/1kg-B --party 2
```

Use direct `party` commands for these separated configurations. The `run`
command expects all three local party files in one configuration directory.
The reference accepts separate A/B configurations with `--config` and `--config-b`.

## Run the complete workflow

```bash
./scripts/run_1kg_workflow.sh
```

The workflow generates inputs, prepares A then B, generates matching keys,
and starts all three parties concurrently with `go run`. Logs are saved under
`output/1kg-logs/`. It waits for all three parties, then reads A/B inputs directly for the reference.
Results are saved under `output/1kg-A/`. It runs the comparison, generates plots,
and summarizes A/B metrics. Reference results use the same prepared A/B inputs, including their
sample overlap.

Override `CONFIG_A`, `CONFIG_B`, `CONFIG_0`, `LOG_DIR`, or `REFERENCE_ENGINE`
(default `python`) through environment variables. With custom directories,
update B's `public_var_list` to point at A's output.

Setup, configuration, individual stages, and output details are the same as the AoU workflow. See the repository [`README.md`](../../../README.md).
