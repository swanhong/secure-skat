# Analysis Scripts

This directory contains analysis helpers plus the Python-first secure-vs-plain
SKAT comparison pipeline.

## SKAT Compare

Primary entrypoint:

```bash
scripts/analysis/skat_compare.py
```

Tiny R helper kept only for the public `SKAT` package reference step:

```bash
scripts/analysis/r_skat_reference.R
```

Internal implementation is split under `scripts/analysis/skat_compare_lib/`:

- `cli.py`: argparse and subcommand dispatch
- `context.py`: dataset/run resolution plus null-model inputs
- `dataset_io.py`: dataset readers and PLINK raw export
- `secure_io.py`: secure output readers
- `compute.py`: secure-compatible SKAT/Burden math
- `reference.py`: `Rscript` bridge to the `SKAT` package
- `reporting.py`: CSV writing and console diagnostics
- `plotting.py`: scatter plot rendering
- `workflow.py`: orchestration of the compare/manual/secure/reference flows

## Environment Setup

Create a repo-local analysis virtualenv and install the committed dependencies:

```bash
python3 -m venv .local/venv/analysis
.local/venv/analysis/bin/pip install -r scripts/analysis/requirements-skat-compare.txt
```

Required tools:

- `python3`
- `plink2`
- `Rscript`

Optional for the `reference`/`compare` subcommands when `--skip-reference` is
not set:

- R package `SKAT`

## Dataset Resolution

When `--dataset` is omitted, `skat_compare.py` resolves the dataset in this
order:

1. `run_metadata.txt` `dataset=` entry
2. Unique local match across `example_data` and `.local/datasets/*` using
   block count plus total variant count inferred from the secure run

The script refuses to compare a run against a dataset whose block count or total
variant count does not match the run metadata/cache hints.

## Main Commands

Auto-resolve the dataset from `run_metadata.txt` and run the full compare flow
without the R package reference:

```bash
python3 scripts/analysis/skat_compare.py compare --run-id ca92 --skip-reference
```

Legacy run without `dataset=` metadata; provide the dataset explicitly:

```bash
python3 scripts/analysis/skat_compare.py compare --run-id e6d9 --dataset example_data --skip-reference
```

Write per-variant debug CSVs:

```bash
python3 scripts/analysis/skat_compare.py compare \
  --run-id ca92 \
  --debug \
  --skip-reference
```

Windowed comparison example:

```bash
python3 scripts/analysis/skat_compare.py compare \
  --run-id ca92 \
  --skip-reference \
  --window-bp 50000 \
  --step-bp 10000 \
  --min-window-variants 2 \
  --window-limit 20 \
  --window-output-tag smoke
```

Run only one stage:

```bash
python3 scripts/analysis/skat_compare.py manual --run-id ca92 --dataset .local/datasets/1000g_all_chr22_anchor50kb_top16
python3 scripts/analysis/skat_compare.py secure --run-id ca92 --dataset .local/datasets/1000g_all_chr22_anchor50kb_top16
python3 scripts/analysis/skat_compare.py reference --run-id ca92 --dataset .local/datasets/1000g_all_chr22_anchor50kb_top16
```

## Shared Arguments

- `--repo-root <path>`
- `--run-id <id>`
- `--dataset <path>`
- `--blocks <spec>`
  - Analysis scope. Defaults to all blocks.
- `--detail-blocks <spec>`
  - Verbose block diagnostics only. Defaults to `1,last`.
- `--debug`
- `--window-bp <int>`
- `--step-bp <int>`
- `--min-window-variants <int>`
- `--window-limit <int>`
- `--window-output-tag <string>`
- `--skip-reference`

Block specs accept comma-separated values such as `all`, `1,last`, or
`1-4,8`.

## Output Notes

Outputs are written under:

```text
.local/tmp/skat_compare/<dataset_tag>/<run_id>/
```

`compare` always writes:

- `summary.csv`
- `block_compare.csv`
- `block_compare_skat_scatter.png`
- `block_compare_burden_scatter.png`

When `--window-bp` is set, it also writes:

- `window_compare_<tag>.csv`
- `window_compare_<tag>_skat_scatter.png`
- `window_compare_<tag>_burden_scatter.png`

When `--debug` is set, it also writes:

- `variant_debug_csv/variant_debug_blockXX.csv`
- `variant_debug_csv/variant_debug_all.csv`

The `summary.csv` `secure` column is always aligned to the selected analysis
blocks. When full-run secure scalars are also present, they are printed in the
console summary as additional diagnostics.
