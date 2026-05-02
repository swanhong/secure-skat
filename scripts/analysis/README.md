# Analysis Scripts

This directory contains analysis helpers plus the simplified plain-vs-secure
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
- `pipeline.py`: dataset/run resolution, file loading, block compare output
- `compute.py`: null-model fitting plus secure-compatible SKAT/Burden math
- `reference.py`: `Rscript` bridge to the `SKAT` package
- `plotting.py`: block scatter plot rendering

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
- `matplotlib` via `requirements-skat-compare.txt`

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

Run only the R package reference:

```bash
python3 scripts/analysis/skat_compare.py reference --run-id ca92 --dataset .local/datasets/1000g_all_chr22_anchor50kb_top16
```

## Shared Arguments

- `--repo-root <path>`
- `--run-id <id>`
- `--dataset <path>`
- `--blocks <spec>`
  - Analysis scope. Defaults to all blocks.
- `--skip-reference`
  - `compare` subcommand only.

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

The `summary.csv` `secure` column is always aligned to the selected analysis
blocks. When full-run secure scalars are also present, they are printed in the
console summary as additional diagnostics.
