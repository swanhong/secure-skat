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
- `pipeline.py`: dataset/run resolution, block-wise compare output, summary stats
- `compute.py`: null-model fitting plus secure-compatible SKAT/Burden math
- `test_plain_modes.py`: experimental compare-only plain variants kept separate from the standard math path
- `reference.py`: block-wise `Rscript` bridge to the `SKAT` package
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

For the null model, `skat_compare.py` uses all covariate columns present in each
party's `cov.txt`. The parties must have the same number of covariate columns.

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

Compare an exact run directory, useful when `run_example.sh --run-base` writes
outside the repo-level `out/` directory:

```bash
python3 scripts/analysis/skat_compare.py compare --run-root datasets/1000g/runs/output_260511_123000_ca92 --skip-reference
```

Legacy run without `dataset=` metadata; provide the dataset explicitly:

```bash
python3 scripts/analysis/skat_compare.py compare --run-id e6d9 --dataset example_data --skip-reference
```

Run compare with the plain test-only local-weight burden variant:

```bash
python3 scripts/analysis/skat_compare.py compare --run-id ca92 --plain-mode local-weight-burden --skip-reference
```

Run only the R package reference:

```bash
python3 scripts/analysis/skat_compare.py reference --run-id ca92 --dataset .local/datasets/1000g_all_chr22_anchor50kb_top16
```

## Shared Arguments

- `--repo-root <path>`
- `--run-id <id>`
  - Run id suffix searched under repo-level `out/`.
- `--run-root <path>`
  - Exact run output directory. Use this instead of `--run-id` for custom run bases.
- `--dataset <path>`
- `--blocks <spec>`
  - Analysis scope. Defaults to all blocks.
- `--skip-reference`
  - `compare` subcommand only.
- `--plain-mode {standard,local-weight-burden}`
  - `compare` subcommand only.
  - `standard` keeps the current pooled/global plain SKAT and burden calculation.
  - `local-weight-burden` keeps plain SKAT unchanged and enables the experimental local burden logic below.
- `--local-weight-mode {direct-total,product-approx}`
  - `compare` subcommand only.
  - Used only when `--plain-mode local-weight-burden`.
  - `direct-total`
    - builds `w_local(p) = 25 * max(1 - p_local, p_local)^24`
    - uses each party's local numerator but the global `2N_total` denominator
    - forms each party's local burden linear term and sums those local terms
  - `product-approx`
    - defines each party's local alt-frequency contribution as `x_p = alt_count_p / (2N_total)`
    - builds an approximate shared weight `w_hat = 25 * product_p (1 - x_p)^24`
    - applies that shared approximate weight to the aggregated score term for the matched test data

## Plain Mode Notes

`--plain-mode local-weight-burden` is a compare-only testing mode. It does not
change the Go secure pipeline, and the R reference path still uses the public
`SKAT` package's standard pooled weighting.

Because the R helper continues to use the public `SKAT` package's standard
pooled weighting, its burden output is expected to differ from the manual/plain
burden output when `--plain-mode local-weight-burden` is enabled.

The `product-approx` submode is a matched-data approximation experiment inside
the compare pipeline. It is useful for testing a separable weight construction,
but it still relies on the current compare script's matched per-variant test
inputs.

## Reference Memory Notes

The R reference helper now runs block-by-block with a shared null model instead
of concatenating every selected block into one giant genotype matrix. This
significantly reduces peak reference memory use when `compare` runs without
`--skip-reference`.

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
- `reference_block_summary.tsv` when the reference step is enabled
- `block_compare_skat_scatter.png`
- `block_compare_burden_scatter.png`

`block_compare.csv` contains per-block plain/manual, secure, and reference Q
values plus pairwise absolute/relative differences.

`summary.csv` contains block-wise error summaries such as max absolute
difference and max relative difference for:

- `plain_vs_reference`
- `secure_vs_reference`
- `plain_vs_secure`

No final aggregate reference scalar is reported; the compare output is centered
on per-block Q values and worst-case block errors.
