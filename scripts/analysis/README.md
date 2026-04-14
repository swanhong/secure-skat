# Analysis Scripts

This directory contains analysis helpers plus the plain-vs-secure SKAT
comparison script.

## Plain SKAT Compare

Entrypoint:

```r
scripts/analysis/compare_plain_skat_example.R
```

The entrypoint is intentionally thin. Most logic now lives in:

- `scripts/analysis/plain_skat_compare/args.R`
- `scripts/analysis/plain_skat_compare/data_io.R`
- `scripts/analysis/plain_skat_compare/secure_io.R`
- `scripts/analysis/plain_skat_compare/skat.R`
- `scripts/analysis/plain_skat_compare/burden.R`
- `scripts/analysis/plain_skat_compare/windows.R`
- `scripts/analysis/plain_skat_compare/reporting.R`
- `scripts/analysis/plain_skat_compare/workflow.R`

## Prerequisites

- `Rscript`
- `plink2`
- Optional: the R `SKAT` package if you want direct package-to-plain comparison

## Test Commands

Parse-check the entrypoint and all submodules:

```bash
Rscript -e 'files <- c("scripts/analysis/compare_plain_skat_example.R", list.files("scripts/analysis/plain_skat_compare", full.names = TRUE)); invisible(lapply(files, parse)); cat("parse_ok\n")'
```

Load the entrypoint without running the full workflow:

```bash
Rscript -e 'source("scripts/analysis/compare_plain_skat_example.R", local = new.env(parent = globalenv())); cat("source_ok\n")'
```

Pick the latest secure run id from `out/output_*`:

```bash
RUN_ID=$(basename "$(ls -dt out/output_* | head -n 1)" | awk -F_ '{print $NF}')
printf '%s\n' "$RUN_ID"
```

Minimal smoke test against block 1 without requiring the `SKAT` package:

```bash
RUN_ID=$(basename "$(ls -dt out/output_* | head -n 1)" | awk -F_ '{print $NF}')
Rscript scripts/analysis/compare_plain_skat_example.R --skip-skat-package . "$RUN_ID" 1
```

Same run, but also write per-block and merged debug CSVs under `.local/tmp/plain_skat_compare/...`:

```bash
RUN_ID=$(basename "$(ls -dt out/output_* | head -n 1)" | awk -F_ '{print $NF}')
Rscript scripts/analysis/compare_plain_skat_example.R --debug --skip-skat-package . "$RUN_ID" 1
```

Windowed comparison example:

```bash
RUN_ID=$(basename "$(ls -dt out/output_* | head -n 1)" | awk -F_ '{print $NF}')
Rscript scripts/analysis/compare_plain_skat_example.R \
  --skip-skat-package \
  --window-bp 50000 \
  --step-bp 10000 \
  --min-window-variants 2 \
  --window-limit 20 \
  --window-output-tag smoke \
  . "$RUN_ID" 1
```

If the R `SKAT` package is installed, run the same comparison with direct
package reference enabled:

```bash
RUN_ID=$(basename "$(ls -dt out/output_* | head -n 1)" | awk -F_ '{print $NF}')
Rscript scripts/analysis/compare_plain_skat_example.R . "$RUN_ID" 1
```

## Arguments

Positional arguments:

1. `repo_root`
2. `run_id`
3. `blocks_to_print`

Supported flags:

- `--dataset <path>`
- `--debug`
- `--skip-skat-package`
- `--window-bp <int>`
- `--step-bp <int>`
- `--min-window-variants <int>`
- `--window-limit <int>`
- `--window-output-tag <string>`

## Output Notes

- Console output prints intermediate comparisons, block-level diagnostics, and
  final SKAT/Burden summaries.
- Debug CSVs are written only when `--debug` is enabled.
- Window CSV and scatter plots are written only when `--window-bp` is provided.
- Cached PLINK `.raw` exports are written under:

```text
.local/tmp/plain_skat_compare/
```
