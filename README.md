# Secure RVAS

Secure RVAS runs privacy-preserving Burden and SKAT tests with MVP (Party 1),
AoU (Party 2), and a helper without genomic inputs (Party 0).

## Step-by-step commands

### 1. Install and configure

On Linux x86-64, install Go, Python 3.11+ (with pip), and R, then run from the repository root:

```bash
bash scripts/setup/install.sh
source scripts/setup/env.sh
go mod vendor
```

Edit each site's `configPrepare.toml` for `run_dir`, input paths, and preparation
settings. `configLocal.toml` contains fixed `{run_dir}` path templates and local
thread settings; it normally needs no edits. MVP, AoU, and the helper each have
an identical `configGlobal.toml` for shared protocol/network settings.
Every stage reads Global, Local, then Prepare before expanding `{run_dir}`.
The helper's prepare config only needs `run_dir`.

Expected output: PLINK 2 at `$HOME/plink2`, R::SKAT under `$HOME/R/library`,
Python packages from `requirements.txt` in the current `python3` environment,
and Go dependencies in `vendor/`. Use the same Python environment when running
the workflow, and load `scripts/setup/env.sh` in each new terminal.
On macOS, install PLINK 2 and R::SKAT separately, and install Python packages with
`python3 -m pip install -r requirements.txt`; the installer above is Linux-only.

### 2. Provide site inputs

Choose one input mode:

- **AoU split test:** localize AoU inputs and split samples/variants into A/B.
  Party 1 represents MVP but contains AoU-derived test data.
- **Independent MVP/AoU:** provide each site's PGEN, phenotype, covariate,
  ancestry, annotation, and gene-panel files at its configured paths; skip the
  split command. MVP input creation is external to this workflow.

For the AoU split test, first place chromosome annotations at
`~/fed_prep_out/chr{chromosome}_annotation.tsv` (see Additional information),
then run in an authorized AoU Workbench environment:

```bash
python3 -m python.datasets.aou.prepare_aou \
  --config config/aou --config-mvp config/mvp
```

Expected output: localized source files in `data/aou/generated/`, MVP test inputs
in `data/mvp/generated/`, AoU test inputs in `data/aou/split/`, and a `split.json`
in each site input root. Matching completed splits are reused; changed
sources/settings or incomplete outputs require regenerating both sites in empty
directories.

### 3. Prepare secure inputs

```bash
go run -mod=vendor secure-rvas.go prepare --config config/mvp --party 1

# On separate machines, copy MVP's public_variants.tsv to AoU and update public_var_list first.
go run -mod=vendor secure-rvas.go prepare --config config/aou --party 2
```

Expected output: `<run_dir>/prepared/<ancestry>/chr<chromosome>/A/` for MVP and
the corresponding `B/` for AoU, with phenotype, covariate, and genotype blocks.
`genes.txt`, `block_sizes.txt`, `pos.txt`, and `public_variants.tsv` are at the
chromosome root, beside `A/` or `B/`.
Party 1 also writes all configured chromosomes to `<run_dir>/public_variants.tsv`;
transfer this single file to Party 2.
Actual blocks live under `<run_dir>/cache/prepared/<ancestry>/chr<chromosome>/`;
the `prepared/` chromosome paths are symlinks to those caches. Matching completed
caches are reused after checking their `config.json`. Changed preparation settings
or old hash-folder caches require a new `run_dir`. Failed builds are discarded
and rebuilt on the next run.
Prepared caches and secure, reference, comparison, and metrics outputs stay under
`run_dir`; input-generation paths are configured separately.
AoU follows that list's gene/column order, zero-fills absent public variants,
and places locally selected variants outside the list in `B/private/`.

AoU's `configPrepare.toml` defaults to:

```toml
public_var_list = "output/mvp/public_variants.tsv"
```

This default assumes MVP's `run_dir` is `output/mvp`.
For separate machines or custom MVP paths, copy this file to AoU (for example,
`<run_dir>/inputs/public_variants.tsv`) and update this setting.
Party 2's `chromosomes` must match the chromosomes in the file.

### 4. Generate shared PRG keys

```bash
go run -mod=vendor secure-rvas.go keygen \
  --config-party0 config/cp0 \
  --config-party1 config/mvp \
  --config-party2 config/aou
```

Expected output: each party's key subset under its configured
`shared_keys_path/<ancestry>/` (default: `<run_dir>/secure/keys/<ancestry>/`).
For separate machines, provision each site's own subset at its configured path.

### 5. Run secure Burden and SKAT

Start these concurrently in three terminals, or on the corresponding machines:

```bash
# Helper
go run -mod=vendor secure-rvas.go party --config config/cp0 --party 0

# MVP
go run -mod=vendor secure-rvas.go party --config config/mvp --party 1

# AoU
go run -mod=vendor secure-rvas.go party --config config/aou --party 2
```

Expected output: secure result TSVs under MVP's `<run_dir>/secure/<ancestry>/`,
including `chr<chromosome>.tsv` and `all_secure_results.tsv`.
Each party writes timing/communication CSVs under its own
`<run_dir>/metrics/<ancestry>/`.

### 6. Run the reference, compare, and plot

For validation in an environment with both sites' prepared inputs available:

```bash
python3 python/analysis/run_reference.py \
  --config config/mvp --config-b config/aou --engine python
python3 python/analysis/compare_secure_to_reference.py --config config/mvp
python3 python/analysis/plot_secure_vs_reference.py --config config/mvp
python3 python/analysis/summarize_metrics.py \
  --config config/mvp --config-b config/aou
```

Expected output under MVP's `run_dir` (default: `output/mvp/`):

```text
reference/<ancestry>/all_r_results.csv
comparison/<ancestry>/all_comparison.csv
comparison/<ancestry>/plots/
```

The final command prints combined MVP/AoU timing, communication, and p-value
accuracy summaries. Use `--engine r` for the R::SKAT reference.

## Complete workflow

After setup, configuration, and annotation/input preparation, run all steps on
one machine:

```bash
# AoU-derived A/B split test (default)
./scripts/run_aou_workflow.sh

# Independent inputs, or reuse an existing split without localization/splitting
INPUT_MODE=sites ./scripts/run_aou_workflow.sh
```

`CONFIG_A`, `CONFIG_B`, and `CONFIG_0` override the three config directories;
`REFERENCE_ENGINE=r` selects R::SKAT. Party logs are written to
`output/aou-logs/partyN.log`; `LOG_DIR` overrides this directory.

## Additional information

- Both sites need the same genome build, biallelic variant IDs in
  `chromosome:position:REF:ALT` form, and corresponding gene IDs.
- Annotations must contain the configured mask columns. The current
  `annotation=pLoF` mask requires an `annotation` column containing `pLoF`;
  `vat_simplify.py` emits `LoF`, so supply the required annotation or explicitly
  configure a matching mask such as `LoF=HC` in both sites.
- `samples_per_cohort=0` uses all eligible local samples. In split mode, sample
  selection and variant-filter settings must match between sites.
- Party 2 requires MVP's public list. `--public-var-list PATH` overrides its
  config setting. Keep received lists outside directories removed by `--clear`;
  `prepare --clear` removes prepared links and derived results while retaining
  cached blocks, source/site inputs, and logs.
- When updating an existing checkout, move `run_dir` from Global/Local to
  `configPrepare.toml` and remove `prepared_cache_dir`. Use `configLocal.toml`;
  legacy `configLocal.PartyN.toml` remains a fallback. Existing external caches
  are not moved or deleted automatically; a new run directory builds fresh
  caches on its own disk. Keep existing input paths if reusing already prepared
  PGEN and table inputs.
- The plaintext reference uses both sites' inputs. It is a validation step in
  an authorized environment, separate from the distributed secure run.
- Secure SKAT uses Wilson–Hilferty; reference outputs include Liu and Davies
  p-values. Timing stages may overlap and should not be summed as wall time.
- An `EOF` or connection reset often follows another party's earlier error;
  inspect all party logs. Reduce `mpc_num_threads` for memory-constrained runs.
- The public-data test is documented in
  [`python/datasets/onekg/README.md`](python/datasets/onekg/README.md).

To create the AoU VAT annotation for the default chromosome-22 test (Bash,
AoU Workbench with `GOOGLE_PROJECT` set):

```bash
annotation_dir="$HOME/fed_prep_out"
mkdir -p "$annotation_dir"
python3 python/datasets/aou/vat_simplify.py \
  --chromosome 22 \
  --pvar <(gsutil -u "$GOOGLE_PROJECT" cat \
    gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/exome/pgen/exome.chr22.pvar) \
  --output "$annotation_dir/chr22_annotation.tsv"
```

Repeat for each configured chromosome. VAT scanning can take a long time;
`prepare_aou` consumes local annotation paths and reuses existing localized files.
Set `vat_annotation` in `config/aou/configPrepare.toml` to change the VAT input
path; `--annotation` overrides it. This is separate from `annotation`, which
points to the normalized input used by secure preparation.
