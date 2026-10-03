# Secure RVAS

Secure RVAS runs privacy-preserving Burden and SKAT rare-variant association tests.
The local 1000 Genomes test workflow is documented separately in [`python/datasets/onekg/README.md`](python/datasets/onekg/README.md).

## Initial setup

Run these commands once after a fresh clone in a Linux x86-64 environment:

```bash
cd "$HOME/secure-skat"
bash scripts/setup/install.sh
source scripts/setup/env.sh
go mod vendor
```

`scripts/setup/install.sh` installs PLINK 2 at `$HOME/plink2` and R::SKAT under `$HOME/R/library`. `go mod vendor` is required because every Go command in this repository uses `-mod=vendor`, while `vendor/` is not stored in Git.

In each new terminal, load the installed paths before running individual commands:

```bash
cd secure-skat
source scripts/setup/env.sh
```

## Configure the AoU run

Review `config/aou/configPrepare.toml`, `configGlobal.toml`, and each
`configLocal.PartyN.toml` before running. `configPrepare.toml` is used only to
create A/B secure-ready inputs; each party run reads the global file and its own
local file. `mpc_num_threads` is the number of independent MPC lanes, not just
a CPU thread count. More lanes increase concurrent memory and communication use.

`prepare` is optional. If Party 1 and Party 2 already have genotype blocks,
phenotype, covariate, gene, and variant-count files in the secure-ready format,
set those paths in their local configs and start with `keygen`.

All stages read the same non-empty autosome array from the configuration:

```toml
chromosomes = [
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
]
```

The current AoU configuration selects high-confidence loss-of-function variants
with MAF at most 0.02:

```toml
masks = ["LoF=HC"]
max_maf = 0.02
```

## Step 0: Create chromosome VAT annotations

`vat_simplify.py` reads the chromosome PVAR from Google Cloud Storage without
first saving a separate local copy. The `--pvar` argument itself accepts a local
file descriptor, so Bash process substitution connects `gsutil cat` to the
Python command.

```bash
cd "$HOME/secure-skat"
source scripts/setup/env.sh

billing_project="${GOOGLE_PROJECT:-${GOOGLE_CLOUD_PROJECT:-}}"
if [[ -z "$billing_project" ]]; then
  echo "Set GOOGLE_PROJECT or GOOGLE_CLOUD_PROJECT" >&2
  exit 1
fi
export GOOGLE_PROJECT="$billing_project"

mkdir -p "$HOME/fed_prep_out"

for chromosome in {1..22}; do
  echo "===== chromosome ${chromosome} ====="

  pvar_uri="gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/exome/pgen/exome.chr${chromosome}.pvar"
  output_path="$HOME/fed_prep_out/chr${chromosome}_annotation.tsv"

  gsutil -u "$billing_project" ls "$pvar_uri" || break

  python3 python/datasets/aou/vat_simplify.py \
    --chromosome "$chromosome" \
    --pvar <(gsutil -u "$billing_project" cat "$pvar_uri") \
    --output "$output_path" || break
done
```

This code generates the configured autosomes sequentially. Each chromosome scans a large portion of the remote VAT and can take a long time.

The default VAT is:
```text
gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz
```

The simplifier applies the following deterministic contract:

- retain PVAR rows whose `FILTER` is `PASS` or `.`, when `FILTER` is present;
- retain biallelic PVAR rows only;
- match VAT rows by chromosome, position, REF, and ALT;
- group transcript rows by `(variant, gene_id)`;
- select one row using MANE Select, canonical transcript, `LoF=HC`, then VAT
  source order as the priority;
- write `variant_key`, gene fields, `LoF`, consequence, and gnomAD/GVS allele
  frequencies.

## Run the complete AoU workflow

After Step 0 has produced every configured chromosome annotation:

```bash
./scripts/run_aou_workflow.sh
```

The script performs the following stages:

```text
0. Localize AoU PGEN, phenotype, and ancestry inputs; normalize VAT annotations
1. Preprocess ancestry-specific A/B secure inputs
2. Generate persistent ancestry-specific shared PRG keys
3. Run secure Burden and SKAT
4. Run the Python or R::SKAT reference
5. Join and compare secure/reference results
6. Generate scatter and Manhattan plots
7. Summarize timing, communication, and accuracy metrics
```

Step 0 inside `run_aou_workflow.sh` is distinct from the VAT extraction above:
it downloads the PGEN triplets and other AoU inputs, then consumes
`$HOME/fed_prep_out/chr{chromosome}_annotation.tsv` to create normalized local
annotations and gene panels.

The workflow's preprocessing command includes `--clear`; rerunning the complete
script replaces the configured `run_dir` before creating new secure inputs.

To select another configuration or the R::SKAT reference engine:

```bash
CONFIG_PATH=config/aou REFERENCE_ENGINE=r ./scripts/run_aou_workflow.sh
```

### Detached execution and monitoring

Use a date-and-time log path and print it before detaching:

```bash
log_path="run-aou-$(date +'%m%d-%H%M').log"
nohup env PYTHONUNBUFFERED=1 ./scripts/run_aou_workflow.sh \
  > "$log_path" 2>&1 < /dev/null &
echo $! > "${log_path%.log}.pid"
echo "log: $log_path"
```

Check the log with:

```bash
tail -f "$log_path"
```

## Run AoU stages individually

All commands below run from the repository root.

### 1. Localize and normalize AoU inputs

```bash
source scripts/setup/env.sh

python3 python/datasets/aou/prepare_aou.py \
  --config config/aou
```

`prepare_aou.py` downloads missing PGEN/PVAR/PSAM, phenotype, and ancestry
files for the chromosomes listed in `config/aou/configGlobal.toml`. Existing localized files are
reused. It normalizes the Step 0 VAT output, computes minor allele frequency
from `gnomad_af`, and writes inputs under `data/aou/generated/`.

### 2. Preprocess secure inputs

```bash
go run -mod=vendor secure-rvas.go prepare \
  --config config/aou \
  --clear
```

Omit `--clear` when the configured `run_dir` must not be replaced.

### Prepare MVP and AoU separately

Use Party 1 for MVP (A) and Party 2 for AoU (B):

```bash
# MVP: use an MVP configuration directory pointing at normalized local inputs.
go run -mod=vendor secure-rvas.go prepare --config config/mvp --party 1

# Copy MVP's public variant lists and set public_var_list in config/aou/configPrepare.toml.
go run -mod=vendor secure-rvas.go prepare --config config/aou --party 2
```

For B, set the input list path in `configPrepare.toml`:

```toml
public_var_list = "/path/mvp-public/chr{chromosome}/public_variants.tsv"
```

`--public-var-list PATH` overrides this setting. Party 2 requires a list;
other prepare modes require no list setting. Each site needs only the global,
prepare, and its own local party configuration. Party 0 needs only the global
and Party 0 local configuration. Generate matching keys once with all three
configuration directories:

```bash
go run -mod=vendor secure-rvas.go keygen --config-party0 config/cp0 --config-party1 config/mvp --config-party2 config/aou
```

Start these commands concurrently on their corresponding sites:

```bash
go run -mod=vendor secure-rvas.go party --config config/cp0 --party 0
go run -mod=vendor secure-rvas.go party --config config/mvp --party 1
go run -mod=vendor secure-rvas.go party --config config/aou --party 2
```

The complete AoU workflow and single-directory examples below describe the
original split/local layout. They require all three Local files and a prepare
configuration without `public_var_list`; use the separate-site commands above
with the current `config/aou`.
The single-directory `keygen --config` and `run --config` commands continue to
support local configurations containing all three parties.

`config/mvp` is the Party 1 template. Update its input paths and column names
for normalized MVP data. `config/aou` is Party 2, and `config/cp0` is Party 0.
The AoU configuration receives copied MVP lists under `data/mvp/public/chrN/`;
change `public_var_list` to their actual location.
Both sites must use inputs on the same genome build with biallelic variant IDs
in `chromosome:position:REF:ALT` form and corresponding gene IDs. Genome-build
alignment is handled during data preparation. REF/ALT and coordinates are
checked against PVAR during extraction.

With `--party 1`, all selected variants become A's public columns and only `A/`
is generated. The ordered list is saved at
`<run_dir>/prepared/<ancestry>/chr<chromosome>/public_variants.tsv`, beside
`block_sizes.txt`. No combined list is created at the run root. It includes empty genes so both
sites retain identical gene/block order.

With `--party 2`, the list fixes gene order and public column order. Public
variants present in the local PGEN are read regardless of local annotation/MAF
selection; absent public variants are zero-filled. Locally selected variants
outside the public list become `B/private`, within the genes listed by MVP.
Only `B/` is generated. A present variant that cannot be extracted as biallelic
causes an error rather than silently becoming an absent variant. Missing calls
retain the existing preprocessing behavior (zero dosage).

For either site, `samples_per_cohort = 0` uses all eligible local samples;
a positive value caps that site's sample count. `shared_rate` and random variant
roles are unused. Party 2 takes gene selection from the public list, and both
sites must configure the same chromosomes. Prepared-cache keys include the
party and public-list content. Use separate run directories for the two sites
and keep the received list outside generated directories when using `--clear`.
`--public-var-list` accepts a `{chromosome}` path template or one combined TSV.
The public list has no sample-specific fields, so one ancestry's lists can be shared.
Lists inside `prepared/` are removed with that directory by `--clear`.
Omitting `--party` retains the original A/B split mode.

### 3. Generate shared PRG keys

```bash
go run -mod=vendor secure-rvas.go keygen \
  --config config/aou
```

`keygen` writes each party's key subset under its configured
`shared_keys_path/<ancestry>/`. It is independent of preprocessing.

### 4. Run secure Burden and SKAT

```bash
go run -mod=vendor secure-rvas.go run \
  --config config/aou
```

The parent only starts the same `party --party 0/1/2` commands that can be run
directly on separate machines. Input loading and all network/MHE/MPC setup happen
inside each party process. `<run_dir>/secure/_SUCCESS` is created only after
every party finishes successfully.

For a distributed run, skip the parent command and start the following on the
three machines at the same time:

```bash
go run -mod=vendor secure-rvas.go party --config config/aou --party 0
go run -mod=vendor secure-rvas.go party --config config/aou --party 1
go run -mod=vendor secure-rvas.go party --config config/aou --party 2
```

An `EOF` or `connection reset by peer` generally means another party or MPC
lane exited first. Inspect the earliest error in the complete log rather than
treating the later network panic as the root cause. On a memory-constrained VM,
retry chromosome 22 with `mpc_num_threads = 2`; use `1` to isolate a
multi-lane-only failure.

### 5. Run the reference

```bash
source scripts/setup/env.sh

python3 python/analysis/run_reference.py \
  --config config/aou \
  --engine python
```

Use `--engine r` for external R::SKAT. `REFERENCE_WORKERS` controls reference
task parallelism and `REFERENCE_BLAS_THREADS` controls BLAS threads per task.

### 6. Compare results

```bash
python3 python/analysis/compare_secure_to_reference.py \
  --config config/aou
```

The comparison joins chromosome, gene, and phenotype identities and retains raw
p-values and absolute errors. Secure Wilson-Hilferty SKAT is compared primarily
with reference SKAT-Liu and R::SKAT's Davies result across every row. The Davies
result follows R::SKAT 2.2.5: rank-one and invalid Davies results use the
modified-Liu fallback. Summaries report both all-row and non-failed R-squared;
empty genes have p-value 1 and are not counted as failures. Non-failed results
exclude only rows where R::SKAT reports `Is_Converged=0`.

### 7. Generate plots and summarize metrics

```bash
python3 python/analysis/plot_secure_vs_reference.py \
  --config config/aou

./scripts/summarize_metrics.sh config/aou
```


## Output layout

The main derived outputs are written below `run_dir`:

```text
<run_dir>/
├── prepared/<ancestry>/chr<chromosome>/
├── secure/<ancestry>/
│   ├── chr<chromosome>.tsv
│   └── all_secure_results.tsv
├── reference/<ancestry>/
├── comparison/<ancestry>/
│   └── plots/
└── metrics/<ancestry>/
```

Timing and communication metrics are separated by ancestry and party. Stage
timings can overlap under parallel execution and should not be added to
reconstruct wall-clock time.

## Repository layout

```text
crypto/                         Lattigo v6 MHE backend
mpc/                            Network and MPC backend
rvas/protocol/                  Secure Burden/SKAT protocol
rvas/workflow/                  secure-rvas orchestration
python/preprocessing/           A/B preprocessing
python/analysis/                Reference, comparison, and plotting
python/datasets/aou/            AoU localization and VAT tools
python/datasets/onekg/          Local public-data preparation
config/aou/                     AoU prepare/global/local configurations
config/1kg/                     Local 1000 Genomes configurations
scripts/                        Workflow, metrics, and setup scripts
data/                           Generated inputs and prepared caches (not tracked)
output/                         Run results (not tracked)
```
