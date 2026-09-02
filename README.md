# Secure RVAS

Secure RVAS runs privacy-preserving Burden and SKAT rare-variant association tests.
The local 1000 Genomes test workflow is documented separately in [`rewrite/testdata/1kgenome/README.md`](rewrite/testdata/1kgenome/README.md).

## Initial setup

Run these commands once after a fresh clone in a Linux x86-64 environment:

```bash
cd "$HOME/secure-skat"
bash setup/install.sh
source setup/env.sh
go mod vendor
```

`setup/install.sh` installs PLINK 2 at `$HOME/plink2` and R::SKAT under `$HOME/R/library`. `go mod vendor` is required because every Go command in this repository uses `-mod=vendor`, while `vendor/` is not stored in Git.

In each new terminal, load the installed paths before running individual commands:

```bash
cd secure-skat
source setup/env.sh
```

## Configure the AoU run

Review `run.aou.conf` before running. Update `run_dir` and the local input paths. `mpc_num_threads` is the number of independent MPC lanes, not just a CPU thread count. More lanes increase concurrent memory and communication use.

All stages read the same non-empty autosome array from the configuration:

```toml
chromosomes = [
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
]
```

## Step 0: Create chromosome VAT annotations

`vat_simplify.py` reads the chromosome PVAR from Google Cloud Storage without
first saving a separate local copy. The `--pvar` argument itself accepts a local
file descriptor, so Bash process substitution connects `gsutil cat` to the
Python command.

```bash
cd "$HOME/secure-skat"
source setup/env.sh

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

  python3 rewrite/testdata/aou/vat_simplify.py \
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
./run_aou_workflow.sh
```

The script performs the following stages:

```text
0. Localize AoU PGEN, phenotype, and ancestry inputs; normalize VAT annotations
1. Preprocess ancestry-specific A/B secure inputs
2. Run secure Burden and SKAT
3. Run the Python or R::SKAT reference
4. Join and compare secure/reference results
5. Generate scatter and Manhattan plots
6. Summarize timing, communication, and accuracy metrics
```

Step 0 inside `run_aou_workflow.sh` is distinct from the VAT extraction above:
it downloads the PGEN triplets and other AoU inputs, then consumes
`$HOME/fed_prep_out/chr{chromosome}_annotation.tsv` to create normalized local
annotations and gene panels.

The workflow's preprocessing command includes `--clear`; rerunning the complete
script replaces the configured `run_dir` before creating new secure inputs.

To select another configuration or the R::SKAT reference engine:

```bash
CONFIG_PATH=run.aou.conf REFERENCE_ENGINE=r ./run_aou_workflow.sh
```

### Detached execution and monitoring

Use a date-and-time log path and print it before detaching:

```bash
log_path="run-aou-$(date +'%m%d-%H%M').log"
nohup env PYTHONUNBUFFERED=1 ./run_aou_workflow.sh \
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
source setup/env.sh

chromosome_values="$(
  python3 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as config_file:
    print(*tomllib.load(config_file)["chromosomes"])
' run.aou.conf
)"
read -r -a chromosomes <<< "$chromosome_values"

python3 rewrite/testdata/aou/prepare_aou.py \
  --chromosome "${chromosomes[@]}"
```

`prepare_aou.py` downloads missing PGEN/PVAR/PSAM, phenotype, and ancestry
files. Existing localized files are reused. It normalizes the Step 0 VAT output,
computes minor allele frequency from `gnomad_af`, and writes inputs under
`rewrite/testdata/aou/generated/`.

### 2. Preprocess secure inputs

```bash
go run -mod=vendor secure-rvas.go prepare \
  --config run.aou.conf \
  --clear
```

Omit `--clear` when the configured `run_dir` must not be replaced.

### 3. Run secure Burden and SKAT

```bash
go run -mod=vendor secure-rvas.go run \
  --config run.aou.conf
```

The parent starts parties 0, 1, and 2 and processes each configured ancestry
and chromosome. `<run_dir>/secure/_SUCCESS` is created only after every party
finishes successfully.

An `EOF` or `connection reset by peer` generally means another party or MPC
lane exited first. Inspect the earliest error in the complete log rather than
treating the later network panic as the root cause. On a memory-constrained VM,
retry chromosome 22 with `mpc_num_threads = 2`; use `1` to isolate a
multi-lane-only failure.

### 4. Run the reference

```bash
source setup/env.sh

python3 rewrite/analysis/run_reference.py \
  --config run.aou.conf \
  --engine python
```

Use `--engine r` for external R::SKAT. `REFERENCE_WORKERS` controls reference
task parallelism and `REFERENCE_BLAS_THREADS` controls BLAS threads per task.

### 5. Compare results

```bash
python3 rewrite/analysis/compare_secure_to_reference.py \
  --config run.aou.conf
```

The comparison joins chromosome, gene, and phenotype identities and retains raw
p-values and absolute errors. Secure Wilson-Hilferty SKAT is compared primarily
with reference SKAT-Liu and with SKAT-Davies when Davies converged.

### 6. Generate plots and summarize metrics

```bash
python3 rewrite/analysis/plot_secure_vs_reference.py \
  --config run.aou.conf

./summarize_metrics.sh run.aou.conf
```

Pass a party ID as the second summary argument to inspect another party:

```bash
./summarize_metrics.sh run.aou.conf 0
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
rewrite/preprocessing/          A/B preprocessing
rewrite/protocol/               Secure Burden/SKAT protocol
rewrite/workflow/               secure-rvas orchestration
rewrite/analysis/               Reference, comparison, and plotting
rewrite/testdata/aou/            AoU localization and VAT tools
rewrite/testdata/1kgenome/       Local public-data test workflow
run.aou.conf                    AoU configuration
run.1kg.conf                    Local 1000 Genomes configuration
```

## Development principles

- Use `-mod=vendor` for every Go build, test, and run command.
- `rewrite/protocol/` directly reuses primitives from `crypto/` and `mpc/`.
- Secure protocol computation remains separate from workflow orchestration and
  reference analysis.
- Do not commit AoU Controlled Tier inputs or VAT-derived outputs.

## License and attribution

The MHE and MPC backend is derived from the Lattigo v6-based SF-GWAS
implementation. See [`LICENSE`](LICENSE) for licensing information.
