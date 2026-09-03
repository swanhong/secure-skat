# Local 1000 Genomes workflow

This workflow runs the same secure Burden/SKAT pipeline as the repository root
README using public 1000 Genomes data and a synthetic phenotype/annotation
fixture. It is intended for local correctness and regression runs, not for
biological interpretation.

Run every command below from the repository root.

## Initial setup

```bash
bash setup/install.sh
source setup/env.sh
go mod vendor
```

## Configure the 1KG run

Review `run.1kg.conf` before running. The checked-in configuration uses
chromosomes 21 and 22, two deterministic test phenotypes, and EUR/AFR/AMR:

```toml
run_dir = "output/pca-test"
chromosomes = [21, 22]
phenotype_columns = ["phenotype1", "phenotype2"]
ancestries = ["EUR", "AFR", "AMR"]
```

Source generation, secure preprocessing, protocol execution, reference
analysis, and metrics summarization all read their run parameters from this
same configuration. To run other autosomes, replace `chromosomes` with the
desired non-empty subset of integers from 1 through 22.

The current fixture selects synthetic high-confidence loss-of-function variants
with MAF at most 0.01:

```toml
masks = ["LoF=HC"]
max_maf = 0.01
```

`binding_ipaddr` is the local listener address. Each
`[servers.partyN].ipaddr` is the address other parties use to reach that party.
The checked-in `127.0.0.1` values assume all three parties run on one host.

## Run the complete 1KG workflow

```bash
./run_1kg_workflow.sh
```

The script performs the following stages:

```text
0. Generate 1000 Genomes source-format test data
1. Preprocess ancestry-specific A/B secure inputs
2. Run secure Burden and SKAT
3. Run the Python or R::SKAT reference
4. Join and compare secure/reference results
5. Generate scatter and Manhattan plots
6. Summarize timing, communication, and accuracy metrics
```

Downloaded and generated source files are reused. By default, preprocessing
does not replace an existing `run_dir`. Allow replacement explicitly with:

```bash
CLEAR_RUN_DIR=1 ./run_1kg_workflow.sh
```

Select external R::SKAT or another configuration with environment variables:

```bash
CONFIG_PATH=run.1kg.conf REFERENCE_ENGINE=r ./run_1kg_workflow.sh
```

### Detached execution and monitoring

```bash
log_path="run-1kg-$(date +'%m%d-%H%M').log"
nohup env PYTHONUNBUFFERED=1 ./run_1kg_workflow.sh \
  > "$log_path" 2>&1 < /dev/null &
echo $! > "${log_path%.log}.pid"
echo "log: $log_path"
```

```bash
tail -f "$log_path"
```

## Run 1KG stages individually

### 1. Generate source-format inputs

```bash
source setup/env.sh

python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --config run.1kg.conf \
  --num-pheno 2
```

The generator reads `chromosomes` from `run.1kg.conf`, downloads missing 1000
Genomes high-coverage phased VCFs, the population panel, and the GENCODE v50
GTF, then writes:

```text
rewrite/testdata/1kgenome/generated/
├── raw/
├── genotype/
│   ├── chr<chromosome>.pgen
│   ├── chr<chromosome>.pvar
│   ├── chr<chromosome>.psam
│   └── chr<chromosome>.afreq
├── gene_panel/chr<chromosome>.tsv
├── annotation/chr<chromosome>.tsv
├── phenotype.csv
├── ancestry_pred.tsv
└── work/
```

The generated data contract is:

- PVAR IDs use `chromosome:position:REF:ALT`;
- the gene panel contains GENCODE v50 protein-coding genes in genomic order;
- every requested chromosome has the same ordered PSAM sample IDs;
- PLINK 2 allele frequencies provide numeric `MAF` values;
- seed-fixed Gaussian phenotypes are named `phenotype1`, `phenotype2`, and so
  on;
- LD-pruned common variants provide 16 PLINK 2 principal components;
- `ancestry_pred.tsv` combines the PCs with 1000 Genomes super-population
  labels in the AoU ancestry-table shape.

The annotation is a synthetic protocol fixture. Variant/gene overlaps and
allele frequencies come from the generated data, but the `LoF` and consequence
labels are synthetic. Approximately one percent of annotation rows receive
`LoF=HC`; the rest receive `LoF=LC`. Generated files are excluded from Git.

### 2. Preprocess secure inputs

```bash
go run -mod=vendor secure-rvas.go prepare \
  --config run.1kg.conf
```

Add `--clear` only when the configured `run_dir` may be replaced.

### 3. Run secure Burden and SKAT

```bash
go run -mod=vendor secure-rvas.go run \
  --config run.1kg.conf
```

The parent starts parties 0, 1, and 2. `<run_dir>/secure/_SUCCESS` is created
only after every configured ancestry and chromosome finishes successfully.

### 4. Run the reference

```bash
source setup/env.sh

python3 rewrite/analysis/run_reference.py \
  --config run.1kg.conf \
  --engine python
```

Use `--engine r` for external R::SKAT.

### 5. Compare results

```bash
python3 rewrite/analysis/compare_secure_to_reference.py \
  --config run.1kg.conf
```

### 6. Generate plots and summarize metrics

```bash
python3 rewrite/analysis/plot_secure_vs_reference.py \
  --config run.1kg.conf

./summarize_metrics.sh run.1kg.conf
```

Pass a party ID as the second summary argument, or override the wrapper's
default `python` interpreter:

```bash
PYTHON_BIN="$PWD/.venv/bin/python" \
  ./summarize_metrics.sh run.1kg.conf 0
```

The Python summary reads `run_dir` and `ancestries` from `run.1kg.conf` and
discovers the corresponding metrics and comparison CSV files automatically.

## Output layout

Protocol and analysis outputs are written below the configured `run_dir`:

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
timings can overlap during parallel execution and should not be added to
reconstruct wall-clock time.
