# SKAT Window Dataset Builder

## 1. Intro

This script builds a secure-skat-friendly party1/party2 window dataset from either a raw VCF or an existing PLINK2 PGEN dataset. It can select rare variants by MAF, group them into fixed genomic windows, materialize per-party PGEN blocks, and emit matching config files.

## 2. File Layout

- `build_pgen_window_dataset.py`: Main entrypoint and pipeline orchestration.
- `argument.py`: CLI arguments and validation.
- `plink.py`: VCF/PGEN input handling and PLINK2 commands.
- `samples.py`: PSAM parsing, phenotype/covariate matching, and party splitting.
- `windowing.py`: Allele-frequency loading and rare-variant window selection.
- `config.py`: Config TOML generation, manifest writing, and optional GCS backup.
- `cloud.py`: Small `gcloud storage` / `gsutil` helpers.
- `utils.py`: Shared path, subprocess, and parsing helpers.
- `prepare_1000g_source.py`: Prepare reusable 1000G PGEN, phenotype, and covariate inputs.
- `make_1000g_synthetic_pheno_table.py`: Commented legacy marker; slated for deletion after migration.

## 3. 1000 Genomes Example

Keep reusable raw 1000G files in one source directory, then build each experiment in its own output directory.

```bash
export G1000_SOURCE="$PWD/datasets/1000g_source"
export G1000_DIR="$PWD/datasets/1000g_all_chr22_test_260511"

mkdir -p "$G1000_SOURCE/raw"
cd "$G1000_SOURCE/raw"

curl -L -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
curl -L -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
curl -L -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

cd -
```

Prepare AoU-style reusable inputs from the 1000G source files.

```bash
python3 scripts/preprocessing/skat/prepare_1000g_source.py \
  --source-dir "$G1000_SOURCE" \
  --chromosome 22 \
  --force
```

Generate local shared keys for the three-process example run.

```bash
mkdir -p "$G1000_DIR/keys"
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_global.bin" 32
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_0_1.bin" 32
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_0_2.bin" 32
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_1_2.bin" 32
```

Build the secure-skat dataset from the prepared PGEN trio and split phenotype/covariate CSV files.

```bash
python3 scripts/preprocessing/skat/build_pgen_window_dataset.py \
  --chromosome 22 \
  --pgen-prefix "$G1000_SOURCE/pgen/1000g.chr22" \
  --raw-dir "$G1000_SOURCE/pgen" \
  --work-dir "$G1000_DIR/work" \
  --pheno-file "$G1000_SOURCE/pheno/phenotype_data.csv" \
  --cov-file "$G1000_SOURCE/pheno/covariate_data.csv" \
  --pheno-col pheno \
  --out-dataset "$G1000_DIR/dataset" \
  --config-template-dir config \
  --config-out-dir "$G1000_DIR/config" \
  --shared-keys-path "$G1000_DIR/keys" \
  --n-samples 120 \
  --maf-threshold 0.05 \
  --window-bp 50000 \
  --min-rare-per-window 2 \
  --max-windows 2 \
  --force
```

Run secure SKAT with all run outputs under `$G1000_DIR/runs`.

```bash
bash run_example.sh \
  --mode skat \
  --dataset "$G1000_DIR/dataset" \
  --config-dir "$G1000_DIR/config" \
  --run-base "$G1000_DIR/runs"
```

After these steps, the main walkthrough directory looks like this.

```text
$G1000_SOURCE/
  raw/         raw VCF, index, and panel
  pgen/        reusable 1000G PGEN trio
  pheno/       reusable global phenotype and covariate tables

$G1000_DIR/
  keys/        generated local shared keys
  work/        preprocessing scratch files
  dataset/     secure-skat input with party1/ and party2/
  config/      generated TOML config files
  runs/        run_example.sh outputs, logs, and protocol caches
```

## 4. AoU Example

```bash
python3 scripts/preprocessing/skat/build_pgen_window_dataset.py \
  --chromosome 22 \
  --pgen-prefix "gs://vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/pgen/exome.chr22" \
  --raw-dir "$HOME/datasets/aou_exome_chr22_test/raw_cache" \
  --work-dir "$HOME/datasets/aou_exome_chr22_test/work" \
  --pheno-file "$HOME/datasets/aou_exome_chr22_test/pheno/phenotype_data.csv" \
  --cov-file "$HOME/datasets/aou_exome_chr22_test/pheno/covariate_data.csv" \
  --pheno-col HDLC_value \
  --out-dataset "$HOME/secure-skat/dataset/aou_exome_chr22_windows" \
  --config-template-dir config \
  --config-out-dir "$HOME/secure-skat/config_aou_exome_chr22_windows" \
  --shared-keys-path "$HOME/datasets/aou_exome_chr22_test/keys" \
  --billing-project "$GOOGLE_PROJECT" \
  --backup-gcs-uri "${WORKSPACE_BUCKET}/dataset/aou_exome_chr22_windows" \
  --n-samples 500 \
  --maf-threshold 0.01 \
  --window-bp 50000 \
  --min-rare-per-window 20 \
  --max-windows 2 \
  --force
```

Split phenotype/covariate mode uses the first column in both CSV files as the sample ID.
Use `--pheno-col <name>` or `--pheno-col-index <1-based-index>` to choose the phenotype.
Use `--cov-cols <names>` or `--cov-col-indices <1-based-indices>` to choose covariates; if both are omitted, every covariate column except the first ID column is used. Covariate values must be numeric.
