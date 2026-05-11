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
- `make_1000g_synthetic_pheno_table.py`: Small helper to create a global 1000G test phenotype/covariate table.

## 3. 1000 Genomes Example

Download the raw 1000 Genomes chr22 VCF and panel file.

```bash
export G1000_DIR="$PWD/datasets/1000g_all_chr22_test_260507"

mkdir -p "$G1000_DIR/source"
cd "$G1000_DIR/source"

curl -L -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
curl -L -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
curl -L -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

cd -
```

Create one global synthetic phenotype/covariate table from the panel. No party split files are needed.

```bash
python3 scripts/preprocessing/skat/make_1000g_synthetic_pheno_table.py \
  --panel-file "$G1000_DIR/source/integrated_call_samples_v3.20130502.ALL.panel" \
  --out "$G1000_DIR/pheno/pheno_cov.tsv"
```

Generate local shared keys for the three-process example run.

```bash
mkdir -p "$G1000_DIR/keys"
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_global.bin" 32
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_0_1.bin" 32
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_0_2.bin" 32
python3 scripts/generateRandomKey.py "$G1000_DIR/keys/shared_key_1_2.bin" 32
```

Build the secure-skat dataset. The script converts VCF to PGEN internally, then creates the party1/party2 split under `$G1000_DIR/dataset`.

```bash
python3 scripts/preprocessing/skat/build_pgen_window_dataset.py \
  --chromosome 22 \
  --vcf "$G1000_DIR/source/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
  --raw-dir "$G1000_DIR/raw_cache" \
  --work-dir "$G1000_DIR/work" \
  --pheno-file "$G1000_DIR/pheno/pheno_cov.tsv" \
  --id-col sample_id \
  --pheno-col pheno \
  --cov-cols superpop_AFR,superpop_AMR,superpop_EAS,superpop_EUR \
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
$G1000_DIR/
  source/      raw VCF and panel
  pheno/       global phenotype/covariate table
  keys/        generated local shared keys
  raw_cache/   VCF-to-PGEN cache
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
  --pheno-file "$HOME/gwas-data-wgs-wb-jaunty-blueberry-8679/pheno/YOUR_PHENO.tsv" \
  --id-col person_id \
  --pheno-col YOUR_TRAIT \
  --cov-cols age,sex,PC1,PC2,PC3,PC4,PC5 \
  --out-dataset "$HOME/secure-skat/dataset/aou_exome_chr22_windows" \
  --config-out-dir "$HOME/secure-skat/config_aou_exome_chr22_windows" \
  --billing-project "$GOOGLE_PROJECT" \
  --backup-gcs-uri "${WORKSPACE_BUCKET}/dataset/aou_exome_chr22_windows" \
  --n-samples 2000 \
  --maf-threshold 0.01 \
  --window-bp 50000 \
  --min-rare-per-window 20 \
  --max-windows 8 \
  --force
```
