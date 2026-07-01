#!/bin/bash
# One-time setup for a fresh AoU workbench node: download plink2 + fetch the phenotype CSV.
# Run once per new node, then use run_fed.sh.
#
#   FED_PHENO_SRC=gs://<your-workspace-bucket>/pheno/<lipids>.csv bash setup_fed.sh
#
# plink2 is a public build; the phenotype CSV lives in your controlled workspace bucket, so its
# path is passed via FED_PHENO_SRC (not hardcoded). GOOGLE_PROJECT is the requester-pays project.
set -e

# --- plink2 (public build) ---
if [ ! -x "$HOME/plink2" ]; then
  echo "=== plink2 ==="
  ( cd "$HOME" && wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip \
      && unzip -o plink2_linux_x86_64_latest.zip plink2 && chmod +x "$HOME/plink2" )
fi
"$HOME/plink2" --version | head -1

# --- phenotype CSV -> ~/fed_prep_in/pheno.csv (controlled bucket, set FED_PHENO_SRC) ---
mkdir -p "$HOME/fed_prep_in"
if [ -n "$FED_PHENO_SRC" ] && [ ! -f "$HOME/fed_prep_in/pheno.csv" ]; then
  echo "=== pheno.csv <- $FED_PHENO_SRC ==="
  gsutil -u "${GOOGLE_PROJECT:?set GOOGLE_PROJECT}" cp "$FED_PHENO_SRC" "$HOME/fed_prep_in/pheno.csv"
fi
[ -f "$HOME/fed_prep_in/pheno.csv" ] && echo "pheno.csv ready" || echo "WARN: no pheno.csv (set FED_PHENO_SRC)"

# --- AoU Controlled-Tier mounts (standard layout) ---
echo "=== AoU mounts ==="
ls ~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/pgen/exome.chr22.pvar \
   ~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/echo_v4_r2.ancestry_preds.tsv

echo "=== setup done -> PLINK2=\$HOME/plink2 ./run_fed.sh ==="
