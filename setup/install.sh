#!/usr/bin/env bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
  echo "This installer requires Linux x86-64." >&2
  exit 1
fi

if [[ -x "$PLINK2" ]]; then
  echo "PLINK 2 already exists at $PLINK2"
else
  download_url="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260818.zip"
  work_dir="$(mktemp -d "${TMPDIR:-/tmp}/secure-rvas-setup.XXXXXX")"
  trap 'rm -rf -- "$work_dir"' EXIT

  echo "Downloading PLINK 2"
  curl --fail --location --retry 3 \
    --output "$work_dir/plink2.zip" \
    "$download_url"
  mkdir "$work_dir/unpacked"
  python3 -m zipfile -e "$work_dir/plink2.zip" "$work_dir/unpacked"
  install -m 0755 "$work_dir/unpacked/plink2" "$PLINK2"
fi

mkdir -p "$R_LIBS_USER"
Rscript - <<'RS'
library_path <- Sys.getenv("R_LIBS_USER")
.libPaths(c(library_path, .libPaths()))

if (!requireNamespace("SKAT", quietly = TRUE)) {
  install.packages(
    "SKAT",
    lib = library_path,
    repos = "https://cloud.r-project.org"
  )
}

cat("SKAT ", as.character(packageVersion("SKAT")), "\n", sep = "")
RS

"$PLINK2" --version
echo "Setup complete"
