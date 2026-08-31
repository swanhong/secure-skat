export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
export LD_LIBRARY_PATH="$(Rscript -e 'cat(R.home("lib"))')${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PLINK2="${PLINK2:-$HOME/plink2}"
export PATH="$HOME:$PATH"
