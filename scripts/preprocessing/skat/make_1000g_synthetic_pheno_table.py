#!/usr/bin/env python3
#
# Legacy script kept only as a deletion marker.
#
# The old make_1000g_synthetic_pheno_table.py workflow has been superseded by:
#
#   scripts/preprocessing/skat/prepare_1000g_source.py
#
# Important helper code was moved there so the current 1000G flow can prepare:
#
#   1. reusable PLINK2 PGEN files from the raw 1000G VCF
#   2. split phenotype_data.csv
#   3. split covariate_data.csv
#
# This file is intentionally commented out and should be deleted once no
# branches/scripts reference the legacy entrypoint.
