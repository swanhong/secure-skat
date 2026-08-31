# Secure RVAS

## Initial setup

Run this once in each new Linux x86-64 environment:

```bash
bash setup/install.sh
```

This installs PLINK 2 at `$HOME/plink2` and R::SKAT under
`$HOME/R/library`. The workflow scripts load these paths automatically. Before
running an individual command in a new terminal, load the same environment:

```bash
source setup/env.sh
```

## Quick command reference

```bash
./run_1kg_workflow.sh
```

The script runs every local 1KG step through ancestry-specific plots. To replace
an existing configured `run_dir`, opt in to preprocessing cleanup:

```bash
CLEAR_RUN_DIR=1 ./run_1kg_workflow.sh
```

The equivalent individual commands are:

```bash
source setup/env.sh

python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome 21 22 \
  --num-pheno 2

go run -mod=vendor secure-rvas.go prepare \
  --config run.1kg.conf

go run -mod=vendor secure-rvas.go run \
  --config run.1kg.conf

python3 rewrite/analysis/run_reference.py \
  --config run.1kg.conf

python3 rewrite/analysis/compare_results.py \
  --config run.1kg.conf

python3 rewrite/analysis/plot_results.py \
  --config run.1kg.conf
```

To rerun preprocessing from a clean `run_dir`, add `--clear`:

```bash
go run -mod=vendor secure-rvas.go prepare \
  --config run.1kg.conf \
  --clear
```

Secure RVAS is a protocol for privacy-preserving Burden and SKAT rare-variant association tests.
It reuses the Lattigo v6-based MHE and MPC primitives in `crypto/` and `mpc/`, while the new preprocessing and secure protocol are implemented under `rewrite/`.

The target workflow runs directly from a local terminal or an All of Us (AoU) Researcher Workbench terminal:

```text
0. Generate local test data from 1000 Genomes
        ↓
1. Preprocess source data into A/B secure inputs
        ↓
2. Run secure Burden and SKAT
        ↓
3. Run R::SKAT on the same secure inputs
        ↓
4. Compare the secure and R results
        ↓
5. Plot the comparison results
```

## Implementation status

| Step | Tool | Status |
|---|---|---|
| 0 | `prepare_1kgenome.py` | Implemented |
| 1 | `secure-rvas prepare` | Implemented for EUR/AFR/AMR local inputs |
| 2 | `secure-rvas run` | Implemented for sequential EUR/AFR/AMR execution |
| 3 | `run_reference.py` | Implemented per ancestry |
| 4 | `compare_results.py` | Implemented per ancestry |
| 5 | `plot_results.py` | Implemented per ancestry |

The local ancestry-aware workflow is implemented end to end.

## Step 0: Generate 1000 Genomes test data

Step 0 prepares source data for local correctness testing. It is not part of
the AoU Workbench workflow.

### Requirements

- Python 3
- `curl`
- PLINK 2 installed by `setup/install.sh`, or available on `PATH`
- Sufficient disk space for the selected 1000 Genomes chromosomes and GENCODE

No additional Python packages are required.

### Basic usage

Generate chromosome 21, chromosome 22, and two phenotypes:

```bash
python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome 21 22 \
  --num-pheno 2
```

`--chromosome` accepts individual chromosomes, comma-separated lists, ranges,
and all autosomes:

```bash
# Chromosomes 1, 2, and 3
python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome 1 2 3 \
  --num-pheno 2

# Chromosomes 1 through 5
python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome 1-5 \
  --num-pheno 2

# Chromosomes 1 through 22
nohup python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome all \
  --num-pheno 2 > prepare-1kg-source.log 2>&1 &
```

`--chromosome 1,2,3` is also supported. Duplicate chromosomes, descending
ranges, and values outside 1–22 are rejected. Existing source downloads and
complete PGEN outputs are reused.

### Source data

The script downloads the following files when they are missing:

- 1000 Genomes high-coverage phased VCFs
- The 1000 Genomes population panel
- The GENCODE v50 GTF

### Outputs

```text
rewrite/testdata/1kgenome/generated/
├── raw/
│   ├── 1kGP_...chr{chromosome}....vcf.gz
│   ├── integrated_call_samples_v3.20130502.ALL.panel
│   └── gencode.v50.annotation.gtf.gz
├── genotype/
│   ├── chr{chromosome}.pgen
│   ├── chr{chromosome}.pvar
│   ├── chr{chromosome}.psam
│   └── chr{chromosome}.afreq
├── gene_panel/
│   └── chr{chromosome}.tsv
├── annotation/
│   └── chr{chromosome}.tsv
├── phenotype.csv
├── ancestry_pred.tsv
└── work/
    ├── phase3.keep
    └── pca/chr{chromosomes}/
        ├── pca.eigenvec
        └── pca.eigenval
```

The generator enforces the following contract:

- PVAR variant IDs use `chromosome:position:REF:ALT`.
- The gene panel contains GENCODE v50 protein-coding genes in genomic order.
- The run fails if the selected chromosomes do not have identical ordered PSAM
  IDs.
- PLINK2 allele frequencies are used to write a numeric `MAF` annotation for
  every variant-gene row.
- `phenotype1`, `phenotype2`, and subsequent columns contain seed-fixed
  Gaussian test phenotypes.
- The requested chromosome PGENs are merged, filtered with `--maf 0.05` and
  `--geno 0.02`, LD-pruned, and used for a 16-component PLINK2 PCA.
- `ancestry_pred.tsv` combines those PC scores with 1000 Genomes super-population
  labels in the five-column AoU ancestry-table format.

The annotation is a synthetic protocol-testing fixture. It uses overlaps
between PVAR variants and GENCODE genes and real frequencies computed from the
generated PGEN, but its `LoF` and consequence labels are synthetic. Approximately
one percent of annotation rows receive `LoF=HC`; the rest receive `LoF=LC`.
These labels must not be interpreted as research annotations. Raw and generated
files under `generated/` are excluded from Git.

## Step 1: Preprocess secure inputs

Step 1 converts the Step 0 source data into A/B secure inputs:

```bash
go run -mod=vendor secure-rvas.go prepare \
  --config run.1kg.conf
```

`secure-rvas prepare` does not download 1000 Genomes data. It consumes the
Step 0 outputs, or user-provided local data that follows the same source-input
contract. The command reads `run.1kg.conf` by default; another configuration can
be selected with `--config`.

The 1KG configuration reads the first 16 values from the `pca_features` array
in `ancestry_pred.tsv` as numeric covariates:

```toml
covariate = "rewrite/testdata/1kgenome/generated/ancestry_pred.tsv"
covariate_id_column = "research_id"
covariate_column = "pca_features"
num_cov = 16

ancestry = "rewrite/testdata/1kgenome/generated/ancestry_pred.tsv"
ancestry_id_column = "research_id"
ancestry_column = "ancestry_pred_other"
ancestries = ["EUR", "AFR", "AMR"]
samples_per_cohort = 0
```

The current 1KG configuration applies all equality masks and then the optional
MAF threshold:

```toml
masks = ["LoF=HC"]
max_maf = 0.01
```

This means `LoF == "HC" AND MAF <= 0.01`. Omitting `max_maf` disables frequency
filtering.

The command produces:

```text
<run_dir>/
├── selected_genes.tsv      # Present when gene selection mode is not `all`
└── prepared/
│   ├── EUR/
│   │   ├── chr21/
│   │   │   ├── A/
│   │   │   │   ├── geno/block.<gene-index>.bin
│   │   │   │   ├── cov.txt
│   │   │   │   └── pheno.txt
│   │   │   ├── B/
│   │   │   │   ├── geno/block.<gene-index>.bin
│   │   │   │   ├── private/block.<gene-index>.bin
│   │   │   │   ├── cov.txt
│   │   │   │   └── pheno.txt
│   │   │   ├── genes.txt
│   │   │   ├── block_sizes.txt
│   │   │   └── pos.txt
│   │   └── chr22/
│   ├── AFR/
│   │   ├── chr21/
│   │   └── chr22/
│   └── AMR/
│       ├── chr21/
│       └── chr22/
```

Phenotypes, PCA covariates, and ancestry labels are loaded once. Each ancestry
uses all eligible samples when `samples_per_cohort = 0`, performs its own seeded
A/B split, and reuses that split across every configured chromosome. Gene and
annotation selection is shared across ancestries.

Gene selection supports `random`, `file`, and `all` modes. In `random` mode,
`per_chromosome` is the number selected independently for each chromosome.

## Step 2: Run secure Burden and SKAT

```bash
go run -mod=vendor secure-rvas.go run \
  --config run.1kg.conf
```

The parent creates temporary shared PRG keys, starts parties 0, 1, and 2, and
waits for all three processes. Each party processes `ancestries` in configuration
order. For every ancestry it loads the matching `prepared/<ancestry>` inputs,
opens an independent network/MPC session, runs its own `SetupNull` once, and then
processes all configured chromosomes sequentially with that ancestry's null
shares. Ancestries and chromosomes are not run in parallel, and no independent
MPC lane is used.

The command writes:

```text
<run_dir>/secure/
├── EUR/
│   ├── chr21.tsv
│   ├── chr22.tsv
│   └── all_secure_results.tsv
├── AFR/
│   └── ...
├── AMR/
│   └── ...
└── _SUCCESS
```

Timing and communication metrics are separated by ancestry:

```text
<run_dir>/metrics/
├── EUR/metrics_party0.csv
├── EUR/metrics_party1.csv
├── EUR/metrics_party2.csv
├── AFR/...
├── AMR/...
└── process_summary.csv
```

Party 1 prints one timing tree per ancestry, headed by `[party1 EUR]`,
`[party1 AFR]`, or `[party1 AMR]`. `process_summary.csv` remains the overall
parent/party process time and peak RSS summary for the complete secure run.

Summarize Party 1 timing, communication, and R² results with:

```bash
./summarize_metrics.sh run.aou.conf
```

The first argument defaults to `run.aou.conf`; the script reads `run_dir` from
that configuration. Pass a second argument to summarize another party, for
example `./summarize_metrics.sh run.aou.conf 0`.

The communication table reports setup and chromosome totals without adding
nested stages, and labels sent-plus-received traffic as that party's total I/O.
When `comparison/<ancestry>/all_comparison.csv` is available, the R² table
reports pooled scores and the worst phenotype-by-chromosome score for each
comparison. Otherwise, the timing and communication tables are still printed
and R² is marked unavailable.

The stage columns are independently measured and may overlap when work runs in
parallel, so they should not be added to reconstruct the chromosome total.

`_SUCCESS` is created only after all parties and chromosomes finish. The result
columns contain the secure Burden p-value and trace-based Wilson-Hilferty SKAT
p-value. Values are stored without clipping, so very small fixed-point or CKKS
noise can place a value slightly above one.

A zero or otherwise invalid kernel moment tuple does not enter the inverse
square-root path. It is evaluated with safe secret-shared inputs and releases
the finite sentinel `z=-9`, which maps to an SKAT p-value of one in float64.

## Step 3: Run the R::SKAT reference

Run R::SKAT on the same A/B prepared inputs consumed by the secure protocol:

```bash
source setup/env.sh

python3 rewrite/analysis/run_reference.py \
  --config run.1kg.conf
```

The wrapper runs every configured ancestry and chromosome in order and writes:

```text
<run_dir>/reference/
├── EUR/
│   ├── chr21.csv
│   ├── chr22.csv
│   └── all_r_results.csv
├── AFR/
│   └── ...
├── AMR/
│   └── ...
└── _SUCCESS
```

For every gene and phenotype it records R Burden, SKAT-Liu, and SKAT-Davies
p-values. `r_skat_davies_converged=1` means Davies converged, `0` means Davies
returned a nonzero failure status, and `NA` means R::SKAT did not run Davies
because no variant remained testable. R returns `p=1` for that degenerate case.
Warnings about monomorphic variants are expected for the small 20-sample
fixture and do not indicate command failure.

## Step 4: Compare secure and R results

```bash
python3 rewrite/analysis/compare_results.py \
  --config run.1kg.conf
```

The comparison requires a one-to-one match on chromosome, gene index, gene ID,
phenotype index, and phenotype name. It compares:

- secure Burden against R Burden;
- secure Wilson-Hilferty SKAT against R SKAT-Liu as the primary approximation
  comparison;
- secure Wilson-Hilferty SKAT against R SKAT-Davies only when Davies converged.

It writes:

```text
<run_dir>/comparison/
├── EUR/all_comparison.csv
├── AFR/all_comparison.csv
├── AMR/all_comparison.csv
└── _SUCCESS
```

The CSV retains all raw p-values and absolute errors. A missing Davies error is
preserved as an empty field rather than converted into an ordinary p-value.

## Step 5: Generate plots

```bash
python3 rewrite/analysis/plot_results.py \
  --config run.1kg.conf
```

Each ancestry receives chromosome-level Burden and SKAT scatter plots and
phenotype-level Manhattan plots across all configured chromosomes:

```text
<run_dir>/comparison/<ancestry>/plots/
├── scatter_burden_pheno<q>_chr<c>.png
├── scatter_skat_liu_pheno<q>_chr<c>.png
├── manhattan_burden_pheno<q>.png
├── manhattan_skat_liu_pheno<q>.png
└── _SUCCESS
```

## Complete local workflow

Run the complete sequence with:

```bash
./run_1kg_workflow.sh
```

Its contents are equivalent to:

```bash
set -euo pipefail

python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome 21 22 \
  --num-pheno 2

go run -mod=vendor secure-rvas.go prepare \
  --config run.1kg.conf

go run -mod=vendor secure-rvas.go run \
  --config run.1kg.conf

python3 rewrite/analysis/run_reference.py \
  --config run.1kg.conf

python3 rewrite/analysis/compare_results.py \
  --config run.1kg.conf

python3 rewrite/analysis/plot_results.py \
  --config run.1kg.conf
```

For a detached Workbench-terminal run, apply `nohup` to the long secure step:

```bash
nohup go run -mod=vendor secure-rvas.go run \
  --config run.1kg.conf > run-1kg.log 2>&1 &
```

Run the R reference, comparison, and plots after `secure/_SUCCESS` appears.

## Previous pooled end-to-end validation

The secure and R workflow was validated on chromosomes 21 and 22 with 20 genes
per chromosome, two phenotypes, and 1,200 samples per cohort, for 80 comparison
rows.

The current run reported:

- Burden mean/maximum absolute difference: `3.34e-6` / `1.74e-5`;
- secure WH versus R-Liu mean/maximum absolute difference: `0.00340` / `0.03733`;
- phenotype/chromosome-level secure WH versus R-Liu R-squared:
  `0.999086`–`0.999903`;
- no nondegenerate secure SKAT result incorrectly fell back to `p=1`.

The `W > R` Hutchinson path was also validated with a 245-variant gene. Its two
secure WH versus R-Liu absolute differences were `0.00206` and `0.00104`.

Repeated Step 0 generation with the same inputs and seeds produced identical
gene-panel, annotation, phenotype, and covariate file hashes. Generated local
outputs live under `output/` and are excluded from Git.

AoU input localization/normalization, workflow-level timing, resume support,
chromosome parallelism, and independent MPC lanes remain later work.

## Repository layout

```text
crypto/                         Lattigo v6 MHE backend
mpc/                            Network and MPC backend
rewrite/preprocessing/          New A/B preprocessing library
rewrite/protocol/               New secure Burden/SKAT protocol
rewrite/workflow/               secure-rvas application workflow
rewrite/analysis/               R reference and result comparison
rewrite/testdata/1kgenome/       Local 1000 Genomes source generator
run.1kg.conf                    Local configuration
run.aou.conf                    Draft AoU configuration; not yet compatible
```

## Development principles

- Use `-mod=vendor` for every Go build, test, and run command.
- `rewrite/protocol/` directly reuses primitives from `crypto/` and `mpc/`.
- Secure protocol computation, including `Finalize` and `Release`, remains
  separate from workflow orchestration.
- The runner uses one MPC session per ancestry and processes ancestries and
  chromosomes sequentially. It does not use chromosome workers or independent
  MPC lanes.
- AoU-specific conversion, workflow-level timing, resume support, and parallel
  MPC lanes remain separate later work.

## License and attribution

The MHE and MPC backend in this repository is derived from the Lattigo v6-based
SF-GWAS implementation. See [`LICENSE`](LICENSE) for licensing information.
