# Secure RVAS

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
```

## Implementation status

| Step | Tool | Status |
|---|---|---|
| 0 | `prepare_1kgenome.py` | Implemented |
| 1 | `secure-rvas prepare` | Implemented for local filesystem inputs |
| 2 | `secure-rvas run` | Implemented and validated with three local parties |
| 3 | `run_reference.py` | Implemented for all configured chromosomes |
| 4 | `compare_results.py` | Implemented with a one-to-one secure/R join |

The complete local terminal workflow is available. Scatter and Manhattan plot
generation is not implemented yet and remains the next analysis step.

## Step 0: Generate 1000 Genomes test data

Step 0 prepares source data for local correctness testing. It is not part of
the AoU Workbench workflow.

### Requirements

- Python 3
- `curl`
- `plink2` available on `PATH`
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
├── covariates.tsv
└── work/phase3.keep
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
- Covariates are AFR, AMR, EAS, and EUR super-population indicators. SAS is the
  reference level to avoid collinearity with the intercept.

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
├── prepared/
│   ├── chr21/
│   │   ├── A/
│   │   │   ├── geno/block.<gene-index>.bin
│   │   │   ├── cov.txt
│   │   │   └── pheno.txt
│   │   ├── B/
│   │   │   ├── geno/block.<gene-index>.bin
│   │   │   ├── private/block.<gene-index>.bin
│   │   │   ├── cov.txt
│   │   │   └── pheno.txt
│   │   ├── genes.txt
│   │   ├── block_sizes.txt
│   │   └── pos.txt
│   └── chr22/
```

Phenotypes and covariates are loaded once. The A/B sample split and row
order are also selected once and reused across every configured
chromosome.

Gene selection supports `random`, `file`, and `all` modes. In `random` mode,
`per_chromosome` is the number selected independently for each chromosome.

## Step 2: Run secure Burden and SKAT

Run three local parties from one parent command:

```bash
go run -mod=vendor secure-rvas.go run \
  --config run.1kg.conf
```

The parent creates temporary shared PRG keys, starts parties 0, 1, and 2, and
waits for all three processes. Each party opens one network and one MPC session,
`SetupNull` runs once for all configured chromosomes, and chromosomes run
sequentially in configuration order. No chromosome worker or independent MPC
lane is used.

The command writes:

```text
<run_dir>/secure/
├── chr21.tsv
├── chr22.tsv
├── all_secure_results.tsv
└── _SUCCESS
```

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
python3 rewrite/analysis/run_reference.py \
  --config run.1kg.conf
```

The wrapper runs every configured chromosome in order and writes:

```text
<run_dir>/reference/
├── chr21.csv
├── chr22.csv
├── all_r_results.csv
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
├── all_comparison.csv
└── _SUCCESS
```

The CSV retains all raw p-values and absolute errors. A missing Davies error is
preserved as an empty field rather than converted into an ordinary p-value.

## Complete local workflow

Use a fresh `run_dir` in `run.1kg.conf` for a new preprocessing run, then execute
the following commands from the repository root:

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
```

For a detached Workbench-terminal run, apply `nohup` to the long secure step:

```bash
nohup go run -mod=vendor secure-rvas.go run \
  --config run.1kg.conf > run-1kg.log 2>&1 &
```

Run the R reference and comparison after `secure/_SUCCESS` appears.

## Current end-to-end validation

The command sequence above was validated from Step 0 through Step 4 on the
current chromosome 21/22 fixture. It produced five genes per chromosome and two
phenotypes, for 20 rows in each combined result file. All three success markers
were created.

The current run reported:

- maximum Burden absolute difference: `3.13e-8`;
- maximum secure WH versus R-Liu absolute difference among nondegenerate rows:
  `0.0388124`;
- 14 degenerate gene-phenotype rows, all with R SKAT and secure SKAT `p=1`;
- no non-finite secure p-values.

Repeated Step 0 generation with the same inputs and seeds produced identical
gene-panel, annotation, phenotype, and covariate file hashes. Generated local
outputs live under `output/` and are excluded from Git.

Per-chromosome secure/plain scatter plots and combined-chromosome Manhattan
plots are the next analysis feature. AoU localization and normalization,
ancestry-specific execution, timing and memory measurement, resume support,
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
- The initial runner uses one MPC session and processes chromosomes
  sequentially. It does not use chromosome workers or independent MPC lanes.
- AoU-specific conversion, timing and memory instrumentation, ancestry splits,
  resume support, and parallel MPC lanes remain separate later work.

## License and attribution

The MHE and MPC backend in this repository is derived from the Lattigo v6-based
SF-GWAS implementation. See [`LICENSE`](LICENSE) for licensing information.
