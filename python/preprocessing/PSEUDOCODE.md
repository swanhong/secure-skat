# Secure RVAS preprocessing

> 상태: local direct workflow 구현 완료. 이 문서는 현재 `secure-rvas prepare`가 사용하는
> Python 경계와 secure-input contract를 설명한다. 통계·representation 계약이 충돌하면
> `.local/summary/ch/02-input-data.md`와 `.local/summary/fed-prep-pseudocode.md`가 우선한다.

## Command boundary

사용자는 저장소 루트에서 다음 명령을 실행한다.

```bash
go run -mod=vendor secure-rvas.go prepare --config run.1kg.conf
```

호출 순서는 다음과 같다.

```text
secure-rvas.go
  → rewrite/workflow.Command
  → workflow.LoadConfig + ValidateConfig
  → workflow.Prepare
  → python3 -m rewrite.preprocessing.prepare
  → prepare_chromosomes
  → prepare_blocks (chromosome별)
```

Go는 config를 JSON request로 encoding하여 Python stdin으로 한 번 전달한다. 이 request를
별도 metadata 파일로 저장하지 않는다. Python module은 end-user config parser가 아니며 TOML
해석과 validation은 `rewrite/workflow/`가 소유한다.

## File ownership

```text
rewrite/preprocessing/
├── PSEUDOCODE.md  # 이 문서
├── model.py       # immutable preprocessing records
├── input.py       # PVAR/PSAM/gene/annotation/sample-table readers
├── selection.py   # random/file gene selection
├── prepare.py     # multi-chromosome request orchestration
├── pipeline.py    # row/variant/role/block transformations
├── plink.py       # PLINK2 extraction and dosage decoding
└── output.py      # secure-input and selected-gene writers
```

`rewrite/workflow/config.go`가 user-facing schema를 소유하고, `prepare.go`가 Python request를
조립한다. `rewrite/testdata/1kgenome/`은 이 library의 입력 형식과 맞는 local source data를
생성할 뿐 preprocessing 계산을 소유하지 않는다.

## Source input contract

```text
PLINK2 prefix:
  .pgen, .pvar, .psam

gene_panel.tsv:
  gene_id, gene_symbol, chromosome, order_index

annotation.tsv:
  variant_key, gene_id, gene_symbol, <annotation columns...>

phenotype.csv:
  configured ID column, configured phenotype columns

covariate.tsv:
  configured ID column, configured covariate columns
```

- `gene_id` is the join identity; `gene_symbol` is display metadata.
- `order_index` is zero-based genomic order.
- PVAR IDs and annotation `variant_key` use the canonical `chrom:pos:ref:alt` form.
- PVAR variants with FILTER other than `PASS` or `.` are excluded.
- If PVAR has no FILTER column, the reader treats every row as `PASS`.
- Phenotype/covariate values `""`, `NA`, `NaN`, `nan` and `.` become NaN and are removed by finite-row
  selection.
- Covariate files do not contain an intercept; the protocol adds it.

## Records

```text
GeneRef:
  gene_id, gene_symbol, chromosome, order_index

VariantRef:
  key, position, filter_value

AnnotationRef:
  variant_key, gene_id, gene_symbol, values

GeneVariants:
  gene, variants[]

GenePlan:
  gene, (variant, role)[]

PhenoCovRows:
  sample_ids[], covariates[n][c-1], phenotypes[n][q]

GeneBlock:
  gene
  public_variants[], private_variants[]
  public_a[nA][mPublic]
  public_b[nB][mPublic]
  private_b[nB][mPrivate]
```

Variant role is one of `shared`, `public_only` or `private`, and belongs to the `(gene, variant)`
occurrence rather than to a global variant key.

## Multi-chromosome orchestration

```text
function prepare_chromosomes(request):
  require at least one chromosome
  read phenotype and covariate tables once
  if file mode, read selected gene TSV once

  for chromosome in config order:
    read PVAR, PSAM, gene panel and annotation

    if first chromosome:
      remember ordered PSAM IDs
      select eligible A/B rows once
    else:
      require identical ordered PSAM IDs

    select gene/variant groups
    prepare_blocks using the same A/B rows

  if mode != all:
    write one selected_genes.tsv
```

Chromosome inputs are not given separate phenotype/covariate selections. They reuse the same in-memory
`rows_a` and `rows_b`, which is why later PSAM order must match the first chromosome.

## Row selection

```text
eligible := PSAM rows that
  exist in phenotype and covariate tables
  and have finite values in every requested column

shuffle eligible once with random.Random(sample_seed)

if samples_per_cohort > 0:
  require at least 2 * samples_per_cohort eligible rows
  A := first n rows
  B := next n rows
else:
  require at least two eligible rows
  A := ceil(N / 2) rows
  B := floor(N / 2) rows
```

Phenotype and covariate column order is exactly the config order. Sample IDs are alignment metadata and
are not written to protocol row files.

## Gene and variant selection

Gene selection modes are:

- `random`: form all non-empty filtered gene groups, sample `per_chromosome` groups with the configured
  seed, then restore genomic order.
- `file`: read requested gene IDs from a gene-panel-format TSV and keep chromosome-local genomic order.
- `all`: keep every chromosome-local gene from the panel, including groups that become empty.

For each annotation occurrence:

```text
keep only if
  its gene is selected
  and its PVAR variant has FILTER PASS or .
  and every equality mask matches
  and max_maf is absent or float(MAF) <= max_maf
```

Different mask columns are AND conditions. Config currently rejects repeated mask columns, so one
column has one accepted exact value. `max_maf` is separate from equality masks; if set, the annotation
must contain `MAF`.

Current implementation note: `select_gene_groups` forwards `max_maf` in `random` mode, which is the
current `run.1kg.conf` path, but omits it when making the final `select_gene_variants` call in `file`
and `all` modes. Those two modes must pass the same argument before they are used with an MAF bound.

A duplicate `(gene_id, variant_key)` annotation occurrence is rejected rather than silently creating
duplicate genotype columns. There is no preprocessing-time 4096-variant cutoff; protocol parameter
capacity remains a later check.

## Role assignment and PLINK extraction

For every selected gene, `assign_roles` shuffles its variants with `random.Random(role_seed)` and assigns
approximately 60% shared, 20% public-only and 20% private roles. The returned order remains PVAR order.

PLINK2 is called separately for A and B:

- A requests shared and public-only variants.
- B requests shared and private variants.
- Rows are reordered to the exact requested sample-ID order.
- ALT hard calls `0`, `1`, `2` become signed `int8`; missing `NA` becomes `0`.
- Fractional dosage is rejected.
- A public-only variant is represented by the real A column and an all-zero B public column.
- B-private columns are written only to `B/private/`.

## Output contract

```text
<run_dir>/
├── selected_genes.tsv                  # random/file modes only
└── prepared/chr<chromosome>/
    ├── A/
    │   ├── geno/block.<g>.bin
    │   ├── cov.txt
    │   └── pheno.txt
    ├── B/
    │   ├── geno/block.<g>.bin
    │   ├── private/block.<g>.bin
    │   ├── cov.txt
    │   └── pheno.txt
    ├── genes.txt
    ├── block_sizes.txt
    └── pos.txt
```

- `block.<g>.bin` is headerless row-major signed `int8`.
- `genes.txt`, `block_sizes.txt` and every block type share the same chromosome-local gene index.
- `block_sizes.txt` is the public variant count.
- `pos.txt` contains public chromosome/position rows in block order.
- A/B row files contain no sample IDs.
- Existing non-empty chromosome output is rejected.
- Partial resume and atomic whole-run output are not implemented.

## Required invariants

1. Each party's genotype, covariate and phenotype row order agrees.
2. Every chromosome reuses the same A/B sample split and X/Y values.
3. A/B public columns agree in order; B public-only columns are zero.
4. Public and private occurrences do not overlap within a gene.
5. Multi-gene variants receive roles independently by `(gene, key)` occurrence.
6. Gene and phenotype order is deterministic for a fixed config and source input.
7. Private variant identity/count is absent from public metadata.

## Deferred preprocessing work

- AoU localization and source normalization
- ancestry-specific sample filtering
- production-scale gene-group extraction instead of one dense chromosome extraction
- duplicate sample-ID policy for production inputs
- fractional dosage policy
- partial resume policy
