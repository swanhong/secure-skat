# Packed score precision diagnostic

This command tests whether the public SKAT score

\[
G^T y_0-(G^T X)\hat\beta
\]

is numerically safe to move from the current SS implementation to across-gene
packed CKKS. It reads an existing `fed_prep.py` output directory and runs in one
local process. It does not start MPC parties or access the original PGEN.

## Workbench commands

Fast plaintext cancellation/conditioning screen for the proposed public-score
migration (run this first):

```bash
FED_OUT=~/fed_prep_out
go run -mod=vendor ./cmd/skat-score-precision \
  --out "$FED_OUT" \
  --public-only \
  --plain-only \
  --top 0 \
  | tee "$FED_OUT/score_precision_plain_public.log"
```

Actual packed PN14 CKKS differential, still without MPC/networking:

```bash
go run -mod=vendor ./cmd/skat-score-precision \
  --out "$FED_OUT" \
  --public-only \
  --ckks PN14QP438 \
  --threads 1 \
  --top 0 \
  | tee "$FED_OUT/score_precision_pn14_public.log"
```

### N14-only precision sweep

The diagnostic also provides two experimental profiles that keep `LogN=14`
and 8192 slots while raising the scale. Their rescaling primes are matched to
the scale; simply increasing `LogDefaultScale` on the original 34-bit chain
would cause scale growth after multiplication and is not an equivalent test.

| Profile | Default scale | Maximum level | Nominal logQP | Purpose |
| --- | ---: | ---: | ---: | --- |
| `PN14QP438` | 34 bits | 9 | 438 | existing baseline |
| `PN14QP431S40` | 40 bits | 7 | 431 | N14 with PN15-like scale |
| `PN14QP436S45` | 45 bits | 6 | 436 | N14 with PN16-like scale |

Run both experimental profiles on the same full gene set:

```bash
go run -mod=vendor ./cmd/skat-score-precision \
  --out "$FED_OUT" \
  --public-only \
  --ckks PN14QP431S40 \
  --threads 1 \
  --top 0 \
  --progress-every 50 \
  --csv "$FED_OUT/score_precision_pn14_s40_public_run1.csv" \
  2>&1 | tee "$FED_OUT/score_precision_pn14_s40_public_run1.log"

go run -mod=vendor ./cmd/skat-score-precision \
  --out "$FED_OUT" \
  --public-only \
  --ckks PN14QP436S45 \
  --threads 1 \
  --top 0 \
  --progress-every 50 \
  --csv "$FED_OUT/score_precision_pn14_s45_public_run1.csv" \
  2>&1 | tee "$FED_OUT/score_precision_pn14_s45_public_run1.log"
```

These are custom, diagnostic-only score-circuit experiments, not new production
parameters; do not copy their names into `configGlobal.toml` or a federated
GWAS run. At depth 3, each profile passes Lattigo's two-party `lambda=128`
statistical-masking bound check for collective refresh. That conservative
budget check is covered by a unit test, but the command above does not execute
a collective refresh or establish end-to-end protocol security. Weight
computation should therefore be evaluated separately as two waves (`base -> w2 -> w4 -> w8`,
refresh, then `w8 -> w16 -> w24`, refresh) before either profile is wired into
the federated path.

For a low-cost CKKS smoke test, add `--max-genes 50`; a limited sample can
reject but cannot accept the design. For a full score-circuit check, run all
genes three times (separate log/CSV names) because CKKS encryption noise is
random. Use the N14-only sweep above before changing the ring dimension; if it
is still borderline, repeat with `--ckks PN15QP880`. Omit
`--public-only` only if the existing B-private CKKS score is also in scope.

The plaintext screen reports SKAT-Q and Burden-L sensitivity, a cancellation
condition number, and an SS-like fixed-point baseline. It can reject an
ill-conditioned CKKS formulation but cannot accept it. The optional CKKS run is
also only a direct-`Enc(beta)` score-circuit differential: it does not reproduce
the secure fixed-point Cholesky result or masked SS-to-CKKS conversion used to
create production `betaRep`. A passing production gate therefore still needs a
controlled run through the real secure beta path.

Both modes use the live code's party-local allele orientation, party-level
subtract-then-aggregate order, non-intercept ridge placement, configured beta
fractional bits, covariate count, and SKAT weights. The stable reference computes
`G^T(y-X beta)` rather than subtracting two large sufficient-statistic terms.
Zero-Q and non-finite genes are counted explicitly and are never silently
removed from a pass denominator.

With `--top 0`, stdout contains cohort sizes plus aggregate error
quantiles/counts and omits the per-gene worst tables. The CSV and any log made
without `--top 0` remain Workbench research outputs; do not export them without
the applicable review. Even the aggregate-only log should follow the project's
normal output-review rules rather than being assumed automatically export-safe.
