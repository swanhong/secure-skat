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

For a low-cost CKKS smoke test, add `--max-genes 50`; a limited sample can
reject but cannot accept the design. For a full score-circuit check, run all
genes three times (separate log/CSV names) because CKKS encryption noise is
random. If PN14 is borderline, repeat with `--ckks PN15QP880`. Omit
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
