# joint-stat

Independent three-party descriptive statistics for figures A/B/C, using the
prepared inputs of the main analysis. It does not run SKAT/Burden or their timers.
Docker/sf-kit integration is not included yet.

## Run

Run the same public options at both sites and the helper:

```sh
go run -mod=vendor secure-rvas.go joint-stat --config config/mvp --party 1
go run -mod=vendor secure-rvas.go joint-stat --config config/aou --party 2
go run -mod=vendor secure-rvas.go joint-stat --config config/cp0 --party 0
```

This single command reads local inputs, computes secure aggregates, releases them
only to parties 1 and 2, and saves both sites' results. No separate statistics
preparation command or count cache is needed. Other parties must be started
separately. Main and joint-stat must not use the same ports concurrently.

Public options (identical at all parties):

- `--carrier-thresholds 20,30,50,100,200,500,1000`
- `--cmac-bins 20,40,80,160,320,640,1280,2560,5120`
- `--batch-cells 1048576`: target input field elements per carrier batch; at most
  1,024 sample rows and at least one row across all genes. This is not a RAM cap.

Use ordinary site preparation with cache format 6 first. Older zero-filled
blocks lack callability and cannot recover it: prepare into a new `run_dir`.
Joint inputs have `A/joint_stat/` or `B/joint_stat/` bit-packed callability files.
Private keys are in `B/private_variants.txt`, one space-separated line per gene
(including empty lines), matching `private_allele_counts.txt`. Other dimensions
and identities come from `pheno.txt`, `genes.txt`, `block_sizes.txt`, and
`public_variants.tsv`. Local and legacy party 0 split outputs omit callability.
Keep prepared files unchanged throughout the run; files are not fingerprinted.

## Definitions

- The gene universe, variants and complete-case analysis samples are those of
  the joint prepared inputs. Local association masks are not substituted.
- Pooled MAF filtering uses ALL-cohort AC and fixed `2 * cohort_N`, matching main.
  At MAF 0.5 no filter is applied. Filter bits remain secret.
- Direction matches main: count REF when pooled analysis ALT AC is **> N**,
  otherwise ALT (including ties, as verified against main's MPC comparator). Missing/unavailable cells contribute neither carrier nor AC;
  the fixed-N orientation rule is not a called-AN frequency estimate.
- **A:** genes meeting each carrier threshold in MVP, AoU, either site
  (`local_any`), or joint. Each sample counts once per gene. Sites must contain
  disjoint participants.
- **B:** MVP/AoU/joint cMAC histograms on the identical genes with >=20 carriers
  in **both** sites. Bins are `[lower, upper)` with an unbounded final bin.
- **C:** distinct observed variants after filtering: shared, MVP-only, AoU-only.
  Multiple gene annotations count once per chromosome. `chr21` and `21` prefixes
  are equated; REF/ALT and indel representations must already be harmonized.
  "Only" refers to prepared inputs, not proof of biological absence elsewhere.
- Annotation is fixed to `HC_pLoF`; joint-stat does not recover variants removed
  upstream or reannotate them. Do not sum ALL with its ancestry subsets or
  combine different MAF thresholds.

## Outputs and timing

Each data site writes `<run_dir>/resource_stats/`:

```text
resource_stats.tsv       # identical aggregate rows at parties 1 and 2
metadata.json            # version, annotation, selected axes, thresholds and batches
metrics_partyN.tsv       # stage seconds, bytes and message counts (helper also writes)
_SUCCESS_partyN          # written after local outputs and metrics succeed
```

TSV columns:

```text
ancestry maf_threshold annotation chromosome gene_set panel series lower upper count status
```

Metrics are separate from main: `local_prepare`, `setup`, `secure_aggregate`,
`release`, `write`, and inclusive `figure_total`. MAFs are computed together, so
stage rows use `maf_threshold=0` to mean all configured MAFs. Secure timing includes
block reads, local contributions and peer waiting. Message counts are application
messages, not network latency rounds. Do not add `figure_total` to its components
or add times across parties. Whole-workflow wall time includes both protocols.

Carrier batches share indicators once and use `SSMultMats` across genes and MAFs,
then one vector comparison for all carrier flags in that batch. Files are read
by sample range. Network and MPC state are reused across chromosomes within each
ancestry; public layout is checked for every chromosome. Fresh domain-separated
PRG streams are derived once per ancestry session and never reset between its
chromosomes. Full AoU/MVP runtime and peak memory still require a server pilot.

## Release boundary and validation

Raw genotype/callability and private variant identities remain local. Gene and
variant counts remain secret until final aggregate release with `RevealSymVec`;
the helper receives no aggregate plaintext. Existing MPC trust assumptions apply.
`released` is not institutional export approval; no automatic cell suppression
is applied. No gene/variant-resolved counts appear in the output table.

```sh
go test -mod=vendor ./rvas/jointstat ./rvas/protocol ./rvas/workflow -count=1 -timeout 120s
python3 -m unittest python.preprocessing.test_joint_stat_metadata python.preprocessing.test_local_prepare python.preprocessing.test_symbolic_alt
```

Tests compare complete TSVs to an independent participant-level oracle, covering
MAF/threshold/bin boundaries, missing/unavailable calls, ties, repeated variants,
private-only/empty genes, batch boundaries and multi-chromosome/ancestry reuse.
