# AoU Workbench runner

`submit_main.sh` submits one coordinator job and returns without waiting. The
coordinator runs the requested chromosomes concurrently on one Batch VM, then
aggregates and plots their results. It does not call the legacy submit path.

```text
submit_main.sh
  -> coordinator
     -> concurrent run_chromosome_task.sh workers
       -> preprocessing
       -> three localhost secure-skat processes
       -> pooled R Burden/Davies
       -> gene_results.csv
     -> aggregation
       -> all_gene_results.csv + summaries + plots
```

## Run configuration

Edit the root `run.conf`, then submit it:

```bash
bash submit_main.sh run.conf
```

`run.conf` separates the two required settings from optional settings:

```bash
annotation_dir=/path/to/chromosome_annotations
data_bits=DATA_BITS
```

Replace `DATA_BITS` with the value validated for the intended cohort size.
Delete or comment out an optional setting to use the default declared in
`submit_main.sh`. Do not add `export`; `run.conf` is loaded only inside the
submit process and does not modify the caller's shell. It uses shell assignment
syntax, so use a trusted local file and do not put spaces around `=`.

`project` defaults first to the Workbench `GOOGLE_CLOUD_PROJECT` value and then
to the active `gcloud` project. `service_account` defaults to the Workbench
`PET_SA_EMAIL` value, the user-specific Terra service account used by the Batch
VM to access controlled resources. Set either value in `run.conf` only when
automatic detection is unavailable or must be overridden.

The phenotype input is comma-separated and defaults to the five configured
lipid columns keyed by `person_id`. The covariate input is tab-separated and
defaults to the AoU ancestry table keyed by `research_id`; its `pca_features`
array is expanded to `PC1` through `PC16`. To use an already-flat covariate
table instead, set `covariate_array_column=` and provide its numeric columns in
`covariate_columns`. Other defaults include `PN14QP436S45`, 30 fractional bits,
50 probes, and seed 42.

Annotation `variant_key` values may be exact PVAR IDs or `POS:REF:ALT`
coordinates. Preprocessing maps coordinate keys to the corresponding PVAR IDs
before gene and variant selection.

`data_bits` is intentionally required. The 60-bit value has been validated
only by the current small tests and is not an AoU full-cohort production
contract.

`chromosomes` defaults to all autosomes. Set a comma-separated subset such as
`chromosomes=21,22` to run and aggregate only those chromosomes.
`max_parallel_chromosomes` defaults to `2` and limits how many chromosome
workers share the coordinator VM at once.

`submit_main.sh` submits one coordinator dsub job without `--wait` and returns
after Batch accepts it. The coordinator localizes shared inputs once, runs the
selected chromosomes concurrently in isolated work directories and disjoint
MPC port ranges, up to `max_parallel_chromosomes` at a time, then aggregates
their results and uploads one result tree. All workers share the coordinator
VM's configured cores and memory. The configured Batch timeout applies to the
complete job. The exact submitted `run.conf` is copied to the run's GCS root.

## Per-chromosome flow

1. Localize one chromosome's PGEN/PVAR/PSAM and annotation.
2. Run preprocessing and write `metadata.json`, `gene_metadata.tsv`, and
   `phenotypes.txt`.
3. Stop with `validation_errors.csv` if any public gene width exceeds 4096.
4. Generate fresh task-local MPC shared keys and run three localhost parties.
5. Run pooled R::SKAT Burden and Davies on the same preprocessing blocks.
6. Join the results and write the chromosome error summary.

Only party A writes `secure_results.csv`. Preprocessing blocks stay on the task
VM and are not copied to final GCS results.

## Outputs

Each successful chromosome is written below `results/`:

```text
chromosomes/chrN/
  secure_results.csv
  r_results.csv
  gene_results.csv
  error_summary.csv
  metadata.json
  party0.log
  party1.log
  party2.log
  run.log
  _SUCCESS
```

The coordinator writes the aggregate below the same result root:

```text
final/
  all_gene_results.csv
  error_summary.csv
  plots/scatter_*.png
  plots/manhattan_*.png
  _SUCCESS
```

Aggregation starts only after every selected chromosome succeeds. A failed
chromosome or public-width validation triggers a best-effort upload of partial
results and logs; the coordinator task's exit status is the authoritative
signal. Temporary genotype and protocol work directories are removed after
each successful chromosome so disk use does not accumulate across the run.

## Reading order

Read the implementation in execution order, not directory order:

1. `submit_main.sh`: package inputs and submit one coordinator job
2. `run_coordinator_task.sh`: shared setup, chromosome loop, and aggregation
3. `run_chromosome_task.sh`: one chromosome's complete flow
4. `rewrite/preprocessing/cli.py`, `prepare()`: original inputs to protocol blocks
5. `rewrite/cmd/secure-skat/run.go`, `runSecure()`: setup once and natural batches
6. `runGeneBatch()` in the same Go file: weights, statistics, finalization, release
7. `run_reference_skat.R`: pooled R Burden and Davies reference
8. `results.py`, `join_results()`, `aggregate_results()`, and plots

`input.go`, key generation, CSV writers, and plotting helpers are support code;
they are not part of the secure statistical algorithm itself.
