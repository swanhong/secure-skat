# Input summary statistics

After Step 0 has created the local PGEN triplets, normalized annotation, and gene
panels, run from the repository root in AoU Workbench:

```bash
python3 -m python.analysis.summarize_inputs --config config/aou
```

This writes `variant_summary.csv`, `gene_summary.csv`, and `dataset_summary.csv`
here. Use `--output-dir PATH` to choose another directory. The utility uses the
existing NumPy and pgenlib dependencies and does not download data.

It includes all samples in the PGEN files and all genes passing the configured
annotation masks and annotation MAF cutoff, plus PASS/`.` and biallelic variant
filters. It does not use ancestry selection, phenotype/covariate availability,
A/B sample selection, random/file gene selection, or shared/private roles.
The configured autosomes must contain the same ordered sample IDs.

- Variant columns: `variant_id,gene,N,MAC,MAF`. `N` counts nonmissing calls;
  `MAC = min(ALT_count, 2*N - ALT_count)` and `MAF = MAC/(2*N)`. All-missing
  variants have `N=0,MAC=0` and an empty MAF. Missing genotypes are not imputed.
- Gene columns: `gene,N_variants,N_samples,total_MAC`. `N_samples` counts samples
  with at least one observed call among the retained gene variants, including
  `0/0`. Genes with no retained variants are omitted; MAC=0 variants are retained.
- Dataset columns: `N_samples_total,N_variants_total,N_genes_total,MAC_total`.
  Samples are counted once across chromosomes. Variant count and MAC are counted
  once per variant ID, even when a variant maps to multiple genes.

These are input-dataset summaries, not phenotype-specific test sample sizes.
Reported MAF is calculated from this cohort; the annotation MAF used for filtering
may come from an external population. Only diploid autosomal hard calls are
supported. CSVs are replaced only after the full calculation succeeds and are
ignored by Git.
