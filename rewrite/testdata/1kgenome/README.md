# Local 1000 Genomes workflow

This workflow runs the secure Burden/SKAT pipeline locally with public 1000 Genomes data and synthetic annotations / phenotypes. The generated fixture is for correctness testing, not biological interpretation.

Run commands from the repository root.

## Download and generate the test data

```bash
source setup/env.sh
python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --config run.1kg.conf \
  --num-pheno 2
```

The generator reads the chromosomes from `run.1kg.conf`. It downloads missing
1000 Genomes VCFs and sample metadata plus GENCODE v50, then writes the local
fixture under:

```text
rewrite/testdata/1kgenome/generated/
```

Downloaded files are reused on later runs.

## Run the complete workflow

```bash
./run_1kg_workflow.sh
```

The workflow also runs the download/generation step automatically, followed by secure preprocessing, secure Burden/SKAT, the reference calculation, comparison, plots, and the metrics summary.

Setup, configuration, individual stages, and output details are the same as the AoU workflow. See the repository [`README.md`](../../../README.md).
