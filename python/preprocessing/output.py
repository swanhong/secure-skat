import csv
from pathlib import Path
import numpy as np
from .model import GeneBlock, GeneRef, PhenoCovRows

def write_selected_genes(
        path: Path,
        genes: tuple[GeneRef, ...],
) -> None:
    with path.open("w", newline="") as file:
        writer = csv.writer(
            file,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(("gene_id", "gene_symbol", "chromosome", "order_index"))
        for gene in genes:
            writer.writerow((gene.gene_id, gene.gene_symbol, gene.chromosome, str(gene.order_index)))

def write_outputs(
    out_dir: Path,
    blocks: tuple[GeneBlock, ...],
    rows_a: PhenoCovRows,
    rows_b: PhenoCovRows,
    party: int = 0,
) -> None:
    if out_dir.exists() and any(out_dir.iterdir()):
        raise ValueError(f"out_dir is not empty: {out_dir}")

    a_geno_dir = out_dir / "A" / "geno"
    b_geno_dir = out_dir / "B" / "geno"
    b_private_dir = out_dir / "B" / "private"

    if party != 2:
        a_geno_dir.mkdir(parents=True, exist_ok=True)
    if party != 1:
        b_geno_dir.mkdir(parents=True, exist_ok=True)
        b_private_dir.mkdir(parents=True, exist_ok=True)
    for index, block in enumerate(blocks):
        if party != 2:
            np.asarray(block.public_a, dtype=np.int8).tofile(a_geno_dir / f"block.{index}.bin")
        if party != 1:
            np.asarray(block.public_b, dtype=np.int8).tofile(b_geno_dir / f"block.{index}.bin")
            np.asarray(block.private_b, dtype=np.int8).tofile(b_private_dir / f"block.{index}.bin")
    for name, rows in (("A", rows_a), ("B", rows_b)):
        if (party == 1 and name == "B") or (party == 2 and name == "A"):
            continue
        np.savetxt(out_dir / name / "cov.txt", rows.covariates, fmt="%.17g", delimiter="\t")
        np.savetxt(out_dir / name / "pheno.txt", rows.phenotypes, fmt="%.17g", delimiter="\t")
    (out_dir / "genes.txt").write_text("".join(f"{block.gene.gene_id}\n" for block in blocks))
    (out_dir / "block_sizes.txt").write_text(
        "".join(f"{len(block.public_variants)}\n" for block in blocks)
    )
    (out_dir / "pos.txt").write_text(
        "".join(
            f"{block.gene.chromosome}\t{variant.position}\n"
            for block in blocks
            for variant in block.public_variants
        )
    )
