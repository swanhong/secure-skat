#!/usr/bin/env python3
"""Federated-private SKAT data PREP: real AoU exome pgen -> per-gene int8 genotype blocks.

Simulates the MVP+AoU two-cohort split from a single AoU pgen: split samples into
cohort A (public-list owner) and B (private), split each gene's variants into shared/public_only/
private. Writes the int8 blocks the Go `skat_fed` mode reads:

    A/geno.<g>.bin    cohort A, public-list variants            [PART A]
    B/geno.<g>.bin    cohort B, public list ALIGNED (public_only cols = 0)   [PART A]
    B/priv.<g>.bin    cohort B, private variants                [PART B]
    {A,B}/cov.txt     real covariates: first 5 AoU ancestry PCs, in geno-row order
    {A,B}/pheno.txt   real phenotype (LDL, inverse-normal), in geno-row order
    manifest.json     per-gene m (public/private), gene list

Block format = row-major n*m int8 (dosage 0/1/2, missing<0 -> Go reads 0); identical to
scripts/plinkBedToBinary.py output and gwas.loadDenseBlocks. Block assembly/alignment + the
plaintext federated==pooled check live in skat_plain_local.py.

DATA ONLY: extracts/shapes genotypes + real PC covariates + real LDL phenotype. No SKAT
computation. Samples = geno ∩ pheno(LDL non-missing) ∩ ancestry(PC); join key person_id==research_id.

Key = chr:pos:ref:alt / GRCh38 / biallelic-only / PASS (see .local/warning.md). Secure SKAT is
n-independent, so N_SUB subsamples samples to keep blocks small; m stays realistic.

    python3 fed_prep.py    # real prep on the workbench (needs plink2 + AoU pgen)
"""
import json
import math
import os
import subprocess

import numpy as np

from skat_plain_local import build_party_blocks

# ---- config (env-overridable; defaults = AoU Workbench Controlled-Tier layout, which AoU
CHR = os.environ.get("FED_CHR", "chr22")
PGEN = os.path.expanduser(os.environ.get("FED_PGEN",
    f"~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/pgen/exome.{CHR}"))
GENCODE = os.path.expanduser(os.environ.get("FED_GENCODE",  # public GENCODE v44 pc-gene coords, committed next to this script
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "gencode_v44_pc_genes.bed")))
ANCESTRY = os.path.expanduser(os.environ.get("FED_ANCESTRY",
    "~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/echo_v4_r2.ancestry_preds.tsv"))
PHENO_CSV = os.path.expanduser(os.environ.get("FED_PHENO", "~/fed_prep_in/pheno.csv"))
PHENO_COL = os.environ.get("FED_PHENO_COL", "inv_LDLC_final_mgdl_6sd_masked")  # LDL, inverse-normal (continuous)
OUT_DIR = os.path.expanduser(os.environ.get("FED_OUT", "~/fed_prep_out"))
N_PCS = 5                    # first N PCs from ancestry_preds pca_features used as covariates (age/sex deferred)
N_SUB = 5000                 # samples per cohort (secure is n-independent; keeps blocks small)
N_GENES = 20
SEED = 71
FRAC_SHARED, FRAC_PUBONLY = 0.6, 0.2   # rest = private
PLINK2 = os.environ.get("PLINK2", "plink2")   # override: PLINK2=/path/to/plink2 python3 fed_prep.py
MAX_ALLELE_LEN = 1000   # --set-all-var-ids cap; long indels keep full chr:pos:ref:alt key (chr22 max=193). Bump if a chrom exceeds this.


def sh(cmd):
    print("  $", cmd)
    subprocess.run(cmd, shell=True, check=True)


def write_int8_block(path, mat):
    """row-major n*m int8 (matches plinkBedToBinary.py + gwas.loadDenseBlocks)."""
    np.asarray(np.rint(mat), dtype=np.int8).tofile(path)


def write_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, out_dir):
    A_blocks, B_aligned, B_priv, _ = build_party_blocks(
        gene_keys, priv_keys, roles, A_geno, B_geno, keycol)
    os.makedirs(f"{out_dir}/A", exist_ok=True)
    os.makedirs(f"{out_dir}/B", exist_ok=True)
    for g in range(len(gene_keys)):
        write_int8_block(f"{out_dir}/A/geno.{g}.bin", A_blocks[g])
        write_int8_block(f"{out_dir}/B/geno.{g}.bin", B_aligned[g])
        write_int8_block(f"{out_dir}/B/priv.{g}.bin", B_priv[g])
    json.dump({"n_genes": len(gene_keys),
               "pub_m": [len(k) for k in gene_keys],
               "priv_m": [b.shape[1] for b in B_priv]},
              open(f"{out_dir}/manifest.json", "w"), indent=2)
    print(f"  wrote blocks + manifest -> {out_dir}")


def load_ancestry_pcs(path, n_pcs):
    """research_id -> [n_pcs PCs] from AoU ancestry_preds.tsv (pca_features = '[v1, ..., v16]')."""
    pcs = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rid, pca = header.index("research_id"), header.index("pca_features")
        for ln in f:
            x = ln.rstrip("\n").split("\t")
            vals = [float(v) for v in x[pca].strip().strip("[]").split(",")]
            if len(vals) < n_pcs:
                raise ValueError(f"pca_features has {len(vals)} PCs < N_PCS={n_pcs}")
            pcs[x[rid]] = vals[:n_pcs]
    return pcs


def write_cov(fam_ids, pcs, out_path):
    """n x n_pcs covariates (ancestry PCs) in geno (.fam) row order. Eligible filter guarantees presence."""
    np.savetxt(out_path, np.asarray([pcs[sid] for sid in fam_ids]))


def load_pheno(path, col):
    """person_id -> phenotype value (csv col); skips blank/NA/non-finite. person_id == research_id."""
    ph = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split(",")
        pid, pc = header.index("person_id"), header.index(col)
        for ln in f:
            x = ln.rstrip("\n").split(",")
            v = x[pc].strip()
            if not v or v.upper() in ("NA", "NAN"):
                continue
            fv = float(v)
            if math.isfinite(fv):  # drop nan/inf so they can't poison the null model
                ph[x[pid]] = fv
    return ph


def write_pheno(fam_ids, pheno, out_path):
    """n-vector phenotype in geno (.fam) row order. Eligible filter guarantees presence."""
    np.savetxt(out_path, np.asarray([pheno[sid] for sid in fam_ids]))


def plink_extract_to_int8(pgen, keep_file, keys_file, n, out_prefix):
    """plink2 (set-var-ids + biallelic + keep + extract) -> .bed -> plinkBedToBinary -> int8;
    returns (n x m matrix, keys, fam_ids). Re-applies --set-all-var-ids on the ORIGINAL pgen so the
    extracted IDs match the keyed .pvar -- avoids rewriting the full keyed pgen. plink2 keeps the
    pgen's variant order, not --extract order, so the caller reorders by the returned keys."""
    here = os.path.dirname(os.path.abspath(__file__))
    sh(f"{PLINK2} --pfile {pgen} --max-alleles 2 --min-alleles 2 "
       f"--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len {MAX_ALLELE_LEN} "
       f"--keep {keep_file} --extract {keys_file} --indiv-sort none --make-bed --out {out_prefix}")
    bim_keys = [ln.split('\t')[1] for ln in open(f"{out_prefix}.bim")]
    fam_ids = [ln.split()[1] for ln in open(f"{out_prefix}.fam")]  # col2 = IID = research_id; geno-row order
    if len(fam_ids) != n:  # plink dropped samples -> g (reshaped n×m) would misalign cov/pheno
        raise SystemExit(f"plink emitted {len(fam_ids)} samples != requested {n} ({out_prefix})")
    m = len(bim_keys)
    sh(f"python3 {here}/../plinkBedToBinary.py {out_prefix}.bed {n} {m} {out_prefix}.bin")
    g = np.fromfile(f"{out_prefix}.bin", dtype=np.int8).reshape(n, m)
    g[g < 0] = 0  # missing -> 0 (matches Go loader)
    return g, bim_keys, fam_ids


def run():
    """Real prep on the workbench. Reads AoU pgen, splits into 2 cohorts, writes genotype blocks."""
    rng = np.random.default_rng(SEED)
    os.makedirs(OUT_DIR, exist_ok=True)

    keyed = f"{OUT_DIR}/{CHR}_keyed"
    sh(f"{PLINK2} --pfile {PGEN} --max-alleles 2 --min-alleles 2 "
       f"--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len {MAX_ALLELE_LEN} "
       f"--make-just-pvar --out {keyed}")

    # (2) gene -> PASS variants (read keyed .pvar; map to GENCODE genes)
    genes = load_gencode_genes(GENCODE, CHR, N_GENES)
    gene_keys_all = scan_pvar_into_genes(f"{keyed}.pvar", genes)  # gene -> ordered [keys] (PASS only)

    # (3) eligible = geno ∩ pheno(LDL) ∩ ancestry(PC); split into cohort A/B; per-gene role split
    pcs = load_ancestry_pcs(ANCESTRY, N_PCS)
    pheno = load_pheno(PHENO_CSV, PHENO_COL)
    psam = [ln.split()[0] for ln in open(f"{PGEN}.psam") if not ln.startswith("#")]  # samples unchanged by keying
    gset = set(psam)
    eligible = [s for s in psam if s in pheno and s in pcs]
    # per-component + pairwise counts so a person_id!=research_id namespace mismatch is visible
    print(f"  geno={len(gset)} pheno={len(pheno)} pc={len(pcs)} | "
          f"geno∩pheno={len(gset & pheno.keys())} geno∩pc={len(gset & pcs.keys())} eligible={len(eligible)}")
    if len(eligible) < 2 * N_SUB:
        raise SystemExit(f"only {len(eligible)} eligible samples < 2*N_SUB={2 * N_SUB} (check person_id==research_id matching)")
    perm = rng.permutation(len(eligible))
    A_ids = [eligible[i] for i in perm[:N_SUB]]
    B_ids = [eligible[i] for i in perm[N_SUB:2 * N_SUB]]
    write_lines(f"{OUT_DIR}/A.keep", A_ids, "#IID"); write_lines(f"{OUT_DIR}/B.keep", B_ids, "#IID")

    gene_keys, priv_keys, roles_all, all_keys = [], [], {}, []
    for keys in gene_keys_all:
        roles, pub, priv = split_roles(rng, keys)
        gene_keys.append(pub); priv_keys.append(priv)
        roles_all.update(roles); all_keys += keys
    A_extract = [k for k in all_keys if roles_all[k] in ("shared", "public_only")]
    B_extract = [k for k in all_keys if roles_all[k] in ("shared", "private")]
    write_lines(f"{OUT_DIR}/A_keys.txt", A_extract)
    write_lines(f"{OUT_DIR}/B_keys.txt", B_extract)

    # (4) extract genotypes (plink2 -> int8), reorder to a single key->column map
    Ag, Ak, Afam = plink_extract_to_int8(PGEN, f"{OUT_DIR}/A.keep", f"{OUT_DIR}/A_keys.txt",
                                         N_SUB, f"{OUT_DIR}/A_geno")
    Bg, Bk, Bfam = plink_extract_to_int8(PGEN, f"{OUT_DIR}/B.keep", f"{OUT_DIR}/B_keys.txt",
                                         N_SUB, f"{OUT_DIR}/B_geno")
    A_geno, B_geno, keycol = merge_cohort_columns(Ag, Ak, Bg, Bk)
    del Ag, Bg  # free pre-merge matrices (RAM-tight workbench)
    # drop any keys plink2 didn't emit (monomorphic/filtered) so build never KeyErrors
    present = set(keycol)
    gene_keys = [[k for k in g if k in present] for g in gene_keys]
    priv_keys = [[k for k in p if k in present] for p in priv_keys]

    # (5) write genotype blocks + real covariates (5 PCs) + real LDL phenotype, all geno-row order
    print("run (real AoU pgen):")
    write_blocks(gene_keys, priv_keys, roles_all, A_geno, B_geno, keycol, OUT_DIR)
    write_cov(Afam, pcs, f"{OUT_DIR}/A/cov.txt")
    write_cov(Bfam, pcs, f"{OUT_DIR}/B/cov.txt")
    write_pheno(Afam, pheno, f"{OUT_DIR}/A/pheno.txt")
    write_pheno(Bfam, pheno, f"{OUT_DIR}/B/pheno.txt")
    print(f"  wrote cov ({N_PCS} PCs) + pheno (LDL) -> {OUT_DIR}/{{A,B}}/")


# ---- helpers ----

def write_lines(path, lines, header=None):
    with open(path, "w") as f:
        if header:
            f.write(header + "\n")
        f.write("\n".join(lines) + "\n")


def split_roles(rng, keys):
    """label each key shared/public_only/private; returns roles dict + (pub_keys, priv_keys)."""
    idx = rng.permutation(len(keys))
    n_sh = int(round(FRAC_SHARED * len(keys)))
    n_pub = int(round(FRAC_PUBONLY * len(keys)))
    roles = {}
    for j, k in enumerate(idx):
        roles[keys[k]] = "shared" if j < n_sh else ("public_only" if j < n_sh + n_pub else "private")
    pub = [k for k in keys if roles[k] in ("shared", "public_only")]
    priv = [k for k in keys if roles[k] == "private"]
    return roles, pub, priv


def load_gencode_genes(bed, chrom, n_genes):
    """gene -> (lo, hi) for n_genes protein-coding genes on chrom, spread across it."""
    genes = []
    for ln in open(bed):
        f = ln.rstrip("\n").split("\t")
        if f[0] != chrom:
            continue
        name = f[3] if len(f) > 3 else f"{f[0]}:{f[1]}"
        genes.append((name, int(f[1]), int(f[2])))
    genes.sort(key=lambda x: x[1])
    step = max(1, len(genes) // n_genes)
    picked = genes[::step][:n_genes]
    return {name: (lo, hi) for name, lo, hi in picked}


def scan_pvar_into_genes(pvar, genes):
    """assign PASS biallelic keys to genes by position; return list of ordered key-lists."""
    buckets = {name: [] for name in genes}
    for ln in open(pvar):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        chrom, pos, vid, ref, alt, flt = f[0], int(f[1]), f[2], f[3], f[4], f[5]
        if flt not in ("PASS", "."):
            continue
        if "," in alt:        # multiallelic already excluded by --max-alleles, belt-and-suspenders
            continue
        for name, (lo, hi) in genes.items():
            if lo <= pos < hi:
                buckets[name].append(vid)  # vid is the chr:pos:ref:alt key set by --set-all-var-ids
                break
    return [v for v in buckets.values() if v]


def merge_cohort_columns(Ag, Ak, Bg, Bk):
    """put A and B genotype matrices into one shared column space keyed by variant key."""
    assert len(set(Ak)) == len(Ak) and len(set(Bk)) == len(Bk), "duplicate variant key in a cohort .bim"
    keys = list(dict.fromkeys(Ak + Bk))
    keycol = {k: i for i, k in enumerate(keys)}
    A = np.zeros((Ag.shape[0], len(keys)), dtype=Ag.dtype)
    B = np.zeros((Bg.shape[0], len(keys)), dtype=Bg.dtype)
    for j, k in enumerate(Ak):
        A[:, keycol[k]] = Ag[:, j]
    for j, k in enumerate(Bk):
        B[:, keycol[k]] = Bg[:, j]
    return A, B, keycol


if __name__ == "__main__":
    run()
