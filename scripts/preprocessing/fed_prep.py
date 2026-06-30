#!/usr/bin/env python3
"""Federated-private SKAT PREP (#4: synthesize two cohorts from a single AoU pgen).

Turns a per-chrom AoU exome pgen into the per-gene int8 genotype blocks the Go `skat_fed`
mode reads, by *simulating* the MVP+AoU scenario from AoU alone: split samples into cohort A
(public-list owner) and B (private), and split each gene's variants into shared / public_only /
private. Outputs (under OUT_DIR):

    A/geno.<g>.bin    cohort A, public-list variants            [PART A]
    B/geno.<g>.bin    cohort B, public list ALIGNED (public_only cols = 0)   [PART A]
    B/priv.<g>.bin    cohort B, private variants                [PART B]
    {A,B}/pheno.txt, {A,B}/cov.txt
    manifest.json     gene list, per-gene m, keys, roles

Block format = row-major n*m int8 (dosage 0/1/2, missing<0 -> Go reads as 0); identical to
scripts/plinkBedToBinary.py output and what gwas.loadDenseBlocks expects.

Genotype extraction = plink2 + plinkBedToBinary (existing, --run only). The NEW logic (gene map,
2-cohort split, B alignment/0-fill, federated-Q check) is pure python and is exercised by
`--selftest` WITHOUT plink2 or AoU data. Run --selftest locally; --run on the workbench.

    python3 fed_prep.py --selftest          # validate logic locally (no plink2/AoU)
    python3 fed_prep.py --run               # real prep on the workbench

Key = chr:pos:ref:alt / GRCh38 / biallelic-only / PASS (see warning.md). Secure SKAT is
n-independent, so N_SUB subsamples samples to keep the dense blocks small; m stays realistic.
"""
import argparse
import json
import os
import subprocess
import numpy as np

# ---- config (edit for the workbench) ----
CHR = "chr22"
PGEN = os.path.expanduser(f"~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/pgen/exome.{CHR}")
GENCODE = os.path.expanduser("~/Projects/mvp-secure/gencode/gencode_v44_pc_genes.bed")
OUT_DIR = os.path.expanduser("~/fed_prep_out")
N_SUB = 5000                 # samples per cohort (secure is n-independent; keeps blocks small)
N_GENES = 20
N_COV = 2
SEED = 71
FRAC_SHARED, FRAC_PUBONLY = 0.6, 0.2   # rest = private
PLINK2 = "plink2"


# ===================== shared plaintext logic (selftest + verify) =====================

def beta_density_weight(count, two_n):
    p = count / two_n
    maf = np.minimum(p, 1 - p)
    return 25.0 * (1 - maf) ** 24


def federated_Q_from_blocks(A_blocks, B_aligned, B_priv, XA, yA, XB, yB):
    """Per-gene Q from the block tensors, mirroring the secure path exactly.
    A_blocks[g]: nA x m_pub ; B_aligned[g]: nB x m_pub (public_only cols already 0) ;
    B_priv[g]: nB x m_privg . Returns list of per-gene Q."""
    nA, nB = len(yA), len(yB)
    two_n = 2.0 * (nA + nB)
    # pooled null model (covariates only)
    XtX = XA.T @ XA + XB.T @ XB
    Xty = XA.T @ yA + XB.T @ yB
    beta = np.linalg.solve(XtX, Xty)
    sigma2 = (yA @ yA + yB @ yB - Xty @ beta) / (nA + nB - XA.shape[1])
    rA, rB = yA - XA @ beta, yB - XB @ beta

    Q = []
    for g in range(len(A_blocks)):
        q = 0.0
        # PART A: public list (A + B aligned, summed column-wise)
        for k in range(A_blocks[g].shape[1]):
            score = A_blocks[g][:, k] @ rA + B_aligned[g][:, k] @ rB
            count = A_blocks[g][:, k].sum() + B_aligned[g][:, k].sum()
            w = beta_density_weight(count, two_n)
            q += w * w * score * score / (2 * sigma2)
        # PART B: private (B only)
        for k in range(B_priv[g].shape[1]):
            score = B_priv[g][:, k] @ rB
            count = B_priv[g][:, k].sum()
            w = beta_density_weight(count, two_n)
            q += w * w * score * score / (2 * sigma2)
        Q.append(float(q))
    return Q


def pooled_Q(union_blocks, XA, yA, XB, yB):
    """Ground truth: per-gene SKAT on the pooled union genotype (0-filled for party-unique).
    union_blocks[g]: (nA+nB) x m_union , rows 0..nA-1 = A, nA.. = B."""
    nA, nB = len(yA), len(yB)
    two_n = 2.0 * (nA + nB)
    X = np.vstack([XA, XB])
    y = np.concatenate([yA, yB])
    XtX = X.T @ X
    Xty = X.T @ y
    beta = np.linalg.solve(XtX, Xty)
    sigma2 = (y @ y - Xty @ beta) / (len(y) - X.shape[1])
    r = y - X @ beta
    Q = []
    for g in range(len(union_blocks)):
        G = union_blocks[g]
        q = 0.0
        for k in range(G.shape[1]):
            score = G[:, k] @ r
            count = G[:, k].sum()
            w = beta_density_weight(count, two_n)
            q += w * w * score * score / (2 * sigma2)
        Q.append(float(q))
    return Q


def build_party_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol):
    """Assemble per-gene blocks from full cohort genotype matrices.
      gene_keys[g]: ordered public-list keys (shared+public_only); priv_keys[g]: private keys (B only)
      roles: dict key -> 'shared'|'public_only'|'private'
      A_geno: nA x ?  ; B_geno: nB x ?  ; keycol[key] -> column index into A_geno/B_geno
    Returns A_blocks, B_aligned, B_priv, union_blocks."""
    nA, nB = A_geno.shape[0], B_geno.shape[0]
    A_blocks, B_aligned, B_priv, union_blocks = [], [], [], []
    for g in range(len(gene_keys)):
        pub = gene_keys[g]
        priv = priv_keys[g]
        # PART A blocks, columns in public-list order
        A_blk = np.zeros((nA, len(pub)), dtype=np.float64)
        B_alg = np.zeros((nB, len(pub)), dtype=np.float64)
        for k, key in enumerate(pub):
            A_blk[:, k] = A_geno[:, keycol[key]]               # A has the whole public list
            if roles[key] == "shared":
                B_alg[:, k] = B_geno[:, keycol[key]]           # shared -> B data; public_only -> stays 0
        # PART B block (B private)
        B_prv = np.zeros((nB, len(priv)), dtype=np.float64)
        for k, key in enumerate(priv):
            B_prv[:, k] = B_geno[:, keycol[key]]
        # pooled union (shared both; public_only A only; private B only)
        union = pub + priv
        U = np.zeros((nA + nB, len(union)), dtype=np.float64)
        for k, key in enumerate(union):
            if roles[key] in ("shared", "public_only"):
                U[:nA, k] = A_geno[:, keycol[key]]
            if roles[key] in ("shared", "private"):
                U[nA:, k] = B_geno[:, keycol[key]]
        A_blocks.append(A_blk); B_aligned.append(B_alg); B_priv.append(B_prv)
        union_blocks.append(U)
    return A_blocks, B_aligned, B_priv, union_blocks


def write_int8_block(path, mat):
    """row-major n*m int8 (matches plinkBedToBinary.py + gwas.loadDenseBlocks)."""
    np.asarray(np.rint(mat), dtype=np.int8).tofile(path)


def write_blocks_and_verify(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, XA, yA, XB, yB, out_dir):
    A_blocks, B_aligned, B_priv, union_blocks = build_party_blocks(
        gene_keys, priv_keys, roles, A_geno, B_geno, keycol)
    # plaintext check BEFORE the secure run: federated (from blocks) == pooled, per gene
    qfed = federated_Q_from_blocks(A_blocks, B_aligned, B_priv, XA, yA, XB, yB)
    qpool = pooled_Q(union_blocks, XA, yA, XB, yB)
    ok = True
    for g, (a, b) in enumerate(zip(qfed, qpool)):
        rel = abs(a - b) / max(abs(b), 1e-12)
        ok &= rel < 1e-9
        print(f"  gene{g:<3} fed={a:14.4f} pooled={b:14.4f} rel={rel:.1e}")
    assert ok, "federated(from blocks) != pooled -- alignment bug"
    if out_dir:
        os.makedirs(f"{out_dir}/A", exist_ok=True)
        os.makedirs(f"{out_dir}/B", exist_ok=True)
        for g in range(len(gene_keys)):
            write_int8_block(f"{out_dir}/A/geno.{g}.bin", A_blocks[g])
            write_int8_block(f"{out_dir}/B/geno.{g}.bin", B_aligned[g])
            write_int8_block(f"{out_dir}/B/priv.{g}.bin", B_priv[g])
        np.savetxt(f"{out_dir}/A/pheno.txt", yA); np.savetxt(f"{out_dir}/A/cov.txt", XA[:, 1:])
        np.savetxt(f"{out_dir}/B/pheno.txt", yB); np.savetxt(f"{out_dir}/B/cov.txt", XB[:, 1:])
        json.dump({"n_genes": len(gene_keys),
                   "pub_m": [len(k) for k in gene_keys],
                   "priv_m": [b.shape[1] for b in B_priv]},
                  open(f"{out_dir}/manifest.json", "w"), indent=2)
        print(f"  wrote blocks + pheno/cov + manifest -> {out_dir}")
    print("  VERIFY OK (federated == pooled)")


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


def synth_cohort(rng, n):
    cov = rng.standard_normal((n, N_COV))
    X = np.hstack([np.ones((n, 1)), cov])
    y = cov @ (0.3 + 0.2 * np.arange(N_COV)) + 1.5 * rng.standard_normal(n)
    return X, y


# ===================== selftest: synthetic, no plink2/AoU =====================

def selftest():
    """Fabricate small synthetic cohorts + genes and check build/align/Q logic."""
    rng = np.random.default_rng(SEED)
    nA, nB, n_genes, m_per_gene = 8, 6, 3, 6
    XA, yA = synth_cohort(rng, nA)
    XB, yB = synth_cohort(rng, nB)

    # one shared genotype pool keyed g<gene>_v<j>; A and B draw independent dosages
    gene_keys, priv_keys, roles_all = [], [], {}
    all_keys = []
    for g in range(n_genes):
        keys = [f"g{g}_v{j}" for j in range(m_per_gene)]
        roles, pub, priv = split_roles(rng, keys)
        gene_keys.append(pub); priv_keys.append(priv)
        roles_all.update(roles); all_keys += keys
    keycol = {k: i for i, k in enumerate(all_keys)}

    def draw(n):
        G = np.zeros((n, len(all_keys)))
        for j in range(len(all_keys)):
            af = 0.05 + 0.4 * rng.random()
            G[:, j] = (rng.random(n) < af).astype(float) + (rng.random(n) < af).astype(float)
        return G
    A_geno, B_geno = draw(nA), draw(nB)

    print("selftest (synthetic, no plink2):")
    write_blocks_and_verify(gene_keys, priv_keys, roles_all, A_geno, B_geno, keycol, XA, yA, XB, yB, out_dir=None)


# ===================== run: real AoU pgen on the workbench =====================

def sh(cmd):
    print("  $", cmd)
    subprocess.run(cmd, shell=True, check=True)


def plink_extract_to_int8(pgen_keyed, keep_file, keys_file, n, out_prefix):
    """plink2 --keep --extract -> .bed -> plinkBedToBinary -> int8; returns (n x m matrix, keys).
    m is read from the output .bim (plink2 may drop variants); plink2 keeps the pgen's variant
    order, not the --extract order, so the caller reorders by the returned keys."""
    here = os.path.dirname(os.path.abspath(__file__))
    sh(f"{PLINK2} --pfile {pgen_keyed} --keep {keep_file} --extract {keys_file} "
       f"--indiv-sort none --make-bed --out {out_prefix}")
    bim_keys = [ln.split('\t')[1] for ln in open(f"{out_prefix}.bim")]
    m = len(bim_keys)
    sh(f"python3 {here}/../plinkBedToBinary.py {out_prefix}.bed {n} {m} {out_prefix}.bin")
    g = np.fromfile(f"{out_prefix}.bin", dtype=np.int8).reshape(n, m).astype(float)
    g[g < 0] = 0  # missing -> 0 (matches Go loader)
    return g, bim_keys


def run():
    """Real prep on the workbench. Reads AoU pgen, synthesizes 2 cohorts, writes blocks."""
    rng = np.random.default_rng(SEED)
    os.makedirs(OUT_DIR, exist_ok=True)

    # (1) biallelic + chr:pos:ref:alt keys (handles the empty AoU ID column)
    keyed = f"{OUT_DIR}/{CHR}_keyed"
    sh(f"{PLINK2} --pfile {PGEN} --max-alleles 2 --min-alleles 2 "
       f"--set-all-var-ids '@:#:$r:$a' --make-pgen --out {keyed}")

    # (2) gene -> PASS variants (read keyed .pvar; map to GENCODE genes)
    genes = load_gencode_genes(GENCODE, CHR, N_GENES)
    gene_keys_all = scan_pvar_into_genes(f"{keyed}.pvar", genes)  # gene -> ordered [keys] (PASS only)

    # (3) synthesize cohorts: sample split + per-gene role split
    psam = [ln.split()[0] for ln in open(f"{keyed}.psam") if not ln.startswith("#")]
    perm = rng.permutation(len(psam))
    A_ids = [psam[i] for i in perm[:N_SUB]]
    B_ids = [psam[i] for i in perm[N_SUB:2 * N_SUB]]
    write_keep(f"{OUT_DIR}/A.keep", A_ids); write_keep(f"{OUT_DIR}/B.keep", B_ids)

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
    Ag, Ak = plink_extract_to_int8(keyed, f"{OUT_DIR}/A.keep", f"{OUT_DIR}/A_keys.txt",
                                   N_SUB, f"{OUT_DIR}/A_geno")
    Bg, Bk = plink_extract_to_int8(keyed, f"{OUT_DIR}/B.keep", f"{OUT_DIR}/B_keys.txt",
                                   N_SUB, f"{OUT_DIR}/B_geno")
    # align both cohorts into one column space (union of A,B extracted keys)
    A_geno, B_geno, keycol = merge_cohort_columns(Ag, Ak, Bg, Bk)
    # drop any keys plink2 didn't emit (monomorphic/filtered) so build never KeyErrors
    present = set(keycol)
    gene_keys = [[k for k in g if k in present] for g in gene_keys]
    priv_keys = [[k for k in p if k in present] for p in priv_keys]

    # (5) cov/pheno (synthetic for #4; oracle uses the same -> validation stays exact)
    XA, yA = synth_cohort(rng, N_SUB)
    XB, yB = synth_cohort(rng, N_SUB)

    print("run (real AoU pgen):")
    write_blocks_and_verify(gene_keys, priv_keys, roles_all, A_geno, B_geno, keycol, XA, yA, XB, yB, OUT_DIR)


# ---- small helpers used only by run() ----

def write_keep(path, ids):
    with open(path, "w") as f:
        f.write("#IID\n" + "\n".join(ids) + "\n")


def write_lines(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


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
    keys = list(dict.fromkeys(Ak + Bk))
    keycol = {k: i for i, k in enumerate(keys)}
    A = np.zeros((Ag.shape[0], len(keys))); B = np.zeros((Bg.shape[0], len(keys)))
    for j, k in enumerate(Ak):
        A[:, keycol[k]] = Ag[:, j]
    for j, k in enumerate(Bk):
        B[:, keycol[k]] = Bg[:, j]
    return A, B, keycol


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true", help="validate logic locally (no plink2/AoU)")
    ap.add_argument("--run", action="store_true", help="real prep on the workbench")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.run:
        run()
    else:
        ap.error("pass --selftest or --run")
