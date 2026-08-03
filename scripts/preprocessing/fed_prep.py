#!/usr/bin/env python3
"""Federated-private SKAT data PREP: real AoU exome pgen -> per-gene int8 genotype blocks.

Splits one AoU pgen into cohort A (public-list owner) + B (private), and each gene's variants into
shared/public_only/private. Writes the int8 blocks the Go `skat_fed` mode reads:

    A/geno.<g>.bin    cohort A, public-list variants                        [PART A]
    B/geno.<g>.bin    cohort B, public list ALIGNED (public_only cols = 0)  [PART A]
    B/priv.<g>.bin    cohort B, private variants                            [PART B]
    {A,B}/cov.txt     covariates = first N_PCS AoU ancestry PCs, geno-row order
    {A,B}/pheno.txt   phenotype (LDL), geno-row order
    manifest.json     per-gene m (public/private)

Block = row-major n*m int8 (dosage 0/1/2, missing<0 -> Go reads 0), same layout as
plinkBedToBinary.py. Samples = geno ∩ pheno(LDL) ∩ ancestry(PC, optional group), join key person_id==research_id.
Key = chr:pos:ref:alt / GRCh38 / biallelic / PASS (see .local/warning.md). Secure SKAT is
n-independent, so N_SUB subsamples to keep blocks small while m stays realistic. Block
assembly/alignment lives in skat_plain_local.py.

    python3 fed_prep.py    # real prep on the workbench (needs plink2 + AoU pgen)
"""
import json
import math
import os
import subprocess
import time

import numpy as np

from skat_plain_local import build_party_blocks

# ---- config (env-overridable; defaults = AoU Workbench Controlled-Tier layout) ----
CHR = os.environ.get("FED_CHR", "chr22")
_V9 = "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel"
PGEN = os.path.expanduser(os.environ.get("FED_PGEN", f"{_V9}/exome/pgen/exome.{CHR}"))
GENCODE = os.path.expanduser(os.environ.get("FED_GENCODE",  # public GENCODE v44 pc-gene coords, committed next to this script
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "gencode_v44_pc_genes.bed")))
ANCESTRY = os.path.expanduser(os.environ.get("FED_ANCESTRY", f"{_V9}/aux/ancestry/ancestry_preds.tsv"))
ANCESTRY_GROUP = os.environ.get("FED_ANCESTRY_GROUP", "").strip().lower()
PHENO_CSV = os.path.expanduser(os.environ.get("FED_PHENO",
    "~/workspace/gwas-data-wgs/pheno/v9_final_lipid_med_corrected_short_read_tot.csv"))
PHENO_COL = os.environ.get("FED_PHENO_COL", "LDLC_final_mgdl_6sd_masked")  # LDL (mg/dL), continuous
OUT_DIR = os.path.expanduser(os.environ.get("FED_OUT", "~/fed_prep_out"))
KEYS_PATH = os.path.expanduser(os.environ.get("FED_KEYS", "~/secure-skat/example_data/keys"))  # MPC PRG seeds (data-independent, reusable)
PORT_BASE = int(os.environ.get("FED_PORT_BASE", "22000"))  # avoid Dataproc/Hadoop ports (8020=HDFS, 8030s=YARN, ...)
CKKS_PARAMS = os.environ.get("FED_CKKS", "PN14QP438")  # PN13QP218 = ~half RAM (slots 4096 > max gene m); on RAM-tight boxes
DATA_BITS = int(os.environ.get("FED_DATABITS", "60"))  # MPC fixed-point total bits; raise if large-n aggregates overflow
FRAC_BITS = int(os.environ.get("FED_FRACBITS", "30"))  # fractional bits (integer range = DATA_BITS-FRAC_BITS)
N_PCS = int(os.environ.get("FED_NPCS", "5"))     # first N PCs from ancestry_preds as covariates (age/sex deferred)
N_SUB = os.environ.get("FED_NSUB", "5000")  # per-cohort samples; int, or "max"/"all" = full eligible split in half
N_GENES = int(os.environ.get("FED_NGENES", "20"))  # genes (spread across chrom); >= chrom total picks all
GENE_LIST = os.environ.get("FED_GENES", "").strip()  # ""=stride-pick N_GENES; "ALL"=every pc gene on CHR; else file of symbols (one per line)
ANNOT = os.path.expanduser(os.environ.get("FED_ANNOT", ""))  # vat_annotate.py table; ""=positional gene assignment, no mask
MASK = os.environ.get("FED_MASK", "pLoF;missenseLC")  # annotation classes to keep; ";"-joined, missenseLC = missense+LC
MAX_MAF = float(os.environ.get("FED_MAXMAF", "0.001"))  # All-by-All serves 1e-4/1e-3/1e-2
AF_SOURCE = os.environ.get("FED_AF", "gnomad")  # gnomad = public + portable (zero-reveal); gvs = AoU-internal
PROBES = int(os.environ.get("FED_PROBES", "0"))    # SKAT p: trace budget; exact if m_pub<=budget, else Hutchinson (0=off)
SEED = 71


def resolve_n_sub(n_eligible):
    """FED_NSUB as an int, or 'max'/'all' -> half the eligible pool (use the full cohort, split A/B)."""
    return n_eligible // 2 if str(N_SUB).strip().lower() in ("max", "all") else int(N_SUB)
FRAC_SHARED, FRAC_PUBONLY = 0.6, 0.2   # rest = private
PLINK2 = os.environ.get("PLINK2", "plink2")   # override: PLINK2=/path/to/plink2 python3 fed_prep.py
MAX_ALLELE_LEN = 1000   # --set-all-var-ids cap so long indels keep a full key. Bump if a chrom exceeds it.


def emit_timing(scope, milliseconds, kind="phase", status="done", count=1):
    """Emit one machine-readable timing record alongside the human timing tree."""
    print(f"[timing] scope={scope} parent=prep party=driver kind={kind} "
          f"status={status} milliseconds={milliseconds:.3f} count={count}", flush=True)


def sh(cmd):
    print("  $", cmd)
    subprocess.run(cmd, shell=True, check=True)


def write_int8_block(path, mat):
    """row-major n*m int8 (matches plinkBedToBinary.py + gwas.loadDenseBlocks)."""
    np.asarray(np.rint(mat), dtype=np.int8).tofile(path)


def write_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, out_dir, gene_names=None):
    A_blocks, B_aligned, B_priv, _ = build_party_blocks(
        gene_keys, priv_keys, roles, A_geno, B_geno, keycol, with_union=False)  # union is fed_compare's job
    os.makedirs(f"{out_dir}/A", exist_ok=True)
    os.makedirs(f"{out_dir}/B", exist_ok=True)
    for g in range(len(gene_keys)):
        write_int8_block(f"{out_dir}/A/geno.{g}.bin", A_blocks[g])
        write_int8_block(f"{out_dir}/B/geno.{g}.bin", B_aligned[g])
        write_int8_block(f"{out_dir}/B/priv.{g}.bin", B_priv[g])
    shared_m = [sum(1 for k in g if roles[k] == "shared") for g in gene_keys]
    pubonly_m = [sum(1 for k in g if roles[k] == "public_only") for g in gene_keys]
    priv_m = [b.shape[1] for b in B_priv]
    json.dump({"ancestry_group": ANCESTRY_GROUP or "all",
               "n_genes": len(gene_keys),
               "gene_symbols": gene_names or [],  # block order; [] on older runs
               "pub_m": [len(k) for k in gene_keys],  # public list = shared + public_only
               "priv_m": priv_m,
               "shared_m": shared_m,        # intersection (both cohorts)
               "pubonly_m": pubonly_m},     # public-list party only
              open(f"{out_dir}/manifest.json", "w"), indent=2)
    print(f"  wrote blocks + manifest -> {out_dir}")
    print(f"  variants: shared(intersection)={sum(shared_m)} public_only={sum(pubonly_m)} "
          f"private={sum(priv_m)} total={sum(shared_m) + sum(pubonly_m) + sum(priv_m)}")


def load_ancestry_pcs(path, n_pcs, ancestry_group=""):
    """Load PCs, optionally keeping one ancestry from the AoU ancestry TSV."""
    ancestry_group = ancestry_group.strip().lower()
    pcs = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rid, pca = header.index("research_id"), header.index("pca_features")
        ancestry = None
        if ancestry_group:
            for column in ("ancestry_pred_other", "ancestry_pred"):
                if column in header:
                    ancestry = header.index(column)
                    break
            if ancestry is None:
                raise ValueError("ancestry TSV has no ancestry_pred_other or ancestry_pred column")
        for ln in f:
            x = ln.rstrip("\n").split("\t")
            if ancestry is not None and x[ancestry].strip().lower() != ancestry_group:
                continue
            vals = [float(v) for v in x[pca].strip().strip("[]").split(",")]
            if len(vals) < n_pcs:
                raise ValueError(f"pca_features has {len(vals)} PCs < N_PCS={n_pcs}")
            pcs[x[rid]] = vals[:n_pcs]
    if ancestry_group and not pcs:
        raise ValueError(f"no ancestry rows matched FED_ANCESTRY_GROUP={ancestry_group}")
    return pcs


def write_cov(fam_ids, pcs, out_path):
    np.savetxt(out_path, np.asarray([pcs[sid] for sid in fam_ids]), delimiter="\t")


def write_config_helpers(gene_keys, chrom, out_dir):
    """Files the Go 'blocks' config needs: geno_block_size_file (m per gene) + snp_position_file
    (chrom<TAB>pos per variant in block order; values unused by SKAT but the file must exist)."""
    sizes = [len(g) for g in gene_keys]
    write_lines(f"{out_dir}/block_sizes.txt", [str(s) for s in sizes])
    cnum = chrom.replace("chr", "")
    with open(f"{out_dir}/pos.txt", "w") as f:
        for keys in gene_keys:
            for k in keys:
                f.write(f"{cnum}\t{k.split(':')[1]}\n")  # key = chr:pos:ref:alt
    print(f"  num_snps={sum(sizes)} -> block_sizes.txt + pos.txt in {out_dir}")
    return sum(sizes)


def write_configs(out_dir, num_snps, n_blocks, keys_path, n_a, n_b):
    """Generate the sfgwas skat_fed config (global + party 0/1/2) with paths into out_dir, so the
    data dims (num_snps/num_inds/num_covs/n_blocks) always match the data. n_blocks = genes that
    actually have variants (< N_GENES when some picked genes are empty). Run: SFGWAS_CONFIG_PATH=<out>/config."""
    cfg = f"{out_dir}/config"
    os.makedirs(cfg, exist_ok=True)
    p1, p2, p3 = PORT_BASE + 20, PORT_BASE + 40, PORT_BASE + 60  # 0-1, 0-2, 1-2 pair bases
    glob = f"""num_main_parties = 2
hub_party_id = 1
debug = false
ckks_params = "{CKKS_PARAMS}"
mpc_num_threads = 2
mpc_field_size = 256
mpc_data_bits = {DATA_BITS}
mpc_frac_bits = {FRAC_BITS}
div_sqrt_max_len = 1000000
mpc_boolean_shares = true
num_inds = [0, {n_a}, {n_b}]
num_snps = {num_snps}
num_covs = {N_PCS}
cov_all_ones = false
geno_file_format = "blocks"
geno_num_blocks = {n_blocks}
binary_pheno = false
private_pid = 2
skat_pvalue_probes = {PROBES}
rotkey_pow2only = true
skip_qc = true
skip_pca = true
use_logistic = false
binding_ipaddr = "0.0.0.0"

[servers.party0]
ipaddr = "127.0.0.1"
ports = {{party1 = "{p1}", party2 = "{p2}"}}
[servers.party1]
ipaddr = "127.0.0.1"
ports = {{party2 = "{p3}"}}
[servers.party2]
ipaddr = "127.0.0.1"
ports = {{}}
"""
    open(f"{cfg}/configGlobal.toml", "w").write(glob)
    for pid, sub in [(0, None), (1, "A"), (2, "B")]:
        loc = (f'shared_keys_path = "{keys_path}"\n'
               f'geno_block_size_file = "{out_dir}/block_sizes.txt"\n'
               f'gene_id_file = "{out_dir}/genes.txt"\n'
               f'output_dir = "{out_dir}/out/party{pid}"\n'
               f'cache_dir = "{out_dir}/cache/party{pid}"\n'
               f'local_num_threads = 4\nmemory_limit = 6000000000\n'
               f'assoc_num_blocks_parallel = 1\n')
        if sub:  # data parties only
            loc += (f'geno_binary_file_prefix = "{out_dir}/{sub}/geno"\n'
                    f'geno_num_blocks = {n_blocks}\n'
                    f'pheno_file = "{out_dir}/{sub}/pheno.txt"\n'
                    f'covar_file = "{out_dir}/{sub}/cov.txt"\n'
                    f'snp_position_file = "{out_dir}/pos.txt"\n')
        if pid == 2:
            loc += f'private_only_prefix = "{out_dir}/B/priv.%d.bin"\n'
        open(f"{cfg}/configLocal.Party{pid}.toml", "w").write(loc)
    print(f"  config -> {cfg}  (run: SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH={cfg})")


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
    run_started = time.perf_counter()
    sample_rng = np.random.default_rng(SEED)
    role_rng = np.random.default_rng(SEED + 1)
    os.makedirs(OUT_DIR, exist_ok=True)
    tmr = {}
    t = time.perf_counter()

    keyed = f"{OUT_DIR}/{CHR}_keyed"
    sh(f"{PLINK2} --pfile {PGEN} --max-alleles 2 --min-alleles 2 "
       f"--set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len {MAX_ALLELE_LEN} "
       f"--make-just-pvar --out {keyed}")
    tmr["keying"] = time.perf_counter() - t; t = time.perf_counter()

    # (2) gene -> PASS variants (read keyed .pvar; map to GENCODE genes)
    genes = load_gencode_genes(GENCODE, CHR, N_GENES)
    annot = load_annotation(ANNOT, MASK, MAX_MAF, AF_SOURCE) if ANNOT else None
    gene_names, gene_keys_all = scan_pvar_into_genes(f"{keyed}.pvar", genes, annot)
    print(f"  genes: {len(genes)} selected -> {len(gene_names)} with variants "
          f"(FED_GENES={GENE_LIST or 'stride'})")

    # (3) eligible = geno ∩ pheno(LDL) ∩ ancestry(PC); split into cohort A/B; per-gene role split
    pcs = load_ancestry_pcs(ANCESTRY, N_PCS, ANCESTRY_GROUP)
    pheno = load_pheno(PHENO_CSV, PHENO_COL)
    psam = [ln.split()[0] for ln in open(f"{PGEN}.psam") if not ln.startswith("#")]  # samples unchanged by keying
    gset = set(psam)
    eligible = [s for s in psam if s in pheno and s in pcs]
    # per-component + pairwise counts so a person_id!=research_id namespace mismatch is visible
    ancestry_label = ANCESTRY_GROUP.upper() if ANCESTRY_GROUP else "all"
    print(f"  geno={len(gset)} pheno={len(pheno)} ancestry[{ancestry_label}]={len(pcs)} | "
          f"geno∩pheno={len(gset & pheno.keys())} geno∩pc={len(gset & pcs.keys())} eligible={len(eligible)}")
    ns = resolve_n_sub(len(eligible))
    if len(eligible) < 2 * ns:
        raise SystemExit(f"only {len(eligible)} eligible samples < 2*n_sub={2 * ns} (check person_id==research_id matching)")
    # 'max'/'all' uses the whole pool; an odd leftover goes to cohort A (A = ns+1).
    a_n = len(eligible) - ns if str(N_SUB).strip().lower() in ("max", "all") else ns
    print(f"  cohorts A={a_n} B={ns} (FED_NSUB={N_SUB}) -> {a_n + ns} of {len(eligible)} eligible used")
    perm = sample_rng.permutation(len(eligible))
    A_ids = [eligible[i] for i in perm[:a_n]]
    B_ids = [eligible[i] for i in perm[a_n:a_n + ns]]
    write_lines(f"{OUT_DIR}/A.keep", A_ids, "#IID"); write_lines(f"{OUT_DIR}/B.keep", B_ids, "#IID")

    gene_keys, priv_keys, roles_all, all_keys = [], [], {}, []
    for keys in gene_keys_all:
        roles, pub, priv = split_roles(role_rng, keys)
        gene_keys.append(pub); priv_keys.append(priv)
        roles_all.update(roles); all_keys += keys
    A_extract = [k for k in all_keys if roles_all[k] in ("shared", "public_only")]
    B_extract = [k for k in all_keys if roles_all[k] in ("shared", "private")]
    write_lines(f"{OUT_DIR}/A_keys.txt", A_extract)
    write_lines(f"{OUT_DIR}/B_keys.txt", B_extract)
    tmr["setup"] = time.perf_counter() - t; t = time.perf_counter()

    # (4) extract genotypes (plink2 -> int8), reorder to a single key->column map
    Ag, Ak, Afam = plink_extract_to_int8(PGEN, f"{OUT_DIR}/A.keep", f"{OUT_DIR}/A_keys.txt",
                                         a_n, f"{OUT_DIR}/A_geno")
    tmr["extract_A"] = time.perf_counter() - t; t = time.perf_counter()
    Bg, Bk, Bfam = plink_extract_to_int8(PGEN, f"{OUT_DIR}/B.keep", f"{OUT_DIR}/B_keys.txt",
                                         ns, f"{OUT_DIR}/B_geno")
    tmr["extract_B"] = time.perf_counter() - t; t = time.perf_counter()
    A_geno, B_geno, keycol = merge_cohort_columns(Ag, Ak, Bg, Bk)
    del Ag, Bg  # free pre-merge matrices (RAM-tight workbench)
    # drop any keys plink2 didn't emit (monomorphic/filtered) so build never KeyErrors
    present = set(keycol)
    gene_keys = [[k for k in g if k in present] for g in gene_keys]
    priv_keys = [[k for k in p if k in present] for p in priv_keys]
    # PART A (public list) and PART B (private) MUST be disjoint per gene, else Q double-counts.
    for pub, priv in zip(gene_keys, priv_keys):
        assert not (set(pub) & set(priv)), "public-list and private variant sets overlap -> Q double-count"

    # (5) write genotype blocks + real covariates (5 PCs) + real LDL phenotype, all geno-row order
    print("run (real AoU pgen):")
    write_blocks(gene_keys, priv_keys, roles_all, A_geno, B_geno, keycol, OUT_DIR, gene_names)
    write_lines(f"{OUT_DIR}/genes.txt", gene_names)  # block index -> symbol; the join key for external tables
    num_snps = write_config_helpers(gene_keys, CHR, OUT_DIR)
    write_configs(OUT_DIR, num_snps, len(gene_keys), KEYS_PATH, len(Afam), len(Bfam))
    write_cov(Afam, pcs, f"{OUT_DIR}/A/cov.txt")
    write_cov(Bfam, pcs, f"{OUT_DIR}/B/cov.txt")
    write_pheno(Afam, pheno, f"{OUT_DIR}/A/pheno.txt")
    write_pheno(Bfam, pheno, f"{OUT_DIR}/B/pheno.txt")
    print(f"  wrote cov ({N_PCS} PCs) + pheno (LDL) -> {OUT_DIR}/{{A,B}}/")
    tmr["write"] = time.perf_counter() - t
    desc = {"keying": "plink2 set-var-ids + biallelic pvar",
            "setup": "gene map + eligible ∩ + role split",
            "extract_A": "plink2 extract cohort A -> int8",
            "extract_B": "plink2 extract cohort B -> int8",
            "write": "blocks + cov + pheno + config"}
    phases = ["keying", "setup", "extract_A", "extract_B", "write"]
    total = time.perf_counter() - run_started
    print("[prep] timing tree (plaintext):")
    for k in phases:
        print(f"  ├─ {k:<12} {tmr[k]:7.1f}s   ({desc[k]})")
    print(f"  └─ {'TOTAL':<12} {total:7.1f}s")
    for k in phases:
        emit_timing(f"prep.{k}", 1000.0 * tmr[k])
    emit_timing("prep.total", 1000.0 * total, kind="total")


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
    """gene -> (lo, hi) for protein-coding genes on chrom, always in genomic order.

    FED_GENES selects which: unset = stride-pick n_genes spread across the chromosome (the historical
    behaviour, positional not biological); "ALL" = every gene on chrom; else a path to a file of gene
    symbols, one per line. Genomic order is enforced in every mode so the block order -- and therefore
    genes.txt, block_sizes.txt and the result vectors -- is reproducible regardless of input order."""
    genes = []
    for ln in open(bed):
        f = ln.rstrip("\n").split("\t")
        if f[0] != chrom:
            continue
        name = f[3] if len(f) > 3 else f"{f[0]}:{f[1]}"
        genes.append((name, int(f[1]), int(f[2])))
    genes.sort(key=lambda x: x[1])
    if GENE_LIST.upper() == "ALL":
        picked = genes
    elif GENE_LIST:
        want = {ln.split("#", 1)[0].strip() for ln in open(GENE_LIST)}
        want.discard("")
        picked = [g for g in genes if g[0] in want]
        missing = want - {g[0] for g in picked}
        if missing:  # a typo'd symbol must not silently shrink the panel
            raise SystemExit(f"FED_GENES={GENE_LIST}: {len(missing)} symbol(s) not in {bed} "
                             f"for {chrom}: {sorted(missing)[:10]}")
    else:
        step = max(1, len(genes) // n_genes)
        picked = genes[::step][:n_genes]
    return {name: (lo, hi) for name, lo, hi in picked}


def load_annotation(path, mask, max_maf, af_source):
    """variant_key -> [gene_symbol, ...], keeping only rows that pass the mask and the MAF cut.

    The table comes from vat_annotate.py, one row per (variant, gene): a variant's consequence is
    gene-specific, so the same variant can be missense for one gene and upstream of its neighbour.
    A "mask" is a SELECTION over annotation classes -- All-by-All serves pLoF, missenseLC,
    pLoF;missenseLC and synonymous, where missenseLC means missense or low-confidence LoF.

    Variants absent from the AF source are KEPT: absence from gnomAD means ultra-rare, and dropping
    them would cut exactly the tail SKAT weights hardest. Counted and reported, not silent."""
    classes = {"missenseLC": ("missense", "LC")}
    wanted = set()
    for part in mask.split(";"):
        part = part.strip()
        wanted.update(classes.get(part, (part,)))
    af_col = {"gnomad": 4, "gvs": 5}[af_source]

    key2genes, n_rows, drop_mask, drop_maf, no_af = {}, 0, 0, 0, 0
    with open(path) as f:
        f.readline()  # header
        for ln in f:
            r = ln.rstrip("\n").split("\t")
            n_rows += 1
            if r[3] not in wanted:
                drop_mask += 1
                continue
            raw = r[af_col].strip()
            if not raw:
                no_af += 1
            elif float(raw) > max_maf:
                drop_maf += 1
                continue
            key2genes.setdefault(r[0], []).append(r[2])  # gene_symbol; gene_id kept in the file
    print(f"  annotation: {n_rows:,} rows -> {len(key2genes):,} variants pass "
          f"mask={mask} maf<={max_maf} af={af_source}")
    print(f"    dropped {drop_mask:,} on annotation, {drop_maf:,} on MAF; "
          f"{no_af:,} kept with no {af_source} AF (treated as rare)")
    return key2genes


def scan_pvar_into_genes(pvar, genes, annot=None):
    """assign PASS biallelic keys to genes; return (symbols, key-lists), both in block order with
    genes that got no variants dropped from both. The symbols are the only record of which gene a
    block is -- nothing downstream carries them, so they must be written out here.

    With `annot`, a variant joins gene G iff the VAT annotated it FOR G with a wanted class. That is
    the assignment All-by-All uses, and the only correct one where genes overlap: variants routinely
    carry annotations for more than one gene, so first-match-by-position would file one under a
    neighbouring gene and then read off that gene's (usually 'other') consequence. Without `annot`,
    fall back to the positional rule."""
    buckets = {name: [] for name in genes}
    n_var = n_placed = n_unannotated = 0
    for ln in open(pvar):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        pos, vid, ref, alt, flt = int(f[1]), f[2], f[3], f[4], f[5]
        if flt not in ("PASS", "."):
            continue
        if "," in alt:        # multiallelic already excluded by --max-alleles, belt-and-suspenders
            continue
        n_var += 1
        if annot is None:
            for name, (lo, hi) in genes.items():
                if lo <= pos < hi:
                    buckets[name].append(vid)  # vid is the key set by --set-all-var-ids
                    break
        else:
            # vat_annotate.py keys on pos:ref:alt -- no contig, because plink2 and the VAT disagree
            # on the 'chr' prefix and that mismatch alone once made the join return nothing.
            hit = annot.get(f"{pos}:{ref}:{alt}")
            if hit is None:
                n_unannotated += 1
                continue
            for name in hit:
                if name in buckets:
                    buckets[name].append(vid)
                    n_placed += 1
    if annot is not None:
        empty = [g for g in genes if not buckets[g]]
        print(f"  assignment: {n_var:,} pvar variants -> {n_placed:,} gene slots "
              f"({n_unannotated:,} carried no wanted annotation)")
        if empty:
            print(f"    {len(empty)} of {len(genes)} selected genes got nothing"
                  f"{' (symbol mismatch?)' if len(empty) > len(genes) // 2 else ''}: "
                  f"{empty[:8]}{' ...' if len(empty) > 8 else ''}")
    kept = [(name, v) for name, v in buckets.items() if v]
    return [name for name, _ in kept], [v for _, v in kept]


def merge_cohort_columns(Ag, Ak, Bg, Bk):
    """put A and B genotype matrices into one shared column space keyed by variant key."""
    assert len(set(Ak)) == len(Ak) and len(set(Bk)) == len(Bk), "duplicate variant key in a cohort .bim"
    keys = list(dict.fromkeys(Ak + Bk))
    keycol = {k: i for i, k in enumerate(keys)}
    A = np.zeros((Ag.shape[0], len(keys)), dtype=Ag.dtype)
    B = np.zeros((Bg.shape[0], len(keys)), dtype=Bg.dtype)
    A[:, np.fromiter((keycol[k] for k in Ak), int, len(Ak))] = Ag
    B[:, np.fromiter((keycol[k] for k in Bk), int, len(Bk))] = Bg
    return A, B, keycol


def check():
    """Dry-run: report sizes (sample N, eligible intersection, per-gene variant counts, A/B split)
    from the real files WITHOUT plink2 extraction or the secure run. Sanity-check before a full run."""
    geno = {ln.split()[0] for ln in open(f"{PGEN}.psam") if not ln.startswith("#")}
    anc = set(load_ancestry_pcs(ANCESTRY, N_PCS, ANCESTRY_GROUP))
    phe = set(load_pheno(PHENO_CSV, PHENO_COL))
    elig = geno & anc & phe
    ancestry_label = ANCESTRY_GROUP.upper() if ANCESTRY_GROUP else "all"
    print(f"samples: geno={len(geno)}  ancestry[{ancestry_label}](PC)={len(anc)}  pheno(LDL)={len(phe)}")
    print(f"  pairwise: geno∩pheno={len(geno & phe)}  geno∩anc={len(geno & anc)}  anc∩pheno={len(anc & phe)}")
    ns = resolve_n_sub(len(elig))
    print(f"  eligible (geno∩anc∩pheno) = {len(elig)}   (n_sub={ns}/cohort, need ≥ 2*n_sub = {2 * ns})")
    if len(elig) < 2 * ns:
        print("  !! eligible < 2*n_sub → lower FED_NSUB, or person_id!=research_id namespace mismatch")

    # Run the SAME assignment the real prep runs, so --check previews what will actually happen.
    # (It used to count positionally and ignore FED_ANNOT entirely, which reported the unfiltered
    # ~600k and made a masked run look like it had changed nothing.)
    genes = load_gencode_genes(GENCODE, CHR, N_GENES)
    print()
    annot = load_annotation(ANNOT, MASK, MAX_MAF, AF_SOURCE) if ANNOT else None
    names, keys = scan_pvar_into_genes(f"{PGEN}.pvar", genes, annot)
    counts = dict(zip(names, (len(k) for k in keys)))
    tot = sum(counts.values())
    how = f"mask={MASK} maf<={MAX_MAF}" if annot else "NO mask (FED_ANNOT unset)"
    print(f"\nvariants (PASS biallelic, {CHR}, {how}): {len(counts)} genes with data, total m={tot:,}")

    sizes = sorted(counts.values())
    if sizes:
        q = lambda p: sizes[min(len(sizes) - 1, int(p * len(sizes)))]
        print(f"  per-gene m: min={sizes[0]} p25={q(.25)} median={q(.5)} p75={q(.75)} max={sizes[-1]}"
              f"   mean={tot / len(sizes):.0f}")
        big = sum(1 for s in sizes if s > PROBES) if PROBES else 0
        if PROBES:
            print(f"  m > FED_PROBES({PROBES}): {big}/{len(sizes)} genes need Hutchinson, "
                  f"{len(sizes) - big} get the exact trace")
    pub_tot = sum(round(FRAC_SHARED * m) + round(FRAC_PUBONLY * m) for m in counts.values())
    priv_tot = tot - pub_tot
    print(f"  => A public-list m={pub_tot:,}   B private m={priv_tot:,}   "
          f"(synthetic {FRAC_SHARED:.0%}/{FRAC_PUBONLY:.0%}/rest split)")

    top = sorted(counts.items(), key=lambda kv: -kv[1])[:15]
    print(f"  largest genes: " + ", ".join(f"{n}={m}" for n, m in top))


def _selfcheck():
    """Gene selection + block/symbol alignment, on the committed GENCODE BED. No AoU data needed:
    python3 scripts/preprocessing/fed_prep.py --selfcheck

    ponytail: temporary scaffold for the FED_GENES / genes.txt change -- delete once the mask work
    lands and a real run has confirmed the block-to-symbol alignment on AoU."""
    import tempfile
    global GENE_LIST
    saved = GENE_LIST
    tmp = tempfile.mkdtemp()

    GENE_LIST = ""
    stride = load_gencode_genes(GENCODE, "chr22", 20)
    assert len(stride) == 20, len(stride)
    assert list(stride)[:3] == ["OR11H1", "FAM246B", "RTL10"], list(stride)[:3]  # historical panel unchanged

    GENE_LIST = "ALL"
    every = load_gencode_genes(GENCODE, "chr22", 20)
    assert len(every) == 447, len(every)
    assert list(every) == sorted(every, key=lambda g: every[g][0]), "ALL must be in genomic order"

    # file order deliberately reversed: output must still be genomic order, else block indices
    # would depend on how someone happened to type the list
    path = f"{tmp}/genes.txt"
    open(path, "w").write("PPARA\nSREBF2\n# a comment\nCSF2RB\n\nIL17RA\n")
    GENE_LIST = path
    picked = load_gencode_genes(GENCODE, "chr22", 20)
    assert list(picked) == ["IL17RA", "CSF2RB", "SREBF2", "PPARA"], list(picked)

    open(path, "w").write("IL17RA\nNOT_A_REAL_GENE\n")
    try:
        load_gencode_genes(GENCODE, "chr22", 20)
        raise AssertionError("a bogus symbol must abort, not silently shrink the panel")
    except SystemExit:
        pass

    # Fixtures below use round, obviously-synthetic coordinates. Anything shaped like a real
    # chr:pos:ref:alt is indistinguishable from an AoU callset variant to a reader, and a variant
    # key is Controlled Tier -- it says that allele exists in the cohort. Keep test data unmistakable.
    G1, G2, G3 = "TESTGENE1", "TESTGENE2", "TESTGENE3"
    genes = {G1: (1000, 2000), G2: (2000, 3000), G3: (3000, 4000)}

    # symbols must stay aligned with key-lists, and empty genes must drop from BOTH
    pvar = f"{tmp}/test.pvar"
    open(pvar, "w").write(
        "#CHROM\tPOS\tID\tREF\tALT\tFILTER\n"
        "chrTEST\t1100\tv1\tG\tA\tPASS\n"      # in G1
        "chrTEST\t3100\tv2\tC\tT\tPASS\n"      # in G3
        "chrTEST\t3200\tv3\tA\tG\tFAIL\n")     # in G3 but filtered out
    names, keys = scan_pvar_into_genes(pvar, genes)
    assert names == [G1, G3], names                      # G2 had no variants -> dropped
    assert len(names) == len(keys), (len(names), len(keys))
    assert [len(k) for k in keys] == [1, 1], keys        # the FAIL row is excluded

    # run() filters again against the keys plink2 actually emitted. That comprehension is what keeps
    # gene_names aligned, so pin it here: length preserved, order preserved, a gene may go empty but
    # must NOT vanish -- if it ever starts vanishing, gene_names has to be filtered in lockstep.
    present = {keys[1][0]}                               # pretend plink2 dropped G1's only variant
    filtered = [[k for k in g if k in present] for g in keys]
    assert len(filtered) == len(names), (len(filtered), len(names))
    assert [len(g) for g in filtered] == [0, 1], filtered
    assert names[1] == G3 and filtered[1] == keys[1], "block->symbol alignment shifted"

    # annotation-driven assignment (option B): the gene comes from the VAT row, not from coordinates
    ann = f"{tmp}/ann.tsv"
    open(ann, "w").write(
        "variant_key\tgene_id\tgene_symbol\tannotation\tgnomad_af\tgvs_af\n"
        # AFs are 0 (passes any cutoff) or 9 (fails any cutoff) so the fixture cannot be mistaken
        # for real frequencies, and the test does not depend on the cutoff's value.
        f"1100:G:A\tENSG_A\t{G1}\tpLoF\t0\t0\n"
        f"3100:C:T\tENSG_B\t{G3}\tmissense\t0\t0\n"
        f"3100:C:T\tENSG_C\t{G2}\tother\t0\t0\n"       # same variant, neighbour gene
        f"3200:A:G\tENSG_B\t{G3}\tsynonymous\t0\t0\n"  # wrong class for our mask
        f"3300:T:C\tENSG_B\t{G3}\tpLoF\t9\t9\n")       # right class, over any cutoff
    a = load_annotation(ann, "pLoF;missenseLC", 0.001, "gnomad")
    assert a == {"1100:G:A": [G1], "3100:C:T": [G3]}, a

    pvar2 = f"{tmp}/t2.pvar"
    open(pvar2, "w").write(
        "#CHROM\tPOS\tID\tREF\tALT\tFILTER\n"
        "chrTEST\t1100\tv1\tG\tA\tPASS\n"
        "chrTEST\t3100\tv2\tC\tT\tPASS\n"
        "chrTEST\t3200\tv3\tA\tG\tPASS\n"
        "chrTEST\t3300\tv4\tT\tC\tPASS\n"
        "chrTEST\t3900\tv5\tG\tC\tPASS\n")  # inside G3's span, but unannotated
    # G3's span contains four of these; positional assignment would take all four, and would also
    # hand 3100 to whichever of G2/G3 sorts first. Annotation mode must not.
    names2, keys2 = scan_pvar_into_genes(pvar2, genes, a)
    assert names2 == [G1, G3], names2
    assert keys2 == [["v1"], ["v2"]], keys2

    GENE_LIST = saved
    print("fed_prep selfcheck OK")


if __name__ == "__main__":
    import sys
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        check() if "--check" in sys.argv else run()
