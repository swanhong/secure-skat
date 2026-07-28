#!/usr/bin/env python3
"""Annotate our exome variants from the All of Us Variant Annotation Table (VAT).

Produces the lookup table that fed_prep.py filters with, reproducing the annotation classes the
All-by-All release uses (github.com/atgu/aou_gwas, utils/annotations.py):

    LoF == 'HC'                          -> pLoF
    LoF == 'LC'                          -> LC
    consequence == 'missense_variant'    -> missense
    consequence == 'synonymous_variant'  -> synonymous
    anything else                        -> other
    (a "mask" such as pLoF;missenseLC is a SELECTION over these, made at filter time)

No CADD/REVEL/SIFT/PolyPhen: All-by-All does not use them either.

Output columns, one row per (variant, gene) -- a variant overlapping two genes gets two rows,
because its consequence is gene-specific:

    variant_key  gene_id  gene_symbol  annotation  gnomad_af  gvs_af

Run:
    GOOGLE_CLOUD_PROJECT=... python3 vat_annotate.py <chr22_keyed.pvar> [out.tsv]

"""
import gzip
import os
import struct
import subprocess
import sys
import time
import zlib
from collections import Counter

VATDIR = os.environ.get(
    "FED_VATDIR",
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/vat")
VAT = f"{VATDIR}/vat_complete.bgz.tsv.gz"
CHROM = os.environ.get("FED_CHR", "chr22")
PROJ = os.environ.get("GOOGLE_CLOUD_PROJECT") or os.environ.get("GOOGLE_PROJECT")
CHUNK = int(os.environ.get("VAT_CHUNK_BYTES", 64_000_000))  # per gsutil range request
PROGRESS_EVERY = 2_000_000  # rows

# 0-based column indexes for the 114-column VAT header
COL_CONTIG, COL_POS, COL_REF, COL_ALT = 2, 3, 4, 5
COL_GVS_AF, COL_GENE_SYMBOL, COL_CONSEQUENCE = 8, 43, 46
COL_GENE_ID, COL_CANONICAL, COL_GNOMAD_AF, COL_LOF = 53, 55, 56, 109
NCOL_MIN = 110

# most severe first; used to collapse a variant's transcripts down to one call
SEVERITY = ["pLoF", "LC", "missense", "synonymous", "other"]
PSEUDO_BIN = 37450  # tabix stashes record COUNTS here, not file offsets -- never read it as one
LINEAR_SHIFT = 14   # linear index granularity: one entry per 16 kb


def gs(args):
    return subprocess.run(["gsutil"] + (["-u", PROJ] if PROJ else []) + args,
                          capture_output=True, check=True).stdout


def chrom_offset(tbi_bytes, chrom, pos=None):
    """Byte offset to start reading at. pos=None -> where the contig begins."""
    d = gzip.decompress(tbi_bytes)
    if d[:4] != b"TBI\x01":
        raise SystemExit(f"not a tabix index (magic={d[:4]!r})")
    n_ref = struct.unpack_from("<i", d, 4)[0]
    l_nm = struct.unpack_from("<i", d, 32)[0]
    names = [n.decode() for n in d[36:36 + l_nm].split(b"\x00") if n]
    if chrom not in names:
        raise SystemExit(f"{chrom} not in index; names seen: {names[:8]}")
    target = names.index(chrom)

    p = 36 + l_nm
    for ref in range(n_ref):
        lo = None
        n_bin = struct.unpack_from("<i", d, p)[0]; p += 4
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", d, p); p += 8
            for _ in range(n_chunk):
                beg = struct.unpack_from("<Q", d, p)[0]; p += 16
                if bin_id != PSEUDO_BIN:
                    off = beg >> 16
                    lo = off if lo is None else min(lo, off)
        n_intv = struct.unpack_from("<i", d, p)[0]; p += 4
        intv = list(struct.unpack_from(f"<{n_intv}Q", d, p)); p += 8 * n_intv
        if ref == target:
            if pos is not None:
                i = min(pos >> LINEAR_SHIFT, n_intv - 1)
                while i >= 0 and (intv[i] >> 16) == 0:
                    i -= 1
                if i < 0:
                    raise SystemExit(f"no linear-index entry at or before {pos:,} on {chrom}")
                return intv[i] >> 16
            for ioff in intv:
                off = ioff >> 16
                if off:
                    lo = off if lo is None else min(lo, off)
            if lo is None:
                raise SystemExit(f"{chrom} has no data chunks in the index")
            return lo
    raise SystemExit("index truncated before reaching the target contig")


def stream_lines(start):
    """Yield decoded lines of the VAT from byte `start`, at bounded memory.

    Reads gsutil's stdout in small blocks and feeds a zlib decompressor incrementally, so nothing
    larger than one block is ever held. (Buffering a whole range request instead needs ~15x its size
    in RAM once decompressed and split -- at a 200 MB range that is several GB, which is what got
    this script OOM-killed on the Dataproc master.)

    BGZF is concatenated gzip members: when one ends, `unused_data` holds the start of the next, so
    a fresh decompressor picks up there. Member and line boundaries both survive chunk edges because
    the decompressor object and the trailing partial line carry across iterations."""
    pos, dec, carry, nbytes = start, zlib.decompressobj(31), "", 0
    while True:
        cmd = ["gsutil"] + (["-u", PROJ] if PROJ else []) + \
              ["cat", "-r", f"{pos}-{pos + CHUNK}", VAT]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                bufsize=1 << 20)
        got = 0
        while True:
            buf = proc.stdout.read(1 << 20)
            if not buf:
                break
            got += len(buf)
            while buf:
                try:
                    out = dec.decompress(buf)
                except zlib.error:
                    proc.kill()
                    return
                if out:
                    lines = (carry + out.decode("utf-8", "replace")).split("\n")
                    carry = lines.pop()
                    yield from lines
                if dec.eof:                      # member finished mid-buffer
                    buf = dec.unused_data
                    dec = zlib.decompressobj(31)
                else:
                    buf = b""
        proc.stdout.close()
        proc.wait()
        nbytes += got
        if got == 0:                             # past the end of the object
            break
        pos += CHUNK + 1
    if carry:
        yield carry


def classify(lof, consequence):
    """The All-by-All annotation classes. consequence may be a comma-joined VEP list."""
    if lof == "HC":
        return "pLoF"
    if lof == "LC":
        return "LC"
    terms = consequence.split(",")
    if "missense_variant" in terms:
        return "missense"
    if "synonymous_variant" in terms:
        return "synonymous"
    return "other"


def main():
    if not PROJ:
        raise SystemExit("set GOOGLE_CLOUD_PROJECT (or GOOGLE_PROJECT) -- buckets are requester-pays")
    if not sys.argv[1:]:
        raise SystemExit(__doc__)
    pvar, out_path = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else
                                   os.path.expanduser(f"~/fed_prep_out/{CHROM}_annotation.tsv"))

    # Our variants define the work: annotate only what we can actually test, and the leftovers
    # become the `unannotated` count rather than a silent loss.
    ours, lo, hi = set(), None, 0
    for ln in open(pvar):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if f[5] in ("PASS", ".") and "," not in f[4]:
            p = int(f[1])
            ours.add(f"{p}:{f[3]}:{f[4]}")  # contig dropped: normalised on both sides
            lo = p if lo is None else min(lo, p)
            hi = max(hi, p)
    if not ours:
        raise SystemExit(f"no PASS biallelic variants in {pvar}")
    print(f"[1/3] {len(ours):,} of our variants, {CHROM}:{lo:,}-{hi:,}")

    print(f"[2/3] streaming VAT from the tabix offset for {CHROM}:{lo:,}")
    off = chrom_offset(gs(["cat", f"{VAT}.tbi"]), CHROM, lo)
    best, af, seen_rows, kept_rows = {}, {}, 0, 0
    tick = time.time()
    for ln in stream_lines(off):
        # Cheap rejection first: only the leading 6 fields are needed to decide, and str.split with
        # maxsplit stops scanning there. Splitting all 114 columns for every row -- the vast
        # majority of which are variants we do not carry -- is what makes this slow.
        head = ln.split("\t", 6)
        if len(head) < 7 or head[COL_CONTIG] != CHROM:
            if seen_rows:
                break            # walked off the end of our contig
            continue             # still in the previous contig's tail
        seen_rows += 1
        if seen_rows % 2_000_000 == 0:
            print(f"      ... {seen_rows / 1e6:.0f}M rows, {kept_rows:,} matched, "
                  f"{time.time() - tick:.0f}s", flush=True)
        pos = int(head[COL_POS])
        if pos > hi:
            break
        key = f"{pos}:{head[COL_REF]}:{head[COL_ALT]}"
        if key not in ours:
            continue
        r = ln.split("\t")       # only now is the full row worth parsing
        if len(r) < NCOL_MIN:
            continue
        kept_rows += 1
        gene = r[COL_GENE_ID] or r[COL_GENE_SYMBOL]
        ann = classify(r[COL_LOF].strip(), r[COL_CONSEQUENCE].strip())
        canon = r[COL_CANONICAL].strip().lower() in ("true", "t", "1")
        # canonical wins outright; among equals, the most severe consequence wins. Nearly every
        # variant we test has a canonical row, so the fallback is an edge case, not the common path.
        rank = (0 if canon else 1, SEVERITY.index(ann))
        cur = best.get((key, gene))
        if cur is None or rank < cur[0]:
            best[(key, gene)] = (rank, ann, r[COL_GENE_SYMBOL])
        af.setdefault(key, (r[COL_GNOMAD_AF].strip(), r[COL_GVS_AF].strip()))
    print(f"      {seen_rows:,} {CHROM} rows scanned, {kept_rows:,} matched our variants "
          f"({time.time() - tick:.0f}s)")

    print(f"[3/3] writing {out_path}")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        f.write("variant_key\tgene_id\tgene_symbol\tannotation\tgnomad_af\tgvs_af\n")
        for (key, gene), (_, ann, sym) in sorted(best.items()):
            g_af, v_af = af.get(key, ("", ""))
            f.write(f"{key}\t{gene}\t{sym}\t{ann}\t{g_af}\t{v_af}\n")

    annotated = {k for k, _ in best}
    miss = len(ours) - len(annotated)
    print(f"\n      annotated {len(annotated):,}/{len(ours):,} "
          f"({100 * len(annotated) / len(ours):.1f}%)   unannotated {miss:,}")
    hist = Counter(v[1] for v in best.values())
    for a in SEVERITY:
        print(f"        {a:<11} {hist.get(a, 0):>8,}")
    multi = len(best) - len(annotated)
    if multi:
        print(f"      {multi:,} extra rows from variants overlapping more than one gene")
    if miss > len(ours) * 0.1:
        print("\n      !! over 10% unannotated -- do not filter on this table until that is explained")


if __name__ == "__main__":
    main()
