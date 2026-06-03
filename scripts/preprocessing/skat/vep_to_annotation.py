#!/usr/bin/env python3
"""Derive a public functional-annotation variant list from a VEP table.

This is the annotation plug-in for the AoU-private hidden-tail selection (paper
Section 7). The fixture/selection step treats the *consequence* annotation as a
public per-position property, so a variant is kept as a functional rare-variant
candidate when its best VEP annotation is protein-coding AND its consequence is
in the allowed loss-of-function / missense set. The same rule is used by
``scripts/preprocessing/build_anchor_window_dataset.py`` (kept in sync here).

Output: one variant ID (``chrom:pos:ref:alt``) per line, suitable for the
``--annotation-file`` option of ``build_1000g_mvp_public_fixture.py``. NOTE: the
emitted IDs must match the dataset's PVAR ``ID`` column convention; when the real
AoU/MVP LoF annotations are available they plug in here without touching the
secure selection logic.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

# Mirrors build_anchor_window_dataset.py ALLOWED_CONSEQUENCES (LoF + missense).
ALLOWED_CONSEQUENCES = {
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "splice_donor_variant",
    "splice_acceptor_variant",
}


def consequence_terms(raw_value: str) -> set[str]:
    return {p.strip() for p in re.split(r"[,&]", raw_value or "") if p.strip()}


def normalize_uploaded_variation(uploaded: str) -> str | None:
    toks = uploaded.split("_", 2)
    if len(toks) != 3 or "/" not in toks[2]:
        return None
    ref, alt = toks[2].split("/", 1)
    if not ref or not alt:
        return None
    return f"{toks[0].removeprefix('chr')}:{toks[1]}:{ref}:{alt}"


def is_functional(row: dict[str, str]) -> bool:
    """A variant row passes the public functional pre-filter."""
    if row.get("BIOTYPE") != "protein_coding":
        return False
    return bool(consequence_terms(row.get("Consequence", "")) & ALLOWED_CONSEQUENCES)


def functional_variants_from_vep(vep_path: Path) -> set[str]:
    """Functional variant IDs (chrom:pos:ref:alt) from a (optionally gzipped) VEP table."""
    opener = gzip.open if str(vep_path).endswith(".gz") else open
    out: set[str] = set()
    with opener(vep_path, "rt") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("##")), delimiter="\t")
        key = reader.fieldnames[0]
        for row in reader:
            vid = normalize_uploaded_variation(row[key])
            if vid is not None and is_functional(row):
                out.add(vid)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vep", required=True, help="VEP table (tab-delimited, optionally .gz)")
    ap.add_argument("--out", required=True, help="Output functional variant-ID list")
    args = ap.parse_args()
    ids = sorted(functional_variants_from_vep(Path(args.vep)))
    Path(args.out).write_text("".join(f"{v}\n" for v in ids))
    print(f"Wrote {len(ids)} functional variant IDs to {args.out}")


if __name__ == "__main__":
    main()
