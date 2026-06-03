#!/usr/bin/env python3
"""Derive a public functional-annotation variant list from a VEP table.

This is the annotation plug-in for the AoU-private hidden-tail selection (paper
Section 7). The fixture/selection step treats the *consequence* annotation as a
public per-position property, so a variant is kept as a functional rare-variant
candidate when its annotation is protein-coding AND its consequence is in the
allowed loss-of-function / missense set.

A variant may have several transcript annotations; ``--rule`` chooses how to
collapse them. Default ``canonical`` picks the single best-ranked annotation
(protein_coding > CANONICAL > allowed > severity) and tests it — matching
``scripts/preprocessing/build_anchor_window_dataset.py``. ``--rule any`` keeps a
variant if ANY transcript is functional (more sensitive). New rules go in
``functional_variants_from_vep``.

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

# Lower = more severe (mirrors build_anchor_window_dataset.py).
CONSEQUENCE_SEVERITY = {
    "transcript_ablation": 0,
    "splice_acceptor_variant": 1,
    "splice_donor_variant": 1,
    "stop_gained": 2,
    "frameshift_variant": 3,
    "stop_lost": 4,
    "start_lost": 5,
    "inframe_insertion": 6,
    "inframe_deletion": 6,
    "missense_variant": 7,
}

# Per-variant rule for collapsing multiple transcript annotations to one decision.
# Default "canonical" matches build_anchor_window_dataset.py (pick the best annotation
# by rank, then test it). Switch with --rule; add new rules in apply_rule().
RULES = ("canonical", "any")


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
    """True if this single transcript row is protein-coding with an allowed consequence."""
    if row.get("BIOTYPE") != "protein_coding":
        return False
    return bool(consequence_terms(row.get("Consequence", "")) & ALLOWED_CONSEQUENCES)


def annotation_rank(row: dict[str, str]) -> tuple[int, int, int, int]:
    """Rank to pick a variant's representative annotation: protein_coding > CANONICAL >
    allowed-consequence > severity (mirrors build_anchor_window_dataset.py)."""
    terms = consequence_terms(row.get("Consequence", ""))
    best_sev = min((CONSEQUENCE_SEVERITY.get(t, 999) for t in terms), default=999)
    return (
        int(row.get("BIOTYPE") == "protein_coding"),
        int(row.get("CANONICAL") == "YES"),
        int(bool(terms & ALLOWED_CONSEQUENCES)),
        -best_sev,
    )


def functional_variants_from_vep(vep_path: Path, rule: str = "canonical") -> set[str]:
    """Functional variant IDs (chrom:pos:ref:alt) from a (optionally gzipped) VEP table.

    rule="canonical": keep the variant iff its best-ranked annotation is functional
                      (canonical/representative transcript; matches the anchor workflow).
    rule="any":       keep the variant iff ANY transcript annotation is functional
                      (more sensitive; e.g. functional in a non-canonical transcript).
    """
    if rule not in RULES:
        raise ValueError(f"rule must be one of {RULES}, got {rule!r}")
    opener = gzip.open if str(vep_path).endswith(".gz") else open
    best: dict[str, dict[str, str]] = {}   # vid -> best-ranked row (for "canonical")
    any_func: set[str] = set()             # vids with at least one functional row (for "any")
    with opener(vep_path, "rt") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("##")), delimiter="\t")
        key = reader.fieldnames[0]
        for row in reader:
            vid = normalize_uploaded_variation(row[key])
            if vid is None:
                continue
            if is_functional(row):
                any_func.add(vid)
            prev = best.get(vid)
            if prev is None or annotation_rank(row) > annotation_rank(prev):
                best[vid] = row
    if rule == "any":
        return set(any_func)
    return {vid for vid, row in best.items() if is_functional(row)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vep", required=True, help="VEP table (tab-delimited, optionally .gz)")
    ap.add_argument("--out", required=True, help="Output functional variant-ID list")
    ap.add_argument("--rule", choices=RULES, default="canonical",
                    help="Per-variant transcript-collapse rule (default: canonical).")
    args = ap.parse_args()
    ids = sorted(functional_variants_from_vep(Path(args.vep), rule=args.rule))
    Path(args.out).write_text("".join(f"{v}\n" for v in ids))
    print(f"Wrote {len(ids)} functional variant IDs ({args.rule} rule) to {args.out}")


if __name__ == "__main__":
    main()
