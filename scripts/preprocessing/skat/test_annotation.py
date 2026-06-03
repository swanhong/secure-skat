#!/usr/bin/env python3
"""Tests for the public functional-annotation pre-filter used in AoU hidden-tail
selection (paper Section 7). Run: python3 scripts/preprocessing/skat/test_annotation.py
(also pytest-compatible)."""
from __future__ import annotations

import tempfile
from pathlib import Path

from build_1000g_mvp_public_fixture import Variant, choose_blocks
from vep_to_annotation import functional_variants_from_vep, normalize_uploaded_variation


def test_vep_functional_derivation():
    rows = [
        "#Uploaded_variation\tConsequence\tBIOTYPE\tCANONICAL",
        "22_100_A/G\tmissense_variant\tprotein_coding\tYES",                    # keep
        "22_200_C/T\tstop_gained\tprotein_coding\tYES",                         # keep (LoF)
        "22_300_G/A\tsynonymous_variant\tprotein_coding\tYES",                  # drop (synonymous)
        "22_400_T/C\tmissense_variant\tlincRNA\tNO",                            # drop (non protein_coding)
        "22_500_A/T\tsplice_donor_variant&intron_variant\tprotein_coding\tYES",  # keep (combo term)
        "chr22_600_A/G\tframeshift_variant\tprotein_coding\tNO",                # keep (chr prefix)
        # multi-transcript X: canonical synonymous (best by rank) + non-canonical missense
        "22_700_A/G\tsynonymous_variant\tprotein_coding\tYES",
        "22_700_A/G\tmissense_variant\tprotein_coding\tNO",
    ]
    with tempfile.NamedTemporaryFile("w", suffix=".tab", delete=False) as f:
        f.write("\n".join(rows) + "\n")
        path = Path(f.name)
    try:
        canon = functional_variants_from_vep(path, rule="canonical")
        anyf = functional_variants_from_vep(path, rule="any")
    finally:
        path.unlink()
    base = {"22:100:A:G", "22:200:C:T", "22:500:A:T", "22:600:A:G"}
    # canonical: X(700) dropped because its best-ranked (canonical) annotation is synonymous.
    assert canon == base, canon
    # any: X(700) kept because a non-canonical transcript is missense.
    assert anyf == base | {"22:700:A:G"}, anyf
    assert normalize_uploaded_variation("22_100_A/G") == "22:100:A:G"
    assert normalize_uploaded_variation("bad") is None


def test_choose_blocks_annotation_gating():
    """An AoU-rare variant enters the hidden tail H only if it passes the public
    annotation pre-filter (and is not in the public MVP list V_M)."""
    def V(p, name):
        return Variant(chrom="22", pos=p, vid=name, ref="A", alt="G")

    variants = {n: V(p, n) for p, n in [(10, "ov"), (20, "mo"), (30, "hA"), (40, "hX")]}
    # total_alleles=256, maf_ub=0.01 => rare iff ALT count in {1,2}.
    ac_a = {"ov": 1, "mo": 50, "hA": 1, "hX": 1}     # ov,hA,hX AoU-rare
    ac_m = {"ov": 2, "mo": 1, "hA": 50, "hX": 50}    # ov,mo MVP-rare (=> V_M); hA,hX not in V_M
    blocks = choose_blocks(variants, ac_a, ac_m, total_alleles=256, maf_ub=0.01,
                           annotation_set={"hA"}, window_bp=50000, num_blocks=1,
                           max_overlap=16, max_mvp_only=16, max_aou_hidden=16)
    b = blocks[0]
    assert sorted(v.vid for v in b.aou_hidden) == ["hA"]   # hX dropped by annotation
    assert sorted(v.vid for v in b.overlap) == ["ov"]
    assert sorted(v.vid for v in b.mvp_only) == ["mo"]


if __name__ == "__main__":
    test_vep_functional_derivation()
    test_choose_blocks_annotation_gating()
    print("all annotation tests passed")
