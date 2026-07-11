#!/usr/bin/env python3
"""Manhattan + secure-vs-plaintext scatter from fed_results.csv (written by fed_compare.py FED_CSV=1).
Reads ONLY aggregate per-gene p-values + public genomic positions — no sample-level data. Plots each
p-value the CSV carries: Burden always; SKAT p when the run computed it (skat_pvalue_probes > 0).

    python3 fed_plot.py [csv] [outdir]
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("fed_results.csv")
outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(csv))

d = np.genfromtxt(csv, delimiter=",", names=True)
chrom = np.atleast_1d(d["chrom"]).astype(int)
pos = np.atleast_1d(d["pos"]).astype(float)
ng = len(chrom)
mlog = lambda p: -np.log10(np.clip(p, 1e-300, 1.0))
order = np.lexsort((pos, chrom))  # genes sorted by (chrom, pos)
thr = 0.05 / ng  # Bonferroni


def manhattan(psec, label, out):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.scatter(np.arange(ng), mlog(psec)[order],
               c=["#1f77b4" if c % 2 == 0 else "#ff7f0e" for c in chrom[order]], s=18)
    ax.axhline(mlog(thr), ls="--", c="red", lw=1, label=f"Bonferroni 0.05/{ng}")
    ax.set_xticks([np.where(chrom[order] == c)[0].mean() for c in np.unique(chrom)])
    ax.set_xticklabels([f"chr{c}" for c in np.unique(chrom)])
    ax.set_ylabel(f"-log10(secure {label} p)")
    ax.set_title(f"Secure federated {label} p-value")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def scatter(psec, pplain, label, out):
    xs, ys = mlog(pplain), mlog(psec)
    ss_tot = np.sum((xs - xs.mean()) ** 2)
    r2 = 1 - np.sum((ys - xs) ** 2) / ss_tot if ss_tot > 0 else float("nan")  # agreement vs y=x
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(xs, ys, s=20, c="#333333")
    hi = max(xs.max(), ys.max()) * 1.05 + 0.1
    ax.plot([0, hi], [0, hi], ls="--", c="red", lw=1, label="y = x")
    ax.text(0.05, 0.95, f"$R^2$ = {r2:.5f}", transform=ax.transAxes, va="top", fontsize=11,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xlabel(f"-log10(plaintext {label} p)")
    ax.set_ylabel(f"-log10(secure {label} p)")
    ax.set_title(f"{label} p-value: secure vs plaintext")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


# One (manhattan, scatter) pair per p-value the CSV carries. Burden is always present; SKAT p appears
# only when the secure run computed it (fed_compare adds the columns when skat_fed_skat_p_out.txt exists).
stats = [("Burden", "burden_p_secure", "burden_p_plain", "burden")]
if "skat_p_secure" in d.dtype.names:
    stats.append(("SKAT", "skat_p_secure", "skat_p_plain", "skat"))

written = []
for label, sec_col, plain_col, tag in stats:
    psec = np.atleast_1d(d[sec_col]).astype(float)
    pplain = np.atleast_1d(d[plain_col]).astype(float)
    manhattan(psec, label, f"{outdir}/manhattan_{tag}.png")
    scatter(psec, pplain, label, f"{outdir}/scatter_{tag}.png")
    written += [f"manhattan_{tag}.png", f"scatter_{tag}.png"]

print(f"wrote {', '.join(written)} in {outdir}  ({ng} genes)")
