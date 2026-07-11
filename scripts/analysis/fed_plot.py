#!/usr/bin/env python3
"""
    run: python3 fed_plot.py [csv] [outdir] 
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
ps = np.atleast_1d(d["burden_p_secure"]).astype(float)
pp = np.atleast_1d(d["burden_p_plain"]).astype(float)
pos = np.atleast_1d(d["pos"]).astype(float)
ng = len(ps)
mlog = lambda p: -np.log10(np.clip(p, 1e-300, 1.0))
thr = 0.05 / ng  # Bonferroni

# --- Manhattan (secure Burden p) --- x = genes sorted by (chrom, pos); color alternates by chrom.
order = np.lexsort((pos, chrom))
fig, ax = plt.subplots(figsize=(10, 4))
ax.scatter(np.arange(ng), mlog(ps)[order],
           c=["#1f77b4" if c % 2 == 0 else "#ff7f0e" for c in chrom[order]], s=18)
ax.axhline(mlog(thr), ls="--", c="red", lw=1, label=f"Bonferroni 0.05/{ng}")
ax.set_xticks([np.where(chrom[order] == c)[0].mean() for c in np.unique(chrom)])
ax.set_xticklabels([f"chr{c}" for c in np.unique(chrom)])
ax.set_ylabel("-log10(secure Burden p)")
ax.set_title("Secure federated Burden p-value")
ax.legend()
fig.tight_layout()
fig.savefig(f"{outdir}/manhattan.png", dpi=150)
plt.close(fig)

# --- Scatter: secure vs plaintext (y=x ⇒ perfect agreement) ---
xs, ys = mlog(pp), mlog(ps)
ss_tot = np.sum((xs - xs.mean()) ** 2)
r2 = 1 - np.sum((ys - xs) ** 2) / ss_tot if ss_tot > 0 else float("nan")  # agreement vs y=x, on -log10 p
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(xs, ys, s=20, c="#333333")
hi = max(xs.max(), ys.max()) * 1.05 + 0.1
ax.plot([0, hi], [0, hi], ls="--", c="red", lw=1, label="y = x")
ax.text(0.05, 0.95, f"$R^2$ = {r2:.5f}", transform=ax.transAxes, va="top", fontsize=11,
        bbox=dict(boxstyle="round", fc="white", ec="0.7"))
ax.set_xlim(0, hi)
ax.set_ylim(0, hi)
ax.set_xlabel("-log10(plaintext Burden p)")
ax.set_ylabel("-log10(secure Burden p)")
ax.set_title("Burden p-value")
ax.legend()
fig.tight_layout()
fig.savefig(f"{outdir}/scatter.png", dpi=150)
plt.close(fig)

print(f"wrote {outdir}/manhattan.png, {outdir}/scatter.png  ({ng} genes)")
