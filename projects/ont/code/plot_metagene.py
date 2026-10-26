#!/usr/bin/env python3
"""Stitch the three computeMatrix outputs (5'UTR, CDS, 3'UTR) into one annotated metagene plot.

Layout: upstream flank | TSS | 5'UTR | start codon | CDS (gene body) | stop codon | 3'UTR | TES | downstream flank
Plots the mean across transcripts (NaN = no signal/coverage ignored) with a +/- SEM band.
"""
import argparse
import gzip
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    with gzip.open(path, "rt") as f:
        hdr = json.loads(f.readline().lstrip("@"))
        rows = [l.rstrip("\n").split("\t")[6:] for l in f]
    first = lambda v: v[0] if isinstance(v, list) else v
    mat = np.array(rows, dtype=float)
    bs = first(hdr["bin size"])
    return mat, first(hdr["upstream"]) // bs, first(hdr["downstream"]) // bs, bs


def smooth(y, k):
    if k <= 1:
        return y
    pad = k // 2
    yp = np.pad(y, pad, mode="edge")
    return np.convolve(yp, np.ones(k) / k, mode="valid")[: len(y)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--utr5", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--utr3", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--title", default="m6A signal across a metagene")
    ap.add_argument("--ylabel", default="mean m6A signal")
    ap.add_argument("--smooth", type=int, default=3, help="moving-average window in bins (1 = none)")
    a = ap.parse_args()

    m5, up5, _, bs = load(a.utr5)
    mc, _, _, _ = load(a.cds)
    m3, _, dn3, _ = load(a.utr3)
    n5, nc, n3 = m5.shape[1] - up5, mc.shape[1], m3.shape[1] - dn3

    full = np.hstack([m5, mc, m3])
    with np.errstate(all="ignore"):
        mean = np.nanmean(full, axis=0)
        nvalid = np.sum(~np.isnan(full), axis=0)
        sem = np.nanstd(full, axis=0, ddof=1) / np.sqrt(np.maximum(nvalid, 1))
    mean, sem = smooth(np.nan_to_num(mean), a.smooth), smooth(np.nan_to_num(sem), a.smooth)

    # boundaries in bin units
    b = np.cumsum([0, up5, n5, nc, n3, dn3])
    x = np.arange(full.shape[1])
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.fill_between(x, mean - sem, mean + sem, color="#4C78A8", alpha=0.25, lw=0)
    ax.plot(x, mean, color="#4C78A8", lw=1.8)
    shade = {1: "#F2F2F2", 3: "#F2F2F2"}  # 5'UTR and 3'UTR shaded; CDS and flanks white
    for i, c in shade.items():
        ax.axvspan(b[i], b[i + 1], color=c, zorder=0)
    for v, lab in [(b[1], "TSS"), (b[2], "start\ncodon"), (b[3], "stop\ncodon"), (b[4], "TES")]:
        ax.axvline(v, color="#888", lw=0.8, ls="--")
    ymax = ax.get_ylim()[1]
    for i, lab in enumerate(["upstream", "5′UTR", "CDS\n(gene body)", "3′UTR", "downstream"]):
        ax.text((b[i] + b[i + 1]) / 2, ymax, lab, ha="center", va="bottom", fontsize=10)
    ax.set_xticks([b[1], b[2], b[3], b[4]])
    ax.set_xticklabels(["TSS", "start codon", "stop codon", "TES"])
    ax.set_xlim(0, full.shape[1] - 1)
    ax.set_ylabel(a.ylabel)
    ax.set_title(f"{a.title} (n = {full.shape[0]:,} transcripts)", pad=34)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out_prefix}.{ext}", dpi=300)

    with open(f"{a.out_prefix}_profile.tsv", "w") as f:
        f.write("bin\tregion\tmean\tsem\tn_transcripts_with_signal\n")
        names = ["upstream"] * up5 + ["5UTR"] * n5 + ["CDS"] * nc + ["3UTR"] * n3 + ["downstream"] * dn3
        for i in range(full.shape[1]):
            f.write(f"{i}\t{names[i]}\t{mean[i]:.6g}\t{sem[i]:.6g}\t{nvalid[i]}\n")
    print(f"bin size {bs} nt; bins: up {up5}, 5UTR {n5}, CDS {nc}, 3UTR {n3}, down {dn3}")


if __name__ == "__main__":
    main()
