#!/usr/bin/env python3
"""Build strand-aware, spliced BED12 files for 5'UTR, CDS and 3'UTR from a GENCODE GTF.

One transcript per gene (Ensembl_canonical protein-coding by default). Only
transcripts with a complete CDS and non-empty 5'UTR, CDS and 3'UTR are kept, so
all three BED files describe the same set of transcripts.
Standard library only.
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

ATTR = re.compile(r'(\w+) "([^"]*)"')


def parse_attrs(s):
    d, tags = {}, []
    for k, v in ATTR.findall(s):
        if k == "tag":
            tags.append(v)
        else:
            d[k] = v
    d["tags"] = tags
    return d


def clip(exons, lo, hi):
    out = []
    for s, e in exons:
        s2, e2 = max(s, lo), min(e, hi)
        if s2 <= e2:
            out.append((s2, e2))
    return out


def length(blocks):
    return sum(e - s + 1 for s, e in blocks)


def bed12(chrom, blocks, strand, name):
    start, end = blocks[0][0] - 1, blocks[-1][1]
    sizes = ",".join(str(e - s + 1) for s, e in blocks) + ","
    starts = ",".join(str(s - 1 - start) for s, e in blocks) + ","
    return "\t".join(map(str, [chrom, start, end, name, 0, strand, start, end, 0, len(blocks), sizes, starts]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--all-protein-coding", action="store_true",
                    help="use every complete protein-coding transcript, not just one canonical per gene")
    ap.add_argument("--min-utr", type=int, default=20, help="minimum spliced UTR length (nt)")
    ap.add_argument("--min-cds", type=int, default=60, help="minimum spliced CDS length (nt)")
    a = ap.parse_args()

    opener = gzip.open if a.gtf.endswith(".gz") else open
    tx = {}
    with opener(a.gtf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("exon", "CDS", "stop_codon"):
                continue
            at = parse_attrs(c[8])
            if at.get("transcript_type") != "protein_coding":
                continue
            if not a.all_protein_coding and "Ensembl_canonical" not in at["tags"]:
                continue
            t = tx.setdefault(at["transcript_id"], dict(
                chrom=c[0], strand=c[6], gene_id=at.get("gene_id", ""), gene_name=at.get("gene_name", ""),
                tags=set(at["tags"]), exons=[], cds=[]))
            span = (int(c[3]), int(c[4]))
            if c[2] == "exon":
                t["exons"].append(span)
            else:  # CDS and stop codon together define the full coding extent
                t["cds"].append(span)

    # one transcript per gene (longest spliced length) unless --all-protein-coding
    best = {}
    for tid, t in tx.items():
        if not t["cds"] or not t["exons"]:
            continue
        if t["tags"] & {"cds_start_NF", "cds_end_NF"}:
            continue
        t["exons"].sort()
        t["len"] = length(t["exons"])
        if a.all_protein_coding:
            best[tid] = t
        elif t["gene_id"] not in best or t["len"] > best[t["gene_id"]][1]["len"]:
            best[t["gene_id"]] = (tid, t)
    items = list(best.items()) if a.all_protein_coding else [v for v in best.values()]

    recs = []
    for tid, t in items:
        ex = t["exons"]
        cmin = min(s for s, e in t["cds"])
        cmax = max(e for s, e in t["cds"])
        left = clip(ex, -1, cmin - 1)
        right = clip(ex, cmax + 1, 10**12)
        cds = clip(ex, cmin, cmax)
        u5, u3 = (left, right) if t["strand"] == "+" else (right, left)
        if length(u5) < a.min_utr or length(u3) < a.min_utr or length(cds) < a.min_cds:
            continue
        recs.append((t["chrom"], ex[0][0], tid, t, u5, cds, u3))
    recs.sort(key=lambda r: (r[0], r[1]))

    import os
    os.makedirs(a.outdir, exist_ok=True)
    outs = {k: open(f"{a.outdir}/{k}.bed", "w") for k in ("utr5", "cds", "utr3")}
    with open(f"{a.outdir}/transcripts.tsv", "w") as info:
        info.write("transcript_id\tgene_id\tgene_name\tchrom\tstrand\tlen_utr5\tlen_cds\tlen_utr3\n")
        for chrom, _, tid, t, u5, cds, u3 in recs:
            name = tid
            outs["utr5"].write(bed12(chrom, u5, t["strand"], name) + "\n")
            outs["cds"].write(bed12(chrom, cds, t["strand"], name) + "\n")
            outs["utr3"].write(bed12(chrom, u3, t["strand"], name) + "\n")
            info.write("\t".join(map(str, [tid, t["gene_id"], t["gene_name"], chrom, t["strand"],
                                           length(u5), length(cds), length(u3)])) + "\n")
    for o in outs.values():
        o.close()
    print(f"kept {len(recs)} transcripts -> {a.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
