#!/usr/bin/env python3
"""Fetch EAS (1000G, n=504) genotypes for candidate variants via remote tabix.

Memory-safe: process one chromosome at a time with a fresh VariantFile, write a
per-chromosome pickle, then release everything. Resumable.

Usage: python3 04_fetch_eas_genotypes.py [chrN ...]
"""
import json, glob, os, sys, time, pickle, gc
import numpy as np
import pandas as pd
import pysam

ROOT = "/home/user/mycobiome_joint_mr"
GDIR = os.path.join(ROOT, "data/gwas_microbe/fungi_genus")
OUTDIR = os.path.join(ROOT, "data/gwas_microbe/eas_chr")
os.makedirs(OUTDIR, exist_ok=True)
VCF = ("https://huggingface.co/datasets/Cloufield/GWASLab-Reference/resolve/main/"
       "1kg/hg19/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz")
GAP, PTHR = 100_000, 1e-5

# ---- collect unique candidate variants ----
recs = {}
for p in sorted(glob.glob(GDIR + "/[0-9]*.json")):
    df = pd.DataFrame(json.load(open(p))["unbinned_variants"])
    df = df[df.pval < PTHR]
    for _, r in df.iterrows():
        recs[(int(r.chrom), int(r.pos), r.ref, r.alt)] = 1
all_chroms = sorted(set(k[0] for k in recs))
arg = [a.replace("chr", "") for a in sys.argv[1:]]
if arg:
    all_chroms = [c for c in all_chroms if str(c) in arg]


def build_regions(pos):
    pos = sorted(pos); regions = []; lo = prev = pos[0]
    for x in pos[1:]:
        if x - prev > GAP:
            regions.append((lo, prev)); lo = x
        prev = x
    regions.append((lo, prev)); return regions


for chrom in all_chroms:
    outp = os.path.join(OUTDIR, f"chr{chrom}.pkl")
    if os.path.exists(outp):
        print("skip chr%d (exists)" % chrom); continue
    keys = [k for k in recs if k[0] == chrom]
    want = {}
    for (c, po, rf, a) in keys:
        want.setdefault(po, set()).add((rf, a))
    regions = build_regions(list(want))
    vcf = pysam.VariantFile(VCF)
    n_samp = len(vcf.header.samples)
    geno = {}
    t0 = time.time()

    def process(lo, hi):
        n = 0
        for rec in vcf.fetch(str(chrom), max(0, lo - 1), hi + 1):
            cands = want.get(rec.pos)
            if not cands:
                continue
            g = np.full(n_samp, -1, dtype=np.int8)
            for j, s in enumerate(rec.samples.values()):
                a, b = s["GT"]
                if a is not None and b is not None:
                    g[j] = a + b
            vr, va = rec.ref, rec.alts[0] if rec.alts else None
            for (rf, a) in cands:
                k = (chrom, rec.pos, rf, a)
                if vr == rf and va == a:
                    geno[k] = g
                elif vr == a and va == rf:
                    gf = g.copy(); gf[g >= 0] = 2 - g[g >= 0]; geno[k] = gf
            n += 1
        return n

    for lo, hi in regions:
        for attempt in range(6):
            try:
                process(lo, hi); break
            except Exception as e:
                try:
                    vcf.close()
                except Exception:
                    pass
                del vcf; gc.collect(); time.sleep(3 + attempt * 2)
                vcf = pysam.VariantFile(VCF)
                if attempt == 5:
                    raise
    with open(outp, "wb") as f:
        pickle.dump({"geno": geno, "n_samp": n_samp}, f)
    print("chr%d: %d/%d variants in %.0fs" % (chrom, len(geno), len(keys), time.time() - t0))
    del vcf, geno, want, regions
    gc.collect()

print("done. per-chr files in", OUTDIR)
