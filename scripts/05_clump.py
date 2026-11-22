#!/usr/bin/env python3
"""Greedy LD clumping (PLINK-style) using cached EAS genotypes.

For each of 115 genus phenotypes and each (r2, kb) parameter set, produce the
set of independent instruments separately at the strict (5e-8) and exploratory
(1e-5) significance thresholds. r2 is orientation-independent, so cached
dosages can be used directly.
"""
import json, glob, os, pickle
import numpy as np
import pandas as pd

ROOT = "/home/user/mycobiome_joint_mr"
GDIR = os.path.join(ROOT, "data/gwas_microbe/fungi_genus")
GENO, N_SAMP = {}, None
for cp in sorted(glob.glob(os.path.join(ROOT, "data/gwas_microbe/eas_chr/chr*.pkl"))):
    c = pickle.load(open(cp, "rb")); GENO.update(c["geno"]); N_SAMP = c["n_samp"]

PARAMS = [("main", 0.1, 500), ("sens_strict", 0.01, 1000),
          ("sens_vstrict", 0.001, 1000), ("sens_loose", 0.1, 250)]
THRESH = {"gw5e8": 5e-8, "expl1e5": 1e-5}


def pairwise_r2(keys):
    X = np.stack([GENO[k].astype(float) for k in keys], axis=0)   # (m, n)
    X = np.where(X < 0, np.nan, X)
    Xc = X - np.nanmean(X, axis=1, keepdims=True)
    Xc = np.where(np.isnan(Xc), 0.0, Xc)       # missing -> 0 (no contribution)
    num = Xc @ Xc.T
    d = np.sqrt((Xc ** 2).sum(axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        R = num / (d[:, None] * d[None, :])
    return np.clip(R, -1, 1) ** 2


def greedy(df, r2thr, kb):
    df = df.sort_values("pval").reset_index(drop=True)
    if len(df) == 0:
        return df
    keys = list(df["key"])
    avail = [k in GENO for k in keys]
    df = df[avail].reset_index(drop=True)
    keys = list(df["key"])
    if len(keys) <= 1:
        return df
    R2 = pairwise_r2(keys)
    keep = np.ones(len(keys), bool)
    for i in range(len(keys)):
        if not keep[i]:
            continue
        for j in range(i + 1, len(keys)):
            if keep[j] and df.chrom[i] == df.chrom[j] and abs(df.pos[i] - df.pos[j]) <= kb * 1000:
                if R2[i, j] >= r2thr:
                    keep[j] = False
    return df[keep].reset_index(drop=True)


rows = []
detail = {}
for p in sorted(glob.glob(GDIR + "/[0-9]*.json")):
    pid = os.path.basename(p).replace(".json", "")
    df = pd.DataFrame(json.load(open(p))["unbinned_variants"])
    df["chrom"] = df.chrom.astype(int); df["pos"] = df.pos.astype(int)
    df["key"] = list(zip(df.chrom, df.pos, df.ref, df.alt))
    for tname, thr in THRESH.items():
        sub = df[df.pval < thr]
        for pname, r2thr, kb in PARAMS:
            cl = greedy(sub, r2thr, kb)
            detail[(pid, tname, pname)] = cl
            for _, r in cl.iterrows():
                rows.append({"pheno_idx": pid, "threshold": tname, "param": pname,
                             "chrom": r.chrom, "pos": r.pos, "ref": r.ref, "alt": r.alt,
                             "rsids": r.rsids, "pval": r.pval, "beta": r.beta,
                             "sebeta": r.sebeta, "maf": r.maf})

out = pd.DataFrame(rows)
out.to_csv(os.path.join(ROOT, "results/instruments_clumped.csv"), index=False)
# summary counts
summ = out.groupby(["threshold", "param"]).size().reset_index(name="total_snps")
nph = out.groupby(["threshold", "param"]).pheno_idx.nunique().reset_index(name="phenos_with_iv")
print(pd.merge(summ, nph).to_string(index=False))
pickle.dump(detail, open(os.path.join(ROOT, "results/instruments_detail.pkl"), "wb"))
print("saved results/instruments_clumped.csv and instruments_detail.pkl")
