#!/usr/bin/env python3
"""Post-clumping verification: strict pairwise-complete r2 among the FINAL
retained instrument sets. Flags any pair (same chrom, within window) whose
pairwise-complete r2 reaches the parameter threshold.
"""
import glob, os, pickle
import numpy as np
import pandas as pd

ROOT = "/home/user/mycobiome_joint_mr"
GENO = {}
for cp in sorted(glob.glob(os.path.join(ROOT, "data/gwas_microbe/eas_chr/chr*.pkl"))):
    GENO.update(pickle.load(open(cp, "rb"))["geno"])
detail = pickle.load(open(os.path.join(ROOT, "results/instruments_detail.pkl"), "rb"))

PARAMS = {"main": (0.1, 500), "sens_strict": (0.01, 1000),
          "sens_vstrict": (0.001, 1000), "sens_loose": (0.1, 250)}


def pair_r2(k1, k2):
    a = GENO[k1].astype(float); b = GENO[k2].astype(float)
    m = (a >= 0) & (b >= 0)
    if m.sum() < 30:
        return np.nan
    x, y = a[m], b[m]
    r = np.corrcoef(x - x.mean(), y - y.mean())[0, 1]
    return r * r


violations = []
max_seen = 0.0
checked_pairs = 0
for (pid, tname, pname), df in detail.items():
    r2thr, kb = PARAMS[pname]
    df = df.reset_index(drop=True)
    keys = list(df["key"])
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if df.chrom[i] != df.chrom[j] or abs(df.pos[i] - df.pos[j]) > kb * 1000:
                continue
            r2 = pair_r2(keys[i], keys[j])
            checked_pairs += 1
            if not np.isnan(r2):
                max_seen = max(max_seen, r2)
                if r2 >= r2thr:
                    violations.append((pid, tname, pname, df.rsids[i], df.rsids[j], round(r2, 3)))

print("instrument sets checked:", len(detail))
print("within-window pairs evaluated:", checked_pairs)
print("maximum pairwise-complete r2 among retained instruments:", round(max_seen, 4))
print("violations (r2 >= threshold):", len(violations))
for v in violations[:20]:
    print("  ", v)
