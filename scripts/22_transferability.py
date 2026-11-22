"""
M1 response: cross-population transferability check. No Japanese mycobiome GWAS
exists, so a direct exposure-effect replication is impossible. We instead quantify
transferability indirectly through the alignment of instrument allele frequencies
between the Han Chinese exposure (maf) and the Japanese BBJ outcome (eaf, allele
aligned at harmonisation), overall and stratified by variant frequency, including
the low-frequency driver variant. Self-contained Python, no R/PLINK.
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy import stats
from mr_core import ROOT

H = pickle.load(open(os.path.join(ROOT, "results/harmonized_forward.pkl"), "rb"))
frames = []
for key in H:
    if key[1] != "expl1e5":
        continue
    d = H[key][H[key].keep]
    frames.append(d[["rsid", "maf", "eaf"]])
a = pd.concat(frames, ignore_index=True).dropna()
a["abs_diff"] = (a.maf - a.eaf).abs()
lo = lambda p: np.log(np.clip(p, 1e-4, 1 - 1e-4) /
                      (1 - np.clip(p, 1e-4, 1 - 1e-4)))

summary = {
    "n_instruments": len(a),
    "pearson_r": stats.pearsonr(a.maf, a.eaf)[0],
    "spearman_rho": stats.spearmanr(a.maf, a.eaf)[0],
    "pearson_r_logit": stats.pearsonr(lo(a.maf), lo(a.eaf))[0],
    "median_abs_diff": a.abs_diff.median(),
    "p90_abs_diff": a.abs_diff.quantile(0.9),
    "pct_diff_gt_0.1": 100 * (a.abs_diff > 0.1).mean(),
    "min_eaf": a.eaf.min(),
    "n_monomorphic_in_JP": int((a.eaf <= 0).sum()),
}
# low-frequency subset (maf < 0.05), where cross-population divergence is largest
rare = a[a.maf < 0.05]
summary.update({
    "n_rare": len(rare),
    "rare_pearson_r": stats.pearsonr(rare.maf, rare.eaf)[0],
    "rare_median_abs_diff": rare.abs_diff.median(),
})
# driver variant
drv = a[a.rsid == "rs111394117"]
if len(drv):
    r0 = drv.iloc[0]
    summary.update(driver_maf=r0.maf, driver_eaf=r0.eaf,
                   driver_abs_diff=r0.abs_diff)

out = pd.DataFrame([summary])
path = os.path.join(ROOT, "results/transferability.csv")
out.to_csv(path, index=False)
a.to_csv(os.path.join(ROOT, "results/transferability_per_snp.csv"), index=False)
for k, v in summary.items():
    print("%-22s %s" % (k, v))
print("\nwrote", path)
