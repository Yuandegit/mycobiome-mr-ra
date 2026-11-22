import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from mr_core import (load_bbj, harmonize, ivw, wald, weighted_median,
                     egger, raps, ROOT)

cl = pd.read_csv(os.path.join(ROOT, "results/instruments_clumped.csv"),
                 dtype={"pheno_idx": str})
bmap = load_bbj()
mp = json.load(open(os.path.join(ROOT, "data/gwas_microbe/fungi_genus/mapping.json")))


def genus(pid):
    s = mp[pid]
    g = s.split("g__")[-1] if "g__" in s else s
    return g.replace("_", " ")


def pval(b, se):
    return 2 * stats.norm.sf(abs(b / se))


results, harmonized_all = [], {}
for thr in ["gw5e8", "expl1e5"]:
    for pid, g in cl[(cl.param == "main") & (cl.threshold == thr)].groupby("pheno_idx"):
        h = harmonize(g, bmap)
        d = h[h.keep].reset_index(drop=True)
        harmonized_all[(pid, thr)] = h
        if len(d) == 0:
            continue
        d["F"] = (d.bx / d.sx)**2
        rec = dict(pheno_idx=pid, genus=genus(pid), threshold=thr,
                   n_snp=len(d), mean_F=d.F.mean(), min_F=d.F.min())
        if len(d) == 1:
            b, se = wald(d); rec.update(method="Wald", b=b, se=se, p=pval(b, se))
        else:
            b, se, Q, dfree = ivw(d)
            rec.update(method="IVW", b=b, se=se, p=pval(b, se),
                       Q=Q, Qp=stats.chi2.sf(Q, dfree))
            wmb, wms = weighted_median(d)
            rec["WM_b"], rec["WM_p"] = wmb, pval(wmb, wms)
            if len(d) >= 3:
                sl, sesl, icpt, Qe = egger(d)
                rec["Egger_b"], rec["Egger_p"] = sl, pval(sl, sesl)
                rec["Egger_int"], rec["Egger_int_p"] = icpt, pval(icpt, 1e3)
                rb, rs = raps(d)
                rec["RAPS_b"], rec["RAPS_p"] = rb, pval(rb, rs)
        results.append(rec)

R = pd.DataFrame(results)
for thr in ["gw5e8", "expl1e5"]:
    sub = R[R.threshold == thr].copy()
    sub["FDR"] = stats.false_discovery_control(sub.p.values)
    R.loc[sub.index, "FDR"] = sub.FDR
R.to_csv(os.path.join(ROOT, "results/forward_mr_results.csv"), index=False)
pickle_h = os.path.join(ROOT, "results/harmonized_forward.pkl")
import pickle
pickle.dump(harmonized_all, open(pickle_h, "wb"))

print("=== STRICT (5e-8) significant IVW/Wald p<0.05 ===")
s = R[(R.threshold == "gw5e8")].sort_values("p")
print(s[["pheno_idx", "genus", "n_snp", "b", "se", "p", "FDR", "mean_F"]].head(15).to_string(index=False))
print("\n=== EXPLORATORY (1e-5) FDR<0.05 count ===")
e = R[R.threshold == "expl1e5"]
print((e.FDR < 0.05).sum(), "of", len(e))
print(e.sort_values("FDR")[["pheno_idx", "genus", "n_snp", "b", "p", "FDR"]].head(12).to_string(index=False))
