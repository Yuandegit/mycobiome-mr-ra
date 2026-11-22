import sys, os, json, pickle
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy import stats
from mr_core import load_bbj, harmonize, ivw, weighted_median, egger, raps, ROOT

cl = pd.read_csv(os.path.join(ROOT, "results/instruments_clumped.csv"),
                 dtype={"pheno_idx": str})
bmap = load_bbj()
HLA = (6, 25_000_000, 35_000_000)


def pval(b, se):
    return 2 * stats.norm.sf(abs(b / se))


rows = []
for thr in ["gw5e8", "expl1e5"]:
    for pid, g in cl[(cl.param == "main") & (cl.threshold == thr)].groupby("pheno_idx"):
        h = harmonize(g, bmap)
        d = h[h.keep].copy()
        d = d[~((d.chrom == HLA[0]) & (d.pos > HLA[1]) & (d.pos < HLA[2]))].reset_index(drop=True)
        if len(d) == 0:
            continue
        if len(d) == 1:
            r = d.iloc[0]; b, se = r.by / r.bx, r.sy / abs(r.bx)
            wm = np.nan
        else:
            b, se, Q, dfree = ivw(d)
            wm = weighted_median(d)[0]
        rows.append(dict(pheno_idx=pid, threshold=thr, n_snp=len(d),
                         b=b, p=pval(b, se), WM_b=wm))

S = pd.DataFrame(rows)
for thr in ["gw5e8", "expl1e5"]:
    ix = S.threshold == thr
    S.loc[ix, "FDR"] = stats.false_discovery_control(S.loc[ix, "p"].values)

S.to_csv(os.path.join(ROOT, "results/sensitivity_noHLA.csv"), index=False)
for thr in ["gw5e8", "expl1e5"]:
    sub = S[S.threshold == thr]
    sig = sub[sub.FDR < 0.05]
    print(f"=== {thr} after removing HLA: {len(sig)} FDR<0.05 of {len(sub)} ===")
    print(sub.sort_values("p").head(8)[["pheno_idx", "n_snp", "b", "p", "FDR", "WM_b"]].to_string(index=False))
