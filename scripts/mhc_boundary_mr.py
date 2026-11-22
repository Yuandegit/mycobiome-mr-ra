import pandas as pd, numpy as np
from scipy import stats

R = "/home/user/mycobiome_joint_mr"
h = pd.read_pickle(R + "/results/harmonized_forward.pkl")

def ivw_one(d):
    d = d[d.keep == True]
    n = len(d)
    if n == 0:
        return np.nan, np.nan, 0
    bx, by, sy = d.bx.values, d.by.values, d.sy.values
    w = bx**2 / sy**2
    if n == 1:
        b = by[0] / bx[0]
        se = sy[0] / abs(bx[0])
    else:
        b = np.sum(w * (by / bx)) / np.sum(w)
        se = np.sqrt(1.0 / np.sum(w))
    p = 2 * stats.norm.sf(abs(b / se))
    return b, p, n

def bh(pvals):
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    fdr = np.empty(m)
    prev = 1.0
    for rank in range(m, 0, -1):
        idx = order[rank - 1]
        val = min(prev, p[idx] * m / rank)
        fdr[idx] = val
        prev = val
    return fdr

scenarios = {
    "full":      None,
    "wide_25_35": (25e6, 35e6),
    "core_28.5_33.5": (28.5e6, 33.5e6),
}
phenos = sorted(set(k[0] for k in h if k[1] == "expl1e5"))
out = pd.DataFrame(index=phenos)

for sname, win in scenarios.items():
    bs, ps, ns = [], [], []
    for ph in phenos:
        d = h[(ph, "expl1e5")].copy()
        if win is not None:
            lo, hi = win
            drop = (d.chrom == 6) & (d.pos >= lo) & (d.pos <= hi)
            d = d[~drop]
        b, p, n = ivw_one(d)
        bs.append(b); ps.append(p); ns.append(n)
    out[sname + "_n"] = ns
    out[sname + "_b"] = bs
    out[sname + "_p"] = ps
    out[sname + "_FDR"] = bh(np.nan_to_num(ps, nan=1.0))

pd.set_option("display.width", 200)
for sname in scenarios:
    f = out[sname + "_FDR"]
    print(f"=== {sname}: FDR<0.05 count={(f<0.05).sum()}, min FDR={f.min():.4f}")

print("\n--- trait 109 (Cutaneotrichosporon) across scenarios ---")
for sname in scenarios:
    print(sname, "n=", out.loc["109", sname+"_n"],
          "b=", round(out.loc["109", sname+"_b"], 4),
          "p=", f"{out.loc['109', sname+'_p']:.3g}",
          "FDR=", f"{out.loc['109', sname+'_FDR']:.3g}")

out.to_csv(R + "/results/mhc_boundary_sensitivity.csv")
print("\nsaved results/mhc_boundary_sensitivity.csv")
