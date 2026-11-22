import sys, os, pickle
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from mr_core import ivw, ROOT
from scipy import stats

H = pickle.load(open(os.path.join(ROOT, "results/harmonized_forward.pkl"), "rb"))
N_SIM = 1000
rng = np.random.default_rng(1)


def global_test(dd, n_sim=N_SIM):
    """MR-PRESSO-style global test: RSS null distribution by simulation."""
    dd = dd.reset_index(drop=True)
    bx, by, sy = dd.bx.values, dd.by.values, dd.sy.values
    b0 = ivw(dd)[0]
    obs_sse = np.sum(((by - b0 * bx) / sy)**2)
    sim_sse = np.array([np.sum((rng.normal(0, sy) / sy)**2) for _ in range(n_sim)])
    return (np.sum(sim_sse >= obs_sse) + 1) / (n_sim + 1)


def studentized_residuals(dd):
    """Leave-one-out studentized residuals (diagnostic, not PRESSO outlier test)."""
    dd = dd.reset_index(drop=True)
    out = []
    for i in range(len(dd)):
        bi = ivw(dd.drop(i))[0]
        ri = (dd.by[i] - bi * dd.bx[i]) / dd.sy[i]
        p = 2 * stats.norm.sf(abs(ri))
        out.append((dd.rsid[i], dd.chrom[i], dd.pos[i], round(ri, 2), round(p, 4)))
    return out


d109 = H[("109", "expl1e5")]; d109 = d109[d109.keep]
gp = global_test(d109)
res = studentized_residuals(d109)
print("=== Cutaneotrichosporon (109) ===")
print("global pleiotropy test p:", gp)
strong = [r for r in res if r[4] < 0.05]
print("studentized-residual outliers (p<0.05):", strong)
print("driver rs111394117 flagged?", any(r[2] == 32187213 for r in res))

glob_sig = total = 0
for (pid, thr), hh in H.items():
    if thr != "expl1e5" or hh.keep.sum() < 4:
        continue
    total += 1
    if global_test(hh[hh.keep]) < 0.05:
        glob_sig += 1
print(f"\nphenotypes with global pleiotropy p<0.05: {glob_sig} of {total}")
