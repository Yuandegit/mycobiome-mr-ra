"""
M3 (batch): run the four pleiotropy-robust estimators and the residual-threshold
trajectory for every one of the 115 genus-level features (exploratory P<1e-5).
Outputs results/robust_estimators_all.csv. Self-contained Python, no R/PLINK.
"""
import sys, os, importlib.util, pickle
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy import stats
from mr_core import ROOT

# load the 18_ module (leading digit prevents a normal import)
_spec = importlib.util.spec_from_file_location(
    "robust18", os.path.join(os.path.dirname(__file__), "18_robust_estimators.py"))
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)

H = pickle.load(open(os.path.join(ROOT, "results/harmonized_forward.pkl"), "rb"))
names = pd.read_csv(os.path.join(ROOT, "results/forward_mr_results.csv"))[
    ["pheno_idx", "genus"]].drop_duplicates("pheno_idx")
name_map = dict(zip(names.pheno_idx.astype(str), names.genus))


def pval(b, se):
    if b is None or np.isnan(b) or np.isnan(se) or se <= 0:
        return np.nan
    return 2 * stats.norm.sf(abs(b / se))


def safe(fn, dd):
    try:
        return fn(dd)
    except Exception:
        return None


rows = []
for key in sorted(H, key=lambda k: int(k[0])):
    pheno, thr = key
    if thr != "expl1e5":
        continue
    d = H[key]
    d = d[d.keep].reset_index(drop=True)
    n = len(d)
    pidx = int(pheno)
    row = {"pheno": pidx, "genus": name_map.get(str(pidx), ""), "n_iv": n}

    mb = safe(R.mode_based, d)
    if mb:
        row.update(mode_b=mb[0], mode_se=mb[1], mode_p=pval(mb[0], mb[1]))

    cm = safe(R.contamination_mixture, d)
    if cm:
        row.update(conmix_b=cm[0], conmix_se=cm[1], conmix_p=pval(cm[0], cm[1]),
                   conmix_psi=cm[2], conmix_ci_lo=cm[3][0], conmix_ci_hi=cm[3][1])

    la = safe(R.mr_lasso, d)
    if la:
        row.update(lasso_b=la[0], lasso_se=la[1], lasso_p=pval(la[0], la[1]),
                   lasso_nvalid=la[2], lasso_npleo=la[3])

    tls = safe(R.weighted_tls, d)
    if tls:
        row.update(tls_b=tls[0], tls_se=tls[1], tls_p=pval(tls[0], tls[1]))

    tt = safe(R.threshold_trajectory, d)
    if tt:
        for t, removed, kept, b, se in tt[1]:
            tag = "thr%s" % str(t).replace(".", "_")
            row[tag + "_b"] = b
            row[tag + "_se"] = se
            row[tag + "_p"] = pval(b, se)
            row[tag + "_removed"] = removed
    rows.append(row)

out = pd.DataFrame(rows)
path = os.path.join(ROOT, "results/robust_estimators_all.csv")
out.to_csv(path, index=False)
print("wrote", path, "rows:", len(out), "cols:", len(out.columns))
