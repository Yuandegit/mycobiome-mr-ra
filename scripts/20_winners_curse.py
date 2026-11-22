"""
M2 response: quantify first-stage winner's curse and report calibrated instrument
strength. For each of the 115 features we report total R^2 under three calibrations:
  (1) raw get_r_from_pn                 -- optimistic upper bound (in-sample selected)
  (2) empirical-Bayes linear shrinkage  -- realistic global calibration
  (3) two-sided truncated conditional MLE -- conservative per-variant lower bound
plus joint (multi-instrument) F statistics. Self-contained Python, no R/PLINK.
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy import stats, optimize
from mr_core import ROOT, r2_from_pn

N_X = 7350
C = stats.norm.ppf(1 - 1e-5 / 2)
H = pickle.load(open(os.path.join(ROOT, "results/harmonized_forward.pkl"), "rb"))


def cond_mle(z, c=C):
    """Signed conditional MLE of the non-centrality under |Z|>c selection."""
    s = np.sign(z)
    z = abs(z)

    def nll(mu):
        sel = 1 - stats.norm.cdf(c - mu) + stats.norm.cdf(-c - mu)
        return -(stats.norm.logpdf(z - mu) - np.log(sel))

    r = optimize.minimize_scalar(nll, bounds=(0, z), method="bounded")
    return s * r.x


# ---- global EB prior variance from all selected instruments (truncation-aware ML)
all_z = []
for key in H:
    if key[1] != "expl1e5":
        continue
    d = H[key]
    all_z.extend((d[d.keep].bx / d[d.keep].sx).tolist())
all_z = np.array(all_z)


def marg_nll(logt):
    t = np.exp(logt)
    sd = np.sqrt(t + 1)
    sel = 2 * stats.norm.sf(C / sd)
    return -np.sum(stats.norm.logpdf(all_z, 0, sd) - np.log(sel))


TAU2 = np.exp(optimize.minimize_scalar(
    marg_nll, bounds=(-4, 3), method="bounded").x)
SHRINK = TAU2 / (TAU2 + 1)


def joint_f(r2_total, k, n=N_X):
    if r2_total <= 0 or r2_total >= 1:
        return np.nan
    return (r2_total / (1 - r2_total)) * ((n - 1 - k) / k)


rows = []
for key in sorted(H, key=lambda k: int(k[0])):
    pheno, thr = key
    if thr != "expl1e5":
        continue
    d = H[key][H[key].keep].reset_index(drop=True)
    k = len(d)
    z = (d.bx / d.sx).values
    r2_raw = r2_from_pn(d.bx, d.sx, N_X)
    # EB: shrink the z, convert the calibrated z back to R2
    z_eb = z * SHRINK
    r2_eb = r2_from_pn(z_eb, np.ones(k), N_X)
    # conditional MLE per variant
    z_cm = np.array([cond_mle(zz) for zz in z])
    r2_cm = r2_from_pn(z_cm, np.ones(k), N_X)
    rows.append({
        "pheno": int(pheno), "n_iv": k,
        "r2_snp_median_raw": np.median(r2_raw),
        "r2_total_raw": r2_raw.sum(),
        "r2_total_eb": r2_eb.sum(),
        "r2_total_condmle": r2_cm.sum(),
        "F_joint_raw": joint_f(r2_raw.sum(), k),
        "F_joint_eb": joint_f(r2_eb.sum(), k),
        "F_joint_condmle": joint_f(r2_cm.sum(), k),
    })

out = pd.DataFrame(rows)
path = os.path.join(ROOT, "results/winners_curse_calibration.csv")
out.to_csv(path, index=False)
print("EB tau2=%.4f shrink=%.4f" % (TAU2, SHRINK))
print("wrote", path, "rows:", len(out))
for c in ["r2_total_raw", "r2_total_eb", "r2_total_condmle"]:
    print("%-16s median %.4f  max %.4f" % (c, out[c].median(), out[c].max()))
for c in ["F_joint_raw", "F_joint_eb", "F_joint_condmle"]:
    print("%-16s median %.2f  min %.2f" % (c, out[c].median(), out[c].min()))
