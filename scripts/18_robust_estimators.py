"""
M3 response: additional pleiotropy-robust estimators and residual-threshold trajectory.
Estimators: weighted mode (MBE), contamination mixture (ConMix), MR-Lasso, weighted TLS (IWTLS).
All implementations are self-contained Python (no R/PLINK).
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy import stats, optimize
from mr_core import ivw, ROOT

H = pickle.load(open(os.path.join(ROOT, "results/harmonized_forward.pkl"), "rb"))


def _arrays(dd):
    dd = dd.reset_index(drop=True)
    return dd.bx.values.astype(float), dd.sx.values.astype(float), \
        dd.by.values.astype(float), dd.sy.values.astype(float)


# ---------------- 1. Weighted mode-based estimator ----------------
def mode_based(dd, n_boot=1000, seed=42):
    bx, sx, by, sy = _arrays(dd)
    theta = by / bx
    w = bx**2 / sy**2
    w = w / w.sum()

    def weighted_mode(th, wt):
        # bandwidth: weighted Silverman rule on the Wald ratios
        med = np.average(th, weights=wt)
        sigma = np.sqrt(np.average((th - med)**2, weights=wt))
        h = 0.9 * min(sigma, stats.iqr(th) / 1.349) * len(th)**(-0.2)
        h = max(h, 1e-3)
        grid = np.linspace(th.min() - 3 * h, th.max() + 3 * h, 2000)
        dens = np.zeros_like(grid)
        for ti, wi in zip(th, wt):
            dens += wi * stats.norm.pdf(grid, ti, h)
        return grid[np.argmax(dens)]

    point = weighted_mode(theta, w)
    rng = np.random.default_rng(seed)
    n = len(theta)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        boots.append(weighted_mode(theta[idx], w[idx] / w[idx].sum()))
    return point, np.std(boots, ddof=1)


# ---------------- 2. Contamination mixture ----------------
def contamination_mixture(dd, n_grid=400):
    bx, sx, by, sy = _arrays(dd)
    theta = by / bx
    se = np.abs(sy / bx)
    lo, hi = theta.min() - 4 * se.max(), theta.max() + 4 * se.max()
    grid = np.linspace(lo, hi, n_grid)
    # uniform component spans the full candidate range
    span = hi - lo

    def neg_loglik(theta0, psi):
        valid = psi * stats.norm.pdf(theta, theta0, se)
        invalid = (1 - psi) / span
        mix = valid + invalid
        return -np.sum(np.log(mix))

    best = None
    for theta0 in grid:
        # profile over mixing proportion psi
        res = optimize.minimize_scalar(
            lambda p: neg_loglik(theta0, np.clip(p, 1e-4, 1 - 1e-4)),
            bounds=(0.05, 0.999), method="bounded")
        if best is None or res.fun < best[1]:
            best = (theta0, res.fun, res.x)
    point, ll, psi = best
    # likelihood-based 95% CI (2*(LL - Lmax) <= 3.84)
    ci = []
    for theta0 in grid:
        l = neg_loglik(theta0, psi)
        if 2 * (l - ll) <= 3.84:
            ci.append(theta0)
    ci = (min(ci), max(ci)) if ci else (np.nan, np.nan)
    return point, (ci[1] - ci[0]) / 2, psi, ci


# ---------------- 3. MR-Lasso ----------------
def mr_lasso(dd):
    """Per-SNP intercepts with an L1 penalty, solved along the LARS-Lasso path.

    The saturated design [bx, I_n] makes cross-validation/BIC model selection
    unstable at small n (it oscillates between 0 and n pleiotropic intercepts).
    We instead use the exact LARS entry order and evaluate the path at the step
    where the causal bx term enters. Variants whose intercepts are already
    non-zero at that point are flagged as pleiotropic; the remaining valid
    instruments are combined by IVW. This is deterministic and reproducible."""
    from sklearn.linear_model import lars_path
    bx, sx, by, sy = _arrays(dd)
    n = len(bx)
    X = np.zeros((n, n + 1))
    X[:, 0] = bx
    for i in range(n):
        X[i, i + 1] = 1.0
    # inverse-variance weighting, then glmnet-style column standardisation
    sw = 1 / sy
    Xw = X * sw[:, None]
    yw = by * sw
    mu = Xw.mean(0)
    sd = Xw.std(0)
    sd[sd == 0] = 1.0
    Xs = (Xw - mu) / sd
    alphas, active, coef = lars_path(Xs, yw, method="lasso")
    # first step at which the causal column (index 0) is active
    step_causal = next((s for s, a in enumerate(active)
                        if 0 in np.atleast_1d(a)), len(active) - 1)
    gamma = coef[1:, step_causal]
    valid = np.abs(gamma) < 1e-6
    d0 = dd.reset_index(drop=True)
    if valid.sum() >= 2:
        b, se, _, _ = ivw(d0[valid])
    elif valid.sum() == 1:
        r = d0[valid].iloc[0]
        b, se = r.by / r.bx, r.sy / abs(r.bx)
    else:
        b, se = np.nan, np.nan
    return b, se, int(valid.sum()), int((~valid).sum())


# ---------------- 4. Weighted TLS (IWTLS, bx and by both measured with error) ----------------
def weighted_tls(dd):
    bx, sx, by, sy = _arrays(dd)

    def objective(theta):
        # orthogonal (total least squares) residual with both SEs
        denom = sy**2 + theta**2 * sx**2
        return np.sum((by - theta * bx)**2 / denom)

    res = optimize.minimize_scalar(objective, method="brent")
    theta = res.x
    # SE via curvature of the profile objective
    grid = np.linspace(theta - abs(theta) * 0.5 - 0.01, theta + abs(theta) * 0.5 + 0.01, 400)
    vals = np.array([objective(t) for t in grid])
    # inverse observed information
    d2 = np.gradient(np.gradient(vals, grid), grid)
    info = d2[np.argmin(abs(grid - theta))]
    se = np.sqrt(1 / info) if info > 0 else np.nan
    return theta, se


# ---------------- 5. Single-pass residual-threshold IVW trajectory ----------------
def threshold_trajectory(dd, thresholds=(2.5, 3.0, 3.5, 4.0)):
    bx0, sx0, by0, sy0 = _arrays(dd)
    b_naive, se_naive, _, _ = ivw(dd)
    r = (by0 - b_naive * bx0) / sy0          # outcome-scaled residual, single pass
    out = []
    for t in thresholds:
        keep = np.abs(r) <= t
        sub = pd.DataFrame({"bx": bx0[keep], "sx": sx0[keep],
                            "by": by0[keep], "sy": sy0[keep]})
        if len(sub) >= 2:
            b, se, _, _ = ivw(sub)
        elif len(sub) == 1:
            b, se = sub.by[0] / sub.bx[0], sub.sy[0] / abs(sub.bx[0])
        else:
            b, se = np.nan, np.nan
        out.append((t, int((~keep).sum()), int(keep.sum()), b, se))
    return b_naive, out


if __name__ == "__main__":
    d = H[("109", "expl1e5")]
    d = d[d.keep]
    print("n instruments:", len(d))
    b0, se0, Q, _ = ivw(d)
    print("IVW: b=%.4f se=%.4f Q=%.2f" % (b0, se0, Q))

    mb, mb_se = mode_based(d)
    print("Mode-based: b=%.4f se=%.4f" % (mb, mb_se))

    cm, cm_se, psi, ci = contamination_mixture(d)
    print("ContamMix: b=%.4f se=%.4f psi=%.3f CI=(%.4f,%.4f)" % (cm, cm_se, psi, ci[0], ci[1]))

    lasso_b, lasso_se, n_valid, n_pleo = mr_lasso(d)
    print("MR-Lasso: b=%.4f se=%.4f validIV=%d pleiotropicIV=%d" % (lasso_b, lasso_se, n_valid, n_pleo))

    tls_b, tls_se = weighted_tls(d)
    print("Weighted TLS: b=%.4f se=%.4f" % (tls_b, tls_se))

    b_naive, traj = threshold_trajectory(d)
    print("\nNaive IVW b=%.4f" % b_naive)
    print("threshold trajectory (thr, removed, kept, IVW b, se):")
    for row in traj:
        print("  %.1f  rem=%2d keep=%2d  b=%.4f se=%.4f" % row)
