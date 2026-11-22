import json, glob, os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/user/mycobiome_joint_mr"
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def load_mapping():
    mp = json.load(open(os.path.join(ROOT, "data/gwas_microbe/fungi_genus/mapping.json")))
    if isinstance(mp, dict) and "phenotypes" in mp:
        mp = mp["phenotypes"]
    out = {}
    for e in mp:
        if isinstance(e, dict):
            out[str(e.get("phenocode", e.get("id"))).zfill(3)] = e
    return out


def load_bbj():
    b = pd.read_csv(os.path.join(ROOT, "data/gwas_disease/bbj_ra_outcome.tsv"), sep="\t")
    return {(r.CHR, r.POS): r for r in b.itertuples()}


def harmonize(inst, bmap, pal_lo=0.42, pal_hi=0.58):
    """Return per-SNP harmonized rows: bx,sx,by,sy,maf,eaf,keep,reason."""
    rows = []
    for r in inst.itertuples():
        b = bmap.get((r.chrom, r.pos))
        rec = dict(rsid=r.rsids, chrom=r.chrom, pos=r.pos, bx=r.beta, sx=r.sebeta,
                   maf=r.maf, keep=False, reason="")
        if b is None:
            rec["reason"] = "no_outcome"
            rows.append(rec); continue
        by, sy = b.BETA, b.SE
        eaf = b.AF_Allele2
        pal = (r.ref, r.alt) in PAL
        if r.alt == b.Allele2 and r.ref == b.Allele1:
            pass
        elif r.alt == b.Allele1 and r.ref == b.Allele2:
            by, eaf = -by, 1 - eaf
        else:
            rec.update(reason="allele_mismatch", by=by, sy=sy, eaf=eaf); rows.append(rec); continue
        if pal and pal_lo < r.maf < pal_hi and pal_lo < eaf < pal_hi:
            rec.update(reason="ambiguous_palindrome", by=by, sy=sy, eaf=eaf); rows.append(rec); continue
        rec.update(by=by, sy=sy, eaf=eaf, keep=True, reason="ok")
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------- estimators ----------------
def ivw(dd):
    bx, by, sy = dd.bx.values, dd.by.values, dd.sy.values
    w = 1 / sy**2
    b = np.sum(w * bx * by) / np.sum(w * bx**2)
    se = np.sqrt(1 / np.sum(w * bx**2))
    # Cochran Q
    Q = np.sum(w * (by - b * bx)**2)
    dfree = len(dd) - 1
    return b, se, Q, dfree


def wald(dd):
    r = dd.iloc[0]
    b = r.by / r.bx
    se = r.sy / abs(r.bx)
    return b, se


def weighted_median(dd):
    bx, by, sy = dd.bx.values, dd.by.values, dd.sy.values
    wb = by / bx
    w = bx**2 / sy**2  # inverse variance of each Wald ratio
    order = np.argsort(wb)
    wb_s, w_s = wb[order], w[order]
    cw = np.cumsum(w_s) / np.sum(w_s)
    i = np.searchsorted(cw, 0.5)
    point = wb_s[i]

    def _wm(bx_, by_, sy_):
        wb_ = by_ / bx_
        w_ = bx_**2 / sy_**2
        o = np.argsort(wb_)
        c = np.cumsum(w_[o]) / np.sum(w_)
        return wb_[o][np.searchsorted(c, 0.5)]

    # Bootstrap SE (TwoSampleMR default); robust when one IV dominates the weights.
    rng = np.random.default_rng(42)
    n = len(bx)
    boots = []
    for _ in range(1000):
        idx = rng.integers(0, n, n)  # one shared index resamples whole observations
        boots.append(_wm(bx[idx], by[idx], sy[idx]))
    return point, np.std(boots, ddof=1)


def egger(dd):
    bx, by, sy = dd.bx.values, dd.by.values, dd.sy.values
    w = 1 / sy**2
    sw = np.sum(w)
    xw = np.sum(w * bx) / sw
    yw = np.sum(w * by) / sw
    num = np.sum(w * (bx - xw) * (by - yw))
    den = np.sum(w * (bx - xw)**2)
    slope = num / den
    intercept = yw - slope * xw
    se_sl = np.sqrt(1 / den)
    Q = np.sum(w * (by - (intercept + slope * bx))**2)
    return slope, se_sl, intercept, Q


def raps(dd):
    # robust adjusted profile score, simplified iterative (Huber loss)
    bx, by, sy = dd.bx.values, dd.by.values, dd.sy.values
    b = ivw(dd)[0]
    for _ in range(100):
        resid = by - b * bx
        s = np.sqrt(np.median(sy**2))
        t = resid / s
        psi = np.where(np.abs(t) <= 1.345, resid, 1.345 * s * np.sign(t))
        nb = b + np.sum((1 / sy**2) * bx * psi) / np.sum((1 / sy**2) * bx**2)
        if abs(nb - b) < 1e-7:
            break
        b = nb
    se = np.sqrt(1 / np.sum((1 / sy**2) * bx**2))
    return b, se


# ---------------- liability scale & Steiger ----------------
def r2_continuous(beta, maf):
    """Point-estimate formula 2 MAF (1-MAF) beta^2; unstable for rare/heavy-tailed
    variants (inflated beta), kept only as a sensitivity comparator."""
    return 2 * maf * (1 - maf) * beta**2


def r2_from_pn(beta, se, n):
    """TwoSampleMR get_r_from_pn: variance explained from the observed
    association strength z; robust to inflated point estimates. Preferred R2."""
    z2 = (beta / se) ** 2
    return z2 / (z2 + n)


def r2_from_lor(lor, af, K, ncase, ncontrol):
    """TwoSampleMR get_r_from_lor (Lee 2012 framework): correct the
    case-enriched sample AF to population AF, then V_G/(V_G+pi^2/3)."""
    from scipy.optimize import brentq
    prop = ncase / (ncase + ncontrol)
    OR = np.exp(lor)

    def eq(f0):
        o0 = f0 / (1 - f0)
        f1 = OR * o0 / (1 + OR * o0)
        return prop * f1 + (1 - prop) * f0 - af
    f0 = brentq(eq, 1e-6, 1 - 1e-6)
    o0 = f0 / (1 - f0)
    f1 = OR * o0 / (1 + OR * o0)
    p = K * f1 + (1 - K) * f0          # population effect allele frequency
    VG = lor**2 * p * (1 - p)
    return VG / (VG + np.pi**2 / 3)


def steiger_r2(bx, sx, maf_x, by, af_y, K_y, ncase, ncontrol,
               n_x, x_binary=False, K_x=None):
    # Exposure R2 from observed association strength (get_r_from_pn); robust
    # for the sparse, heavy-tailed CLR mycobiome traits.
    r2x = r2_from_pn(bx, sx, n_x)
    r2y = r2_from_lor(by, af_y, K_y, ncase, ncontrol)
    return r2x, r2y
