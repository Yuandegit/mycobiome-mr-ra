import json, sys
import numpy as np
import pandas as pd
import scipy.stats as st
sys.path.insert(0, 'scripts')
import mr_core as mc

ROOT = '/home/user/mycobiome_joint_mr'
N_EXP = 7350

# ---- mapping: scale + genus name ----
mapping = json.load(open(f'{ROOT}/data/gwas_microbe/fungi_genus/mapping.json'))
def genus_of(v):
    for part in v.split('.'):
        if part.startswith('g__'):
            return part[3:]
    return v
meta = {k: ('Continuous' if v.startswith('Continuous') else 'Binary', genus_of(v))
        for k, v in mapping.items()}

harm = pd.read_pickle(f'{ROOT}/results/harmonized_forward.pkl')
fwd = pd.read_csv(f'{ROOT}/results/forward_mr_results.csv')
fwd = fwd[fwd.threshold == 'expl1e5'].set_index('pheno_idx')

rows = []
for idx in range(115):
    key = (f'{idx:03d}', 'expl1e5')
    k = harm[key]
    k = k[k.keep].reset_index(drop=True)
    scale, genus = meta[f'{idx:03d}']
    # per-SNP F and method-B R2
    zx = k.bx / k.sx
    F = zx**2
    r2 = mc.r2_from_pn(k.bx.values, k.sx.values, N_EXP)
    r = fwd.loc[idx]
    # weighted median bootstrap (correct SE)
    wm_b, wm_se = mc.weighted_median(k)
    wm_p = 2 * (1 - st.norm.cdf(abs(wm_b / wm_se)))
    # Egger intercept and p from a fresh egger fit (TwoSampleMR R-lm
    # convention: residual-scale sigma from Egger Q / (n-2), normal p)
    eg_slope, eg_se_sl, eg_int, eg_Q = mc.egger(k)
    n_k = len(k)
    w = 1 / k.sy.values**2
    sw = np.sum(w)
    xw = np.sum(w * k.bx.values) / sw
    den = np.sum(w * (k.bx.values - xw)**2)
    fe_int_se = np.sqrt(1/sw + xw**2/den)          # fixed scale=1
    if n_k > 2:
        sigma = np.sqrt(eg_Q / (n_k - 2))          # R lm residual scale
    else:
        sigma = 1.0
    eg_int_se = fe_int_se * sigma
    eg_int_p = 2 * (1 - st.norm.cdf(abs(eg_int / eg_int_se)))
    rows.append(dict(
        id=idx, genus=genus, scale=scale, n_iv=len(k),
        F_med=np.median(F), F_min=np.min(F),
        r2_per_med=np.median(r2), r2_total=np.sum(r2),
        ivw_b=r.b, ivw_se=r.se, ivw_p=r.p, fdr=r.FDR,
        wm_b=wm_b, wm_p=wm_p,
        eg_int=eg_int, eg_int_p=eg_int_p,
        Q=r.Q, Qp=r.Qp, mde80=2.8*r.se,
    ))

sup = pd.DataFrame(rows)
sup.to_csv(f'{ROOT}/results/suppl_115_traits.csv', index=False)

# ---- Table 1: aggregate by scale ----
def agg_block(d):
    return dict(
        n=len(d),
        iv_med=d.n_iv.median(), iv_q1=d.n_iv.quantile(.25), iv_q3=d.n_iv.quantile(.75),
        F_med=d.F_med.median(), F_min=d.F_min.min(), F_max=d.F_med.max(),
        r2_per=d.r2_per_med.median(), r2_total=d.r2_total.median(),
        mde_med=d.mde80.median(), mde_q1=d.mde80.quantile(.25), mde_q3=d.mde80.quantile(.75),
        n_nom=(d.ivw_p < .05).sum(), n_fdr=(d.fdr < .05).sum(),
    )
t1 = []
for label, d in [('Continuous', sup[sup.scale=='Continuous']),
                 ('Binary', sup[sup.scale=='Binary']),
                 ('All', sup)]:
    a = agg_block(d); a['scale'] = label; t1.append(a)
t1 = pd.DataFrame(t1).set_index('scale')
t1.to_csv(f'{ROOT}/results/table1_aggregate.csv')
print(t1.to_string())
print('\nsuppl rows:', len(sup))
