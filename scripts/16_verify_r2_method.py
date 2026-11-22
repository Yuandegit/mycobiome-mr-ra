import sys; sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd, json
from mr_core import r2_continuous, r2_from_lor

N_EXP = 7350
NCASE, NCONTROL = 5348, 173268
K = 0.01
mapping = json.load(open('data/gwas_microbe/fungi_genus/mapping.json'))

d = pd.read_pickle('results/harmonized_forward.pkl')
rows = []
for (idx, thr), df in d.items():
    if thr != 'expl1e5':
        continue
    is_bin = mapping[str(idx).zfill(3)].startswith('Binary')
    k = df[df.keep]
    for _, r in k.iterrows():
        z = r.bx / r.sx
        r2x_B = z**2 / (z**2 + N_EXP)          # get_r_from_pn (z/p based)
        r2x_A = r2_continuous(r.bx, r.maf)     # 2 MAF (1-MAF) beta^2
        r2y = r2_from_lor(r.by, r.eaf, K, NCASE, NCONTROL)
        rows.append((idx, is_bin, r.rsid, r2x_A, r2x_B, r2y,
                     r2x_A > r2y, r2x_B > r2y))

v = pd.DataFrame(rows, columns=['idx','binary','rsid','r2x_A','r2x_B','r2y',
                                'fwd_A','fwd_B'])
print('total kept IVs:', len(v))
print('--- Method A (point estimate, current) ---')
print('  forward (r2x>r2y):', v.fwd_A.sum(), ' reverse:', (~v.fwd_A).sum())
print('--- Method B (z/p based, TwoSampleMR get_r_from_pn) ---')
print('  forward (r2x>r2y):', v.fwd_B.sum(), ' reverse:', (~v.fwd_B).sum())

# continuous only
c = v[~v.binary]
print('\ncontinuous IVs:', len(c), ' method-B reverse:', (~c.fwd_B).sum())
b = v[v.binary]
print('binary IVs:', len(b), ' method-B reverse:', (~b.fwd_B).sum())

# driver SNP
drv = v[v.rsid == 'rs111394117']
print('\nDriver rs111394117:')
for _, r in drv.iterrows():
    print(f"  r2x_A={r.r2x_A:.4f}  r2x_B={r.r2x_B:.4f}  r2y(liability)={r.r2y:.4f}")
    print(f"  method A forward={r.fwd_A}  method B forward={r.fwd_B}")

# how many distinct phenotypes flip
flips = v[~v.fwd_B].idx.unique()
print('\nphenotypes with >=1 method-B reverse IV:', len(flips))
v.to_csv('results/r2_method_comparison.csv', index=False)
print('saved results/r2_method_comparison.csv')
