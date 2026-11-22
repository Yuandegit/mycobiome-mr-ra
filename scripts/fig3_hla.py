import sys; sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from plot_style import *

d = pd.read_pickle('results/harmonized_forward.pkl')
df = d[('109', 'expl1e8')] if ('109','expl1e8') in d else d[('109','expl1e5')]
k = df[df.keep].copy().reset_index(drop=True)
k['name'] = [r if isinstance(r, str) else f"chr{c}:{p}" for r, c, p in zip(k.rsid, k.chrom, k.pos)]
k['wald'] = k.by / k.bx
k['wse'] = k.sy / np.abs(k.bx)
k['w'] = 1.0 / (k.sy**2 / k.bx**2)          # IVW weight
driver = k.name == 'rs111394117'

# leave-one-out IVW
def ivw(sub):
    b = np.sum(sub.w * sub.wald) / np.sum(sub.w)
    se = np.sqrt(1.0 / np.sum(sub.w))
    return b, se
loo_b, loo_se = [], []
for i in range(len(k)):
    b, se = ivw(k.drop(i))
    loo_b.append(b); loo_se.append(se)
k['loo'] = loo_b; k['loo_se'] = loo_se

fig = plt.figure(figsize=(7.2, 6.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.42, wspace=0.38)

# ---- Panel A: per-SNP Wald forest ----
axA = fig.add_subplot(gs[0, 0])
order = np.arange(len(k))[::-1]
for i, (_, row) in enumerate(k.iterrows()):
    yp = order[i]
    isd = row['name'] == 'rs111394117'
    col = C_RED if isd else C_BLUE
    lo, hi = row.wald - 1.96*row.wse, row.wald + 1.96*row.wse
    axA.plot([lo, hi], [yp, yp], lw=1.0, color=col)
    axA.scatter(row.wald, yp, s=18 if isd else 12, c=col, zorder=5)
# pooled IVW
pooled, pooled_se = ivw(k)
axA.axvline(pooled, ls='--', lw=1.0, color=C_NAVY, alpha=0.8)
axA.axvline(0, lw=0.7, color='#444444')
axA.set_yticks(order); axA.set_yticklabels(k.name, fontsize=5.6)
axA.set_xlabel('Per-SNP Wald ratio (95% CI)')
axA.set_title('A  Per-SNP Wald ratios (22 IVs); dashed = pooled IVW', loc='left', fontweight='bold')
axA.grid(axis='y', visible=False)
axA.set_xlim(-2.2, 1.6)

# ---- Panel B: leave-one-out ----
axB = fig.add_subplot(gs[0, 1])
orderB = np.arange(len(k))[::-1]
for i, (_, row) in enumerate(k.iterrows()):
    yp = orderB[i]
    isd = row['name'] == 'rs111394117'
    col = C_RED if isd else C_GREEN
    lo, hi = row.loo - 1.96*row.loo_se, row.loo + 1.96*row.loo_se
    axB.plot([lo, hi], [yp, yp], lw=1.0 if isd else 0.7, color=col)
    axB.scatter(row.loo, yp, s=18 if isd else 12, c=col, zorder=5)
axB.axvline(pooled, ls='--', lw=1.0, color=C_NAVY, alpha=0.8)
axB.axvline(0, lw=0.7, color='#444444')
axB.set_yticks(orderB); axB.set_yticklabels(k.name, fontsize=5.6)
axB.set_xlabel('IVW estimate after excluding the SNP (95% CI)')
axB.set_title('B  Leave-one-out; excluding rs111394117 → β=+0.012 (null)', loc='left', fontweight='bold')
axB.grid(axis='y', visible=False)
axB.set_xlim(-0.22, 0.08)

# ---- Panel C: funnel ----
axC = fig.add_subplot(gs[1, 0])
prec = 1.0/k.wse
axC.scatter(k.wald[~driver], prec[~driver], s=16, c=C_BLUE, zorder=4)
axC.scatter(k.wald[driver], prec[driver], s=40, c=C_RED, edgecolors='white', linewidth=0.5, zorder=5)
axC.axvline(pooled, ls='--', lw=1.0, color=C_NAVY)
axC.axvline(0, lw=0.7, color='#444444')
axC.annotate('rs111394117', xy=(k.wald[driver].iloc[0], prec[driver].iloc[0]),
             xytext=(-1.7, 3.2), fontsize=6.5, color=C_RED,
             arrowprops=dict(arrowstyle='->', lw=0.7, color=C_RED))
axC.set_xlabel('Wald ratio'); axC.set_ylabel('Precision (1 / SE)')
axC.set_title('C  Funnel plot (asymmetry driven by the HLA SNP)', loc='left', fontweight='bold')
axC.set_xlim(-2.2, 1.6)

# ---- Panel D: HLA region schematic ----
axD = fig.add_subplot(gs[1, 1])
axD.set_xlim(24, 36); axD.set_ylim(0, 10); axD.axis('off')
# HLA region band
axD.add_patch(Rectangle((25, 6.4), 10, 1.2, facecolor=C_GREY, alpha=0.25, edgecolor=C_NAVY, lw=0.8))
axD.text(30, 7.9, 'HLA complex (chr6: 25–35 Mb)', ha='center', fontsize=7, color=C_NAVY)
# NOTCH4 gene
axD.add_patch(Rectangle((31.6, 4.4), 0.55, 0.9, facecolor=C_PURPLE, alpha=0.8, edgecolor='none'))
axD.text(31.88, 4.05, 'NOTCH4', ha='center', fontsize=6.5, color=C_PURPLE)
# driver SNP
sx = 32.187
axD.scatter(sx, 6.4, s=46, c=C_RED, zorder=6, edgecolors='white', linewidth=0.5)
axD.annotate('rs111394117\n(chr6:32,187,213)', xy=(sx, 6.5), xytext=(32.6, 8.6),
             fontsize=6.5, color=C_RED, arrowprops=dict(arrowstyle='->', lw=0.7, color=C_RED))
# LD proxies
axD.text(30, 2.75, 'JPT LD (1000G phase 3, n=104)\n66 proxies r² ≥ 0.8 (4 at r²=1.0); 1,394 total\nmarker of an extended haplotype; NOTCH4 is a marker,\nnot an established causal gene (clumping used EAS, n=504)',
         ha='center', va='top', fontsize=6.2, color='#333333', linespacing=1.5)
axD.set_title('D  Driver in NOTCH4 intron (HLA class III)', loc='left', fontweight='bold')

fig.savefig('figures/Figure3_hla_dissection.pdf')
fig.savefig('figures/Figure3_hla_dissection.png')
print('saved Figure3; pooled IVW =', round(pooled,4), '; LOO driver =', round(k[driver].loo.iloc[0],4))
