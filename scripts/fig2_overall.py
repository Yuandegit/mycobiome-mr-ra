import sys; sys.path.insert(0, 'scripts')
import numpy as np, pandas as pd, json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from plot_style import *

r = pd.read_csv('results/forward_mr_results.csv')
mapping = json.load(open('data/gwas_microbe/fungi_genus/mapping.json'))
def ptype(idx):
    return 'Continuous' if mapping[str(idx).zfill(3)].startswith('Continuous') else 'Binary'

expl = r[r.threshold == 'expl1e5'].copy()
expl['ptype'] = [ptype(i) for i in expl.pheno_idx]
strict = r[r.threshold == 'gw5e8'].copy()

fig = plt.figure(figsize=(10.6, 3.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.18], wspace=0.62)


def volcano(ax, sub, tag, n_family):
    fdr_sig = sub.FDR < 0.05
    nom = (sub.p < 0.05) & (~fdr_sig)
    other = ~fdr_sig & ~nom
    for mask, col, s in [(other, C_GREY, 16), (nom, C_ORANGE, 24), (fdr_sig, C_RED, 52)]:
        ax.scatter(sub.b[mask], -np.log10(sub.p[mask]), s=s, c=col,
                   alpha=0.85 if col != C_GREY else 0.7,
                   edgecolors='white' if col == C_RED else 'none',
                   linewidth=0.5, zorder=5 if col == C_RED else 3)
    ax.axhline(-np.log10(0.05), ls='--', lw=0.7, color=C_NAVY, alpha=0.7)
    ax.axhline(-np.log10(0.05 / n_family), ls=':', lw=0.8, color=C_NAVY, alpha=0.9)
    ax.axvline(0, lw=0.7, color='#444444', alpha=0.6)
    ax.set_xlabel('IVW β (%s scale)' % tag)
    ax.set_ylabel(r'$-\log_{10}(P)$')
    return fdr_sig


# Panel A1: continuous CLR (58 features)
axA1 = fig.add_subplot(gs[0, 0])
cont = expl[expl.ptype == 'Continuous']
sig1 = volcano(axA1, cont, 'continuous CLR', len(cont))
axA1.set_title('A1  Continuous CLR (n=%d)' % len(cont), loc='left', fontweight='bold')
hit = cont[sig1]
if len(hit):
    h = hit.iloc[0]
    axA1.annotate('Cutaneotrichosporon\nHLA artifact (Fig. 3)',
                  xy=(h.b, -np.log10(h.p)), xytext=(h.b - 0.02, 13.2),
                  fontsize=6.0, color=C_RED, ha='left', va='top',
                  arrowprops=dict(arrowstyle='->', lw=0.7, color=C_RED))

# Panel A2: binary log-odds (57 features)
axA2 = fig.add_subplot(gs[0, 1])
bina = expl[expl.ptype == 'Binary']
volcano(axA2, bina, 'binary log-odds', len(bina))
axA2.set_title('A2  Binary detection (n=%d)' % len(bina), loc='left', fontweight='bold')

# Panel B: strict forest (point-and-interval)
axB = fig.add_subplot(gs[0, 2])
hits = strict[strict.p < 0.05].sort_values('b')
labels = list(hits.genus)
yy = np.arange(len(hits))[::-1]
for i, (_, row) in enumerate(hits.iterrows()):
    yp = yy[i]; lo, hi = row.b - 1.96*row.se, row.b + 1.96*row.se
    axB.plot([lo, hi], [yp, yp], lw=1.4, color=C_BLUE)
    axB.plot([lo, lo], [yp-0.12, yp+0.12], lw=1.0, color=C_BLUE)
    axB.plot([hi, hi], [yp-0.12, yp+0.12], lw=1.0, color=C_BLUE)
    axB.scatter(row.b, yp, s=26, c=C_BLUE, zorder=5)
axB.axvline(0, lw=0.8, color='#444444', alpha=0.7)
axB.set_yticks(yy); axB.set_yticklabels(labels)
axB.set_xlabel('Single-SNP Wald ratio (95% CI)')
axB.set_title('B  Strict single-IV hits (all FDR 0.34–0.35)', loc='left', fontweight='bold')
axB.grid(axis='y', visible=False)

handles = [Line2D([0],[0], marker='o', color='none', markerfacecolor=C_GREY, markersize=7, label='P ≥ 0.05'),
           Line2D([0],[0], marker='o', color='none', markerfacecolor=C_ORANGE, markersize=7, label='Nominal P<0.05'),
           Line2D([0],[0], marker='o', color='none', markerfacecolor=C_RED, markersize=8, label='FDR<0.05 (artifact)')]
axA2.legend(handles=handles, loc='upper right', frameon=True, framealpha=1.0,
            facecolor='white', edgecolor='#999999', fontsize=6.5,
            handletextpad=0.35, borderaxespad=0.4, borderpad=0.5, labelspacing=0.4)

fig.subplots_adjust(left=0.06, right=0.98, bottom=0.20, top=0.87)
fig.savefig('figures/Figure2_overall_mr.pdf')
fig.savefig('figures/Figure2_overall_mr.png')
print('saved Figure2')
