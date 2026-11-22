import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Publication style (Frontiers)
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial','DejaVu Sans','Helvetica'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.linewidth': 0.8,
    'axes.grid': True,
    'grid.linewidth': 0.4,
    'grid.alpha': 0.35,
    'grid.color': '#999999',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# Color-blind friendly palette
C_NAVY   = '#2c3e50'
C_BLUE   = '#3b7dd8'
C_RED    = '#d64545'
C_ORANGE = '#e8833a'
C_GREEN  = '#2e9e6b'
C_GREY   = '#9aa3ad'
C_PURPLE = '#7d5ba6'
