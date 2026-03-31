import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Modo servidor (sin ventana)
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

def generate_quality_panel(gff_path, output_png="VALENCIA_Quality_Report.png"):
    data = []
    try:
        with open(gff_path, 'r', encoding='latin-1') as f:
            for line in f:
                if line.startswith('#') or '\t' not in line: continue
                fields = line.split('\t')
                if len(fields) < 9: continue
                if fields[2].lower() in ['mrna', 'transcript']:
                    attr = fields[8]
                    tx = re.search(r'transcripts_edit_distance=([0-9.]+)', attr, re.I)
                    pr = re.search(r'proteins_edit_distance=([0-9.]+)', attr, re.I)
                    cds = re.search(r'cds_edit_distance=([0-9.]+)', attr, re.I)
                    if tx and pr and cds:
                        t_v, p_v, c_v = float(tx.group(1)), float(pr.group(1)), float(cds.group(1))
                        data.append({'AED_tx': t_v, 'AED_pr': p_v, 'AED_cds': c_v, 'Delta': abs(p_v - c_v)})
        
        df = pd.DataFrame(data)
        if df.empty:
            print("⚠️ No AED metrics found to plot.")
            return

        plt.rcParams['savefig.dpi'] = 900 
        sns.set_theme(style="white")
        fig = plt.figure(figsize=(32, 12))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.4)

        # PLOT A: Structure vs Function
        inner_gs = gs[0].subgridspec(4, 4, hspace=0.12, wspace=0.12)
        ax_main = fig.add_subplot(inner_gs[1:, :-1]); ax_hx = fig.add_subplot(inner_gs[0, :-1], sharex=ax_main); ax_hy = fig.add_subplot(inner_gs[1:, -1], sharey=ax_main)
        sc = ax_main.scatter(df['AED_tx'], df['AED_pr'], c=df['Delta'], s=6, cmap="YlOrRd", vmin=0, vmax=1, alpha=0.7, rasterized=True)
        ax_hx.hist(df['AED_tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1); ax_hy.hist(df['AED_pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        ax_main.axvline(0.3, color='red', linestyle='--', alpha=0.5); ax_main.axhline(0.1, color='red', linestyle='--', alpha=0.5)
        ax_main.set_xlabel('AED Transcripts (Structure)'); ax_main.set_ylabel('AED Proteins (Function)')
        fig.colorbar(sc, cax=make_axes_locatable(ax_main).append_axes("bottom", size="4%", pad=0.8), orientation='horizontal').set_label('Δ AED (Discrepancy)')

        # PLOT B: Consistency Scatter
        ax_corr = fig.add_subplot(gs[1])
        sns.scatterplot(data=df, x='AED_cds', y='AED_pr', alpha=0.15, s=5, color='#34495e', ax=ax_corr, rasterized=True)
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identity (X=Y)')
        ax_corr.set_title('B. Correlation: CDS vs Protein', fontsize=16, fontweight='bold')
        ax_corr.set_xlabel('AED CDS'); ax_corr.set_ylabel('AED Proteins'); ax_corr.legend()

        # PLOT C: Hexbin Density
        ax_dens = fig.add_subplot(gs[2])
        hb = ax_dens.hexbin(df['AED_cds'], df['AED_pr'], gridsize=40, cmap="magma", mincnt=1, bins='log')
        ax_dens.set_title('C. Point Density (CDS vs Protein)', fontsize=16, fontweight='bold')
        ax_dens.set_xlabel('AED CDS'); ax_dens.set_ylabel('AED Proteins')
        fig.colorbar(hb, cax=make_axes_locatable(ax_dens).append_axes("bottom", size="4%", pad=0.8), orientation='horizontal').set_label('Log10(Gene Count)')

        fig.suptitle(f'VALENCIA Quality Analysis (n={len(df)})', fontsize=24, fontweight='bold', y=0.98)
        plt.savefig(output_png, bbox_inches='tight')
        plt.close()
        print(f"✅ Report saved: {output_png}")

    except Exception as e:
        print(f"❌ Plotting error: {e}")