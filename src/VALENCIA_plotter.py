import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

def generate_quality_panel(gff_path, output_png):
    print(f"📊 Generating Improved Quality Report from: {gff_path}")
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
                        data.append({
                            'tx': t_v, 'pr': p_v, 'cds': c_v, 
                            'Delta': abs(p_v - c_v)
                        })
        
        df = pd.DataFrame(data)
        if df.empty: return

        plt.rcParams['savefig.dpi'] = 900 
        sns.set_theme(style="white")
        fig = plt.figure(figsize=(36, 12))
        
        # Layout: 3 columnas
        gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.35)

        # --- 1. PLOT A: AJUSTE DE TAMAÑO ---
        # Usamos hspace y wspace muy bajos para que los histogramas "peguen" con el centro
        inner_gs = gs[0].subgridspec(4, 4, hspace=0.05, wspace=0.05)
        ax_main = fig.add_subplot(inner_gs[1:, :-1])
        ax_hist_x = fig.add_subplot(inner_gs[0, :-1], sharex=ax_main)
        ax_hist_y = fig.add_subplot(inner_gs[1:, -1], sharey=ax_main)

        sc = ax_main.scatter(df['tx'], df['pr'], c=df['Delta'], 
                            s=4, cmap="YlOrRd", vmin=0, vmax=1, alpha=0.7, rasterized=True)

        ax_hist_x.hist(df['tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1)
        # El histograma lateral ahora ocupa exactamente el mismo espacio vertical que ax_main
        ax_hist_y.hist(df['pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        
        ax_hist_x.tick_params(labelbottom=False); ax_hist_y.tick_params(labelleft=False)
        ax_main.axvline(0.1, color='red', linestyle='--', alpha=0.5)
        ax_main.axhline(0.1, color='red', linestyle='--', alpha=0.5)
        ax_main.set_xlabel('lev_edit_distance transcripts', fontweight='bold')
        ax_main.set_ylabel('lev_edit_distance proteins', fontweight='bold')

        divider = make_axes_locatable(ax_main)
        cax = divider.append_axes("bottom", size="3.5%", pad=0.8)
        fig.colorbar(sc, cax=cax, orientation='horizontal').set_label('Δ lev_edit_distance (Discrepancy)')

        # --- 2. PLOT B: CORRELACIÓN ---
        ax_corr = fig.add_subplot(gs[1])
        sns.scatterplot(data=df, x='cds', y='pr', alpha=0.15, s=5, color='#34495e', ax=ax_corr, rasterized=True)
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Internal consistency (X=Y)')
        ax_corr.fill_between([0, 1], [0, 1], [1, 1], color='red', alpha=0.04, label='Potential frameshift zone')
        ax_corr.set_title('B. Correlation: CDS vs Protein', fontsize=18, fontweight='bold')
        ax_corr.set_xlabel('CDS Distance', fontweight='bold'); ax_corr.set_ylabel('Protein Distance', fontweight='bold')
        ax_corr.legend(loc='upper left')

        # --- 3. PLOT C: HISTOGRAMA DE DISTRIBUCIÓN ---
        ax_dist = fig.add_subplot(gs[2])
        # Graficamos la distribución de la diferencia (Protein - CDS)
        # Un histograma con línea de densidad (KDE) es mucho más limpio
        sns.histplot(df['pr'] - df['cds'], bins=50, kde=True, color='#2E86C1', ax=ax_dist, edgecolor='white')
        
        ax_dist.set_title('C. Distribution of Editing Difference', fontsize=18, fontweight='bold')
        ax_dist.set_xlabel('Difference (Protein Distance - CDS Distance)', fontweight='bold')
        ax_dist.set_ylabel('Frequency (Number of Genes)', fontweight='bold')
        ax_dist.axvline(0, color='red', linestyle=':', label='Zero Difference')
        ax_dist.grid(True, linestyle=':', alpha=0.4)
        ax_dist.legend()

        fig.suptitle(f'VALENCIA Annotation quality analysis (n={len(df)})', fontsize=26, fontweight='bold', y=0.98)
        plt.savefig(output_png, bbox_inches='tight')
        plt.close()
        print(f"✅ Panel generated: {output_png}")

    except Exception as e:
        print(f"❌ Error: {e}")