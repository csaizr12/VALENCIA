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
    print(f"📊 Generando Panel Maestro (Colorbar Fina) desde: {gff_path}")
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
                        data.append({'tx': t_v, 'pr': p_v, 'cds': c_v, 'Delta': abs(p_v - c_v)})
        
        df = pd.DataFrame(data)
        if df.empty: return

        plt.rcParams['savefig.dpi'] = 900 
        sns.set_theme(style="white")
        fig = plt.figure(figsize=(36, 14))
        
        gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.4)

       # --- 1. PLOT A: DESPLAZAMIENTO DE LEYENDA (COLORBAR) ---
        # Aumentamos ligeramente el hspace final para dar aire
        inner_gs = gs[0].subgridspec(6, 4, 
                                     height_ratios=[1.2, 1, 1, 1, 1, 0.2], 
                                     hspace=0.05, 
                                     wspace=0.05)
        
        ax_main = fig.add_subplot(inner_gs[1:5, :-1])
        ax_hist_x = fig.add_subplot(inner_gs[0, :-1], sharex=ax_main)
        ax_hist_y = fig.add_subplot(inner_gs[1:5, -1], sharey=ax_main)

        # Dibujado con Magma
        sc = ax_main.scatter(df['tx'], df['pr'], c=df['Delta'], 
                            s=5, cmap="magma", vmin=0, vmax=1, alpha=0.7, rasterized=True)

        ax_hist_x.hist(df['tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1)
        ax_hist_y.hist(df['pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        
        # Etiquetas de los ejes
        ax_hist_x.set_ylabel('Nb. transcripts', fontsize=10)
        ax_hist_y.set_xlabel('Nb. transcripts', fontsize=10)
        ax_hist_x.tick_params(labelbottom=False, bottom=False)
        ax_hist_y.tick_params(labelleft=False, left=False)
        
        # --- EL TRUCO ESTÁ AQUÍ ---
        # Añadimos un padding mayor al label del eje X para que no choque con la barra
        ax_main.set_xlabel('lev_edit_distance transcripts', fontweight='bold', fontsize=12, labelpad=15)
        ax_main.set_ylabel('lev_edit_distance proteins', fontweight='bold', fontsize=12)

        # Ubicación de la Colorbar con espacio extra
        cax_a = fig.add_subplot(inner_gs[5, :-1])
        # Usamos pad en el label de la colorbar para que respire
        cbar = fig.colorbar(sc, cax=cax_a, orientation='horizontal', aspect=50)
        cbar.set_label('Δ lev_edit_distance (Discrepancy)', fontweight='bold', fontsize=10, labelpad=10)
        cbar.ax.tick_params(labelsize=8)

        # Leyenda lateral (AED transcripts/proteins)
        ax_hist_y.legend(handles=[Patch(facecolor='#45a049', label='AED transcripts'),
                                 Patch(facecolor='#e91e63', label='AED proteins')],
                         loc='upper left', bbox_to_anchor=(1.05, 1.0), frameon=True)
        # --- PLOTS B y C (Igual que antes) ---
        ax_corr = fig.add_subplot(gs[1])
        sns.scatterplot(data=df, x='cds', y='pr', alpha=0.15, s=6, color='#34495e', ax=ax_corr, rasterized=True)
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identity (X=Y)')
        ax_corr.fill_between([0, 1], [0, 1], [1, 1], color='red', alpha=0.04, label='Potential frameshift zone')
        ax_corr.set_title('B. Correlation: CDS vs Protein', fontsize=18, fontweight='bold')
        ax_corr.legend(loc='upper left')

        ax_dist = fig.add_subplot(gs[2])
        sns.histplot(df['pr'] - df['cds'], bins=60, kde=True, color='#2E86C1', ax=ax_dist, edgecolor='white')
        ax_dist.set_title('C. Distribution of Editing Difference', fontsize=18, fontweight='bold')
        ax_dist.axvline(0, color='red', linestyle=':', label='Zero Difference')
        ax_dist.legend()

        fig.suptitle(f'VALENCIA Annotation Quality Analysis (n={len(df)})', fontsize=28, fontweight='bold', y=0.98)
        
        plt.savefig(output_png, bbox_inches='tight')
        plt.close()
        print(f"✅ Panel generado con colorbar fina: {output_png}")

    except Exception as e:
        print(f"❌ Error: {e}")