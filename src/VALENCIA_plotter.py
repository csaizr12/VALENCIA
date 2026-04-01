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
    data = []
    try:
        # --- 1. EXTRACCIÓN DE DATOS ---
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
        if df.empty: 
            print(f"⚠️ No se encontraron datos en {gff_path}")
            return

        # --- 2. CONFIGURACIÓN DE FIGURA ---
        plt.rcParams['savefig.dpi'] = 900 
        sns.set_theme(style="white")
        fig = plt.figure(figsize=(36, 14))
        
        # Grid principal para los 3 paneles
        gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.4)

        # --- 3. PANEL A: ESTRUCTURA VS FUNCIÓN ---
        # 7 filas para control total de espacios
        inner_gs = gs[0].subgridspec(7, 4, 
                                     height_ratios=[1.2, 1, 1, 1, 1, 0.6, 0.2], 
                                     hspace=0.05, 
                                     wspace=0.05)
        
        ax_main = fig.add_subplot(inner_gs[1:5, :-1])
        ax_hist_x = fig.add_subplot(inner_gs[0, :-1], sharex=ax_main)
        ax_hist_y = fig.add_subplot(inner_gs[1:5, -1], sharey=ax_main)

        # Scatter con Magma
        sc = ax_main.scatter(df['tx'], df['pr'], c=df['Delta'], 
                            s=5, cmap="magma", vmin=0, vmax=1, alpha=0.7, rasterized=True)

        # Histogramas marginales
        ax_hist_x.hist(df['tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1)
        ax_hist_y.hist(df['pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        
        ax_hist_x.set_ylabel('Nb. transcripts', fontsize=10)
        ax_hist_y.set_xlabel('Nb. transcripts', fontsize=10)
        ax_hist_x.tick_params(labelbottom=False, bottom=False)
        ax_hist_y.tick_params(labelleft=False, left=False)

        ax_main.set_xlabel('Lev_edit_distance transcripts', fontweight='bold', fontsize=12, labelpad=25)
        ax_main.set_ylabel('Lev_edit_distance proteins', fontweight='bold', fontsize=12)

        # Colorbar fina
        cax_a = fig.add_subplot(inner_gs[6, :-1])
        cbar = fig.colorbar(sc, cax=cax_a, orientation='horizontal', aspect=50)
        cbar.set_label('Δ Lev_edit_distance', fontweight='bold', fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        # LEYENDA PANEL A (Corregida y movida arriba a la derecha del panel rosa)
        leg_elements = [
            Patch(facecolor='#45a049', label='Lev_edit_distance transcripts'),
            Patch(facecolor='#e91e63', label='Lev_edit_distance proteins')
        ]
        ax_hist_y.legend(handles=leg_elements, loc='upper left', 
                         bbox_to_anchor=(1.05, 1.15), frameon=True, fontsize=11)

        # --- 4. PANEL B: CORRELACIÓN ---
        ax_corr = fig.add_subplot(gs[1])
        sns.scatterplot(data=df, x='cds', y='pr', alpha=0.15, s=6, color='#34495e', ax=ax_corr, rasterized=True)
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identity (X=Y)')
        ax_corr.fill_between([0, 1], [0, 1], [1, 1], color='red', alpha=0.04, label='Potential frameshift zone')
        
        ax_corr.set_title('Correlation: CDS vs protein', fontsize=18, fontweight='bold', pad=20)
        ax_corr.set_xlabel('Lev_edit_distance CDS (Nucleotide level)', fontweight='bold', fontsize=12)
        ax_corr.set_ylabel('Lev_edit_distance Proteins (Amino acid level)', fontweight='bold', fontsize=12)
        ax_corr.legend(loc='upper left', fontsize=11)
        ax_corr.grid(True, linestyle=':', alpha=0.4)

        # --- 5. PANEL C: DISTRIBUCIÓN ---
        ax_dist = fig.add_subplot(gs[2])
        sns.histplot(df['pr'] - df['cds'], bins=60, kde=True, color='#2E86C1', ax=ax_dist, edgecolor='white')
        
        ax_dist.set_title('Distribution of editing difference', fontsize=18, fontweight='bold', pad=20)
        ax_dist.set_xlabel('Lev_edit_distance difference (protein - CDS)', fontweight='bold', fontsize=12)
        ax_dist.set_ylabel('Number of genes (Frequency)', fontweight='bold', fontsize=12)
        
        ax_dist.axvline(0, color='red', linestyle=':', label='Zero difference')
        ax_dist.grid(True, linestyle=':', alpha=0.4)
        ax_dist.legend(loc='upper right', fontsize=11)

        # Título General
        fig.suptitle(f'VALENCIA Annotation Quality Analysis (n={len(df)})', fontsize=28, fontweight='bold', y=0.98)
        
        # --- 6. GUARDADO ESTABLE (Sin bbox_inches='tight' para evitar el error) ---
        # Dejamos márgenes manuales para que nada se corte
        plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.9, hspace=0.3, wspace=0.3)
        
        plt.savefig(output_png)
        plt.close()
        print(f"✅ Panel generado correctamente: {output_png}")

    except Exception as e:
        print(f"❌ Error crítico en generate_quality_panel: {e}")