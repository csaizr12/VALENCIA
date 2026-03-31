import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Necesario para ejecución en pipelines/servidores
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

def generate_quality_panel(gff_path, output_png):
    """
    Genera el reporte de calidad VALENCIA con Paneles A y B detallados.
    Incluye detección de frameshifts y marginales de alta resolución.
    """
    print(f"🧬 Processing GFF for Quality Report: {gff_path}")
    data = []
    
    try:
        with open(gff_path, 'r', encoding='latin-1') as f:
            for line in f:
                if line.startswith('#') or '\t' not in line: continue
                fields = line.split('\t')
                if len(fields) < 9: continue
                
                if fields[2].lower() in ['mrna', 'transcript']:
                    attr = fields[8]
                    # Extracción robusta (insensible a mayúsculas)
                    tx = re.search(r'transcripts_edit_distance=([0-9.]+)', attr, re.I)
                    pr = re.search(r'proteins_edit_distance=([0-9.]+)', attr, re.I)
                    cds = re.search(r'cds_edit_distance=([0-9.]+)', attr, re.I)
                    
                    if tx and pr and cds:
                        t_v, p_v, c_v = float(tx.group(1)), float(pr.group(1)), float(cds.group(1))
                        data.append({
                            'AED_tx': t_v, 'AED_pr': p_v, 'AED_cds': c_v, 
                            'Delta': abs(p_v - c_v)
                        })
        
        df = pd.DataFrame(data)
        if df.empty:
            print("⚠️ No valid AED metrics found in GFF.")
            return

        # --- CONFIGURACIÓN ESTÉTICA ---
        plt.rcParams['savefig.dpi'] = 900 
        sns.set_theme(style="white")
        fig = plt.figure(figsize=(28, 14))
        
        # Layout: 2 columnas principales
        gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.45)

        # --- 1. PLOT A: STRUCTURE VS FUNCTION (PUNTOS + MARGINALES) ---
        inner_gs = gs[0].subgridspec(4, 4, hspace=0.12, wspace=0.12)
        ax_main = fig.add_subplot(inner_gs[1:, :-1])
        ax_hist_x = fig.add_subplot(inner_gs[0, :-1], sharex=ax_main)
        ax_hist_y = fig.add_subplot(inner_gs[1:, -1], sharey=ax_main)

        # Scatter principal
        sc = ax_main.scatter(df['AED_tx'], df['AED_pr'], c=df['Delta'], 
                            s=4, cmap="YlOrRd", vmin=0, vmax=1, 
                            edgecolors='none', alpha=0.7, rasterized=True)

        # Histogramas Marginales
        ax_hist_x.hist(df['AED_tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1)
        ax_hist_y.hist(df['AED_pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        
        ax_hist_x.set_ylabel('Nb. Transcripts', fontsize=10)
        ax_hist_y.set_xlabel('Nb. Transcripts', fontsize=10)
        ax_hist_x.tick_params(labelbottom=False); ax_hist_y.tick_params(labelleft=False)

        # Umbrales y Ejes
        ax_main.axvline(0.3, color='red', linestyle='--', alpha=0.5, linewidth=1.2)
        ax_main.axhline(0.1, color='red', linestyle='--', alpha=0.5, linewidth=1.2)
        ax_main.set_xlim(-0.01, 1.01); ax_main.set_ylim(-0.01, 1.01)
        ax_main.set_xlabel('AED Transcripts (Structure Match)', fontweight='bold', fontsize=14)
        ax_main.set_ylabel('AED Proteins (Function Match)', fontweight='bold', fontsize=14)

        # Colorbar Plot A
        divider = make_axes_locatable(ax_main)
        cax = divider.append_axes("bottom", size="3.5%", pad=0.85)
        cbar = fig.colorbar(sc, cax=cax, orientation='horizontal')
        cbar.set_label('Δ AED (Max Discrepancy CDS vs Protein | 1.0 = Potential Frameshift)', fontweight='bold', fontsize=11)

        # Leyenda Plot A
        leg_elements = [Patch(facecolor='#45a049', label='AED transcripts'),
                        Patch(facecolor='#e91e63', label='AED proteins')]
        ax_hist_y.legend(handles=leg_elements, loc='upper left', bbox_to_anchor=(1.1, 1.0), frameon=True, fontsize=11)

        # --- 2. PLOT B: INTERNAL CONSISTENCY (CORRELATION) ---
        ax_corr = fig.add_subplot(gs[1])
        sns.scatterplot(data=df, x='AED_cds', y='AED_pr', alpha=0.15, s=5, 
                        color='#34495e', ax=ax_corr, edgecolor=None, rasterized=True)
        
        # Línea de Identidad y Zona de Frameshift
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', linewidth=1.5, label='Internal Consistency (X=Y)')
        ax_corr.fill_between([0, 1], [0, 1], [1, 1], color='red', alpha=0.04, label='Potential Frameshift Zone')
        
        ax_corr.set_xlabel('AED CDS (Nucleotide level)', fontweight='bold', fontsize=14)
        ax_corr.set_ylabel('AED Proteins (Amino acid level)', fontweight='bold', fontsize=14)
        ax_corr.set_title('B. Correlation: CDS vs Protein Distance', fontsize=18, fontweight='bold')
        ax_corr.set_xlim(-0.01, 1.01); ax_corr.set_ylim(-0.01, 1.01)
        ax_corr.legend(loc='upper left', fontsize=12)
        ax_corr.grid(True, linestyle=':', alpha=0.4)

      
        # TÍTULO GLOBAL
        fig.suptitle(f'VALENCIA Annotation Quality Analysis (n={len(df)})', fontsize=26, fontweight='bold', y=0.98)
        
        # GUARDADO
        plt.savefig(output_png, bbox_inches='tight')
        plt.close()
        print(f"✅ Master Quality Report generated successfully at: {output_png}")

    except Exception as e:
        print(f"❌ Error in Plotter: {e}")