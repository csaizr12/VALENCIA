import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

def generate_quality_panel(gff_path, output_png="VALENCIA_Quality_Report_Definitivo.png"):
    """
    Genera el panel de calidad VALENCIA de 3 gráficos a 900 DPI.
    Soluciona la alineación del histograma rosa y usa un color de alto contraste.
    """
    print(f"📊 Generando Panel Definitivo (Simetría + Alto Contraste) desde: {gff_path}")
    data = []
    try:
        # 1. Extracción de Métricas robusta con regex
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

        # --- CONFIGURACIÓN DE ESTILO Y RESOLUCIÓN ---
        # 900 DPI para calidad máxima en TFG
        plt.rcParams['savefig.dpi'] = 900 
        sns.set_theme(style="white")
        fig = plt.figure(figsize=(36, 12))
        
        # Layout: 3 columnas principales
        gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.35)

        # --- PLOT A: SIMETRÍA PERFECTA Y COLOR 'PLASMA' (Alto Contraste) ---
        # 5 filas para aislar la Colorbar y que no estire el histograma rosa
        inner_gs = gs[0].subgridspec(5, 4, hspace=0.0, wspace=0.0)
        
        # Gráfico Principal (Cuadro grande) - Filas 1 a 4
        ax_main = fig.add_subplot(inner_gs[1:4, :-1])
        # Histograma X (Verde) - Fila 0
        ax_hist_x = fig.add_subplot(inner_gs[0, :-1], sharex=ax_main)
        # Histograma Y (Rosa) - Filas 1 a 4 (ALINEADO CON AX_MAIN)
        ax_hist_y = fig.add_subplot(inner_gs[1:4, -1], sharey=ax_main)

        # DIBUJADO CON NUEVO COLOR 'PLASMA' (Máxima visibilidad de degradado)
        sc = ax_main.scatter(df['tx'], df['pr'], c=df['Delta'], 
                            s=4, cmap="plasma", vmin=0, vmax=1, alpha=0.7, rasterized=True)

        ax_hist_x.hist(df['tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1)
        ax_hist_y.hist(df['pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        
        # Limpieza de ticks interiores
        ax_hist_x.tick_params(labelbottom=False, bottom=False)
        ax_hist_y.tick_params(labelleft=False, left=False)
        
        ax_main.axvline(0.1, color='red', linestyle='--', alpha=0.5)
        ax_main.axhline(0.1, color='red', linestyle='--', alpha=0.5)
        ax_main.set_xlabel('lev_edit_distance transcripts', fontweight='bold')
        ax_main.set_ylabel('lev_edit_distance proteins', fontweight='bold')

        # UBICACIÓN INDEPENDIENTE DE LA COLORBAR: En la última fila libre (fila 4)
        cax_a = fig.add_subplot(inner_gs[4, :-1])
        fig.colorbar(sc, cax=cax_a, orientation='horizontal').set_label('Δ lev_edit_distance (Discrepancy)')

        # Leyenda derecha para Plot A
        leg_a = [Patch(facecolor='#45a049', label='AED transcripts'),
                 Patch(facecolor='#e91e63', label='AED proteins')]
        ax_hist_y.legend(handles=leg_a, loc='upper left', bbox_to_anchor=(1.1, 1.0), frameon=True)

        # --- PLOT B y C (Se mantienen igual) ---
        ax_corr = fig.add_subplot(gs[1])
        sns.scatterplot(data=df, x='cds', y='pr', alpha=0.15, s=5, color='#34495e', ax=ax_corr, rasterized=True)
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identity (X=Y)')
        ax_corr.fill_between([0, 1], [0, 1], [1, 1], color='red', alpha=0.04, label='Potential frameshift zone')
        
        ax_corr.set_title('B. Correlation: CDS vs Protein', fontsize=18, fontweight='bold')
        ax_corr.set_xlabel('AED CDS (Nucleotide level)', fontweight='bold'); ax_corr.set_ylabel('AED Proteins (Amino acid level)', fontweight='bold')
        ax_corr.legend(loc='upper left')

        ax_dist = fig.add_subplot(gs[2])
        sns.histplot(df['pr'] - df['cds'], bins=50, kde=True, color='#2E86C1', ax=ax_dist, edgecolor='white')
        ax_dist.set_title('C. Distribution of Editing Difference', fontsize=18, fontweight='bold')
        ax_dist.set_xlabel('Difference (Protein - CDS)', fontweight='bold')
        ax_dist.set_ylabel('Frequency (Number of Genes)', fontweight='bold')
        ax_dist.axvline(0, color='red', linestyle=':', label='Zero Difference')
        ax_dist.grid(True, linestyle=':', alpha=0.4)
        ax_dist.legend()

        fig.suptitle(f'VALENCIA Annotation quality analysis (n={len(df)})', fontsize=26, fontweight='bold', y=0.98)
        
        # GUARDADO FINAL (Asegura PNG a 900 DPI)
        plt.savefig(output_png, bbox_inches='tight')
        plt.close()
        print(f"✅ Master Quality Report generado con éxito a 900 DPI: {output_png}")

    except Exception as e:
        print(f"❌ Error: {e}")