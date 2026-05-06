import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib.patches import Patch
from pathlib import Path

def generate_quality_panel(gff_path, output_folder):
    data = []
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        # --- 1. EXTRACCIÓN DE DATOS ---
        fp_pattern = r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'
        with open(gff_path, 'r', encoding='latin-1') as f:
            for line in f:
                if line.startswith('#') or '\t' not in line: continue
                fields = line.split('\t')
                if fields[2].lower() in ['mrna', 'transcript']:
                    attr = fields[8]
                    tx = re.search(r'transcripts_edit_distance=' + fp_pattern, attr, re.I)
                    pr = re.search(r'proteins_edit_distance=' + fp_pattern, attr, re.I)
                    cds = re.search(r'cds_edit_distance=' + fp_pattern, attr, re.I)
                    
                    if tx and pr and cds:
                        data.append({
                            'tx': float(tx.group(1)), 
                            'pr': float(pr.group(1)), 
                            'cds': float(cds.group(1)), 
                            'diff': float(pr.group(1)) - float(cds.group(1)),
                            'Delta': abs(float(pr.group(1)) - float(cds.group(1)))
                        })
        df = pd.DataFrame(data)
        if df.empty:
            print(f"No se encontraron datos en {gff_path}")
            return

        n_samples = len(df)
        main_title = f"VALENCIA annotation quality analysis (n={n_samples})"
        sns.set_theme(style="whitegrid")

        # --- PANEL A: ESTRUCTURA VS FUNCIÓN ---
        fig_a = plt.figure(figsize=(14, 12))
        gs = fig_a.add_gridspec(7, 4, height_ratios=[1.2, 1, 1, 1, 1, 0.6, 0.2], hspace=0.1, wspace=0.1)
        ax_main = fig_a.add_subplot(gs[1:5, :-1])
        ax_hist_x = fig_a.add_subplot(gs[0, :-1], sharex=ax_main)
        ax_hist_y = fig_a.add_subplot(gs[1:5, -1], sharey=ax_main)
        
        sc = ax_main.scatter(df['tx'], df['pr'], c=df['Delta'], s=5, cmap="magma", vmin=0, vmax=1, alpha=0.7)
        ax_hist_x.hist(df['tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1)
        ax_hist_y.hist(df['pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        
        fig_a.suptitle(main_title, fontsize=20, fontweight='bold')
        ax_main.set_xlabel('Lev_edit_distance transcripts', fontweight='bold', fontsize=12)
        ax_main.set_ylabel('Lev_edit_distance proteins', fontweight='bold', fontsize=12)
        
        cax = fig_a.add_subplot(gs[6, :-1])
        fig_a.colorbar(sc, cax=cax, orientation='horizontal', label='Δ Lev_edit_distance')
        
        plt.savefig(output_dir / "VALENCIA_panel_A_structure.svg", format='svg', bbox_inches='tight')
        plt.close(fig_a)

        # --- PANEL B: CORRELACIÓN ---
        plt.figure(figsize=(10, 10))
        ax_corr = sns.scatterplot(data=df, x='cds', y='pr', alpha=0.2, s=15, color='#34495e')
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identity (X=Y)')
        
        plt.title(f"{main_title}\nCorrelation: CDS vs protein", fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Lev_edit_distance CDS (Nucleotides)', fontweight='bold')
        plt.ylabel('Lev_edit_distance proteins (Amino acids)', fontweight='bold')
        plt.legend()
        
        plt.savefig(output_dir / "VALENCIA_panel_B_correlation.svg", format='svg', bbox_inches='tight')
        plt.close()

        # --- PANEL C: DISTRIBUCIÓN ---
        plt.figure(figsize=(10, 10))
        ax_dist = sns.histplot(df['diff'], bins=100, kde=True, color='#2E86C1', edgecolor='white')
        
        # Ajuste de zoom para ver la distribución central
        ax_dist.set_xlim(-0.05, 0.15) 
        ax_dist.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero difference')
        
        plt.title(f"{main_title}\nDistribution of editing difference", fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Difference (Protein dist. - CDS dist.)', fontweight='bold')
        plt.ylabel('Number of transcripts', fontweight='bold')
        plt.legend()
        
        plt.savefig(output_dir / "VALENCIA_panel_C_distribution.svg", format='svg', bbox_inches='tight')
        plt.close()

    except Exception as e:
        print(f"Error crítico: {e}")