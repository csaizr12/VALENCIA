import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib.patches import Patch
from pathlib import Path

# Updated signature to accept description and species
def generate_quality_panel(gff_path, output_folder, description=None, species="Unknown species"):
    data = []
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set default description if none is provided via terminal
    if not description:
        description = Path(gff_path).stem

        # --- 1. DATA EXTRACTION ---
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
            print(f"No data found in {gff_path}")
            return

        n_samples = len(df)
        sns.set_theme(style="whitegrid")

        # --- PANEL A: STRUCTURE VS FUNCTION ---
        fig_a = plt.figure(figsize=(14, 12))
        gs = fig_a.add_gridspec(7, 4, height_ratios=[1.2, 1, 1, 1, 1, 0.6, 0.2], hspace=0.15, wspace=0.15)
        
        ax_main = fig_a.add_subplot(gs[1:5, :-1])
        ax_hist_x = fig_a.add_subplot(gs[0, :-1], sharex=ax_main)
        ax_hist_y = fig_a.add_subplot(gs[1:5, -1], sharey=ax_main)
        
        sc = ax_main.scatter(df['tx'], df['pr'], c=df['Delta'], s=5, cmap="magma", vmin=0, vmax=1, alpha=0.7, rasterized=True)
        
        ax_hist_x.hist(df['tx'], bins=80, color='#45a049', edgecolor='black', linewidth=0.1)
        ax_hist_x.set_ylabel('Nb. transcripts', fontweight='bold', fontsize=10) 
        
        ax_hist_y.hist(df['pr'], bins=80, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)
        ax_hist_y.set_xlabel('Nb. proteins', fontweight='bold', fontsize=10)
        
        # --- DYNAMIC TITLES FOR PANEL A ---
        fig_a.suptitle(f"VALENCIA Quality Analysis: {species} (n={n_samples})", 
                       fontsize=18, fontweight='bold', y=0.98)
        
        plt.figtext(0.5, 0.94, description, ha='center', fontsize=14, style='italic', color='#2c3e50')
        plt.figtext(0.5, 0.91, "Transcript-protein edit distance correlation", ha='center', fontsize=12, color='gray')

        ax_main.set_xlabel('Lev_edit_distance transcripts', fontweight='bold', fontsize=12)
        ax_main.set_ylabel('Lev_edit_distance proteins', fontweight='bold', fontsize=12)
        
        ax_hist_x.tick_params(labelbottom=False, bottom=False)
        ax_hist_y.tick_params(labelleft=False, left=False)
        
        cax = fig_a.add_subplot(gs[6, :-1])
        fig_a.colorbar(sc, cax=cax, orientation='horizontal', label='Δ Lev_edit_distance')
        
        leg_elements = [Patch(facecolor='#45a049', label='Transcripts'), Patch(facecolor='#e91e63', label='Proteins')]
        ax_hist_y.legend(handles=leg_elements, loc='upper left', bbox_to_anchor=(1.05, 1.2), frameon=True)
        
        plt.savefig(output_dir / f"VALENCIA_quality_correlation_{description.replace(' ', '_')}.svg", format='svg', bbox_inches='tight')
        plt.close(fig_a)
       # --- PANEL B: CORRELATION (CDS vs Protein) ---
        plt.figure(figsize=(10, 10))
        ax_corr = sns.scatterplot(data=df, x='cds', y='pr', alpha=0.2, s=15, color='#34495e', rasterized=True)
        ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identity (X=Y)')
        
        # dynamic titles for Panel B
        plt.suptitle(f"VALENCIA Quality Analysis: {species} (n={n_samples})", 
                     fontweight='bold', fontsize=15, y=0.98)
        plt.title(f"{description}\nCorrelation: CDS vs Protein", 
                  fontsize=12, pad=15, style='italic', color='#2c3e50')
        
        plt.xlabel('Lev_edit_distance CDS (Nucleotides)', fontweight='bold')
        plt.ylabel('Lev_edit_distance proteins (Amino acids)', fontweight='bold')
        plt.legend()
        
        # Save B with dynamic name
        plt.savefig(output_dir / f"VALENCIA_CDS_protein_correlation_{description.replace(' ', '_')}.svg", 
                    format='svg', bbox_inches='tight')
        plt.close()

        # --- PANEL C: DISTRIBUTION (Absolute difference) ---
        plt.figure(figsize=(10, 10))
        absolute_diff = df['diff'].abs()
        
        ax_dist = sns.histplot(absolute_diff, bins=100, kde=True, color='#2E86C1', 
                               edgecolor='white', rasterized=True)
        
        ax_dist.set_xlim(0, 0.15) 
        ax_dist.axvline(0, color='red', linestyle='--', linewidth=2, label='Match (Diff = 0)')
        
        # dynamic titles for Panel C
        plt.suptitle(f"VALENCIA Quality Analysis: {species} (n={n_samples})", 
                     fontweight='bold', fontsize=15, y=0.98)
        plt.title(f"{description}\nDistribution of absolute editing difference", 
                  fontsize=12, pad=15, style='italic', color='#2c3e50')
        
        plt.xlabel('|Protein dist. - CDS dist.| (Absolute difference)', fontweight='bold')
        plt.ylabel('Number of transcripts', fontweight='bold')
        plt.legend()
        
        # Save C with dynamic name
        plt.savefig(output_dir / f"VALENCIA_edit_distance_distribution_{description.replace(' ', '_')}.svg", 
                    format='svg', bbox_inches='tight', dpi=300)
        plt.close()