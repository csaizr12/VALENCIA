import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import sys
import numpy as np
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

def extract_metrics(gff_path):
    """Extracts AED metrics from GFF mRNA features."""
    data = []
    print(f"🧬 Processing GFF file: {gff_path}")
    try:
        with open(gff_path, 'r', encoding='latin-1') as f:
            for line in f:
                if line.startswith('#') or '\t' not in line: continue
                fields = line.split('\t')
                if fields[2].lower() in ['mrna', 'transcript']:
                    attr = fields[8]
                    # Case-insensitive search for distance attributes
                    tx = re.search(r'transcripts_edit_distance=([0-9.]+)', attr, re.I)
                    pr = re.search(r'proteins_edit_distance=([0-9.]+)', attr, re.I)
                    cds = re.search(r'cds_edit_distance=([0-9.]+)', attr, re.I)
                    
                    if tx and pr and cds:
                        t_v, p_v, c_v = float(tx.group(1)), float(pr.group(1)), float(cds.group(1))
                        # Delta: Absolute difference between Nucleotide (CDS) and Amino Acid (Protein) distances
                        data.append({
                            'AED_tx': t_v, 
                            'AED_pr': p_v, 
                            'AED_cds': c_v, 
                            'Delta': abs(p_v - c_v)
                        })
    except Exception as e:
        print(f"❌ Error: {e}")
    return pd.DataFrame(data)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 script.py <input.gff>")
        sys.exit(1)

    df = extract_metrics(sys.argv[1])
    if df.empty:
        print("⚠️ No valid AED metrics found in GFF.")
        return

    # Set professional style
    sns.set_theme(style="white")
    fig = plt.figure(figsize=(26, 12))
    # Layout: 2 main columns (A and B)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.4)

    # --- 1. PLOT A: STRUCTURE VS FUNCTION (HEXBIN + MARGINALS) ---
    inner_gs = gs[0].subgridspec(4, 4, hspace=0.1, wspace=0.1)
    
    ax_main = fig.add_subplot(inner_gs[1:, :-1])
    ax_hist_x = fig.add_subplot(inner_gs[0, :-1], sharex=ax_main)
    ax_hist_y = fig.add_subplot(inner_gs[1:, -1], sharey=ax_main)

    # Hexbin showing maximum discrepancy (Delta)
    hb = ax_main.hexbin(df['AED_tx'], df['AED_pr'], C=df['Delta'], 
                        reduce_C_function=np.max, gridsize=40, 
                        cmap="YlOrRd", mincnt=1, edgecolors='none', vmin=0, vmax=1)

    # Marginal Histograms (Green for X, Crimson for Y)
    ax_hist_x.hist(df['AED_tx'], bins=60, color='#45a049', edgecolor='black', linewidth=0.2)
    ax_hist_y.hist(df['AED_pr'], bins=60, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.2)

    # Marginal labels
    ax_hist_x.set_ylabel('Nb. Transcripts', fontsize=9)
    ax_hist_y.set_xlabel('Nb. Transcripts', fontsize=9)
    ax_hist_x.tick_params(labelbottom=False)
    ax_hist_y.tick_params(labelleft=False)

    # Quality Thresholds (0.3 and 0.1)
    ax_main.axvline(0.3, color='red', linestyle='--', alpha=0.6)
    ax_main.axhline(0.1, color='red', linestyle='--', alpha=0.6)
    ax_main.set_xlim(-0.02, 1.02); ax_main.set_ylim(-0.02, 1.02)
    ax_main.set_xlabel('AED Transcripts (Structure match)', fontweight='bold', fontsize=11)
    ax_main.set_ylabel('AED Proteins (Function match)', fontweight='bold', fontsize=11)

    # REFINED COLORBAR (Positioned exactly below the main plot)
    divider = make_axes_locatable(ax_main)
    cax = divider.append_axes("bottom", size="4%", pad=0.7)
    cbar = fig.colorbar(hb, cax=cax, orientation='horizontal')
    cbar.set_label('Δ AED (Max discrepancy CDS vs Protein | 1.0 = Frameshift)', fontweight='bold', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # EXTERNAL LEGEND (Positioned to the right of the marginal Y histogram)
    leg_elements = [Patch(facecolor='#45a049', label='AED transcripts'),
                    Patch(facecolor='#e91e63', label='AED proteins')]
    ax_hist_y.legend(handles=leg_elements, loc='upper left', bbox_to_anchor=(1.05, 1.0), frameon=True)

    # --- 2. PLOT B: INTERNAL CONSISTENCY (CDS VS PROTEIN CORRELATION) ---
    ax_corr = fig.add_subplot(gs[1])
    sns.scatterplot(data=df, x='AED_cds', y='AED_pr', alpha=0.2, s=10, color='#34495e', ax=ax_corr, edgecolor=None)
    
    # Identity Line (X=Y)
    ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', linewidth=2, label='Internal consistency (X=Y)')
    # Frameshift Zone Highlight
    ax_corr.fill_between([0, 1], [0, 1], [1, 1], color='red', alpha=0.05, label='Potential frameshift zone')
    
    ax_corr.set_xlabel('AED CDS (Nucleotide level)', fontweight='bold', fontsize=11)
    ax_corr.set_ylabel('AED Proteins (Amino acid level)', fontweight='bold', fontsize=11)
    ax_corr.set_title('Correlation: CDS edit distance vs Protein edit distance', fontsize=14, fontweight='bold')
    ax_corr.set_xlim(-0.02, 1.02); ax_corr.set_ylim(-0.02, 1.02)
    ax_corr.legend(loc='upper left', frameon=True)
    ax_corr.grid(True, linestyle=':', alpha=0.5)

    # Main Title
    fig.suptitle(f'VALENCIA Annotation quality panel (n={len(df)})', fontsize=20, fontweight='bold', y=0.98)
    
    # Save the figure
    output_file = "AED_VALENCIA_quality_panel_.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()