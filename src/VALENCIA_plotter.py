import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np  
from matplotlib.patches import Patch
from pathlib import Path

def generate_quality_panel(gff_path, output_folder, description=None, species=None):
    data = []
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. DATA EXTRACTION (Now safely outside the IF) ---
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

    # Clean infinite and null values that cause matplotlib broadcasting errors
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['tx', 'pr', 'cds'])
    if df.empty:
        print(f"No valid data left in {gff_path} after removing NaNs/Infs")
        return

    n_samples = len(df)
    sns.set_theme(style="whitegrid")

    # Determine internal description for filenames if None
    file_desc = description if description else Path(gff_path).stem

    # --- HELPER: DYNAMIC TITLE LOGIC ---
    def set_dynamic_titles(fig_or_plt, n, spec, desc, metric_text):
        # Line 1: Main Title
        main_t = f"VALENCIA annotation quality analysis (n={n})"
        if hasattr(fig_or_plt, 'suptitle'):
            fig_or_plt.suptitle(main_t, fontsize=16, fontweight='bold', y=0.98)
        else:
            plt.suptitle(main_t, fontsize=16, fontweight='bold', y=0.98)

        # Line 2: Combined Info (Species | Description)
        info_parts = [i for i in [spec, desc] if i]
        if info_parts:
            combined = " | ".join(info_parts)
            plt.figtext(0.5, 0.94, combined, ha='center', fontsize=13, style='italic', color='#2c3e50')
            y_metric = 0.91
        else:
            # If no info, move metric up to match image_929c78.png
            y_metric = 0.94

        # Line 3: Metric
        plt.figtext(0.5, y_metric, metric_text, ha='center', fontsize=11, color='gray')

    # --- PANEL A: STRUCTURE VS FUNCTION ---
    fig_a = plt.figure(figsize=(14, 12))
    gs = fig_a.add_gridspec(7, 4, height_ratios=[1.2, 1, 1, 1, 1, 0.6, 0.2], hspace=0.15, wspace=0.15)
    
    ax_main = fig_a.add_subplot(gs[1:5, :-1])
    ax_hist_x = fig_a.add_subplot(gs[0, :-1], sharex=ax_main)
    ax_hist_y = fig_a.add_subplot(gs[1:5, -1], sharey=ax_main)
    
    # Create explicit bin edges to bypass the uniform-bin bug in Python 3.14t
    bins_x = np.linspace(df['tx'].min(), df['tx'].max(), 81) if df['tx'].min() != df['tx'].max() else 80
    bins_y = np.linspace(df['pr'].min(), df['pr'].max(), 81) if df['pr'].min() != df['pr'].max() else 80

    sc = ax_main.scatter(df['tx'], df['pr'], c=df['Delta'], s=5, cmap="magma", vmin=0, vmax=1, alpha=0.7, rasterized=True)
    ax_hist_x.hist(df['tx'], bins=bins_x, color='#45a049', edgecolor='black', linewidth=0.1)
    ax_hist_y.hist(df['pr'], bins=bins_y, color='#e91e63', orientation='horizontal', edgecolor='black', linewidth=0.1)

    set_dynamic_titles(fig_a, n_samples, species, description, "Transcript-protein edit distance correlation")

    ax_main.set_xlabel('Lev_edit_distance transcripts', fontweight='bold', fontsize=12)
    ax_main.set_ylabel('Lev_edit_distance proteins', fontweight='bold', fontsize=12)
    ax_hist_x.tick_params(labelbottom=False, bottom=False)
    ax_hist_y.tick_params(labelleft=False, left=False)
    
    cax = fig_a.add_subplot(gs[6, :-1])
    fig_a.colorbar(sc, cax=cax, orientation='horizontal', label='Δ Lev_edit_distance')
    
    leg_elements = [Patch(facecolor='#45a049', label='Transcripts'), Patch(facecolor='#e91e63', label='Proteins')]
    ax_hist_y.legend(handles=leg_elements, loc='upper left', bbox_to_anchor=(1.05, 1.2), frameon=True)
    
    plt.savefig(output_dir / f"VALENCIA_quality_correlation_{file_desc.replace(' ', '_')}.svg", format='svg', bbox_inches='tight')
    plt.close(fig_a)

    # --- PANEL B: CORRELATION ---
    plt.figure(figsize=(10, 10))
    ax_corr = sns.scatterplot(data=df, x='cds', y='pr', alpha=0.2, s=15, color='#34495e', rasterized=True)
    ax_corr.plot([0, 1], [0, 1], color='red', linestyle='--', label='Identity (X=Y)')
    
    set_dynamic_titles(plt, n_samples, species, description, "Correlation: CDS vs Protein")
    plt.subplots_adjust(top=0.88)
    
    plt.xlabel('Lev_edit_distance CDS (Nucleotides)', fontweight='bold')
    plt.ylabel('Lev_edit_distance proteins (Amino acids)', fontweight='bold')
    plt.legend()
    plt.savefig(output_dir / f"VALENCIA_CDS_protein_correlation_{file_desc.replace(' ', '_')}.svg", format='svg', bbox_inches='tight')
    plt.close()

    # --- PANEL C: DISTRIBUTION ---
    plt.figure(figsize=(10, 10))
    absolute_diff = df['diff'].abs()
    ax_dist = sns.histplot(absolute_diff, bins=100, kde=True, color='#2E86C1', edgecolor='white', rasterized=True)
    ax_dist.set_xlim(0, 0.15) 
    ax_dist.axvline(0, color='red', linestyle='--', linewidth=2, label='Match (Diff = 0)')
    
    set_dynamic_titles(plt, n_samples, species, description, "Distribution of absolute editing difference")
    plt.subplots_adjust(top=0.88)
    
    plt.xlabel('|Protein dist. - CDS dist.| (Absolute difference)', fontweight='bold')
    plt.ylabel('Number of transcripts', fontweight='bold')
    plt.legend()
    plt.savefig(output_dir / f"VALENCIA_edit_distance_distribution_{file_desc.replace(' ', '_')}.svg", format='svg', bbox_inches='tight', dpi=300)
    plt.close()