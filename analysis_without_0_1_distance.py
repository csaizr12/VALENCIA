import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

def get_attr(attr_line, key):
    pattern = rf'{key}\s*[=:]\s*([^;]+)'
    match = re.search(pattern, attr_line, re.IGNORECASE)
    if match:
        val = match.group(1).strip().replace('"', '')
        return val if val.upper() != "NA" else "NA"
    return "NA"

def generate_filtered_class_code_plots(gff_files):
    output_dir = "results_plots"
    os.makedirs(output_dir, exist_ok=True)

    legend_mapping = {
        'complete': 'Complete (=)',
        'SubsequencesTarget': 'Subsequences target (c)',
        'SubsequencesReferences': 'Subsequences references (k)',
        'TotalIntronsRetention': 'Total introns retention (m)',
        'PartialIntronRetention': 'Partial intron retention (n)',
        'PotentialIsoform': 'Potential isoform (j)',
        'PartialExonOverlap': 'Partial exon overlap (o)',
        'RetainedIntronSingleExon': 'Retained intron single exon (e)',
        'Unknown': 'Unknown intergenic (u)',
    }

    data = {}

    for full_path in gff_files:
        if not os.path.exists(full_path):
            continue

        pipeline_name = os.path.basename(os.path.dirname(full_path))
        species_raw = os.path.basename(os.path.dirname(os.path.dirname(full_path)))
        
        species_name = species_raw.replace('_dataset_test', '').replace('_datset_test', '').replace('_', ' ')
        if not species_name:
            species_name = "Unknown_Species"
            
        if species_name not in data:
            data[species_name] = {}
        if pipeline_name not in data[species_name]:
            data[species_name][pipeline_name] = {}

        try:
            with open(full_path, 'r', errors='ignore') as f:
                for line in f:
                    if '\tmRNA\t' not in line: 
                        continue

                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue

                    attr = parts[8]
                    dist_str = get_attr(attr, 'transcripts_edit_distance')
                    
                    try:
                        dist = float(dist_str)
                    except ValueError:
                        continue

                    # Filter: Only keep transcripts strictly between 0 and 1
                    if 0.0 < dist < 1.0:
                        c = get_attr(attr, 'transcripts_class_code')
                        if c == "NA":
                            c = get_attr(attr, 'class_code')
                        
                        data[species_name][pipeline_name][c] = data[species_name][pipeline_name].get(c, 0) + 1

        except Exception as e:
            print(f"  [Error] Processing {full_path}: {e}")

    distinct_colors = [
        '#0f4c81', '#ffb347', '#1d8348', '#f1948a', '#6c3483', 
        '#bb8fce', '#117a65', '#73c6b6', '#9a7d0a', '#f7dc6f',
        '#5dade2', '#2e4053', '#aeb6bf', '#edbb99', '#f5b041'
    ]

    for species_name, pipelines_data in data.items():
        flat_data = []
        for pipe, counts in pipelines_data.items():
            for code, count in counts.items():
                clean_code = legend_mapping.get(code, code)
                flat_data.append({'Pipeline': pipe, 'Class Code': clean_code, 'Count': count})

        if not flat_data:
            print(f"  -> No intermediate data found for {species_name}. Skipping plot.")
            continue

        df = pd.DataFrame(flat_data)
        df_pivot = df.pivot(index='Pipeline', columns='Class Code', values='Count').fillna(0)
        
        def custom_sort(pipe_name):
            base = pipe_name.replace('test_UTR_', 'test_')
            is_utr = 1 if 'UTR' in pipe_name else 0
            return (base, is_utr)

        ordered_pipelines = sorted(df_pivot.index, key=custom_sort)
        df_pivot = df_pivot.reindex(ordered_pipelines)
        
        totals = df_pivot.sum(axis=1)
        
        df_pivot = df_pivot.div(totals, axis=0) * 100
        
        df_pivot.index = [f"{pipe} (n={int(total)})" for pipe, total in zip(df_pivot.index, totals)]

        fig, ax = plt.subplots(figsize=(14, 12))
        
        used_colors = distinct_colors[:len(df_pivot.columns)]
        df_pivot.plot(kind='barh', stacked=True, color=used_colors, ax=ax)
        
        title_text = f"Distribution of transcript class codes ({species_name})\nFiltered range: 0.0 < dist < 1.0"
        x_label_text = "Proportion of total transcripts (%)"
        legend_title_text = "VALENCIA transcript class codes"

        ax.set_title(title_text, fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel(x_label_text, fontsize=12, labelpad=10)
        ax.set_ylabel('Annotation pipeline', fontsize=12, labelpad=10)
        
        ax.set_xlim(0, 100)
        ax.set_xticks(np.arange(0, 101, 10))
        ax.set_xticklabels([str(x) for x in np.arange(0, 101, 10)], fontsize=9)
        
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.tick_params(which='minor', length=4, color='#7f8c8d')
        ax.tick_params(which='major', length=7, color='#2c3e50')
        
        ax.grid(axis='x', which='major', linestyle='--', alpha=0.6)
        ax.legend(title=legend_title_text, bbox_to_anchor=(1.02, 1), loc='upper left')
        
        plt.tight_layout()
        
        clean_name = species_name.lower().replace(' ', '_')
        output_filename = f"{output_dir}/{clean_name}_transcript_class_codes_mountain.svg"
        plt.savefig(output_filename, format='svg', dpi=300)
        plt.close()
        
        print(f"Generated plot: {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 script.py <file1.gff3> <file2.gff3> ...")
        sys.exit(1)

    input_files = sys.argv[1:]
    generate_filtered_class_code_plots(input_files)