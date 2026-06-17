import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

def generate_class_code_plots(gff_files):
    output_dir = "dataset_test/results_plots"
    os.makedirs(output_dir, exist_ok=True)

    tc_pattern = re.compile(r'transcripts_class_code[ =]"?([^";\s]+)"?', re.IGNORECASE)
    pc_pattern = re.compile(r'proteins_class_code[ =]"?([^";\s]+)"?', re.IGNORECASE)
    fallback_pattern = re.compile(r'class_code[ =]"?([^";\s]+)"?', re.IGNORECASE)

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
        'NA': 'Without evidence (NA)'
    }

    data = {
        'transcript': {},
        'protein': {}
    }

    for full_path in gff_files:
        if not os.path.exists(full_path):
            continue

        pipeline_name = os.path.basename(os.path.dirname(full_path))
        species_raw = os.path.basename(os.path.dirname(os.path.dirname(full_path)))
        
        species_name = species_raw.replace('_dataset_test', '').replace('_datset_test', '').replace('_', ' ')
        if not species_name:
            species_name = "Unknown_Species"

        for mode in ['transcript', 'protein']:
            if species_name not in data[mode]:
                data[mode][species_name] = {}
            if pipeline_name not in data[mode][species_name]:
                data[mode][species_name][pipeline_name] = {}

        try:
            with open(full_path, 'r', errors='ignore') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue

                    attr = parts[8]

                    tc_match = tc_pattern.search(attr)
                    if tc_match:
                        c = tc_match.group(1).strip()
                        data['transcript'][species_name][pipeline_name][c] = data['transcript'][species_name][pipeline_name].get(c, 0) + 1
                    else:
                        fb_match = fallback_pattern.search(attr)
                        if fb_match and '\tmRNA\t' in line:
                            c = fb_match.group(1).strip()
                            data['transcript'][species_name][pipeline_name][c] = data['transcript'][species_name][pipeline_name].get(c, 0) + 1

                    pc_match = pc_pattern.search(attr)
                    if pc_match:
                        c = pc_match.group(1).strip()
                        data['protein'][species_name][pipeline_name][c] = data['protein'][species_name][pipeline_name].get(c, 0) + 1

        except Exception as e:
            print(f"  [Error] Processing {full_path}: {e}")

    distinct_colors = [
        '#0f4c81', '#ffb347', '#1d8348', '#f1948a', '#6c3483', 
        '#bb8fce', '#117a65', '#73c6b6', '#9a7d0a', '#f7dc6f',
        '#5dade2', '#2e4053', '#aeb6bf', '#edbb99', '#f5b041'
    ]

    for mode in ['transcript', 'protein']:
        for species_name, pipelines_data in data[mode].items():
            flat_data = []
            for pipe, counts in pipelines_data.items():
                for code, count in counts.items():
                    clean_code = legend_mapping.get(code, code)
                    flat_data.append({'Pipeline': pipe, 'Class Code': clean_code, 'Count': count})

            if not flat_data:
                continue

            df = pd.DataFrame(flat_data)
            df_pivot = df.pivot(index='Pipeline', columns='Class Code', values='Count').fillna(0)
            
            def custom_sort(pipe_name):
                base = pipe_name.replace('test_UTR_', 'test_')
                is_utr = 1 if 'UTR' in pipe_name else 0
                return (base, is_utr)

            ordered_pipelines = sorted(df_pivot.index, key=custom_sort)
            df_pivot = df_pivot.reindex(ordered_pipelines)
            
            pipeline_totals = df_pivot.sum(axis=1).astype(int)
            df_pivot.index = [f"{pipe} (N={pipeline_totals[pipe]})" for pipe in df_pivot.index]
            
            df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0) * 100

            fig, ax = plt.subplots(figsize=(14, 12))
            df_pivot.plot(kind='barh', stacked=True, color=distinct_colors, ax=ax)
            
            title_text = f"Distribution of {mode} class codes ({species_name})"
            x_label_text = f"Proportion of total {mode}s (%)"
            legend_title_text = f"VALENCIA {mode} class codes"

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
            output_filename = f"{output_dir}/{clean_name}_{mode}_class_codes_comparison.svg"
            plt.savefig(output_filename, format='svg', dpi=300)
            plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 script.py <file1.gff3> <file2.gff3> ...")
        sys.exit(1)

    input_files = sys.argv[1:]
    generate_class_code_plots(input_files)