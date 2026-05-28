import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def extract_chrom_stats_for_distance_1(base_path, annotation_mapping, species_name, output_dir, mode='transcript'):
    """
    Parses GFF3 files to compute absolute counts and relative percentages 
    of distance="1" elements per chromosome, exporting a CSV and an SVG heatmap.
    Appends the total element count to the pipeline names on the Y-axis.
    """
    if not os.path.exists(base_path):
        print(f"Error: Path '{base_path}' not found.")
        return

    print(f"\n Analyzing Chromosomes for {species_name} ({mode.upper()} - Distance = 1 Analysis)...")
    
    distance_pattern = re.compile(r'distance[ =]"?1"?')
    gene_id_pattern = re.compile(r'gene_id[ =]"?([^";\s]+)"?')
    species_results = []
    pipeline_totals = {}

    for technical_folder, official_name in annotation_mapping.items():
        folder_path = os.path.join(base_path, technical_folder)
        if not os.path.exists(folder_path):
            continue
            
        gff_files = [f for f in os.listdir(folder_path) if f.endswith('.gff3')]
        
        for file_name in gff_files:
            full_path = os.path.join(folder_path, file_name)
            chrom_total_elements = {}
            chrom_distance1_elements = {}
            chrom_genes = {}
            has_distance_data = False
            
            try:
                with open(full_path, 'r', errors='ignore') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        parts = line.split('\t')
                        if len(parts) < 9:
                            continue
                        
                        if mode == 'protein' and parts[2].lower() not in ['cds', 'protein', 'mrna']:
                            continue
                        elif mode == 'transcript' and parts[2].lower() not in ['transcript', 'mrna', 'exon']:
                            continue
                            
                        chrom = parts[0].strip()
                        attributes = parts[8]
                        
                        chrom_total_elements[chrom] = chrom_total_elements.get(chrom, 0) + 1
                        
                        if distance_pattern.search(attributes):
                            has_distance_data = True
                            chrom_distance1_elements[chrom] = chrom_distance1_elements.get(chrom, 0) + 1
                            
                            gene_match = gene_id_pattern.search(attributes)
                            if gene_match:
                                gene_id = gene_match.group(1)
                                if chrom not in chrom_genes:
                                    chrom_genes[chrom] = set()
                                chrom_genes[chrom].add(gene_id)
                
                if has_distance_data:
                    # Calculate genome-wide absolute total for this specific pipeline
                    pipeline_totals[official_name] = sum(chrom_total_elements.values())
                    
                    for chrom in sorted(chrom_total_elements.keys()):
                        total_el_chrom = chrom_total_elements[chrom]
                        dist1_el_chrom = chrom_distance1_elements.get(chrom, 0)
                        
                        percentage_dist1 = (dist1_el_chrom / total_el_chrom) * 100 if total_el_chrom > 0 else 0.0
                        total_genes = len(chrom_genes.get(chrom, set()))
                        
                        species_results.append({
                            'Pipeline': official_name,
                            'Chromosome': chrom,
                            f'Total {mode.capitalize()}s (Dist 1)': dist1_el_chrom,
                            f'Total Chromosome {mode.capitalize()}s': total_el_chrom,
                            'Percentage (Dist 1 %)': round(percentage_dist1, 2),
                            'Total Unique Genes': total_genes
                        })
                    break
                    
            except Exception as e:
                print(f"Error processing {file_name}: {e}")

    if species_results:
        df = pd.DataFrame(species_results)
        clean_name = species_name.lower().replace(' ', '_')
        
        if mode == 'protein':
            cbar_label = '% of Proteins with distance = 1 (Relative to chromosome total)'
            plot_title = f"Chromosomal Distribution of distance=1 Proteins (%)\n({species_name})"
            output_csv = f"{output_dir}/{clean_name}_chrom_protein_distance1_summary.csv"
            output_svg = f"{output_dir}/{clean_name}_chrom_protein_distance1_percentage_distribution.svg"
        else:
            cbar_label = '% of Transcripts with distance = 1 (Relative to chromosome total)'
            plot_title = f"Chromosomal Distribution of distance=1 Transcripts (%)\n({species_name})"
            output_csv = f"{output_dir}/{clean_name}_chrom_transcript_distance1_summary.csv"
            output_svg = f"{output_dir}/{clean_name}_chrom_transcript_distance1_percentage_distribution.svg"

        df.to_csv(output_csv, index=False)
        print(f"Table saved to: {output_csv}")
        
        try:
            # SOLUCIÓN ERROR 1: Se usa pivot_table con aggfunc='mean' para evitar colapsos por claves duplicadas (como GEMOMA)
            df_plot = df.pivot_table(index='Pipeline', columns='Chromosome', values='Percentage (Dist 1 %)', aggfunc='mean').fillna(0)
            
            ordered_pipelines = [name for name in annotation_mapping.values() if name in df_plot.index]
            seen = set()
            unique_pipelines = [x for x in ordered_pipelines if not (x in seen or seen.add(x))]
            
            # Reindex utilizing the predefined strict pipeline order
            df_plot = df_plot.reindex(unique_pipelines[::-1])
            
            # Map new dynamic labels containing the absolute total counts: Pipeline (N=X)
            new_labels = [f"{idx} (N={pipeline_totals.get(idx, 0)})" for idx in df_plot.index]
            
            fig, ax = plt.subplots(figsize=(13, 9))
            
            # SOLUCIÓN ERROR 2: Eliminados los argumentos inválidos 'edgecolors' y 'linewidths' de matshow
            cax = ax.matshow(df_plot, cmap='twilight_shifted', aspect='auto')

            # Dibujar las líneas de rejilla de forma nativa y compatible con matshow
            ax.set_xticks(np.arange(-0.5, len(df_plot.columns), 1), minor=True)
            ax.set_yticks(np.arange(-0.5, len(df_plot.index), 1), minor=True)
            ax.grid(which='minor', color='#b0b0b0', linestyle='-', linewidth=0.4)

            cbar = fig.colorbar(cax, pad=0.02)
            cbar.set_label(cbar_label, fontsize=11, labelpad=10)
            cbar.ax.tick_params(labelsize=9)
            
            ax.set_xticks(range(len(df_plot.columns)))
            ax.set_xticklabels(df_plot.columns, rotation=45, ha='left', fontsize=10)
            
            # Apply the new dynamically computed labels to the Y axis
            ax.set_yticks(range(len(df_plot.index)))
            ax.set_yticklabels(new_labels, fontsize=10)
            ax.xaxis.set_ticks_position('bottom')
            
            ax.set_title(plot_title, fontsize=13, fontweight='bold', pad=20)
            ax.set_xlabel('Chromosomes', fontsize=11, labelpad=10)
            ax.set_ylabel('Annotation pipeline', fontsize=11, labelpad=10)
            
            plt.tight_layout()
            plt.savefig(output_svg, format='svg', dpi=300)
            plt.close()
            print(f"SVG percentage plot successfully saved to: {output_svg}")
            
        except Exception as graph_err:
            print(f"Error plotting SVG: {graph_err}")
    else:
        print(f" No records found for {species_name} in {mode} mode.")

if __name__ == "__main__":
    
    output_dir = "results_plots"
    os.makedirs(output_dir, exist_ok=True)
    input_base_dir = "dataset_test"

    species_list = [
        {"folder": "Arabidopsis_thaliana_dataset_test", "name": "Arabidopsis thaliana", "ref": "ARAPORT11"},
        {"folder": "Nicotiana_benthamiana_dataset_test", "name": "Nicotiana benthamiana", "ref": "GuoTTv1"},
        {"folder": "Oryza_sativa_dataset_test", "name": "Oryza sativa", "ref": "AGIS1v1"},
        {"folder": "Solanum_lycopersicum_datset_test", "name": "Solanum lycopersicum", "ref": "ITAG5.0"}
    ]

    for species in species_list:
        path = os.path.join(input_base_dir, species["folder"])
        
        if not os.path.exists(path):
            print(f"Warning: Data path '{path}' does not exist. Skipping {species['name']}.")
            continue
            
        mapping = {
            'test_BRAT': 'BRAKER1', 'test_UTR_BRAT': 'BRAKER1 + UTR',
            'test_BRAP': 'BRAKER2', 'test_UTR_BRAP': 'BRAKER2 + UTR',
            'test_BRATP': 'BRAKER3', 'test_UTR_BRATP': 'BRAKER3 + UTR',
            'test_EGX': 'Eukaryotic Genome Annotation Pipeline', 'test_UTR_EGX': 'Eukaryotic Genome Annotation Pipeline + UTR',
            'test_EVI': 'EviAnn', 'test_UTR_EVI': 'EviAnn + UTR',
            'test_EVO': 'ANNEVO', 'test_UTR_EVO': 'ANNEVO + UTR',
            'test_GEM': 'GEMOMA', 'test_UTR_GEM': 'GEMOMA',
            'test_HEL': 'HELIXER', 'test_UTR_HEL': 'HELIXER + UTR',
            'test_MAK': 'MAKER', 'test_UTR_MAK': 'MAKER + UTR',
            'test_REF': f'{species["ref"]}', 'test_UTR_REF': f'{species["ref"]} + UTR'
        }
        
        # 1. Run Analysis for Structural Transcripts
        extract_chrom_stats_for_distance_1(path, mapping, species["name"], output_dir, mode='transcript')
        
        # 2. Run Analysis for Functional Proteins
        extract_chrom_stats_for_distance_1(path, mapping, species["name"], output_dir, mode='protein')