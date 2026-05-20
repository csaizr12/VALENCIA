import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def extract_chrom_stats_for_distance_1(base_path, annotation_mapping, species_name, output_dir):
    """
    Parses GFF3 files to compute absolute counts and relative percentages 
    of distance="1" transcripts per chromosome, exporting a CSV and an SVG heatmap.
    """
    if not os.path.exists(base_path):
        print(f"Error: Path '{base_path}' not found.")
        return

    print(f"\n Analyzing Chromosomes for {species_name} (Distance = 1 Analysis)...")
    
    distance_pattern = re.compile(r'distance[ =]"?1"?')
    gene_id_pattern = re.compile(r'gene_id[ =]"?([^";\s]+)"?')
    species_results = []

    for technical_folder, official_name in annotation_mapping.items():
        folder_path = os.path.join(base_path, technical_folder)
        if not os.path.exists(folder_path):
            continue
            
        gff_files = [f for f in os.listdir(folder_path) if f.endswith('.gff3')]
        
        for file_name in gff_files:
            full_path = os.path.join(folder_path, file_name)
            chrom_total_transcripts = {}
            chrom_distance1_transcripts = {}
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
                        
                        chrom = parts[0].strip()
                        attributes = parts[8]
                        
                        chrom_total_transcripts[chrom] = chrom_total_transcripts.get(chrom, 0) + 1
                        
                        if distance_pattern.search(attributes):
                            has_distance_data = True
                            chrom_distance1_transcripts[chrom] = chrom_distance1_transcripts.get(chrom, 0) + 1
                            
                            gene_match = gene_id_pattern.search(attributes)
                            if gene_match:
                                gene_id = gene_match.group(1)
                                if chrom not in chrom_genes:
                                    chrom_genes[chrom] = set()
                                chrom_genes[chrom].add(gene_id)
                
                if has_distance_data:
                    for chrom in sorted(chrom_total_transcripts.keys()):
                        total_tx_chrom = chrom_total_transcripts[chrom]
                        dist1_tx_chrom = chrom_distance1_transcripts.get(chrom, 0)
                        
                        percentage_dist1 = (dist1_tx_chrom / total_tx_chrom) * 100 if total_tx_chrom > 0 else 0.0
                        total_genes = len(chrom_genes.get(chrom, set()))
                        
                        species_results.append({
                            'Pipeline': official_name,
                            'Chromosome': chrom,
                            'Total Transcripts (Dist 1)': dist1_tx_chrom,
                            'Total Chromosome Transcripts': total_tx_chrom,
                            'Percentage Transcripts (Dist 1 %)': round(percentage_dist1, 2),
                            'Total Unique Genes': total_genes
                        })
                    break
                    
            except Exception as e:
                print(f"Error processing {file_name}: {e}")

    if species_results:
        df = pd.DataFrame(species_results)
        
        # CSV Export Execution inside results_plots folder
        clean_name = species_name.lower().replace(' ', '_')
        output_csv = f"{output_dir}/{clean_name}_chrom_distance1_summary.csv"
        df.to_csv(output_csv, index=False)
        print(f"Table saved to: {output_csv}")
        
        # SVG Heatmap Generation inside results_plots folder
        try:
            df_plot = df.pivot(index='Pipeline', columns='Chromosome', values='Percentage Transcripts (Dist 1 %)').fillna(0)
            
            ordered_pipelines = [name for name in annotation_mapping.values() if name in df_plot.index]
            seen = set()
            unique_pipelines = [x for x in ordered_pipelines if not (x in seen or seen.add(x))]
            df_plot = df_plot.reindex(unique_pipelines[::-1])
            
            fig, ax = plt.subplots(figsize=(12, 8))
            cax = ax.matshow(df_plot, cmap='YlGnBu', aspect='auto')
            
            cbar = fig.colorbar(cax, pad=0.02)
            cbar.set_label('% of Transcripts with Distance = 1 (Relative to Chromosome Total)', fontsize=11, labelpad=10)
            
            ax.set_xticks(range(len(df_plot.columns)))
            ax.set_xticklabels(df_plot.columns, rotation=45, ha='left', fontsize=10)
            ax.set_yticks(range(len(df_plot.index)))
            ax.set_yticklabels(df_plot.index, fontsize=10)
            ax.xaxis.set_ticks_position('bottom')
            
            ax.set_title(f"Chromosomal Distribution of Distance=1 Transcripts (%)\n({species_name})", fontsize=13, fontweight='bold', pad=20)
            ax.set_xlabel('Chromosomes / Scaffolds', fontsize=11, labelpad=10)
            ax.set_ylabel('Annotation pipeline', fontsize=11, labelpad=10)
            
            plt.tight_layout()
            
            output_svg = f"{output_dir}/{clean_name}_chrom_distance1_percentage_distribution.svg"
            plt.savefig(output_svg, format='svg', dpi=300)
            plt.close()
            print(f"SVG Percentage Plot successfully saved to: {output_svg}")
            
        except Exception as graph_err:
            print(f"Error plotting SVG: {graph_err}")
    else:
        print(f" No records found for {species_name}.")

if __name__ == "__main__":
    
    # Target directory configuration matching the output structure
    output_dir = "results_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    # Path configuration for input datasets located in results_dataset_test
    input_base_dir = "results_dataset_test"

    species_list = [
        {"folder": "Arabidopsis_thaliana_dataset_test", "name": "Arabidopsis thaliana", "ref": "ARAPORT11"},
        {"folder": "Nicotiana_benthamiana_dataset_test", "name": "Nicotiana benthamiana", "ref": "GuoTTv1"},
        {"folder": "Oryza_sativa_dataset_test", "name": "Oryza sativa", "ref": "AGIS1v1"},
        {"folder": "Solanum_lycopersicum_datset_test", "name": "Solanum lycopersicum", "ref": "ITAG5.0"}
    ]

    for species in species_list:
        # Construct path relative to TFG execution root folder
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
            'test_GEM': 'GEMOMA', 'test_UTR_GEM': 'GEMOMA + UTR',
            'test_HEL': 'HELIXER', 'test_UTR_HEL': 'HELIXER + UTR',
            'test_MAK': 'MAKER', 'test_UTR_MAK': 'MAKER + UTR',
            'test_REF': f'{species["ref"]}', 'test_UTR_REF': f'{species["ref"]} + UTR'
        }
        
        extract_chrom_stats_for_distance_1(path, mapping, species["name"], output_dir)