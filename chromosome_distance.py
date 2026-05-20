import os
import re
import pandas as pd
import matplotlib.pyplot as plt

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
            cbar_label = '% of Proteins with Distance = 1 (Relative to Chromosome Total)'
            plot_title = f"Chromosomal Distribution of Distance=1 Proteins (%)\n({species_name})"
            output_csv = f"{output_dir}/{clean_name}_chrom_protein_distance1_summary.csv"
            output_svg = f"{output_dir}/{clean_name}_chrom_protein_distance1_percentage_distribution.svg"
        else:
            cbar_label = '% of Transcripts with Distance = 1 (Relative to Chromosome Total)'
            plot_title = f"Chromosomal Distribution of Distance=1 Transcripts (%)\n({species_name})"
            output_csv = f"{output_dir}/{clean_name}_chrom_transcript_distance1_summary.csv"
            output_svg = f"{output_dir}/{clean_name}_chrom_transcript_distance1_percentage_distribution.svg"

        df.to_csv(output_csv, index=False)
        print(f"Table saved to: {output_csv}")
        
        try:
            df_plot = df.pivot(index='Pipeline', columns='Chromosome', values='Percentage (Dist 1 %)').fillna(0)
            
            ordered_pipelines = [name for name in annotation_mapping.values() if name in df_plot.index]
            seen = set()
            unique_pipelines = [x for x in ordered_pipelines if not (x in seen or seen.add(x))]
            
            # Reindex utilizing the predefined strict pipeline order
            df_plot = df_plot.reindex(unique_pipelines[::-1])
            
            # Map new dynamic labels containing the absolute total counts: Pipeline (N=X)
            new_labels = [f"{idx} (N={pipeline_totals.get(idx, 0)})" for idx in df_plot.index]
            
            fig, ax = plt.subplots(figsize=(13, 9))
            
            # PALETTE CHANGE: Using 'viridis' (high-contrast academic standard)
            cax = ax.matshow(df_plot, cmap='viridis', aspect='auto')
            
            cbar = fig.colorbar(cax, pad=0.02)
            cbar.set_label(cbar_label, fontsize=11, labelpad=10)
            
            ax.set_xticks(range(len(df_plot.columns)))
            ax.set_xticklabels(df_plot.columns, rotation=45, ha='left', fontsize=10)
            
            # Apply the new dynamically computed labels to the Y axis
            ax.set_yticks(range(len(df_plot.index)))
            ax.set_yticklabels(new_labels, fontsize=10)
            ax.xaxis.set_ticks_position('bottom')
            
            ax.set_title(plot_title, fontsize=13, fontweight='bold', pad=20)
            ax.set_xlabel('Chromosomes / Scaffolds', fontsize=11, labelpad=10)
            ax.set_ylabel('Annotation pipeline', fontsize=11, labelpad=10)
            
            plt.tight_layout()
            plt.savefig(output_svg, format='svg', dpi=300)
            plt.close()
            print(f"SVG