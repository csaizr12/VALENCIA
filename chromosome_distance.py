import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def generate_plots_from_files(archivos_gff):
    output_dir = "results_plots"
    os.makedirs(output_dir, exist_ok=True)

    distance_pattern = re.compile(r'distance[ =]"?1"?')
    gene_id_pattern = re.compile(r'gene_id[ =]"?([^";\s]+)"?')

    # Process each evidence transcript and protein separately
    for mode in ['transcript', 'protein']:
        results = []
        pipeline_totals = {}

        for full_path in archivos_gff:
            if not os.path.exists(full_path):
                continue
                
            # Obtein pipeline and species names from the file path
            pipeline_name = os.path.basename(os.path.dirname(full_path))
            species_raw = os.path.basename(os.path.dirname(os.path.dirname(full_path)))
            
            # Clean species name for display and keys
            species_name = species_raw.replace('_dataset_test', '').replace('_datset_test', '').replace('_', ' ')
            if not species_name:
                species_name = "Unknown_Species"

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
            except Exception as e:
                print(f"  Processing error {full_path}: {e}")
                continue

            # Keep results only if we found distance=1 data for this file
            if has_distance_data:
                # Use species_name and pipeline_name as keys for totals
                key = (species_name, pipeline_name)
                pipeline_totals[key] = sum(chrom_total_elements.values())
                
                for chrom in sorted(chrom_total_elements.keys()):
                    total_el_chrom = chrom_total_elements[chrom]
                    dist1_el_chrom = chrom_distance1_elements.get(chrom, 0)
                    
                    percentage_dist1 = (dist1_el_chrom / total_el_chrom) * 100 if total_el_chrom > 0 else 0.0
                    total_genes = len(chrom_genes.get(chrom, set()))
                    
                    results.append({
                        'Species': species_name,
                        'Pipeline': pipeline_name,
                        'Chromosome': chrom,
                        'Percentage (Dist 1 %)': round(percentage_dist1, 2),
                        'Total Unique Genes': total_genes
                    })

        df = pd.DataFrame(results)

        # Generating plots for each species
        for species in df['Species'].unique():
            df_species = df[df['Species'] == species]
            clean_name = species.lower().replace(' ', '_')
            
            if mode == 'protein':
                cbar_label = '% of Proteins with distance = 1 (relative to chromosome total)'
                plot_title = f"Chromosomal distribution of edit distance = 1 | Proteins (%)\n({species})"
                output_svg = f"{output_dir}/{clean_name}_chrom_protein_distance1_percentage_distribution.svg"
            else:
                cbar_label = '% of Transcripts with distance = 1 (relative to chromosome total)'
                plot_title = f"Chromosomal distribution of edit distance = 1 | Transcripts (%)\n({species})"
                output_svg = f"{output_dir}/{clean_name}_chrom_transcript_distance1_percentage_distribution.svg"
            
            try:
                df_plot = df_species.pivot_table(index='Pipeline', columns='Chromosome', values='Percentage (Dist 1 %)', aggfunc='mean').fillna(0)
                
                # Reorder pipelines to have a consistent order across species
                ordered_pipelines = sorted(df_plot.index, reverse=True)
                df_plot = df_plot.reindex(ordered_pipelines)
                
                new_labels = [f"{idx} (N={pipeline_totals.get((species, idx), 0)})" for idx in df_plot.index]
                
                fig, ax = plt.subplots(figsize=(14, 10))
                discrete_heatmap_cmap = plt.get_cmap('Spectral_r', 20)            
                cax = ax.matshow(df_plot, cmap=discrete_heatmap_cmap, aspect='auto', vmin=0, vmax=100)

                ax.set_xticks(np.arange(-0.5, len(df_plot.columns), 1), minor=True)
                ax.set_yticks(np.arange(-0.5, len(df_plot.index), 1), minor=True)
                ax.grid(which='minor', color='#b0b0b0', linestyle='-', linewidth=0.4)

                cbar = fig.colorbar(cax, pad=0.02, ticks=np.arange(0, 101, 5))
                cbar.set_label(cbar_label, fontsize=11, labelpad=10)
                cbar.ax.tick_params(labelsize=9)
                
                ax.set_xticks(range(len(df_plot.columns)))
                ax.set_xticklabels(df_plot.columns, rotation=90, ha='right', fontsize=9)
                
                ax.set_yticks(range(len(df_plot.index)))
                ax.set_yticklabels(new_labels, fontsize=10)
                ax.xaxis.set_ticks_position('bottom')
                
                ax.set_title(plot_title, fontsize=13, fontweight='bold', pad=20)
                ax.set_xlabel('Chromosomes', fontsize=11, labelpad=10)
                ax.set_ylabel('Annotation pipeline', fontsize=11, labelpad=10)
                
                plt.subplots_adjust(bottom=0.25)
                plt.tight_layout()
                
                plt.savefig(output_svg, format='svg', dpi=300)
                plt.close()
                print(f" Success: {output_svg}")
                
            except Exception as graph_err:
                print(f" Generation error {species}: {graph_err}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python3 script.py <archivo1.gff3> <archivo2.gff3> ...")
        sys.exit(1)

    archivos_gff = sys.argv[1:]
    generate_plots_from_files(archivos_gff)