import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def extract_chrom_stats_for_distance_1(base_path, annotation_mapping, species_name):
    """
    Scans GFF3 files, filters records with distance="1", counts the
    total number of transcripts and unique genes per chromosome, and
    generates both a summary CSV and an elegant SVG heatmap.
    """
    if not os.path.exists(base_path):
        return

    print(f"\n Analyzing Chromosomes for {species_name} (Distance = 1)...")
    
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
            chrom_transcripts = {}
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
                        
                        attributes = parts[8]
                        
                        if distance_pattern.search(attributes):
                            has_distance_data = True
                            chrom = parts[0].strip()
                            
                            chrom_transcripts[chrom] = chrom_transcripts.get(chrom, 0) + 1
                            
                            gene_match = gene_id_pattern.search(attributes)
                            if gene_match:
                                gene_id = gene_match.group(1)
                                if chrom not in chrom_genes:
                                    chrom_genes[chrom] = set()
                                chrom_genes[chrom].add(gene_id)
                
                if has_distance_data:
                    for chrom in sorted(chrom_transcripts.keys()):
                        total_tx = chrom_transcripts[chrom]
                        total_genes = len(chrom_genes.get(chrom, set()))
                        
                        species_results.append({
                            'Pipeline': official_name,
                            'Chromosome': chrom,
                            'Total Transcripts (Dist 1)': total_tx,
                            'Total Unique Genes': total_genes
                        })
                    break
                    
            except Exception as e:
                print(f"Error processing {file_name}: {e}")

    if species_results:
        df = pd.DataFrame(species_results)
        
        # 1. Guardar la tabla en CSV
        clean_name = species_name.lower().replace(' ', '_')
        output_csv = f"{clean_name}_chrom_distance1_summary.csv"
        df.to_csv(output_csv, index=False)
        print(f"Table saved to: {output_csv}")
        
        # 2. Generar el gráfico SVG (Heatmap de distribución)
        try:
            # Pivotamos los datos para tener los Pipelines en las filas y los Cromosomas en las columnas
            df_plot = df.pivot(index='Pipeline', columns='Chromosome', values='Total Transcripts (Dist 1)').fillna(0)
            
            # Asegurar el orden correcto inverso en el eje Y para los pipelines mapeados
            ordered_pipelines = [name for name in annotation_mapping.values() if name in df_plot.index]
            seen = set()
            unique_pipelines = [x for x in ordered_pipelines if not (x in seen or seen.add(x))]
            df_plot = df_plot.reindex(unique_pipelines[::-1])
            
            # Dibujar la figura
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Usamos matshow para crear el mapa de calor con una paleta elegante ('YlGnBu' o 'Viridis')
            cax = ax.matshow(df_plot, cmap='YlGnBu', aspect='auto')
            
            # Añadir barra de color (leyenda de intensidad)
            cbar = fig.colorbar(cax, pad=0.02)
            cbar.set_label('Transcript Count (Distance = 1)', fontsize=11, labelpad=10)
            
            # Configurar las etiquetas de los ejes
            ax.set_xticks(range(len(df_plot.columns)))
            ax.set_xticklabels(df_plot.columns, rotation=45, ha='left', fontsize=10)
            ax.set_yticks(range(len(df_plot.index)))
            ax.set_yticklabels(df_plot.index, fontsize=10)
            
            # Mover los ticks de las X a la parte inferior (Matplotlib por defecto los pone arriba en matshow)
            ax.xaxis.set_ticks_position('bottom')
            
            # Títulos y formato académico
            ax.set_title(f"Chromosomal distribution of distance=1 transcripts\n({species_name})", fontsize=13, fontweight='bold', pad=20)
            ax.set_xlabel('Chromosomes / Scaffolds', fontsize=11, labelpad=10)
            ax.set_ylabel('Annotation pipeline', fontsize=11, labelpad=10)
            
            plt.tight_layout()
            
            # Guardar el gráfico como vector SVG listo para la memoria
            output_svg = f"{clean_name}_chrom_distance1_distribution.svg"
            plt.savefig(output_svg, format='svg', dpi=300)
            plt.close()
            print(f"SVG Plot successfully saved to: {output_svg}")
            
        except Exception as graph_err:
            print(f"Error plotting SVG: {graph_err}")
            
    else:
        print(f" No records found for {species_name}.")

if __name__ == "__main__":
    
    def get_valid_path(folder):
        if os.path.exists(folder): return folder
        if os.path.exists(f"../{folder}"): return f"../{folder}"
        return None

    species_list = [
        {"folder": "Arabidopsis_thaliana_dataset_test", "name": "Arabidopsis thaliana", "ref": "ARAPORT11"},
        {"folder": "Nicotiana_benthamiana_dataset_test", "name": "Nicotiana benthamiana", "ref": "GuoTTv1"},
        {"folder": "Oryza_sativa_dataset_test", "name": "Oryza sativa", "ref": "AGIS1v1"},
        {"folder": "Solanum_lycopersicum_datset_test", "name": "Solanum lycopersicum", "ref": "ITAG5.0"}
    ]

    for species in species_list:
        path = get_valid_path(species["folder"])
        if not path:
            continue
            
        mapping = {
            'test_BRAT': 'BRAKER1', 'test_UTR_BRAT': 'BRAKER1 + UTR',
            'test_BRAP': 'BRAKER2', 'test_UTR_BRAP': 'BRAKER2 + UTR',
            'test_BRATP': 'BRAKER3', 'test_UTR_BRATP': 'BRAKER3 + UTR',
            'test_EGX': 'Eukaryotic Genome Annotation Pipeline', 'test_UTR_EGX': 'Eukaryotic Genome Annotation Pipeline + UTR',
            'test_EVI': 'Evidence-based', 'test_UTR_EVI': 'Evidence-based + UTR',
            'test_EVO': 'ANNEVO', 'test_UTR_EVO': 'ANNEVO + UTR',
            'test_GEM': 'Genome-scale Metabolic Model', 'test_GS': 'Genome-scale Metabolic Model',
            'test_UTR_GEM': 'Genome-scale Metabolic Model + UTR', 'test_UTR_GS': 'Genome-scale Metabolic Model + UTR',
            'test_HEL': 'HELIXER', 'test_UTR_HEL': 'HELIXER + UTR',
            'test_MAK': 'MAKER', 'test_UTR_MAK': 'MAKER + UTR',
            'test_REF': f'{species["ref"]}', 'test_UTR_REF': f'{species["ref"]} + UTR'
        }
        
        extract_chrom_stats_for_distance_1(path, mapping, species["name"])