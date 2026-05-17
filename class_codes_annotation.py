import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def plot_species_class_codes(folder_name, species_name, output_filename):
    """
    Scans the subfolders of a specific dataset, parses the GFF3 files,
    and generates a single horizontal stacked bar plot for that species.
    """
    # COMPROBACIÓN AUTOMÁTICA DE RUTAS: Busca de forma local o un nivel arriba
    if os.path.exists(folder_name):
        base_path = folder_name
    elif os.path.exists(f"../{folder_name}"):
        base_path = f"../{folder_name}"
    else:
        print(f"Error: The base path '{folder_name}' could not be found locally or in '../'. Skipping {species_name}.")
        return

    # Diccionario de referencias específicas por especie
    species_references = {
        "Arabidopsis thaliana": "ARAPORT11",
        "Nicotiana benthamiana": "GuoTTv1",
        "Oryza sativa": "AGIS1v1",
        "Solanum lycopersicum": "ITAG5.0"
    }
    
    # Obtenemos la referencia correspondiente o dejamos "REFERENCE" por defecto
    ref_name = species_references.get(species_name, "REFERENCE")

    # Mapeo estricto de las filas del gráfico con placeholders dinámicos para la referencia
    annotation_mapping = {
        'test_BRAT': 'BRAKER1', 
        'test_UTR_BRAT': 'BRAKER1 + UTR',
        'test_BRAP': 'BRAKER2',
        'test_UTR_BRAP': 'BRAKER2 + UTR',
        'test_BRATP': 'BRAKER3',
        'test_UTR_BRATP': 'BRAKER3 + UTR',
        'test_EGX': 'Eukaryotic Genome Annotation Pipeline',
        'test_UTR_EGX': 'Eukaryotic Genome Annotation Pipeline + UTR',
        'test_EVI': 'Evidence-based',
        'test_UTR_EVI': 'Evidence-based + UTR',
        'test_EVO': 'ANNEVO',
        'test_UTR_EVO': 'ANNEVO + UTR',
        'test_GEM': 'Genome-scale Metabolic Model',
        'test_GS': 'Genome-scale Metabolic Model',
        'test_UTR_GEM': 'Genome-scale Metabolic Model + UTR',
        'test_UTR_GS': 'Genome-scale Metabolic Model + UTR',
        'test_HEL': 'HELIXER',
        'test_UTR_HEL': 'HELIXER + UTR',
        'test_MAK': 'MAKER',
        'test_UTR_MAK': 'MAKER + UTR',
        'test_REF': f'{ref_name}',          
        'test_UTR_REF': f'{ref_name} + UTR' 
    }

    legend_mapping = {
        'complete': 'Complete (=)',
        'SubsequencesTarget': 'Subsequences target (c)',
        'PotentialIsoform': 'Potential isoform (j)',
        'PartialExonOverlap': 'Partial exon overlap (o)',
        'TotalIntronsRetention': 'Total introns retention (m)',
        'RetainedIntronSingleExon': 'Retained intron single exon (e)',
        'Unknown': 'Unknown intergenic (u)',
        'ExonicOverlapOppStran': 'Exonic overlap Opp. (x)',
        'SubsequencesReferences': 'Subsequences references (k)',
        'NA': 'Not annotated (NA)'
    }

    class_code_pattern = re.compile(r'class_code[ =]"?([^";\s]+)"?')
    all_data = []

    print(f"\nProcessing dataset for: {species_name} (Reference: {ref_name})...")

    # Process each technical folder in the defined order
    for technical_folder, official_name in annotation_mapping.items():
        folder_path = os.path.join(base_path, technical_folder)
        
        if not os.path.exists(folder_path):
            continue # Silently skip missing folders (handled automatically per species)
            
        gff_files = [f for f in os.listdir(folder_path) if f.endswith('.gff3')]
        counts = {}
        file_detected = False
        
        for file_name in gff_files:
            full_path = os.path.join(folder_path, file_name)
            contains_class_code = False
            
            try:
                with open(full_path, 'r', errors='ignore') as f:
                    for _ in range(100):
                        line = f.readline()
                        if not line:
                            break
                        if 'class_code' in line:
                            contains_class_code = True
                            break
                
                if contains_class_code:
                    file_detected = True
                    with open(full_path, 'r', errors='ignore') as f:
                        for line in f:
                            if line.startswith('#') or not line.strip():
                                continue
                            parts = line.split('\t')
                            if len(parts) < 9:
                                continue
                            
                            match = class_code_pattern.search(parts[8])
                            if match:
                                code = match.group(1).strip()
                                counts[code] = counts.get(code, 0) + 1
                    break # Valid file found, stop scanning this subfolder
                    
            except Exception as e:
                print(f"Error reading {file_name}: {e}")
        
        if file_detected and counts:
            for code, count in counts.items():
                clean_code_name = legend_mapping.get(code, code)
                all_data.append({
                    'Annotation': official_name, 
                    'Class Code': clean_code_name, 
                    'Count': count
                })

    df = pd.DataFrame(all_data)
    if df.empty:
        print(f"No valid data collected for {species_name}.")
        return

    # Data structuring and sorting
    df_pivot = df.pivot(index='Annotation', columns='Class Code', values='Count').fillna(0)
    ordered_names = [name for name in annotation_mapping.values() if name in df_pivot.index]
    
    # Remove duplicates from ordered_names while preserving order
    seen = set()
    unique_ordered_names = [x for x in ordered_names if not (x in seen or seen.add(x))]
    df_pivot = df_pivot.reindex(unique_ordered_names[::-1])
    
    # High-contrast color palette
    distinct_colors = [
        '#0f4c81', '#ffb347', '#1d8348', '#f1948a', '#6c3483', 
        '#bb8fce', '#117a65', '#73c6b6', '#9a7d0a', '#f7dc6f',
        '#5dade2', '#2e4053', '#aeb6bf', '#edbb99', '#f5b041'
    ]
    
    fig, ax = plt.subplots(figsize=(14, 12))
    df_pivot.plot(kind='barh', stacked=True, color=distinct_colors, ax=ax)
    
    # Dynamic title application
    ax.set_title(f"Distribution of Class Codes by Annotation Type ({species_name})", fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Total transcript count', fontsize=12, labelpad=10)
    ax.set_ylabel('Annotation pipeline', fontsize=12, labelpad=10)
    ax.legend(title='GFF3 class codes', bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_filename, format='svg', dpi=300)
    plt.close() # Close plot to free memory


if __name__ == "__main__":
    
    # 1. ARABIDOPSIS THALIANA (Usa ARAPORT11)
    plot_species_class_codes(
        folder_name="Arabidopsis_thaliana_dataset_test",
        species_name="Arabidopsis thaliana",
        output_filename="athaliana_class_codes_comparison.svg"
    )

    # 2. NICOTIANA BENTHAMIANA (Usa GuoTTv1)
    plot_species_class_codes(
        folder_name="Nicotiana_benthamiana_dataset_test",
        species_name="Nicotiana benthamiana",
        output_filename="nbenthamiana_class_codes_comparison.svg"
    )

    # 3. ORYZA SATIVA (Usa AGIS1v1)
    plot_species_class_codes(
        folder_name="Oryza_sativa_dataset_test",
        species_name="Oryza sativa",
        output_filename="osativa_class_codes_comparison.svg"
    )

    # 4. SOLANUM LYCOPERSICUM (Usa ITAG5.0)
    plot_species_class_codes(
        folder_name="Solanum_lycopersicum_datset_test",
        species_name="Solanum lycopersicum",
        output_filename="slycopersicum_class_codes_comparison.svg"
    )