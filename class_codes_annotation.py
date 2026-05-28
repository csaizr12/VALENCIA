import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def plot_valencia_class_codes(folder_name, species_name, output_filename, evaluation_type='transcript'):
    """
    Scans the subfolders of a specific dataset, parses the GFF3 files,
    and generates a horizontal stacked bar plot using a single unified logic.
    """
    # PATH AUTO-CHECK: Configured to find inputs inside 'results_dataset_test' from TFG root
    possible_paths = [
        folder_name,
        os.path.join("dataset_test", folder_name),
        os.path.join("VALENCIA", folder_name),
        os.path.join("VALENCIA", "dataset_test", folder_name)
    ]
    
    base_path = None
    for path in possible_paths:
        if os.path.exists(path):
            base_path = path
            break

    if not base_path:
        print(f"Error: Base path for '{folder_name}' could not be resolved. Skipping {species_name} [{evaluation_type}].")
        return

    # Species-specific reference mapping dictionary
    species_references = {
        "Arabidopsis thaliana": "ARAPORT11",
        "Nicotiana benthamiana": "GuoTTv1",
        "Oryza sativa": "AGIS1v1",
        "Solanum lycopersicum": "ITAG5.0"
    }
    
    ref_name = species_references.get(species_name, "REFERENCE")

    # Strict row mapping with dynamic reference name placeholders
    annotation_mapping = {
        'test_BRAT': 'BRAKER1', 
        'test_UTR_BRAT': 'BRAKER1 + UTR',
        'test_BRAP': 'BRAKER2',
        'test_UTR_BRAP': 'BRAKER2 + UTR',
        'test_BRATP': 'BRAKER3',
        'test_UTR_BRATP': 'BRAKER3 + UTR',
        'test_EGX': 'Eukaryotic Genome Annotation Pipeline',
        'test_UTR_EGX': 'Eukaryotic Genome Annotation Pipeline + UTR',
        'test_EVI': 'EviAnn',
        'test_UTR_EVI': 'EviAnn + UTR',
        'test_EVO': 'ANNEVO',
        'test_UTR_EVO': 'ANNEVO + UTR',
        'test_GEM': 'GEMOMA',
        'test_UTR_GEM': 'GEMOMA + UTR',
        'test_HEL': 'HELIXER',
        'test_UTR_HEL': 'HELIXER + UTR',
        'test_MAK': 'MAKER',
        'test_UTR_MAK': 'MAKER + UTR',
        'test_REF': f'{ref_name}',          
        'test_UTR_REF': f'{ref_name} + UTR' 
    }

    # DYNAMIC METADATA AND LEGEND CONFIGURATION BASED ON EVALUATION TYPE
    if evaluation_type == 'protein':
        legend_mapping = {
            'complete': 'Complete Structural Match (=)',
            'SubsequencesTarget': 'Subsequences target (c)',
            'PotentialIsoform': 'Potential isoform (j)',
            'PartialExonOverlap': 'Partial exon overlap (o)',
            'TotalIntronsRetention': 'Total introns retention (m)',
            'RetainedIntronSingleExon': 'Retained intron single exon (e)',
            'PartialIntronRetention': 'Partial intron retention (n)',
            'SubsequencesReferences': 'Subsequences references (k)',
        }
        title_text = f"Distribution of protein class codes ({species_name})"
        x_label_text = "Total protein count"
        legend_title_text = "VALENCIA protein class codes"
    else:
        legend_mapping = {
            'complete': 'Complete Structural Match (=)',
            'SubsequencesTarget': 'Subsequences target (c)',
            'PotentialIsoform': 'Potential isoform (j)',
            'PartialExonOverlap': 'Partial exon overlap (o)',
            'TotalIntronsRetention': 'Total introns retention (m)',
            'RetainedIntronSingleExon': 'Retained intron single exon (e)',
            'PartialIntronRetention': 'Partial intron retention (n)',
            'SubsequencesReferences': 'Subsequences references (k)',
        }
        title_text = f"Distribution of transcript class codes ({species_name})"
        x_label_text = "Total transcript count"
        legend_title_text = "GFF3 transcript class codes"

    class_code_pattern = re.compile(r'class_code[ =]"?([^";\s]+)"?')
    all_data = []

    print(f"Running VALENCIA [{evaluation_type.upper()}] pipeline for: {species_name}...")

    # Core GFF3 file parsing algorithm
    for technical_folder, official_name in annotation_mapping.items():
        folder_path = os.path.join(base_path, technical_folder)
        
        if not os.path.exists(folder_path):
            continue 
            
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
                    break 
                    
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
        print(f"No valid data collected for {species_name} [{evaluation_type}].")
        return

    # Data pivoting and reverse indexing (identical sorting for both modes)
    df_pivot = df.pivot(index='Annotation', columns='Class Code', values='Count').fillna(0)
    ordered_names = [name for name in annotation_mapping.values() if name in df_pivot.index]
    
    seen = set()
    unique_ordered_names = [x for x in ordered_names if not (x in seen or seen.add(x))]
    df_pivot = df_pivot.reindex(unique_ordered_names[::-1])
    
    # High-contrast color palette configuration
    distinct_colors = [
        '#0f4c81', '#ffb347', '#1d8348', '#f1948a', '#6c3483', 
        '#bb8fce', '#117a65', '#73c6b6', '#9a7d0a', '#f7dc6f',
        '#5dade2', '#2e4053', '#aeb6bf', '#edbb99', '#f5b041'
    ]
    
    fig, ax = plt.subplots(figsize=(14, 12))
    df_pivot.plot(kind='barh', stacked=True, color=distinct_colors, ax=ax)
    
    # Plot formatting and styling execution
    ax.set_title(title_text, fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel(x_label_text, fontsize=12, labelpad=10)
    ax.set_ylabel('Annotation pipeline', fontsize=12, labelpad=10)
    ax.legend(title=legend_title_text, bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_filename, format='svg', dpi=300)
    plt.close()


if __name__ == "__main__":
    
    # Target directory for the generated SVG figures
    output_dir = "results_plots"
    os.makedirs(output_dir, exist_ok=True)

    # 1. TRANSCRIPT EVALUATION
    plot_valencia_class_codes(
        folder_name="Arabidopsis_thaliana_dataset_test",
        species_name="Arabidopsis thaliana",
        output_filename=f"{output_dir}/athaliana_transcript_class_codes_comparison.svg",
        evaluation_type='transcript'
    )

    plot_valencia_class_codes(
        folder_name="Nicotiana_benthamiana_dataset_test",
        species_name="Nicotiana benthamiana",
        output_filename=f"{output_dir}/nbenthamiana_transcript_class_codes_comparison.svg",
        evaluation_type='transcript'
    )

    plot_valencia_class_codes(
        folder_name="Oryza_sativa_dataset_test",
        species_name="Oryza sativa",
        output_filename=f"{output_dir}/osativa_class_codes_comparison.svg",
        evaluation_type='transcript'
    )

    plot_valencia_class_codes(
        folder_name="Solanum_lycopersicum_datset_test",
        species_name="Solanum lycopersicum",
        output_filename=f"{output_dir}/slycopersicum_transcript_class_codes_comparison.svg",
        evaluation_type='transcript'
    )



    # 2. PROTEIN EVALUATION
    plot_valencia_class_codes(
        folder_name="Arabidopsis_thaliana_dataset_test",
        species_name="Arabidopsis thaliana",
        output_filename=f"{output_dir}/athaliana_protein_class_codes_comparison.svg",
        evaluation_type='protein'
    )

    plot_valencia_class_codes(
        folder_name="Nicotiana_benthamiana_dataset_test",
        species_name="Nicotiana benthamiana",
        output_filename=f"{output_dir}/nbenthamiana_protein_class_codes_comparison.svg",
        evaluation_type='protein'
    )

    plot_valencia_class_codes(
        folder_name="Oryza_sativa_dataset_test",
        species_name="Oryza sativa",
        output_filename=f"{output_dir}/osativa_protein_class_codes_comparison.svg",
        evaluation_type='protein'
    )

    plot_valencia_class_codes(
        folder_name="Solanum_lycopersicum_datset_test",
        species_name="Solanum lycopersicum",
        output_filename=f"{output_dir}/slycopersicum_protein_class_codes_comparison.svg",
        evaluation_type='protein'
    )