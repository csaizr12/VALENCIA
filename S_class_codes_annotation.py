import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def detect_and_parse_gff(folder_path):
    """
    Scans a folder, looks for .gff3 files, and parses the first one
    that contains the 'class_code' attribute.
    """
    class_code_pattern = re.compile(r'class_code[ =]"?([^";\s]+)"?')
    
    if not os.path.exists(folder_path):
        return None
        
    gff_files = [f for f in os.listdir(folder_path) if f.endswith('.gff3')]
    
    for file_name in gff_files:
        full_path = os.path.join(folder_path, file_name)
        counts = {}
        contains_class_code = False
        
        try:
            with open(full_path, 'r', errors='ignore') as f:
                # Fast header inspection (checks the first 100 lines)
                for _ in range(100):
                    line = f.readline()
                    if not line:
                        break
                    if 'class_code' in line:
                        contains_class_code = True
                        break
            
            # If the file contains class_codes, process it entirely
            if contains_class_code:
                with open(full_path, 'r', errors='ignore') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        parts = line.split('\t')
                        if len(parts) < 9:
                            continue
                        
                        match = class_code_pattern.search(parts[8])
                        if match:
                            code = match.group(1).strip() # Strip spaces to avoid mapping mismatches
                            counts[code] = counts.get(code, 0) + 1
                return counts # Return counts and stop scanning this folder
                
        except Exception as e:
            print(f"   Error reading {file_name}: {e}")
            
    return None

def plot_valencia_class_codes(base_path, annotation_mapping, output_filename="slycopersicum_class_codes_comparison.svg"):
    """
    Iterates through folders defined in the mapping dictionary and generates
    a horizontal stacked bar plot using official annotation names saved as SVG
    in the current results folder.
    """
    all_data = []

    # FIXED: The keys now match EXACTLY the values parsed from your GFF3 files
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

    # Iterate following the exact order defined in the dictionary
    for technical_folder, official_name in annotation_mapping.items():
        folder_path = os.path.join(base_path, technical_folder)
        
        if not os.path.exists(folder_path):
            print(f" Skipping: Folder '{folder_path}' does not exist")
            continue
            
        print(f"Searching for valid GFF3 in: {technical_folder}...")
        counts = detect_and_parse_gff(folder_path)
        
        if counts:
            for code, count in counts.items():
                # Map the code to its clean long description
                clean_code_name = legend_mapping.get(code, code)
                all_data.append({
                    'Annotation': official_name, 
                    'Class Code': clean_code_name, 
                    'Count': count
                })
        else:
            print(f"No .gff3 file with 'class_code' found in {technical_folder}")

    df = pd.DataFrame(all_data)
    
    if df.empty:
        print("No valid data collected for plotting. Check your folder paths or file content.")
        return

    # Pivot data to structure horizontal stacked bars
    df_pivot = df.pivot(index='Annotation', columns='Class Code', values='Count').fillna(0)
    
    # Reindex the rows to preserve the exact order of the dictionary
    ordered_names = [name for name in annotation_mapping.values() if name in df_pivot.index]
    df_pivot = df_pivot.reindex(ordered_names[::-1])
    
    # Custom high-contrast color palette for publication (Alternating Contrast Palette)
    distinct_colors = [
        '#0f4c81', '#ffb347', '#1d8348', '#f1948a', '#6c3483', 
        '#bb8fce', '#117a65', '#73c6b6', '#9a7d0a', '#f7dc6f',
        '#5dade2', '#2e4053', '#aeb6bf', '#edbb99', '#f5b041'
    ]
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    df_pivot.plot(kind='barh', stacked=True, color=distinct_colors, ax=ax)
    
    # Updated title for Solanum lycopersicum
    ax.set_title('Distribution of class codes by annotation type (Solanum lycopersicum)', fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Total transcript count', fontsize=12, labelpad=10)
    ax.set_ylabel('Annotation pipeline', fontsize=12, labelpad=10)
    ax.legend(title='GFF3 class codes', bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_filename, format='svg', dpi=300)
    plt.show()

if __name__ == "__main__":
    
    # Path updated to point to your tomato folder (incorporating the technical typo 'datset')
    my_base_path = "../Solanum_lycopersicum_datset_test"
    
    my_annotations = {
        # Standard Pipelines
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
        'test_GS': 'Genome-scale Metabolic Model',         # Match for your test_GS directory
        'test_UTR_GS': 'Genome-scale Metabolic Model + UTR', # Kept for consistency if needed later
        'test_HEL': 'HELIXER',
        'test_UTR_HEL': 'HELIXER + UTR',
        'test_MAK': 'MAKER',
        'test_UTR_MAK': 'MAKER + UTR',
        'test_REF': 'ARAPORT11',
        'test_UTR_REF': 'ARAPORT11 + UTR'
    }
    
    plot_valencia_class_codes(
        base_path=my_base_path,
        annotation_mapping=my_annotations,
        output_filename="slycopersicum_class_codes_comparison.svg"
    )