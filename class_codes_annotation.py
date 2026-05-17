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
                print(f"   ✨ Correct file detected: {file_name}")
                with open(full_path, 'r', errors='ignore') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        parts = line.split('\t')
                        if len(parts) < 9:
                            continue
                        
                        match = class_code_pattern.search(parts[8])
                        if match:
                            code = match.group(1)
                            counts[code] = counts.get(code, 0) + 1
                return counts # Return counts and stop scanning this folder
                
        except Exception as e:
            print(f"   ⚠️ Error reading {file_name}: {e}")
            
    return None

def plot_valencia_class_codes(base_path, annotation_mapping, output_img="athaliana_class_codes_comparison.svg"):
    """
    Iterates through folders defined in the mapping dictionary and generates
    a horizontal stacked bar plot using official annotation names saved as SVG.
    """
    all_data = []

    # Dictionary to rename the class codes in the final legend
    legend_mapping = {
        '=': 'Complete Match (=)',
        'c': 'Contained (c)',
        'j': 'Alternative Splicing (j)',
        'o': 'Other Exon Overlap (o)',
        'p': 'Possible Polymerase Run-on (p)',
        'u': 'Unknown Intergenic', # Updated as requested
        'i': 'Intronic (i)',
        'x': 'Exonic Overlap Opp. Strand (x)',
        's': 'Intronic Opp. Strand (s)',
        'r': 'Repeat (r)',
        'k': 'Frg. Match Opp. Strand (k)'
    }

    # Iterate following the exact order defined in the dictionary
    for technical_folder, official_name in annotation_mapping.items():
        folder_path = os.path.join(base_path, technical_folder)
        
        if not os.path.exists(folder_path):
            print(f"⚠️ Skipping: Folder {technical_folder} does not exist")
            continue
            
        print(f"Searching for valid GFF3 in: {technical_folder}...")
        counts = detect_and_parse_gff(folder_path)
        
        if counts:
            for code, count in counts.items():
                # If the code is in our legend_mapping, use the long name, otherwise keep the code
                clean_code_name = legend_mapping.get(code, f"Code {code}")
                all_data.append({
                    'Annotation': official_name, 
                    'Class Code': clean_code_name, 
                    'Count': count
                })
        else:
            print(f"   ❌ No .gff3 file with 'class_code' found in {technical_folder}")

    df = pd.DataFrame(all_data)
    
    if df.empty:
        print("❌ No valid data collected for plotting. Check your folder names or file content.")
        return

    # Pivot data to structure horizontal stacked bars
    df_pivot = df.pivot(index='Annotation', columns='Class Code', values='Count').fillna(0)
    
    # Reindex the rows to preserve the exact order of the dictionary.
    ordered_names = [name for name in annotation_mapping.values() if name in df_pivot.index]
    df_pivot = df_pivot.reindex(ordered_names[::-1])
    
    # Generate the plot
    plt.figure(figsize=(14, 8))
    
    # CHANGED COLORS: Using 'viridis' colormap for a modern, high-contrast, publication-ready style
    # Other good options: 'plasma', 'inferno', 'Set2', 'Dark2'
    ax = df_pivot.plot(kind='barh', stacked=True, figsize=(14, 8), colormap='viridis')
    
    # Visual adjustments
    plt.title('Distribution of Class Codes by Annotation Type (Arabidopsis thaliana)', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Total Transcript Count', fontsize=12, labelpad=10)
    plt.ylabel('Annotation Pipeline', fontsize=12, labelpad=10)
    plt.legend(title='GFF3 Class Codes', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Export explicitly as high-quality vector SVG format
    plt.savefig(output_img, format='svg', dpi=300)
    plt.show()
    print(f"\n🎉 Plot successfully saved to: {output_img}!")


# =====================================================================
# SCRIPT CONFIGURATION & EXECUTION FOR ATHALIANA
# =====================================================================
if __name__ == "__main__":
    
    # 1. Target path pointing to the Arabidopsis dataset directory
    my_base_path = "./Arabidopsis_thaliana_dataset_test"
    
    # 2. Definitive mapping and strict ordering for your TFG
    my_annotations = {
        'test_BRAP': 'BRAKER2',
        'test_BRAT': 'BRAKER1', 
        'test_BRATP': 'BRAKER3', 
        'test_EGX': 'Eukaryotic Genome Annotation Pipeline',       
        'test_EVI': 'Evidence-based',
        'test_EVO': 'ANNEVO', 
        'test_GEM': 'Genome-scale Metabolic Model',
        'test_HEL': 'HELIXER',
        'test_MAK': 'MAKER',
        'test_REF': 'ARAPORT11'     
    }
    
    # Run the function
    plot_valencia_class_codes(
        base_path=my_base_path,
        annotation_mapping=my_annotations,
        output_img="athaliana_class_codes_comparison.svg"
    )