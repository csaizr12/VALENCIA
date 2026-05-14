import os
import re

# Function to filter pseudogenes from a gff file and write the clean gff and a log file with the filtered pseudogenes
def filter_pseudogenes(gff_file, outbase):
    clean_gff = os.path.join(outbase, "annotation_target_no_pseudogenes.gff")
    log_file = os.path.join(outbase, "filter_pseudogenes.log")

    os.makedirs(outbase, exist_ok=True)

    count_pseudogenes = 0
    # open gff file, clean gff file, and log file
    with open(gff_file, "r") as gff_input, open(clean_gff, "w") as gff_output, open(log_file, "w") as log_output:
        log_output.write("Filtering pseudogenes from {}\n".format(gff_file))
        
        for line in gff_input:
            if line.startswith("#"):
                gff_output.write(line)
                continue
            
            line_lower = line.lower()
            if "pseudogene" in line_lower or "pseudo=true" in line_lower:
                match = re.search(r'ID=([^;]+)', line)
                feature_id = match.group(1) if match else "NA"
                
                log_output.write(f"Excluding pseudogene: {feature_id}\n")
                count_pseudogenes += 1
            else:
                gff_output.write(line)
        
        log_output.write(f"Total pseudogenes filtered: {count_pseudogenes}\n")
    
    print(f"Filtering complete. {count_pseudogenes} pseudogenes removed. Clean GFF saved to {clean_gff}. Log saved to {log_file}.")
    
    return clean_gff