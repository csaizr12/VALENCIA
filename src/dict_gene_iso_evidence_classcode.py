import os

# parse a reference map file to update a gene-isoform dictionary
# with evidence types a class codes 
def add_refmap_info(gene_isoform_dict, refmap_path):
    # identify evidence type
    fname = os.path.basename(refmap_path)
    if 'transcripts_evidence' in fname:
        evidence_type = 'transcripts'
    elif 'proteins_evidence' in fname:
        evidence_type = 'proteins'
    else:
        evidence_type = 'unknown'
        
    # process the file
    with open(refmap_path, "r") as f:
        for line in f:
            if line.startswith("ref_gene"):
                continue
            fields = line.strip().split('\t')
            class_code = fields[2]
            id_list = fields[3]
            # parse IDs and update dictionary
            parts = id_list.split("|")
            gen_id = parts[0].strip()
            iso_list = parts[1].strip()
            isoforms = [i.strip() for i in iso_list.split(",")]
            # check if the Gene ID exists in the main dictionary
            if gen_id not in gene_isoform_dict:
                    continue
            # if the specific Isoform ID exists within that Gene´s entry:
            for isoform_id in isoforms:
                    if isoform_id in gene_isoform_dict[gen_id]:
                     gene_isoform_dict[gen_id][isoform_id].update({evidence_type: class_code})
            
    return gene_isoform_dict


        





        

