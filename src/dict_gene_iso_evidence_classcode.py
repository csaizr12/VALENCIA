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
            
            gene_id, isoform_list = fields[3].split('|', 1)

            target_gene = gene_isoform_dict.get(gene_id.strip())
            
            if target_gene:
                for isoform in isoform_list.split(','):
                    iso_id = isoform.strip()
                    if iso_id in target_gene:
                        target_gene[iso_id].update({evidence_type: class_code})

    return gene_isoform_dict


        





        

