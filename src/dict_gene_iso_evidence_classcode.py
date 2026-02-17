import os

# define a function that updates a gene-isform dictionary
# using a reference map file
def add_refmap_info(gene_isoform_dict, refmap_path):
    # get only the filename part from the full reference map path
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
            # look up the gene_id in the gene_isoform_dict; if not found, get None.
            target_gene = gene_isoform_dict.get(gene_id.strip(), None)
            # if the gene is no present in the dictionary, skip and conitnue
            if target_gene is None:
                continue
            if target_gene:
                for isoform in isoform_list.split(','):
                    iso_id = isoform.strip()
                    if iso_id in target_gene:
                        target_gene[iso_id].update({evidence_type: class_code})

    return gene_isoform_dict


        

1



        

