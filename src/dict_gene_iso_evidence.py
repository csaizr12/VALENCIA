def add_refmap_info(gene_isoform_dict, refmap_path):
    if 'trasncripts' in refmap_path:
        evidence_type = 'transcripts'
    elif 'proteins' in refmap_path:
        evidence_type = 'proteins'
    else:
        evidence_type = 'unknown'

    with open(refmap_path, "r") as f:
        for line in f:
            if line.startswith("ref_gene"):
                continue
            fields = line.strip().split('\t')
            class_code = fields[2]
            id_list = fields[3]

            gen_id, isoform_list = id_list.split("|")

            gen_id = gen_id.strip()
            isoforms = [i.strip() for i in isoform_list.split(",") ]
            if gen_id in gene_isoform_dict:
                    continue
            for isoform_id in isoforms:
                    gene_isoform_dict[gen_id][isoform_id]['class_code'] = class_code
                    gene_isoform_dict[gen_id][isoform_id]['evidence_type'] = evidence_type
         
    return gene_isoform_dict


        





        

