import os
def add_refmap_info(gene_isoform_dict, refmap_path):
    fname = os.path.basename(refmap_path).lower()
    if 'trasncripts' in fname:
        evidence_type = 'transcripts'
    elif 'proteins' in fname:
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

            parts = id_list.split("|")
            if len(parts) < 2:
                 continue
            gen_id = parts[0].strip()
            iso_list = parts[1].strip()

            isoforms = [i.strip() for i in iso_list.split(",")]
            if gen_id not in gene_isoform_dict:
                    continue
            for isoform_id in isoforms:
                    if isoform_id in gene_isoform_dict[gen_id]:
                     gene_isoform_dict[gen_id][isoform_id]['class_code'] = class_code
                     gene_isoform_dict[gen_id][isoform_id]['evidence_type'] = evidence_type
         
    return gene_isoform_dict


        





        

