import re

def get_gene_isoform_dict_from_target_annotation(target_annotation):
    gene_isoform_dict = {}

    # first parse to get genes
    with open(target_annotation, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attributes = fields[8]
            id_match = re.search(r'gene_id=([^;]+)', attributes)
            parent_match = re.search(r'Parent=([^;]+)', attributes)

            if id_match and not parent_match:
                gene_id = id_match.group(1)
                gene_isoform_dict.setdefault(gene_id, {})

    # second parse to get isoforms and assign them to genes
    with open(target_annotation, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attributes = fields[8]
            # obtein the gene_id and isoform_id from the attributes
            id_match = re.search(r'gene_id=([^;]+)', attributes) 
            parent_match = re.search(r'Parent=([^;]+)', attributes)
            if not id_match or not parent_match:
                continue
            isoform_id = id_match.group(1)
            parent_id = parent_match.group(1) if parent_match else None

            if parent_id in gene_isoform_dict:
                gene_isoform_dict[parent_id].setdefault(isoform_id, {})

    return gene_isoform_dict
