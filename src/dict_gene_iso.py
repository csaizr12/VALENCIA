import re

def get_gene_isoform_dict_from_target_annotation(target_annotation):
    gene_isoform_dict = {}

    # parse the target annotation file and fill the gene
    with open(target_annotation, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            features = fields[2]
            attributes = fields[8]
            # consider only gene features to fill the gene_isoform_dict
            if features == "gene":
                gene_id_match = re.search(r'ID=([^;]+)', attributes)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    gene_isoform_dict.setdefault(gene_id, {})
    
    # parse the target annotation file again to fill the isoforms for each gene
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
            if not id_match:
                continue
            isoform_id = id_match.group(1)
            parent_id = parent_match.group(1) if parent_match else None
            # if the parent_id is in the gene_isoform_dict, add the isoform_id to the corresponding gene
            if parent_id in gene_id:
                gene_isoform_dict[parent_id].setdefault(isoform_id, {})
                continue
            # if the isoform_id is not in the gene_isoform_dict, check if it starts with any of the gene_ids and add it to the corresponding gene
            for gen in gene_isoform_dict:
                if isoform_id != gen and isoform_id.startswith(gen):
                    gene_isoform_dict[gen].setdefault(isoform_id, {})
                    break

    return gene_isoform_dict
            