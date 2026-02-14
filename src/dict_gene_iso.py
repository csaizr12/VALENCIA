import re

# parse a genomic annotation file to create a  
# dictionary structure of genes and their corresponding isoforms
def get_gene_isoform_dict_from_target_annotation(target_annotation):
    gene_isoform_dict = {}
    # file processing
    with open(target_annotation, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            # obtein the id and parent match from the attributes
            id_match = re.search(r'ID=([^;]+)', attributes)
            parent_match = re.search(r'Parent=([^;]+)', attributes)
            # identify biological role
            if id_match:
                gene_id = id_match.group(1)
                # if there is no parent, the ID belongs to a gene
                if not parent_match:
                    if gene_id not in gene_isoform_dict:
                        gene_isoform_dict[gene_id] = {}
                # if there is a parent, the ID belongs to a isoform and the Parent belongs to a gen        
                elif parent_match:
                     isoform_id = id_match.group(1)
                     gene_id = parent_match.group(1)
                     # only add the isoform if the parent gen existe
                     if gene_id in gene_isoform_dict:
                         if isoform_id not in gene_isoform_dict[gene_id]:
                           gene_isoform_dict[gene_id][isoform_id] = {}
    return gene_isoform_dict
            