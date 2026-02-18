import os
import re


CLASS_CODE_TRANSLATION = {"=": "complete", "c": "SubsequencesTarget", "k":"SubsequencesReferences",
                          "m":"TotalIntronsRetention", "n":"PartialIntronRetention", "j":"PotentialIsoform",
                          "o": "PartialExonOverlap", "e":"RetainedIntronSingleExon"}

# parse a genomic annotation file to create a  
# dictionary structure of genes and their corresponding isoforms
def get_gene_isoform_dict_from_target_annotation(target_annotation):
    gene_isoform_dict = {}
    # file processing
    for line in target_annotation:
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
            class_code = CLASS_CODE_TRANSLATION[fields[2]]
            id_field = fields[1].strip()
            for gene in fields[3].split(","):
                gene_id, iso_id = gene.split('|')
            # look up the gene_id in the gene_isoform_dict; if not found, get None.
                target_gene = gene_isoform_dict.get(gene_id.strip(), None)
            # if the gene is no present in the dictionary, skip and conitnue
                if target_gene is None: 
                    continue
                if target_gene:
                    if iso_id in target_gene:
                        if id_field.startswith("MSTRG"):
                            target_gene[iso_id].setdefault(evidence_type, {})
                            target_gene[iso_id][evidence_type].update({"class_code": class_code, "match_sequence": id_field,})
                        if id_field.startswith("PAC"):
                            target_gene[iso_id].setdefault(evidence_type, {})
                            target_gene[iso_id][evidence_type].update({"class_code": class_code, "target_id": id_field})
    return gene_isoform_dict
