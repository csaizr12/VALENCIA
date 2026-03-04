import sys
import re

# this function adds features to the original gff file based on the gene-isoform dictionary only for mRNA lines
# , and writes the output in a new gff file
def add_features_to_gff(outbase, gff_file, gene_isoform_dict):
 with open(gff_file, "r") as gff_input:
    with open(outbase / "Athaliana_447_Araport11.gene_exons_with_evidence_features.gff3", "w") as gff_output:
        # write a comment about the command used to run the script
        gff_output.write("##CMD:  {}\n".format(" ".join(sys.argv)))
        for line in gff_input:
            # if the line is a comment, write it as is
            if line.startswith("#"):
                    gff_output.write(line)
                    continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            type = fields[2]
            # if the line is not an mRNA, write it as is
            if type != "mRNA":
                gff_output.write(line)
            # if the line is an mRNA, we look for the gene and isoform in the gene_isoform_dict
            else:
                id_match = re.search(r'ID=([^;]+)', attributes)
                parent_match = re.search(r'Parent=([^;]+)', attributes)
                isoform_id = id_match.group(1)
                gene_id = parent_match.group(1)
                # look up the gene_id in the gene_isoform_dict; if not found, get None.
                target_gene = gene_isoform_dict.get(gene_id.strip(), None)
                # get the features for the isoform
                features = target_gene[isoform_id]
                evidence_info = []
                # for each evidence type, we add the transcript evidence and the protein evidence, class code and edit distance to the attributes or NA if not available
                for evidence_type, evidence_features in features.items():
                        if evidence_type == "transcripts_evidence":
                            evidence_match = evidence_features.get("match_sequence", 'NA')
                            class_code = evidence_features.get("class_code", "NA")
                            edit_distance = evidence_features.get("edit_distance", "NA")
                            evidence_info.append(f"{evidence_type}_match_sequence={evidence_match};{evidence_type}_class_code={class_code};{evidence_type}_edit_distance={edit_distance}")
                        elif evidence_type == "proteins_evidence":
                            evidence_match = evidence_features.get("match_sequence", 'NA')
                            class_code = evidence_features.get("class_code", "NA")
                            edit_distance = evidence_features.get("edit_distance", "NA")
                            evidence_info.append(f"{evidence_type}_match_sequence={evidence_match};{evidence_type}_class_code={class_code};{evidence_type}_edit_distance={edit_distance}")
                # if we have evidence info, we add it to the attributes; if not, we add evidence_info=NA
                if evidence_info:
                    new_attributes = attributes + ";" + ";".join(evidence_info)
                else:
                    new_attributes = attributes + ";evidence_info=NA"
                # replace the attributes field with the new attributes
                fields[8] = new_attributes
                gff_output.write("\t".join(fields) + "\n")
               
           

