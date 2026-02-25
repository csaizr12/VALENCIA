# necesito: archivo de arabidopsis araport 11 (gff), diccionario gen iso ()
import re

def add_features_to_gff(outbase, gff_file, gene_isoform_dict):
 with open(gff_file, "r") as gff_input:
    with open(outbase / "Athaliana_447_Araport11.gene_exons_edit_distance.gff3", "w") as gff_output:
        for line in gff_input:
            if line.startswith("#"):
                gff_output.write(line)
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            type = fields[2]
            if type != "mRNA":
                gff_output.write(line)
                continue
            else:
                id_match = re.search(r'ID=([^;]+)', attributes)
                parent_match = re.search(r'Parent=([^;]+)', attributes)
                if id_match and parent_match:
                    isoform_id = id_match.group(1)
                    gene_id = parent_match.group(1)
                    target_gene = gene_isoform_dict.get(gene_id.strip(), None)
                    if target_gene and isoform_id in target_gene:
                        features = target_gene[isoform_id]
                        evidence_info = []
                        for evidence_type, evidence_data in features.items():
                            class_code = evidence_data.get("class_code", "NA")
                            edit_distance = evidence_data.get("edit_distance", "NA")
                            evidence_info.append(f"{evidence_type}_class_code={class_code};{evidence_type}_edit_distance={edit_distance}")
                        new_attributes = attributes + ";" + ";".join(evidence_info)
                        fields[8] = new_attributes
                        gff_output.write("\t".join(fields) + "\n")
                    else:
                         gff_output.write(line)
        
    # save the new gff file in outbase with a name similar to the original but with added features
# abrimos el archivo arabidopsis Araport 11
    # abrimos el arhivo arabiopsis Araport 11 para editar 
        # lo dejamos como esta pero si estamos en la columna de RNAm 
            # añdimos la informacion del dict gen iso
# quiero que me guardes el resultado en outbase con un nombre similar