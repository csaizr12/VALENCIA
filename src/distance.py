from Bio import SeqIO
from Levenshtein import distance

def edit_distance(parsed_evidence_gene_isoform_dict, transcript_target_fasta, transcript_evidence_fasta,
                  protein_target_fasta, protein_evidence_fasta):
    # index fasta files
    records_transcript_target = SeqIO.index(transcript_target_fasta, "fasta")
    records_transcript_evidence = SeqIO.index(transcript_evidence_fasta, "fasta")
    records_protein_target = SeqIO.index(protein_target_fasta, "fasta")
    records_protein_evidence = SeqIO.index(protein_evidence_fasta, "fasta")

    # run gene_isoform_dict  
    for target_gene_id, evidence_found in parsed_evidence_gene_isoform_dict.items():
        for target_isoform_id, features in evidence_found.items():
            if not evidence_found:
                continue
            # obtein evidence type and mach_id
            print(target_gene_id, features)
            evidence_type = list(features.keys())[0]
            matching_evidence_id = features[evidence_type]['match_sequence']
            # if target_isoform_id and match_in in records, charge sequences
            print(evidence_type)
            if evidence_type == "transcripts":
                seq_target = str(records_transcript_target[target_isoform_id].seq)
                seq_evidence = str(records_transcript_evidence[matching_evidence_id].seq)
            elif evidence_type == "proteins":
                seq_target = str(records_protein_target[target_isoform_id].seq)
                seq_evidence = str(records_protein_evidence[matching_evidence_id].seq)
            # obtain edit distance with target and evidence
            edit_distance = distance(seq_target, seq_evidence)
            # save it in a gene_isoform_dict
            parsed_evidence_gene_isoform_dict[target_gene_id][target_isoform_id][evidence_type]["edit_distance"] = edit_distance

    
    print(parsed_evidence_gene_isoform_dict)

    return parsed_evidence_gene_isoform_dict


    