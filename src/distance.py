from Bio import SeqIO
from Levenshtein import distance

# function to calculate edit distance between target and evidence sequence, and add it to the gene_isoform_dict
def edit_distance(parsed_evidence_gene_isoform_dict, transcript_target_fasta, transcript_evidence_fasta,
                  protein_target_fasta, protein_evidence_fasta, CDS_target_fasta, CDS_evidence_fasta):
    # index fasta files
    records_transcript_target = SeqIO.index(transcript_target_fasta, "fasta")
    records_transcript_evidence = SeqIO.index(transcript_evidence_fasta, "fasta")
    records_protein_target = SeqIO.index(protein_target_fasta, "fasta")
    records_protein_evidence = SeqIO.index(protein_evidence_fasta, "fasta")
    records_CDS_target = SeqIO.index(CDS_target_fasta, "fasta")
    records_CDS_evidence = SeqIO.index(CDS_evidence_fasta, "fasta")


    # run gene_isoform_dict  
    for target_gene_id, evidence_found in parsed_evidence_gene_isoform_dict.items():
        for target_isoform_id, features in evidence_found.items():
            if not features:
                continue
            # obtein evidence type and mach_id
            evidence_types = list(features.keys())
            for evidence_type in evidence_types:
                matching_evidence_id = features[evidence_type]['match_sequence']
                 # if target_isoform_id and match_in in records, charge sequences
                if evidence_type == "transcripts":
                    seq_target = str(records_transcript_target[target_isoform_id].seq)
                    seq_evidence = str(records_transcript_evidence[matching_evidence_id].seq)
                elif evidence_type == "proteins":
                    seq_target = str(records_protein_target[target_isoform_id].seq)
                    seq_evidence = str(records_protein_evidence[matching_evidence_id].seq)
                elif evidence_type == "CDS":
                    seq_target = str(records_CDS_target[target_isoform_id].seq)
                    seq_evidence = str(records_CDS_evidence[matching_evidence_id].seq)
                # obtain edit distance with target and evidence
                len_target = len(seq_target)
                len_evidence = len(seq_evidence)
                longest_len = max(len_target, len_evidence)
                edit_distance = distance(seq_target, seq_evidence)
                shared_bases = longest_len - edit_distance
                #Adapted aed calc from ingenannot
                sensitivity = shared_bases /(shared_bases + (len_evidence - shared_bases))
                specificity = shared_bases / (shared_bases + (len_target - shared_bases))
                print(sensitivity, specificity, len_target, len_evidence, edit_distance, shared_bases)
                lev_edit_distance = 1 - ((sensitivity + specificity) /2)
                # save it in a gene_isoform_dict
                parsed_evidence_gene_isoform_dict[target_gene_id][target_isoform_id][evidence_type]["edit_distance"] = lev_edit_distance

    return parsed_evidence_gene_isoform_dict


    