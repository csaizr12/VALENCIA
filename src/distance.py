from Bio import SeqIO
from Levenshtein import distance

def edit_distance(gene_isoform_dict, transcript_target_fasta, transcript_evidence_fasta,
             protein_target_fasta, protein_evidence_fasta):
    # index fasta files
    records_transcript_target = SeqIO.index(transcript_target_fasta, "fasta")
    records_transcript_evidence = SeqIO.index(transcript_evidence_fasta, "fasta")
    records_protein_target = SeqIO.index(protein_target_fasta, "fasta")
    records_protein_evidence = SeqIO.index(protein_evidence_fasta, "fasta")

    # run gene_isoform_dict  
    for gene_id, isoforms in gene_isoform_dict.items():
        for iso_id, info in isoforms.items():
            if not info:
                continue
            # obtein evidence type and mach_id
            evidence_type = list(info.keys())[0]
            match_id = info[evidence_type]['match_sequence']
            # if iso_id and match_in in records, charge sequences
            if evidence_type == "transcripts":
                if iso_id not in records_transcript_target or match_id not in records_transcript_evidence:
                    continue
                seq_target = str(records_transcript_target[iso_id].seq)
                seq_evidence = str(records_transcript_evidence[match_id].seq)
            elif evidence_type == "protein":
                 if iso_id not in records_protein_target or match_id not in records_protein_evidence:
                    continue
                 seq_target = str(records_protein_target[iso_id].seq)
                 seq_evidence = str(records_protein_evidence[match_id].seq)
            else:
                continue
            # obtein edit distance with target and evidence
            edit_distance = distance(seq_target, seq_evidence)
            # save it in a gene_isoform_dict
            gene_isoform_dict[gene_id][iso_id]["edit_distance"] = edit_distance
    print(gene_isoform_dict)

    return gene_isoform_dict

    