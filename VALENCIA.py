# import 
import argparse

from pathlib import Path

from src.gffread import run_gffread

#function generate parse_arguments:
def parse_arguments():
    desc = 'script to evaluated annotation with biological evidence'
    parser = argparse.ArgumentParser(description=desc)
    help_transcriptome = '(Required) Annotation with transcriptome evidence'
    parser.add_argument('--transcriptome_evidence','-t', 
                        type=str, help=help_transcriptome, required=True)
    help_protein = '(Required) Annotation with protein evidence'
    parser.add_argument('--protein_evidence', '-p',
                         type=str, help=help_protein, required=True)
    help_genome = '(Required) Genome assembly in fasta format'
    parser.add_argument('--genome_assembly', '-g',
                         type=str, help=help_genome, required=True)
    help_target = '(Required) Target annotation to be evaluated'
    parser.add_argument('--annotation_target', '-x',
                         type=str, help=help_target, required=True)
    help_outbase = '(Required) Outbase'
    parser.add_argument('--outbase','-o', 
                        type=str, help=help_outbase, required=True)
    return parser.parse_args()

# function to get arguments and return dictionary
def get_arguments():
    parser = parse_arguments()
    return {'transcripts_evidence': Path(parser.transcriptome_evidence), 
            'proteins_evidence': Path(parser.protein_evidence), 
            'genome_assembly': Path(parser.genome_assembly), 
            'annotation_target': Path(parser.annotation_target), 
            'outbase': Path(parser.outbase)}

# function main
def main():
    args = get_arguments()
    outbase = args["outbase"]
    if not outbase.exists():
        outbase.mkdir(parents=True, exist_ok=True)
    #generate_sequences for evidence
    for option, path in args.items():
        if "evidence" in option:
            kinds = ["transcripts", "proteins"]
        elif option == "annotation_target":
            kinds = ["annotation_target"]
        results = run_gffread(outbase, args["genome_assembly"],
                              args["annotation_target"], path, kinds=kinds)
        if results["returncode"] != 0:
              print("Error in {}: {}".format(option, results["log_msg"]))        

# run main function
if __name__ == '__main__':
    main()