# import 
import argparse

from pathlib import Path

#function generate parse_arguments:
def parse_arguments():
    desc = 'script to evaluated annotation with biological evidence'
    # parser argument 
    parser = argparse.ArgumentParser(description=desc)
    # generate arguments for transcriptome
    help_transcriptome = '(Requiered) Annotation with transcriptome evidence'
    # add argument for transcriptome
    parser.add_argument('--transcriptome', '-t', type=str, help=help_transcriptome, required=True)
    # generate arguments for protein evidence
    help_protein = '(Required) Annotation with protein evidence'
    # add argument for protein evidence
    parser.add_argument('--protein', '-p', type=str, help=help_protein, required=True)
    # add argument for genome assembly fasta format
    help_genome = '(Required) Genome assembly in fasta format'
    # add argument for genome assembly
    parser.add_argument('--genome', '-g', type=str, help=help_genome, required=True)
    # add argument for target annotation
    help_target = '(Required) Target annotation to be evaluated'
    # add argument for target annotation
    parser.add_argument('--target', '-x', type=str, help=help_target, required=True)
     # generate arguments for outbase
    help_outbase = '(Required) Outbase'
    # add argument for outbase
    parser.add_argument('--outbase','-o', type=str, help=help_outbase, required=True)
    # return parser arguments
    return parser.parse_args()

# function to get arguments and return diccionary
def get_arguments():
    # parse arguments
    parser = parse_arguments()
    # return diccionary with arguments
    return {'transcriptome_evidence': Path(parser.transcriptome),
            'protein_evidence': Path(parser.protein), 'genome_assembly': Path(parser.genome), 'target_annotation': Path(parser.target), 'outbase': Path(parser.outbase)}

# function main
def main():
    # get arguments
    args = get_arguments()
    # print arguments
    print(args)

# run main function
if __name__ == '__main__':
    main()